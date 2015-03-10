# Author: Noel Dawe
import ROOT

import os
import math
import subprocess

import goodruns
import externaltools

# rootpy imports
from rootpy.plotting import Hist
from rootpy.tree.filtering import EventFilter, EventFilterList
from rootpy.tree import Tree, TreeChain, TreeModel, TreeBuffer
from rootpy.extern.argparse import ArgumentParser
from rootpy.io import root_open
from rootpy import stl
# Show stack traces for ROOT warning messages
#from rootpy import log
#import logging
#log["/ROOT"].show_stack(min_level=logging.WARNING)

# local imports
from higgstautau import hepmc
from higgstautau import tautools
from higgstautau import eventshapes
from higgstautau import datasets
from higgstautau import utils
from higgstautau.batch import ATLASStudent
from higgstautau.units import GeV
from higgstautau.mixins import *
from higgstautau.filters import *
from higgstautau.hadhad import branches as hhbranches
from higgstautau.hadhad.objects import define_objects
from higgstautau.hadhad.models import *
from higgstautau.hadhad.filters import *
from higgstautau import httcpana
from higgstautau import mass
from higgstautau.mass import is_MET_bisecting
from higgstautau.overlap import TauJetOverlapRemoval
from higgstautau.embedding import (
    EmbeddingPileupPatch, EmbeddingIsolation, EmbeddingCorrections)
from higgstautau.systematics import Systematics
from higgstautau.met import METRecalculation
from higgstautau.jetcalibration import JetCalibration
from higgstautau.tauspinner import EmbeddingTauSpinner, TauSpinner
from higgstautau.patches import ElectronIDpatch, TauIDpatch
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.trigger.efficiency import TauTriggerEfficiency
from higgstautau.trigger.emulation import (
    TauTriggerEmulation, update_trigger_trees)
#from higgstautau.trigger.matching import (
#    TauTriggerMatchIndex, TauTriggerMatchThreshold)
from higgstautau.pileup import (
    PileupTemplates, PileupReweight, get_pileup_reweighting_tool,
    averageIntPerXingPatch, PileupScale)
from higgstautau.rand import RandomRunNumber, RandomSeed
from higgstautau import log; log = log[__name__]

from higgstautau.hadhad.bdtscore import BDTScore

class hhskim(ATLASStudent):

    def __init__(self, options, **kwargs):
        super(hhskim, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--local', action='store_true', default=False)
        parser.add_argument('--syst-terms', default=None)
        parser.add_argument('--student-verbose', action='store_true', default=False)
        parser.add_argument('--student-very-verbose', action='store_true', default=False)
        parser.add_argument('--redo-selection', action='store_true', default=False)
        parser.add_argument('--nominal-values', action='store_true', default=False)
        args = parser.parse_args(options)
        self.args = args
        if args.syst_terms is not None:
            args.syst_terms = set([
                eval('Systematics.%s' % term) for term in
                args.syst_terms.split(',')])
        if args.local:
            def merge(inputs, output, metadata):
                # merge output trees
                root_output = output + '.root'
                log.info("merging output trees")
                subprocess.call(['hadd', root_output] + inputs)
                if metadata.datatype == datasets.DATA:
                    # merge GRLs
                    log.info("merging GRL fragments")
                    grl = goodruns.GRL()
                    for input in inputs:
                        grl |= goodruns.GRL('%s:/lumi' % input)
                    grl.save('%s:/lumi' % root_output)

            hhskim.merge = staticmethod(merge)

    def work(self):
        # get argument values
        local = self.args.local
        syst_terms = self.args.syst_terms
        datatype = self.metadata.datatype
        year = self.metadata.year
        verbose = self.args.student_verbose
        very_verbose = self.args.student_very_verbose
        redo_selection = True #self.args.redo_selection
        nominal_values = self.args.nominal_values

        # get the dataset name
        dsname = os.getenv('INPUT_DATASET_NAME', None)
        if dsname is None:
            # attempt to guess dsname from dirname
            if self.files:
                dsname = os.path.basename(os.path.dirname(self.files[0]))

        # is this a signal sample?
        # if so we will also keep some truth information in the output below
        is_signal = datatype == datasets.MC and (
            '_VBFH' in dsname or
            '_ggH' in dsname or
            '_ZH' in dsname or
            '_WH' in dsname or
            '_ttH' in dsname)
        log.info("DATASET: {0}".format(dsname))
        log.info("IS SIGNAL: {0}".format(is_signal))

        # is this an inclusive signal sample for overlap studies?
        is_inclusive_signal = is_signal and '_inclusive' in dsname

        # is this a BCH-fixed sample? (temporary)
        is_bch_sample = 'r5470_r4540_p1344' in dsname
        if is_bch_sample:
            log.warning("this is a BCH-fixed r5470 sample")

        # onfilechange will contain a list of functions to be called as the
        # chain rolls over to each new file
        onfilechange = []
        count_funcs = {}

        if datatype != datasets.DATA:
            # count the weighted number of events
            if local:
                def mc_weight_count(event):
                    return event.hh_mc_weight
            else:
                def mc_weight_count(event):
                    return event.mc_event_weight

            count_funcs = {
                'mc_weight': mc_weight_count,
            }

        # three instances of the pileup reweighting tool are created to write
        # out the nominal, high and low pileup weights
        pileup_tool = None
        pileup_tool_high = None
        pileup_tool_low = None

        if local:
            # local means running on the skims, the output of this script
            # running on the grid
            if datatype == datasets.DATA:
                # merge the GRL fragments
                merged_grl = goodruns.GRL()

                def update_grl(student, grl, name, file, tree):
                    grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())

                onfilechange.append((update_grl, (self, merged_grl,)))

            if datatype == datasets.DATA:
                merged_cutflow = Hist(1, 0, 1, name='cutflow', type='D')
            else:
                merged_cutflow = Hist(2, 0, 2, name='cutflow', type='D')

            def update_cutflow(student, cutflow, name, file, tree):
                # record a cut-flow
                year = student.metadata.year
                datatype = student.metadata.datatype
                cutflow[1].value += file.cutflow_event[1].value
                if datatype != datasets.DATA:
                    cutflow[2].value += file.cutflow_event_mc_weight[1].value

            onfilechange.append((update_cutflow, (self, merged_cutflow,)))

        else:
            # get pileup reweighting tool
            pileup_tool = get_pileup_reweighting_tool(
                year=year,
                use_defaults=True)
            pileup_tool_high = get_pileup_reweighting_tool(
                year=year,
                use_defaults=True,
                systematic='high')
            pileup_tool_low = get_pileup_reweighting_tool(
                year=year,
                use_defaults=True,
                systematic='low')

            if datatype not in (datasets.EMBED, datasets.MCEMBED):
                # merge TrigConfTrees
                metadirname = '%sMeta' % self.metadata.treename
                trigconfchain = ROOT.TChain('%s/TrigConfTree' % metadirname)
                map(trigconfchain.Add, self.files)
                metadir = self.output.mkdir(metadirname)
                metadir.cd()
                trigconfchain.Merge(self.output, -1, 'fast keep')
                self.output.cd()

            if datatype == datasets.DATA:
                # merge GRL XML strings
                merged_grl = goodruns.GRL()
                for fname in self.files:
                    with root_open(fname) as f:
                        for key in f.Lumi.keys():
                            merged_grl |= goodruns.GRL(
                                str(key.ReadObj().GetString()),
                                from_string=True)
                lumi_dir = self.output.mkdir('Lumi')
                lumi_dir.cd()
                xml_string= ROOT.TObjString(merged_grl.str())
                xml_string.Write(self.metadata.treename)
                self.output.cd()

        self.output.cd()

        # create the output tree
        model = get_model(datatype, dsname,
                          prefix=None if local else 'hh_',
                          is_inclusive_signal=is_inclusive_signal)
        log.info("Output Model:\n\n{0}\n\n".format(model))
        outtree = Tree(name=self.metadata.treename,
                       model=model)

        if local:
            tree = outtree
        else:
            tree = outtree.define_object(name='tree', prefix='hh_')

        tree.define_object(name='tau', prefix='tau_')
        tree.define_object(name='tau1', prefix='tau1_')
        tree.define_object(name='tau2', prefix='tau2_')
        tree.define_object(name='truetau1', prefix='truetau1_')
        tree.define_object(name='truetau2', prefix='truetau2_')
        tree.define_object(name='jet1', prefix='jet1_')
        tree.define_object(name='jet2', prefix='jet2_')
        tree.define_object(name='jet3', prefix='jet3_')

        mmc_objects = [
            tree.define_object(name='mmc0', prefix='mmc0_'),
            tree.define_object(name='mmc1', prefix='mmc1_'),
            tree.define_object(name='mmc2', prefix='mmc2_'),
        ]

        for mmc_obj in mmc_objects:
            mmc_obj.define_object(name='resonance', prefix='resonance_')
            mmc_obj.define_object(name='tau1', prefix='tau1_')
            mmc_obj.define_object(name='tau2', prefix='tau2_')
            mmc_obj.define_object(name='nu1', prefix='nu1_')
            mmc_obj.define_object(name='nu2', prefix='nu2_')
            mmc_obj.define_object(name='vistau1', prefix='vistau1_')
            mmc_obj.define_object(name='vistau2', prefix='vistau2_')



        trigger_emulation = TauTriggerEmulation(
            year=year,
            passthrough=local or datatype != datasets.MC or year > 2011,
            count_funcs=count_funcs)

        if not trigger_emulation.passthrough:
            onfilechange.append(
                (update_trigger_trees, (self, trigger_emulation,)))

        trigger_config = None

        if datatype not in (datasets.EMBED, datasets.MCEMBED):
            # trigger config tool to read trigger info in the ntuples
            trigger_config = get_trigger_config()
            # update the trigger config maps on every file change
            onfilechange.append((update_trigger_config, (trigger_config,)))

        # define the list of event filters
        if local and syst_terms is None and not redo_selection:
            event_filters = None
        else:
            tau_ntrack_recounted_use_ntup = False
            if year > 2011:
                # peek at first tree to determine if the extended number of
                # tracks is already stored
                with root_open(self.files[0]) as test_file:
                    test_tree = test_file.Get(self.metadata.treename)
                    tau_ntrack_recounted_use_ntup = (
                        'tau_out_track_n_extended' in test_tree)

            event_filters = EventFilterList([
                GRLFilter(
                    self.grl,
                    passthrough=(
                        local or (
                            datatype not in (datasets.DATA, datasets.EMBED))),
                    count_funcs=count_funcs),
                CoreFlags(
                    passthrough=local,
                    count_funcs=count_funcs),
                EmbeddingPileupPatch(
                    passthrough=(
                        local or year > 2011 or datatype != datasets.EMBED),
                    count_funcs=count_funcs),
                averageIntPerXingPatch(
                    passthrough=(
                        local or year < 2012 or datatype != datasets.MC),
                    count_funcs=count_funcs),
                PileupTemplates(
                    year=year,
                    passthrough=(
                        local or is_bch_sample or datatype not in (
                            datasets.MC, datasets.MCEMBED)),
                    count_funcs=count_funcs),
                RandomSeed(
                    datatype=datatype,
                    count_funcs=count_funcs),
                BCHSampleRunNumber(
                    passthrough=not is_bch_sample,
                    count_funcs=count_funcs),
                RandomRunNumber(
                    tree=tree,
                    datatype=datatype,
                    pileup_tool=pileup_tool,
                    passthrough=local,
                    count_funcs=count_funcs),
                trigger_emulation,
                Triggers(
                    year=year,
                    tree=tree,
                    datatype=datatype,
                    passthrough=datatype in (datasets.EMBED, datasets.MCEMBED),
                    count_funcs=count_funcs),
                PileupReweight(
                    year=year,
                    tool=pileup_tool,
                    tool_high=pileup_tool_high,
                    tool_low=pileup_tool_low,
                    tree=tree,
                    passthrough=(
                        local or (
                            datatype not in (datasets.MC, datasets.MCEMBED))),
                    count_funcs=count_funcs),
                PriVertex(
                    passthrough=local,
                    count_funcs=count_funcs),
                LArError(
                    passthrough=local,
                    count_funcs=count_funcs),
                TileError(
                    passthrough=local,
                    count_funcs=count_funcs),
                TileTrips(
                    passthrough=(
                        local or datatype in (datasets.MC, datasets.MCEMBED)),
                    count_funcs=count_funcs),
                JetCopy(
                    tree=tree,
                    passthrough=local,
                    count_funcs=count_funcs),
                # IMPORTANT!
                # JetCalibration MUST COME BEFORE ANYTHING THAT REFERS TO
                # jet.fourvect since jet.fourvect IS CACHED!
                JetCalibration(
                    datatype=datatype,
                    year=year,
                    verbose=very_verbose,
                    passthrough=local or nominal_values,
                    count_funcs=count_funcs),
                # in situ TES shift for 2012 data
                TauEnergyShift(
                    passthrough=(
                        local or datatype != datasets.DATA
                        or year < 2012 or nominal_values),
                    count_funcs=count_funcs),
                # truth matching must come before systematics due to
                # TES_TRUE/FAKE
                TruthMatching(
                    passthrough=datatype == datasets.DATA,
                    count_funcs=count_funcs),
                NvtxJets(
                    tree=tree,
                    count_funcs=count_funcs),
                # PUT THE SYSTEMATICS "FILTER" BEFORE
                # ANY FILTERS THAT REFER TO OBJECTS
                # BUT AFTER CALIBRATIONS
                # Systematics must also come before anything that refers to
                # thing.fourvect since fourvect is cached!
                Systematics(
                    terms=syst_terms,
                    year=year,
                    datatype=datatype,
                    tree=tree,
                    verbose=verbose,
                    passthrough=not syst_terms,
                    count_funcs=count_funcs),
                JetIsPileup(
                    passthrough=(
                        local or year < 2012 or
                        datatype not in (datasets.MC, datasets.MCEMBED)),
                    count_funcs=count_funcs),
                LArHole(
                    tree=tree,
                    passthrough=year > 2011,
                    count_funcs=count_funcs),
                JetCleaning(
                    datatype=datatype,
                    year=year,
                    count_funcs=count_funcs),
                ElectronVeto(
                    count_funcs=count_funcs),
                MuonVeto(
                    year=year,
                    count_funcs=count_funcs),
                TauPT(2,
                    thresh=20 * GeV,
                    count_funcs=count_funcs),
                TauHasTrack(2,
                    count_funcs=count_funcs),
                TauEta(2,
                    count_funcs=count_funcs),
                TauElectronVeto(2,
                    count_funcs=count_funcs),
                TauMuonVeto(2,
                    count_funcs=count_funcs),
                TauAuthor(2,
                    count_funcs=count_funcs),
                TauCrack(2,
                    count_funcs=count_funcs),
                TauLArHole(2,
                    tree=tree,
                    passthrough=year > 2011,
                    count_funcs=count_funcs),
                # before selecting the leading and subleading taus
                # be sure to only consider good candidates
                TauIDMedium(2,
                    count_funcs=count_funcs),
                #TauTriggerMatchIndex(
                #    config=trigger_config,
                #    year=year,
                #    datatype=datatype,
                #    passthrough=datatype == datasets.EMBED,
                #    count_funcs=count_funcs),
                # Select two leading taus at this point
                # 25 and 35 for data
                # 20 and 30 for MC to leave room for TES uncertainty
                TauLeadSublead(
                    lead=(
                        35 * GeV if datatype == datasets.DATA or local
                        else 30 * GeV),
                    sublead=(
                        25 * GeV if datatype == datasets.DATA or local
                        else 20 * GeV),
                    count_funcs=count_funcs),
                # taus are sorted (in decreasing order) by pT from here on
                TauIDSelection(
                    tree=tree,
                    count_funcs=count_funcs),
                TaudR(3.2,
                    count_funcs=count_funcs),
                #TauTriggerMatchThreshold(
                #    datatype=datatype,
                #    tree=tree,
                #    count_funcs=count_funcs),
                TauTriggerEfficiency(
                    year=year,
                    datatype=datatype,
                    tree=tree,
                    tes_systematic=self.args.syst_terms and (
                        Systematics.TES_TERMS & self.args.syst_terms),
                    passthrough=datatype == datasets.DATA,
                    count_funcs=count_funcs),
                PileupScale(
                    tree=tree,
                    year=year,
                    datatype=datatype,
                    passthrough=local,
                    count_funcs=count_funcs),
                TauIDScaleFactors(
                    year=year,
                    passthrough=datatype == datasets.DATA,
                    count_funcs=count_funcs),
                TauFakeRateScaleFactors(
                    year=year,
                    datatype=datatype,
                    tree=tree,
                    tes_up=(self.args.syst_terms is not None and
                        (Systematics.TES_FAKE_TOTAL_UP in self.args.syst_terms or
                         Systematics.TES_FAKE_FINAL_UP in self.args.syst_terms)),
                    tes_down=(self.args.syst_terms is not None and
                        (Systematics.TES_FAKE_TOTAL_DOWN in self.args.syst_terms or
                         Systematics.TES_FAKE_FINAL_DOWN in self.args.syst_terms)),
                    passthrough=datatype in (datasets.DATA, datasets.EMBED),
                    count_funcs=count_funcs),
                HiggsPT(
                    year=year,
                    tree=tree,
                    passthrough=not is_signal or local,
                    count_funcs=count_funcs),
                TauTrackRecounting(
                    year=year,
                    use_ntup_value=tau_ntrack_recounted_use_ntup,
                    passthrough=local,
                    count_funcs=count_funcs),
                MCWeight(
                    datatype=datatype,
                    tree=tree,
                    passthrough=local or datatype == datasets.DATA,
                    count_funcs=count_funcs),
                EmbeddingIsolation(
                    tree=tree,
                    passthrough=(
                        local or year < 2012 or
                        datatype not in (datasets.EMBED, datasets.MCEMBED)),
                    count_funcs=count_funcs),
                EmbeddingCorrections(
                    tree=tree,
                    year=year,
                    passthrough=(
                        local or
                        datatype not in (datasets.EMBED, datasets.MCEMBED)),
                    count_funcs=count_funcs),
                TauSpinner(
                     year=year,
                     tree=tree,
                     passthrough=(
                        local or
                        not (is_signal or
                             datatype in (datasets.EMBED, datasets.MCEMBED))),
                     count_funcs=count_funcs),
                # put MET recalculation after tau selection but before tau-jet
                # overlap removal and jet selection because of the RefAntiTau
                # MET correction
                METRecalculation(
                    terms=syst_terms,
                    year=year,
                    tree=tree,
                    refantitau=not nominal_values,
                    verbose=verbose,
                    very_verbose=very_verbose,
                    count_funcs=count_funcs),
                TauJetOverlapRemoval(
                    count_funcs=count_funcs),
                JetPreselection(
                    count_funcs=count_funcs),
                NonIsolatedJet(
                    tree=tree,
                    count_funcs=count_funcs),
                JetSelection(
                    year=year,
                    count_funcs=count_funcs),
                RecoJetTrueTauMatching(
                    passthrough=datatype == datasets.DATA or local,
                    count_funcs=count_funcs),
                BCHCleaning(
                    tree=tree,
                    passthrough=year == 2011 or local,
                    datatype=datatype,
                    count_funcs=count_funcs),
                ClassifyInclusiveHiggsSample(
                    tree=tree,
                    passthrough=not is_inclusive_signal,
                    count_funcs=count_funcs),
            ])

            # set the event filters
            self.filters['event'] = event_filters

        # peek at first tree to determine which branches to exclude
        with root_open(self.files[0]) as test_file:
            test_tree = test_file.Get(self.metadata.treename)
            ignore_branches = test_tree.glob(
                hhbranches.REMOVE,
                exclude=hhbranches.KEEP)
            ignore_branches_output = test_tree.glob(
                hhbranches.REMOVE_OUTPUT,
                exclude=hhbranches.KEEP_OUTPUT)

        # initialize the TreeChain of all input files
        chain = TreeChain(
            self.metadata.treename,
            files=self.files,
            ignore_branches=ignore_branches,
            events=self.events,
            onfilechange=onfilechange,
            filters=event_filters,
            cache=True,
            cache_size=50000000,
            learn_entries=100)

        if local:
            copied = [
                'EventNumber',
            ]

            hh_buffer = TreeBuffer()
            buffer = TreeBuffer()
            for name, value in chain._buffer.items():
                if name.startswith('hh_'):
                    hh_buffer[name[3:]] = value
                elif name in copied:
                    buffer[name] = value
            outtree.set_buffer(
                hh_buffer,
                create_branches=False,
                visible=True)
            outtree.set_buffer(
                buffer,
                create_branches=True,
                visible=False)

        else:
            # additional decorations on existing objects
            if year > 2011 and datatype in (datasets.MC, datasets.MCEMBED):
                class Decorations(TreeModel):
                    jet_ispileup = stl.vector('bool')

                chain.set_buffer(Decorations(), create_branches=True)

            # include the branches in the input chain in the output tree
            # set branches to be removed in ignore_branches
            outtree.set_buffer(
                chain._buffer,
                ignore_branches=ignore_branches + ignore_branches_output,
                create_branches=True,
                ignore_duplicates=True,
                transfer_objects=True,
                visible=False)

        # define tree objects
        define_objects(chain, year)

        # create the BDT classifier
        bdt = BDTScore(year=year)

        # create the HttCP
        cp = httcpana.HttCP()

        # create the MMC
        mmc = mass.MMC(year=year)
        # report which packages have been loaded
        externaltools.report()

        self.output.cd()

        # The main event loop
        # the event filters above are automatically run for each event and only
        # the surviving events are looped on
        for event in chain:

            if local and syst_terms is None and not redo_selection:
                outtree.Fill()
                continue

            # sort taus and jets in decreasing order by pT
            event.taus.sort(key=lambda tau: tau.pt, reverse=True)
            event.jets.sort(key=lambda jet: jet.pt, reverse=True)

            # tau1 is the leading tau
            # tau2 is the subleading tau
            tau1, tau2 = event.taus
            jets = list(event.jets)
            jet1, jet2, jet3 = None, None, None
            beta = None

            # leading track
            tau1_lead_track_tlv = ROOT.TLorentzVector()
            tau1_lead_track_qoverp = -1E10
            tau1_lead_track_d0 = -1E10
            tau1_lead_track_z0 = -1E10

            tau1_lead_track_tlv, tau1_lead_track_qoverp, tau1_lead_track_d0, tau1_lead_track_z0 = cp.find_lead_track ( tau1 )
            tau1.lead_track_pt = tau1_lead_track_tlv.Pt()
            tau1.lead_track_eta = tau1_lead_track_tlv.Eta()
            tau1.lead_track_phi = tau1_lead_track_tlv.Phi()
            tau1_lead_track_qoverp = tau1.lead_track_qoverp
            tau1_lead_track_d0 = tau1.lead_track_d0
            tau1_lead_track_z0 = tau1.lead_track_z0

            tau2_lead_track_tlv = ROOT.TLorentzVector()
            tau2_lead_track_qoverp = -1E10
            tau2_lead_track_d0 = -1E10
            tau2_lead_track_z0 = -1E10

            tau2_lead_track_tlv, tau2_lead_track_qoverp, tau2_lead_track_d0, tau2_lead_track_z0 = cp.find_lead_track ( tau2 )
            tau2.lead_track_pt = tau2_lead_track_tlv.Pt()
            tau2.lead_track_eta = tau2_lead_track_tlv.Eta()
            tau2.lead_track_phi = tau2_lead_track_tlv.Phi()
            tau2_lead_track_qoverp = tau2.lead_track_qoverp
            tau2_lead_track_d0 = tau2.lead_track_d0
            tau2_lead_track_z0 = tau2.lead_track_z0


            if len(jets) >= 2:
                jet1, jet2 = jets[:2]

                # determine boost of system
                # determine jet CoM frame
                beta = (jet1.fourvect + jet2.fourvect).BoostVector()
                tree.jet_beta.copy_from(beta)

                jet1.fourvect_boosted.copy_from(jet1.fourvect)
                jet2.fourvect_boosted.copy_from(jet2.fourvect)
                jet1.fourvect_boosted.Boost(beta * -1)
                jet2.fourvect_boosted.Boost(beta * -1)

                tau1.fourvect_boosted.copy_from(tau1.fourvect)
                tau2.fourvect_boosted.copy_from(tau2.fourvect)
                tau1.fourvect_boosted.Boost(beta * -1)
                tau2.fourvect_boosted.Boost(beta * -1)

                tau1.min_dr_jet = min(
                    tau1.fourvect.DeltaR(jet1.fourvect),
                    tau1.fourvect.DeltaR(jet2.fourvect))
                tau2.min_dr_jet = min(
                    tau2.fourvect.DeltaR(jet1.fourvect),
                    tau2.fourvect.DeltaR(jet2.fourvect))

                #sphericity, aplanarity = eventshapes.sphericity_aplanarity(
                #    [tau1.fourvect,
                #     tau2.fourvect,
                #     jet1.fourvect,
                #     jet2.fourvect])

                # sphericity
                #tree.sphericity = sphericity
                # aplanarity
                #tree.aplanarity = aplanarity

                #sphericity_boosted, aplanarity_boosted = eventshapes.sphericity_aplanarity(
                #    [tau1.fourvect_boosted,
                #     tau2.fourvect_boosted,
                #     jet1.fourvect_boosted,
                #     jet2.fourvect_boosted])

                # sphericity
                #tree.sphericity_boosted = sphericity_boosted
                # aplanarity
                #tree.aplanarity_boosted = aplanarity_boosted

                # tau centrality (degree to which they are between the two jets)
                tau1.centrality = eventshapes.eta_centrality(
                    tau1.fourvect.Eta(),
                    jet1.fourvect.Eta(),
                    jet2.fourvect.Eta())

                tau2.centrality = eventshapes.eta_centrality(
                    tau2.fourvect.Eta(),
                    jet1.fourvect.Eta(),
                    jet2.fourvect.Eta())

                # boosted tau centrality
                tau1.centrality_boosted = eventshapes.eta_centrality(
                    tau1.fourvect_boosted.Eta(),
                    jet1.fourvect_boosted.Eta(),
                    jet2.fourvect_boosted.Eta())

                tau2.centrality_boosted = eventshapes.eta_centrality(
                    tau2.fourvect_boosted.Eta(),
                    jet1.fourvect_boosted.Eta(),
                    jet2.fourvect_boosted.Eta())


                # 3rd leading jet
                if len(jets) >= 3:
                    jet3 = jets[2]
                    jet3.fourvect_boosted.copy_from(jet3.fourvect)
                    jet3.fourvect_boosted.Boost(beta * -1)

            elif len(jets) == 1:
                jet1 = jets[0]

                tau1.min_dr_jet = tau1.fourvect.DeltaR(jet1.fourvect)
                tau2.min_dr_jet = tau2.fourvect.DeltaR(jet1.fourvect)

                #sphericity, aplanarity = eventshapes.sphericity_aplanarity(
                #    [tau1.fourvect,
                #     tau2.fourvect,
                #     jet1.fourvect])

                # sphericity
                #tree.sphericity = sphericity
                # aplanarity
                #tree.aplanarity = aplanarity

            RecoJetBlock.set(tree, jet1, jet2, jet3, local=local)

            # mass of ditau + leading jet system
            if jet1 is not None:
                tree.mass_tau1_tau2_jet1 = (
                    tau1.fourvect + tau2.fourvect + jet1.fourvect).M()

            # full sphericity and aplanarity
            #sphericity_full, aplanarity_full = eventshapes.sphericity_aplanarity(
            #    [tau1.fourvect, tau2.fourvect] + [jet.fourvect for jet in jets])

            #tree.sphericity_full = sphericity_full
            #tree.aplanarity_full = aplanarity_full

            #####################################
            # number of tracks from PV minus taus
            #####################################
            ntrack_pv = 0
            ntrack_nontau_pv = 0
            privtx = ROOT.TVector3()
            for vxp in event.vertices:
                # primary vertex
                if vxp.type == 1:
                    ntrack_pv = vxp.nTracks
                    ntrack_nontau_pv = ntrack_pv - tau1.numTrack - tau2.numTrack
                    break
                if vpx.type == 1 and vxp.nTrack >= 4:
                    privtx.SetXYZ(vxp.x,vxp.y,vxp.z)
                    break
            tree.ntrack_pv = ntrack_pv
            tree.ntrack_nontau_pv = ntrack_nontau_pv

            #########################
            # MET variables
            #########################
            METx = event.MET.etx
            METy = event.MET.ety
            MET = event.MET.et
            MET_vect = Vector2(METx, METy)
            MET_4vect = LorentzVector()
            MET_4vect.SetPxPyPzE(METx, METy, 0., MET)
            MET_4vect_boosted = LorentzVector()
            MET_4vect_boosted.copy_from(MET_4vect)
            if beta is not None:
                MET_4vect_boosted.Boost(beta * -1)

            tree.MET_et = MET
            tree.MET_etx = METx
            tree.MET_ety = METy
            tree.MET_phi = event.MET.phi
            dPhi_tau1_tau2 = abs(tau1.fourvect.DeltaPhi(tau2.fourvect))
            dPhi_tau1_MET = abs(tau1.fourvect.DeltaPhi(MET_4vect))
            dPhi_tau2_MET = abs(tau2.fourvect.DeltaPhi(MET_4vect))
            tree.dPhi_tau1_tau2 = dPhi_tau1_tau2
            tree.dPhi_tau1_MET = dPhi_tau1_MET
            tree.dPhi_tau2_MET = dPhi_tau2_MET
            tree.dPhi_min_tau_MET = min(dPhi_tau1_MET, dPhi_tau2_MET)
            tree.MET_bisecting = is_MET_bisecting(
                dPhi_tau1_tau2,
                dPhi_tau1_MET,
                dPhi_tau2_MET)

            sumET = event.MET.sumet
            tree.MET_sumet = sumET
            if sumET != 0:
                tree.MET_sig = ((2. * MET / GeV) /
                    (utils.sign(sumET) * sqrt(abs(sumET / GeV))))
            else:
                tree.MET_sig = -1.

            tree.MET_centrality = eventshapes.phi_centrality(
                tau1.fourvect,
                tau2.fourvect,
                MET_vect)
            tree.MET_centrality_boosted = eventshapes.phi_centrality(
                tau1.fourvect_boosted,
                tau2.fourvect_boosted,
                MET_4vect_boosted)

            tree.number_of_good_vertices = len(event.vertices)

            ##########################
            # Jet and sum pt variables
            ##########################
            tree.numJets = len(event.jets)

            # sum pT with only the two leading jets
            tree.sum_pt = sum(
                [tau1.pt, tau2.pt] +
                [jet.pt for jet in jets[:2]])

            # sum pT with all selected jets
            tree.sum_pt_full = sum(
                [tau1.pt, tau2.pt] +
                [jet.pt for jet in jets])

            # vector sum pT with two leading jets and MET
            tree.vector_sum_pt = sum(
                [tau1.fourvect, tau2.fourvect] +
                [jet.fourvect for jet in jets[:2]] +
                [MET_4vect]).Pt()

            # vector sum pT with all selected jets and MET
            tree.vector_sum_pt_full = sum(
                [tau1.fourvect, tau2.fourvect] +
                [jet.fourvect for jet in jets] +
                [MET_4vect]).Pt()

            # resonance pT
            tree.resonance_pt = sum(
                [tau1.fourvect, tau2.fourvect, MET_4vect]).Pt()

            #############################
            # tau <-> vertex association
            #############################
            tree.tau_same_vertex = (
                tau1.privtx_x == tau2.privtx_x and
                tau1.privtx_y == tau2.privtx_y and
                tau1.privtx_z == tau2.privtx_z)

            tau1.vertex_prob = ROOT.TMath.Prob(
                tau1.privtx_chiSquared,
                int(tau1.privtx_numberDoF))

            tau2.vertex_prob = ROOT.TMath.Prob(
                tau2.privtx_chiSquared,
                int(tau2.privtx_numberDoF))


            ##########################
            # MMC Mass
            ##########################

            mmc_result = mmc.mass(
                tau1, tau2,
                METx, METy, sumET,
                njets=len(event.jets),
                EventNumber=event.EventNumber)

            mmc1_mass = 0.0
            for mmc_method, mmc_object in enumerate(mmc_objects):
                mmc_mass, mmc_resonance, mmc_met, \
                    mmc_tau1, mmc_tau2, \
                    mmc_nu1, mmc_nu2, \
                    mmc_vistau1, mmc_vistau2 = mmc_result[mmc_method]

                if verbose:
                    log.info("MMC (method %d): %f" % (mmc_method, mmc_mass))
                if mmc_method == 1:
                    mmc1_mass = mmc_mass

                mmc_object.mass = mmc_mass
                mmc_object.MET_et = mmc_met.Mod()
                mmc_object.MET_etx = mmc_met.X()
                mmc_object.MET_ety = mmc_met.Y()
                mmc_object.MET_phi = math.pi - mmc_met.Phi()
                if mmc_mass > 0:
                    FourMomentum.set(mmc_object.resonance, mmc_resonance)
                    FourMomentum.set(mmc_object.tau1, mmc_tau1)
                    FourMomentum.set(mmc_object.tau2, mmc_tau2)
                    FourMomentum.set(mmc_object.nu1, mmc_nu1)
                    FourMomentum.set(mmc_object.nu2, mmc_nu2)
                    FourMomentum.set(mmc_object.vistau1, mmc_vistau1)
                    FourMomentum.set(mmc_object.vistau2, mmc_vistau2)


            ############################
            # collinear and visible mass
            ############################
            vis_mass, collin_mass, tau1_x, tau2_x = mass.collinearmass(
                tau1, tau2, METx, METy)

            tree.mass_vis_tau1_tau2 = vis_mass
            tree.mass_collinear_tau1_tau2 = collin_mass
            tau1.collinear_momentum_fraction = tau1_x
            tau2.collinear_momentum_fraction = tau2_x

            ###########################
            # Match jets to VBF partons
            ###########################
            #if datatype == datasets.MC and 'VBF' in dsname and year == 2011:
            #    # get partons (already sorted by eta in hepmc) FIXME!!!
            #    parton1, parton2 = hepmc.get_VBF_partons(event)
            #    tree.mass_true_quark1_quark2 = (parton1.fourvect + parton2.fourvect).M()
            #    # order here needs to be revised since jets are no longer
            #    # sorted by eta but instead by pT
            #    PartonBlock.set(tree, parton1, parton2)
            #    if len(jets) >= 2:
            #        jet1, jet2 = jets[:2]
            #        for i, jet in zip((1, 2), (jet1, jet2)):
            #            for parton in (parton1, parton2):
            #                if utils.dR(jet.eta, jet.phi, parton.eta, parton.phi) < .8:
            #                    setattr(tree, 'jet%i_matched' % i, True)


            ############################
            # MVA BDTScore
            ############################
            # calculate BDTScore

            inputs_vbf = [mmc1_mass,
                          tree.dEta_jets,
                          tree.eta_product_jets,
                          tree.mass_jet1_jet2,
                          tree.tau1_centrality,
                          tree.tau2_centrality,
                          tree.dR_tau1_tau2,
                          tree.MET_centrality,
                          tree.vector_sum_pt]

            inputs_boosted = [mmc1_mass,
                              tree.dR_tau1_tau2,
                              tree.MET_centrality,
                              tree.sum_pt_full,
                              tree.tau_pt_ratio]

            tree.BDTScore_vbf_0 = bdt.get_score(event.EventNumber, 0, 'VBF', inputs_vbf)
            tree.BDTScore_boosted_0 = bdt.get_score(event.EventNumber, 0, 'Boosted', inputs_boosted)
            tree.BDTScore_vbf_1 = bdt.get_score(event.EventNumber, 1, 'VBF', inputs_vbf)
            tree.BDTScore_boosted_1 = bdt.get_score(event.EventNumber, 1, 'Boosted', inputs_boosted)

            ############################
            # CP angle variables
            ############################
            # calculate Acoplanarity angles

            # for pi-pi decay
            tree.Acoplanarity_IP = cp.Acoplanarity_IP( tau1, tau2, event.vertices )

            # for pi-rho , rho-rho decay : leadtrk
            tree.Acoplanarity_tau1_IP_tau2_rho_leadtrk = \
                cp.Acoplanarity_tau1_IP_tau2_rho_leadtrk( tau1, tau2, event.vertices )
            tree.Acoplanarity_tau2_IP_tau1_rho_leadtrk = \
                cp.Acoplanarity_tau2_IP_tau1_rho_leadtrk( tau1, tau2, event.vertices )
            tree.Acoplanarity_rho_leadtrk = cp.Acoplanarity_rho_leadtrk( tau1, tau2 )

            if year > 2011:
                # for pi-rho , rho-rho decay : cluster
                tree.Acoplanarity_tau1_IP_tau2_rho_cluster = \
                    cp.Acoplanarity_tau1_IP_tau2_rho_cluster( tau1, tau2, event.vertices )
                tree.Acoplanarity_tau2_IP_tau1_rho_cluster = \
                    cp.Acoplanarity_tau2_IP_tau1_rho_cluster( tau1, tau2, event.vertices )
                tree.Acoplanarity_rho_cluster = cp.Acoplanarity_rho_cluster( tau1, tau2 )

            # use mmc output information
            for mmc_method, mmc_object in enumerate(mmc_objects):
                if mmc_method == 0 :
                    tree.Acoplanarity_mmc0 = \
                        cp.Acoplanarity_mmc(mmc_object.vistau1,
                                            mmc_object.nu1,
                                            mmc_object.vistau2,
                                            mmc_object.nu2)
                    tree.Acoplanarity_mmc0_visboost = \
                        cp.Acoplanarity_mmc_visboost(mmc_object.vistau1,
                                                     mmc_object.nu1,
                                                     mmc_object.vistau2,
                                                     mmc_object.nu2)

                if mmc_method == 1 :
                    tree.Acoplanarity_mmc1 = \
                        cp.Acoplanarity_mmc(mmc_object.vistau1,
                                            mmc_object.nu1,
                                            mmc_object.vistau2,
                                            mmc_object.nu2)
                    tree.Acoplanarity_mmc1_visboost = \
                        cp.Acoplanarity_mmc_visboost(mmc_object.vistau1,
                                                     mmc_object.nu1,
                                                     mmc_object.vistau2,
                                                     mmc_object.nu2)

                if mmc_method == 2 :
                    tree.Acoplanarity_mmc2 = \
                        cp.Acoplanarity_mmc(mmc_object.vistau1,
                                            mmc_object.nu1,
                                            mmc_object.vistau2,
                                            mmc_object.nu2)
                    tree.Acoplanarity_mmc2_visboost = \
                        cp.Acoplanarity_mmc_visboost(mmc_object.vistau1,
                                                     mmc_object.nu1,
                                                     mmc_object.vistau2,
                                                     mmc_object.nu2)

            # Fill the tau block
            # This must come after the RecoJetBlock is filled since
            # that sets the jet_beta for boosting the taus
            RecoTauBlock.set(event, tree, datatype, tau1, tau2, local=local)
            if datatype != datasets.DATA:
                TrueTauBlock.set(tree, tau1, tau2)

            # fill the output tree
            outtree.Fill(reset=True)

        externaltools.report()

        # flush any baskets remaining in memory to disk
        self.output.cd()
        outtree.FlushBaskets()
        outtree.Write()

        if local:
            if datatype == datasets.DATA:
                xml_string = ROOT.TObjString(merged_grl.str())
                xml_string.Write('lumi')
            merged_cutflow.Write()
