import ROOT
import math

from argparse import ArgumentParser

from rootpy.tree.filtering import *
from rootpy.tree import Tree, TreeBuffer, TreeChain
from rootpy.plotting import Hist
from rootpy.io import open as ropen

from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from atlastools.filtering import GRLFilter
from atlastools.batch import ATLASStudent

from higgstautau.mixins import *
from higgstautau import hepmc
from higgstautau import tautools
from higgstautau.models import *
from higgstautau.hadhad.models import *
from higgstautau import eventshapes
from higgstautau import eventview
from higgstautau.filters import *
from higgstautau.hadhad.filters import *
from higgstautau.hadhad.categories import *
from higgstautau import mass
#from higgstautau.mass.ditaumass import HAD1P, HAD3P
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.systematics import Systematics
from higgstautau.jetcalibration import JetCalibration
from higgstautau.overlap import TauJetOverlapRemoval

from goodruns import GRL
import subprocess

from externaltools import TauFakeRates
from ROOT import TauFakeRates as TFR

#ROOT.gErrorIgnoreLevel = ROOT.kFatal
YEAR = 2011
VERBOSE = False


class SkimExtraModel(TreeModel):

    number_of_good_vertices = IntCol()


class TauIDProcessor(ATLASStudent):

    def work(self):

        # trigger config tool to read trigger info in the ntuples
        trigger_config = get_trigger_config()

        OutputModel = (RecoTauBlock + EventVariables + SkimExtraModel +
            TrueTauBlock)

        onfilechange = []
        # update the trigger config maps on every file change
        onfilechange.append((update_trigger_config, (trigger_config,)))

        cutflow = Hist(2, 0, 2, name='cutflow', type='D')

        # initialize the TreeChain of all input files (each containing one tree named self.metadata.treename)
        chain = TreeChain(self.metadata.treename,
                         files=self.files,
                         events=self.events,
                         cache=True,
                         cache_size=10000000,
                         learn_entries=30,
                         onfilechange=onfilechange)

        # create output tree
        self.output.cd()
        tree = Tree(name='higgstautauhh', model=OutputModel)

        copied_variables = ['actualIntPerXing',
                            'averageIntPerXing',
                            'RunNumber',
                            'EventNumber',
                            'lbn']

        tree.set_buffer(
                chain.buffer,
                branches=copied_variables,
                create_branches=True,
                visible=False)
        chain.always_read(copied_variables)

        # set the event filters
        event_filters = EventFilterList([
            #Triggers(
            #    datatype=self.metadata.datatype,
            #    year=YEAR,
            #    skim=False),
            PriVertex(),
            LArError(),
            LArHole(datatype=self.metadata.datatype),
            JetCleaning(
                datatype=self.metadata.datatype,
                year=YEAR),
            TauAuthor(1),
            TauHasTrack(1),
            TauPT(1, thresh=25 * GeV),
            TauEta(1),
            TauCrack(1),
            TauLArHole(1),
            #TauTriggerMatch(
            #    config=trigger_config,
            #    year=YEAR,
            #    datatype=self.metadata.datatype,
            #    skim=False,
            #    tree=tree,
            #    min_taus=1),
        ])

        self.filters['event'] = event_filters
        chain.filters += event_filters

        # define tree collections
        chain.define_collection(name="taus", prefix="tau_", size="tau_n", mix=TauFourMomentum)
        chain.define_collection(name="taus_EF", prefix="trig_EF_tau_",
                                size="trig_EF_tau_n", mix=TauFourMomentum)

        # jet_* etc. is AntiKt4LCTopo_* in tau-perf D3PDs
        chain.define_collection(name="jets", prefix="jet_", size="jet_n", mix=FourMomentum)
        chain.define_collection(name="truetaus", prefix="trueTau_", size="trueTau_n", mix=MCTauFourMomentum)
        chain.define_collection(name="mc", prefix="mc_", size="mc_n", mix=MCParticle)
        chain.define_collection(name="muons", prefix="mu_staco_", size="mu_staco_n")
        chain.define_collection(name="electrons", prefix="el_", size="el_n")
        chain.define_collection(name="vertices", prefix="vxp_", size="vxp_n")

        from externaltools import PileupReweighting
        from ROOT import Root
        # Initialize the pileup reweighting tool
        pileup_tool = Root.TPileupReweighting()
        if YEAR == 2011:
            pileup_tool.AddConfigFile(PileupReweighting.get_resource('mc11b_defaults.prw.root'))
            pileup_tool.AddLumiCalcFile('lumi/2011/hadhad/ilumicalc_histograms_None_178044-191933.root')
        elif YEAR == 2012:
            pileup_tool.AddConfigFile(PileupReweighting.get_resource('mc12a_defaults.prw.root'))
            pileup_tool.SetDataScaleFactors(1./1.11)
            pileup_tool.AddLumiCalcFile('lumi/2012/hadhad/ilumicalc_histograms_None_200841-205113.root')
        else:
            raise ValueError('No pileup reweighting defined for year %d' %
                    YEAR)
        # discard unrepresented data (with mu not simulated in MC)
        pileup_tool.SetUnrepresentedDataAction(2)
        pileup_tool.Initialize()

        # entering the main event loop...
        for event in chain:
            tree.reset()

            event.vertices.select(vertex_selection)
            tree.number_of_good_vertices = len(event.vertices)

            # match only with visible true taus
            event.truetaus.select(lambda tau: tau.vis_Et > 10 * GeV and abs(tau.vis_eta) < 2.5)

            if len(event.truetaus) == 1:
                true_tau = event.truetaus[0]
                TrueTauBlock.set(tree, 1, true_tau)
            else:
                continue

            # Truth-matching
            matched_reco = None
            reco_index = true_tau.tauAssoc_index
            tau = event.taus.getitem(reco_index)
            if tau in event.taus:
                matched_reco = tau
            else:
                continue

            tree.MET = event.MET_RefFinal_BDTMedium_et

            # fill tau block
            RecoTauBlock.set(event, tree, matched_reco, None)

            # set the event weight
            tree.pileup_weight = pileup_tool.GetCombinedWeight(event.RunNumber,
                                                               event.mc_channel_number,
                                                               event.averageIntPerXing)
            tree.mc_weight = event.mc_event_weight
            tree.Fill(reset=True)

        self.output.cd()
        tree.FlushBaskets()
        tree.Write()
        total_events = event_filters[0].total
        cutflow[0] = total_events
        cutflow[1] = total_events
        cutflow.Write()
