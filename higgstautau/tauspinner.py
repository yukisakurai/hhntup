"""
See instructions here:
https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG4/software/common/TauSpinner/trunk/README.rst

Based on original implementation by Daniele:
https://twiki.cern.ch/twiki/bin/view/Main/DanieleZanziTauSpinner
"""
from rootpy.tree.filtering import EventFilter
from . import log; log = log[__name__]

import ROOT

class EmbeddingTauSpinner(EventFilter):
    """
    This class applies Tau Spinner
    """
    def __init__(self, year, tree, passthrough=False, **kwargs):

        super(EmbeddingTauSpinner, self).__init__(
            passthrough=passthrough,
            **kwargs)

        if not passthrough:
            import TauSpinnerTool
            self.tree = tree
            self.d3pd_mc = ROOT.TauSpinnerHelpers.D3PD_MC()
            if year == 2012:
                self.tool = ROOT.TauSpinnerHelpers.HSG4.build8TeV()
            elif year == 2011:
                self.tool = ROOT.TauSpinnerHelpers.HSG4.build7TeV()
            else:
                raise ValueError(
                    "No TauSpinner config for year {0}".format(year))

    def passes(self, event):
        # If passthrough is True (not embedding) then this
        # method is never called
        self.d3pd_mc.set_addresses(
           event.mc_pt,
           event.mc_eta,
           event.mc_phi,
           event.mc_m,
           event.mc_pdgId,
           event.mc_status,
           event.mc_child_index)
        status = self.tool.read_event(self.d3pd_mc)
        self.tree.embedding_spin_weight = self.tool.get_spin_weight()

        return True

class TauSpinner(EventFilter):
    """
    This class applies Tau Spinner
    """
    def __init__(self, year, tree, passthrough=False, **kwargs):

        super(TauSpinner, self).__init__(
            passthrough=passthrough,
            **kwargs)

        if not passthrough:
            import TauSpinnerTool
            self.tree = tree
            self.d3pd_mc = ROOT.TauSpinnerHelpers.D3PD_MC()
            if year == 2012:
                self.tool = ROOT.TauSpinnerHelpers.HSG4.build8TeV()
            elif year == 2011:
                self.tool = ROOT.TauSpinnerHelpers.HSG4.build7TeV()
            else:
                raise ValueError(
                    "No TauSpinner config for year {0}".format(year))
            from externaltools import HtautauCPAna
            self.cptool = ROOT.HtautauCPAna()

    def passes(self, event):
        # If passthrough is True (not signal) then this
        # method is never called
        self.tree.true_pi0_minus_pt.clear()
        self.tree.true_pi0_minus_eta.clear()
        self.tree.true_pi0_minus_phi.clear()
        self.tree.true_pi0_minus_m.clear()
        self.tree.true_prong_minus_pt.clear()
        self.tree.true_prong_minus_eta.clear()
        self.tree.true_prong_minus_phi.clear()
        self.tree.true_prong_minus_m.clear()
        self.tree.true_pi0_plus_pt.clear()
        self.tree.true_pi0_plus_eta.clear()
        self.tree.true_pi0_plus_phi.clear()
        self.tree.true_pi0_plus_m.clear()
        self.tree.true_prong_plus_pt.clear()
        self.tree.true_prong_plus_eta.clear()
        self.tree.true_prong_plus_phi.clear()
        self.tree.true_prong_plus_m.clear()

        self.d3pd_mc.set_addresses(
           event.mc_pt,
           event.mc_eta,
           event.mc_phi,
           event.mc_m,
           event.mc_pdgId,
           event.mc_status,
           event.mc_child_index)

        status = self.tool.read_event(self.d3pd_mc)
        self.tree.embedding_spin_weight = self.tool.get_spin_weight()

        status = self.tool.read_event(self.d3pd_mc)
        self.tree.cp_odd_weight = self.tool.get_pseudoscalar_weight()


        status = self.tool.read_event(self.d3pd_mc)
        self.tree.cp_even_weight = self.tool.get_scalar_weight()

        pcles = self.tool.get_decay_particles()

        if pcles.X_idx==-1 or pcles.tau_idx==-1 or pcles.tau2_idx==-1:
            return True

        if event.mc_child_index[pcles.tau_idx].size()==0 or event.mc_child_index[pcles.tau2_idx].size()==0:
            return True

        self.tree.true_tau_minus_decay = ROOT.TauSpinnerHelpers.get_decay_mode(pcles.tau,pcles.tau_daughters);
        self.tree.true_tau_plus_decay = ROOT.TauSpinnerHelpers.get_decay_mode(pcles.tau2,pcles.tau2_daughters);
        dau_vistau_minus_size = pcles.tau_daughters.size();
        dau_vistau_plus_size = pcles.tau2_daughters.size();

        tlv_vis_minus = ROOT.TLorentzVector();    tlv_vis_plus = ROOT.TLorentzVector()
        tlv_nu_minus = ROOT.TLorentzVector();     tlv_nu_plus = ROOT.TLorentzVector()
        tlv_vis_minus_prong = ROOT.TLorentzVector();  tlv_vis_plus_prong = ROOT.TLorentzVector()
        tlv_vis_minus_pi0 = ROOT.TLorentzVector();    tlv_vis_plus_pi0 = ROOT.TLorentzVector()

        PV_xyz =  ROOT.TVector3( event.mc_vx_x[pcles.tau_idx],
                                 event.mc_vx_y[pcles.tau_idx],
                                 event.mc_vx_z[pcles.tau_idx] )
        vis_minus_idx = event.mc_child_index[pcles.tau_idx][0]
        vis_plus_idx = event.mc_child_index[pcles.tau2_idx][0]
        vis_minus_xyz = ROOT.TVector3( event.mc_vx_x[vis_minus_idx],
                                       event.mc_vx_y[vis_minus_idx],
                                       event.mc_vx_z[vis_minus_idx] )
        vis_plus_xyz = ROOT.TVector3( event.mc_vx_x[vis_plus_idx],
                                      event.mc_vx_y[vis_plus_idx],
                                      event.mc_vx_z[vis_plus_idx] )

        for i in range(0,pcles.tau_daughters.size()):
            daughters = ROOT.TauSpinnerHelpers.simp2tlv(pcles.tau_daughters[i])
            pdgId = pcles.tau_daughters[i].pdgid()

            if ROOT.TMath.Abs(pdgId)==16:
                tlv_nu_minus += daughters
            elif ROOT.TMath.Abs(pdgId)==111 or ROOT.TMath.Abs(pdgId)==310 or ROOT.TMath.Abs(pdgId)==130:
                tlv_vis_minus += daughters
                tlv_vis_minus_pi0 += daughters
                self.tree.true_pi0_minus_pt.push_back(daughters.Pt())
                self.tree.true_pi0_minus_eta.push_back(daughters.Eta())
                self.tree.true_pi0_minus_phi.push_back(daughters.Phi())
                self.tree.true_pi0_minus_m.push_back(daughters.M())
            else:
                tlv_vis_minus += daughters
                tlv_vis_minus_prong += daughters
                self.tree.true_prong_minus_pt.push_back(daughters.Pt())
                self.tree.true_prong_minus_eta.push_back(daughters.Eta())
                self.tree.true_prong_minus_phi.push_back(daughters.Phi())
                self.tree.true_prong_minus_m.push_back(daughters.M())

        for i in range(0,pcles.tau2_daughters.size()):
            daughters = ROOT.TauSpinnerHelpers.simp2tlv(pcles.tau2_daughters[i])
            pdgId = pcles.tau2_daughters[i].pdgid()

            if ROOT.TMath.Abs(pdgId)==16:
                tlv_nu_plus += daughters
            elif ROOT.TMath.Abs(pdgId)==111 or ROOT.TMath.Abs(pdgId)==310 or ROOT.TMath.Abs(pdgId)==130:
                tlv_vis_plus += daughters
                tlv_vis_plus_pi0 += daughters
                self.tree.true_pi0_plus_pt.push_back(daughters.Pt())
                self.tree.true_pi0_plus_eta.push_back(daughters.Eta())
                self.tree.true_pi0_plus_phi.push_back(daughters.Phi())
                self.tree.true_pi0_plus_m.push_back(daughters.M())
            else:
                tlv_vis_plus += daughters
                tlv_vis_plus_prong += daughters
                self.tree.true_prong_plus_pt.push_back(daughters.Pt())
                self.tree.true_prong_plus_eta.push_back(daughters.Eta())
                self.tree.true_prong_plus_phi.push_back(daughters.Phi())
                self.tree.true_prong_plus_m.push_back(daughters.M())

        self.tree.true_tau_minus_pt  = event.mc_pt[pcles.tau_idx]
        self.tree.true_tau_minus_eta = event.mc_eta[pcles.tau_idx]
        self.tree.true_tau_minus_phi = event.mc_phi[pcles.tau_idx]
        self.tree.true_tau_minus_m   = event.mc_m[pcles.tau_idx]
        self.tree.true_tau_plus_pt  = event.mc_pt[pcles.tau2_idx]
        self.tree.true_tau_plus_eta = event.mc_eta[pcles.tau2_idx]
        self.tree.true_tau_plus_phi = event.mc_phi[pcles.tau2_idx]
        self.tree.true_tau_plus_m   = event.mc_m[pcles.tau2_idx]

        self.tree.true_vistau_minus_pt  = tlv_vis_minus.Pt()
        self.tree.true_vistau_minus_eta = tlv_vis_minus.Eta()
        self.tree.true_vistau_minus_phi = tlv_vis_minus.Phi()
        self.tree.true_vistau_minus_m   = tlv_vis_minus.M()
        self.tree.true_vistau_plus_pt  = tlv_vis_plus.Pt()
        self.tree.true_vistau_plus_eta = tlv_vis_plus.Eta()
        self.tree.true_vistau_plus_phi = tlv_vis_plus.Phi()
        self.tree.true_vistau_plus_m   = tlv_vis_plus.M()

        self.tree.true_nu_minus_pt  = tlv_nu_minus.Pt()
        self.tree.true_nu_minus_eta = tlv_nu_minus.Eta()
        self.tree.true_nu_minus_phi = tlv_nu_minus.Phi()
        self.tree.true_nu_minus_m   = tlv_nu_minus.M()
        self.tree.true_nu_plus_pt  = tlv_nu_plus.Pt()
        self.tree.true_nu_plus_eta = tlv_nu_plus.Eta()
        self.tree.true_nu_plus_phi = tlv_nu_plus.Phi()
        self.tree.true_nu_plus_m   = tlv_nu_plus.M()

        self.tree.true_Acoplanarity = self.cptool.Acoplanarity(tlv_vis_minus,tlv_nu_minus,tlv_vis_plus,tlv_nu_plus)
        self.tree.true_Acollinearity = self.cptool.Acollinearity(tlv_vis_minus,tlv_nu_minus,tlv_vis_plus,tlv_nu_plus)

        # calculate true pca
        vis_minus_pca_xyz = self.cptool.linear_pca_xyz3D(tlv_vis_minus,PV_xyz,vis_minus_xyz)
        vis_plus_pca_xyz = self.cptool.linear_pca_xyz3D(tlv_vis_plus,PV_xyz,vis_plus_xyz)

        # for pi-pi decay mode
        self.tree.true_Acoplanarity_IP = \
            self.cptool.Acoplanarity_IP(tlv_vis_minus,tlv_vis_plus,PV_xyz,vis_minus_pca_xyz,vis_plus_pca_xyz)
        self.tree.true_Acoplanarity_IP_lab = \
            self.cptool.Acoplanarity_IP(tlv_vis_minus,tlv_vis_plus,PV_xyz,vis_minus_pca_xyz,vis_plus_pca_xyz)
        self.tree.true_Acoplanarity_IP_psi = \
            self.cptool.Acoplanarity_IP(tlv_vis_minus,tlv_vis_plus,PV_xyz,vis_minus_pca_xyz,vis_plus_pca_xyz)

        # for pi-rho decay mode
        if not tlv_vis_minus_pi0.Pt()==0 and not tlv_vis_plus_pi0.Pt()==0:

            if not tlv_vis_minus_pi0.Pt()==0 and tlv_vis_plus_pi0.Pt()==0:
                self.tree.true_Acoplanarity_IP_rho = \
                    self.cptool.Acoplanarity_IP_rho( \
                    tlv_vis_plus,PV_xyz,vis_plus_pca_xyz,tlv_vis_minus_prong,tlv_vis_minus_pi0,True)

            else:
                self.tree.true_Acoplanarity_IP_rho = \
                    self.cptool.Acoplanarity_IP_rho( \
                    tlv_vis_minus,PV_xyz,vis_minus_pca_xyz,tlv_vis_plus_prong,tlv_vis_plus_pi0,True)

        # for rho-rho decay mode
        if not tlv_vis_minus_pi0.Pt()==0 and not tlv_vis_plus_pi0.Pt()==0:
            self.tree.true_Acoplanarity_rho = \
                self.cptool.Acoplanarity_rho( \
                tlv_vis_minus_prong,tlv_vis_minus_pi0,tlv_vis_plus_prong,tlv_vis_plus_pi0,True)

        self.tree.true_Acoplanarity_visboost = self.cptool.Acoplanarity_visboost(tlv_vis_minus,tlv_nu_minus,tlv_vis_plus,tlv_nu_plus)

        return True

