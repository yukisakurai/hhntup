from rootpy import log as rlog
from . import log; log = log[__name__]
import ROOT

from units import GeV
import datasets

import os
from math import sqrt
from externaltools import HtautauCPAna


class HttCP(object):

    INITED = False

    def __init__(self):

        if HttCP.INITED:
            raise RuntimeError("do not create more than one HtautauCPAna")
        HttCP.INITED = True

        self.tool = ROOT.HtautauCPAna()
        self.pi_mass = 139.57018
        self.pi0_mass = 134.9766
        self.rho_mass = 775.4


    def find_privtx(self, vertices):

        privtx_xyz = ROOT.TVector3()
        for vxp in vertices:
            if vxp.type == 1 and vxp.nTracks >= 4:
                privtx_xyz.SetXYZ(vxp.x,vxp.y,vxp.z)
                break

        return privtx_xyz

    def find_lead_track(self, tau):

        # find lead track
        tau_lead_track_tlv = ROOT.TLorentzVector()
        tau_lead_track_d0 = 0.0
        tau_lead_track_z0 = 0.0
        tau_lead_track_qoverp = 0.0
        for i in range(tau.track_pt.size()):
            if tau_lead_track_tlv.Pt() < tau.track_pt[i]:
                tau_lead_track_tlv.SetPtEtaPhiM(
                    tau.track_pt[i],
                    tau.track_eta[i],
                    tau.track_phi[i],
                    self.pi_mass)
                tau_lead_track_d0 = tau.track_d0[i]
                tau_lead_track_z0 = tau.track_z0[i]
                tau_lead_track_qoverp = tau.track_qoverp[i]

        return tau_lead_track_tlv, tau_lead_track_qoverp, tau_lead_track_d0, tau_lead_track_z0

    def Acoplanarity_IP(self, tau1, tau2, vertices):

        privtx_xyz = self.find_privtx( vertices )

        tau1_lead_track_tlv = ROOT.TLorentzVector()
        tau1_lead_track_qoverp = 0.0
        tau1_lead_track_d0 = 0.0
        tau1_lead_track_z0 = 0.0

        tau1_lead_track_tlv, \
            tau1_lead_track_qoverp, \
            tau1_lead_track_d0, \
            tau1_lead_track_z0 = self.find_lead_track ( tau1 )


        tau2_lead_track_tlv = ROOT.TLorentzVector()
        tau2_lead_track_qoverp = 0.0
        tau2_lead_track_d0 = 0.0
        tau2_lead_track_z0 = 0.0

        tau2_lead_track_tlv, \
            tau2_lead_track_qoverp, \
            tau2_lead_track_d0, \
            tau2_lead_track_z0 = self.find_lead_track ( tau2 )

        tau1_lead_track_pca_xyz = self.tool.pca_xyz3D (
            privtx_xyz,
            tau1_lead_track_qoverp,tau1_lead_track_d0,tau1_lead_track_z0,
            tau1_lead_track_tlv.Phi(),tau1_lead_track_tlv.Theta() )

        tau2_lead_track_pca_xyz = self.tool.pca_xyz3D (
            privtx_xyz,
            tau2_lead_track_qoverp,tau2_lead_track_d0,tau2_lead_track_z0,
            tau2_lead_track_tlv.Phi(),tau2_lead_track_tlv.Theta() )

        return self.tool.Acoplanarity_IP(
            tau1_lead_track_tlv, tau2_lead_track_tlv,
            privtx_xyz, tau1_lead_track_pca_xyz,tau2_lead_track_pca_xyz )

    def Acoplanarity_rho_leadtrk(self, tau1, tau2):

        tau1_tlv = ROOT.TLorentzVector()
        tau1_tlv.SetPtEtaPhiM(
            tau1.fourvect.Pt(),
            tau1.fourvect.Eta(),
            tau1.fourvect.Phi(),
            self.rho_mass)

        tau2_tlv= ROOT.TLorentzVector()
        tau2_tlv.SetPtEtaPhiM(
            tau2.fourvect.Pt(),
            tau2.fourvect.Eta(),
            tau2.fourvect.Phi(),
            self.rho_mass)

        tau1_lead_track_tlv = ROOT.TLorentzVector()
        tau1_lead_track_tlv, \
            tau1_lead_track_qoverp, \
            tau1_lead_track_d0, \
            tau1_lead_track_z0 = self.find_lead_track ( tau1 )

        tau2_lead_track_tlv = ROOT.TLorentzVector()
        tau2_lead_track_tlv, \
            tau2_lead_track_qoverp, \
            tau2_lead_track_d0, \
            tau2_lead_track_z0 = self.find_lead_track ( tau2 )

        tau1_pi0_tlv = tau1_tlv - tau1_lead_track_tlv;
        tau2_pi0_tlv = tau2_tlv - tau2_lead_track_tlv;

        return self.tool.Acoplanarity_rho(
            tau1_lead_track_tlv, tau1_pi0_tlv,
            tau2_lead_track_tlv, tau2_pi0_tlv)

    def Acoplanarity_tau1_IP_tau2_rho_leadtrk(self, tau1, tau2, vertices):

        privtx_xyz = self.find_privtx( vertices )

        tau1_lead_track_tlv = ROOT.TLorentzVector()
        tau1_lead_track_tlv, \
            tau1_lead_track_qoverp, \
            tau1_lead_track_d0, \
            tau1_lead_track_z0 = self.find_lead_track ( tau1 )

        tau1_lead_track_pca_xyz = self.tool.pca_xyz3D (
            privtx_xyz,
            tau1_lead_track_qoverp,tau1_lead_track_d0,tau1_lead_track_z0,
            tau1_lead_track_tlv.Phi(),tau1_lead_track_tlv.Theta() )

        tau2_tlv= ROOT.TLorentzVector()
        tau2_tlv.SetPtEtaPhiM(
            tau2.fourvect.Pt(),
            tau2.fourvect.Eta(),
            tau2.fourvect.Phi(),
            self.rho_mass)

        tau2_lead_track_tlv = ROOT.TLorentzVector()
        tau2_lead_track_tlv, \
            tau2_lead_track_qoverp, \
            tau2_lead_track_d0, \
            tau2_lead_track_z0 = self.find_lead_track ( tau2 )

        tau2_pi0_tlv = tau2_tlv - tau2_lead_track_tlv;

        return self.tool.Acoplanarity_IP_rho(
            tau1_lead_track_tlv, privtx_xyz, tau1_lead_track_pca_xyz,
            tau2_lead_track_tlv, tau2_pi0_tlv)

    def Acoplanarity_tau2_IP_tau1_rho_leadtrk(self, tau1, tau2, vertices):

        privtx_xyz = self.find_privtx( vertices )

        tau2_lead_track_tlv = ROOT.TLorentzVector()
        tau2_lead_track_tlv, \
            tau2_lead_track_qoverp, \
            tau2_lead_track_d0, \
            tau2_lead_track_z0 = self.find_lead_track ( tau2 )

        tau2_lead_track_pca_xyz = self.tool.pca_xyz3D (
            privtx_xyz,
            tau2_lead_track_qoverp,tau2_lead_track_d0,tau2_lead_track_z0,
            tau2_lead_track_tlv.Phi(),tau2_lead_track_tlv.Theta() )

        tau1_tlv= ROOT.TLorentzVector()
        tau1_tlv.SetPtEtaPhiM(
            tau1.fourvect.Pt(),
            tau1.fourvect.Eta(),
            tau1.fourvect.Phi(),
            self.rho_mass)

        tau1_lead_track_tlv = ROOT.TLorentzVector()
        tau1_lead_track_tlv, \
            tau1_lead_track_qoverp, \
            tau1_lead_track_d0, \
            tau1_lead_track_z0 = self.find_lead_track ( tau1 )

        tau1_pi0_tlv = tau1_tlv - tau1_lead_track_tlv;

        return self.tool.Acoplanarity_IP_rho(
            tau2_lead_track_tlv, privtx_xyz, tau2_lead_track_pca_xyz,
            tau1_lead_track_tlv, tau1_pi0_tlv)

    def Acoplanarity_rho_cluster(self, tau1, tau2):


        tau1_vis_tlv= ROOT.TLorentzVector()
        tau1_vis_tlv.SetPtEtaPhiM(
            tau1.pi0_vistau_pt,
            tau1.pi0_vistau_eta,
            tau1.pi0_vistau_phi,
            self.pi_mass)

        tau1_cl1_tlv= ROOT.TLorentzVector()
        tau1_cl1_tlv.SetPtEtaPhiM(
            tau1.pi0_cl1_pt,
            tau1.pi0_cl1_eta,
            tau1.pi0_cl1_phi,
            self.pi0_mass)

        tau1_cl2_tlv= ROOT.TLorentzVector()
        tau1_cl2_tlv.SetPtEtaPhiM(
            tau1.pi0_cl2_pt,
            tau1.pi0_cl2_eta,
            tau1.pi0_cl2_phi,
            self.pi0_mass)

        tau1_pi0_tlv = ROOT.TLorentzVector()
        tau1_pi0_tlv = tau1_cl1_tlv + tau1_cl2_tlv
        tau1_prong_tlv = ROOT.TLorentzVector()
        tau1_prong_tlv = tau1_vis_tlv - tau1_pi0_tlv

        tau2_vis_tlv= ROOT.TLorentzVector()
        tau2_vis_tlv.SetPtEtaPhiM(
            tau2.pi0_vistau_pt,
            tau2.pi0_vistau_eta,
            tau2.pi0_vistau_phi,
            self.pi_mass)

        tau2_cl1_tlv= ROOT.TLorentzVector()
        tau2_cl1_tlv.SetPtEtaPhiM(
            tau2.pi0_cl1_pt,
            tau2.pi0_cl1_eta,
            tau2.pi0_cl1_phi,
            self.pi0_mass)

        tau2_cl2_tlv= ROOT.TLorentzVector()
        tau2_cl2_tlv.SetPtEtaPhiM(
            tau2.pi0_cl2_pt,
            tau2.pi0_cl2_eta,
            tau2.pi0_cl2_phi,
            self.pi0_mass)

        tau2_pi0_tlv = ROOT.TLorentzVector()
        tau2_pi0_tlv = tau2_cl1_tlv + tau2_cl2_tlv
        tau2_prong_tlv = ROOT.TLorentzVector()
        tau2_prong_tlv = tau2_vis_tlv - tau2_pi0_tlv

        return self.tool.Acoplanarity_rho(
            tau1_prong_tlv, tau1_pi0_tlv,
            tau2_prong_tlv, tau2_pi0_tlv)

    def Acoplanarity_tau1_IP_tau2_rho_cluster(self, tau1, tau2, vertices):

        privtx_xyz = self.find_privtx( vertices )

        tau1_lead_track_tlv = ROOT.TLorentzVector()
        tau1_lead_track_tlv, \
            tau1_lead_track_qoverp, \
            tau1_lead_track_d0, \
            tau1_lead_track_z0 = self.find_lead_track ( tau1 )

        tau1_lead_track_pca_xyz = self.tool.pca_xyz3D (
            privtx_xyz,
            tau1_lead_track_qoverp,tau1_lead_track_d0,tau1_lead_track_z0,
            tau1_lead_track_tlv.Phi(),tau1_lead_track_tlv.Theta() )

        tau2_vis_tlv= ROOT.TLorentzVector()
        tau2_vis_tlv.SetPtEtaPhiM(
            tau2.pi0_vistau_pt,
            tau2.pi0_vistau_eta,
            tau2.pi0_vistau_phi,
            self.pi_mass)

        tau2_cl1_tlv= ROOT.TLorentzVector()
        tau2_cl1_tlv.SetPtEtaPhiM(
            tau2.pi0_cl1_pt,
            tau2.pi0_cl1_eta,
            tau2.pi0_cl1_phi,
            self.pi0_mass)

        tau2_cl2_tlv= ROOT.TLorentzVector()
        tau2_cl2_tlv.SetPtEtaPhiM(
            tau2.pi0_cl2_pt,
            tau2.pi0_cl2_eta,
            tau2.pi0_cl2_phi,
            self.pi0_mass)

        tau2_pi0_tlv = ROOT.TLorentzVector()
        tau2_pi0_tlv = tau2_cl1_tlv + tau2_cl2_tlv
        tau2_prong_tlv = ROOT.TLorentzVector()
        tau2_prong_tlv = tau2_vis_tlv - tau2_pi0_tlv

        return self.tool.Acoplanarity_IP_rho(
            tau1_lead_track_tlv, privtx_xyz, tau1_lead_track_pca_xyz,
            tau2_prong_tlv, tau2_pi0_tlv)

    def Acoplanarity_tau2_IP_tau1_rho_cluster(self, tau1, tau2, vertices):

        privtx_xyz = self.find_privtx( vertices )

        tau2_lead_track_tlv = ROOT.TLorentzVector()
        tau2_lead_track_tlv, \
            tau2_lead_track_qoverp, \
            tau2_lead_track_d0, \
            tau2_lead_track_z0 = self.find_lead_track ( tau2 )

        tau2_lead_track_pca_xyz = self.tool.pca_xyz3D (
            privtx_xyz,
            tau2_lead_track_qoverp,tau2_lead_track_d0,tau2_lead_track_z0,
            tau2_lead_track_tlv.Phi(),tau2_lead_track_tlv.Theta() )

        tau1_vis_tlv= ROOT.TLorentzVector()
        tau1_vis_tlv.SetPtEtaPhiM(
            tau1.pi0_vistau_pt,
            tau1.pi0_vistau_eta,
            tau1.pi0_vistau_phi,
            self.pi_mass)

        tau1_cl1_tlv= ROOT.TLorentzVector()
        tau1_cl1_tlv.SetPtEtaPhiM(
            tau1.pi0_cl1_pt,
            tau1.pi0_cl1_eta,
            tau1.pi0_cl1_phi,
            self.pi0_mass)

        tau1_cl2_tlv= ROOT.TLorentzVector()
        tau1_cl2_tlv.SetPtEtaPhiM(
            tau1.pi0_cl2_pt,
            tau1.pi0_cl2_eta,
            tau1.pi0_cl2_phi,
            self.pi0_mass)

        tau1_pi0_tlv = ROOT.TLorentzVector()
        tau1_pi0_tlv = tau1_cl1_tlv + tau1_cl2_tlv
        tau1_prong_tlv = ROOT.TLorentzVector()
        tau1_prong_tlv = tau1_vis_tlv - tau1_pi0_tlv

        return self.tool.Acoplanarity_IP_rho(
            tau2_lead_track_tlv, privtx_xyz, tau2_lead_track_pca_xyz,
            tau1_prong_tlv, tau1_pi0_tlv)

    def Acoplanarity_mmc(self, vistau1, nu1, vistau2, nu2):

        vistau1_tlv = ROOT.TLorentzVector();
        vistau1_tlv.SetPtEtaPhiM(vistau1.pt,vistau1.eta,vistau1.phi,vistau1.m)
        vistau2_tlv = ROOT.TLorentzVector();
        vistau2_tlv.SetPtEtaPhiM(vistau2.pt,vistau2.eta,vistau2.phi,vistau2.m)

        nu1_tlv = ROOT.TLorentzVector();
        nu1_tlv.SetPtEtaPhiM(nu1.pt,nu1.eta,nu1.phi,nu1.m)
        nu2_tlv = ROOT.TLorentzVector();
        nu2_tlv.SetPtEtaPhiM(nu2.pt,nu2.eta,nu2.phi,nu2.m)

        return self.tool.Acoplanarity(vistau1_tlv,nu1_tlv,vistau2_tlv,nu2_tlv)

    def Acoplanarity_mmc_visboost(self, vistau1, nu1, vistau2, nu2):

        vistau1_tlv = ROOT.TLorentzVector();
        vistau1_tlv.SetPtEtaPhiM(vistau1.pt,vistau1.eta,vistau1.phi,vistau1.m)
        vistau2_tlv = ROOT.TLorentzVector();
        vistau2_tlv.SetPtEtaPhiM(vistau2.pt,vistau2.eta,vistau2.phi,vistau2.m)

        nu1_tlv = ROOT.TLorentzVector();
        nu1_tlv.SetPtEtaPhiM(nu1.pt,nu1.eta,nu1.phi,nu1.m)
        nu2_tlv = ROOT.TLorentzVector();
        nu2_tlv.SetPtEtaPhiM(nu2.pt,nu2.eta,nu2.phi,nu2.m)

        return self.tool.Acoplanarity_visboost(vistau1_tlv,nu1_tlv,vistau2_tlv,nu2_tlv)
