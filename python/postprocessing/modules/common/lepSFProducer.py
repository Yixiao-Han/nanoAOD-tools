import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class lepSFProducer(Module):
    def __init__(self, muonSelectionTag, electronSelectionTag):
        self.runA = False

	if muonSelectionTag=="LooseWP_2016":
            mu_f=["2016_Mu_Trg.root",
                  "2016_Mu_ID.root",
                  "2016_Mu_Iso.root"]
            mu_h = ["IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio",
                    "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio",
                    "LooseISO_LooseID_pt_eta/pt_abseta_ratio"]
        elif muonSelectionTag=="LooseWP_2017":
            mu_f=["2017_Mu_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root",
                  "2017_Mu_RunBCDEF_mc_ID.root",
                  "2017_Mu_RunBCDEF_mc_ISO.root"]
            mu_h = ["IsoMu27_PtEtaBins/pt_abseta_ratio",
                    "NUM_LooseID_DEN_genTracks_pt_abseta",
                    "NUM_LooseRelIso_DEN_LooseID_pt_abseta"]
        elif muonSelectionTag=="LooseWP_2018":
            mu_f=["2018_Mu_EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root",
                  "2018_Mu_EfficienciesStudies_2018_rootfiles_RunABCD_SF_ID.root",
                  "2018_Mu_EfficienciesStudies_2018_rootfiles_RunABCD_SF_ISO.root"]
            mu_h = ["IsoMu24_PtEtaBins/pt_abseta_ratio",
                    "NUM_LooseID_DEN_TrackerMuons_pt_abseta",
                    "NUM_LooseRelIso_DEN_LooseID_pt_abseta"]
        elif muonSelectionTag=="LooseWP_2018_runA":
            self.runA = True
            mu_f=["2018_Mu_EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root",
                  "2018_Mu_EfficienciesStudies_2018_rootfiles_RunABCD_SF_ID.root",
                  "2018_Mu_EfficienciesStudies_2018_rootfiles_RunABCD_SF_ISO.root"]
            mu_h = ["IsoMu24_PtEtaBins/pt_abseta_ratio",
                    "NUM_LooseID_DEN_TrackerMuons_pt_abseta",
                    "NUM_LooseRelIso_DEN_LooseID_pt_abseta"]

        if electronSelectionTag=="GPMVA90_2016":
            el_f = ["2016_El_EGM2D_eleGSF.root",
                    "2016_El_EGM2D_eleMVA90.root"]
            el_h = ["EGamma_SF2D", "EGamma_SF2D"]
        elif electronSelectionTag=="GPMVA90_2017":
            el_f = ["2017_El_MVA90.root",
                    "2017_El_egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root",
                    "2017_El_egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root"]
            el_h = ["EGamma_SF2D", "EGamma_SF2D", "EGamma_SF2D"]
        elif electronSelectionTag=="GPMVA90_2018":
            el_f = ["2018_El_MVA90.root",
                    "2018_El_egammaEffi.txt_EGM2D_updatedAll.root"]
            el_h = ["EGamma_SF2D", "EGamma_SF2D"]

        mu_f = ["%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in mu_f]
        el_f = ["%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in el_f]

        self.mu_f = ROOT.std.vector(str)(len(mu_f))
        self.mu_h = ROOT.std.vector(str)(len(mu_f))
        for i in range(len(mu_f)): self.mu_f[i] = mu_f[i]; self.mu_h[i] = mu_h[i];
        self.el_f = ROOT.std.vector(str)(len(el_f))
        self.el_h = ROOT.std.vector(str)(len(el_f))
        for i in range(len(el_f)): self.el_f[i] = el_f[i]; self.el_h[i] = el_h[i];
        if "/LeptonEfficiencyCorrector_cc.so" not in ROOT.gSystem.GetLibraries():
            print "Load C++ Worker"
            ROOT.gROOT.ProcessLine(".L %s/src/PhysicsTools/NanoAODTools/python/postprocessing/helpers/LeptonEfficiencyCorrector.cc+" % os.environ['CMSSW_BASE'])
    def beginJob(self):
        self._worker_mu = ROOT.LeptonEfficiencyCorrectorCppWorker(self.mu_f,self.mu_h)
        if not self.runA:
	    self._worker_el = ROOT.LeptonEfficiencyCorrectorCppWorker(self.el_f,self.el_h)

    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

	if self.runA:
	    self.out.branch("Muon_SF_runA"       , "F", lenVar="nMuon"    )
	    self.out.branch("Muon_SFErr_runA"    , "F", lenVar="nMuon"    )
        else:
	    self.out.branch("Muon_SF"       , "F", lenVar="nMuon"    )
	    self.out.branch("Muon_SFErr"    , "F", lenVar="nMuon"    )
	    self.out.branch("Electron_SF"   , "F", lenVar="nElectron")
	    self.out.branch("Electron_SFErr", "F", lenVar="nElectron")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons     = Collection(event, "Muon")

        sf_mu = [ self._worker_mu.getSF(mu.pdgId,mu.pt,mu.eta) for mu in muons ]
        sferr_mu = [ self._worker_mu.getSFErr(mu.pdgId,mu.pt,mu.eta) for mu in muons ]

	if self.runA:
	    self.out.fillBranch("Muon_SF_runA"       , sf_mu)
	    self.out.fillBranch("Muon_SFErr_runA"    , sferr_mu)
	else:
            electrons = Collection(event, "Electron")
            sf_el = [ self._worker_el.getSF(el.pdgId,el.pt,el.eta) for el in electrons ]
            sferr_el = [ self._worker_el.getSFErr(el.pdgId,el.pt,el.eta) for el in electrons ]
	    self.out.fillBranch("Muon_SF"       , sf_mu)
	    self.out.fillBranch("Muon_SFErr"    , sferr_mu)
	    self.out.fillBranch("Electron_SF"   , sf_el)
	    self.out.fillBranch("Electron_SFErr", sferr_el)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

lepSF_2016 = lambda : lepSFProducer( "LooseWP_2016", "GPMVA90_2016")
lepSF_2017 = lambda : lepSFProducer( "LooseWP_2017", "GPMVA90_2017")
lepSF_2018_runA = lambda : lepSFProducer( "LooseWP_2018_runA", "GPMVA90_2018")
lepSF_2018 = lambda : lepSFProducer( "LooseWP_2018", "GPMVA90_2018")
