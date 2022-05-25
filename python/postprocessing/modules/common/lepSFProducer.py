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
            mu_f=["Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root",
                  "Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root",
                  "Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root"]
            mu_h = ["NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt",
                    "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                    "NUM_LooseRelIso_DEN_LooseID_abseta_pt"]
        elif muonSelectionTag=="LooseWP_2016_APV":
            mu_f=["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root",
                  "Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root",
                  "Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root"]
            mu_h = ["NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt",
                    "Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root",
                    "Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root"]
        elif muonSelectionTag=="LooseWP_2017":
            mu_f=["Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root",
                  "Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root",
                  "Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root"]
            mu_h = ["NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt",
                    "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                    "NUM_LooseRelIso_DEN_LooseID_abseta_pt"]
        elif muonSelectionTag=="LooseWP_2018":
            mu_f=["Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root",
                  "Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root",
                  "Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root"]
            mu_h = ["NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt",
                    "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                    "NUM_LooseRelIso_DEN_LooseID_abseta_pt"]
        # elif muonSelectionTag=="LooseWP_2018_runA":
        #     self.runA = True
        #     mu_f=["Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root",
        #           "Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root",
        #           "Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root"]
        #     mu_h = ["NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt",
        #             "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
        #             "NUM_LooseRelIso_DEN_LooseID_abseta_pt"]

        if electronSelectionTag=="GPMVA90_2016":
            el_f = ["2016_El_EGM2D_eleGSF.root",
                    "2016_El_EGM2D_eleMVA90.root"]
            el_h = ["EGamma_SF2D", "EGamma_SF2D"]
        elif electronSelectionTag=="GPMVA90_2017":
            el_f = ["egammaEffi.txt_EGM2D_MVA90iso_UL17.root",
                    #"egammaEffi_ptAbove20.txt_EGM2D_UL2017.root",
                    #"egammaEffi_ptBelow20.txt_EGM2D_UL2017.root"
                    ]
            el_h = ["EGamma_SF2D", 
            #"EGamma_SF2D", 
            #"EGamma_SF2D"
            ]
        elif electronSelectionTag=="GPMVA90_2018":
            el_f = ["egammaEffi.txt_Ele_wp90iso_EGM2D.root",
                    #"2018_El_egammaEffi.txt_EGM2D_updatedAll.root"
                    ]
            el_h = ["EGamma_SF2D", 
            #"EGamma_SF2D"
            ]

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
	self._worker_el = ROOT.LeptonEfficiencyCorrectorCppWorker(self.el_f,self.el_h)

    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

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

#	if self.runA:
#	self.out.fillBranch("Muon_SF_runA"       , sf_mu)
#	self.out.fillBranch("Muon_SFErr_runA"    , sferr_mu)
	#else:
        electrons = Collection(event, "Electron")
        sf_el = [ self._worker_el.getSF(el.pdgId,el.pt,el.eta) for el in electrons ]
        sferr_el = [ self._worker_el.getSFErr(el.pdgId,el.pt,el.eta) for el in electrons ]
	self.out.fillBranch("Muon_SF"       , sf_mu)
	self.out.fillBranch("Muon_SFErr"    , sferr_mu)
	self.out.fillBranch("Electron_SF"   , sf_el)
	self.out.fillBranch("Electron_SFErr", sferr_el)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

lepSF_2016_APV = lambda : lepSFProducer( "LooseWP_2016_APV", "GPMVA90_2016_APV")
lepSF_2016 = lambda : lepSFProducer( "LooseWP_2016", "GPMVA90_2016")
lepSF_2017 = lambda : lepSFProducer( "LooseWP_2017", "GPMVA90_2017")
#lepSF_2018_runA = lambda : lepSFProducer( "LooseWP_2018_runA", "GPMVA90_2018")
lepSF_2018 = lambda : lepSFProducer( "LooseWP_2018", "GPMVA90_2018")
