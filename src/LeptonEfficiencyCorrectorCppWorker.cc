#include "PhysicsTools/NanoAODTools/interface/LeptonEfficiencyCorrectorCppWorker.h"

#include <iostream>
#include <string>
#include <vector>
#include <TH2.h>
#include <TFile.h>

LeptonEfficiencyCorrectorCppWorker:: LeptonEfficiencyCorrectorCppWorker(std::vector<std::string> files, std::vector<std::string> histos) {
  effmaps_.clear();
  if(files.size()!=histos.size()) {
    std::cout << "ERROR! There should be one histogram per input file! Returning 0 as SF." << std::endl;
    return;
  }

  for(int i=0; i<(int)files.size();++i) {
    std::cout << " --> opening this file : " << files[i].c_str() << std::endl;
    TFile *f = TFile::Open(files[i].c_str(),"read");
    if(!f) {
      std::cout << "WARNING! File " << files[i] << " cannot be opened. Skipping this scale factor " << std::endl;
      continue;
    }
    TH2F *hist = (TH2F*)(f->Get(histos[i].c_str()))->Clone(("eff_"+histos[i]).c_str());
    hist->SetDirectory(0);
    if(!hist) {
      std::cout << "ERROR! Histogram " << histos[i] << " not in file " << files[i] << ". Not considering this SF. " << std::endl;
      continue;
    } else {
      std::cout << "Loading histogram " << histos[i] << " from file " << files[i] << "... " << std::endl;
    }
    std::cout << "---- doing a push_back" << std::endl;
    effmaps_.push_back(hist);
    std::cout << "---- closing the file" << std::endl;
    f->Close();
    std::cout << "---- file closed" << std::endl;
  }
  std::cout << "---- the end of the init function " << std::endl;
}

void LeptonEfficiencyCorrectorCppWorker::setLeptons(int nLep, int *lepPdgId, float *lepPt, float *lepEta) {
  nLep_ = nLep; Lep_pdgId_ = lepPdgId; Lep_pt_ = lepPt; Lep_eta_ = lepEta;
}

float LeptonEfficiencyCorrectorCppWorker::getSF(int pdgid, float pt, float eta) {
  float out=1.;
  float x = abs(pdgid)==13 ? pt : eta;
  float y = abs(pdgid)==13 ? fabs(eta) : pt;
  for(std::vector<TH2F*>::iterator hist=effmaps_.begin(); hist<effmaps_.end(); ++hist) {
    WeightCalculatorFromHistogram wc(*hist);
    out *= wc.getWeight(x,y);
  }
  return out;
}

float LeptonEfficiencyCorrectorCppWorker::getSFErr(int pdgid, float pt, float eta) {
  float out=1.;
  float x = pt;
  float y = abs(pdgid)==13 ? fabs(eta) : eta;
  for(std::vector<TH2F*>::iterator hist=effmaps_.begin(); hist<effmaps_.end(); ++hist) {
    WeightCalculatorFromHistogram wc(*hist);
    out *= wc.getWeightErr(x,y);
  }
  return out;
}

const std::vector<float> & LeptonEfficiencyCorrectorCppWorker::run() {
  ret_.clear();
  for (int iL = 0, nL = nLep_; iL < nL; ++iL) {
    ret_.push_back(getSF((Lep_pdgId_)[iL], (Lep_pt_)[iL], (Lep_eta_)[iL]));
  }
  return ret_;
}

