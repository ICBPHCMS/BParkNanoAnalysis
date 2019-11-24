/*
A macro to produce input variables for BDT
For data (--isMC 0), both KEE and KMuMu variables are stored in the output
For MC (--isMC 1), the output consists of either KEE (--isEE 1) or KMuMu (--isEE 0) variables
----------------
g++ -Wall -o RDF_MC_Kll_jobs `root-config --cflags --glibs ` RDF_MC_Kll_jobs.cpp
----------------
./RDF_MC_Kll_jobs --JOBid (-1,1,2...) --inList nanoAOD_list.txt --outFile /path/to/output_file.root --isMC (0,1) --isEE (0,1) --isResonant (0,1) --testFile /path/to/input_file.root
*/

#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TStopwatch.h"
#include <vector>
#include <string>
#include <regex>

#include "interface/functions.h"
//using namespace ROOT::VecOps;


int main(int argc, char **argv){

  TStopwatch t;
  t.Start();

  if(argc < 2) {
    std::cout << " Missing arguments " << std::endl;
    return -1;
  }

  std::string JOBid = "-1";
  std::string inList = "-1";
  std::string outFName = "-1";
  int isMC = 0;
  int isEE = 0;
  int isResonant = 0;
  std::string testFile = "-1";

  parseInputs(argc, argv, JOBid, inList, outFName, isMC, isEE, isResonant, testFile);

  if(inList != "-1" && JOBid == "-1"){
    std::cout << " running a test => careful: splitting file based but missing JOBid and output file path " << std::endl;
  }

  if(inList == "-1" && testFile == "-1"){
    std::cout << " by default " << std::endl;
  }

  std::cout << " inList = " << inList << " JOBid = " << JOBid << " outFName = " << outFName << " testFile = " << testFile
	    << " isMC = " << isMC << " isEE = " << isEE << " isResonant = " << isResonant << "\n" << std::endl;

  std::vector<std::string> inputFileList;
  
  if(inList != "-1"){
    std::string rootFileName;
    std::ifstream inFileLong;
    inFileLong.open(inList.c_str(), std::ios::in);
    while(!inFileLong.eof()){
      inFileLong >> rootFileName;
      if( inFileLong.eof() ) break;
      inputFileList.push_back(rootFileName);
    }
  }
  else if(testFile != "-1"){
    inputFileList.push_back(testFile);
  }
  else{
    if(isMC){
    std::string base = "/eos/cms//store/group/cmst3/group/bpark/BParkingNANO_2019Oct25/";
	if(isEE) 
inputFileList.push_back(base + "BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee/191025_125913/0000/BParkNANO_mc_2019Oct25_1.root");
	else 
inputFileList.push_back(base + "BuToKJpsi_ToMuMu_probefilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_ToMuMu/191025_125744/0000/BParkNANO_mc_2019Oct25_1.root");
    }
    else{
    std::string base = "/eos/cms//store/group/cmst3/group/bpark/BParkingNANO_2019Oct21/";
    inputFileList.push_back(base + "ParkingBPH2/crab_data_Run2018A_part2/191021_131326/0000/BParkNANO_data_2019Oct21_993.root");
    /*
    //A2
    inputFileList.push_back(base + "ParkingBPH2/crab_data_Run2018A_part2/191021_131326/0000/BParkNANO_data_2019Oct21_*.root");
    inputFileList.push_back(base + "ParkingBPH2/crab_data_Run2018A_part2/191021_131326/0001/BParkNANO_data_2019Oct21_*.root");
    //A3
    inputFileList.push_back(base + "ParkingBPH3/crab_data_Run2018A_part3/191021_131508/0000/BParkNANO_data_2019Oct21_*.root");
    inputFileList.push_back(base + "ParkingBPH3/crab_data_Run2018A_part3/191021_131508/0001/BParkNANO_data_2019Oct21_*.root");
    */
    }
  }//defaut
 
  inList = "-1";
  for(auto ij : inputFileList) std::cout << " file = " << ij << std::endl;
 
  
  ROOT::RDataFrame d("Events", inputFileList);
  
  auto n0 = d.Define("lumi_All", "luminosityBlock")
             .Define("eventToT_All", "event")
             .Define("run_All", "run");

  auto n = n0;
             
  if(!isMC){
    n = n0.Define("nBToKEE_All", "nBToKEE")
     .Define("KEE_fit_mass_All", "BToKEE_fit_mass")
     .Define("KEE_fit_massErr_All", "BToKEE_fit_massErr")
     .Define("KEE_cos2D_All", "BToKEE_cos2D")
     .Define("KEE_vtxProb_All", "BToKEE_svprob")
     .Define("KEE_l_xy_unc_All", "BToKEE_l_xy_unc")
     .Define("KEE_l_xyS_All", "BToKEE_l_xy/BToKEE_l_xy_unc")
     .Define("KEE_pT_All", "BToKEE_fit_pt")
     .Define("KEE_eta_All", "BToKEE_fit_eta")
     .Define("KEE_phi_All", "BToKEE_fit_phi")
     .Define("KEE_mll_fullfit_All", "BToKEE_mll_fullfit")
     .Define("KEE_k_eta_All", "BToKEE_fit_k_eta")
     .Define("KEE_k_phi_All", "BToKEE_fit_k_phi")
     .Define("KEE_k_pt_All", "BToKEE_fit_k_pt")
     .Define("KEE_k_dz_All", "Take(ProbeTracks_dz, BToKEE_kIdx)")
     .Define("KEE_k_dzE_All", "Take(ProbeTracks_dzS, BToKEE_kIdx)")
     .Define("KEE_k_dxy_All", "Take(ProbeTracks_dxy, BToKEE_kIdx)")
     .Define("KEE_k_dxyE_All", "Take(ProbeTracks_dxyS, BToKEE_kIdx)")
     .Define("KEE_k_isPck_All", "Take(ProbeTracks_isPacked, BToKEE_kIdx)")//int
     .Define("e1_eta_All", "BToKEE_fit_l1_eta")
     .Define("e2_eta_All", "BToKEE_fit_l2_eta")
     .Define("e1_phi_All", "BToKEE_fit_l1_phi")
     .Define("e2_phi_All", "BToKEE_fit_l2_phi")
     .Define("e1_pt_All", "BToKEE_fit_l1_pt")
     .Define("e2_pt_All", "BToKEE_fit_l2_pt")
     .Define("e1_dz_All", "Take(Electron_dz, BToKEE_l1Idx)")
     .Define("e2_dz_All", "Take(Electron_dz, BToKEE_l2Idx)")
     .Define("e1_dzE_All", "Take(Electron_dzErr, BToKEE_l1Idx)")
     .Define("e2_dzE_All", "Take(Electron_dzErr, BToKEE_l2Idx)")
     .Define("e1_dxy_All", "Take(Electron_dxy, BToKEE_l1Idx)")
     .Define("e2_dxy_All", "Take(Electron_dxy, BToKEE_l2Idx)")
     .Define("e1_dxyE_All", "Take(Electron_dxyErr, BToKEE_l1Idx)")
     .Define("e2_dxyE_All", "Take(Electron_dxyErr, BToKEE_l2Idx)")
     .Define("e1_isConvVeto_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_convVeto, BToKEE_l1Idx)")//bool
     .Define("e2_isConvVeto_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_convVeto, BToKEE_l2Idx)")//bool
     .Define("e1_seedID_All", "Take(Electron_unBiased, BToKEE_l1Idx)")
     .Define("e2_seedID_All", "Take(Electron_unBiased, BToKEE_l2Idx)")
     .Define("e1_mvaID_All", "Take(Electron_mvaId, BToKEE_l1Idx)")
     .Define("e2_mvaID_All", "Take(Electron_mvaId, BToKEE_l2Idx)")
     .Define("e1_isPF_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l1Idx)")//bool
     .Define("e2_isPF_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l2Idx)")//bool
     .Define("e1_isLowPt_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l1Idx)")//bool
     .Define("e2_isLowPt_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l2Idx)")//bool     
     .Define("e1_isPFoverlap_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l1Idx)")//bool
     .Define("e2_isPFoverlap_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l2Idx)")//bool
     .Define("Idx_sel", bkgSB, { "nBToKEE_All", "e1_isPFoverlap_All", "e2_isPFoverlap_All",  "KEE_fit_mass_All"} )
    
     .Define("nBToKMM_All", "nBToKMuMu")
     .Define("KMM_fit_mass_All", "BToKMuMu_fit_mass")
     .Define("KMM_fit_massErr_All", "BToKMuMu_fit_massErr")
     .Define("KMM_cos2D_All", "BToKMuMu_cos2D")
     .Define("KMM_vtxProb_All", "BToKMuMu_svprob")
     .Define("KMM_l_xy_unc_All", "BToKMuMu_l_xy_unc")
     .Define("KMM_l_xyS_All", "BToKMuMu_l_xy/BToKMuMu_l_xy_unc")
     .Define("KMM_pT_All", "BToKMuMu_fit_pt")
     .Define("KMM_eta_All", "BToKMuMu_fit_eta")
     .Define("KMM_phi_All", "BToKMuMu_fit_phi")
     .Define("KMM_mll_fullfit_All", "BToKMuMu_mll_fullfit")
     .Define("KMM_k_eta_All", "BToKMuMu_fit_k_eta")
     .Define("KMM_k_phi_All", "BToKMuMu_fit_k_phi")
     .Define("KMM_k_pt_All", "BToKMuMu_fit_k_pt")
     .Define("KMM_k_dz_All", "Take(ProbeTracks_dz, BToKMuMu_kIdx)")
     .Define("KMM_dzE_All", "Take(ProbeTracks_dzS, BToKMuMu_kIdx)")
     .Define("KMM_k_dxy_All", "Take(ProbeTracks_dxy, BToKMuMu_kIdx)")
     .Define("KMM_k_dxyE_All", "Take(ProbeTracks_dxyS, BToKMuMu_kIdx)")
     .Define("KMM_k_isPck_All", "Take(ProbeTracks_isPacked, BToKMuMu_kIdx)")//int  
     .Define("m1_eta_All", "BToKMuMu_fit_l1_eta")
     .Define("m2_eta_All", "BToKMuMu_fit_l2_eta")
     .Define("m1_phi_All", "BToKMuMu_fit_l1_phi")
     .Define("m2_phi_All", "BToKMuMu_fit_l2_phi")
     .Define("m1_pt_All", "BToKMuMu_fit_l1_pt")
     .Define("m2_pt_All", "BToKMuMu_fit_l2_pt")
     .Define("m1_dz_All", "Take(Muon_dz, BToKMuMu_l1Idx)")
     .Define("m2_dz_All", "Take(Muon_dz, BToKMuMu_l2Idx)")
     .Define("m1_dzE_All", "Take(Muon_dzErr, BToKMuMu_l1Idx)")
     .Define("m2_dzE_All", "Take(Muon_dzErr, BToKMuMu_l2Idx)")
     .Define("m1_dxy_All", "Take(Muon_dxy, BToKMuMu_l1Idx)")
     .Define("m2_dxy_All", "Take(Muon_dxy, BToKMuMu_l2Idx)")
     .Define("m1_dxyE_All", "Take(Muon_dxyErr, BToKMuMu_l1Idx)")
     .Define("m2_dxyE_All", "Take(Muon_dxyErr, BToKMuMu_l2Idx)")
     .Define("m1_mvaID_All", "Take(Muon_mvaId, BToKMuMu_l1Idx)")
     .Define("m2_mvaID_All", "Take(Muon_mvaId, BToKMuMu_l2Idx)")
     .Define("m1_isPF_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l1Idx)")//bool
     .Define("m2_isPF_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l2Idx)")//bool
     .Define("m1_isTriggering_All", "Take(Muon_isTriggering, BToKMuMu_l1Idx)")//int
     .Define("m2_isTriggering_All", "Take(Muon_isTriggering, BToKMuMu_l2Idx)")//int
     .Define("KMM_nExtraTrg_All", "(nTriggerMuon - m1_isTriggering_All - m2_isTriggering_All)");
  
    std::cout << " totN start = " << *(d.Count())
    	      << " valid KEE triplets = " << *(n.Filter("nBToKEE_All > 0").Count())
    	      << " valid KMM triplets = " << *(n.Filter("nBToKMM_All > 0").Count()) << std::endl;
  }
  
  
  if(isMC && isEE){
    n = n0.Define("isEE_All", Form("%d",isEE))
     .Define("isResonant_All", Form("%d",isResonant))
     .Define("nBToKEE_All", "nBToKEE")
     .Define("KEE_fit_mass_All", "BToKEE_fit_mass")
     .Define("KEE_fit_massErr_All", "BToKEE_fit_massErr")
     .Define("KEE_cos2D_All", "BToKEE_cos2D")
     .Define("KEE_vtxProb_All", "BToKEE_svprob")
     .Define("KEE_l_xy_unc_All", "BToKEE_l_xy_unc")
     .Define("KEE_l_xyS_All", "BToKEE_l_xy/BToKEE_l_xy_unc")
     .Define("KEE_pT_All", "BToKEE_fit_pt")
     .Define("KEE_eta_All", "BToKEE_fit_eta")
     .Define("KEE_phi_All", "BToKEE_fit_phi")
     .Define("KEE_mll_fullfit_All", "BToKEE_mll_fullfit")
     .Define("KEE_k_eta_All", "BToKEE_fit_k_eta")
     .Define("KEE_k_phi_All", "BToKEE_fit_k_phi")
     .Define("KEE_k_pt_All", "BToKEE_fit_k_pt")
     .Define("KEE_k_dz_All", "Take(ProbeTracks_dz, BToKEE_kIdx)")
     .Define("KEE_k_dzE_All", "Take(ProbeTracks_dzS, BToKEE_kIdx)")
     .Define("KEE_k_dxy_All", "Take(ProbeTracks_dxy, BToKEE_kIdx)")
     .Define("KEE_k_dxyE_All", "Take(ProbeTracks_dxyS, BToKEE_kIdx)")
     .Define("KEE_k_isPck_All", "Take(ProbeTracks_isPacked, BToKEE_kIdx)")//int
     .Define("e1_eta_All", "BToKEE_fit_l1_eta")
     .Define("e2_eta_All", "BToKEE_fit_l2_eta")
     .Define("e1_phi_All", "BToKEE_fit_l1_phi")
     .Define("e2_phi_All", "BToKEE_fit_l2_phi")
     .Define("e1_pt_All", "BToKEE_fit_l1_pt")
     .Define("e2_pt_All", "BToKEE_fit_l2_pt")
     .Define("e1_dz_All", "Take(Electron_dz, BToKEE_l1Idx)")
     .Define("e2_dz_All", "Take(Electron_dz, BToKEE_l2Idx)")
     .Define("e1_dzE_All", "Take(Electron_dzErr, BToKEE_l1Idx)")
     .Define("e2_dzE_All", "Take(Electron_dzErr, BToKEE_l2Idx)")
     .Define("e1_dxy_All", "Take(Electron_dxy, BToKEE_l1Idx)")
     .Define("e2_dxy_All", "Take(Electron_dxy, BToKEE_l2Idx)")
     .Define("e1_dxyE_All", "Take(Electron_dxyErr, BToKEE_l1Idx)")
     .Define("e2_dxyE_All", "Take(Electron_dxyErr, BToKEE_l2Idx)")
     .Define("e1_isConvVeto_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_convVeto, BToKEE_l1Idx)")//bool
     .Define("e2_isConvVeto_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_convVeto, BToKEE_l2Idx)")//bool
     .Define("e1_seedID_All", "Take(Electron_unBiased, BToKEE_l1Idx)")
     .Define("e2_seedID_All", "Take(Electron_unBiased, BToKEE_l2Idx)")
     .Define("e1_mvaID_All", "Take(Electron_mvaId, BToKEE_l1Idx)")
     .Define("e2_mvaID_All", "Take(Electron_mvaId, BToKEE_l2Idx)")
     .Define("e1_isPF_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l1Idx)")//bool
     .Define("e2_isPF_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l2Idx)")//bool
     .Define("e1_isLowPt_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l1Idx)")//bool
     .Define("e2_isLowPt_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l2Idx)")//bool     
     .Define("e1_isPFoverlap_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l1Idx)")//bool
     .Define("e2_isPFoverlap_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l2Idx)");//bool

    std::cout << " totN start = " << *(d.Count())
	      << " valid KEE triplets = " << *(n.Filter("nBToKEE_All > 0").Count());
  }

  
  if(isMC && !isEE){
    n = n0.Define("isEE_All", Form("%d",isEE))
     .Define("isResonant_All", Form("%d",isResonant))
     .Define("nBToKMM_All", "nBToKMuMu")
     .Define("KMM_fit_mass_All", "BToKMuMu_fit_mass")
     .Define("KMM_fit_massErr_All", "BToKMuMu_fit_massErr")
     .Define("KMM_cos2D_All", "BToKMuMu_cos2D")
     .Define("KMM_vtxProb_All", "BToKMuMu_svprob")
     .Define("KMM_l_xy_unc_All", "BToKMuMu_l_xy_unc")
     .Define("KMM_l_xyS_All", "BToKMuMu_l_xy/BToKMuMu_l_xy_unc")
     .Define("KMM_pT_All", "BToKMuMu_fit_pt")
     .Define("KMM_eta_All", "BToKMuMu_fit_eta")
     .Define("KMM_phi_All", "BToKMuMu_fit_phi")
     .Define("KMM_mll_fullfit_All", "BToKMuMu_mll_fullfit")
     .Define("KMM_k_eta_All", "BToKMuMu_fit_k_eta")
     .Define("KMM_k_phi_All", "BToKMuMu_fit_k_phi")
     .Define("KMM_k_pt_All", "BToKMuMu_fit_k_pt")
     .Define("KMM_k_dz_All", "Take(ProbeTracks_dz, BToKMuMu_kIdx)")
     .Define("KMM_dzE_All", "Take(ProbeTracks_dzS, BToKMuMu_kIdx)")
     .Define("KMM_k_dxy_All", "Take(ProbeTracks_dxy, BToKMuMu_kIdx)")
     .Define("KMM_k_dxyE_All", "Take(ProbeTracks_dxyS, BToKMuMu_kIdx)")
     .Define("KMM_k_isPck_All", "Take(ProbeTracks_isPacked, BToKMuMu_kIdx)")//int  
     .Define("m1_eta_All", "BToKMuMu_fit_l1_eta")
     .Define("m2_eta_All", "BToKMuMu_fit_l2_eta")
     .Define("m1_phi_All", "BToKMuMu_fit_l1_phi")
     .Define("m2_phi_All", "BToKMuMu_fit_l2_phi")
     .Define("m1_pt_All", "BToKMuMu_fit_l1_pt")
     .Define("m2_pt_All", "BToKMuMu_fit_l2_pt")
     .Define("m1_dz_All", "Take(Muon_dz, BToKMuMu_l1Idx)")
     .Define("m2_dz_All", "Take(Muon_dz, BToKMuMu_l2Idx)")
     .Define("m1_dzE_All", "Take(Muon_dzErr, BToKMuMu_l1Idx)")
     .Define("m2_dzE_All", "Take(Muon_dzErr, BToKMuMu_l2Idx)")
     .Define("m1_dxy_All", "Take(Muon_dxy, BToKMuMu_l1Idx)")
     .Define("m2_dxy_All", "Take(Muon_dxy, BToKMuMu_l2Idx)")
     .Define("m1_dxyE_All", "Take(Muon_dxyErr, BToKMuMu_l1Idx)")
     .Define("m2_dxyE_All", "Take(Muon_dxyErr, BToKMuMu_l2Idx)")
     .Define("m1_mvaID_All", "Take(Muon_mvaId, BToKMuMu_l1Idx)")
     .Define("m2_mvaID_All", "Take(Muon_mvaId, BToKMuMu_l2Idx)")
     .Define("m1_isPF_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l1Idx)")//bool
     .Define("m2_isPF_All", "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l2Idx)")//bool
     .Define("m1_isTriggering_All", "Take(Muon_isTriggering, BToKMuMu_l1Idx)")//int
     .Define("m2_isTriggering_All", "Take(Muon_isTriggering, BToKMuMu_l2Idx)")//int
     .Define("KMM_nExtraTrg_All", "(nTriggerMuon - m1_isTriggering_All - m2_isTriggering_All)");

    std::cout << " totN start = " << *(d.Count())
            << " valid KMM triplets = " << *(n.Filter("nBToKMM_All > 0").Count()) << std::endl;
  }  

  
  auto tree_All = n;
    
  if(!isMC){

    tree_All = n.Define("good_KEE_All", KEEcut, {"nBToKEE_All", "e1_pt_All", "e2_pt_All", "KEE_k_pt_All",
	  "KEE_cos2D_All", "KEE_vtxProb_All", "KEE_l_xyS_All", "KEE_l_xy_unc_All", "KEE_pT_All", "KEE_eta_All",
	  "e1_isPFoverlap_All", "e2_isPFoverlap_All"})
            .Define("KEE_rkVtx_All", flagReverseRank, {"nBToKEE_All", "good_KEE_All", "KEE_vtxProb_All"})
            .Define("good_KMM_All", KMMcut, {"nBToKMM_All", "m1_pt_All", "m2_pt_All", "KMM_k_pt_All", "KMM_nExtraTrg_All",
	   "KMM_cos2D_All", "KMM_vtxProb_All", "KMM_l_xyS_All", "KMM_l_xy_unc_All", "KMM_pT_All", "KMM_eta_All"})
            .Define("KMM_rkVtx_All", flagReverseRank, {"nBToKMM_All", "good_KMM_All", "KMM_vtxProb_All"});

  }

  if(isMC){

    if(isEE){  

      tree_All = n.Define("good_KEE_All", KEEcut, {"nBToKEE_All", "e1_pt_All", "e2_pt_All", "KEE_k_pt_All",
	    "KEE_cos2D_All", "KEE_vtxProb_All", "KEE_l_xyS_All", "KEE_l_xy_unc_All", "KEE_pT_All", "KEE_eta_All",
	    "e1_isPFoverlap_All", "e2_isPFoverlap_All"})
              .Define("KEE_rkVtx_All", flagReverseRank, {"nBToKEE_All", "good_KEE_All", "KEE_vtxProb_All"});

    }
    else{

      tree_All = n.Define("good_KMM_All", KMMcut, {"nBToKMM_All", "m1_pt_All", "m2_pt_All", "KMM_k_pt_All", "KMM_nExtraTrg_All",
	    "KMM_cos2D_All", "KMM_vtxProb_All", "KMM_l_xyS_All", "KMM_l_xy_unc_All", "KMM_pT_All", "KMM_eta_All"})
	      .Define("KMM_rkVtx_All", flagReverseRank, {"nBToKMM_All", "good_KMM_All", "KMM_vtxProb_All"});

    }

  }


  if(isMC){
        
    std::string l1_eta = isEE ? "e1_eta_All" : "m1_eta_All";
    std::string l2_eta = isEE ? "e2_eta_All" : "m2_eta_All";
    std::string k_eta = isEE ? "KEE_k_eta_All" : "KMM_k_eta_All";
    std::string l1_phi = isEE ? "e1_phi_All" : "m1_phi_All";
    std::string l2_phi = isEE ? "e2_phi_All" : "m2_phi_All";
    std::string k_phi = isEE ? "KEE_k_phi_All" : "KMM_k_phi_All";
    std::string n_Btriplet = isEE ? "nBToKEE_All" : "nBToKMM_All";

    // 443 = JPsi    521 = B+  321 = K 
    std::string B_l1_genParent = isEE ? "Take(Electron_genPartFlav, BToKEE_l1Idx)" : "Take(Muon_genPartFlav, BToKMuMu_l1Idx)";
    std::string B_l2_genParent = isEE ? "Take(Electron_genPartFlav, BToKEE_l2Idx)" : "Take(Muon_genPartFlav, BToKMuMu_l2Idx)";
    std::string B_k_genParent = isEE ? "Take(ProbeTracks_genPartFlav, BToKEE_kIdx)" : "Take(ProbeTracks_genPartFlav, BToKMuMu_kIdx)";
    std::string GenPart_l1_idx = isEE ? "Take(Electron_genPartIdx, BToKEE_l1Idx)" : "Take(Muon_genPartIdx, BToKMuMu_l1Idx)";
    std::string GenPart_l2_idx = isEE ? "Take(Electron_genPartIdx, BToKEE_l2Idx)" : "Take(Muon_genPartIdx, BToKMuMu_l2Idx)";
    std::string GenPart_k_idx = isEE ? "Take(ProbeTracks_genPartIdx, BToKEE_kIdx)" : "Take(ProbeTracks_genPartIdx, BToKMuMu_kIdx)";
    
    auto mc  = tree_All.Define("B_l1_genParent_GenAll", B_l1_genParent.c_str())
      .Define("B_l2_genParent_GenAll", B_l2_genParent.c_str())
      .Define("B_k_genParent_GenAll", B_k_genParent.c_str())
      .Define("GenPart_l1_idx_GenAll", GenPart_l1_idx.c_str())
      .Define("GenPart_l2_idx_GenAll", GenPart_l2_idx.c_str())
      .Define("GenPart_k_idx_GenAll", GenPart_k_idx.c_str())
      .Define("GenPart_l1_pdgId_GenAll", "Take(GenPart_pdgId, GenPart_l1_idx_GenAll)")
      .Define("GenPart_l2_pdgId_GenAll", "Take(GenPart_pdgId, GenPart_l2_idx_GenAll)")
      .Define("GenPart_k_pdgId_GenAll", "Take(GenPart_pdgId, GenPart_k_idx_GenAll)")
      .Define("GenMothPart_l1_idx_GenAll", "Take(GenPart_genPartIdxMother, GenPart_l1_idx_GenAll)")
      .Define("GenMothPart_l2_idx_GenAll", "Take(GenPart_genPartIdxMother, GenPart_l2_idx_GenAll)")
      .Define("GenMothPart_k_idx_GenAll", "Take(GenPart_genPartIdxMother, GenPart_k_idx_GenAll)")
      .Define("GenMothPart_l1_pdgId_GenAll", "Take(GenPart_pdgId, GenMothPart_l1_idx_GenAll)")
      .Define("GenMothPart_l2_pdgId_GenAll", "Take(GenPart_pdgId, GenMothPart_l2_idx_GenAll)")
      .Define("GenMothPart_k_pdgId_GenAll", "Take(GenPart_pdgId, GenMothPart_k_idx_GenAll)")
      .Define("GenGMothPart_l1_idx_GenAll", "Take(GenPart_genPartIdxMother, GenMothPart_l1_idx_GenAll)")
      .Define("GenGMothPart_l2_idx_GenAll", "Take(GenPart_genPartIdxMother, GenMothPart_l2_idx_GenAll)")
      .Define("GenGMothPart_k_idx_GenAll", "Take(GenPart_genPartIdxMother, GenMothPart_k_idx_GenAll)")
      .Define("GenGMothPart_l1_pdgId_GenAll", "Take(GenPart_pdgId, GenGMothPart_l1_idx_GenAll)") 
      .Define("GenGMothPart_l2_pdgId_GenAll", "Take(GenPart_pdgId, GenGMothPart_l2_idx_GenAll)") 
      .Define("GenGMothPart_k_pdgId_GenAll", "Take(GenPart_pdgId, GenGMothPart_k_idx_GenAll)")
      .Define("GenPart_l1_eta_GenAll", "Take(GenPart_eta, GenPart_l1_idx_GenAll)")
      .Define("GenPart_l2_eta_GenAll", "Take(GenPart_eta, GenPart_l2_idx_GenAll)")
      .Define("GenPart_k_eta_GenAll", "Take(GenPart_eta, GenPart_k_idx_GenAll)")
      .Define("GenPart_l1_phi_GenAll", "Take(GenPart_phi, GenPart_l1_idx_GenAll)")
      .Define("GenPart_l2_phi_GenAll", "Take(GenPart_phi, GenPart_l2_idx_GenAll)")
      .Define("GenPart_k_phi_GenAll", "Take(GenPart_phi, GenPart_k_idx_GenAll)")
      .Define("l1_eta_all", l1_eta.c_str())
      .Define("l2_eta_all", l2_eta.c_str())
      .Define("k_eta_all", k_eta.c_str())
      .Define("l1_phi_all", l1_phi.c_str())
      .Define("l2_phi_all", l2_phi.c_str())
      .Define("k_phi_all", k_phi.c_str())
      .Define("n_Btriplet_all", n_Btriplet.c_str());


    auto mcGenMatched = mc.Define("isGenMatched_All", flagGenMatchExt, {"isResonant_All",
          "GenPart_l1_pdgId_GenAll", "GenPart_l2_pdgId_GenAll", "GenPart_k_pdgId_GenAll",
          "GenMothPart_l1_pdgId_GenAll", "GenMothPart_l2_pdgId_GenAll", "GenMothPart_k_pdgId_GenAll",
          "GenGMothPart_l1_pdgId_GenAll", "GenGMothPart_l2_pdgId_GenAll", "GenGMothPart_k_pdgId_GenAll"})
                          .Define("dRwithGen_All", computedR, {"isGenMatched_All",
				 "GenPart_l1_eta_GenAll", "GenPart_l2_eta_GenAll", "GenPart_k_eta_GenAll",
				 "l1_eta_all", "l2_eta_all", "k_eta_all",
				 "GenPart_l1_phi_GenAll", "GenPart_l2_phi_GenAll", "GenPart_k_phi_GenAll",
				 "l1_phi_all", "l2_phi_all", "k_phi_all"})
                          .Define("rank_gendR_All", flagRank, {"n_Btriplet_all", "isGenMatched_All", "dRwithGen_All"});
                          
    auto mcGenMatched_i = mcGenMatched;

    if(isEE)
        mcGenMatched_i = mcGenMatched.Define("Idx_sel", si3sg, { "nBToKEE_All", "e1_isPFoverlap_All", 
                                             "e2_isPFoverlap_All",  "KEE_fit_mass_All", "rank_gendR_All"} );

    //give eta, phi, l1, l2, k e gen
    std::cout << "\n computed dR for genmatched " << std::endl;

    if(outFName == "-1")
      mcGenMatched_i.Snapshot("newtree", Form("output_RDF_Kll_isMC%d_isEE%d_BDT.root", isMC, isEE), "\\b([^ ]*)(_All)|\\b([^ ]*)(_sel)");
    else 
      mcGenMatched_i.Snapshot("newtree", outFName.c_str(), "\\b([^ ]*)(_All)|\\b([^ ]*)(_sel)");
    //all triplets with Mc gen matched info

  }//isMC
  else{

    if(outFName == "-1")
      tree_All.Snapshot("newtree", Form("output_RDF_Kll_isMC%d_BDT.root", isMC), "\\b([^ ]*)(_All)|\\b([^ ]*)(_sel)");
    else
      tree_All.Snapshot("newtree", outFName.c_str(), "\\b([^ ]*)(_All)|\\b([^ ]*)(_sel)");
  }

  
  t.Stop(); 
  t.Print();

} 
