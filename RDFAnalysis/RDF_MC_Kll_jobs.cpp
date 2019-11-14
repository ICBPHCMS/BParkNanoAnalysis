//g++ -Wall -o RDF_MC_Kll_jobs `root-config --cflags --glibs ` RDF_MC_Kll_jobs.cpp 

//./RDF_MC_Kll_jobs --JOBid (-1, 1,2..) --inList input_list.txt --outFile outFileName --isMC (0, 1) --isEE (0, 1) --isResonant (0, 1) --testFile /path/to/input_file.root

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <string>
#include "TStopwatch.h"
#include <vector>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include "interface/functions.h"


#include <regex>

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
    std::cout << " configuration ERROR => splitting file based but missing JOBid and output file path " << std::endl;
    return -1;
  }

  if(inList == "-1" && testFile == "-1"){
    std::cout << " by default " << std::endl;
  }

  std::cout << " inList = " << inList << " JOBid = " << JOBid << " outFName = " << outFName << " testFile = " << testFile
	    << " isMC = " << isMC << " isEE = " << isEE << "isResonant = " << isResonant << "\n" << std::endl;

  std::vector<std::string> inputFileList;
  
  if(inList != "-1"){
    std::string rootFileName;
    std::ifstream inFileLong;
    inFileLong.open(inList.c_str(), std::ios::in);
    while(!inFileLong.eof()){
      if(inFileLong >> rootFileName){
	inputFileList.push_back(rootFileName);
	//	ch.Add(rootFileName.c_str());
	std::cout << " adding " << rootFileName << std::endl;
      }
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


  //  for(auto ij : inputFileList) std::cout << " file = " << ij << std::endl;

  //  return 10;
  //  std::cout << " isMC = " << isMC << " isEE = " << isEE << std::endl;
  //  std::cout << " doCutBased = " << doCutBased << " isResonant = " << isResonant << std::endl;

  //ROOT::RDataFrame df(ch);
  ROOT::RDataFrame d("Events", inputFileList);

  std::string nB = isEE ? "nBToKEE" : "nBToKMuMu";
  // std::string isElectronChannel = isEE ? "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 1)" : 
  //                                        "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)";
  //  std::string event = "ROOT::VecOps::RVec<int> (nBtriplet_All, event)";
  std::string l1Trg = isEE ? "ROOT::VecOps::RVec<int> (nBtriplet_All, -1)" : "Take(Muon_isTriggering, BToKMuMu_l1Idx)";
  std::string l2Trg = isEE ? "ROOT::VecOps::RVec<int> (nBtriplet_All, -1)" : "Take(Muon_isTriggering, BToKMuMu_l2Idx)";
  // for ele as nTriggerMuon flattened over triplets
  std::string nTrg = isEE ? "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, nTriggerMuon)" : "(nTriggerMuon - B_l1_isTriggering_All - B_l2_isTriggering_All)"; 
  // std::string B_l1_pT = isEE ? "Take(Electron_pt, BToKEE_l1Idx)" : "Take(Muon_pt, BToKMuMu_l1Idx)";
  // std::string B_l2_pT = isEE ? "Take(Electron_pt, BToKEE_l2Idx)" : "Take(Muon_pt, BToKMuMu_l2Idx)";
  // std::string B_k_pT = isEE ? "Take(ProbeTracks_pt, BToKEE_kIdx)" : "Take(ProbeTracks_pt, BToKMuMu_kIdx)";
  std::string B_fit_mass = isEE ? "BToKEE_fit_mass" : "BToKMuMu_fit_mass";
  std::string B_fit_massErr = isEE ? "BToKEE_fit_massErr" : "BToKMuMu_fit_massErr";
  std::string B_cos2D = isEE ? "BToKEE_cos2D" : "BToKMuMu_cos2D";
  std::string B_vtxProb = isEE ? "BToKEE_svprob" : "BToKMuMu_svprob";
  std::string B_pT = isEE ? "BToKEE_fit_pt" : "BToKMuMu_fit_pt";
  std::string B_eta = isEE ? "BToKEE_fit_eta" : "BToKMuMu_fit_eta";
  std::string B_phi = isEE ? "BToKEE_fit_phi" : "BToKMuMu_fit_phi";
  //std::string B_mll_llfit = isEE ? "BToKEE_mll_llfit" : "BToKMuMu_mll_llfit";
  std::string B_mll_fullfit = isEE ? "BToKEE_mll_fullfit" : "BToKMuMu_mll_fullfit";
  std::string B_l_xy_unc = isEE ? "BToKEE_l_xy_unc" : "BToKMuMu_l_xy_unc";
  std::string B_l_xyS = isEE ? "BToKEE_l_xy/BToKEE_l_xy_unc" : "BToKMuMu_l_xy/BToKMuMu_l_xy_unc";
  // can add further requirements on lepton ID
  std::string B_l1_eta = isEE ? "BToKEE_fit_l1_eta" : "BToKMuMu_fit_l1_eta";
  std::string B_l2_eta = isEE ? "BToKEE_fit_l2_eta" : "BToKMuMu_fit_l2_eta";
  std::string B_k_eta = isEE ? "BToKEE_fit_k_eta" : "BToKMuMu_fit_k_eta";
  std::string B_l1_phi = isEE ? "BToKEE_fit_l1_phi" : "BToKMuMu_fit_l1_phi";
  std::string B_l2_phi = isEE ? "BToKEE_fit_l2_phi" : "BToKMuMu_fit_l2_phi";
  std::string B_k_phi = isEE ? "BToKEE_fit_k_phi" : "BToKMuMu_fit_k_phi";
  std::string B_l1_pt = isEE ? "BToKEE_fit_l1_pt" : "BToKMuMu_fit_l1_pt";
  std::string B_l2_pt = isEE ? "BToKEE_fit_l2_pt" : "BToKMuMu_fit_l2_pt";
  std::string B_k_pt = isEE ? "BToKEE_fit_k_pt" : "BToKMuMu_fit_k_pt";

  std::string B_l1_dz = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_dz, BToKEE_l1Idx)":
                               "(ROOT::VecOps::RVec<float>) Take(Muon_dz, BToKMuMu_l1Idx)";
  std::string B_l2_dz = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_dz, BToKEE_l2Idx)":
                               "(ROOT::VecOps::RVec<float>) Take(Muon_dz, BToKMuMu_l2Idx)";
  std::string B_k_dz = isEE ? "(ROOT::VecOps::RVec<float>) Take(ProbeTracks_dz, BToKEE_kIdx)":
                               "(ROOT::VecOps::RVec<float>) Take(ProbeTracks_dz, BToKMuMu_kIdx)";
  std::string B_l1_dzE = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_dzErr, BToKEE_l1Idx)":
                               "(ROOT::VecOps::RVec<float>) Take(Muon_dzErr, BToKMuMu_l1Idx)";
  std::string B_l2_dzE = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_dzErr, BToKEE_l2Idx)":
                               "(ROOT::VecOps::RVec<float>) Take(Muon_dzErr, BToKMuMu_l2Idx)";
  std::string B_k_dzE = isEE ? "(ROOT::VecOps::RVec<float>) Take(ProbeTracks_dzS, BToKEE_kIdx)":
                               "(ROOT::VecOps::RVec<float>) Take(ProbeTracks_dzS, BToKMuMu_kIdx)";
  std::string B_l1_dxy = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_dxy, BToKEE_l1Idx)":
                               "(ROOT::VecOps::RVec<float>) Take(Muon_dxy, BToKMuMu_l1Idx)";
  std::string B_l2_dxy = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_dxy, BToKEE_l2Idx)":
                               "(ROOT::VecOps::RVec<float>) Take(Muon_dxy, BToKMuMu_l2Idx)";
  std::string B_k_dxy = isEE ? "(ROOT::VecOps::RVec<float>) Take(ProbeTracks_dxy, BToKEE_kIdx)":
                               "(ROOT::VecOps::RVec<float>) Take(ProbeTracks_dxy, BToKMuMu_kIdx)";
  std::string B_l1_dxyE = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_dxyErr, BToKEE_l1Idx)":
                               "(ROOT::VecOps::RVec<float>) Take(Muon_dxyErr, BToKMuMu_l1Idx)";
  std::string B_l2_dxyE = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_dxyErr, BToKEE_l2Idx)":
                               "(ROOT::VecOps::RVec<float>) Take(Muon_dxyErr, BToKMuMu_l2Idx)";
  std::string B_k_dxyE = isEE ? "(ROOT::VecOps::RVec<float>) Take(ProbeTracks_dxyS, BToKEE_kIdx)":
                               "(ROOT::VecOps::RVec<float>) Take(ProbeTracks_dxyS, BToKMuMu_kIdx)";
  //ID
  std::string B_l1_isConvVeto = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_convVeto, BToKEE_l1Idx)" : 
                                       "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)";
  std::string B_l2_isConvVeto = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_convVeto, BToKEE_l2Idx)" : 
                                       "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)";
  std::string B_l1_seedID = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_unBiased, BToKEE_l1Idx)" : 
                                   "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)";
  std::string B_l2_seedID = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_unBiased, BToKEE_l2Idx)" : 
                                   "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)";
  std::string B_l1_mvaID = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_mvaId, BToKEE_l1Idx)" : 
                                  "(ROOT::VecOps::RVec<float>) Take(Muon_mvaId, BToKMuMu_l1Idx)";
  std::string B_l2_mvaID = isEE ? "(ROOT::VecOps::RVec<float>) Take(Electron_mvaId, BToKEE_l2Idx)" : 
                                  "(ROOT::VecOps::RVec<float>) Take(Muon_mvaId, BToKMuMu_l2Idx)";

  std::string B_l1_isPF = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l1Idx)" : 
                                 "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l1Idx)"; 
  std::string B_l2_isPF = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l2Idx)" : 
                                 "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l2Idx)";
  std::string B_k_isPck = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(ProbeTracks_isPacked, BToKEE_kIdx)" : 
                                 "(ROOT::VecOps::RVec<unsigned int>) Take(ProbeTracks_isPacked, BToKMuMu_kIdx)";

  // std::string B_l1_isLowPt = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l1Idx)" : 
  //                                "(ROOT::VecOps::RVec<unsigned int>) (nBtriplet_All, 0)"; 
  // std::string B_l2_isLowPt = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l2Idx)" : 
  //                                "(ROOT::VecOps::RVec<unsigned int>) (nBtriplet_All, 0)";
  std::string B_l1_isPFoverlap = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l1Idx)" : 
                                        "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)"; 
  std::string B_l2_isPFoverlap = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l2Idx)" : 
                                        "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)"; 
  std::string weights = "ROOT::VecOps::RVec<float>(nBtriplet_All, 1.)";


  auto n = d.Define("lumi_All", "luminosityBlock")
    .Define("eventToT_All", "event")
    .Define("run_All", "run")

    .Define("isResonant_All", (std::string(argv[4])).c_str())
    .Define("isEE_All", (std::string(argv[2])).c_str())

    .Define("nBtriplet_All", nB.c_str())
    .Define("B_l1_isTriggering_All", l1Trg.c_str())
    .Define("B_l2_isTriggering_All", l2Trg.c_str())
    .Define("nExtraTrg_All", nTrg.c_str())

    .Define("B_fit_mass_All", B_fit_mass.c_str())
    .Define("B_fit_massErr_All", B_fit_massErr.c_str())
    .Define("B_cos2D_All", B_cos2D.c_str())
    .Define("B_vtxProb_All", B_vtxProb.c_str())
    .Define("B_l_xy_unc_All", B_l_xy_unc.c_str())
    .Define("B_l_xyS_All", B_l_xyS.c_str())
    .Define("B_pT_All", B_pT.c_str())
    .Define("B_eta_All", B_eta.c_str())
    .Define("B_phi_All", B_phi.c_str())
    .Define("B_mll_fullfit_All", B_mll_fullfit.c_str())
    //    .Define("B_mll_llfit_All", B_mll_llfit.c_str())

    .Define("B_l1_eta_All", B_l1_eta.c_str())
    .Define("B_l2_eta_All", B_l2_eta.c_str())
    .Define("B_k_eta_All", B_k_eta.c_str())
    .Define("B_l1_phi_All", B_l1_phi.c_str())
    .Define("B_l2_phi_All", B_l2_phi.c_str())
    .Define("B_k_phi_All", B_k_phi.c_str())
    .Define("B_l1_pt_All", B_l1_pt.c_str())
    .Define("B_l2_pt_All", B_l2_pt.c_str())
    .Define("B_k_pt_All", B_k_pt.c_str())
  
    .Define("B_l1_dz_All", B_l1_dz.c_str())
    .Define("B_l2_dz_All", B_l2_dz.c_str())
    .Define("B_k_dz_All", B_k_dz.c_str())
    .Define("B_l1_dzE_All", B_l1_dzE.c_str())
    .Define("B_l2_dzE_All", B_l2_dzE.c_str())
    .Define("B_k_dzE_All", B_k_dzE.c_str())
    .Define("B_l1_dxy_All", B_l1_dxy.c_str())
    .Define("B_l2_dxy_All", B_l2_dxy.c_str())
    .Define("B_k_dxy_All", B_k_dxy.c_str())
    .Define("B_l1_dxyE_All", B_l1_dxyE.c_str())
    .Define("B_l2_dxyE_All", B_l2_dxyE.c_str())
    .Define("B_k_dxyE_All", B_k_dxyE.c_str())
         
    .Define("B_l1_isConvVeto_All", B_l1_isConvVeto.c_str())
    .Define("B_l2_isConvVeto_All", B_l2_isConvVeto.c_str())
    .Define("B_l1_seedID_All", B_l1_seedID.c_str())
    .Define("B_l2_seedID_All", B_l2_seedID.c_str())
    .Define("B_l1_mvaID_All", B_l1_mvaID.c_str())
    .Define("B_l2_mvaID_All", B_l2_mvaID.c_str())

    .Define("B_l1_isPF_All", B_l1_isPF.c_str())
    .Define("B_l2_isPF_All", B_l2_isPF.c_str())
    .Define("B_k_isPck_All", B_k_isPck.c_str())
    // .Define("B_l1_isLowPt_All", B_l1_isLowPt.c_str())
    // .Define("B_l2_isLowPt_All", B_l2_isLowPt.c_str())
    .Define("B_l1_isPFoverlap_All", B_l1_isPFoverlap.c_str())
    .Define("B_l2_isPFoverlap_All", B_l2_isPFoverlap.c_str())
    .Define("weights_All", weights.c_str());

  

  //define a branch flagging the good candidates and filter events with 0
  auto tree_All = n.Define("cutBase_goodB_All", cutBased, {"nBtriplet_All", "B_l1_pt_All", "B_l2_pt_All", "B_k_pt_All", "nExtraTrg_All", 
	"B_cos2D_All", "B_vtxProb_All", "B_l_xyS_All", "B_l_xy_unc_All", "B_pT_All", "B_eta_All", 
	"B_l1_isPFoverlap_All", "B_l2_isPFoverlap_All"})
    .Define("rankVtx_All", flagReverseRank, {"nBtriplet_All", "cutBase_goodB_All", "B_vtxProb_All"});
  // //auto selected = n2.Filter("Any(good_B == true)", "goodB");


  /*
  std::vector<std::string> listColumns_All = {"B_fit_mass_All", "cutBase_goodB_All", "rankVtx_All", 
					      "B_l1_pT_All", "B_l2_pT_All", "B_k_pT_All", "B_l_xyS_All", 
					      "B_cos2D_All", "B_vtxProb_All", "B_pT_All", "B_eta_All",
					      "B_mll_fullfit_All", "B_mll_llfit_All", "B_l_xy_unc_All", 
					      "nBtriplet_All", "eventToT_All", 
					      "B_l1_eta_All", "B_l2_eta_All", "B_k_eta_All", "B_l1_phi_All",
					      "B_l2_phi_All", "B_k_phi_All", "B_l1_pt_All", "B_l2_pt_All", "B_k_pt_All",
					      "B_l1_isPF_All", "B_l2_isPF_All", 
					      "B_l1_isLowPt_All", "B_l2_isLowPt_All", 
					      "B_l1_isPFoverlap_All", "B_l2_isPFoverlap_All"};
  */

  if(isMC){
    // listColumns_All.push_back("isGenMatched_All");
    // listColumns_All.push_back("dRwithGen_All");
    // listColumns_All.push_back("rank_gendR_All");
        
    // 443 = JPsi    521 = B+  321 = K 
    std::string B_l1_genParent = isEE ? "Take(Electron_genPartFlav, BToKEE_l1Idx)" : "Take(Muon_genPartFlav, BToKMuMu_l1Idx)";
    std::string B_l2_genParent = isEE ? "Take(Electron_genPartFlav, BToKEE_l2Idx)" : "Take(Muon_genPartFlav, BToKMuMu_l2Idx)";
    std::string B_k_genParent = isEE ? "Take(ProbeTracks_genPartFlav, BToKEE_kIdx)" : "Take(ProbeTracks_genPartFlav, BToKMuMu_kIdx)";
    std::string GenPart_l1_idx = isEE ? "Take(Electron_genPartIdx, BToKEE_l1Idx)" : "Take(Muon_genPartIdx, BToKMuMu_l1Idx)";
    std::string GenPart_l2_idx = isEE ? "Take(Electron_genPartIdx, BToKEE_l2Idx)" : "Take(Muon_genPartIdx, BToKMuMu_l2Idx)";
    std::string GenPart_k_idx = isEE ? "Take(ProbeTracks_genPartIdx, BToKEE_kIdx)" : "Take(ProbeTracks_genPartIdx, BToKMuMu_kIdx)";
    
    auto mc = tree_All.Define("B_l1_genParent_GenAll", B_l1_genParent.c_str())
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
      .Define("GenPart_k_phi_GenAll", "Take(GenPart_phi, GenPart_k_idx_GenAll)");


    std::cout << " flagGenMatched" << std::endl;
    auto mcGenMatched = mc.Define("isGenMatched_All", flagGenMatchExt, {"isResonant_All", 
	  "GenPart_l1_pdgId_GenAll", "GenPart_l2_pdgId_GenAll", "GenPart_k_pdgId_GenAll", 
	  "GenMothPart_l1_pdgId_GenAll", "GenMothPart_l2_pdgId_GenAll", "GenMothPart_k_pdgId_GenAll", 
	  "GenGMothPart_l1_pdgId_GenAll", "GenGMothPart_l2_pdgId_GenAll", "GenGMothPart_k_pdgId_GenAll"})
      .Define("dRwithGen_All", computedR, {"isGenMatched_All", 
	    "GenPart_l1_eta_GenAll", "GenPart_l2_eta_GenAll", "GenPart_k_eta_GenAll", 
	    "B_l1_eta_All", "B_l2_eta_All", "B_k_eta_All", 
	    "GenPart_l1_phi_GenAll", "GenPart_l2_phi_GenAll", "GenPart_k_phi_GenAll", 
	    "B_l1_phi_All", "B_l2_phi_All", "B_k_phi_All"})
      .Define("rank_gendR_All", flagRank, {"nBtriplet_All", "isGenMatched_All", "dRwithGen_All"});

    //give eta, phi, l1, l2, k e gen
    std::cout << " computed dR for genmatched " << std::endl;

    //    mcGenMatched.Snapshot("newtree", Form("newfile_isMC%d_isEE%d_All.root", isMC, isEE), listColumns_All);
    //    mcGenMatched.Snapshot("newtree", Form("newfile_isMC%d_isEE%d_All.root", isMC, isEE), "\\b([^ ]*)(_All)");

    if(outFName == "-1")
      mcGenMatched.Snapshot("newtree", Form("output_RDF_Kll_isMC%d_isEE%d_Job.root", isMC, isEE), "\\b([^ ]*)(_All)");
    else 
      mcGenMatched.Snapshot("newtree", outFName.c_str(), "\\b([^ ]*)(_All)");
    //all triplets with Mc gen matched info

  }//isMC
  else{

    /*
    if(doCutBased){
      std::vector<std::string> listColumns;
      auto filtered = apply_CutBased(isMC, tree_All, listColumns);
      filtered.Snapshot("newtree", Form("newfile_isMC%d_isEE%d_CB.root", isMC, isEE), listColumns);
    }
    */

    if(outFName == "-1")
      tree_All.Snapshot("newtree", Form("output_RDF_Kll_isMC%d_isEE%d_Job.root", isMC, isEE), "\\b([^ ]*)(_All)");
    else
      tree_All.Snapshot("newtree", outFName.c_str(), "\\b([^ ]*)(_All)");
  }


  std::cout << " totN start = " << *(d.Count()) 
	    << " valid triplets = " << *(n.Filter("nBtriplet_All > 0").Count()) << std::endl;
  
  // your code goes here 
  t.Stop(); 
  t.Print();


}




// next steps:
// - add branch for MC matching
// configure inputs/ selections/ and type: MC vs DATA ee vs mumu
// save in parallel trees for BDT training 
// ...



//some links
//https://root.cern/doc/master/classROOT_1_1RDF_1_1RInterface.html#a233b7723e498967f4340705d2c4db7f8
//https://root.cern.ch/doc/master/namespaceROOT_1_1VecOps.html#a7dcd060b97f6c82621ba0d8f376ad195
//https://root.cern.ch/doc/master/df004__cutFlowReport_8C.html
//https://github.com/arizzi/nail/blob/master/vbfHmumuAna.py

