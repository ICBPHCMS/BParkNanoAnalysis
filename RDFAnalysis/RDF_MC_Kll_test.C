#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <string>
#include "TStopwatch.h"
#include <vector>


//flag triplet passing the selection criteria
ROOT::VecOps::RVec<bool> IsGood(unsigned int nB,
				ROOT::VecOps::RVec<float>& pT1, ROOT::VecOps::RVec<float>& pT2, ROOT::VecOps::RVec<float>& pTk,
				ROOT::VecOps::RVec<unsigned int> nTrg, ROOT::VecOps::RVec<float>& cos2D, ROOT::VecOps::RVec<float>& vtxP,
				ROOT::VecOps::RVec<float>& disp, ROOT::VecOps::RVec<float>& dispU, ROOT::VecOps::RVec<float>& pT, ROOT::VecOps::RVec<float>& eta, 
				ROOT::VecOps::RVec<bool>& l1_ok, ROOT::VecOps::RVec<bool>& l2_ok) {

  ROOT::VecOps::RVec<bool> goodB(nB, false);
  for (auto ij=0; ij<nB; ++ij){
    if(pT1[ij] > 1. && pT2[ij] > 0.5 && pTk[ij] > 0.8 &&
       nTrg[ij] > 0 && cos2D[ij] > 0.99 && vtxP[ij] > 0.1 && disp[ij] > 2 && dispU[ij] != 0 && 
       pT[ij] > 3. && std::abs(eta[ij]) < 2.4 && l1_ok[ij] == true && l2_ok[ij] == true)
      goodB[ij]  = true;
  }

  return goodB;
}


//save for each good triplet the rank - choice is now wrt vtxCL - 
//useful to plot just 1 triplet per event 
ROOT::VecOps::RVec<int> Rank(unsigned int nB, 
			     ROOT::VecOps::RVec<bool>& goodB, ROOT::VecOps::RVec<float>& vtxP){

  auto sortIndices = Argsort(vtxP);
  ROOT::VecOps::RVec<int> rank;
  rank.resize(nB);
  int nRank = 0;
  for (auto ij=0; ij<nB; ++ij){

    if(goodB[sortIndices[ij]]){
      rank[sortIndices[ij]] = nRank;
      ++nRank;
    }
    else{
      rank[sortIndices[ij]] = -1;
    }
  }
  return rank; 
}



using namespace ROOT::VecOps;

//root RDF_MC_Kll_test.C'()' 
void RDF_MC_Kll_test(int isMC, int isEE){

  TStopwatch t; 
  t.Start(); 

  std::string inputFileList = "/eos/cms//store/group/cmst3/group/bpark/BParkingNANO_2019Oct14/";
  if(isMC && !isEE) inputFileList += "BuToKJpsi_ToMuMu_probefilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_ToMuMu/191014_075537/0000/BParkNANO_mc_2019Oct14_*.root";
  else if(isMC && isEE) inputFileList += "BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee/191014_075238/0000/BParkNANO_mc_2019Oct14_*.root";


  std::cout << " isMC = " << isMC << "isEE = " << isEE << std::endl;
  std::cout << "inputFileList = " << inputFileList << std::endl;

  ROOT::RDataFrame d("Events", inputFileList.c_str());


  std::string isElectronChannel = isEE ? "1" : "0";
  std::string nB = isEE ? "nBToKEE" : "nBToKMuMu";
  std::string l1Trg = isEE ? "-1" : "Take(Muon_isTriggering, BToKMuMu_l1Idx)";
  std::string l2Trg = isEE ? "-1" : "Take(Muon_isTriggering, BToKMuMu_l2Idx)";
  std::string nTrg = isEE ? "RVec<unsigned int> (nBtriplet, nTriggerMuon)" : "(nTriggerMuon - B_l1_isTriggering - B_l2_isTriggering)"; // for ele as nTriggerMuon flattened over triplets
  std::string B_l1_pT = isEE ? "Take(Electron_pt, BToKEE_l1Idx)" : "Take(Muon_pt, BToKMuMu_l1Idx)";
  std::string B_l2_pT = isEE ? "Take(Electron_pt, BToKEE_l2Idx)" : "Take(Muon_pt, BToKMuMu_l2Idx)";
  std::string B_k_pT = isEE ? "Take(ProbeTracks_pt, BToKEE_kIdx)" : "Take(ProbeTracks_pt, BToKMuMu_kIdx)";
  std::string B_fit_mass = isEE ? "BToKEE_fit_mass" : "BToKMuMu_fit_mass";
  std::string B_cos2D = isEE ? "BToKEE_cos2D" : "BToKMuMu_cos2D";
  std::string B_vtxProb = isEE ? "BToKEE_svprob" : "BToKMuMu_svprob";
  std::string B_pT = isEE ? "BToKEE_fit_pt" : "BToKMuMu_fit_pt";
  std::string B_eta = isEE ? "BToKEE_fit_eta" : "BToKMuMu_fit_eta";
  std::string B_mll_llfit = isEE ? "BToKEE_mll_llfit" : "BToKMuMu_mll_llfit";
  std::string B_mll_fullfit = isEE ? "BToKEE_mll_fullfit" : "BToKMuMu_mll_fullfit";
  std::string B_l_xy_unc = isEE ? "BToKEE_l_xy_unc" : "BToKMuMu_l_xy_unc";
  std::string B_l_xyS = isEE ? "BToKEE_l_xy/BToKEE_l_xy_unc" : "BToKMuMu_l_xy/BToKMuMu_l_xy_unc";
  std::string B_l1_isGood = isEE ? "(RVec<bool>) (Take(Electron_isPFoverlap, BToKEE_l1Idx) == 0)" : "RVec<bool> (nBtriplet, true)"; // can add further requirements
  std::string B_l2_isGood = isEE ? "(RVec<bool>) (Take(Electron_isPFoverlap, BToKEE_l2Idx) == 0)" : "RVec<bool> (nBtriplet, true)"; // can add further requirements


  auto n = d.Define("lumi", "luminosityBlock")
    .Define("isEE", isElectronChannel.c_str())
    .Define("nBtriplet", nB.c_str())
    .Define("B_l1_isTriggering", l1Trg.c_str())
    .Define("B_l2_isTriggering", l2Trg.c_str())
    .Define("nExtraTrg", nTrg.c_str())
    .Define("B_l1_pT", B_l1_pT.c_str())
    .Define("B_l2_pT", B_l2_pT.c_str())
    .Define("B_k_pT", B_k_pT.c_str())
    .Define("B_fit_mass", B_fit_mass.c_str())
    .Define("B_cos2D", B_cos2D.c_str())
    .Define("B_vtxProb", B_vtxProb.c_str())
    .Define("B_pT", B_pT.c_str())
    .Define("B_eta", B_eta.c_str())
    .Define("B_mll_llfit", B_mll_llfit.c_str())
    .Define("B_mll_fullfit", B_mll_fullfit.c_str())
    .Define("B_l_xy_unc", B_l_xy_unc.c_str())
    .Define("B_l_xyS", B_l_xyS.c_str())
    .Define("B_l1_isGood", B_l1_isGood.c_str())
    .Define("B_l2_isGood", B_l2_isGood.c_str());


  //  n.Foreach(InvertValue, {"B_l1_isGood"});

  auto n2 = n.Define("good_B", IsGood, {"nBtriplet", "B_l1_pT", "B_l2_pT", "B_k_pT", "nExtraTrg", "B_cos2D", "B_vtxProb", "B_l_xyS", "B_l_xy_unc", "B_pT", "B_eta", 
	"B_l1_isGood", "B_l2_isGood"})
    .Define("rankVtx", Rank, {"nBtriplet", "good_B", "B_vtxProb"});
  

  auto selected = n2.Filter("Any(good_B == true)", "goodB");


  std::vector<string> listColumns = {"B_fit_mass", "good_B", "rankVtx", 
				     "B_l1_pT", "B_l2_pT", "B_k_pT", "B_l_xyS", "B_cos2D", "B_vtxProb", "B_pT", "B_eta",
				     "B_mll_fullfit", "B_mll_llfit", "B_l_xy_unc", 
				     "nBtriplet", "event", "B_l1_isGood", "B_l2_isGood"};

  //for(auto ij : listColumns) std::cout << ij << std::endl;

  if(isMC){

    // 443 = JPsi    521 = B+
    std::string B_l1_genParent = isEE ? "Take(Electron_genPartFlav, BToKEE_l1Idx)" : "Take(Muon_genPartFlav, BToKMuMu_l1Idx)";
    std::string B_l2_genParent = isEE ? "Take(Electron_genPartFlav, BToKEE_l2Idx)" : "Take(Muon_genPartFlav, BToKMuMu_l2Idx)";
    std::string B_k_genParent = isEE ? "Take(ProbeTracks_genPartFlav, BToKEE_kIdx)" : "Take(ProbeTracks_genPartFlav, BToKMuMu_kIdx)";
    
    auto mc = n2.Define("B_l1_genParent", B_l1_genParent.c_str())
      .Define("B_l2_genParent", B_l2_genParent.c_str())
      .Define("B_k_genParent", B_k_genParent.c_str());
    
    listColumns.push_back("B_l1_genParent");
    listColumns.push_back("B_l2_genParent");
    listColumns.push_back("B_k_genParent");

    mc.Snapshot("newtree", Form("newfile_isMC%d_isEE%d.root", isMC, isEE), listColumns);
  }

  else{


  /*add all variables needed for BDT*/
  /*maybe ?    "event"       "luminosityBlock"       "run"     */

    selected.Snapshot("newtree", Form("newfile_isMC%d_isEE%d.root", isMC, isEE), listColumns);
  }

  // your code goes here 
  t.Stop(); 
  t.Print();


  // if(isMC){
  //   //5279.26 +/- 0.17 
  //   auto s = n2.Filter("Any(B_fit_mass > )")
  //   n2.SnapShot("signal", "MC_signalTree.root")

  // }  
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

