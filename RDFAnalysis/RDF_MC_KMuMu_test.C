#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>


//flag triplet passing the selection criteria
ROOT::VecOps::RVec<bool> IsGood(unsigned int nB, 
				ROOT::VecOps::RVec<float>& pT1, ROOT::VecOps::RVec<float>& pT2, ROOT::VecOps::RVec<float>& pTk, 
 				ROOT::VecOps::RVec<unsigned int>& nTrg, ROOT::VecOps::RVec<float>& cos2D, ROOT::VecOps::RVec<float>& vtxP, 
 				ROOT::VecOps::RVec<float>& disp, ROOT::VecOps::RVec<float>& dispU, ROOT::VecOps::RVec<float>& pT, ROOT::VecOps::RVec<float>& eta) { 

  ROOT::VecOps::RVec<bool> goodB(nB, false);
  for (auto ij=0; ij<nB; ++ij){
    if (pT1[ij] > 1. && pT2[ij] > 0.5 && pTk[ij] > 0.8 &&
     	nTrg[ij] > 0 && cos2D[ij] > 0.99 && vtxP[ij] > 0.1 && disp[ij] > 2 && dispU[ij] != 0 && pT[ij] > 3. && std::abs(eta[ij]) < 2.4)
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


void RDF_MC_KMuMu_test(){

  ROOT::RDataFrame d("Events", "/eos/cms//store/group/cmst3/group/bpark/BParkingNANO_2019Oct14/BuToKJpsi_ToMuMu_probefilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_ToMuMu/191014_075537/0000/BParkNANO_mc_2019Oct14_*.root");


  auto n = d.Define("B_l1_isTriggering", "Take(Muon_isTriggering, BToKMuMu_l1Idx)")
    .Define("B_l2_isTriggering", "Take(Muon_isTriggering, BToKMuMu_l2Idx)")
    .Define("nExtraTrg", "nTriggerMuon - B_l1_isTriggering - B_l2_isTriggering")
    .Define("B_l1_pT", "Take(Muon_pt, BToKMuMu_l1Idx)")
    .Define("B_l2_pT", "Take(Muon_pt, BToKMuMu_l2Idx)")
    .Define("B_k_pT", "Take(ProbeTracks_pt, BToKMuMu_kIdx)")
    .Define("B_fit_mass", "BToKMuMu_fit_mass")
    .Define("B_cos2D", "BToKMuMu_cos2D")
    .Define("B_vtxProb", "BToKMuMu_svprob")
    .Define("B_pT", "BToKMuMu_fit_pt")
    .Define("B_eta", "BToKMuMu_fit_eta")
    .Define("B_mll_llfit", "BToKMuMu_mll_llfit")
    .Define("B_mll_fullfit", "BToKMuMu_mll_fullfit")
    .Define("B_l_xy_unc", "BToKMuMu_l_xy_unc")
    .Define("B_l_xyS", "BToKMuMu_l_xy/BToKMuMu_l_xy_unc");

  auto n2 = n.Define("good_B", IsGood, {"nBToKMuMu", "B_l1_pT", "B_l2_pT", "B_k_pT", "nExtraTrg", "B_cos2D", "B_vtxProb", "B_l_xyS", "B_l_xy_unc", "B_pT", "B_eta"})
    .Define("rankVtx", Rank, {"nBToKMuMu", "good_B", "BToKMuMu_svprob"});


  auto g = n2.Filter("Any(good_B == true)", "goodB");

  g.Snapshot("newtree", "newfile.root", {"B_fit_mass", "good_B", "rankVtx", "B_l1_pT", "B_l2_pT", "B_k_pT", "B_l_xyS", "B_cos2D", "B_vtxProb", "B_pT", "B_eta", 
             "B_mll_fullfit", "B_mll_llfit", "B_l_xy_unc"});

  
}




// next steps:
// - add branch for MC matching
// configure inputs/ selections/ and type: MC vs DATA ee vs mumu
// save in parallel trees for BDT training 
// ...
