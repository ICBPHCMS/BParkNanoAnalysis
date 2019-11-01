//g++ -Wall -o RDF_MC_Kll_test `root-config --cflags --glibs ` RDF_MC_Kll_test.cpp 

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <string>
#include "TStopwatch.h"
#include <vector>


//flag triplet passing the selection criteria

ROOT::VecOps::RVec<int> IsGood(unsigned int nB,
			       ROOT::VecOps::RVec<float>& pT1, ROOT::VecOps::RVec<float>& pT2, 
			       ROOT::VecOps::RVec<float>& pTk, ROOT::VecOps::RVec<unsigned int> nTrg, 
			       ROOT::VecOps::RVec<float>& cos2D, ROOT::VecOps::RVec<float>& vtxP,
			       ROOT::VecOps::RVec<float>& disp, ROOT::VecOps::RVec<float>& dispU, 
			       ROOT::VecOps::RVec<float>& pT, ROOT::VecOps::RVec<float>& eta, 
			       ROOT::VecOps::RVec<unsigned int>& l1_PFoverlap, ROOT::VecOps::RVec<unsigned int>& l2_PFoverlap) {
  
  ROOT::VecOps::RVec<int> goodB(nB, 0);
  for (auto ij=0; ij<nB; ++ij){
    if(pT1[ij] > 1. && pT2[ij] > 0.5 && pTk[ij] > 0.8 &&
       nTrg[ij] > 0 && cos2D[ij] > 0.99 && vtxP[ij] > 0.1 && disp[ij] > 2 && dispU[ij] != 0 && 
       pT[ij] > 3. && std::abs(eta[ij]) < 2.4 && l1_PFoverlap[ij] == 0 && l2_PFoverlap[ij] == 0)
      goodB[ij] = 1;
  }

  return goodB;
}


//save indices of good candidate - to be used later for filtering 
ROOT::VecOps::RVec<unsigned int> FilterGood(ROOT::VecOps::RVec<int>& good_triples){

  ROOT::VecOps::RVec<unsigned int> goodB;
  for (auto ij=0; ij<good_triples.size(); ++ij){
    if(good_triples[ij] == 1) 
      goodB.push_back(ij);
  }

  return goodB;
}


//save for each good triplet the rank - choice is now wrt vtxCL - 
//useful to plot just 1 triplet per event 
ROOT::VecOps::RVec<int> Rank(unsigned int nB, 
			     ROOT::VecOps::RVec<int>& goodB, ROOT::VecOps::RVec<float>& vtxP){

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



//save for each good triplet the rank - choice is now wrt vtxCL - 
//useful to plot just 1 triplet per event 
ROOT::VecOps::RVec<int> Rankv2(ROOT::VecOps::RVec<float>& vtxP){

  auto sortIndices = Argsort(vtxP);
  ROOT::VecOps::RVec<int> rank;
  auto totN = vtxP.size();
  rank.resize(totN);
  int nRank = 0;
  for (auto ij=0; ij<totN; ++ij){
    rank[sortIndices[ij]] = nRank;
    ++nRank;
  }
  return rank; 
}


ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> 
apply_CutBased(int isMC, ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>& dt_All, std::vector<std::string>& listC){
  //get vector with indices of good triplets
  // do not use lepton ID => to be updated
  auto n2v2 = dt_All.Define("idx_goodB", FilterGood, {"good_B_All"});
  
  //filter branches: only save good triplets
  auto filteresCandidates = n2v2.Define("B_fit_mass", "Take(B_fit_mass_All, idx_goodB)")
    .Define("B_l1_pT", "Take(B_l1_pT_All, idx_goodB)")
    .Define("B_l2_pT", "Take(B_l2_pT_All, idx_goodB)")
    .Define("B_k_pT", "Take(B_k_pT_All, idx_goodB)")
    .Define("B_l_xyS", "Take(B_l_xyS_All, idx_goodB)")
    .Define("B_l_xy_unc", "Take(B_l_xy_unc_All, idx_goodB)")
    .Define("B_cos2D", "Take(B_cos2D_All, idx_goodB)")
    .Define("B_vtxProb", "Take(B_vtxProb_All, idx_goodB)")
    .Define("B_pT", "Take(B_pT_All, idx_goodB)")
    .Define("B_eta", "Take(B_eta_All, idx_goodB)")
    .Define("B_mll_fullfit", "Take(B_mll_fullfit_All, idx_goodB)")
    .Define("B_mll_llfit", "Take(B_mll_llfit_All, idx_goodB)")
    .Define("B_l1_isPF", "Take(B_l1_isPF_All, idx_goodB)")
    .Define("B_l2_isPF", "Take(B_l2_isPF_All, idx_goodB)")
    .Define("B_l1_isPFoverlap", "Take(B_l1_isPFoverlap_All, idx_goodB)")
    .Define("B_l2_isPFoverlap", "Take(B_l2_isPFoverlap_All, idx_goodB)")
    .Define("rankVtx", Rankv2, {"B_vtxProb"})  // should be as .Define("rankVtx", "Take(rankVtx_All, idx_goodB)")
    .Define("nBtriplet", "(unsigned int) idx_goodB.size()")
    .Define("eventToT", "Take(eventToT_All, idx_goodB)");
  
  listC = {"B_fit_mass", "idx_goodB", "rankVtx", 
	   "B_l1_pT", "B_l2_pT", "B_k_pT", "B_l_xyS", 
	   "B_cos2D", "B_vtxProb", "B_pT", "B_eta",
	   "B_mll_fullfit", "B_mll_llfit", "B_l_xy_unc", 
	   "nBtriplet", "eventToT", "B_l1_isPF", "B_l2_isPF", 
	   "B_l1_isPFoverlap", "B_l2_isPFoverlap"};
  
  return filteresCandidates;
};


//using namespace ROOT::VecOps;

int main(int argc, char **argv){

  TStopwatch t;
  t.Start();

  int isMC = atoi(argv[1]);
  int isEE = atoi(argv[2]);
  int doCutBased = atoi(argv[3]);
  int isResonant = atoi(argv[4]);

  std::string inputFileList = "/eos/cms//store/group/cmst3/group/bpark/BParkingNANO_2019Oct14/";
  if(isMC && !isEE) inputFileList += "BuToKJpsi_ToMuMu_probefilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_ToMuMu/191014_075537/0000/BParkNANO_mc_2019Oct14_*.root";
  else if(isMC && isEE) inputFileList += "BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee/191014_075238/0000/BParkNANO_mc_2019Oct14_*.root";


  std::cout << " isMC = " << isMC << "isEE = " << isEE << std::endl;
  std::cout << "inputFileList = " << inputFileList << std::endl;


  std::cout << " doCutBased = " << doCutBased << " isResonant = " << isResonant << std::endl;

  ROOT::RDataFrame d("Events", inputFileList.c_str());

  std::string nB = isEE ? "nBToKEE" : "nBToKMuMu";
  std::string isElectronChannel = isEE ? "ROOT::VecOps::RVec<int> (nBtriplet_All, 1)" : "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)";
  std::string event = "ROOT::VecOps::RVec<int> (nBtriplet_All, event)";
  std::string l1Trg = isEE ? "ROOT::VecOps::RVec<int> (nBtriplet_All, -1)" : "Take(Muon_isTriggering, BToKMuMu_l1Idx)";
  std::string l2Trg = isEE ? "ROOT::VecOps::RVec<int> (nBtriplet_All, -1)" : "Take(Muon_isTriggering, BToKMuMu_l2Idx)";
  // for ele as nTriggerMuon flattened over triplets
  std::string nTrg = isEE ? "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, nTriggerMuon)" : "(nTriggerMuon - B_l1_isTriggering_All - B_l2_isTriggering_All)"; 
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
  // can add further requirements on lepton ID
  std::string B_l1_isPF = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l1Idx)" : 
                                 "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l1Idx)"; 
  std::string B_l2_isPF = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l2Idx)" : 
                                 "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l2Idx)";
  std::string B_l1_isPFoverlap = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l1Idx)" : 
                                        "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)"; 
  std::string B_l2_isPFoverlap = isEE ? "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l2Idx)" : 
                                        "ROOT::VecOps::RVec<unsigned int> (nBtriplet_All, 0)"; 


  auto n = d.Define("lumi", "luminosityBlock")
    .Define("nBtriplet_All", nB.c_str())
    .Define("isEE_All", isElectronChannel.c_str())
    .Define("eventToT_All", event.c_str())
    .Define("B_l1_isTriggering_All", l1Trg.c_str())
    .Define("B_l2_isTriggering_All", l2Trg.c_str())
    .Define("nExtraTrg_All", nTrg.c_str())
    .Define("B_l1_pT_All", B_l1_pT.c_str())
    .Define("B_l2_pT_All", B_l2_pT.c_str())
    .Define("B_k_pT_All", B_k_pT.c_str())
    .Define("B_fit_mass_All", B_fit_mass.c_str())
    .Define("B_cos2D_All", B_cos2D.c_str())
    .Define("B_vtxProb_All", B_vtxProb.c_str())
    .Define("B_pT_All", B_pT.c_str())
    .Define("B_eta_All", B_eta.c_str())
    .Define("B_mll_llfit_All", B_mll_llfit.c_str())
    .Define("B_mll_fullfit_All", B_mll_fullfit.c_str())
    .Define("B_l_xy_unc_All", B_l_xy_unc.c_str())
    .Define("B_l_xyS_All", B_l_xyS.c_str())
    .Define("B_l1_isPF_All", B_l1_isPF.c_str())
    .Define("B_l2_isPF_All", B_l2_isPF.c_str())
    .Define("B_l1_isPFoverlap_All", B_l1_isPFoverlap.c_str())
    .Define("B_l2_isPFoverlap_All", B_l2_isPFoverlap.c_str());



  /////////////////////////////
  //flag the best triplet matched to gen
  //  ROOT::VecOps::RVec<int> findGenMatch = [&isResonant](ROOT::VecOps::RVec<int>& l1_pdgId, ROOT::VecOps::RVec<int>& l2_pdgId, 
  auto findGenMatch = [&isResonant](ROOT::VecOps::RVec<int>& l1_pdgId, ROOT::VecOps::RVec<int>& l2_pdgId, 
				    ROOT::VecOps::RVec<int>& k_pdgId, ROOT::VecOps::RVec<int>& Ml1_pdgId, 
				    ROOT::VecOps::RVec<int>& Ml2_pdgId, ROOT::VecOps::RVec<int>& Mk_pdgId, 
				    ROOT::VecOps::RVec<int>& GMl1_pdgId, ROOT::VecOps::RVec<int>& GMl2_pdgId, 
				    ROOT::VecOps::RVec<int>& GMk_pdgId){

    auto totN = l1_pdgId.size();
    ROOT::VecOps::RVec<int> matched(totN, 0);
    
    for (auto ij=0; ij<totN; ++ij){

    // std::cout << " l1_pdgId[ij] = " << l1_pdgId[ij] << " l2_pdgId[ij] " << l2_pdgId[ij] << " k_pdgId[ij]  = " << k_pdgId[ij] 
    // 	      << " Ml1_pdgId[ij] = " << Ml1_pdgId[ij] << " Ml2_pdgId[ij] " << Ml2_pdgId[ij] << " Mk_pdgId[ij]  = " << Mk_pdgId[ij]       
    // 	      << " GMl1_pdgId[ij] = " << GMl1_pdgId[ij] << " GMl2_pdgId[ij] " << GMl2_pdgId[ij] << std::endl;

    if(l1_pdgId[ij] == -1 || l2_pdgId[ij] == -1 || k_pdgId[ij] == -1) continue;
    if(l1_pdgId[ij] != -1. * l2_pdgId[ij] || std::abs(k_pdgId[ij]) != 321) continue;

    if(isResonant){
      if( Ml1_pdgId[ij] == Ml2_pdgId[ij] && Ml1_pdgId[ij] == 443 && 
	  GMl2_pdgId[ij] == GMl1_pdgId[ij] && GMl1_pdgId[ij] == Mk_pdgId[ij] && 
	  std::abs(GMl1_pdgId[ij]) == 521) 
	matched[ij] = 1;
    }
    else if (!isResonant){
      if(Ml1_pdgId[ij] == Ml2_pdgId[ij] && Ml2_pdgId[ij] == Mk_pdgId[ij] && 
	 std::abs(Ml1_pdgId[ij]) == 521) 
	matched[ij] = 1;
    }
    }
    return matched;
  };

  //define a branch flagging the good candidates and filter events with 0
  auto tree_All = n.Define("good_B_All", IsGood, {"nBtriplet_All", "B_l1_pT_All", "B_l2_pT_All", "B_k_pT_All", "nExtraTrg_All", 
	"B_cos2D_All", "B_vtxProb_All", "B_l_xyS_All", "B_l_xy_unc_All", "B_pT_All", "B_eta_All", 
	"B_l1_isPFoverlap_All", "B_l2_isPFoverlap_All"})
    .Define("rankVtx_All", Rank, {"nBtriplet_All", "good_B_All", "B_vtxProb_All"});
  // //auto selected = n2.Filter("Any(good_B == true)", "goodB");



  std::vector<std::string> listColumns_All = {"B_fit_mass_All", "good_B_All", "rankVtx_All", 
					      "B_l1_pT_All", "B_l2_pT_All", "B_k_pT_All", "B_l_xyS_All", 
					      "B_cos2D_All", "B_vtxProb_All", "B_pT_All", "B_eta_All",
					      "B_mll_fullfit_All", "B_mll_llfit_All", "B_l_xy_unc_All", 
					      "nBtriplet_All", "eventToT_All", "B_l1_isPF_All", "B_l2_isPF_All", 
					      "B_l1_isPFoverlap_All", "B_l2_isPFoverlap_All"};






  if(isMC){
    // listColumns_All.push_back("B_l1_genParent_All");
    // listColumns_All.push_back("B_l2_genParent_All");
    // listColumns_All.push_back("B_k_genParent_All");
    // listColumns_All.push_back("GenPart_l1_idx_All");    
    // listColumns_All.push_back("GenPart_l2_idx_All");
    // listColumns_All.push_back("GenPart_k_idx_All");
    // listColumns_All.push_back("GenPart_l1_pdgId_All");    
    // listColumns_All.push_back("GenPart_l2_pdgId_All");
    // listColumns_All.push_back("GenPart_k_pdgId_All");
    // listColumns_All.push_back("GenMothPart_l1_idx_All");
    // listColumns_All.push_back("GenMothPart_l2_idx_All");
    // listColumns_All.push_back("GenMothPart_k_idx_All");
    // listColumns_All.push_back("GenMothPart_l1_pdgId_All");
    // listColumns_All.push_back("GenMothPart_l2_pdgId_All");
    // listColumns_All.push_back("GenMothPart_k_pdgId_All");
    // listColumns_All.push_back("GenGMothPart_l1_idx_All");
    // listColumns_All.push_back("GenGMothPart_l2_idx_All");
    // listColumns_All.push_back("GenGMothPart_k_idx_All");
    // listColumns_All.push_back("GenGMothPart_l1_pdgId_All");
    // listColumns_All.push_back("GenGMothPart_l2_pdgId_All");
    // listColumns_All.push_back("GenGMothPart_k_pdgId_All");
    listColumns_All.push_back("isGenMatched_All");
    
    
    // 443 = JPsi    521 = B+  321 = K 
    std::string B_l1_genParent = isEE ? "Take(Electron_genPartFlav, BToKEE_l1Idx)" : "Take(Muon_genPartFlav, BToKMuMu_l1Idx)";
    std::string B_l2_genParent = isEE ? "Take(Electron_genPartFlav, BToKEE_l2Idx)" : "Take(Muon_genPartFlav, BToKMuMu_l2Idx)";
    std::string B_k_genParent = isEE ? "Take(ProbeTracks_genPartFlav, BToKEE_kIdx)" : "Take(ProbeTracks_genPartFlav, BToKMuMu_kIdx)";
    std::string GenPart_l1_idx = isEE ? "Take(Electron_genPartIdx, BToKEE_l1Idx)" : "Take(Muon_genPartIdx, BToKMuMu_l1Idx)";
    std::string GenPart_l2_idx = isEE ? "Take(Electron_genPartIdx, BToKEE_l2Idx)" : "Take(Muon_genPartIdx, BToKMuMu_l2Idx)";
    std::string GenPart_k_idx = isEE ? "Take(ProbeTracks_genPartIdx, BToKEE_kIdx)" : "Take(ProbeTracks_genPartIdx, BToKMuMu_kIdx)";
    
    auto getPdg = [] (ROOT::VecOps::RVec<int>& idx, ROOT::VecOps::RVec<int>& pdgL){
      int totN = idx.size();
      ROOT::VecOps::RVec<int> pdgs(totN, -1);
      for(auto ij=0; ij<totN; ++ij)
	if(idx[ij] != -1) pdgs[ij] = pdgL[idx[ij]];
      return pdgs;
    };

    auto mc = tree_All.Define("B_l1_genParent_All", B_l1_genParent.c_str())
      .Define("B_l2_genParent_All", B_l2_genParent.c_str())
      .Define("B_k_genParent_All", B_k_genParent.c_str())
      .Define("GenPart_l1_idx_All", GenPart_l1_idx.c_str())
      .Define("GenPart_l2_idx_All", GenPart_l2_idx.c_str())
      .Define("GenPart_k_idx_All", GenPart_k_idx.c_str())
      .Define("GenPart_l1_pdgId_All", getPdg, {"GenPart_l1_idx_All", "GenPart_pdgId"})
      .Define("GenPart_l2_pdgId_All", getPdg, {"GenPart_l2_idx_All", "GenPart_pdgId"})
      .Define("GenPart_k_pdgId_All", getPdg, {"GenPart_k_idx_All", "GenPart_pdgId"})
      .Define("GenMothPart_l1_idx_All", getPdg, {"GenPart_l1_idx_All", "GenPart_genPartIdxMother"})
      .Define("GenMothPart_l2_idx_All", getPdg, {"GenPart_l2_idx_All", "GenPart_genPartIdxMother"})
      .Define("GenMothPart_k_idx_All", getPdg, {"GenPart_k_idx_All", "GenPart_genPartIdxMother"})
      .Define("GenMothPart_l1_pdgId_All", getPdg, {"GenMothPart_l1_idx_All", "GenPart_pdgId"})
      .Define("GenMothPart_l2_pdgId_All", getPdg, {"GenMothPart_l2_idx_All", "GenPart_pdgId"})
      .Define("GenMothPart_k_pdgId_All", getPdg, {"GenMothPart_k_idx_All", "GenPart_pdgId"})
      .Define("GenGMothPart_l1_idx_All", getPdg, {"GenMothPart_l1_idx_All", "GenPart_genPartIdxMother"})
      .Define("GenGMothPart_l2_idx_All", getPdg, {"GenMothPart_l2_idx_All", "GenPart_genPartIdxMother"})
      .Define("GenGMothPart_k_idx_All", getPdg, {"GenMothPart_k_idx_All", "GenPart_genPartIdxMother"})
      .Define("GenGMothPart_l1_pdgId_All", getPdg, {"GenGMothPart_l1_idx_All", "GenPart_pdgId"}) 
      .Define("GenGMothPart_l2_pdgId_All", getPdg, {"GenGMothPart_l2_idx_All", "GenPart_pdgId"}) 
      .Define("GenGMothPart_k_pdgId_All", getPdg, {"GenGMothPart_k_idx_All", "GenPart_pdgId"});

    //    filteredTriplets
    /*
      .Define("B_l1_genParent", "Take(B_l1_genParent_All, idx_goodB)")
      .Define("B_l2_genParent", "Take(B_l2_genParent_All, idx_goodB)")
      .Define("B_k_genParent", "Take(B_k_genParent_All, idx_goodB)");
    */

    std::cout << " prima di mcGenMatched " << std::endl;
    auto mcGenMatched = mc.Define("isGenMatched_All", findGenMatch, {"GenPart_l1_pdgId_All", "GenPart_l2_pdgId_All", 
     	  "GenPart_k_pdgId_All", "GenMothPart_l1_pdgId_All", "GenMothPart_l2_pdgId_All", "GenMothPart_k_pdgId_All", 
	  "GenGMothPart_l1_pdgId_All", "GenGMothPart_l2_pdgId_All", "GenGMothPart_k_pdgId_All"});


    std::cout << " now saving " << std::endl;

    if(doCutBased){
      std::vector<std::string> listColumns;
      auto filtered_mc = apply_CutBased(isMC, mcGenMatched, listColumns);

      auto filtered_mc_Gen = filtered_mc.Define("B_l1_genParent", "Take(B_l1_genParent_All, idx_goodB)")
	.Define("B_l2_genParent", "Take(B_l2_genParent_All, idx_goodB)")
	.Define("B_k_genParent", "Take(B_k_genParent_All, idx_goodB)")
	.Define("isGenMatched", "Take(isGenMatched_All, idx_goodB)");

      filtered_mc_Gen.Snapshot("newtree", Form("newfile_isMC%d_isEE%d_CB.root", isMC, isEE), listColumns);
      //MC cut based con info sul gen matched
    }
    else  mcGenMatched.Snapshot("newtree", Form("newfile_isMC%d_isEE%d_All.root", isMC, isEE), listColumns_All);
    //all triplets with Mc gen matched info
  }//isMC
  else{

    if(doCutBased){
      std::vector<std::string> listColumns;
      auto filtered = apply_CutBased(isMC, tree_All, listColumns);
      filtered.Snapshot("newtree", Form("newfile_isMC%d_isEE%d_CB.root", isMC, isEE), listColumns);
    }
    else tree_All.Snapshot("newtree", Form("newfile_isMC%d_isEE%d_All.root", isMC, isEE), listColumns_All);
  }
  
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

