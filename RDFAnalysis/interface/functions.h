#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <string>
#include "TStopwatch.h"
#include <vector>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>



void parseInputs(int argc, char **argv,  
		 std::string& JOBid, std::string& inList, std::string& outF,
		 int& isMC, int& isEE, int& isResonant, 
		 std::string& testFile){

  for(int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--JOBid") {
      if(i + 1 < argc) {
        JOBid = argv[i+1];
        break;
      }
      else {
	std::cerr << " --JOBid option requires one argument " << std::endl;
        return;
      }
    }
  }
  for(int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--inList") {
      if (i + 1 < argc) {
        inList = argv[i+1];
        break;
      }
      else {
	std::cerr << " --inList option requires one argument " << std::endl;
        return;
      }
    }
  }
  for(int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--outFile") {
      if (i + 1 < argc) {
        outF = argv[i+1];
        break;
      }
      else {
	std::cerr << " --outFile option requires one argument " << std::endl;
        return;
      }
    }
  }
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isMC") {
      if (i + 1 < argc) {
        isMC = atoi(argv[i+1]);
        break;
      }
      else {
	std::cerr << " --isMC option requires one argument " << std::endl;
        return;
      }
    }
  }
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isEE") {
      if (i + 1 < argc) {
	isEE = atoi(argv[i+1]);
	break;
      }
      else {
	std::cerr << " --isEE option requires one argument " << std::endl;
	return;
      }
    }
  }
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isResonant") {
      if (i + 1 < argc) {
	isResonant = atoi(argv[i+1]);
	break;
      }
      else {
	std::cerr << " --isEE option requires one argument " << std::endl;
	return;
      }
    }
  }
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--testFile") {
      if (i + 1 < argc) {
	testFile = argv[i+1];
	break;
      } 
      else {
	std::cerr << " --testFile option requires one argument " << std::endl;
	return;
      }
    }
  }

  return;
}


int CountNperEvent(ROOT::VecOps::RVec<unsigned int>& goodIdxs){
  return int(goodIdxs.size());
}


//IsGood
ROOT::VecOps::RVec<int> cutBased(unsigned int nB,
				 ROOT::VecOps::RVec<float>& pT1, ROOT::VecOps::RVec<float>& pT2,
				 ROOT::VecOps::RVec<float>& pTk, ROOT::VecOps::RVec<unsigned int> nTrg,
				 ROOT::VecOps::RVec<float>& cos2D, ROOT::VecOps::RVec<float>& vtxP,
				 ROOT::VecOps::RVec<float>& disp, ROOT::VecOps::RVec<float>& dispU,
				 ROOT::VecOps::RVec<float>& pT, ROOT::VecOps::RVec<float>& eta,
				 ROOT::VecOps::RVec<unsigned int>& l1_PFoverlap, ROOT::VecOps::RVec<unsigned int>& l2_PFoverlap) {
  
  ROOT::VecOps::RVec<int> goodB(nB, 0);
  for(unsigned int ij=0; ij<nB; ++ij){
    if(pT1[ij] > 1. && pT2[ij] > 0.5 && pTk[ij] > 0.8 &&
       nTrg[ij] > 0 && cos2D[ij] > 0.99 && vtxP[ij] > 0.1 && disp[ij] > 2 && dispU[ij] != 0 &&
       pT[ij] > 3. && std::abs(eta[ij]) < 2.4 && l1_PFoverlap[ij] == 0 && l2_PFoverlap[ij] == 0)
      goodB[ij] = 1;
  }
  return goodB;
}


//KEE IsGood
ROOT::VecOps::RVec<int> KEEcut(unsigned int nB,
			       ROOT::VecOps::RVec<float>& pT1, ROOT::VecOps::RVec<float>& pT2,
			       ROOT::VecOps::RVec<float>& pTk,
			       ROOT::VecOps::RVec<float>& cos2D, ROOT::VecOps::RVec<float>& vtxP,
			       ROOT::VecOps::RVec<float>& disp, ROOT::VecOps::RVec<float>& dispU,
			       ROOT::VecOps::RVec<float>& pT, ROOT::VecOps::RVec<float>& eta,
			       ROOT::VecOps::RVec<unsigned int>& l1_PFoverlap, ROOT::VecOps::RVec<unsigned int>& l2_PFoverlap) {

  ROOT::VecOps::RVec<int> goodB(nB, 0);
  for(unsigned int ij=0; ij<nB; ++ij){
    if(pT1[ij] > 1. && pT2[ij] > 0.5 && pTk[ij] > 0.8 &&
       cos2D[ij] > 0.99 && vtxP[ij] > 0.1 && disp[ij] > 2 && dispU[ij] != 0 &&
       pT[ij] > 3. && std::abs(eta[ij]) < 2.4 && l1_PFoverlap[ij] == 0 && l2_PFoverlap[ij] == 0)
      goodB[ij] = 1;
  }
  return goodB;
}


//KMM IsGood
ROOT::VecOps::RVec<int> KMMcut(unsigned int nB,
			       ROOT::VecOps::RVec<float>& pT1, ROOT::VecOps::RVec<float>& pT2,
			       ROOT::VecOps::RVec<float>& pTk, ROOT::VecOps::RVec<unsigned int> nTrg,
			       ROOT::VecOps::RVec<float>& cos2D, ROOT::VecOps::RVec<float>& vtxP,
			       ROOT::VecOps::RVec<float>& disp, ROOT::VecOps::RVec<float>& dispU,
			       ROOT::VecOps::RVec<float>& pT, ROOT::VecOps::RVec<float>& eta) {

  ROOT::VecOps::RVec<int> goodB(nB, 0);
  for(unsigned int ij=0; ij<nB; ++ij){
    if(pT1[ij] > 1. && pT2[ij] > 0.5 && pTk[ij] > 0.8 &&
       nTrg[ij] > 0 && cos2D[ij] > 0.99 && vtxP[ij] > 0.1 && disp[ij] > 2 && dispU[ij] != 0 &&
       pT[ij] > 3. && std::abs(eta[ij]) < 2.4)
      goodB[ij] = 1;
  }
  return goodB;
}


//FilterGood
ROOT::VecOps::RVec<unsigned int> selectGoodIdx(ROOT::VecOps::RVec<int>& goodIdxs){

  ROOT::VecOps::RVec<unsigned int> goodB;
  for(unsigned int ij=0; ij<goodIdxs.size(); ++ij){
    if(goodIdxs[ij] == 1)
      goodB.push_back(ij);
  }
  return goodB;
}


ROOT::VecOps::RVec<unsigned int> selectGoodIdx(ROOT::VecOps::RVec<unsigned int>& goodIdxs){

  ROOT::VecOps::RVec<unsigned int> goodB;
  for(unsigned int ij=0; ij<goodIdxs.size(); ++ij){
    if(goodIdxs[ij] == 1)
      goodB.push_back(ij);
  }
  return goodB;
}




ROOT::VecOps::RVec<unsigned int> selectMllIdx (ROOT::VecOps::RVec<float>& ll_mass, 
					       ROOT::VecOps::RVec<unsigned int>& isEle){
  
  ROOT::VecOps::RVec<unsigned int> goodB;
  unsigned int totN = ll_mass.size();
  float min = -99;
  float max = -99;
  if(totN > 0 && isEle[0] == 0){
    min =  3.0172; //3.0964 - 3.* 0.0264;                                     
    max = 3.1756; //3.0964 + 3.* 0.0264;                                      
  }
  else if(totN > 0 && isEle[0] == 1){
    min =  2.9771; //3.0956 - 3.* 0.0395;                                     
    max = 3.2141; //3.0956 + 3.* 0.0395;                                      
  }
  for(unsigned int ij=0; ij<totN; ++ij){
    if(min < ll_mass[ij] && ll_mass[ij] < max)
      goodB.push_back(ij);
  }
  return goodB;
}
  

ROOT::VecOps::RVec<unsigned int> flagGenMatchExt(int& isResonant,
					ROOT::VecOps::RVec<int>& l1_pdgId, ROOT::VecOps::RVec<int>& l2_pdgId,
					ROOT::VecOps::RVec<int>& k_pdgId, ROOT::VecOps::RVec<int>& Ml1_pdgId,
					ROOT::VecOps::RVec<int>& Ml2_pdgId, ROOT::VecOps::RVec<int>& Mk_pdgId,
					ROOT::VecOps::RVec<int>& GMl1_pdgId, ROOT::VecOps::RVec<int>& GMl2_pdgId,
					ROOT::VecOps::RVec<int>& GMk_pdgId){
  auto totN = l1_pdgId.size();
  ROOT::VecOps::RVec<unsigned int> matched(totN, 0);
  for(unsigned int ij=0; ij<totN; ++ij){
    // std::cout << " l1_pdgId[ij] = " << l1_pdgId[ij] << " l2_pdgId[ij] " << l2_pdgId[ij] << " k_pdgId[ij]  = " << k_pdgId[ij]
    //        << " Ml1_pdgId[ij] = " << Ml1_pdgId[ij] << " Ml2_pdgId[ij] " << Ml2_pdgId[ij] << " Mk_pdgId[ij]  = " << Mk_pdgId[ij]
    //        << " GMl1_pdgId[ij] = " << GMl1_pdgId[ij] << " GMl2_pdgId[ij] " << GMl2_pdgId[ij] << std::endl;    

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
}


ROOT::VecOps::RVec<int> getPdg (ROOT::VecOps::RVec<int>& idx, ROOT::VecOps::RVec<int>& pdgL){
  int totN = idx.size();
  ROOT::VecOps::RVec<int> pdgs(totN, -1);
  for(auto ij=0; ij<totN; ++ij)
    if(idx[ij] != -1) pdgs[ij] = pdgL[idx[ij]];
  return pdgs;
}



ROOT::VecOps::RVec<float> computedR(ROOT::VecOps::RVec<unsigned int>& matchedIdxs,
				    ROOT::VecOps::RVec<float>& gen_eta1, ROOT::VecOps::RVec<float>& gen_eta2, ROOT::VecOps::RVec<float>& gen_eta3, 
				    ROOT::VecOps::RVec<float>& reco_eta1, ROOT::VecOps::RVec<float>& reco_eta2, ROOT::VecOps::RVec<float>& reco_eta3,
				    ROOT::VecOps::RVec<float>& gen_phi1, ROOT::VecOps::RVec<float>& gen_phi2, ROOT::VecOps::RVec<float>& gen_phi3, 
				    ROOT::VecOps::RVec<float>& reco_phi1, ROOT::VecOps::RVec<float>& reco_phi2, ROOT::VecOps::RVec<float>& reco_phi3){


  auto idx = ROOT::VecOps::Nonzero(matchedIdxs);

  auto genEta1 = Take(gen_eta1, idx);
  auto genEta2 = Take(gen_eta2, idx);
  auto genEta3 = Take(gen_eta3, idx); 

  auto genPhi1 = Take(gen_phi1, idx);
  auto genPhi2 = Take(gen_phi2, idx); 
  auto genPhi3 = Take(gen_phi3, idx);

  auto recoEta1 = Take(reco_eta1, idx);
  auto recoEta2 = Take(reco_eta2, idx);
  auto recoEta3 = Take(reco_eta3, idx);

  auto recoPhi1 = Take(reco_phi1, idx);
  auto recoPhi2 = Take(reco_phi2, idx);
  auto recoPhi3 = Take(reco_phi3, idx);

  ROOT::VecOps::RVec<float> dR_1 = ROOT::VecOps::DeltaR(genEta1, recoEta1, genPhi1, recoPhi1);
  ROOT::VecOps::RVec<float> dR_2 = ROOT::VecOps::DeltaR(genEta2, recoEta2, genPhi2, recoPhi2);
  ROOT::VecOps::RVec<float> dR_3 = ROOT::VecOps::DeltaR(genEta3, recoEta3, genPhi3, recoPhi3);


  ROOT::VecOps::RVec<float> sum(gen_eta1.size(), -1);
  int selected = 0;
  for(unsigned int ij=0; ij<matchedIdxs.size(); ++ij){
    if(matchedIdxs[ij] == 1){
      sum[ij] = dR_1[selected] + dR_2[selected] + dR_3[selected];
      ++selected;
    }
    else sum[ij] = -1;
  }
  return sum;
}



ROOT::VecOps::RVec<int> flagRank(unsigned int nB,
				 ROOT::VecOps::RVec<unsigned int>& goodB, 
				 ROOT::VecOps::RVec<float>& vtxP){

  auto sortIndices = Argsort(vtxP);
  ROOT::VecOps::RVec<int> rank;
  rank.resize(nB);
  int nRank = 0;
  for(unsigned int ij=0; ij<nB; ++ij){

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


//bigger to smaller
ROOT::VecOps::RVec<int> flagReverseRank(unsigned int nB,
					ROOT::VecOps::RVec<int>& goodB, 
					ROOT::VecOps::RVec<float>& vtxP){

  auto sortIndices = Argsort(vtxP);
  auto rev_sortIndices = Reverse(sortIndices);

  ROOT::VecOps::RVec<int> rank;
  rank.resize(nB);
  int nRank = 0;
  for(unsigned int ij=0; ij<nB; ++ij){

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


//Rankv2
ROOT::VecOps::RVec<int> flagRank_forselected(ROOT::VecOps::RVec<float>& vtxP){

  auto sortIndices = Argsort(vtxP);
  ROOT::VecOps::RVec<int> rank;
  auto totN = vtxP.size();
  rank.resize(totN);
  int nRank = 0;
  for(unsigned int ij=0; ij<totN; ++ij){
    rank[sortIndices[ij]] = nRank;
    ++nRank;
  }
  return rank;
}


//auto v1_intersect_v2 = Intersect(v1, v2);

/*
ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>
  apply_CutBased(int isMC, ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>& dt_All, std::vector<std::string>& listC){
  //get vector with indices of good triplets
  // do not use lepton ID => to be updated                                                                                                                                          
  auto n2v2 = dt_All.Define("idx_goodB", selectGoodIdx, {"cutBase_goodB_All"});
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
    .Define("rankVtx", "Take(rankVtx_All, idx_goodB)")                                                                   
    .Define("nBtriplet", "(unsigned int) idx_goodB.size()")
    .Define("eventToT", "Take(eventToT_All, idx_goodB)")
    .Define("weights", "Take(weights_All, idx_goodB)");
  listC = {"B_fit_mass", "idx_goodB", "rankVtx",
           "B_l1_pT", "B_l2_pT", "B_k_pT", "B_l_xyS",
           "B_cos2D", "B_vtxProb", "B_pT", "B_eta",
           "B_mll_fullfit", "B_mll_llfit", "B_l_xy_unc",
           "nBtriplet", "eventToT", "B_l1_isPF", "B_l2_isPF",
           "B_l1_isPFoverlap", "B_l2_isPFoverlap"};
  return filteresCandidates;
}
*/
