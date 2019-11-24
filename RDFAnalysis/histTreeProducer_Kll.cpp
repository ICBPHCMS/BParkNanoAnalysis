/*
----------------
g++ -Wall -o histTreeProducer_Kll `root-config --cflags --glibs ` histTreeProducer_Kll.cpp
----------------
./histTreeProducer_Kll --JOBid (-1,1,2...) --inList nanoAOD_list.txt --outFile /path/to/output_file.root --isMC (0,1) --isEE (0,1) --isResonant (0,1) --testFile /path/to/input_file.root
----------------
Add --tree 1 to store the output as tree
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
//#include <regex>

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
  int isTree = 0;
  std::string testFile = "-1";

  parseInputs(argc, argv, JOBid, inList, outFName, isMC, isEE, isResonant, testFile);

  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--tree") {
      if (i + 1 < argc) {
	isTree = atoi(argv[i+1]);
	break;
      }
      else {
	std::cerr << " --tree option requires one argument " << std::endl;
	return 1;
      }
    }
  }

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
    std::string base = "/vols/cms/vc1116/BParking/BDT_variables/";
	if(isEE) 
inputFileList.push_back(base + "tree_BParkingNANO_2019Oct25_MC_BuToKee/*_1.root");
	else 
inputFileList.push_back(base + "");
    }
    else{
    std::string base = "/vols/cms/vc1116/BParking/BDT_variables/";
    inputFileList.push_back(base + "tree_BParkingNANO_2019Oct21_Run2018A_part2/*_1.root");
    }
  }//defaut
 
  inList = "-1";
  for(auto ij : inputFileList) std::cout << " file = " << ij << std::endl;
 
  ROOT::RDataFrame d("newtree", inputFileList);
 
  
  auto n = d.Define("lumi_all", "lumi_All")
    .Define("eventToT_all", "eventToT_All")
    .Define("run_all", "run_All")
    .Define("KEE_fit_mass_sel", "Take(KEE_fit_mass_All, Idx_sel)")
    .Define("KEE_fit_massErr_sel", "Take(KEE_fit_massErr_All, Idx_sel)")
    .Define("KEE_cos2D_sel", "Take(KEE_cos2D_All, Idx_sel)")
    .Define("KEE_vtxProb_sel", "Take(KEE_vtxProb_All, Idx_sel)")
    .Define("KEE_l_xy_unc_sel", "Take(KEE_l_xy_unc_All, Idx_sel)")
    .Define("KEE_l_xyS_sel", "Take(KEE_l_xyS_All, Idx_sel)")
    .Define("KEE_pT_sel", "Take(KEE_pT_All, Idx_sel)")
    .Define("KEE_eta_sel", "Take(KEE_eta_All, Idx_sel)")
    .Define("KEE_phi_sel",  "Take(KEE_phi_All, Idx_sel)")
    .Define("KEE_mll_fullfit_sel", "Take(KEE_mll_fullfit_All, Idx_sel)")
    .Define("KEE_k_eta_sel", "Take(KEE_k_eta_All, Idx_sel)")
    .Define("KEE_k_phi_sel", "Take(KEE_k_phi_All, Idx_sel)")
    .Define("KEE_k_pt_sel",  "Take(KEE_k_pt_All, Idx_sel)")
    .Define("KEE_k_dz_sel", "Take(KEE_k_dz_All, Idx_sel)")
    .Define("KEE_k_dzE_sel", "Take(KEE_k_dzE_All, Idx_sel)")
    .Define("KEE_k_dxy_sel", "Take(KEE_k_dxy_All, Idx_sel)")
    .Define("KEE_k_dxyE_sel", "Take(KEE_k_dxyE_All, Idx_sel)")
    .Define("KEE_k_isPck_sel", "Take(KEE_k_isPck_All, Idx_sel)")
    .Define("e1_eta_sel", "Take(e1_eta_All, Idx_sel)")
    .Define("e2_eta_sel", "Take(e2_eta_All, Idx_sel)")
    .Define("e1_phi_sel", "Take(e1_phi_All, Idx_sel)")
    .Define("e2_phi_sel", "Take(e2_phi_All, Idx_sel)")
    .Define("e1_pt_sel", "Take(e1_pt_All, Idx_sel)")
    .Define("e2_pt_sel", "Take(e2_pt_All, Idx_sel)")
    .Define("e1_dz_sel", "Take(e1_dz_All, Idx_sel)")
    .Define("e2_dz_sel", "Take(e2_dz_All, Idx_sel)")
    .Define("e1_dzE_sel", "Take(e1_dzE_All, Idx_sel)")
    .Define("e2_dzE_sel", "Take(e2_dzE_All, Idx_sel)")
    .Define("e1_dxy_sel", "Take(e1_dxy_All, Idx_sel)")
    .Define("e2_dxy_sel", "Take(e2_dxy_All, Idx_sel)")
    .Define("e1_dxyE_sel", "Take(e1_dxyE_All, Idx_sel)")
    .Define("e2_dxyE_sel", "Take(e2_dxyE_All, Idx_sel)")
    .Define("e1_isConvVeto_sel", "Take(e1_isConvVeto_All, Idx_sel)")
    .Define("e2_isConvVeto_sel", "Take(e2_isConvVeto_All, Idx_sel)")
    .Define("e1_seedID_sel", "Take(e1_seedID_All, Idx_sel)")
    .Define("e2_seedID_sel", "Take(e2_seedID_All, Idx_sel)")
    .Define("e1_mvaID_sel", "Take(e1_mvaID_All, Idx_sel)")
    .Define("e2_mvaID_sel", "Take(e2_mvaID_All, Idx_sel)")
    .Define("e1_isPF_sel", "Take(e1_isPF_All, Idx_sel)")
    .Define("e2_isPF_sel", "Take(e2_isPF_All, Idx_sel)")
    .Define("e1_isLowPt_sel", "Take(e1_isLowPt_All, Idx_sel)")
    .Define("e2_isLowPt_sel", "Take(e2_isLowPt_All, Idx_sel)");
  
 
  if(!isTree){
    auto h_KEE_fit_mass = n.Histo1D( {"KEE_fit_mass", "", 75, 4.5, 6.0}, "KEE_fit_mass_sel");
    auto h_KEE_fit_massErr = n.Histo1D( {"KEE_fit_massErr", "", 50, 0., 1.}, "KEE_fit_massErr_sel");
    auto h_KEE_cos2D = n.Histo1D( {"KEE_cos2D", "", 60, -0.1, 1.1}, "KEE_cos2D_sel");
    auto h_KEE_vtxProb = n.Histo1D( {"KEE_vtxProb", "", 50, 0., 1.}, "KEE_vtxProb_sel");
    auto h_KEE_l_xy_unc = n.Histo1D( {"KEE_l_xy_unc", "", 30, 0., 0.6}, "KEE_l_xy_unc_sel");
    auto h_KEE_l_xyS = n.Histo1D( {"KEE_l_xyS", "", 700, 0., 70.}, "KEE_l_xyS_sel");
    auto h_KEE_pT = n.Histo1D( {"KEE_pT", "", 700, 0., 70.}, "KEE_pT_sel"); 
    auto h_KEE_eta = n.Histo1D( {"KEE_eta", "", 80, -4., 4.}, "KEE_eta_sel");   
    auto h_KEE_phi = n.Histo1D( {"KEE_phi", "", 80, -4., 4.}, "KEE_phi_sel");  
    auto h_KEE_mll_fullfit = n.Histo1D( {"KEE_mll_fullfit", "", 300, 0., 6.}, "KEE_mll_fullfit_sel");   
    auto h_KEE_k_eta = n.Histo1D( {"KEE_k_eta", "", 60, -3., 3.}, "KEE_k_eta_sel"); 
    auto h_KEE_k_phi = n.Histo1D( {"KEE_k_phi", "", 80, -4., 4.}, "KEE_k_phi_sel");   
    auto h_KEE_k_pt = n.Histo1D( {"KEE_k_pt", "", 250, 0., 25.}, "KEE_k_pt_sel");  
    auto h_KEE_k_dz = n.Histo1D( {"KEE_k_dz", "", 300, -15., 15.}, "KEE_k_dz_sel");   
    auto h_KEE_k_dzE = n.Histo1D( {"KEE_k_dzE", "", 400, -2000., 2000.}, "KEE_k_dzE_sel");  
    auto h_KEE_k_dxy = n.Histo1D( {"KEE_k_dxy", "", 20, -1., 1.}, "KEE_k_dxy_sel");  
    auto h_KEE_k_dxyE = n.Histo1D( {"KEE_k_dxyE", "", 800, -40., 40.}, "KEE_k_dxyE_sel");  
    auto h_KEE_k_isPck = n.Histo1D( {"KEE_k_isPck", "", 2, 0., 2.}, "KEE_k_isPck_sel");
    auto h_e1_eta = n.Histo1D( {"e1_eta", "", 60, -3., 3.}, "e1_eta_sel");  
    auto h_e2_eta = n.Histo1D( {"e2_eta", "", 60, -3., 3.}, "e2_eta_sel");  
    auto h_e1_phi = n.Histo1D( {"e1_phi", "", 80, -4., 4.}, "e1_phi_sel");  
    auto h_e2_phi = n.Histo1D( {"e2_phi", "", 80, -4., 4.}, "e2_phi_sel");  
    auto h_e1_pt = n.Histo1D( {"e1_pt", "", 400, 0., 40.}, "e1_pt_sel");
    auto h_e2_pt = n.Histo1D( {"e2_pt", "", 250, 0., 25.}, "e2_pt_sel");  
    auto h_e1_dz = n.Histo1D( {"e1_dz", "", 300, -15., 15.}, "e1_dz_sel");  
    auto h_e2_dz = n.Histo1D( {"e2_dz", "", 300, -15., 15.}, "e2_dz_sel");  
    auto h_e1_dzE = n.Histo1D( {"e1_dzE", "", 50, 0., 5.}, "e1_dzE_sel");  
    auto h_e2_dzE = n.Histo1D( {"e2_dzE", "", 40, 0., 4.}, "e2_dzE_sel");  
    auto h_e1_dxy = n.Histo1D( {"e1_dxy", "", 100, -5., 5.}, "e1_dxy_sel");   
    auto h_e2_dxy = n.Histo1D( {"e2_dxy", "", 180, -9., 9.}, "e2_dxy_sel");  
    auto h_e1_dxyE = n.Histo1D( {"e1_dxyE", "", 40, 0., 4.}, "e1_dxyE_sel");  
    auto h_e2_dxyE = n.Histo1D( {"e2_dxyE", "", 80, 0., 8.}, "e2_dxyE_sel");
    auto h_e1_isConvVeto = n.Histo1D( {"e1_isConvVeto", "", 2, 0., 2.}, "e1_isConvVeto_sel");  
    auto h_e2_isConvVeto = n.Histo1D( {"e2_isConvVeto", "", 2, 0., 2.}, "e2_isConvVeto_sel");  
    auto h_e1_seedID = n.Histo1D( {"e1_seedID", "", 150, 0., 15.}, "e1_seedID_sel");  
    auto h_e2_seedID = n.Histo1D( {"e2_seedID", "", 200, -5., 15.}, "e2_seedID_sel");  
    auto h_e1_mvaID = n.Histo1D( {"e1_mvaID", "", 200, -10., 10.}, "e1_mvaID_sel");  
    auto h_e2_mvaID = n.Histo1D( {"e2_mvaID", "", 200, -10., 10.}, "e2_mvaID_sel");
    auto h_e1_isPF = n.Histo1D( {"e1_isPF", "", 2, 0., 2.}, "e1_isPF_sel");  
    auto h_e2_isPF = n.Histo1D( {"e2_isPF", "", 2, 0., 2.}, "e2_isPF_sel");  
    auto h_e1_isLowPt = n.Histo1D( {"e1_isLowPt", "", 2, 0., 2.}, "e1_isLowPt_sel");  
    auto h_e2_isLowPt = n.Histo1D( {"e2_isLowPt", "", 2, 0., 2.}, "e2_isLowPt_sel");
  
    std::string outName = Form("histos_RDF_Kll_isMC%d_isEE%d.root", isMC, isEE);
    if(outFName != "-1") outName = outFName;  
  
    TFile outHistos(outName.c_str(), "recreate");

    outHistos.cd();

    h_KEE_fit_mass->Write();
    h_KEE_fit_massErr->Write();
    h_KEE_cos2D->Write();
    h_KEE_vtxProb->Write();
    h_KEE_l_xy_unc->Write();
    h_KEE_l_xyS->Write();  
    h_KEE_pT->Write();
    h_KEE_eta->Write();
    h_KEE_phi->Write();
    h_KEE_mll_fullfit->Write();
    h_KEE_k_eta->Write();  
    h_KEE_k_phi->Write();  
    h_KEE_k_pt->Write();
    h_KEE_k_dz->Write();  
    h_KEE_k_dzE->Write();
    h_KEE_k_dxy->Write();
    h_KEE_k_dxyE->Write();
    h_KEE_k_isPck->Write();
    h_e1_eta->Write();  
    h_e2_eta->Write();  
    h_e1_phi->Write();  
    h_e2_phi->Write();
    h_e1_pt->Write();  
    h_e2_pt->Write();
    h_e1_dz->Write();  
    h_e2_dz->Write();  
    h_e1_dzE->Write();  
    h_e2_dzE->Write();
    h_e1_dxy->Write();  
    h_e2_dxy->Write();  
    h_e1_dxyE->Write();  
    h_e2_dxyE->Write();
    h_e1_isConvVeto->Write();  
    h_e2_isConvVeto->Write();
    h_e1_seedID->Write();  
    h_e2_seedID->Write();
    h_e1_mvaID->Write();  
    h_e2_mvaID->Write();
    h_e1_isPF->Write();  
    h_e2_isPF->Write();
    h_e1_isLowPt->Write();   
    h_e2_isLowPt->Write();  
  
    outHistos.Close();

  }
  else{

    if(outFName == "-1")
      n.Snapshot("selTree", Form("tree_RDF_Kll_isMC%d_isEE%d.root", isMC, isEE), "\\b([^ ]*)(_sel)|\\b([^ ]*)(_all)");
    else
      n.Snapshot("selTree", outFName.c_str(), "\\b([^ ]*)(_sel)|\\b([^ ]*)(_all)");

  }
  
  t.Stop(); 
  t.Print();

}  
 
