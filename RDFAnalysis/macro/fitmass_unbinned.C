//plot and fit BToKll 

//example to run
//root fitBmass_unbinned.C'(isEle, "file.root")'



#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TChain.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooArgusBG.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "TText.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

#include "../src/GBRMath.cc"
#include "../src/RooDoubleCBFast.cc"

#include <ROOT/RVec.hxx>

using namespace RooFit;
using namespace ROOT::VecOps;

void fitmass_unbinned(int isMC, int isEleFinalState, std::string inFile, 
		      int isB, int isJPsiBin, int isAll, int is2PF, int is2LowPt, int only1perEvt){

  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::cout << " isMC = " << isMC << " isEle = " << isEleFinalState << " isB (or ll) " << isB
	    << " isJPsiBin = " << isJPsiBin << " isAll = " << isAll 
	    << " is2PF = " << is2PF << " is2LowPt = " << is2LowPt << " only1perEvt = " << only1perEvt << std::endl;

  float minMll =  3.018;
  float maxMll = 3.174; 
  if(isEleFinalState){
    minMll =  2.972;
    maxMll = 3.224; 
  }

  RooWorkspace w("w");    
  RooRealVar x("x", "", 4.5, 6.);
  if(!isB) x.setRange(2., 4.);
  RooDataSet data ("data", "data", RooArgSet(x)); //, Cut("isGenMatched_All ==1"));
  

  TFile* inF = TFile::Open(inFile.c_str());

  
  TTree* tree = (TTree*)inF->Get("newtree");
  
  std::vector<float>* B_fit_mass_All = 0;
  std::vector<float>* B_mll_llfit_All = 0;
  std::vector<unsigned int>* isGenMatched_All = 0;
  std::vector<unsigned int>* B_l1_isPF_All = 0;
  std::vector<unsigned int>* B_l2_isPF_All = 0;
  std::vector<unsigned int>* B_l1_isPFoverlap_All = 0;
  std::vector<unsigned int>* B_l2_isPFoverlap_All = 0;
  std::vector<unsigned int>* cutBase_goodB_All = 0;
  std::vector<float>* dRwithGen_All = 0;
  std::vector<float>* rankVtx_All = 0;

  tree->SetBranchStatus("*", 0);

  if(isB)  tree->SetBranchStatus("B_fit_mass_All", 1);        tree->SetBranchAddress("B_fit_mass_All", &B_fit_mass_All);
  if(!isB || isJPsiBin) tree->SetBranchStatus("B_mll_llfit_All", 1);       tree->SetBranchAddress("B_mll_llfit_All", &B_mll_llfit_All);

  if(isAll || is2PF || is2LowPt){
    tree->SetBranchStatus("B_l1_isPF_All", 1);        tree->SetBranchAddress("B_l1_isPF_All", &B_l1_isPF_All);
    tree->SetBranchStatus("B_l2_isPF_All", 1);        tree->SetBranchAddress("B_l2_isPF_All", &B_l2_isPF_All);
    tree->SetBranchStatus("B_l1_isPFoverlap_All", 1); tree->SetBranchAddress("B_l1_isPFoverlap_All", &B_l1_isPFoverlap_All);
    tree->SetBranchStatus("B_l2_isPFoverlap_All", 1); tree->SetBranchAddress("B_l2_isPFoverlap_All", &B_l2_isPFoverlap_All);
  }
  if(isMC){
    tree->SetBranchStatus("dRwithGen_All", 1);             tree->SetBranchAddress("dRwithGen_All", &dRwithGen_All);
    tree->SetBranchStatus("isGenMatched_All", 1);          tree->SetBranchAddress("isGenMatched_All", &isGenMatched_All);
  }
  if(!isMC){
    tree->SetBranchStatus("cutBase_goodB_All", 1);           tree->SetBranchAddress("cutBase_goodB_All", &cutBase_goodB_All);
    tree->SetBranchStatus("rankVtx_All", 1);                 tree->SetBranchAddress("rankVtx_All", &rankVtx_All);
  }


  auto totN = tree->GetEntriesFast();
  RVec<float> rankV; 
  RVec<float> xVec; 

  int counter = 0;
  int counterH = 0;
  for(auto iEvt=0; iEvt<totN; ++iEvt){
    tree->GetEntry(iEvt);
    
    for(int ij=0; ij<B_fit_mass_All->size(); ++ij){

      if(isMC && isGenMatched_All->at(ij) != 1) continue;
      if(!isMC && cutBase_goodB_All->at(ij) != 1) continue;
      if(isJPsiBin && (B_mll_llfit_All->at(ij) < minMll || B_mll_llfit_All->at(ij) > maxMll)) continue;
      if(B_l1_isPFoverlap_All->at(ij) == 1 || B_l2_isPFoverlap_All->at(ij) == 1) continue;

      if(isAll){
	xVec.push_back(ij);
	if(isMC) rankV.push_back(dRwithGen_All->at(ij));
	else rankV.push_back(rankVtx_All->at(ij));
      }
      if(is2PF && B_l1_isPF_All->at(ij) == 1 && B_l2_isPF_All->at(ij) == 1){
	xVec.push_back(ij);
	if(isMC) rankV.push_back(dRwithGen_All->at(ij));
	else rankV.push_back(rankVtx_All->at(ij));
      }
      if(is2LowPt && B_l1_isPF_All->at(ij) == 0 && B_l2_isPF_All->at(ij) == 0){
	xVec.push_back(ij);
	if(isMC) rankV.push_back(dRwithGen_All->at(ij));
	else rankV.push_back(rankVtx_All->at(ij));
      }
    }//loop over branch

    if(rankV.size() > 0){
   
      auto idx = ArgMin(rankV);

      //find the best dR matched
      if(only1perEvt){
	x = isB ? B_fit_mass_All->at(xVec[idx]) : B_mll_llfit_All->at(xVec[idx]);
	data.add(RooArgSet(x));
	++counter;
	//if(B_fit_mass_All->at(xVec[idx]) > 5.6) ++counterH;
      }
      else{
	for(auto ip : xVec){
	  x = isB ? B_fit_mass_All->at(ip) : B_mll_llfit_All->at(ip);
	  data.add(RooArgSet(x));
	  ++counter;
	  //if(B_fit_mass_All->at(ip) > 5.6) ++counterH;
	}
      }
    }
    xVec.clear();
    rankV.clear();
  }
  std::cout << " counter = " << counter << " counterH = " << counterH << std::endl;  


  w.import(x);
  w.factory("nBkg[100, 0, 5.e3]"); 
  w.factory("nSig[1.e5, 0.0, 1.e6]");      
  
  //ele double CB
  if(isEleFinalState){
    if(isB)
w.factory("RooDoubleCBFast::smodel(x,mu[5.3,4.5,6],sigma[0.04,0.001,0.08], aL[1.5, 0, 5.], nL[2, 0, 10], aR[1.5, 0, 5.], nR[10, 0, 15])"); 
    else
w.factory("RooDoubleCBFast::smodel(x,mu[3.1,2.,4],sigma[0.04,0.02,0.08], aL[1.5, 0, 5.], nL[2, 0, 10], aR[1.5, 0, 5.], nR[10, 0, 15])");
  }

  //muon double CB
  if(!isEleFinalState){ 
    if(isB)
w.factory("RooDoubleCBFast::smodel(x,mu[5.3,4.5,6],sigma[0.02,0.001,0.05], aL[1.5, 0, 3.], nL[2, 0, 5], aR[1.5, 0, 5.], nR[10, 5, 15])");
    else
w.factory("RooDoubleCBFast::smodel(x,mu[3.1,2.,4],sigma[0.02,0.001,0.05], aL[1.5, 0, 3.], nL[2, 0, 5], aR[1.5, 0, 5.], nR[10, 5, 15])");
  }

  w.factory("Exponential::bmodel(x,tau[-2,-3,0])");

  w.factory("SUM::model(nSig * smodel)");
  RooAbsPdf * model = w.pdf("model");
    
  w.Print();



  //fit JPsi
  RooFitResult * rJ = model->fitTo(data, Minimizer("Minuit2"),Save(true));
  std::cout << "\n  >>>>>> rJ fit status = " << rJ->status() << std::endl;

  RooPlot * plotJ = w.var("x")->frame();
  if(isEleFinalState){
    if(isB) plotJ->SetXTitle("K(JPsi)ee mass (GeV)");
    else plotJ->SetXTitle("(JPsi)ee mass (GeV)");
  }
  else{
    if(isB) plotJ->SetXTitle("K(JPsi)#mu#mu mass (GeV)");
    else plotJ->SetXTitle("(JPsi)#mu#mu mass (GeV)");
  }
  plotJ->SetTitle("");
  if(isB) plotJ->SetAxisRange(4.5,6);
  else plotJ->SetAxisRange(2.,4);
  data.plotOn(plotJ, RooFit::Binning(200));
  model->plotOn(plotJ);
  //  model->plotOn(plotJ, Components("bmodel"),LineStyle(kDashed), LineColor(kBlue));
  model->plotOn(plotJ, Components("smodel"),LineColor(kRed+1));
  model->paramOn(plotJ, RooFit::Layout(0.2,0.4,0.9),RooFit::Format("NEA",AutoPrecision(1)));
  plotJ->getAttLine()->SetLineColorAlpha(kWhite, 0.2);
  plotJ->getAttText()->SetTextSize(0.03);
  plotJ->getAttText()->SetTextFont(42);

  float chi2_J = plotJ->chiSquare(6);

  std::cout << " chi2_J  = " << chi2_J << std::endl;
  
  RooRealVar* parS_J = (RooRealVar*) rJ->floatParsFinal().find("nSig");

  auto nEv_postFitJ = parS_J->getValV();
  auto nEvError_postFitJ = parS_J->getError();

  std::cout << " **** JPsi => N events signal = \t " << parS_J->getValV() << " error = " << parS_J->getError() 
	    << std::endl;


  RooRealVar* parMeanJ = (RooRealVar*) rJ->floatParsFinal().find("mu");
  RooRealVar* parSigmaJ = (RooRealVar*) rJ->floatParsFinal().find("sigma");
  
  float meanValJ = parMeanJ->getValV();
  float sigmaValJ = parSigmaJ->getValV();

  std::cout << "\n  parMean = " << parMeanJ->getValV() << " parSigma = " << parSigmaJ->getValV() << std::endl;
    
  w.var("x")->setRange("signalRange", meanValJ - 3.*sigmaValJ, meanValJ + 3.*sigmaValJ);


    TCanvas * cc = new TCanvas();
    cc->SetLogy(0);
    plotJ->Draw();
    TLatex tL;
    tL.SetNDC();
    tL.SetTextSize(0.04);
    tL.SetTextFont(42);
    tL.DrawLatex(0.65,0.8, Form("chi2 = %.1f ",chi2_J));
    TLatex tL2;
    tL2.SetNDC();
    tL2.SetTextSize(0.05);
    tL2.SetTextFont(42);
    std::string category = "All";
    if(is2PF) category = "PF-PF";
    if(is2LowPt) category = "LowPt-LowPt";
    tL2.DrawLatex(0.65,0.7, category.c_str());
    TLatex tL3;
    tL3.SetNDC();
    tL3.SetTextSize(0.05);
    tL3.SetTextFont(42);
    std::string type = "q^{2} = full";
    if(isJPsiBin) type = "q^{2} = JPsi";
    tL3.DrawLatex(0.65,0.6, type.c_str());
    TLatex tL4;
    tL4.SetNDC();
    tL4.SetTextSize(0.05);
    tL4.SetTextFont(42);
    std::string nGen = "all triplets";
    if(only1perEvt) nGen = "best triplet";
    tL4.DrawLatex(0.65,0.5, nGen.c_str());


    std::string outName = "plots/Bmass_MC/";
    if(!isMC) outName = "plots/Bmass_DATA/";
    if(isB){
      if(isEleFinalState) outName += "Kee";
      else outName += "Kmm";
      if(isJPsiBin) outName += "_JPsibin";
    }
    else{
      if(isEleFinalState) outName += "ee";
      else outName += "mm";
    }
    if(isAll){
      outName += "_All";
    }
    if(is2PF){
      outName += "_2PF";
    }
    if(is2LowPt){
      outName += "_2LowPt";
    }
    if(only1perEvt) outName += "_bestdRgenMatched";

    cc->Print((outName+".png").c_str(), "png");
    gPad->SetLogy();
    cc->Print((outName+"_log.png").c_str(), "png");

 
}  



  /*  
  RVec<float> valW;
  ROOT::RDataFrame d("newtree", inFile.c_str());
  auto d2 = d.Define("isGenMatched_idx",  [](const RVec<unsigned int>& matchedIdxs){return Nonzero(matchedIdxs);}, {"isGenMatched_All"});
  auto dd2 = d2.Define("B_fit_mass_genMatched", "Take(B_fit_mass_All, isGenMatched_idx)");
  //  auto loadData = [&valW](RVec<float>& val){return Concatenate(valW, val); };
  auto loadData = [&valW](RVec<float>& val){std::cout << "val.size() = " << val.size() << std::endl; };
  void ddd2 = dd2.Foreach(loadData, {"B_fit_mass_genMatched"});
  
  for(auto ij:valW){
    x = ij;
    data_BgenMatched.add(RooArgSet(x));
  }
  */
