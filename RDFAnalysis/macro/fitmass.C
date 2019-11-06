//plot and fit BToKll 

//example to run
//root fitBmass_fromHistos.C'(1, "/vols/cms/amartell/BParking/data_allStat/outMassHistos_Kee_tightSel_KeePrefitMass.root")'
//root fitBmass_fromHistos.C'(0, "/vols/cms/amartell/BParking/data_allStat/outMassHistos_Kmumu_tightSel.root")'


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


using namespace RooFit;


void fitmass(int isEleFinalState, std::string inFile, int isB, int isJPsiBin){

  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::cout << " isEleFinalState = " << isEleFinalState << std::endl;


  TFile* inF = TFile::Open(inFile.c_str());

  TH1F* h_mass;
  //  TH1F* h_JPsimass;

  if(isB && isJPsiBin == 0) h_mass = (TH1F*)inF->Get("B_mass")->Clone("h_mass");  
  else if(isB == 0) h_mass = (TH1F*)inF->Get("JPsi_mass")->Clone("h_mass");  
  else if(isB == 1 && isJPsiBin == 1) h_mass = (TH1F*)inF->Get("B_mass_JPsibin")->Clone("h_mass");  

  

  RooWorkspace w("w");    

  if(isB)  w.factory("x[4.5, 6.]");  
  else w.factory("x[2., 4.]");  

  // w.factory("nbackground[10000, 0, 100000]");   
  // w.factory("nbackgroundR[10, 0, 100]");   
  // w.factory("nsignal[100, 0.0, 10000]");

  w.factory("nBkg[100, 0, 5.e3]"); 
  w.factory("nSig[1.e5, 0.0, 1.e6]");      

  
  //  w.factory("Gaussian::B_smodel1(xB,muB[5.3,4.5,6],sigmaB[0.05,0,0.1])");
  //  w.factory("RooCBShape::B_smodel(xB,muB[5.3,4.5,6],sigmaB[0.05,0,2.],alphaB[0.01], nB[2])");


  //ele double CB
  //JPsi
  //parMean = 3.09558 parSigma = 0.039532
  //  [#1] INFO:Eval -- RooRealVar::setRange(xJPsi) new range named 'signalRange' created with bounds [2.97698,3.21417]
  if(!isB) w.factory("RooDoubleCBFast::smodel(x,mu[3.1,2.,4],sigma[0.04,0.001,0.08], aL[1.5, 0, 5.], nL[2, 0, 10], aR[1.5, 0, 5.], nR[10, 0, 15])");
  else w.factory("RooDoubleCBFast::smodel(x,mu[5.3,4.5,6],sigma[0.04,0.001,0.08], aL[1.5, 0, 5.], nL[2, 0, 10], aR[1.5, 0, 5.], nR[10, 0, 15])");


  //muon double CB
  //JPsi
  // parMean = 3.09641 parSigma = 0.0263931
  //  [#1] INFO:Eval -- RooRealVar::setRange(xJPsi) new range named 'signalRange' created with bounds [3.01723,3.17559]
  if(!isB) w.factory("RooDoubleCBFast::smodel(x,mu[3.1,2.,4],sigma[0.02,0.001,0.05], aL[1.5, 0, 3.], nL[2, 0, 5], aR[1.5, 0, 5.], nR[10, 5, 15])");
  else w.factory("RooDoubleCBFast::smodel(x,mu[5.3,4.5,6],sigma[0.02,0.001,0.05], aL[1.5, 0, 3.], nL[2, 0, 5], aR[1.5, 0, 5.], nR[10, 5, 15])");


  w.factory("Exponential::bmodel(x,tau[-2,-3,0])");
  //  w.factory("Exponential::JPsi_bmodel(xJPsi,tauJ[-0.5,-3,0])");

  //argus
  //  w.factory("ArgusBG::B_argus(mes,5.291,argpar[-20,-100,-1])");
  //w.factory("ArgusBG::JPsi_bmodel(xJPsi,muJ[3.,2.,4],cJ[5, 0, 100])");

  //  w.factory("FCONV::smodel(x,smodel1,smodel2)");

  // RooAbsPdf * B_smodel = w.pdf("B_smodel2");

  // RooAbsPdf * B_bmodel = w.pdf("B_bmodel1");
  // RooAbsPdf * JPsi_bmodel = w.pdf("JPsi_bmodel1");



  //  w.factory("SUM::model(nBkg * bmodel, nSig * smodel)");
  w.factory("SUM::model(nSig * smodel)");
  RooAbsPdf * model = w.pdf("model");
  // w.factory("SUM::JPsi_model(nBkg * JPsi_bmodel, nS * JPsi_smodel)");
  // RooAbsPdf * JPsi_model = w.pdf("JPsi_model");
    
  RooDataHist hMass("hMass", "hMass", *w.var("x"), Import(*(h_mass)));

  w.Print();

  //fit JPsi
  RooFitResult * rJ = model->fitTo(hMass, Minimizer("Minuit2"),Save(true));
  std::cout << "\n  >>>>>> rJ fit status = " << rJ->status() << std::endl;

  RooPlot * plotJ = w.var("x")->frame();
  if(isEleFinalState){
    if(isB && isJPsiBin) plotJ->SetXTitle("JPsibin - K(JPsi)ee mass (GeV)");
    else if(isB) plotJ->SetXTitle("K(JPsi)ee mass (GeV)");
    else plotJ->SetXTitle("(JPsi)ee mass (GeV)");
  }
  else{
    if(isB && isJPsiBin) plotJ->SetXTitle("JPsibin - K(JPsi)#mu#mu mass (GeV)");
    else if(isB) plotJ->SetXTitle("K(JPsi)#mu#mu mass (GeV)");
    else plotJ->SetXTitle("(JPsi)#mu#mu mass (GeV)");
  }
  plotJ->SetTitle("");
  if(isB) plotJ->SetAxisRange(4.,6);
  else plotJ->SetAxisRange(2.,4);
  hMass.plotOn(plotJ);
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
  //  RooRealVar* parB_J = (RooRealVar*) rJ->floatParsFinal().find("nBkg");
  auto nEv_postFitJ = parS_J->getValV();
  auto nEvError_postFitJ = parS_J->getError();
  // auto nBkg_postFitJ = parB_J->getValV();
  // auto nBkgError_postFitJ = parB_J->getError();

  std::cout << " **** JPsi => N events signal = \t " << parS_J->getValV() << " error = " << parS_J->getError() 
    /*<< " bkg events = " << parB_J->getValV() << " error = " << parB_J->getError() */
	    << std::endl;


  RooRealVar* parMeanJ = (RooRealVar*) rJ->floatParsFinal().find("mu");
  RooRealVar* parSigmaJ = (RooRealVar*) rJ->floatParsFinal().find("sigma");
  
  float meanValJ = parMeanJ->getValV();
  float sigmaValJ = parSigmaJ->getValV();

  std::cout << "\n  parMean = " << parMeanJ->getValV() << " parSigma = " << parSigmaJ->getValV() << std::endl;
    
  w.var("x")->setRange("signalRange", meanValJ - 3.*sigmaValJ, meanValJ + 3.*sigmaValJ);
  /*
  RooAbsReal* bkgIntegralJ = w.pdf("bmodel")->createIntegral(*(w.var("x")), NormSet(*(w.var("x"))), Range("signalRange")) ;
  std::cout << "\n  bkgIntegral = " << bkgIntegralJ->getVal() << std::endl;
  auto nBkgInt_postFitJ = bkgIntegralJ->getVal() * nBkg_postFitJ;
  auto nBkgIntError_postFitJ = nBkgError_postFitJ * bkgIntegralJ->getVal();
  */

    TCanvas * cc = new TCanvas();
    cc->SetLogy(0);
    plotJ->Draw();

    TLatex tL;
    tL.SetNDC();
    tL.SetTextSize(0.04);
    tL.SetTextFont(42);
    tL.DrawLatex(0.65,0.8, Form("chi2 = %.1f ",chi2_J));
    // TLatex tL2;
    // tL2.SetNDC();
    // tL2.SetTextSize(0.05);
    // tL2.SetTextFont(42);
    // tL2.DrawLatex(0.65,0.85, Form("B %.1f +/- %1.f",nBkgInt_postFitJ, nBkgIntError_postFitJ));


    std::string outName = "plots/Bmass_MC/";
    if(isB){
      if(isEleFinalState) outName += "Kee";
      else outName += "Kmm";
      if(isJPsiBin) outName += "_JPsibin";
    }
    else{
      if(isEleFinalState) outName += "ee";
      else outName += "mm";
    }


    cc->Print((outName+".png").c_str(), "png");
    gPad->SetLogy();
    cc->Print((outName+"_log.png").c_str(), "png");

    // if(isEleFinalState && isB){
    //   cc->Print(("plots/Bmass_MC/Kee_%s.png",h_mass->GetName()), "png");
    //   gPad->SetLogy();
    //   cc->Print(Form("plots/Bmass_MC/Kee_%s_log.png",h_mass->GetName()), "png");
    // }
    // else if(isEleFinalState && !isB){
    //   cc->Print(Form("plots/Bmass_MC/ee_%s.png",h_mass->GetName()), "png");
    //   gPad->SetLogy();
    //   cc->Print(Form("plots/Bmass_MC/ee_%s_log.png",h_mass->GetName()), "png");
    // }
    // else if(!isEleFinalState && !isB){
    //   cc->Print(Form("plots/Bmass_MC/mumu_%s.png",h_mass->GetName()), "png");
    //   gPad->SetLogy();
    //   cc->Print(Form("plots/Bmass_MC/mumu_%s_log.png",h_mass->GetName()), "png");
    // }
    // else if(!isEleFinalState && isB){
    //   cc->Print(Form("plots/Bmass_MC/Kmumu_%s.png",h_mass->GetName()), "png");
    //   gPad->SetLogy();
    //   cc->Print(Form("plots/Bmass_MC/Kmumu_%s_log.png",h_mass->GetName()), "png");
    // }

  
}  
