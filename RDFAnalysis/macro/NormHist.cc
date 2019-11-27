/*
A macro to normalize and save to a file all histograms from an input file
root NormHist.cc'("/vols/cms/vc1116/BParking/BDT_variables/SB_HIST_BParkingNANO_2019Oct21_Run2018A_part2.root")'
https://root.cern.ch/doc/master/loopdir11_8C.html
*/

void NormHist(std::string inFile) {
   
   TFile* inF = TFile::Open(inFile.c_str());
   
   std::string outFile = inFile;
   outFile = std::regex_replace(outFile, std::regex("\\.root"), "_NORM.root");
   
   TFile outF(outFile.c_str(), "recreate");
   outF.cd();
   
   for(auto k : *inF->GetListOfKeys()) {
      TKey *key = static_cast<TKey*>(k);
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      TH1 *h = key->ReadObject<TH1>();
      double_t scale = 1/h->Integral();
      h->Scale(scale);
      h->SetOption("HIST"); 
      h->Write();
   }
   
   outF.Close();

}
