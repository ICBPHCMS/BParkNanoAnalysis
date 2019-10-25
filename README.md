# BParkNanoAnalysis with RDataFrame



## Getting started 

```shell
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
git cms-init
git clone https://github.com/amartelli/BParkNanoAnalysis
```


## Run the analyzer
need ROOT version > 6.18
reads the input trees and dumps a reduced tree 
rejecting all the bad triplets plus
adding a column with the rank of the triplet in the event - to plot 1 triplet per event - 
Compiled and non macros
```shell
cd RDFAnalysis
source BParkNanoAnalysis/RDFAnalysis/settings.sh
root RDF_MC_KMuMu_test.C'((int)isMC, (int)isEE)'
g++ -Wall -o RDF_MC_Kll_test `root-config --cflags --glibs ` RDF_MC_Kll_test.cpp 
./RDF_MC_Kll_test (int)isMC, (int)isEE 
```

## Next, to come
