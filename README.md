# BParkNanoAnalysis with RDataFrame



## Getting started (need ROOT version > 6.18)

```shell
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
git clone https://github.com/amartelli/BParkNanoAnalysis
source BParkNanoAnalysis/RDFAnalysis/settings.sh
```


## Run the analyzer
reads the input trees and dumps a reduced tree rejecting the events without good triplets plus
a column flagging the good triplet plus
a column with the rank of the triplet in the event - to plot 1 triplet per event - 

```shell
cd RDFAnalysis
root RDF_MC_KMuMu_test.C
```

## Next, to come
