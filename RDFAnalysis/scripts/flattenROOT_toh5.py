#!/usr/bin/python                                                                                                                                          

import ROOT
import sys
from root_numpy import array2tree

import numpy as np
import pandas
#these 2 would allow to write directly in root... 
#but problem with version installed
#import root_pandas
#from root_pandas import readwrite

fileInputName = sys.argv[1]
#fileInputName = "../newfile_isMC1_isEE0.root"
print ("fileInputName = ", fileInputName)


df = ROOT.ROOT.RDataFrame("newtree", fileInputName)
npy = df.AsNumpy()

pd = pandas.DataFrame(npy)

print ("pre flat = ", pd)

#https://gist.github.com/rainsunny/a4fd8f760194dc55da67b0fd2384395e
def expand(df, expand_column):
    """
    Expand df on expand_column
    df: pandas.DataFrame
    expand_column: the column name to expand
    """
    lens = [ len(item) for item in  df[expand_column]]
    d = {}
    d[expand_column] = np.concatenate(df[expand_column].values)
    for col in df.columns.values:
        if col != expand_column:
            d[col] = np.repeat(df[col].values, lens)
    return pandas.DataFrame(d)


flat = expand(pd, "eventToT")

#this works in case (by hand)
#lens = [len(item) for item in pd['eventToT']]
#print("expected events = ", lens)
#flat = pandas.DataFrame( {"eventToT" : np.repeat(pd['eventToT'].values,lens), 
#                          "B_fit_mass" : np.concatenate(pd['B_fit_mass'].values)})


print ("post flat = ", flat)

#outFileName = (fileInputName+"_test.h5")
flat.to_hdf("testFile.h5", key='df', mode='w')

#end_np = flat.to_numpy()
#end_np = flat.rename_axis().values
#end_np = flat.values()
#tree = array2tree(end_np)
#tree.Print()
