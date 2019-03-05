import glob, os, sys, subprocess
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

MASH="/home-1/bsolomo9@jhu.edu/software/mash-Linux64-v2.1/mash"
inputFile=sys.argv[1]
outName=sys.argv[2]
#suffix=".fa.msh"

with open(inputFile) as myFile:
    inputMatrix=json.load(myFile)

Z = linkage(distMatrix)

plt.figure()
dendrogram(Z)
plt.show()
plt.savefig(outName)
