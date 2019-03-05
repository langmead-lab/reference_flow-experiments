import glob, os, sys, subprocess
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

MASH="/home-1/bsolomo9@jhu.edu/software/mash-Linux64-v2.1/mash"
dataDir=sys.argv[1]
outName=sys.argv[2]
#suffix=".fa.msh"
suffix="-1M.msh"

mashList = glob.glob(dataDir+"*"+suffix)

#Build the 1D condensed data matrix (Upper triangle of distance matrix, in order 0,1 to 0,N, 1,2, ...
# Build an NxN empty array
distMatrix = []

# Build distance matrix
for i1,v1 in enumerate(mashList):
	for i2, v2 in enumerate(mashList):
		#pcall = "{} {} {}".format("echo", "hi", "blah")
		#proc = subprocess.Popen(["echo", "test"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if i2 < i1:
			proc = subprocess.Popen([MASH, "dist", v1, v2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			output, err = proc.communicate()
			spout = output.split("\t")
			dist = float(spout[2])
			distMatrix.append(dist)

print len(distMatrix)

Z = linkage(distMatrix)

plt.figure()
dendrogram(Z)
plt.show()
plt.savefig(outName)
