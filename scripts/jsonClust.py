import glob, os, sys, subprocess
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

MASH="/home-1/bsolomo9@jhu.edu/software/mash-Linux64-v2.1/mash"
inputFile=sys.argv[1]
outName=sys.argv[2]
#suffix=".fa.msh"

#Import NxN square matrix
#with open(inputFile) as myFile:
#    inputMatrix=json.load(myFile)

inputMatrix = [ [0, 0.1, 0.3],
				[0.1, 0, 0.2],
				[0.3, 0.2, 0] ]

inputLabels= ["A", "B", "C"]

maxDist=0.2
k= 2

#Build the 1D condensed data matrix (Upper triangle of distance matrix, in order 0,1 to 0,N, 1,2, ...
distMatrix = []

# Build distance matrix
for i1,v1 in enumerate(inputMatrix):
    for i2, v2 in enumerate(v1):
        #pcall = "{} {} {}".format("echo", "hi", "blah")
        #proc = subprocess.Popen(["echo", "test"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if i2 > i1:
			distMatrix.append(v2)

Z = linkage(distMatrix, method='ward')
print Z

clusters = fcluster(Z, maxDist, criterion='distance')
#clusters = fcluster(Z,k, criterion='maxclust') 
print clusters


plt.figure()
dendrogram(
    Z,
    truncate_mode='lastp',  # show only the last p merged clusters
    p=20,  # show only the last p merged clusters
    leaf_rotation=90.,
    leaf_font_size=12.,
    show_contracted=True,  # to get a distribution impression in truncated branches
)

#dendrogram(Z)
plt.show()
plt.savefig(outName)
