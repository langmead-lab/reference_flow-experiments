import glob, os, sys, subprocess
#from scipy.cluster.hierarchy import dendrogram, linkage
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
import numpy as np
#import gzip
import json

jsonPrefix=sys.argv[1] #gzipped!
outName=sys.argv[2]
minFilter=0.4

with open(jsonPrefix+"_hapA.json") as myFile:
	matrixHapA=json.load(myFile)
with open(jsonPrefix+"_hapB.json") as myFile:
	matrixHapB=json.load(myFile)
totalVar = 0.0 #kept
total= 0.0 #total
# Build an NxN empty array
distMatrix = []
numSamples = len(matrixHapA[0])
for i in range(numSamples):
	distMatrix.append([0.0]*numSamples)

# Compute similarity between NxN samples for each of M variants
# There's a lot of variants, consider subsampling vari (M) based on total expression
for vari in range(len(matrixHapA)):
	count=0.0
	total+=1
	rowA=matrixHapA[vari]
	rowB=matrixHapB[vari]
	for sample in rowA:
		if sample:
			count+=1
	for sample in rowB:
		if sample:
			count+=1

	if count / (2.0*numSamples) >= minFilter:
		totalVar+=1.0
		temp=[0]*numSamples
		for sami in range(numSamples):
			for samj in range(sami+1, numSamples):
				if matrixHapA[vari][sami]==matrixHapA[vari][samj]:
					distMatrix[sami][samj]+=1.0
					distMatrix[samj][sami]+=1.0
				if matrixHapB[vari][sami]==matrixHapB[vari][samj]:
					distMatrix[sami][samj]+=1.0
					distMatrix[samj][sami]+=1.0
		
#Normalize matrix to 0-1 DISTANCE instead of similarity
#Max similarity is 2*total variants
normalizedSim=np.array(distMatrix)
normalizedSim=normalizedSim / (totalVar*2)

#Set diagonals to sim=1
for i in range(numSamples):
	normalizedSim[i][i]=1.0

#Swap sim and dif
sim2diff=np.array(1.0)
normalizedDist=np.subtract(sim2diff,normalizedSim)

print "Total Variants Kept: " + str(totalVar)
print "Total Variants: " + str(total)
# Double output to be 'safe'
#print normalizedDist

with open(outName, 'w') as outFile:
    json.dump(normalizedDist, outFile)


