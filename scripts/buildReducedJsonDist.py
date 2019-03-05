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
minFilter=0.2

with open(jsonPrefix+"_hapA.json") as myFile:
	matrixHapA=json.load(myFile)
with open(jsonPrefix+"_hapB.json") as myFile:
	matrixHapB=json.load(myFile)


outMatrixHapA=[]
outMatrixHapB=[]

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
		outMatrixHapA.append(rowA)
		outMatrixHapB.append(rowB)
	
	if total % 1000 == 0:
		print "Processed {} variants, kept {} variants.".format(str(total), totalVar)

	
onHapA = outName+"_hapA.json"
onHapB = outName+"_hapB.json"
with open(onHapA, 'w') as outFile:
    json.dump(outMatrixHapA, outFile)

with open(onHapB, 'w') as outFile:
    json.dump(outMatrixHapB, outFile)
