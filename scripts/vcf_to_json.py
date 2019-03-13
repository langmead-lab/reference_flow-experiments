import glob, os, sys, subprocess
#from scipy.cluster.hierarchy import dendrogram, linkage
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
import numpy as np
#import gzip
import json

vcf=sys.argv[1] #not gzipped!
outName=sys.argv[2]

# Build an NxN empty array
distMatrix = []
matrix=False
length = 0
head_l = 9


matrixHapA=[]
matrixHapB=[]
count = 0

with open(vcf) as myFile:
	for line in myFile:
		line = line.strip()
		if matrix==True:
			count+=1
			if count%1000 == 0:
				print "Processed " + str(count) + " lines."
			temp_hapA = []
			temp_hapB = []
			row = line.split("\t")
			for val in row[9:]:
				#tuple-ify
				sval = val.split("|")
				if len(sval)!=2:
					print "ERROR ON LINE: " + line
					break
				temp_hapA.append(bool(int(sval[0])))
				temp_hapB.append(bool(int(sval[1])))
				#temp.append( (sval[0], sval[1]) )
			#distMatrix.append(temp)
			matrixHapA.append(temp_hapA)
			matrixHapB.append(temp_hapB)
		if matrix==False and line[0]=="#" and line[1]!="#": #This should be the 'header' file which labels the rows and the like
			header = line.split("\t")
			length = len(header)-head_l
			matrix = True


onHapA = outName+"_hapA.json"
onHapB = outName+"_hapB.json"
with open(onHapA, 'w') as outFile:
	json.dump(matrixHapA, outFile)

with open(onHapB, 'w') as outFile:
	json.dump(matrixHapB, outFile)
	
'''
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

'''
