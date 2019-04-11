import seaborn as sns
import pandas as pd
import json
import glob, os, sys, subprocess
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np 

from sklearn.datasets import load_iris, load_digits
from sklearn.model_selection import train_test_split


jsonPrefix=sys.argv[1] #/net/langmead-bigmem-ib.bluecrab.cluster/storage/bsolomo9/json_vcf/reduced_min02_chr21_hg38 (Leave off suffix)

# CSV has same order as json (json and csv both read line by line from vcf file)
csv_file=sys.argv[2] #~/1kgenomes_sampleinfo_parsed.csv 
outFig=sys.argv[3]


with open(jsonPrefix+"_hapA.json") as myFile:
    matrixHapA=json.load(myFile)
with open(jsonPrefix+"_hapB.json") as myFile:
    matrixHapB=json.load(myFile)

# Build tissue dictionary
count=0
popDic={}
superPop_labels=["EAS","EUR","AFR","AMR","SAS"]
superPop=[]
with open(csv_file) as myFile:
	for line in myFile:
		line=line.strip()
		spline=line.split(',')
		sample=spline[1]
		k=spline[3]
		
		if k=="CHB" or k=="JPT" or k=="CHS" or k=="CDX" or k=="KHV":
			superPop.append(0)
		elif k=="CEU" or k=="TSI" or k=="FIN" or k=="GBR" or k=="IBS":
			superPop.append(1)
		elif k=="YRI" or k=="LWK" or k=="GWD" or k=="MSL" or k=="ESN" or k=="ASW" or k=="ACB":
			superPop.append(2)
		elif k=="MXL" or k=="PUR" or k=="CLM" or k=="PEL":
			superPop.append(3)
		elif k=="GIH" or k=="PJL" or k=="BEB" or k=="STU" or k=="ITU":
			superPop.append(4)

			
dataArray_labels=[]
for i in range(len(matrixHapA[0])):
	dataArray_labels.append("var{}".format(i))

print(len(matrixHapA))
print(len(matrixHapB[0]))
print(len(dataArray_labels))
print(superPop)
print(superPop_labels)

onek_df = pd.DataFrame(matrixHapA, columns=dataArray_labels)
onek_df['superpop'] = pd.Series(superPop).map(dict(zip(range(len(superPop_labels)),superPop_labels)))
sns_plot = sns.pairplot(onek_df, hue='species')
sns_plot.savefig(outFig)


#print(iris.data)
#print(iris.feature_names)
#print(iris.target)
#print(iris.target_names)

#iris_df = pd.DataFrame(iris.data, columns=iris.feature_names)
#iris_df['species'] = pd.Series(iris.target).map(dict(zip(range(3),iris.target_names)))
#sns.pairplot(iris_df, hue='species');
