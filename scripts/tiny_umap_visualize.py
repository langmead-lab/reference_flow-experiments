import seaborn as sns
import pandas as pd
import json
import glob, os, sys, subprocess
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np 
import umap

from sklearn.datasets import load_iris, load_digits
from sklearn.model_selection import train_test_split


jsonPrefix=sys.argv[1] #/net/langmead-bigmem-ib.bluecrab.cluster/storage/bsolomo9/json_vcf/reduced_min02_chr21_hg38 (Leave off suffix)

# CSV has same order as json (json and csv both read line by line from vcf file)
csv_file=sys.argv[2] #~/1kgenomes_sampleinfo_parsed.csv 
outFig=sys.argv[3]


print("Loading jsons")
with open(jsonPrefix+"_hapA.json") as myFile:
    matrixHapA=json.load(myFile)
#with open(jsonPrefix+"_hapB.json") as myFile:
#    matrixHapB=json.load(myFile)

print("Loaded jsons")


# int conversion from bool / reduce matrix
tiny=100
tinyMatrix=[]
for i in range(tiny):
	temp=[]
	for si in range(len(matrixHapA[i])): #its all the same but why not
		if matrixHapA[i][si]:
			temp.append(1)
		else:
			temp.append(0)
	tinyMatrix.append(temp)
		
# change matrix so that samples have arrays of variants rather then variants have arrays of samples

corrTinyMatrix=np.swapaxes(tinyMatrix,0,1)

print(len(corrTinyMatrix))
print(len(corrTinyMatrix[0]))


# Build tissue dictionary
count=0
popDic={}
superPop_labels=["EAS","EUR","AFR","AMR","SAS"]
superPop=[]
print("Building superpop array")
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
print("Built superpop array")
			
print("Building data array lebels")
dataArray_labels=[]
for i in range(len(corrTinyMatrix[0])):
	dataArray_labels.append("var{}".format(i))
print("Built data array labels")

print("Reduce data")

reducer = umap.UMAP()
embedding = reducer.fit_transform(corrTinyMatrix)
embedding.shape

#print(len(matrixHapA))
#print(len(matrixHapB[0]))
print(len(dataArray_labels))
print(superPop)
print(superPop_labels)


print("pandas dataframing")
onek_df = pd.DataFrame(corrTinyMatrix, columns=dataArray_labels)
print("Adding data labels")
#onek_df['superpop'] = pd.Series(superPop).map(dict(zip(range(len(superPop_labels)),superPop_labels)))
print("Plotting")
um_plot=plt.scatter(embedding[:, 0], embedding[:, 1], c=[sns.color_palette()[x] for x in superPop])
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP projection of the DOUBLETINY dataset', fontsize=24);

#sns_plot = sns.pairplot(onek_df, hue='superpop')
print("Saving figure")
um_plot.savefig(outFig)


#print(iris.data)
#print(iris.feature_names)
#print(iris.target)
#print(iris.target_names)

#iris_df = pd.DataFrame(iris.data, columns=iris.feature_names)
#iris_df['species'] = pd.Series(iris.target).map(dict(zip(range(3),iris.target_names)))
#sns.pairplot(iris_df, hue='species');
