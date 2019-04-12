import seaborn as sns
import pandas as pd

import glob, os, sys, subprocess
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np 

from sklearn.datasets import load_iris, load_digits
from sklearn.model_selection import train_test_split


sketchFile=sys.argv[1] #/net/langmead-bigmem-ib.bluecrab.cluster/storage/bsolomo9/dash_registers_03_11_2019.txt
csv_file=sys.argv[2] #~/1kgenomes_sampleinfo_parsed.csv 

# Build tissue dictionary
count=0
popDic={}
with open(csv_file) as myFile:
    for line in myFile:
        line=line.strip()
        spline=line.split(',')
        sample=spline[1]
        pop=spline[3]
        popDic[sample]=pop

superPop_labels=["EAS","EUR","AFR","AMR","SAS"]
superPop=[]
dataArray=[]
with open(sketchFile) as myFile:
	for line in myFile:
		line = line.strip()
		if line:
			spline=line.split(" [")
			name=spline[0].split("/")[-1].split("_")[0]
			tempArray=spline[1].split(",")
			
			outArray=[]
			for i in tempArray:
				i=i.strip()
				pi = i.split("]")[0]
				outArray.append(int(pi))
		
			k=popDic[name]
	
			dataArray.append(outArray)
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
for i in range(len(dataArray[0])):
	dataArray_labels.append("reg{}".format(i))

print(len(dataArray))
print(dataArray_labels)
print(superPop)
print(superPop_labels)

#iris_df = pd.DataFrame(iris.data, columns=iris.feature_names)
#iris_df['species'] = pd.Series(iris.target).map(dict(zip(range(3),iris.target_names)))
#sns.pairplot(iris_df, hue='species');

#print(iris.data)
#print(iris.feature_names)
#print(iris.target)
#print(iris.target_names)

#iris_df = pd.DataFrame(iris.data, columns=iris.feature_names)
#iris_df['species'] = pd.Series(iris.target).map(dict(zip(range(3),iris.target_names)))
#sns.pairplot(iris_df, hue='species');
