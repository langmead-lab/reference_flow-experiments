import re, sys
import os, glob


csv_file = sys.argv[1]
cluster_file = sys.argv[2]
numClust = 10


class clustStats:
	def __init__(self):
		self.popCount={"CHB":0, "JPT":0, "CHS":0, "CDX":0,"KHV":0,"CEU":0,"TSI":0,"FIN":0,"GBR":0,"IBS":0,"YRI":0,"LWK":0,"GWD":0,"MSL":0,"ESN":0,"ASW":0,"ACB":0,"MXL":0,"PUR":0,"CLM":0,"PEL":0,"GIH":0,"PJL":0,"BEB":0,"STU":0,"ITU":0}

	def getPop(self):
		print self.popCount

	def getSuperPop(self):
		spDic={"EAS":0,"EUR":0,"AFR":0,"AMR":0,"SAS":0}
		for k in self.popCount.keys():
			if k=="CHB" or k=="JPT" or k=="CHS" or k=="CDX" or k=="KHV":
				spDic["EAS"]+=self.popCount[k]
			elif k=="CEU" or k=="TSI" or k=="FIN" or k=="GBR" or k=="IBS":
				spDic["EUR"]+=self.popCount[k]
			elif k=="YRI" or k=="LWK" or k=="GWD" or k=="MSL" or k=="ESN" or k=="ASW" or k=="ACB":
				spDic["AFR"]+=self.popCount[k]
			elif k=="MXL" or k=="PUR" or k=="CLM" or k=="PEL":
				spDic["AMR"]+=self.popCount[k]
			elif k=="GIH" or k=="PJL" or k=="BEB" or k=="STU" or k=="ITU":
				spDic["SAS"]+=self.popCount[k]
		
		print spDic	

count=0
popDic={}
with open(csv_file) as myFile:
	for line in myFile:
		line=line.strip()
		spline=line.split(',')
		sample=spline[1]
		pop=spline[3]
		popDic[sample]=pop


countArray=[]
for i in range(numClust):
	countArray.append(clustStats())

with open(cluster_file) as myFile:
	for line in myFile:
		line=line.strip()
		if line:
			spline = line.split(": ")
			sample = spline[0]
			clust = int(spline[1])-1
			pop=popDic[sample]
			#print pop
			#print clust
			#print pop
			#print countArray[clust]
			countArray[clust].popCount[pop]+=1

for i in countArray:
	i.getPop()

for i in countArray:
	i.getSuperPop()
