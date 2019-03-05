import glob, os, sys, subprocess
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

# Build an NxN empty array
heatMap = []
for i in range(len(mashList)):
	tempList = [-1] * len(mashList)
	heatMap.append(tempList)


for i1,v1 in enumerate(mashList):
	for i2, v2 in enumerate(mashList):
		#pcall = "{} {} {}".format("echo", "hi", "blah")
		#proc = subprocess.Popen(["echo", "test"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if i2 <= i1:
			proc = subprocess.Popen([MASH, "dist", v1, v2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			output, err = proc.communicate()
			spout = output.split("\t")
			dist = float(spout[2])
			heatMap[i1][i2] = dist
			heatMap[i2][i1] = dist

fig, ax = plt.subplots()
im=ax.pcolor(np.array(heatMap), cmap='RdBu')
ax=plt.gca()
plt.colorbar(im, ax=ax)
plt.show()
fig.savefig(outName)
#print heatMap
