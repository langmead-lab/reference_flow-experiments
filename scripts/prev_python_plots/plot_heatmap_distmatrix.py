import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import re, sys
import os, glob

#from plot_from_exps import read_true_pos
inputFile=sys.argv[1]
outName=sys.argv[2]

inputMatrix = np.load(inputFile)

fig, ax = plt.subplots()
im=ax.pcolor(inputMatrix, cmap='RdBu')
ax=plt.gca()
plt.colorbar(im, ax=ax)
plt.show()
fig.savefig(outName)


#log_fn = ''

#list_tp = read_true_pos(log_fn)
#list_s1, list_s2 = read_sample(log_fn)
#df = pd.DataFrame(index=list(set(list_s1)), columns=list(set(list_s2)))
#for i, t in enumerate(list_tp):
#    df.loc[list_s1[i], list_s2[i]] = t

#: TODO
#df1 = df.drop(index='HG00628').dropna(axis=1).fillna(0)
#sns.heatmap(df1, cmap='YlGnBu')
#plt.show()
