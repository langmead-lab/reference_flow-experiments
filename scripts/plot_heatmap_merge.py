import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from plot_from_exps import read_true_pos

def read_sample(log_fn):
    list_sample1 = []
    list_sample2 = []
    with open(log_fn, 'r') as f:
        for line in f:
            if line.startswith('sample'):
                sample1 = line.split('=')[1].split(',')[0].rstrip()
                sample2 = line.split('=')[1].split(',')[1].rstrip()
                list_sample1.append(sample1)
                list_sample2.append(sample2)
    return list_sample1, list_sample2

log_fn = ''

list_tp = read_true_pos(log_fn)
list_s1, list_s2 = read_sample(log_fn)
df = pd.DataFrame(index=list(set(list_s1)), columns=list(set(list_s2)))
for i, t in enumerate(list_tp):
    df.loc[list_s1[i], list_s2[i]] = t

#: TODO
df1 = df.drop(index='HG00628').dropna(axis=1).fillna(0)
sns.heatmap(df1, cmap='YlGnBu')
plt.show()