# import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
# from plot_from_exps import read_true_pos

def read_true_pos(log_fn):
    list_tp = []
    with open(log_fn, 'r') as f:
        for line in f:
            #: report sensitivity
            if line.startswith('sensitivity_all'):
                tp = int(line.split()[3][1:])
                list_tp.append(tp)
    return list_tp

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

log_fn = 'NA12878-h37maj-q10-100dip-merged_random_partial.log'
dist_fn = '../reduced_min02_chr21_hg38_dist.csv'
samples_fn = 'better_list.txt'

list_tp = read_true_pos(log_fn)
list_s1, list_s2 = read_sample(log_fn)
df = pd.DataFrame(index=list(set(list_s1)), columns=list(set(list_s2)))
for i, t in enumerate(list_tp):
    df.loc[list_s1[i], list_s2[i]] = t

list_samples_for_dist = pd.read_csv(samples_fn, header=None)
df_dist = pd.read_csv(dist_fn, header=None)
df_dist.index = list_samples_for_dist[0]
df_dist.columns = list_samples_for_dist[0]
dict_NA12878 = {}
for i, e in enumerate(df_dist.loc['NA12878']):
    dict_NA12878[df_dist.index[i]] = e

#: TODO
df1 = df.drop(index='HG01131').dropna(axis=1).fillna(0)
df1['sort_index'] = [dict_NA12878[i] for i in df1.index]
df1 = df1.sort_values(by='sort_index')
sns.heatmap(df1.iloc[:-1,:-1], cmap='YlGnBu')
# sns.heatmap(df1, cmap='YlGnBu')
plt.show()