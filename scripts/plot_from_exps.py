import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s', '--samples',
        help='txt file recording the order of samples used in experiments'
    )
    parser.add_argument(
        '-l', '--log',
        help='log file for experiments'
    )
    parser.add_argument(
        '-p', '--ped',
        help='ped file recoding population and family info'
    )
    parser.add_argument(
        '-sp', '--superpopulation',
        help='superpopulation table'
    )
    return parser.parse_args()

def read_true_pos(log_fn):
    list_tp = []
    with open(log_fn, 'r') as f:
        for line in f:
            if line.startswith('sensitivity_all'):
                tp = int(line.split()[3][1:])
                list_tp.append(tp)
    return list_tp

def read_samples(sample_fn):
    list_samples = []
    with open(sample_fn, 'r') as f:
        for line in f:
            list_samples.append(line.rstrip())
    return list_samples

def read_ped(ped_fn):
    '''
    Reads the ped file and returns a dict where keys are indivs and 
    values are corresponding populations
    '''
    ped_df = pd.read_csv(ped_fn, sep='\t')
    dd = ped_df[['Individual ID', 'Population']]

    popd = {}
    for i in range(dd.shape[0]):
        key = dd['Individual ID'][i]
        value = dd['Population'][i]
        popd[key] = value

    return popd

if __name__ == '__main__':
    args = parse_args()
    sample_fn = args.samples
    log_fn = args.log
    ped_fn = args.ped
    superpop_fn = args.superpopulation
    
    TP_OFFSET = 19719900 # num of TP in the first pass
    TP_FLOOR  = 19908912 # num of TP using purely hg19
    TP_CEIL   = 19929495 # num of TP using personalized genomes

    list_tp = read_true_pos(log_fn)
    list_samples = read_samples(sample_fn)
    data_df = pd.DataFrame(list_tp)
    data_df.columns = ['Num TP']
    data_df.index = list_samples
    dict_pop = read_ped(ped_fn)
    list_pop = [dict_pop[data_df.index[i]] for i in range(data_df.shape[0])]
    data_df['Population'] = list_pop
    data_df['Num TP'] += TP_OFFSET

    list_tp = [TP_FLOOR] + list(data_df['Num TP']) + [TP_CEIL]
    list_sample = ['hg19'] + list(data_df.index) + ['personalized']

    superpop_df = pd.read_csv(superpop_fn, sep='\t', header=None, index_col=0)
    superpop_groups = superpop_df.groupby(2)
    dict_superpop = {}
    for n, g in superpop_groups:
        for i in superpop_groups.groups[n]:
            dict_superpop[i] = n
    dict_pop['hg19'] = 'hg19'
    dict_pop['personalized'] = 'CEU'
    dict_superpop['hg19'] = 'hg19'

    dict_color = {'AFR': 'C0', 'AMR': 'C1', 'CEU': 'C2', 'EAS': 'C3', 'EUR': 'C4', 'SAS': 'C5', 'hg19': 'C6'}

    list_tp.sort()
    x = np.arange(len(list_tp))
    bar = plt.bar(x, list_tp)
    plt.ylim(bottom=1.1*min(list_tp)-0.1*max(list_tp), top=max(list_tp))
    existed_superpop = []
    for i in range(len(list_sample)):
        pop = dict_pop[list_sample[i]]
        superpop = dict_superpop[pop]
        if superpop not in existed_superpop:
            existed_superpop.append(superpop)
            bar[i].set_label(superpop)
        color = dict_color[superpop]
        bar[i].set_color(color)
    plt.legend()
    plt.title('Aligned to h37maj and then 100 random indivs; mapq_th=10')
    plt.show()
    plt.clf()