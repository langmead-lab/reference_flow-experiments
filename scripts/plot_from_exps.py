import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

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
    
    ''' NA12878 '''
    # TP_OFFSET = 0
    # TP_OFFSET = 19719900 # num of TP in the first pass h37maj
    # TP_H37    = 19908912 # num of TP using purely h37
    # TP_H37MAJ = 19918604 # num of TP using h37maj (snvs+indels)
    # TP_PER    = 19929495 # num of TP using personalized genomes
    
    ''' HG03518 '''
    # TP_OFFSET = 0 # num of TP in the first pass h37maj
    TP_OFFSET = 19706590 # num of TP in the first pass h37maj
    TP_H37    = 19903321 # num of TP using purely h37
    TP_H37MAJ = 19913067 # num of TP using h37maj (snvs+indels)
    TP_PER    = 19929234 # num of TP using personalized genomes
    
    NUM_READS = 20000000
    

    list_tp = read_true_pos(log_fn)
    list_samples = read_samples(sample_fn)
    data_df = pd.DataFrame(list_tp)
    data_df.columns = ['Num TP (Second Pass)']
    data_df.index = list_samples
    dict_pop = read_ped(ped_fn)
    list_pop = [dict_pop[data_df.index[i]] for i in range(data_df.shape[0])]
    data_df['Population'] = list_pop
    data_df['Num TP'] = data_df['Num TP (Second Pass)'] + TP_OFFSET

    #: special items people would like to compare
    data_df.loc['h37'] = [TP_H37-TP_OFFSET, 'h37', TP_H37]
    data_df.loc['personalized'] = [TP_PER-TP_OFFSET, 'personalized', TP_PER]
    data_df.loc['h37maj'] = [TP_H37MAJ-TP_OFFSET, 'h37maj', TP_H37MAJ]
    
    data_df['Sensitivity'] = data_df['Num TP'] / NUM_READS

    superpop_df = pd.read_csv(superpop_fn, sep='\t', header=None, index_col=0)
    superpop_groups = superpop_df.groupby(2)
    dict_superpop = {}
    for n, g in superpop_groups:
        for i in superpop_groups.groups[n]:
            dict_superpop[i] = n
    dict_pop['h37'] = 'h37'
    dict_pop['h37maj'] = 'h37maj'
    dict_pop['personalized'] = 'personalized'
    dict_superpop['h37'] = 'h37'
    dict_superpop['h37maj'] = 'h37maj'
    dict_superpop['personalized'] = 'personalized'
    dict_color = {'AFR': 'C0', 'AMR': 'C1', 'CEU': 'C2', 'EAS': 'C3', 'EUR': 'C4', 'SAS': 'C5', 'h37': 'C6', 'h37maj': 'C7', 'personalized': 'C8'}

    data_df = data_df.sort_values(by='Sensitivity')
    list_superpop = []
    for i in data_df['Population']:
        list_superpop.append(dict_superpop[i])
    list_superpop
    data_df['Superpopulation'] = list_superpop

    #: sort categories by the median of each superpop
    order_by_superpop = list(data_df.groupby('Superpopulation').median().sort_values('Sensitivity').index)

    ''' Bar plot '''
    # x = np.arange(data_df.shape[0])
    # bar = plt.bar(x, data_df['Sensitivity'])
    # plt.ylim(bottom=1.1*min(data_df['Sensitivity'])-0.1*max(data_df['Sensitivity']), top=max(data_df['Sensitivity']))
    # existed_superpop = []
    # for i in range(data_df.shape[0]):
    #     pop = dict_pop[data_df.index[i]]
    #     superpop = dict_superpop[pop]
    #     if superpop not in existed_superpop:
    #         existed_superpop.append(superpop)
    #         bar[i].set_label(superpop)
    #     color = dict_color[superpop]
    #     bar[i].set_color(color)
    # plt.legend()
    # plt.ylabel('Sensitivity')
    # plt.title('First pass: h37maj; second pass: 100 random indivs; mapq_th=10 (NA12878)')
    # # plt.title('First pass: h37maj; second pass: 100 random indivs (hap); mapq_th=10 (NA12878)')
    # # plt.title('Aligned against 100 random indivs (NA12878)')
    # # plt.title('Aligned against 100 random indivs (HG03518)')
    # # plt.title('First pass: h37maj; second pass: 100 random indivs; mapq_th=10 (HG03518)')
    # plt.show()
    # plt.clf()

    ''' Box plot '''
    # sns.boxplot(x='Population', y='Sensitivity', hue='Superpopulation', data=data_df)
    # sns.swarmplot(x='Population', y='Sensitivity', hue='Superpopulation', data=data_df)
    sns.boxplot(x='Superpopulation', y='Sensitivity', data=data_df, order=order_by_superpop)
    sns.swarmplot(x='Superpopulation', y='Sensitivity', data=data_df, order=order_by_superpop)
    # plt.title('First pass: h37maj; second pass: 100 random indivs; mapq_th=10 (NA12878)')
    # plt.title('First pass: h37maj; second pass: 100 random indivs (hap); mapq_th=10 (NA12878)')
    # plt.title('Aligned against 100 random indivs (NA12878)')
    # plt.title('Aligned against 100 random indivs (HG03518)')
    plt.title('First pass: h37maj; second pass: 100 random indivs; mapq_th=10 (HG03518)')
    plt.show()
    plt.clf()

    # sns.barplot(x='Population', y='Sensitivity', hue='Superpopulation', data=data_df)
    # plt.ylim(bottom=1.1*min(data_df['Sensitivity'])-0.1*max(data_df['Sensitivity']), top=max(data_df['Sensitivity']))
    # plt.title('First pass: h37maj; second pass: 100 random indivs; mapq_th=10 (NA12878)')
    # plt.show()
    # plt.clf()