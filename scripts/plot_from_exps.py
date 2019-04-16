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
        '--tsv_exp',
        help='csv specifying all experiments of interest'
    )
    parser.add_argument(
        '-p', '--ped',
        help='ped file recoding population and family info'
    )
    parser.add_argument(
        '-sp', '--superpopulation',
        help='superpopulation table'
    )
    parser.add_argument(
        '-c', '--cluster',
        help='txt file specifying unsupervised clusters'
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

def read_indivs(sample_fn):
    list_indivs = []
    with open(sample_fn, 'r') as f:
        for line in f:
            list_indivs.append(line.rstrip())
    return list_indivs

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
    tsv_exp_fn = args.tsv_exp
    ped_fn = args.ped
    cluster_fn = args.cluster
    superpop_fn = args.superpopulation

    #: constants
    NUM_READS = 20000000
    dict_NA12878 = {'h37': 19908912, 'h37maj': 19918604, 'personalized': 19929495}
    dict_HG03518 = {'h37': 19903321, 'h37maj': 19913067, 'personalized': 19929234}
    dict_samples = {'NA12878': dict_NA12878, 'HG03518': dict_HG03518}
    LABEL_GROUPS = ['Population', 'Superpopulation', 'Cluster']

    dict_pop = read_ped(ped_fn)

    df_exp_spec = pd.read_csv(tsv_exp_fn, sep='\t')
    list_exp = []
    for name_exp in df_exp_spec['filename']:
        name_exp = name_exp.split('/')[1]
        list_exp.append(name_exp[:name_exp.rfind('.log')])
    list_indivs = read_indivs(sample_fn)
    df_data = pd.DataFrame(columns=list_exp, index=list_indivs)

    #: adds population info
    list_pop = [dict_pop[df_data.index[i]] for i in range(df_data.shape[0])]
    df_data.insert(df_data.shape[1], 'Population', list_pop)

    #: adds superpopulation info
    superpop_df = pd.read_csv(superpop_fn, sep='\t', header=None, index_col=0)
    superpop_groups = superpop_df.groupby(2)
    dict_superpop = {}
    for n, g in superpop_groups:
        for i in superpop_groups.groups[n]:
            dict_superpop[i] = n
    #: add a superpop column
    list_superpop = []
    for i in df_data['Population']:
        list_superpop.append(dict_superpop[i])
    df_data.insert(df_data.shape[1], 'Superpopulation', list_superpop)

    #: process clusrterting data
    df_cluster = pd.read_csv(cluster_fn, sep=': ', header=None)
    dict_cluster = {}
    for i, name in enumerate(df_cluster[0]):
        dict_cluster[name] = df_cluster[1][i]
    c_label = []
    for name in df_data.index:
        c_label.append(dict_cluster[name])
    df_data.insert(df_data.shape[1], 'Cluster', c_label)

    list_included_samples = []

    for id_exp in range(df_exp_spec.shape[0]):
        sample = df_exp_spec['sample_id'][id_exp]
        #: add 'h37', 'h37maj', 'per' numbers for a new sample
        if sample not in list_included_samples:
            l_id = []
            
            TPR_H37 = dict_samples[sample]['h37'] / NUM_READS
            l_h37 = [TPR_H37] * (len(df_data.columns) - len(LABEL_GROUPS))
            for i in range(len(LABEL_GROUPS)):
                l_h37.append('h37')
            l_id.append(sample+'-h37')
            
            TPR_H37MAJ = dict_samples[sample]['h37maj'] / NUM_READS
            l_h37maj = [TPR_H37MAJ] * (len(df_data.columns) - len(LABEL_GROUPS))
            for i in range(len(LABEL_GROUPS)):
                l_h37maj.append('h37maj')
            l_id.append(sample+'-h37maj')
            
            TPR_PER = dict_samples[sample]['personalized'] / NUM_READS
            l_per = [TPR_PER] * (len(df_data.columns) - len(LABEL_GROUPS))
            for i in range(len(LABEL_GROUPS)):
                l_per.append('personalized')
            l_id.append(sample+'-personalized')

            ll = [l_h37, l_h37maj, l_per]
            df_ll = pd.DataFrame(ll, columns=list(df_data.columns), index=l_id)
            df_data = df_data.append(df_ll)
            list_included_samples.append(sample)

        TP_OFFSET = int(df_exp_spec['1pass-truepos'][id_exp])
        log_fn = df_exp_spec['filename'][id_exp]
        list_tp = read_true_pos(log_fn)
        list_tpr = [(i + TP_OFFSET) / NUM_READS for i in list_tp]
        df_data[list_exp[id_exp]][0:len(list_tp)] = list_tpr
        

    plot_target = list_exp[2]

    #: sort categories by the median of each cluster
    order_by_cluster = list(df_data.groupby('Cluster').median().sort_values(plot_target).index)
    #: sort categories by the median of each superpop
    order_by_superpop = list(df_data.groupby('Superpopulation').median().sort_values(plot_target).index)

    ''' Bar plot '''
    # x = np.arange(data_df.shape[0])
    # bar = plt.bar(x, data_df['Sensitivity'])
    # plt.ylim(bottom=1.1*min(data_df['Sensitivity'])-0.1*max(data_df['Sensitivity']), top=max(data_df['Sensitivity']))
    # existed_superpop = []
    # dict_color = {'AFR': 'C0', 'AMR': 'C1', 'CEU': 'C2', 'EAS': 'C3', 'EUR': 'C4', 'SAS': 'C5', 'h37': 'C6', 'h37maj': 'C7', 'personalized': 'C8'}
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
    #: by population
    # sns.boxplot(x='Population', y='Sensitivity', hue='Superpopulation', data=data_df)
    # sns.swarmplot(x='Population', y='Sensitivity', hue='Superpopulation', data=data_df)
    
    #: by superpopulation
    sns.boxplot(x='Superpopulation', y=plot_target, data=df_data.iloc[:-len(LABEL_GROUPS)], order=order_by_superpop)
    sns.swarmplot(x='Superpopulation', y=plot_target, data=df_data.iloc[:-len(LABEL_GROUPS)], order=order_by_superpop)
    
    #: by cluster
    # sns.boxplot(x='Cluster', y='Sensitivity', data=data_df, order=order_by_cluster)
    # sns.swarmplot(x='Cluster', y='Sensitivity', data=data_df, order=order_by_cluster)
    plt.title(plot_target)
    plt.ylabel('Sensitivity')
    plt.show()
    plt.clf()

    # sns.barplot(x='Population', y='Sensitivity', hue='Superpopulation', data=data_df)
    # plt.ylim(bottom=1.1*min(data_df['Sensitivity'])-0.1*max(data_df['Sensitivity']), top=max(data_df['Sensitivity']))
    # plt.title('First pass: h37maj; second pass: 100 random indivs; mapq_th=10 (NA12878)')
    # plt.show()
    # plt.clf()

    #: save data to a tsv
    df_data.to_csv('experiments.tsv', sep='\t')