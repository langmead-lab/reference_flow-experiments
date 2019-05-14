import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
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

# def read_truepos_and_sample(log_fn):
#     list_tp = []
#     list_sample = []
#     with open(log_fn, 'r') as f:
#         for line in f:
#             if line.startswith('sample'):
#                 samp = line.split('=')[1].rstrip()
#                 list_sample.append(samp)
#             #: report sensitivity
#             if line.startswith('sensitivity_all'):
#                 tp = int(line.split()[3][1:])
#                 list_tp.append(tp)
#     return list_tp, list_sample

def read_unaligned(log_fn):
    list_unaligned = []
    with open(log_fn, 'r') as f:
        for line in f:
            #: report sensitivity
            if line.startswith('Unaligned:'):
                un = int(line.split()[1])
                list_unaligned.append(un)
            # if line.startswith('unaligned            ='):
            #     un = int(line.split()[3][1:])
            #     list_unaligned.append(un)
            # else:
            #     print ('F')
    return list_unaligned

def read_true_pos(log_fn):
    list_tp = []
    with open(log_fn, 'r') as f:
        for line in f:
            #: report sensitivity
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

    PLOT_UNALIGNED = False
    #: constants
    if PLOT_UNALIGNED:
        NUM_READS = 1
    else:
        NUM_READS = 20000000
    if PLOT_UNALIGNED:
        dict_NA12878 = {'linear': 10196, 'major': 6939, 'personalized': 3849}#, 'major-personalized': 2655}
    else:
        dict_NA12878 = {'linear': 19908912, 'major': 19918604, 'personalized': 19929501}#, 'major-personalized': 19932867}
    dict_HG03518 = {'linear': 19903321, 'major': 19913067, 'personalized': 19929234}#, 'major-personalized': 20000000}
    dict_samples = {'NA12878': dict_NA12878, 'HG03518': dict_HG03518}
    LABEL_GROUPS = ['Population', 'Superpopulation', 'Cluster', 'Merge']

    dict_pop = read_ped(ped_fn)

    df_exp_spec = pd.read_csv(tsv_exp_fn, sep='\t')
    list_exp = []
    for name_exp in df_exp_spec['filename']:
        name_exp = name_exp.split('/')[1]
        list_exp.append(name_exp[:name_exp.rfind('.log')])
    list_indivs = read_indivs(sample_fn)
    num_single_indiv = 0
    for i in list_indivs:
        if i.count(',') == 0:
            num_single_indiv += 1
    df_data = pd.DataFrame(columns=list_exp, index=list_indivs)

    #: adds population info
    list_pop = []
    for i, indiv in enumerate(df_data.index):
        if i < num_single_indiv:
            list_pop.append(dict_pop[indiv])
        else:
            pop1 = dict_pop[indiv.split(',')[0]]
            pop2 = dict_pop[indiv.split(',')[1]]
            if pop1 <= pop2:
                list_pop.append(pop1 + ',' + pop2)
            else:
                list_pop.append(pop2 + ',' + pop1)
    df_data.insert(df_data.shape[1], 'Population', list_pop)

    #: adds superpopulation info
    superpop_df = pd.read_csv(superpop_fn, sep='\t', header=None, index_col=0)
    superpop_groups = superpop_df.groupby(2)
    dict_superpop = {}
    for n, g in superpop_groups:
        for i in superpop_groups.groups[n]:
            dict_superpop[i] = n
    #: add a superpop column
    list_superpop = [dict_superpop[df_data['Population'][i]] for i in range(num_single_indiv)]
    list_superpop = []
    for i, indiv in enumerate(df_data.index):
        if i < num_single_indiv:
            list_superpop.append(dict_superpop[df_data['Population'][i]])
        else:
            superpop1 = dict_superpop[df_data['Population'][i].split(',')[0]]
            superpop2 = dict_superpop[df_data['Population'][i].split(',')[1]]
            if superpop1 <= superpop2:
                list_superpop.append(superpop1 + ',' + superpop2)
            else:
                list_superpop.append(superpop2 + ',' + superpop1)
    df_data.insert(df_data.shape[1], 'Superpopulation', list_superpop)

    #: process clusrterting data
    df_cluster = pd.read_csv(cluster_fn, sep=': ', header=None)
    dict_cluster = {}
    for i, name in enumerate(df_cluster[0]):
        dict_cluster[name] = str(df_cluster[1][i])
    list_cluster = []
    for i, indiv in enumerate(df_data.index):
        if i < num_single_indiv:
            list_cluster.append(dict_cluster[indiv])
        else:
            cluster1 = dict_cluster[indiv.split(',')[0]]
            cluster2 = dict_cluster[indiv.split(',')[1]]
            if cluster1 <= cluster2:
                list_cluster.append(cluster1 + ',' + cluster2)
            else:
                list_cluster.append(cluster2 + ',' + cluster1)
    df_data.insert(df_data.shape[1], 'Cluster', list_cluster)

    #: adds merging info
    list_merge = ['major-diploid-mergeAS' if i.count(',') > 0 else 'major-diploid' for i in df_data.index]
    df_data.insert(df_data.shape[1], 'Merge', list_merge)

    list_included_samples = []

    for id_exp in range(df_exp_spec.shape[0]):
        sample = df_exp_spec['sample_id'][id_exp]
        #: add 'linear', 'major', 'per' numbers for a new sample
        if sample not in list_included_samples:
            l_id = []
            
            TPR_H37 = dict_samples[sample]['linear'] / NUM_READS
            l_h37 = [TPR_H37] * (len(df_data.columns) - len(LABEL_GROUPS))
            for i in range(len(LABEL_GROUPS)):
                l_h37.append('linear')
            l_id.append(sample+'-linear')
            
            TPR_H37MAJ = dict_samples[sample]['major'] / NUM_READS
            l_h37maj = [TPR_H37MAJ] * (len(df_data.columns) - len(LABEL_GROUPS))
            for i in range(len(LABEL_GROUPS)):
                l_h37maj.append('major')
            l_id.append(sample+'-major')
            
            TPR_PER = dict_samples[sample]['personalized'] / NUM_READS
            l_per = [TPR_PER] * (len(df_data.columns) - len(LABEL_GROUPS))
            for i in range(len(LABEL_GROUPS)):
                l_per.append('personalized')
            l_id.append(sample+'-personalized')

            # TPR_PER2 = dict_samples[sample]['major-personalized'] / NUM_READS
            # l_per2 = [TPR_PER2] * (len(df_data.columns) - len(LABEL_GROUPS))
            # for i in range(len(LABEL_GROUPS)):
            #     l_per2.append('major-personalized')
            # l_id.append(sample+'-major-personalized')

            # ll = [l_h37, l_h37maj, l_per, l_per2]
            special_items = [l_h37, l_h37maj, l_per]
            df_special_items = pd.DataFrame(special_items, columns=list(df_data.columns), index=l_id)
            df_data = df_data.append(df_special_items)
            list_included_samples.append(sample)

        TP_OFFSET = int(df_exp_spec['1pass-truepos'][id_exp])
        log_fn = df_exp_spec['filename'][id_exp]
        
        if PLOT_UNALIGNED:
            list_tp = read_unaligned(log_fn)
            list_tpr = [(i) / NUM_READS for i in list_tp]
        else:
            list_tp = read_true_pos(log_fn)
            list_tpr = [(i + TP_OFFSET) / NUM_READS for i in list_tp]

        # df_data[list_exp[id_exp]][0:len(list_tp)] = list_tpr
        if df_exp_spec['2pass-merge'][id_exp] == 'False':
            df_data[list_exp[id_exp]][0:num_single_indiv] = list_tpr
        else:
            df_data[list_exp[id_exp]][num_single_indiv:len(list_indivs)] = list_tpr


    df_data['all'] = df_data[list_exp[0]].combine_first(df_data[list_exp[1]])

    #: select plot target
    # plot_target = list_exp[2]
    plot_target = 'all'

    #: sort categories by the median of each cluster
    order_by_cluster = list(df_data.groupby('Cluster').median().sort_values(plot_target).index)
    #: sort categories by the median of each superpop
    order_by_superpop = list(df_data.groupby('Superpopulation').median().sort_values(plot_target).index)
    #: sort categories by the median of each merging method
    order_by_merge = list(df_data.groupby('Merge').median().sort_values(plot_target).index)

    df_data['sort_by_superpop'] = [order_by_superpop.index(i) for i in df_data['Superpopulation']]
    df2 = df_data[:-len(LABEL_GROUPS)].sort_values(by='sort_by_superpop', ascending=True)

    PLOT_MERGE_ONLY = False
    #: plot personalized and merging results
    if PLOT_MERGE_ONLY:
        df2 = df_data[num_single_indiv:-3].sort_values(by='sort_by_superpop', ascending=True)
        df2 = df2.append(df_data.iloc[-1])
        # df2 = df_data[num_single_indiv:-4].sort_values(by='sort_by_superpop', ascending=True)
        # df2 = df2.drop(['NA12878-linear', 'NA12878-major', 'personalized'])
        # order_by_merge = list(df2.groupby('Merge').median().sort_values(plot_target).index)
    PLOT_WITHOUT_MERGE = False
    if PLOT_WITHOUT_MERGE:
        df3 = df_data.iloc[:num_single_indiv]
        df3 = df3.append(df_data.iloc[-3:])
        order_by_superpop = list(df3.groupby('Superpopulation').median().sort_values(plot_target).index)
    #: TODO
    # LABEL_GROUPS -= 2

    params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
        #  'xtick.labelsize':'x-large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)

    ''' Bar plot '''
    # x = np.arange(data_df.shape[0])
    # bar = plt.bar(x, data_df['Sensitivity'])
    # plt.ylim(bottom=1.1*min(data_df['Sensitivity'])-0.1*max(data_df['Sensitivity']), top=max(data_df['Sensitivity']))
    # existed_superpop = []
    # dict_color = {'AFR': 'C0', 'AMR': 'C1', 'CEU': 'C2', 'EAS': 'C3', 'EUR': 'C4', 'SAS': 'C5', 'linear': 'C6', 'major': 'C7', 'personalized': 'C8'}
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
    # plt.title('First pass: major; second pass: 100 random indivs; mapq_th=10 (NA12878)')
    # # plt.title('First pass: major; second pass: 100 random indivs (hap); mapq_th=10 (NA12878)')
    # # plt.title('Aligned against 100 random indivs (NA12878)')
    # # plt.title('Aligned against 100 random indivs (HG03518)')
    # # plt.title('First pass: major; second pass: 100 random indivs; mapq_th=10 (HG03518)')
    # plt.show()
    # plt.clf()

    ''' Box plot '''
    #: by population
    # sns.boxplot(x='Population', y='Sensitivity', hue='Superpopulation', data=data_df)
    # sns.swarmplot(x='Population', y='Sensitivity', hue='Superpopulation', data=data_df)
    
    #: by superpopulation
    # sns.boxplot(x='Superpopulation', y=plot_target, data=df_data.iloc[:-len(LABEL_GROUPS)], order=order_by_superpop)
    # sns.swarmplot(x='Superpopulation', y=plot_target, data=df_data.iloc[:-len(LABEL_GROUPS)], order=order_by_superpop)
    # sns.boxplot(x='Superpopulation', y=plot_target, data=df_data, order=order_by_superpop)
    # sns.swarmplot(x='Superpopulation', y=plot_target, data=df_data, order=order_by_superpop)
    
    #: by cluster
    # sns.boxplot(x='Cluster', y='Sensitivity', data=data_df, order=order_by_cluster)
    # sns.swarmplot(x='Cluster', y='Sensitivity', data=data_df, order=order_by_cluster)
    
    #: by merge
    if PLOT_UNALIGNED:
        sns.boxplot(x='Merge', y='all', data=df_data.iloc[-len(special_items):], order=order_by_merge)
        sns.boxplot(x='Merge', y='all', hue='Superpopulation', data=df_data.iloc[num_single_indiv:-len(special_items)], order=order_by_merge)
        sns.boxplot(x='Merge', y='all', hue='Superpopulation', data=df_data.iloc[:num_single_indiv], order=order_by_merge)
    elif PLOT_MERGE_ONLY:
        sns.boxplot(x='Superpopulation', y='all', data=df2)
        # sns.boxplot(x='Merge', y='all', data=df_data.iloc[-1:], order=order_by_merge)
        # # sns.boxplot(x='Merge', y='all', data=df_data.iloc[-2:], order=order_by_merge)
        # sns.boxplot(x='Merge', y='all', hue='Superpopulation', data=df2.iloc[:-2], order=order_by_merge)
    elif PLOT_WITHOUT_MERGE:
        sns.boxplot(x='Superpopulation', y='all', data=df3, order=order_by_superpop)
    else:
        sns.boxplot(x='Merge', y='all', data=df_data.iloc[-len(special_items):], order=order_by_merge)
        sns.boxplot(x='Merge', y='all', hue='Superpopulation', data=df2.iloc[num_single_indiv:-len(special_items)], order=order_by_merge)
        sns.boxplot(x='Merge', y='all', hue='Superpopulation', data=df2.iloc[:num_single_indiv], order=order_by_merge)
    
    # sns.boxplot(x='Merge', y='all', hue='Superpopulation', data=df_data.iloc[num_single_indiv:-len(special_items)], order=order_by_merge)
    # sns.boxplot(x='Merge', y='all', hue='Superpopulation', data=df_data.iloc[:num_single_indiv], order=order_by_merge)
    # sns.swarmplot(x='Merge', y='NA12878-major-q10-100dip-merged_random_partial', data=df_data.iloc[num_single_indiv:-len(special_items)], order=order_by_merge)
    # sns.swarmplot(x='Merge', y='NA12878-major-q10-100dip', data=df_data.iloc[:num_single_indiv], order=order_by_merge)

    
    # plt.title(plot_target)
    plt.xlabel('')
    if PLOT_UNALIGNED:
        plt.ylabel('#Unaligned')
    else:
        plt.ylabel('Sensitivity')
    # plt.gca().get_legend().remove()
    # plt.legend(loc='lower right', ncol=2)
    plt.show()
    plt.clf()

    # sns.barplot(x='Population', y='Sensitivity', hue='Superpopulation', data=data_df)
    # plt.ylim(bottom=1.1*min(data_df['Sensitivity'])-0.1*max(data_df['Sensitivity']), top=max(data_df['Sensitivity']))
    # plt.title('First pass: major; second pass: 100 random indivs; mapq_th=10 (NA12878)')
    # plt.show()
    # plt.clf()

    #: save data to a tsv
    # df_data.to_csv('experiments.tsv', sep='\t')