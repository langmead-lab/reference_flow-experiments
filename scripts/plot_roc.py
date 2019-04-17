import pandas as pd
import matplotlib.pyplot as plt

def get_tpr_fdr(df, n_var):
    g_nv = df.groupby('numvar').groups
    num_total = g_nv[n_var].shape[0]
    
    g_q = df.groupby('mapq').groups
    list_tp = []
    list_fd = []
    for qual in g_q:
        cnt = df.iloc[g_q[qual] & g_nv[n_var]]['dist'].value_counts()
        if len(list_tp) > 0:
            num_tp = list_tp[-1]
        else:
            num_tp = 0
        num_tp_in_group = 0
        for i in range(0, 11):
            if i in cnt:
                num_tp_in_group += cnt[i]
        # print (num_tp_in_group)
        try:
            if len(list_fd) > 0:
                num_fd = sum(cnt) - cnt[-3] - num_tp_in_group + list_fd[-1]
            else:
                num_fd = sum(cnt) - cnt[-3] - num_tp_in_group
        except:
            if len(list_fd) > 0:
                num_fd = sum(cnt) - num_tp_in_group + list_fd[-1]
            else:
                num_fd = sum(cnt) - num_tp_in_group

        list_tp.append(num_tp + num_tp_in_group)
        list_fd.append(num_fd)
    
    # tpr = num_tp / num_total
    # fdr = num_fd / num_total

    list_tpr = [i / num_total for i in list_tp]
    list_fdr = [i / num_total for i in list_fd]
    
    print (list_tpr)
    print (list_fdr)
    return list_tpr, list_fdr

df_ref = pd.read_pickle('10M-ref.sam-stats.pkl')
df_h37maj = pd.read_pickle('10M-h37maj.sam-stats.pkl')
df_per = pd.read_pickle('10M-per.sam-stats.pkl')

list_num_var = [0, 1, 2, 3]
plt.figure(1)
for i, num_var in enumerate(list_num_var):
    # num_var = 2
    plt.subplot(2, 2, i+1)
    list_tpr_ref, list_fdr_ref = get_tpr_fdr(df_ref, num_var)
    list_tpr_h37maj, list_fdr_h37maj = get_tpr_fdr(df_h37maj, num_var)
    list_tpr_per, list_fdr_per = get_tpr_fdr(df_per, num_var)

    plt.plot(list_fdr_ref, list_tpr_ref, marker='.', label='ref')
    plt.plot(list_fdr_h37maj, list_tpr_h37maj, marker='x', label='h37maj')
    plt.plot(list_fdr_per, list_tpr_per, marker='+', label='per (diploid)')
    plt.xscale('log')
    plt.ylabel('TPR (sensitivity)')
    plt.xlabel('FDR (1-precision)')
    plt.legend()
    plt.title('Chr21, number of vars = {0}'.format(num_var))
plt.show()
plt.clf()