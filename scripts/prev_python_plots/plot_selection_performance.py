import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

list_col = ['num_highq_reads', 'sensitivity_all', 'sensitivity_highq', 'sensitivity_lowq']
list_no_chr = [str(i) for i in list(range(1,23))]
list_chr = ['chr'+str(i) for i in list(range(1,23))]
df = pd.DataFrame(index=list_chr, columns=list_col)

dict_nhr = {}
with open('num_highq_wgs.txt', 'r') as f:
    ch = None
    for line in f:
        if ch == None and line.startswith('chr'):
            ch = line.rstrip()
        elif ch != None:
            dict_nhr[ch] = int(line.rstrip())
            ch = None
df['num_highq_reads'] = [dict_nhr[i] / 200000 for i in df.index]

dict_sa = {}
with open('sensitivity_all.txt', 'r') as f:
    ch = None
    for line in f:
        if ch == None and line.startswith('chr'):
            ch = line.rstrip()
        elif ch != None:
            tp = int(line.split()[3][1:])
            total = int(line.split()[5][:-1])
            dict_sa[ch] = float(tp / total)
            ch = None
df['sensitivity_all'] = [dict_sa[i] for i in df.index]

dict_sh = {}
with open('sensitivity_high.txt', 'r') as f:
    ch = None
    for line in f:
        if ch == None and line.startswith('chr'):
            ch = line.rstrip()
        elif ch != None:
            tp = int(line.split()[3][1:])
            total = int(line.split()[5][:-1])
            dict_sh[ch] = float(tp / total)
            ch = None
df['sensitivity_highq'] = [dict_sh[i] for i in df.index]

dict_sl = {}
with open('sensitivity_low.txt', 'r') as f:
    ch = None
    for line in f:
        if ch == None and line.startswith('chr'):
            ch = line.rstrip()
        elif ch != None:
            tp = int(line.split()[3][1:])
            total = int(line.split()[5][:-1])
            dict_sl[ch] = float(tp / total)
            ch = None
df['sensitivity_lowq'] = [dict_sl[i] for i in df.index]

plt.rcParams.update({'font.size': 20})

# plt.plot(list_no_chr, df['num_highq_reads'], label='fraction of q>=10 reads')
plt.plot(list_no_chr, df['sensitivity_highq'], label='MAPQ>=10')
plt.plot(list_no_chr, df['sensitivity_lowq'], label='MAPQ<10')
plt.plot(list_no_chr, df['sensitivity_all'], label='all')
# plt.yscale('log')
plt.ylabel('Sensitivity')
plt.xlabel('Chromosome')
plt.legend()
plt.tight_layout()
plt.show()
plt.clf()


# fig, ax1 = plt.subplots()
# ax1.bar(df.index, df['num_highq_reads'])
# ax1.set_ylabel('Number of reads')

# ax2 = ax1.twinx()
# ax2.plot(df.index, df['sensitivity_all'], label='all')
# ax2.plot(df.index, df['sensitivity_highq'], label='highq')
# ax2.plot(df.index, df['sensitivity_lowq'], label='lowq')
# ax2.set_yscale('log')
# ax2.set_ylabel('Sensitivity')

# plt.xlabel('Chromosome')
# plt.legend()
# plt.show()
# plt.clf()