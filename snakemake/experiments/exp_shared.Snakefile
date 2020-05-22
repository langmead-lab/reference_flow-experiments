'''
Miscellaneous rules for experiments
'''
def organize_accuracy(fn_input, fn_output):
    with open(fn_input, 'r') as f:
        list_tp = []
        list_all = []
        for line in f:
            if line.count('sensitivity_all') > 0:
                line = line.split()
                list_tp.append(int(line[3][1:]))
                list_all.append(int(line[5][:-1]))
    f_out = open(fn_output, 'w')
    f_out.write('{0}\n{1}\n{2}\n'.format(sum(list_tp), sum(list_all), sum(list_tp) / sum(list_all)))
    return

def summarize_allelc_bias(fn_input, fn_output):
    df = pd.read_csv(str(fn_input), sep='\t')
    len_raw = len(df)
    #: filter by BIAS_MIN_READ_COUNT
    df_f = df[(df['REF_COUNT'] + df['ALT_COUNT']) >= BIAS_MIN_READ_COUNT]
    len_filtered = len(df_f)

    ref = df_f['REF_COUNT']
    alt = df_f['ALT_COUNT']
    ref_to_alt = sum(ref) / sum(alt)
    num_refbias = len(df_f[ref / (ref+alt) >= 0.5 + BIAS_TAIL_THRDS])
    num_altbias = len(df_f[ref / (ref+alt) <= 0.5 - BIAS_TAIL_THRDS])

    with open(str(fn_output), 'w') as f:
        # Fraction of sites passing the filter
        f.write(str(len_filtered / len_raw) + '\n')
        # Sum of REFs
        f.write(str(sum(ref)) + '\n')
        # Sum of ALTs
        f.write(str(sum(alt)) + '\n')
        # REF-to-ALT
        f.write(str(ref_to_alt) + '\n')
        # Num of REF biases
        f.write(str(num_refbias) + '\n')
        # Num of ALT biases
        f.write(str(num_altbias) + '\n')
    return

def retrieve_indiv_variants(indiv, fn_input, fn_output):
    df_out = pd.DataFrame()

    list_rcd = []
    list_snp = []
    list_indel = []
    for i, fn in enumerate(fn_input):
        if fn.count(INDIV[i]) <= 0:
            exit(1)
        with open(fn, 'r') as f:
            for line in f:
                line = line.split('\t')
                if line[0] == 'SN' and line[2] == 'number of records:':
                    num_rcd = int(line[3])
                elif line[0] == 'SN' and line[2] == 'number of SNPs:':
                    num_snp = int(line[3])
                elif line[0] == 'SN' and line[2] == 'number of indels:':
                    num_indel = int(line[3])
                    break
        list_rcd.append(num_rcd)
        list_snp.append(num_snp)
        list_indel.append(num_indel)

    dict_indiv_to_pop = build_dict_indiv_to_pop(FAMILY)
    dict_pop_to_spop = build_dict_pop_to_spop(SPOP)

    df_out['INDIV'] = INDIV
    df_out['SUPER_POP'] = [dict_pop_to_spop[dict_indiv_to_pop[s]] for s in INDIV]
    df_out['POP'] = [dict_indiv_to_pop[s] for s in INDIV]
    df_out['RECORD'] = list_rcd
    df_out['NUM_SNP'] = list_snp
    df_out['NUM_INDEL'] = list_indel

    df_out.to_csv(str(fn_output), sep='\t', index=None)

    return

