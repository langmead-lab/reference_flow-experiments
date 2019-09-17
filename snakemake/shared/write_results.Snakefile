'''
Summarize results as a TSV
'''
rule write_as_tsv:
    input:
        acc = expand(os.path.join(DIR_RESULTS,
            '{INDIV}' + '-{0}-h37maj-{1}-{2}.acc'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            # '{INDIV}-' + CHROM + '-h37maj-' + POP_DIRNAME + '.acc'),
            INDIV = INDIV),
        check_onepass_acc = expand(os.path.join(DIR, '{INDIV}_onepass_acc.done'), INDIV = INDIV)
    output:
        tsv = os.path.join(DIR_RESULTS, 'all.tsv'),
        check = touch(temp(os.path.join(DIR, 'accuracy.done')))
    run:
        dict_indiv_to_pop = build_dict_indiv_to_pop(FAMILY)
        dict_pop_to_spop = build_dict_pop_to_spop(SPOP)
        
        df = pd.DataFrame()
        list_fn = os.listdir(DIR_RESULTS)
        list_fn = [os.path.join(DIR_RESULTS, fn) for fn in list_fn]
        list_exp = []
        list_indiv = []
        list_method = []
        list_tp = []
        list_all = []
        list_sensitivity = []
        for fn in list_fn:
            if fn.endswith('.acc'):
                with open(fn, 'r') as f:
                    for i, line in enumerate(f):
                        if i == 0:
                            list_tp.append(int(line))
                        elif i == 1:
                            list_all.append(int(line))
                fn = fn[: fn.rfind('.acc')]
                # list_exp.append(os.path.basename(fn))
                bn = os.path.basename(fn).split('-')
                indiv = bn[0]
                method = '-'.join(bn[1:])
                list_indiv.append(indiv)
                list_method.append(method)
        # df['Experiments'] = list_exp
        df['Inidvidual'] = list_indiv
        df['Super Population'] = [dict_pop_to_spop[dict_indiv_to_pop[s]] for s in list_indiv]
        df['Method'] = list_method
        df['True Positive'] = list_tp
        df['Num Reads'] = list_all
        df.to_csv(output.tsv, sep='\t', index=None)
