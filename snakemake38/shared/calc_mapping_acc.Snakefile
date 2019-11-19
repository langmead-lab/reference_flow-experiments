'''
Calculate mapping accuracy for all settings
'''
rule calc_grc_accuracy:
    input:
        sam = os.path.join(DIR_FIRST_PASS, CHROM + '-grc.sam'),
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, CHROM + '-grc.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-grc.acc')
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} > {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule calc_major_accuracy_gnomad:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-liftover-sorted.sam'.format(CHROM)),
        # sam = os.path.join(DIR_FIRST_PASS, CHROM + '-major.sam'),
        # var_genome = PREFIX_MAJOR + '.var',
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad.acc_log'.format(CHROM)),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + 'chr{}-major-gnomad.acc'.format(CHROM))
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} > {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule calc_major_accuracy_onekg:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-liftover-sorted.sam'.format(CHROM)),
        # sam = os.path.join(DIR_FIRST_PASS, CHROM + '-major.sam'),
        # var_genome = PREFIX_MAJOR + '.var',
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg.acc_log'.format(CHROM)),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + 'chr{}-major-1kg.acc'.format(CHROM))
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} > {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule calc_per_accuracy:
    input:
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
        sam = os.path.join(DIR_FIRST_PASS,
            '{}-per-merged-liftover-sorted.sam'.format(CHROM))
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, CHROM + '-per.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-per.acc')
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} >> \
        {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule calc_per_haptohap_accuracy:
    input:
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
        sam = os.path.join(DIR_FIRST_PASS,
            '{}-per_h2h-merged-liftover-sorted.sam'.format(CHROM))
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, CHROM + '-per_h2h.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-per_h2h.acc')
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} >> \
        {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule calc_refflow_accuracy_gnomad:
    input:
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
        sam = os.path.join(DIR_SECOND_PASS, 'chr{}-refflow-{}-{}-liftover-gnomad-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, 'chr{}-refflow-{}-{}-gnomad.acc_log'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        acc = os.path.join(DIR_RESULTS,
            '{INDIV}' + '-chr{0}-major-{1}-{2}-gnomad.acc'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} >> \
        {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule calc_refflow_accuracy_onekg:
    input:
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
        sam = os.path.join(DIR_SECOND_PASS, 'chr{}-refflow-{}-{}-liftover-1kg-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, 'chr{}-refflow-{}-{}-1kg.acc_log'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        acc = os.path.join(DIR_RESULTS,
            '{INDIV}' + '-chr{0}-major-{1}-{2}-1kg.acc'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} >> \
        {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule check_standard_accuracy:
    input:
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + 'chr{}-major-gnomad.acc'.format(CHROM)),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + 'chr{}-major-1kg.acc'.format(CHROM)),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-grc.acc'),
            INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-per.acc'),
        #     INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-per_h2h.acc'),
        #     INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'standard_acc.done')))

rule check_refflow_accuracy:
    input:
        gnomad = expand(os.path.join(DIR_RESULTS,
            '{INDIV}' + '-chr{0}-major-{1}-{2}-gnomad.acc'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV),
        onekg = expand(os.path.join(DIR_RESULTS,
            '{INDIV}' + '-chr{0}-major-{1}-{2}-1kg.acc'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV),
        # onepass = expand(
        #     os.path.join(DIR_RESULTS,
        #     '{INDIV}-' + CHROM + '-{GROUP}-' + POP_GENOME_SUFFIX + '.acc'),
        #     INDIV = INDIV, GROUP = GROUP)
    output:
        touch(temp(os.path.join(DIR, 'refflow_acc.done')))

'''
Summarize results as a TSV
'''
rule check_mapping_acc_and_write_as_tsv:
    input:
        acc_standard = os.path.join(DIR, 'standard_acc.done'),
        acc_refflow  = os.path.join(DIR, 'refflow_acc.done'),
        # acc_vg = os.path.join(DIR, 'vg_acc.done')
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
        df.to_csv(output.tsv, sep='\t', index=None, float_format = '%.4f')
