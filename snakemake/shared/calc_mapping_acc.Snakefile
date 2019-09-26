'''
Calculate mapping accuracy for all settings
'''
rule calc_per_accuracy:
    input:
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
        merge_paths = os.path.join(DIR_FIRST_PASS, '{}-per.merge_paths'.format(CHROM))
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, CHROM + '-per.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-per.acc')
    run:
        with open(input.merge_paths, 'r') as f:
            for line in f:
                fn = line.rstrip()
                shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
                -c {CHROM} -g {input.gold} -p 2 -vr {input.var_reads} \
                -vs {input.var_reads} -n {fn} >> {output.acc_log}')
        organize_accuracy(output.acc_log, output.acc)

#: ---- OLD MERGING
# rule personalize_merge_and_calc_accuracy:
#     input:
#         samA = DIR_FIRST_PASS + CHROM + '-per_hapA.sam',
#         samB = DIR_FIRST_PASS + CHROM + '-per_hapB.sam',
#         gold = PREFIX_PER + '_1.sam',
#         var_reads = PREFIX_PER + '.var',
#     output:
#         merge_path = os.path.join(DIR_FIRST_PASS, CHROM + '-per-merged.paths'),
#         acc_log = os.path.join(DIR_FIRST_PASS, CHROM + '-per.acc_log'),
#         acc = os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-per.acc')
#     run:
#         shell('{PYTHON} {DIR_SCRIPTS}/merge_sam.py \
#             -n1 {input.samA} -id1 hapA \
#             -n2 {input.samB} -id2 hapB \
#             -rs 0 -l {output.merge_path};')
#         with open(output.merge_path, 'r') as f:
#             for line in f:
#                 fn = line.rstrip()
#                 shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
#                 -c {CHROM} -g {input.gold} -p 2 -vr {input.var_reads} \
#                 -vs {input.var_reads} -n {fn} >> {output.acc_log}')
#         organize_accuracy(output.acc_log, output.acc)
#: ----

rule calc_grc_accuracy:
    input:
        sam = os.path.join(DIR_FIRST_PASS, CHROM + '-grch37.sam'),
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, CHROM + '-grch37.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-grch37.acc')
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} > {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule calc_major_accuracy:
    input:
        sam = os.path.join(DIR_FIRST_PASS, CHROM + '-h37maj.sam'),
        var_genome = PREFIX_MAJOR + '.var',
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, CHROM + '-h37maj.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-h37maj.acc')
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -vs {input.var_genome} -n {input.sam} > {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule check_standard_accuracy:
    input:
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-h37maj.acc'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-per.acc'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-grch37.acc'),
            INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'standard_acc.done')))

rule calc_refflow_firstpass_accuracy:
    input:
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
        first_sam = os.path.join(DIR_FIRST_PASS, CHROM + 
            '-h37maj-mapqgeq' + ALN_MAPQ_THRSD + '.sam'),
        first_var_genome = PREFIX_MAJOR + '.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, '{0}-{1}-{2}-firstpass.acc_log'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    shell:
        '{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -vs {input.first_var_genome} -n {input.first_sam} >> \
        {output.acc_log}'

rule calc_refflow_accuracy:
    input:
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
        # first_sam = os.path.join(DIR_FIRST_PASS, CHROM + 
        #     '-h37maj-mapqgeq' + ALN_MAPQ_THRSD + '.sam'),
        first_var_genome = PREFIX_MAJOR + '.var',
        acc_log = os.path.join(DIR_FIRST_PASS, '{0}-{1}-{2}-firstpass.acc_log'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        second_sam_path = os.path.join(DIR_SECOND_PASS, '{0}-h37maj-{1}-{2}.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        second_group = os.path.join(DIR_SECOND_PASS, '{0}-h37maj-{1}-{2}.ids'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        second_var_genome = expand(
            os.path.join(
            DIR_POP_GENOME_BLOCK,
            CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.var'),
            # CHROM + '_superpop_{GROUP}_thrds' + str(POP_THRSD) +
            # '_S' + str(POP_STOCHASTIC) + '_b' + str(POP_BLOCK_SIZE) + 
            # '_ld' + str(POP_USE_LD) + '.var'),
            INDIV = INDIV, GROUP = GROUP
        )
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, '{0}-{1}-{2}.acc_log'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    run:
        # if os.path.exists(output.acc_log):
        #     shell('rm {output.acc_log};')
        # shell(
        #     '{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        #     -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        #     -vs {input.first_var_genome} -n {input.first_sam} >> \
        #     {output.acc_log};')
        shell('cp {input.acc_log} {output.acc_log};')
        list_group = []
        with open(input.second_group) as f:
            for line in f:
                list_group.append(line.rstrip())
        list_fn = []
        with open(input.second_sam_path) as f:
            for line in f:
                list_fn.append(line.rstrip())
        for i, fn in enumerate(list_fn):
            if i == 0:
                second_var_genome = input.first_var_genome
            else:
                second_var_genome = os.path.join(
                    DIR_POP_GENOME_BLOCK,
                    CHROM + '_superpop_' + list_group[i] + '_' + POP_DIRNAME + '.var'
                )
            shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
            -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
            -vs {second_var_genome} -n {fn} >> \
            {output.acc_log};')

rule sum_refflow_accuracy:
    input:
        acc_log = os.path.join(DIR_FIRST_PASS, '{0}-{1}-{2}.acc_log'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        acc = os.path.join(DIR_RESULTS,
            '{INDIV}' + '-{0}-h37maj-{1}-{2}.acc'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    run:
        organize_accuracy(input.acc_log, output.acc)

'''
Summarize results as a TSV
'''
rule check_mapping_acc_and_write_as_tsv:
    input:
        acc = expand(os.path.join(DIR_RESULTS,
            '{INDIV}' + '-{0}-h37maj-{1}-{2}.acc'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV),
        acc_standard = os.path.join(DIR, 'standard_acc.done')
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
