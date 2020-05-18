'''
Calculate mapping accuracy for all settings
'''
CALC_ACC_MEM_USAGE_GB = 100
rule merge_gold_sam:
    input:
        first = PREFIX_PER + '_1.sam',
        second = PREFIX_PER + '_2.sam'
    output:
        PREFIX_PER + '.sam'
    shell:
        'cat {input.first} > {output};'
        'cat {input.second} >> {output}'

rule calc_grc_accuracy:
    input:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC.sam'),
        gold = PREFIX_PER + '.sam',
        var_reads = PREFIX_PER + '-per.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-grc.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + EXP_LABEL + '-grc.acc')
    resources:
        mem_gb = CALC_ACC_MEM_USAGE_GB
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS_EXP}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} > {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule get_major_sam_sorted:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted.bam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted.sam')
    shell:
        '{SAMTOOLS} view -h -o {output} {input}'

rule calc_major_accuracy:
    input:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted.sam'),
        gold = PREFIX_PER + '.sam',
        var_reads = PREFIX_PER + '-per.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + EXP_LABEL + '-major.acc')
    resources:
        mem_gb = CALC_ACC_MEM_USAGE_GB
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS_EXP}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} > {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule get_per_sam_sorted:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted.bam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted.sam')
    shell:
        '{SAMTOOLS} view -h -o {output} {input}'

rule calc_per_accuracy:
    input:
        gold = PREFIX_PER + '.sam',
        var_reads = PREFIX_PER + '-per.var',
        sam = os.path.join(DIR_FIRST_PASS,
            EXP_LABEL + '-per-merged-liftover-sorted.sam')
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + EXP_LABEL + '-per.acc')
    resources:
        mem_gb = CALC_ACC_MEM_USAGE_GB
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS_EXP}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} >> \
        {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule calc_refflow_accuracy:
    input:
        gold = PREFIX_PER + '.sam',
        var_reads = PREFIX_PER + '-per.var',
        sam = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL +
            '-refflow-{}-{}-liftover.sam'.format(ALN_MAPQ_THRSD, POP_DIRNAME)
            )
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-refflow-{}-{}.acc_log'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        acc = os.path.join(DIR_RESULTS,
            '{INDIV}-' + EXP_LABEL + '-major-{}-{}.acc'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    resources:
        mem_gb = CALC_ACC_MEM_USAGE_GB
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS_EXP}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} >> \
        {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule calc_popspecific_onepass_accuarcy:
    input:
        gold = PREFIX_PER + '.sam',
        var_reads = PREFIX_PER + '-per.var',
        sam = os.path.join(DIR_FIRST_PASS,
            EXP_LABEL + '-{GROUP}-' + POP_DIRNAME +'-liftover-sorted.sam')
    output:
        acc_log = os.path.join(DIR_FIRST_PASS,
            EXP_LABEL + '-{GROUP}-' + POP_DIRNAME +'-liftover-sorted.acc_log'),
        acc = os.path.join(DIR_RESULTS,
            '{INDIV}-' + EXP_LABEL + '-{GROUP}-' + POP_DIRNAME + '.acc'),
    resources:
        mem_gb = CALC_ACC_MEM_USAGE_GB
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS_EXP}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} >> \
        {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule get_vg_sam_sorted:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}-sorted.bam'.format(ALLELE_FREQ_FOR_VG))
    output:
        temp(os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}-sorted.sam'.format(ALLELE_FREQ_FOR_VG)))
    shell:
        '{SAMTOOLS} view -h -o {output} {input}'

rule calc_vg_accuracy:
    input:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}-sorted.sam'.format(ALLELE_FREQ_FOR_VG)),
        gold = PREFIX_PER + '.sam',
        var_reads = PREFIX_PER + '-per.var'
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}.acc_log'.format(ALLELE_FREQ_FOR_VG)),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + EXP_LABEL + '-vg_{}.acc'.format(ALLELE_FREQ_FOR_VG))
    resources:
        mem_gb = CALC_ACC_MEM_USAGE_GB
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS_EXP}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} >> \
        {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule check_vg_accuracy:
    input:
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + EXP_LABEL + '-vg_{}.acc'.format(ALLELE_FREQ_FOR_VG)),
            INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'vg_acc.done')))

rule check_standard_accuracy:
    input:
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + EXP_LABEL + '-major.acc'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + EXP_LABEL + '-grc.acc'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + EXP_LABEL + '-per.acc'),
            INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'standard_acc.done')))

rule check_refflow_accuracy:
    input:
        expand(os.path.join(DIR_RESULTS,
            '{INDIV}-' + EXP_LABEL + '-major-{}-{}.acc'.format(
                ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV),
#         expand(
#             os.path.join(DIR_RESULTS,
#             '{INDIV}-' + EXP_LABEL + '-{GROUP}-' + POP_DIRNAME + '.acc'),
#             INDIV = INDIV, GROUP = GROUP)
    output:
        touch(temp(os.path.join(DIR, 'refflow_acc.done')))

'''
Summarize results as a TSV
'''
rule check_mapping_acc_and_write_as_tsv:
    input:
        acc_standard = os.path.join(DIR, 'standard_acc.done'),
        acc_refflow  = os.path.join(DIR, 'refflow_acc.done'),
        acc_vg = os.path.join(DIR, 'vg_acc.done')
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
