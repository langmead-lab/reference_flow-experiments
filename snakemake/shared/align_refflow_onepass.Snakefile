rule refflow_onepass_calc_accuracy:
    input:
        sam = os.path.join(DIR_FIRST_PASS, CHROM + '-{GROUP}-' + POP_DIRNAME +'.sam'),
        var_genome = os.path.join(
            DIR_POP_GENOME_BLOCK,
            CHROM + '_superpop_{GROUP}_thrds' + str(POP_THRSD) +
            '_S' + str(POP_STOCHASTIC) + '_b' + str(POP_BLOCK_SIZE) + 
            '_ld' + str(POP_USE_LD) + '.var'
        ),
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, CHROM + '-{GROUP}-' + POP_GENOME_SUFFIX + '.acc_log'),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-{GROUP}-' + POP_GENOME_SUFFIX + '.acc')
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -vs {input.var_genome} -n {input.sam} > {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule check_refflow_onepass_accuracy:
    input:
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + CHROM + '-{GROUP}-' + POP_GENOME_SUFFIX + '.acc'),
            INDIV = INDIV, GROUP = GROUP)
    output:
        touch(temp(os.path.join(DIR, 'rf_onepass_acc.done')))

