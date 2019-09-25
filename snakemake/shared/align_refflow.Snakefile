'''
Rules for two-pass refflow
'''
rule refflow_separate_fistpass_results:
    input:
        sam = os.path.join(DIR_FIRST_PASS, CHROM + '-h37maj.sam')
    output:
        highq = os.path.join(DIR_FIRST_PASS, CHROM + 
            '-h37maj-mapqgeq' + ALN_MAPQ_THRSD + '.sam'),
        lowq = os.path.join(DIR_FIRST_PASS, CHROM + 
            '-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.sam'),
        lowq_reads = os.path.join(DIR_FIRST_PASS, CHROM +
            '-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.fq')
    shell:
        'awk -v var="{ALN_MAPQ_THRSD}" \
        \'{{ if ($5 >= var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.highq};'
        'awk -v var={ALN_MAPQ_THRSD} \
        \'{{ if ($5 < var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.lowq};'
        'samtools fastq {output.lowq} > {output.lowq_reads}'

rule refflow_align_secondpass:
    input:
        reads = os.path.join(DIR_FIRST_PASS, CHROM +
            '-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.fq'),
        idx1 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.1.bt2',
        idx2 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.2.bt2',
        idx3 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.3.bt2',
        idx4 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.4.bt2',
        idx5 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.1.bt2',
        idx6 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.2.bt2'
    params:
        index = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX
    output:
        # sam = os.path.join(DIR_SECOND_PASS, CHROM + '-h37maj-' + ALN_MAPQ_THRSD + '-{GROUP}-' + POP_DIRNAME + '.sam')
        sam = PREFIX_SECOND_PASS + '.sam'
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {THREADS} -x {params.index} -U {input.reads} -S {output.sam};'

rule refflow_merge_secondpass:
    input:
        sam = expand(
            PREFIX_SECOND_PASS + '.sam',
            INDIV = INDIV, GROUP = GROUP)
    output:
        path = os.path.join(DIR_SECOND_PASS, '{0}-h37maj-{1}-{2}.paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        id = os.path.join(DIR_SECOND_PASS, '{0}-h37maj-{1}-{2}.ids'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        merge_paths = os.path.join(DIR_SECOND_PASS, '{0}-h37maj-{1}-{2}.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    params:
        prefix = os.path.join(DIR_SECOND_PASS, '2ndpass')
    run:
        indiv = wildcards.INDIV
        dir_2p = os.path.join(DIR, 'experiments/' + indiv + '/' + POP_DIRNAME)
        maj_lowq = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/' + CHROM +
            '-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.sam')

        shell('echo {maj_lowq} > {output.path};')
        shell('echo "h37maj" > {output.id};')
        for g in GROUP:
            fn = os.path.join(dir_2p, '{0}-h37maj-{1}-{2}-{3}.sam'.format(CHROM, ALN_MAPQ_THRSD, g, POP_DIRNAME))
            shell('ls {fn} >> {output.path};')
            shell('echo {g} >> {output.id};')
        shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
            -ids {output.id} -rs {RAND_SEED} -p {params.prefix} \
            -l {output.merge_paths};')
        # shell('{PYTHON} {DIR_SCRIPTS}/merge_sam.py \
        #     -ns {output.path} -ids {output.id} \
        #     -rs 0 -l {output.merge_paths}')
        
rule check_secondpass:
    input:
        sam = expand(
            os.path.join(DIR_SECOND_PASS, '{0}-h37maj-{1}-{2}.paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV, GROUP = GROUP
        ),
        merge = expand(
            os.path.join(DIR_SECOND_PASS, '{0}-h37maj-{1}-{2}.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV, GROUP = GROUP
        )
    output:
        touch(temp(os.path.join(DIR, 'refflow_secondpass.done')))

''' Old merge methods '''
# rule refflow_align_and_merge_secondpass:
#     input:
#         reads = os.path.join(DIR_FIRST_PASS, CHROM +
#             '-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.fq'),
#         idx = expand(
#             DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.{idx}.bt2',
#             idx = IDX_ITEMS, GROUP = GROUP)
#     params:
#         prefix = os.path.join(DIR_SECOND_PASS, '2ndpass')
#     #     index = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX
#     output:
#         merge_log = os.path.join(DIR_SECOND_PASS, 'merged_log'),
#         merge_id = os.path.join(DIR_SECOND_PASS, 'merge_id'),
#         merge_path = os.path.join(DIR_SECOND_PASS, 'merge_path')
#     run:
#         maj_lowq = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/' + CHROM +
#             '-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.sam')
#         shell('echo {maj_lowq} > {output.merge_path};')
#         shell('echo "h37maj" > {output.merge_id};')
#         for g in GROUP:
#             shell('echo {g} >> {output.merge_id};')
#             index = DIR_POP_GENOME_BLOCK_IDX + \
#                 CHROM + '_superpop_' + g + \
#                 '_thrds{}_S{}_b{}_ld{}'.format(POP_THRSD, POP_STOCHASTIC, POP_BLOCK_SIZE, POP_USE_LD)
#             shell(
#             'bowtie2 --reorder --threads {THREADS} -x {index} \
#             -U {input.reads} | \
#             {PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.merge_path} \
#             -ids {output.merge_id} -rs {RAND_SEED} -p {params.prefix} \
#             -l {output.merge_log} -id {g};')
#             shell('tail -1 {output.merge_log} >> {output.merge_path};')
