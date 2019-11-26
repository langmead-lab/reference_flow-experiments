''' Align to refflow genomes in one-pass '''
# rule refflow_align_onepass:
#     input:
#         reads = PREFIX_PER + '_1.fq',
#         idx1 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.1.bt2',
#         idx2 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.2.bt2',
#         idx3 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.3.bt2',
#         idx4 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.4.bt2',
#         idx5 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.1.bt2',
#         idx6 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.2.bt2'
#     params:
#         index = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX
#     output:
#         sam = os.path.join(DIR_FIRST_PASS, CHROM + '-{GROUP}-' + POP_DIRNAME +'.sam')
#     threads: THREADS
#     shell:
#         'bowtie2 --threads {threads} -x {params.index} -U {input.reads} -S {output.sam};'

'''
Rules for two-pass refflow
'''
rule refflow_separate_fistpass_results_gnomad:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad.sam'.format(CHROM))
    output:
        highq = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-gnomad-mapqgeq{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lowq = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-gnomad-mapqlt{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lowq_reads = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-gnomad-mapqlt{}.fq'.format(CHROM, ALN_MAPQ_THRSD))
    shell:
        'awk -v var="{ALN_MAPQ_THRSD}" \
        \'{{ if ($5 >= var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.highq};'
        'awk -v var={ALN_MAPQ_THRSD} \
        \'{{ if ($5 < var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.lowq};'
        'samtools fastq {output.lowq} > {output.lowq_reads}'

rule refflow_separate_fistpass_results_onekg:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg.sam'.format(CHROM))
    output:
        highq = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-1kg-mapqgeq{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lowq = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-1kg-mapqlt{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lowq_reads = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-1kg-mapqlt{}.fq'.format(CHROM, ALN_MAPQ_THRSD))
    shell:
        'awk -v var="{ALN_MAPQ_THRSD}" \
        \'{{ if ($5 >= var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.highq};'
        'awk -v var={ALN_MAPQ_THRSD} \
        \'{{ if ($5 < var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.lowq};'
        'samtools fastq {output.lowq} > {output.lowq_reads}'

rule refflow_align_secondpass_gnomad:
    input:
        reads = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-gnomad-mapqlt{}.fq'.format(CHROM, ALN_MAPQ_THRSD)),
        idx1 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.1.bt2',
        idx2 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.2.bt2',
        idx3 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.3.bt2',
        idx4 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.4.bt2',
        idx5 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.1.bt2',
        idx6 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.2.bt2'
    params:
        index = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX
    output:
        sam = PREFIX_SECOND_PASS + '-gnomad.sam'
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {threads} -x {params.index} -U {input.reads} -S {output.sam};'

rule refflow_align_secondpass_onekg:
    input:
        reads = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-1kg-mapqlt{}.fq'.format(CHROM, ALN_MAPQ_THRSD)),
        idx1 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.1.bt2',
        idx2 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.2.bt2',
        idx3 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.3.bt2',
        idx4 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.4.bt2',
        idx5 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.1.bt2',
        idx6 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.2.bt2'
    params:
        index = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX
    output:
        sam = PREFIX_SECOND_PASS + '-1kg.sam'
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {threads} -x {params.index} -U {input.reads} -S {output.sam};'

rule refflow_merge_secondpass_gnomad:
    input:
        sam = expand(
            PREFIX_SECOND_PASS + '-gnomad.sam',
            INDIV = INDIV, GROUP = GROUP),
        lowq = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-gnomad-mapqlt{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
    output:
        path = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-gnomad.paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        id = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-gnomad.ids'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        merge_paths = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-gnomad.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    params:
        prefix = os.path.join(DIR_SECOND_PASS, '2ndpass-gnomad')
    run:
        indiv = wildcards.INDIV
        dir_2p = os.path.join(DIR, 'experiments/' + indiv + '/' + POP_DIRNAME)
        # maj_lowq = os.path.join(DIR, 'experiments/' + indiv + '/' + CHROM +
        #     '-major-mapqlt' + ALN_MAPQ_THRSD + '.sam')

        shell('echo {input.lowq} > {output.path};')
        # shell('echo {maj_lowq} > {output.path};')
        shell('echo "major" > {output.id};')
        for g in GROUP:
            fn = os.path.join(dir_2p, 'chr{0}-major-{1}-{2}-{3}-gnomad.sam'.format(CHROM, ALN_MAPQ_THRSD, g, POP_DIRNAME))
            shell('ls {fn} >> {output.path};')
            shell('echo {g} >> {output.id};')
        shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
            -ids {output.id} -rs {RAND_SEED} -p {params.prefix} \
            -l {output.merge_paths};')

rule refflow_merge_secondpass_onekg:
    input:
        sam = expand(
            PREFIX_SECOND_PASS + '-1kg.sam',
            INDIV = INDIV, GROUP = GROUP),
        lowq = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-1kg-mapqlt{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
    output:
        path = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-1kg.paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        id = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-1kg.ids'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        merge_paths = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-1kg.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    params:
        prefix = os.path.join(DIR_SECOND_PASS, '2ndpass-1kg')
    run:
        indiv = wildcards.INDIV
        dir_2p = os.path.join(DIR, 'experiments/' + indiv + '/' + POP_DIRNAME)
        # maj_lowq = os.path.join(DIR, 'experiments/' + indiv + '/' + CHROM +
        #     '-major-mapqlt' + ALN_MAPQ_THRSD + '.sam')

        shell('echo {input.lowq} > {output.path};')
        # shell('echo {maj_lowq} > {output.path};')
        shell('echo "major" > {output.id};')
        for g in GROUP:
            fn = os.path.join(dir_2p, 'chr{0}-major-{1}-{2}-{3}-1kg.sam'.format(CHROM, ALN_MAPQ_THRSD, g, POP_DIRNAME))
            shell('ls {fn} >> {output.path};')
            shell('echo {g} >> {output.id};')
        shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
            -ids {output.id} -rs {RAND_SEED} -p {params.prefix} \
            -l {output.merge_paths};')

rule check_secondpass:
    input:
        # sam = expand(
        #     os.path.join(DIR_SECOND_PASS, '{0}-major-{1}-{2}.paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        #     INDIV = INDIV, GROUP = GROUP
        # ),
        # gnomad = expand(
        #     os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-gnomad.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        #     INDIV = INDIV, GROUP = GROUP
        # ),
        onekg = expand(
            os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-1kg.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV, GROUP = GROUP
        )
    output:
        touch(temp(os.path.join(DIR, 'refflow_secondpass.done')))
