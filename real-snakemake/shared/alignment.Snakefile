''' Reference flow: first pass '''
rule align_to_ref:
    input:
        reads1 = READS1,
        idx = expand(
            os.path.join(DIR, 'grch37/wg.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        index = os.path.join(DIR, 'grch37/indexes/wg')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-GRCh37.sam')
    threads: THREADS
    shell:
        'bowtie2 --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam}'

rule align_to_major:
    input:
        reads1 = READS1,
        idx = expand(
            os.path.join(DIR, 'major/wg_h37maj.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        index = os.path.join(DIR, 'major/indexes/wg_h37maj')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-h37maj.sam')
    threads: THREADS
    shell:
        'bowtie2 --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam}'

DIR_PER_IDX = os.path.join(DIR, 'personalized/{INDIV}/indexes/')
''' Personalized '''
rule align_to_perA:
    input:
        reads1 = READS1,
        idxA = expand(
            os.path.join(DIR_PER_IDX, 'wg_perA.{IDX_ITEMS}.bt2'),
            IDX_ITEMS = IDX_ITEMS, INDIV = INDIV)
    params:
        indexA = os.path.join(DIR_PER_IDX, 'wg_perA')
    output:
        samA = os.path.join(DIR_FIRST_PASS, 'wg-per_hapA.sam')
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {THREADS} -x {params.indexA} -U {input.reads1} -S {output.samA};'

rule align_to_perB:
    input:
        reads1 = READS1,
        idxB = expand(
            os.path.join(DIR_PER_IDX, 'wg_perB.{IDX_ITEMS}.bt2'),
            IDX_ITEMS = IDX_ITEMS, INDIV = INDIV)
    params:
        indexB = os.path.join(DIR_PER_IDX, 'wg_perB')
    output:
        samB = os.path.join(DIR_FIRST_PASS, 'wg-per_hapB.sam')
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {THREADS} -x {params.indexB} -U {input.reads1} -S {output.samB};'

rule merge_per:
    input:
        samA = os.path.join(DIR_FIRST_PASS, 'wg-per_hapA.sam'),
        samB = os.path.join(DIR_FIRST_PASS, 'wg-per_hapB.sam')
    output:
        path = os.path.join(DIR_FIRST_PASS, 'wg-per.paths'),
        id = os.path.join(DIR_FIRST_PASS, 'wg-per.ids'),
        merge_paths = os.path.join(DIR_FIRST_PASS, 'wg-per.merge_paths'),
        samA = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapA.sam'),
        samB = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapB.sam')
    params:
        os.path.join(DIR_FIRST_PASS, 'wg-per-merged')
    run:
        #: prepare ids
        for h in ['hapA', 'hapB']:
            shell('echo {h} >> {output.id};')
        #: prepare paths
        shell('ls {input.samA} >> {output.path};')
        shell('ls {input.samB} >> {output.path};')
        #: merge_incremental
        shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
            -ids {output.id} -rs {RAND_SEED} -p {params} \
            -l {output.merge_paths};')

''' Refflow '''
rule refflow_separate_fistpass_results:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-h37maj.sam')
    output:
        highq = os.path.join(DIR_FIRST_PASS,
            'wg-h37maj-mapqgeq' + ALN_MAPQ_THRSD + '.sam'),
        lowq = os.path.join(DIR_FIRST_PASS,
            'wg-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.sam'),
        lowq_reads = os.path.join(DIR_FIRST_PASS,
            'wg-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.fq')
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
        reads = os.path.join(DIR_FIRST_PASS,
            'wg-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.fq'),
        idx1 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.1.bt2'),
        idx2 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.2.bt2'),
        idx3 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.3.bt2'),
        idx4 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.4.bt2'),
        idx5 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.1.bt2'),
        idx6 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.2.bt2')
    params:
        index = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX)
    output:
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
        path = os.path.join(DIR_SECOND_PASS, 'wg-h37maj-{}-{}.paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        id = os.path.join(DIR_SECOND_PASS, 'wg-h37maj-{}-{}.ids'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        merge_paths = os.path.join(DIR_SECOND_PASS, 'wg-h37maj-{}-{}.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    params:
        prefix = os.path.join(DIR_SECOND_PASS, '2ndpass')
    run:
        # indiv = wildcards.INDIV
        dir_2p = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME)
        maj_lowq = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/wg-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.sam')

        shell('echo {maj_lowq} > {output.path};')
        shell('echo "h37maj" > {output.id};')
        for g in GROUP:
            fn = os.path.join(dir_2p, 'wg-h37maj-{}-{}-{}.sam'.format(ALN_MAPQ_THRSD, g, POP_DIRNAME))
            shell('ls {fn} >> {output.path};')
            shell('echo {g} >> {output.id};')
        shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
            -ids {output.id} -rs {RAND_SEED} -p {params.prefix} \
            -l {output.merge_paths};')

rule check_alignment_standard:
    input:
        grc = expand(
            os.path.join(DIR_FIRST_PASS, 'wg-GRCh37.sam'),
            INDIV = INDIV),
        major = expand(
            os.path.join(DIR_FIRST_PASS, 'wg-h37maj.sam'),
            INDIV = INDIV),
        perA = expand(
            os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapA.sam'),
            # os.path.join(DIR_FIRST_PASS, 'wg-per_hapA.sam'),
            INDIV = INDIV),
        perB = expand(
            os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapB.sam'),
            # os.path.join(DIR_FIRST_PASS, 'wg-per_hapB.sam'),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'alignment_standard.done')))

rule check_alignment_refflow:
    input:
        merge_paths = expand(
            os.path.join(DIR_SECOND_PASS, 'wg-h37maj-{}-{}.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'alignment_refflow.done')))

# rule check_secondpass:
#     input:
#         merge_paths = expand(
#             os.path.join(DIR_SECOND_PASS, 'wg-h37maj-{}-{}.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
#             INDIV = INDIV)
#     output:
#         touch(temp(os.path.join(DIR, 'secondpass.done')))
