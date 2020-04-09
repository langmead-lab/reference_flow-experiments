''' Reference flow: first pass '''
rule align_to_ref:
    input:
        reads1 = READS1,
        reads2 = READS2,
        idx = expand(
            os.path.join(DIR, 'grc/wg.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        index = os.path.join(DIR, 'grc/wg')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-GRC.sam')
    threads: THREADS
    shell:
        'bowtie2 --threads {THREADS} -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam}'
        # 'bowtie2 --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam}'

DIR_PER_IDX = os.path.join(DIR, 'personalized/{INDIV}/indexes/')
''' Personalized '''
rule align_to_perA:
    input:
        reads1 = READS1,
        reads2 = READS2,
        idx = expand(
            os.path.join(DIR_PER_IDX, 'wg-perA.{IDX_ITEMS}.bt2'),
            IDX_ITEMS = IDX_ITEMS, INDIV = INDIV)
    params:
        index = os.path.join(DIR_PER_IDX, 'wg-perA')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-per_hapA.sam')
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {THREADS} -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam}'
        # 'bowtie2 --reorder --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam}'

rule align_to_perB:
    input:
        reads1 = READS1,
        reads2 = READS2,
        idx = expand(
            os.path.join(DIR_PER_IDX, 'wg-perB.{IDX_ITEMS}.bt2'),
            IDX_ITEMS = IDX_ITEMS, INDIV = INDIV)
    params:
        index = os.path.join(DIR_PER_IDX, 'wg-perB')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-per_hapB.sam')
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {THREADS} -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam};'
        # 'bowtie2 --reorder --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam};'

rule merge_per:
    input:
        samA = os.path.join(DIR_FIRST_PASS, 'wg-per_hapA.sam'),
        samB = os.path.join(DIR_FIRST_PASS, 'wg-per_hapB.sam')
    output:
        path = os.path.join(DIR_FIRST_PASS, 'wg-per.paths'),
        label = os.path.join(DIR_FIRST_PASS, 'wg-per.ids'),
        merge_paths = os.path.join(DIR_FIRST_PASS, 'wg-per.merge_paths'),
        samA = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapA.sam'),
        samB = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapB.sam')
    params:
        os.path.join(DIR_FIRST_PASS, 'wg-per-merged')
    run:
        #: prepare ids
        for h in ['hapA', 'hapB']:
            shell('echo {h} >> {output.label};')
        #: prepare paths
        shell('ls {input.samA} >> {output.path};')
        shell('ls {input.samB} >> {output.path};')
        #: merge_incremental
        shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
            -ids {output.label} -rs {RAND_SEED} -p {params} \
            -l {output.merge_paths} --paired-end')

rule check_alignment_grc_and_per:
    input:
        grc = expand(
            os.path.join(DIR_FIRST_PASS, 'wg-GRC.sam'),
            INDIV = INDIV),
        perA = expand(
            os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapA.sam'),
            INDIV = INDIV),
        perB = expand(
            os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapB.sam'),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'alignment_grc_per.done')))
