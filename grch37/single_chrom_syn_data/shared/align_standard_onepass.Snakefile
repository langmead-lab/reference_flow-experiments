'''
Rules for one-pass alignment methods including grch37, major and personalized.

Grch37 and major are straightforward one-pass alignments, but the personalized
alignment method all are aligned to both hapA and hapB, and a merging step
is performed to generate the "best" outputs.

Checkpoint:
    temp(os.path.join(DIR, 'standard_onepass.done'))
'''
rule align_to_major:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        idx = expand(PREFIX_MAJOR_IDX + '.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS)
    params:
        PREFIX_MAJOR_IDX
    output:
        sam = os.path.join(DIR_FIRST_PASS, CHROM + '-h37maj.sam')
    threads: THREADS
    shell:
        'bowtie2 --threads {THREADS} -x {params} -U {input.reads1} -S {output.sam}'

rule align_to_grc:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        idx = expand(os.path.join(
            DIR_GRCH37_IDX, CHROM + '_grch37.{i}.bt2'), i = IDX_ITEMS)
    params:
        DIR_GRCH37_IDX + CHROM + '_grch37'
    output:
        sam = os.path.join(DIR_FIRST_PASS, CHROM + '-grch37.sam')
    threads: THREADS
    shell:
        'bowtie2 --threads {THREADS} -x {params} -U {input.reads1} -S {output.sam}'

# rule check_firstpass:
#     input:
#         sam = expand(
#             os.path.join(DIR_FIRST_PASS, CHROM + '-h37maj.sam'),
#             INDIV = INDIV),
#         highq = expand(
#             os.path.join(DIR_FIRST_PASS, CHROM + 
#             '-h37maj-mapqgeq' + ALN_MAPQ_THRSD + '.sam'),
#             INDIV = INDIV),
#         lowq_reads = expand(
#             os.path.join(DIR_FIRST_PASS, CHROM +
#             '-h37maj-mapqlt' + ALN_MAPQ_THRSD + '.fq'),
#             INDIV = INDIV)
#     output:
#         touch(temp(os.path.join(DIR, 'firstpass.done')))

rule align_to_per_haploid_setting:
    input:
        readsA1 = PREFIX_PER + '_hapA_1.fq',
        readsB1 = PREFIX_PER + '_hapB_1.fq',
        idxA = expand(DIR_PER_IDX + CHROM + '-per_hapA.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS, INDIV = INDIV),
        idxB = expand(DIR_PER_IDX + CHROM + '-per_hapB.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS, INDIV = INDIV)
    params:
        indexA = DIR_PER_IDX + CHROM + '-per_hapA',
        indexB = DIR_PER_IDX + CHROM + '-per_hapB'
    output:
        samA = os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapA_haploid.sam'),
        samB = os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapB_haploid.sam')
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {THREADS} -x {params.indexA} -U {input.readsA1} -S {output.samA};'
        'bowtie2 --reorder --threads {THREADS} -x {params.indexB} -U {input.readsB1} -S {output.samB}'

rule align_to_per:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        idxA = expand(DIR_PER_IDX + CHROM + '-per_hapA.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS, INDIV = INDIV),
        idxB = expand(DIR_PER_IDX + CHROM + '-per_hapB.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS, INDIV = INDIV)
    params:
        indexA = DIR_PER_IDX + CHROM + '-per_hapA',
        indexB = DIR_PER_IDX + CHROM + '-per_hapB'
    output:
        samA = os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapA.sam'),
        samB = os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapB.sam')
        # samA = DIR_FIRST_PASS + CHROM + '-per_hapA.sam',
        # samB = DIR_FIRST_PASS + CHROM + '-per_hapB.sam'
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {THREADS} -x {params.indexA} -U {input.reads1} -S {output.samA};'
        'bowtie2 --reorder --threads {THREADS} -x {params.indexB} -U {input.reads1} -S {output.samB}'

rule merge_per:
    input:
        samA = os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapA.sam'),
        samB = os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapB.sam')
        # samA = DIR_FIRST_PASS + CHROM + '-per_hapA.sam',
        # samB = DIR_FIRST_PASS + CHROM + '-per_hapB.sam'
    output:
        path = os.path.join(DIR_FIRST_PASS, '{}-per.paths'.format(CHROM)),
        id = os.path.join(DIR_FIRST_PASS, '{}-per.ids'.format(CHROM)),
        merge_paths = os.path.join(DIR_FIRST_PASS, '{}-per.merge_paths'.format(CHROM)),
        samA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB.sam'.format(CHROM))
    params:
        os.path.join(DIR_FIRST_PASS, '{}-per-merged'.format(CHROM))
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

rule check_standard_onepass:
    input:
        maj = expand(
            os.path.join(DIR_FIRST_PASS, CHROM + '-h37maj.sam'),
            INDIV = INDIV),
        grc = expand(
            os.path.join(DIR_FIRST_PASS, CHROM + '-grch37.sam'),
            INDIV = INDIV),
        samA = expand(
            os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA.sam'.format(CHROM)),
            # os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapA.sam'),
            INDIV = INDIV),
        samB = expand(
            os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB.sam'.format(CHROM)),
            # os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapB.sam'),
            INDIV = INDIV),
        samAh = expand(
            os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapA_haploid.sam'),
            INDIV = INDIV),
        samBh = expand(
            os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapB_haploid.sam'),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'standard_onepass.done')))
