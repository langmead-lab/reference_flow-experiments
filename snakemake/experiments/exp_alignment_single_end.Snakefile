''' Reference flow: first pass '''
rule align_to_ref:
    input:
        reads1 = READS1,
        reads2 = READS2,
        idx = expand(
            os.path.join(DIR, 'grc/indexes/' + EXP_LABEL + '.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        index = os.path.join(DIR, 'grc/indexes/' + EXP_LABEL)
    output:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC.sam')
    threads: THREADS
    shell:
        'bowtie2 --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam}'

DIR_PER_IDX = os.path.join(DIR, 'personalized/{INDIV}/indexes/')
PER_LABEL = ['A', 'B']
''' Personalized '''
rule align_to_personalized:
    input:
        reads1 = READS1,
        reads2 = READS2,
        idx1 = os.path.join(DIR_PER_IDX, EXP_LABEL + '-per{PER_LABEL}.1.bt2'),
        idx2 = os.path.join(DIR_PER_IDX, EXP_LABEL + '-per{PER_LABEL}.2.bt2'),
        idx3 = os.path.join(DIR_PER_IDX, EXP_LABEL + '-per{PER_LABEL}.3.bt2'),
        idx4 = os.path.join(DIR_PER_IDX, EXP_LABEL + '-per{PER_LABEL}.4.bt2'),
        idxr1 = os.path.join(DIR_PER_IDX, EXP_LABEL + '-per{PER_LABEL}.rev.1.bt2'),
        idxr2 = os.path.join(DIR_PER_IDX, EXP_LABEL + '-per{PER_LABEL}.rev.2.bt2')
    params:
        index = os.path.join(DIR_PER_IDX, EXP_LABEL + '-per{PER_LABEL}')
    output:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per_hap{PER_LABEL}.sam')
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam}'

rule merge_per:
    input:
        samA = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per_hapA.sam'),
        samB = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per_hapB.sam')
    output:
        path = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per.paths'),
        label = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per.ids'),
        merge_paths = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per.merge_paths'),
        samA = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapA.sam'),
        samB = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapB.sam')
    params:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged')
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
            -l {output.merge_paths}')

rule check_alignment_grc_and_per:
    input:
        grc = expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC.sam'),
            INDIV = INDIV),
        perA = expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapA.sam'),
            INDIV = INDIV),
        perB = expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapB.sam'),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'alignment_grc_per.done')))
