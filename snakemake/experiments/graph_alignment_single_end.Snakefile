'''Perform `vg map` in single-end mode.
'''
rule vg_map:
    input:
        vg = os.path.join(DIR_VG, EXP_LABEL + '-{}.vg'.format(GRAPH_AF_THRSD)),
        xg = os.path.join(DIR_VG, EXP_LABEL + '-{}.xg'.format(GRAPH_AF_THRSD)),
        gcsa = os.path.join(DIR_VG, EXP_LABEL + '-{}.gcsa'.format(GRAPH_AF_THRSD)),
        reads1 = READS1
    output:
        bam = os.path.join(
            DIR_FIRST_PASS,
            EXP_LABEL + '-vg_{}.bam'.format(GRAPH_AF_THRSD)
        ),
        log_map = os.path.join(
            DIR_VG, '{INDIV}' + EXP_LABEL + '-vg_map.time_log'
        )
    threads:
        THREADS
    params:
        prefix = os.path.join(DIR_VG, EXP_LABEL + '-{}'.format(GRAPH_AF_THRSD))
    shell:
        '{TIME} -v -o {output.log_map} {VG} map -d {params.prefix} -f {input.reads1} '
        '-t {threads} --surject-to bam > {output.bam}'

'''Align using HISAT2'''
rule ht_align:
    input:
        idx = expand(
            os.path.join(DIR_HT2, 'indexes/ht2-{}.'.format(GRAPH_AF_THRSD) + '{idx}.ht2'),
            idx=[1, 2, 3, 4, 5, 6, 7, 8]),
        reads1 = READS1
    output:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-ht2_{}.sam'.format(GRAPH_AF_THRSD))
    params:
        os.path.join(DIR_HT2, 'indexes/ht2-{}'.format(GRAPH_AF_THRSD))
    threads:
        THREADS
    shell:
        'hisat2 -p {threads} -x {params} -U {input.reads1} -k 10 --no-spliced-alignment '
        '--no-temp-splicesite | {SAMTOOLS} view -h -F 256 > {output.sam}'
        # 'hisat2 -p {threads} -x {params} -U {input.reads1} -S {output.sam} -k 10 --no-spliced-alignment --no-temp-splicesite'
