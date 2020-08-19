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

