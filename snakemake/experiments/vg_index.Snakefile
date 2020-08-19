'''
VG: build graph, index graph, alignment
'''
rule vg_construct:
    input:
        vcfgz = os.path.join(DIR_VG, EXP_LABEL + '-{}.vcf.gz'.format(GRAPH_AF_THRSD)),
        tbi = os.path.join(DIR_VG, EXP_LABEL + '-{}.vcf.gz.tbi'.format(GRAPH_AF_THRSD))
    output:
        vg = os.path.join(DIR_VG, EXP_LABEL + '-{}.vg'.format(GRAPH_AF_THRSD)),
        log_construct = os.path.join(DIR_VG, EXP_LABEL + '-vg_construct.time_log')
    threads:
        THREADS
    shell:
        '{TIME} -v -o {output.log_construct} {VG} construct '
        '-v {input.vcfgz} -r {GENOME} -t {threads} > {output.vg}'

rule vg_index:
    input:
        vg = os.path.join(DIR_VG, EXP_LABEL + '-{}.vg'.format(GRAPH_AF_THRSD))
    output:
        xg = os.path.join(DIR_VG, EXP_LABEL + '-{}.xg'.format(GRAPH_AF_THRSD)),
        gcsa = os.path.join(DIR_VG, EXP_LABEL + '-{}.gcsa'.format(GRAPH_AF_THRSD)),
        log_index = os.path.join(DIR_VG, EXP_LABEL + '-vg_index.time_log')
    params:
        k = 16,
        tmp = os.path.join(DIR_VG, 'tmp')
    threads:
        THREADS
    shell:
        'mkdir -p {params.tmp};'
        '{TIME} -v -o {output.log_index} {VG} index -x {output.xg} '
        '-g {output.gcsa} -k {params.k} -t {threads} -b {params.tmp} {input.vg}'

# rule sort_vg:
#     input:
#         os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}.bam'.format(GRAPH_AF_THRSD))
#     output:
#         os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}-sorted.bam'.format(GRAPH_AF_THRSD))
#     threads: THREADS
#     shell:
#         '{SAMTOOLS} sort -@ {threads} -o {output} {input}'

rule check_vg_index:
    input:
        vg = os.path.join(DIR_VG, EXP_LABEL + '-{}.vg'.format(GRAPH_AF_THRSD)),
        xg = os.path.join(DIR_VG, EXP_LABEL + '-{}.xg'.format(GRAPH_AF_THRSD)),
        gcsa = os.path.join(DIR_VG, EXP_LABEL + '-{}.gcsa'.format(GRAPH_AF_THRSD))
    output:
        touch(temp(os.path.join(DIR, 'vg_index.done')))

rule check_vg_map:
    input:
        expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}.bam'.format(GRAPH_AF_THRSD)),
            INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}-sorted.bam'.format(GRAPH_AF_THRSD)),
        #     INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'vg_map.done')))
