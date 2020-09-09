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

'''
Align reads using HISAT2
'''

DIR_HT2 = os.path.join(DIR, 'hisat2')

rule get_ht_snps:
    input:
        vcfgz = os.path.join(DIR, EXP_LABEL + '-{}.vcf.gz'.format(GRAPH_AF_THRSD))
    output:
        vcf = temp(os.path.join(DIR_HT2, EXP_LABEL + '_filtered-{}.vcf'.format(GRAPH_AF_THRSD))),
        snp = os.path.join(DIR_HT2, EXP_LABEL + '_filtered-{}.ht2.snp'.format(GRAPH_AF_THRSD))
    shell:
        'bgzip -cd {input.vcfgz} > {output.vcf};'
        '{PYTHON} {DIR_SCRIPTS_EXP}/convert_vcf_to_hisat2_snp.py -v {output.vcf} -s {output.snp}'

rule build_ht_index:
    input:
       snp = os.path.join(DIR_HT2, EXP_LABEL + '_filtered-{}.ht2.snp'.format(GRAPH_AF_THRSD))
    output:
        expand(
            os.path.join(DIR_HT2, 'indexes/ht2-{}.'.format(GRAPH_AF_THRSD) + '{idx}.ht2'),
            idx=[1, 2, 3, 4, 5, 6, 7, 8])
    params:
        os.path.join(DIR_HT2, 'indexes/ht2-{}'.format(GRAPH_AF_THRSD))
    threads:
        THREADS
    shell:
        'hisat2-build -p {threads} --snp {input.snp} -f {GENOME} {params}'

rule check_ht_align:
    input:
        expand(
            os.path.join(
                DIR_FIRST_PASS,
                EXP_LABEL + '-ht2_{}.sam'.format(GRAPH_AF_THRSD)
            ),
            INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'ht_align.done')))
