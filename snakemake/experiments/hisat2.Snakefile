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
