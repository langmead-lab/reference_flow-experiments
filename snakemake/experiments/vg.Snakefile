'''
VG: build graph, index graph, alignment
'''
rule filter_vcf_by_freq:
    input:
        vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf')
    output:
        vcfgz = os.path.join(DIR_VG, EXP_LABEL + '-{}.vcf.gz'.format(ALLELE_FREQ_FOR_VG))
    shell:
        '{BCFTOOLS} view -O z -m2 -M2 -V mnps,other -q {ALLELE_FREQ_FOR_VG} {input.vcf} > {output.vcfgz};'

rule build_vcf_index:
    input:
        vcfgz = os.path.join(DIR_VG, EXP_LABEL + '-{}.vcf.gz'.format(ALLELE_FREQ_FOR_VG))
    output:
        csi = os.path.join(DIR_VG, EXP_LABEL + '-{}.vcf.gz.csi'.format(ALLELE_FREQ_FOR_VG)),
        tbi = os.path.join(DIR_VG, EXP_LABEL + '-{}.vcf.gz.tbi'.format(ALLELE_FREQ_FOR_VG))
    shell:
        'tabix -p vcf {input.vcfgz};'
        '{BCFTOOLS} index {input.vcfgz};'

rule vg_construct:
    input:
        vcfgz = os.path.join(DIR_VG, EXP_LABEL + '-{}.vcf.gz'.format(ALLELE_FREQ_FOR_VG)),
        tbi = os.path.join(DIR_VG, EXP_LABEL + '-{}.vcf.gz.tbi'.format(ALLELE_FREQ_FOR_VG))
    output:
        vg = os.path.join(DIR_VG, EXP_LABEL + '-{}.vg'.format(ALLELE_FREQ_FOR_VG)),
        log_construct = os.path.join(DIR_VG, EXP_LABEL + '-vg_construct.time_log')
    threads:
        THREADS
    shell:
        '{TIME} -v -o {output.log_construct} {VG} construct '
        '-v {input.vcfgz} -r {GENOME} -t {threads} > {output.vg}'

rule vg_index:
    input:
        vg = os.path.join(DIR_VG, EXP_LABEL + '-{}.vg'.format(ALLELE_FREQ_FOR_VG))
    output:
        xg = os.path.join(DIR_VG, EXP_LABEL + '-{}.xg'.format(ALLELE_FREQ_FOR_VG)),
        gcsa = os.path.join(DIR_VG, EXP_LABEL + '-{}.gcsa'.format(ALLELE_FREQ_FOR_VG)),
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

rule vg_map:
    input:
        vg = os.path.join(DIR_VG, EXP_LABEL + '-{}.vg'.format(ALLELE_FREQ_FOR_VG)),
        xg = os.path.join(DIR_VG, EXP_LABEL + '-{}.xg'.format(ALLELE_FREQ_FOR_VG)),
        gcsa = os.path.join(DIR_VG, EXP_LABEL + '-{}.gcsa'.format(ALLELE_FREQ_FOR_VG)),
        reads1 = READS1,
        reads2 = READS2
    output:
        bam = os.path.join(
            DIR_FIRST_PASS,
            EXP_LABEL + '-vg_{}.bam'.format(ALLELE_FREQ_FOR_VG)
        ),
        log_map = os.path.join(
            DIR_VG, '{INDIV}' + EXP_LABEL + '-vg_map.time_log'
        )
    threads:
        THREADS
    params:
        prefix = os.path.join(DIR_VG, EXP_LABEL + '-{}'.format(ALLELE_FREQ_FOR_VG))
    shell:
        '{TIME} -v -o {output.log_map} {VG} map -d {params.prefix} -f {input.reads1} -f {input.reads2} '
        '-t {threads} --surject-to bam > {output.bam}'

rule sort_vg:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}.bam'.format(ALLELE_FREQ_FOR_VG))
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}-sorted.bam'.format(ALLELE_FREQ_FOR_VG))
    threads: THREADS
    shell:
        '{SAMTOOLS} sort -@ {threads} -o {output} {input}'

rule check_vg_index:
    input:
        vg = os.path.join(DIR_VG, EXP_LABEL + '-{}.vg'.format(ALLELE_FREQ_FOR_VG)),
        xg = os.path.join(DIR_VG, EXP_LABEL + '-{}.xg'.format(ALLELE_FREQ_FOR_VG)),
        gcsa = os.path.join(DIR_VG, EXP_LABEL + '-{}.gcsa'.format(ALLELE_FREQ_FOR_VG))
    output:
        touch(temp(os.path.join(DIR, 'vg_index.done')))

rule check_vg_map:
    input:
        expand(
            os.path.join(
                DIR_FIRST_PASS,
                EXP_LABEL + '-vg_{}.bam'.format(ALLELE_FREQ_FOR_VG)
            ),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-vg_{}-sorted.bam'.format(ALLELE_FREQ_FOR_VG)),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'vg_map.done')))
