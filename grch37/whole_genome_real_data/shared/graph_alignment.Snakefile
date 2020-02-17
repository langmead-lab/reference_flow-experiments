'''
VG: build graph, index graph, alignment
'''
rule filter_vcf_by_freq:
    input:
        vcf = os.path.join(DIR, 'wg_filtered.vcf')
    output:
        vcfgz = os.path.join(DIR_VG, 'wg-{}.vcf.gz'.format(ALLELE_FREQ_FOR_VG)),
    shell:
        '{BCFTOOLS} view -O z -m2 -M2 -V mnps,other -q {ALLELE_FREQ_FOR_VG} {input.vcf} > {output.vcfgz};'

rule build_vcf_index:
    input:
        vcfgz = os.path.join(DIR_VG, 'wg-{}.vcf.gz'.format(ALLELE_FREQ_FOR_VG))
    output:
        csi = os.path.join(DIR_VG, 'wg-{}.vcf.gz.csi'.format(ALLELE_FREQ_FOR_VG)),
        tbi = os.path.join(DIR_VG, 'wg-{}.vcf.gz.tbi'.format(ALLELE_FREQ_FOR_VG))
    shell:
        'tabix -p vcf {input.vcfgz};'
        '{BCFTOOLS} index {input.vcfgz};'

rule vg_construct:
    input:
        vcfgz = os.path.join(DIR_VG, 'wg-{}.vcf.gz'.format(ALLELE_FREQ_FOR_VG)),
        tbi = os.path.join(DIR_VG, 'wg-{}.vcf.gz.tbi'.format(ALLELE_FREQ_FOR_VG))
    output:
        vg = os.path.join(DIR_VG, 'wg-{}.vg'.format(ALLELE_FREQ_FOR_VG)),
        log_construct = os.path.join(DIR_VG, 'wg-vg_construct.time_log')
    threads:
        THREADS
    shell:
        '{TIME} -v -o {output.log_construct} {VG} construct '
        '-v {input.vcfgz} -r {GENOME} -t {THREADS} > {output.vg}'

rule vg_index:
    input:
        vg = os.path.join(DIR_VG, 'wg-{}.vg'.format(ALLELE_FREQ_FOR_VG)),
    output:
        xg = os.path.join(DIR_VG, 'wg-{}.xg'.format(ALLELE_FREQ_FOR_VG)),
        gcsa = os.path.join(DIR_VG, 'wg-{}.gcsa'.format(ALLELE_FREQ_FOR_VG)),
        log_index = os.path.join(DIR_VG, 'wg-vg_index.time_log')
    params:
        k = 16,
        tmp = os.path.join(DIR_VG, 'tmp')
    threads:
        THREADS
    shell:
        'mkdir -p {params.tmp};'
        '{TIME} -v -o {output.log_index} {VG} index -x {output.xg} '
        '-g {output.gcsa} -k {params.k} -t {THREADS} -b {params.tmp} {input.vg}'

rule vg_map:
    input:
        vg = os.path.join(DIR_VG, 'wg-{}.vg'.format(ALLELE_FREQ_FOR_VG)),
        xg = os.path.join(DIR_VG, 'wg-{}.xg'.format(ALLELE_FREQ_FOR_VG)),
        gcsa = os.path.join(DIR_VG, 'wg-{}.gcsa'.format(ALLELE_FREQ_FOR_VG)),
        reads1 = READS1
    output:
        bam = os.path.join(DIR_FIRST_PASS,
            'wg-vg_{}.bam'.format(ALLELE_FREQ_FOR_VG)),
        log_map = os.path.join(DIR_VG, '{INDIV}-wg-vg_map.time_log')
    threads:
        THREADS
    params:
        prefix = os.path.join(DIR_VG, 'wg-{}'.format(ALLELE_FREQ_FOR_VG))
    shell:
        '{TIME} -v -o {output.log_map} {VG} map -d {params.prefix} -f {input.reads1} '
        '-t {THREADS} --surject-to bam > {output.bam}'

rule check_vg_index:
    input:
        vg = os.path.join(DIR_VG, 'wg-{}.vg'.format(ALLELE_FREQ_FOR_VG)),
        xg = os.path.join(DIR_VG, 'wg-{}.xg'.format(ALLELE_FREQ_FOR_VG)),
        gcsa = os.path.join(DIR_VG, 'wg-{}.gcsa'.format(ALLELE_FREQ_FOR_VG)),
    output:
        touch(temp(os.path.join(DIR, 'vg_index.done')))

rule check_vg_map:
    input:
        bam = expand(os.path.join(DIR_FIRST_PASS,
            'wg-vg_{}.bam'.format(ALLELE_FREQ_FOR_VG)),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'vg_map.done')))
