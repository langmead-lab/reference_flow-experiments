'''
VG: build graph, index graph, alignment
'''
rule filter_vcf_by_freq:
    input:
        vcf = VCF
    output:
        vcfgz = os.path.join(DIR_VG, '{}-{}.vcf.gz'.format(CHROM, ALLELE_FREQ_FOR_VG)),
    shell:
        '{BCFTOOLS} view -O z -m2 -M2 -V mnps,other -q {ALLELE_FREQ_FOR_VG} {input.vcf} > {output.vcfgz};'

rule build_vcf_index:
    input:
        vcfgz = os.path.join(DIR_VG, '{}-{}.vcf.gz'.format(CHROM, ALLELE_FREQ_FOR_VG)),
    output:
        csi = os.path.join(DIR_VG, '{}-{}.vcf.gz.csi'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        tbi = os.path.join(DIR_VG, '{}-{}.vcf.gz.tbi'.format(CHROM, ALLELE_FREQ_FOR_VG))
    shell:
        'tabix -p vcf {input.vcfgz};'
        '{BCFTOOLS} index {input.vcfgz};'

rule vg_construct:
    input:
        vcfgz = os.path.join(DIR_VG, '{}-{}.vcf.gz'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        tbi = os.path.join(DIR_VG, '{}-{}.vcf.gz.tbi'.format(CHROM, ALLELE_FREQ_FOR_VG))
    output:
        vg = os.path.join(DIR_VG, '{}-{}.vg'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        log_construct = os.path.join(DIR_VG, '{}-vg_construct.time_log'.format(CHROM))
    threads:
        THREADS
    shell:
        '{TIME} -v -o {output.log_construct} {VG} construct '
        '-v {input.vcfgz} -r {GENOME} -t {THREADS} > {output.vg}'

rule vg_index:
    input:
        vg = os.path.join(DIR_VG, '{}-{}.vg'.format(CHROM, ALLELE_FREQ_FOR_VG)),
    output:
        xg = os.path.join(DIR_VG, '{}-{}.xg'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        gcsa = os.path.join(DIR_VG, '{}-{}.gcsa'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        # xg = os.path.join(DIR_VG, '{}.xg'.format(CHROM)),
        # gcsa = os.path.join(DIR_VG, '{}.gcsa'.format(CHROM)),
        log_index = os.path.join(DIR_VG, '{}-vg_index.time_log'.format(CHROM))
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
        vg = os.path.join(DIR_VG, '{}-{}.vg'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        xg = os.path.join(DIR_VG, '{}-{}.xg'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        gcsa = os.path.join(DIR_VG, '{}-{}.gcsa'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        reads1 = PREFIX_PER + '_1.fq',
    output:
        bam = os.path.join(DIR_FIRST_PASS,
            '{}-vg_{}.bam'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        log_map = os.path.join(DIR_VG, '{INDIV}' + '-{}-vg_map.time_log'.format(CHROM))
    threads:
        THREADS
    params:
        prefix = os.path.join(DIR_VG, '{}-{}'.format(CHROM, ALLELE_FREQ_FOR_VG))
    shell:
        '{TIME} -v -o {output.log_map} {VG} map -d {params.prefix} -f {input.reads1} '
        '-t {THREADS} --surject-to bam > {output.bam}'

rule sort_vg:
    input:
        os.path.join(DIR_FIRST_PASS, '{}-vg_{}.bam'.format(CHROM, ALLELE_FREQ_FOR_VG))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-vg_{}-sorted.sam'.format(CHROM, ALLELE_FREQ_FOR_VG))
    threads: THREADS
    shell:
        'samtools sort -@ {THREADS} -o {output} -O sam {input}'

rule check_vg_index:
    input:
        vg = os.path.join(DIR_VG, '{}-{}.vg'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        xg = os.path.join(DIR_VG, '{}-{}.xg'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        gcsa = os.path.join(DIR_VG, '{}-{}.gcsa'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        # vg = os.path.join(DIR_VG, '{}.vg'.format(CHROM)),
        # xg = os.path.join(DIR_VG, '{}.xg'.format(CHROM)),
        # gcsa = os.path.join(DIR_VG, '{}.gcsa'.format(CHROM))
    output:
        touch(temp(os.path.join(DIR, 'vg_index.done')))

rule check_vg_map:
    input:
        expand(
            os.path.join(DIR_FIRST_PASS,
            '{}-vg_{}.bam'.format(CHROM, ALLELE_FREQ_FOR_VG)),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, '{}-vg_{}-sorted.sam'.format(CHROM, ALLELE_FREQ_FOR_VG)),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'vg_map.done')))

rule calc_vg_accuracy:
    input:
        sam = os.path.join(DIR_FIRST_PASS, '{}-vg_{}-sorted.sam'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        gold = PREFIX_PER + '_1.sam',
        var_reads = PREFIX_PER + '.var',
    output:
        acc_log = os.path.join(DIR_FIRST_PASS, '{}-vg_{}.acc_log'.format(CHROM, ALLELE_FREQ_FOR_VG)),
        acc = os.path.join(DIR_RESULTS, '{INDIV}-' + '{}-vg_{}.acc'.format(CHROM, ALLELE_FREQ_FOR_VG))
    run:
        shell('{PYTHON} -O {DIR_SCRIPTS}/analyze_diploid_indels.py \
        -c {CHROM} -g {input.gold} -p 0 -vr {input.var_reads} \
        -n {input.sam} > {output.acc_log};')
        organize_accuracy(output.acc_log, output.acc)

rule check_vg_accuracy:
    input:
        expand(
            os.path.join(DIR_RESULTS, '{INDIV}-' + '{}-vg_{}.acc'.format(CHROM, ALLELE_FREQ_FOR_VG)),
            INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'vg_acc.done')))

