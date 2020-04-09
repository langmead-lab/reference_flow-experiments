# Get HETs for individuals
rule get_het:
    input:
        vcf = os.path.join(DIR, 'wg_filtered.vcf')
    output:
        het = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf')
    shell:
        '{BCFTOOLS} view -s {INDIV} -i "AC>0" -m2 -M2 {input.vcf} |'
        '{PYTHON} {DIR_SCRIPTS}/remove_het_overlapping_indel.py |'
        '{BCFTOOLS} view -g het > {output.het}'

rule liftover_lift_major:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-major.sam'),
        lft = os.path.join(DIR_MAJOR, 'wg-major.lft')
    output:
        temp(os.path.join(DIR_FIRST_PASS, 'wg-major-liftover.sam'))
    params:
        os.path.join(DIR_FIRST_PASS, 'wg-major-liftover')
    threads: THREADS
    shell:
        '{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params} -t {threads}'

rule liftover_serialize_per:
    input:
        vcf = os.path.join(DIR_PER, '{INDIV}-wg-per.vcf'),
        length_map = LENGTH_MAP
    output:
        lftA = os.path.join(DIR_PER, 'wg-perA.lft'),
        lftB = os.path.join(DIR_PER, 'wg-perB.lft')
    params:
        A = os.path.join(DIR_PER, 'wg-perA'),
        B = os.path.join(DIR_PER, 'wg-perB')
    shell:
        '{LIFTOVER} serialize -v {input.vcf} -p {params.A} -g 0 -s {wildcards.INDIV} -k {input.length_map};'
        '{LIFTOVER} serialize -v {input.vcf} -p {params.B} -g 1 -s {wildcards.INDIV} -k {input.length_map};'

rule liftover_lift_perA:
    input:
        samA = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapA.sam'),
        lftA = os.path.join(DIR_PER, 'wg-perA.lft'),
        vcf = os.path.join(DIR_PER, '{INDIV}-wg-per.vcf')
    output:
        A = temp(os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapA-liftover.sam')),
    params:
        A = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapA-liftover'),
    threads: THREADS
    shell:
        '{LIFTOVER} lift -a {input.samA} -l {input.lftA} -p {params.A} -t {threads}'

rule liftover_lift_perB:
    input:
        samB = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapB.sam'),
        lftB = os.path.join(DIR_PER, 'wg-perB.lft'),
        vcf = os.path.join(DIR_PER, '{INDIV}-wg-per.vcf')
    output:
        B = temp(os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapB-liftover.sam'))
    params:
        B = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapB-liftover')
    threads: THREADS
    shell:
        '{LIFTOVER} lift -a {input.samB} -l {input.lftB} -p {params.B} -t {threads}'

rule merge_per_allinone:
    input:
        A = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapA-liftover.sam'),
        B = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapB-liftover.sam')
    output:
        temp(os.path.join(DIR_FIRST_PASS, 'wg-per-merged-liftover.sam'.format(CHROM)))
    shell:
        'cp {input.A} {output};'
        'grep -hv "^@" {input.B} >> {output}'

rule sort_grc:
    input:
        os.path.join(DIR_FIRST_PASS, 'wg-GRC.sam')
    output:
        os.path.join(DIR_FIRST_PASS, 'wg-GRC-sorted.bam')
    threads: THREADS
    run:
        shell('samtools sort -@ {threads} -o {output} -O BAM {input}')

rule sort_major:
    input:
        os.path.join(DIR_FIRST_PASS, 'wg-major-liftover.sam')
    output:
        os.path.join(DIR_FIRST_PASS, 'wg-major-liftover-sorted.bam')
    threads: THREADS
    run:
        shell('samtools sort -@ {threads} -o {output} -O BAM {input}')
        
rule sort_per:
    input:
        os.path.join(DIR_FIRST_PASS, 'wg-per-merged-liftover.sam')
    output:
        os.path.join(DIR_FIRST_PASS, 'wg-per-merged-liftover-sorted.bam')
    threads: THREADS
    run:
        shell('samtools sort -@ {threads} -o {output} -O BAM {input}')

'''
Find reads overlapping interesting regions (HETs here).
Using bedtools improves speed.
'''
rule filter_reads_overlapping_het_grc:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        # bam = os.path.join(DIR_FIRST_PASS, 'wg-GRC.sam')
        bam = os.path.join(DIR_FIRST_PASS, 'wg-GRC-sorted.bam')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-GRC-{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | samtools view -h >> {output.sam}'
        
rule filter_reads_overlapping_het_major:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        bam = os.path.join(DIR_FIRST_PASS, 'wg-major-liftover-sorted.bam')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-major-{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | samtools view -h >> {output.sam}'

rule filter_reads_overlapping_het_refflow:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        bam = os.path.join(DIR_SECOND_PASS,
            'wg-refflow-{}-{}-liftover-sorted.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        # sam = os.path.join(DIR_FIRST_PASS, 'wg-major-{INDIV}_het_no_overlaps.sam')
        sam = os.path.join(DIR_SECOND_PASS,
            'wg-refflow-{}-{}-'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | samtools view -h >> {output.sam}'

rule filter_reads_overlapping_het_per:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        # bam = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-liftover.sam')
        bam = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-liftover-sorted.bam')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-per-{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | samtools view -h >> {output.sam}'

DIR_BIAS_EVAL = '/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/bias_inspector/src'
rule calc_allelic_bias_grc:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-GRC-{INDIV}_het_no_overlaps.sam'),
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
    output:
        os.path.join(DIR_FIRST_PASS, 'wg-GRC-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_BIAS_EVAL}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

rule calc_allelic_bias_major:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-major-{INDIV}_het_no_overlaps.sam'),
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
    output:
        os.path.join(DIR_FIRST_PASS, 'wg-major-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_BIAS_EVAL}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

rule calc_allelic_bias_refflow:
    input:
        sam = os.path.join(DIR_SECOND_PASS,
            'wg-refflow-{}-{}-'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '{INDIV}_het_no_overlaps.sam'),
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
    output:
        os.path.join(DIR_FIRST_PASS, 'wg-refflow-{}-{}'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_BIAS_EVAL}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

rule calc_allelic_bias_per:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-per-{INDIV}_het_no_overlaps.sam'),
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
    output:
        os.path.join(DIR_FIRST_PASS, 'wg-per-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_BIAS_EVAL}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

rule check_bias_exp:
    input:
        expand(
            os.path.join(DIR_FIRST_PASS, 'wg-GRC-{INDIV}.allele.bias'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, 'wg-major-{INDIV}.allele.bias'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, 'wg-refflow-{}-{}'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '-{INDIV}.allele.bias'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, 'wg-per-{INDIV}.allele.bias'),
            INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'bias_exp.done')))
