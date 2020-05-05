# Get HETs for individuals
rule get_het:
    input:
        vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf')
    output:
        het = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf')
    shell:
        '{BCFTOOLS} view -s {INDIV} -i "AC>0" -m2 -M2 {input.vcf} |'
        '{PYTHON} {DIR_SCRIPTS_EXP}/remove_het_overlapping_indel.py |'
        '{BCFTOOLS} view -g het > {output.het}'

rule liftover_lift_major:
    input:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major.sam'),
        lft = os.path.join(DIR_MAJOR, EXP_LABEL + '-major.lft')
    output:
        temp(os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover.sam'))
    params:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover')
    threads: THREADS
    shell:
        '{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params} -t {threads}'

rule liftover_serialize_per:
    input:
        vcf = os.path.join(DIR_PER, EXP_LABEL + '-per.vcf'),
        length_map = LENGTH_MAP
    output:
        lftA = os.path.join(DIR_PER, EXP_LABEL + '-perA.lft'),
        lftB = os.path.join(DIR_PER, EXP_LABEL + '-perB.lft')
    params:
        A = os.path.join(DIR_PER, EXP_LABEL + '-perA'),
        B = os.path.join(DIR_PER, EXP_LABEL + '-perB')
    shell:
        '{LIFTOVER} serialize -v {input.vcf} -p {params.A} -g 0 -s {wildcards.INDIV} -k {input.length_map};'
        '{LIFTOVER} serialize -v {input.vcf} -p {params.B} -g 1 -s {wildcards.INDIV} -k {input.length_map};'

rule liftover_lift_perA:
    input:
        samA = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapA.sam'),
        lftA = os.path.join(DIR_PER, EXP_LABEL + '-perA.lft')
    output:
        A = temp(os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapA-liftover.sam')),
    params:
        A = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapA-liftover'),
    threads: THREADS
    shell:
        '{LIFTOVER} lift -a {input.samA} -l {input.lftA} -p {params.A} -t {threads}'

rule liftover_lift_perB:
    input:
        samB = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapB.sam'),
        lftB = os.path.join(DIR_PER, EXP_LABEL + '-perB.lft')
    output:
        B = temp(os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapB-liftover.sam'))
    params:
        B = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapB-liftover')
    threads: THREADS
    shell:
        '{LIFTOVER} lift -a {input.samB} -l {input.lftB} -p {params.B} -t {threads}'

rule merge_per_allinone:
    input:
        A = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapA-liftover.sam'),
        B = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapB-liftover.sam')
    output:
        temp(os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover.sam'.format(CHROM)))
    shell:
        'cp {input.A} {output};'
        'grep -hv "^@" {input.B} >> {output}'

rule sort_grc:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC.sam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-sorted.bam')
    threads: THREADS
    shell:
        'samtools sort -@ {threads} -o {output} -O BAM {input}'

rule sort_major:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover.sam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted.bam')
    threads: THREADS
    shell:
        'samtools sort -@ {threads} -o {output} -O BAM {input}'
        
rule sort_per:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover.sam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted.bam')
    threads: THREADS
    shell:
        'samtools sort -@ {threads} -o {output} -O BAM {input}'

'''
Find reads overlapping interesting regions (HETs here).
Using bedtools improves speed.
'''
# vcf_het = '/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/bias_inspector/variant_analysis/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer-norm-het_no_overlaps.vcf'
# vcf_het = '/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/refflow-exp/snakemake/SRR622457/personalized/new/compare/shared_het.vcf'
vcf_het = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
rule filter_reads_overlapping_het_grc:
    input:
        # vcf = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
        vcf = vcf_het,
        bam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-sorted.bam')
    output:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | samtools view -h >> {output.sam}'
        
rule filter_reads_overlapping_het_major:
    input:
        # vcf = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
        vcf = vcf_het,
        bam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted.bam')
    output:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | samtools view -h >> {output.sam}'

rule filter_reads_overlapping_het_refflow:
    input:
        # vcf = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
        vcf = vcf_het,
        bam = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        # sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-{INDIV}_het_no_overlaps.sam')
        sam = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | samtools view -h >> {output.sam}'

rule filter_reads_overlapping_het_per:
    input:
        # vcf = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
        vcf = vcf_het,
        bam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted.bam')
    output:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | samtools view -h >> {output.sam}'

DIR_BIAS_EVAL = '/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/bias_inspector/src'
rule calc_allelic_bias_grc:
    input:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-{INDIV}_het_no_overlaps.sam'),
        #vcf = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
        vcf = vcf_het,
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_BIAS_EVAL}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

rule calc_allelic_bias_major:
    input:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-{INDIV}_het_no_overlaps.sam'),
        #vcf = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
        vcf = vcf_het,
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_BIAS_EVAL}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

rule calc_allelic_bias_refflow:
    input:
        sam = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '{INDIV}_het_no_overlaps.sam'),
        #vcf = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
        vcf = vcf_het,
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-refflow-{}-{}'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_BIAS_EVAL}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

rule calc_allelic_bias_per:
    input:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-{INDIV}_het_no_overlaps.sam'),
        #vcf = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
        vcf = vcf_het,
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_BIAS_EVAL}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

rule check_bias_exp:
    input:
        expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-{INDIV}.allele.bias'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-{INDIV}.allele.bias'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-refflow-{}-{}'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '-{INDIV}.allele.bias'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-{INDIV}.allele.bias'),
            INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'bias_exp.done')))
