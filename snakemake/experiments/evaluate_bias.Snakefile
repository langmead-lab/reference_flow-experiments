# Get HETs for individuals
rule get_het:
    input:
        vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf')
    output:
        het = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf')
    shell:
        '{BCFTOOLS} view -s {wildcards.INDIV} -i "AC>0" -m2 -M2 {input.vcf} |'
        '{PYTHON} {DIR_SCRIPTS_EXP}/remove_het_overlapping_indel.py |'
        '{BCFTOOLS} view -g het > {output.het}'

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
        '{LEVIOSAM} serialize -v {input.vcf} -p {params.A} -g 0 -s {wildcards.INDIV} -k {input.length_map};'
        '{LEVIOSAM} serialize -v {input.vcf} -p {params.B} -g 1 -s {wildcards.INDIV} -k {input.length_map};'

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
        '{LEVIOSAM} lift -a {input.sam} -l {input.lft} -p {params} -t {threads}'

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
        '{LEVIOSAM} lift -a {input.samA} -l {input.lftA} -p {params.A} -t {threads}'

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
        '{LEVIOSAM} lift -a {input.samB} -l {input.lftB} -p {params.B} -t {threads}'

rule merge_per_allinone:
    input:
        A = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapA-liftover.sam'),
        B = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-hapB-liftover.sam')
    output:
        temp(os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover.sam'.format(CHROM)))
    shell:
        'cp {input.A} {output};'
        'grep -hv "^@" {input.B} >> {output}'


'''
1.  Find reads overlapping interesting regions (HETs here).
2.  Calculate biases.
3.  Summarize allelic bias results
'''
vcf_het = os.path.join(DIR, EXP_LABEL + '_{INDIV}_het_no_overlaps.vcf'),
# vcf_het = '/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/bias_inspector/variant_analysis/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer-norm-het_no_overlaps.vcf'
# vcf_het = '/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/refflow-exp/snakemake/SRR622457/personalized/new/compare/shared_het.vcf'

BIAS_OBJECTS = [
    EXP_LABEL + '-GRC',
    EXP_LABEL + '-major-liftover',
    EXP_LABEL + '-per-merged-liftover',
    # EXP_LABEL + '-vg-from_randflow_ld',
    EXP_LABEL + '-ht2_{}'.format(GRAPH_AF_THRSD),
    EXP_LABEL + '-vg_{}'.format(GRAPH_AF_THRSD)
]
rule sam_to_bam:
    input:
        os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}.sam')
    output:
        os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}.bam')
    wildcard_constraints:
        BIAS_OBJECTS='^((?!vg).)*$'
    shell:
        '{SAMTOOLS} view -bh -o {output} {input}'

ruleorder: sort_refflow > sort > sam_to_bam
rule sort:
    input:
        # os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}.bam')
        os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}.sam')
    output:
        os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}-sorted.bam')
    threads: THREADS
    shell:
        '{SAMTOOLS} sort -@ {threads} -o {output} -O BAM {input}'

ruleorder: filter_reads_overlapping_het_refflow > filter_reads_overlapping_het
rule filter_reads_overlapping_het:
    input:
        vcf = vcf_het,
        bam = os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}-sorted.bam')
    output:
        sam = os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}-{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | {SAMTOOLS} view -h >> {output.sam}'

rule filter_reads_overlapping_het_refflow:
    input:
        vcf = vcf_het,
        bam = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        sam = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '{INDIV}_het_no_overlaps.sam')
    shell:
        '{BEDTOOLS} intersect -a {input.bam} -b {input.vcf} | {SAMTOOLS} view -h >> {output.sam}'

ruleorder: calc_allelic_bias_refflow > calc_allelic_bias
rule calc_allelic_bias:
    input:
        sam = os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}-{INDIV}_het_no_overlaps.sam'),
        vcf = vcf_het,
    output:
        os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_SCRIPTS_EXP}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

rule calc_allelic_bias_refflow:
    input:
        sam = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '{INDIV}_het_no_overlaps.sam'),
        vcf = vcf_het,
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-refflow-{}-{}'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '-{INDIV}.allele.bias')
    shell:
        '{PYTHON} {DIR_SCRIPTS_EXP}/evaluate_allele.py -v {input.vcf} -s {input.sam} -f {GENOME} -o {output}'

ruleorder: summarize_refflow > summarize
rule summarize:
    input:
        os.path.join(DIR_FIRST_PASS, '{BIAS_OBJECTS}-{INDIV}.allele.bias')
    output:
        os.path.join(DIR_RESULTS_BIAS, '{BIAS_OBJECTS}-{INDIV}.summary.bias')
    run:
        summarize_allelc_bias(input, output)

rule summarize_refflow:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-refflow-{}-{}'.format(ALN_MAPQ_THRSD, POP_DIRNAME) + '-{INDIV}.allele.bias')
    output:
        os.path.join(DIR_RESULTS_BIAS, EXP_LABEL + '-' + POP_DIRNAME + '-{INDIV}.summary.bias')
    run:
        summarize_allelc_bias(input, output)

rule check_refbias_and_write_to_tsv:
    input:
        expand(
            os.path.join(DIR_RESULTS_BIAS, '{BIAS_OBJECTS}-{INDIV}.summary.bias'),
            BIAS_OBJECTS = BIAS_OBJECTS, INDIV = INDIV),
        expand(
           os.path.join(DIR_RESULTS_BIAS, EXP_LABEL + '-' + POP_DIRNAME + '-{INDIV}.summary.bias'),
            INDIV = INDIV),
    output:
        tsv = os.path.join(DIR_RESULTS_BIAS, 'bias.tsv'),
        done = touch(temp(os.path.join(DIR, 'bias_exp.done')))
    run:
        dict_indiv_to_pop = build_dict_indiv_to_pop(FAMILY)
        dict_pop_to_spop = build_dict_pop_to_spop(SPOP)

        tmp = os.listdir(DIR_RESULTS_BIAS)
        list_fn = []
        for t in tmp:
            if t.endswith('.bias'):
                list_fn.append(t)
        list_fn = [os.path.join(DIR_RESULTS_BIAS, fn) for fn in list_fn]
        list_indiv = []
        list_method = []
        list_pass_rate = []
        list_sum_ref = []
        list_sum_alt = []
        list_ref_to_alt = []
        list_num_refbias = []
        list_num_altbias = []
        for fn in list_fn:
            basename = os.path.basename(fn).split('-')
            list_indiv.append(basename[-1].split('.')[0])
            method = '-'.join(basename[:-1])
            list_method.append(method.split('-')[1])
            # basename = os.path.basename(fn).split('-')
            # list_indiv.append(basename[0])
            # method = '-'.join(basename[1:])
            # list_method.append(method[: method.rfind('.bias')])
            # list_method.append(basename[1][:basename[1].rfind('.bias')])
            with open(fn, 'r') as f:
                list_pass_rate.append(f.readline().rstrip())
                list_sum_ref.append(f.readline().rstrip())
                list_sum_alt.append(f.readline().rstrip())
                list_ref_to_alt.append(f.readline().rstrip())
                list_num_refbias.append(f.readline().rstrip())
                list_num_altbias.append(f.readline().rstrip())
        df = pd.DataFrame()
        df['INDIV'] = list_indiv
        df['Super_Population'] = [dict_pop_to_spop[dict_indiv_to_pop[s]] for s in list_indiv]
        df['Method'] = list_method
        df['Pass_Rate({})'.format(BIAS_MIN_READ_COUNT)] = list_pass_rate
        df['SUM_REF'] = list_sum_ref
        df['SUM_ALT'] = list_sum_alt
        df['REF_TO_ALT'] = list_ref_to_alt
        df['NUM_REFBIAS({})'.format(0.5+BIAS_TAIL_THRDS)] = list_num_refbias
        df['NUM_ALTBIAS({})'.format(0.5-BIAS_TAIL_THRDS)] = list_num_altbias
        df = df.sort_values('INDIV')
        df.to_csv(output.tsv, sep = '\t', index = None, float_format = '%.4f')
