''' Prepare HET sites from the 1KG VCF '''
#: SNV-only
rule build_individual_het:
    input:
        vcf = PHASED_VCF_F
    output:
        os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}.vcf')
    shell:
        '{BCFTOOLS} view -s {wildcards.INDIV} {input.vcf} | \
         {BCFTOOLS} view -i "AC>0" -v snps -g het -m2 -M2 > {output}'

#: HETs including SNVs and INDELs
rule get_het_with_indel:
    input:
        vcf = PHASED_VCF_F
    output:
        het = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_withindel.vcf')
    shell:
        '{BCFTOOLS} view -s {wildcards.INDIV} {input.vcf} | {BCFTOOLS} view -i "AC>0" -g het -m2 -M2 > {output.het}'

#: Remove INDELS and SNVs overlapping INDELs
rule get_het_with_indel_processed:
    input:
        het = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_withindel.vcf')
    output:
        het = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf')
    shell:
        'cat {input.het} | python {DIR_SCRIPTS}/remove_het_overlapping_indel.py > {output.het}'

''' Allelic bias measurement '''
rule calc_major_bias_gnomad:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-liftover-sorted.sam'.format(CHROM))
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'major-gnomad-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'major-gnomad-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'major-gnomad-refbias.txt')
    shell:
        'echo "major-gnomad" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
            -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_major_bias_onekg:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-liftover-sorted.sam'.format(CHROM))
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'major-1kg-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'major-1kg-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'major-1kg-refbias.txt')
    shell:
        'echo "major-1kg" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
            -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_grc_bias:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-grc-sorted.sam'.format(CHROM))
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'grc-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'grc-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'grc-refbias.txt')
    shell:
        'echo "grc" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
           -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_per_bias:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_FIRST_PASS,
            'chr{}-per-merged-liftover-sorted.sam'.format(CHROM))
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'per-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'per-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'per-refbias.txt')
    shell:
        'echo "per" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
            -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_per_bias_haploid_setting:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_FIRST_PASS,
            '{}-per_h2h-merged-liftover-sorted.sam'.format(CHROM))
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'per_h2h-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'per_h2h-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'per_h2h-refbias.txt')
    shell:
        'echo "per_h2h" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
            -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_refflow_bias_gnomad:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_SECOND_PASS,
            'chr{}-refflow-{}-{}-liftover-gnomad-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        list_path = os.path.join(DIR_SECOND_PASS, 'refflow-{}-gnomad.paths'.format(POP_DIRNAME)),
        list_id = os.path.join(DIR_SECOND_PASS, 'refflow-{}-gnomad.ids'.format(POP_DIRNAME)),
        bias = os.path.join(DIR_SECOND_PASS, 'refflow-{}-gnomad.txt'.format(POP_DIRNAME))
    shell:
        'echo "refflow" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
           -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_refflow_bias_onekg:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_SECOND_PASS,
            'chr{}-refflow-{}-{}-liftover-1kg-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        list_path = os.path.join(DIR_SECOND_PASS, 'refflow-{}-1kg.paths'.format(POP_DIRNAME)),
        list_id = os.path.join(DIR_SECOND_PASS, 'refflow-{}-1kg.ids'.format(POP_DIRNAME)),
        bias = os.path.join(DIR_SECOND_PASS, 'refflow-{}-1kg.txt'.format(POP_DIRNAME))
    shell:
        'echo "refflow" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
           -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_vg_bias:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_FIRST_PASS, '{}-vg_{}-sorted.sam'.format(CHROM, ALLELE_FREQ_FOR_VG))
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.paths'.format(ALLELE_FREQ_FOR_VG)),
        list_id = os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.ids'.format(ALLELE_FREQ_FOR_VG)),
        bias = os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.txt'.format(ALLELE_FREQ_FOR_VG))
    params:
        id = 'vg_{}'.format(ALLELE_FREQ_FOR_VG)
    shell:
        'echo {params.id} > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
           -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule summarize_grc:
    input:
        os.path.join(DIR_FIRST_PASS,
            'grc-refbias.txt')
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-grc.bias')
    run:
        summarize_allelc_bias(input, output)

rule summarize_major_gnomad:
    input:
        os.path.join(DIR_FIRST_PASS,
            'major-gnomad-refbias.txt')
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-major-gnomad.bias')
    run:
        summarize_allelc_bias(input, output)

rule summarize_major_onekg:
    input:
        os.path.join(DIR_FIRST_PASS,
            'major-1kg-refbias.txt')
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-major-1kg.bias')
    run:
        summarize_allelc_bias(input, output)

rule summarize_per:
    input:
        os.path.join(DIR_FIRST_PASS,
            'per-refbias.txt')
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-per.bias')
    run:
        summarize_allelc_bias(input, output)

rule summarize_per_haploid_setting:
    input:
        os.path.join(DIR_FIRST_PASS,
            'per_h2h-refbias.txt')
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-per_h2h.bias')
    run:
        summarize_allelc_bias(input, output)

rule summarize_refflow_gnomad:
    input:
        os.path.join(DIR_SECOND_PASS,
            'refflow-{}-gnomad.txt'.format(POP_DIRNAME))
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-' + POP_DIRNAME + '-gnomad.bias')
    run:
        summarize_allelc_bias(input, output)

rule summarize_refflow_onekg:
    input:
        os.path.join(DIR_SECOND_PASS,
            'refflow-{}-1kg.txt'.format(POP_DIRNAME))
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-' + POP_DIRNAME + '-1kg.bias')
    run:
        summarize_allelc_bias(input, output)

rule summarize_vg:
    input:
        os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.txt'.format(ALLELE_FREQ_FOR_VG))
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-vg_' + '{}.bias'.format(ALLELE_FREQ_FOR_VG))
    run:
        summarize_allelc_bias(input, output)

rule check_refbias_and_write_to_tsv:
    input:
        expand(
            os.path.join(DIR_RESULTS_BIAS, '{INDIV}-grc.bias'),
            INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_RESULTS_BIAS, '{INDIV}-major-gnomad.bias'),
        #     INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS_BIAS, '{INDIV}-major-1kg.bias'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS_BIAS, '{INDIV}-per.bias'),
            INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_RESULTS_BIAS, '{INDIV}-per_h2h.bias'),
        #     INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_RESULTS_BIAS, '{INDIV}-' + POP_DIRNAME + '-gnomad.bias'),
        #     INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS_BIAS, '{INDIV}-' + POP_DIRNAME + '-1kg.bias'),
            INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_RESULTS_BIAS, '{INDIV}-vg_' + '{}.bias'.format(ALLELE_FREQ_FOR_VG)),
        #     INDIV = INDIV),
    output:
        tsv = os.path.join(DIR_RESULTS_BIAS, 'bias.tsv'),
        done = touch(temp(os.path.join(DIR, 'refbias.done')))
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
            list_indiv.append(basename[0])
            method = '-'.join(basename[1:])
            list_method.append(method[: method.rfind('.bias')])
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

rule find_strongly_biased_reads_refflow:
    input:
        list_path = os.path.join(DIR_SECOND_PASS, 'refflow-{}.paths'.format(POP_DIRNAME)),
        list_id = os.path.join(DIR_SECOND_PASS, 'refflow-{}.ids'.format(POP_DIRNAME)),
        bias = os.path.join(DIR_SECOND_PASS, 'refflow-{}.txt'.format(POP_DIRNAME)),
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = PREFIX_PER + '.vcf'
    output:
        tsv = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + '{}-above{}_or_below{}.reads.tsv'.format(
                POP_DIRNAME,
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )),
        reads = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + '{}-above{}_or_below{}.reads'.format(
                POP_DIRNAME,
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            ))
    run:
        range_bias = '0-{},{}-1'.format(0.5 - BIAS_TAIL_THRDS, 0.5 + BIAS_TAIL_THRDS)
        shell('{PYTHON} {DIR_SCRIPTS}/refbias/find_reads_given_HET.py \
            -s {input.list_path} -v {input.vcf} -f {input.bias} -id {input.list_id} -m {output.reads} -r {range_bias} -vc {input.vcf_complete} -indiv {wildcards.INDIV}')

rule find_unbiased_reads_per:
    input:
        list_path = os.path.join(DIR_FIRST_PASS, 'per-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'per-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'per-refbias.txt'),
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = PREFIX_PER + '.vcf'
    output:
        tsv = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per-between_45_55.reads.tsv'
            ),
        reads = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per-between_45_55.reads'
            )
    run:
        range_bias = '0.45-0.55'
        shell('{PYTHON} {DIR_SCRIPTS}/refbias/find_reads_given_HET.py \
            -s {input.list_path} -v {input.vcf} -f {input.bias} -id {input.list_id} -m {output.reads} -r {range_bias} --sample 0.01 -vc {input.vcf_complete} -indiv {wildcards.INDIV}')

rule find_strongly_biased_reads_per:
    input:
        list_path = os.path.join(DIR_FIRST_PASS, 'per-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'per-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'per-refbias.txt'),
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = PREFIX_PER + '.vcf'
    output:
        tsv = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per-above{}_or_below{}.reads.tsv'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )),
        reads = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per-above{}_or_below{}.reads'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            ))
    run:
        range_bias = '0-{},{}-1'.format(0.5 - BIAS_TAIL_THRDS, 0.5 + BIAS_TAIL_THRDS)
        shell('{PYTHON} {DIR_SCRIPTS}/refbias/find_reads_given_HET.py \
            -s {input.list_path} -v {input.vcf} -f {input.bias} -id {input.list_id} -m {output.reads} -r {range_bias} -vc {input.vcf_complete} -indiv {wildcards.INDIV}')

rule find_strongly_biased_reads_per_haploid:
    input:
        list_path = os.path.join(DIR_FIRST_PASS, 'per_h2h-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'per_h2h-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'per_h2h-refbias.txt'),
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = PREFIX_PER + '.vcf'
    output:
        tsv = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per_h2h-above{}_or_below{}.reads.tsv'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )),
        reads = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per_h2h-above{}_or_below{}.reads'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            ))
    run:
        range_bias = '0-{},{}-1'.format(0.5 - BIAS_TAIL_THRDS, 0.5 + BIAS_TAIL_THRDS)
        shell('{PYTHON} {DIR_SCRIPTS}/refbias/find_reads_given_HET.py \
            -s {input.list_path} -v {input.vcf} -f {input.bias} -id {input.list_id} -m {output.reads} -r {range_bias} -vc {input.vcf_complete} -indiv {wildcards.INDIV}')

rule find_strongly_biased_reads_grc:
    input:
        list_path = os.path.join(DIR_FIRST_PASS, 'grc-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'grc-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'grc-refbias.txt'),
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = PREFIX_PER + '.vcf'
    output:
        tsv = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'grc-above{}_or_below{}.reads.tsv'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )),
        reads = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'grc-above{}_or_below{}.reads'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            ))
    run:
        range_bias = '0-{},{}-1'.format(0.5 - BIAS_TAIL_THRDS, 0.5 + BIAS_TAIL_THRDS)
        shell('{PYTHON} {DIR_SCRIPTS}/refbias/find_reads_given_HET.py \
            -s {input.list_path} -v {input.vcf} -f {input.bias} -id {input.list_id} -m {output.reads} -r {range_bias} -vc {input.vcf_complete} -indiv {wildcards.INDIV}')

rule find_strongly_biased_reads_major:
    input:
        list_path = os.path.join(DIR_FIRST_PASS, 'major-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'major-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'major-refbias.txt'),
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = PREFIX_PER + '.vcf'
    output:
        tsv = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'major-above{}_or_below{}.reads.tsv'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )),
        reads = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'major-above{}_or_below{}.reads'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            ))
    run:
        range_bias = '0-{},{}-1'.format(0.5 - BIAS_TAIL_THRDS, 0.5 + BIAS_TAIL_THRDS)
        shell('{PYTHON} {DIR_SCRIPTS}/refbias/find_reads_given_HET.py \
            -s {input.list_path} -v {input.vcf} -f {input.bias} -id {input.list_id} -m {output.reads} -r {range_bias} -vc {input.vcf_complete} -indiv {wildcards.INDIV}')

rule find_strongly_biased_reads_vg:
    input:
        list_path = os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.paths'.format(ALLELE_FREQ_FOR_VG)),
        list_id = os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.ids'.format(ALLELE_FREQ_FOR_VG)),
        bias = os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.txt'.format(ALLELE_FREQ_FOR_VG)),
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = PREFIX_PER + '.vcf'
    output:
        tsv = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'vg_{}-above{}_or_below{}.reads.tsv'.format(
                ALLELE_FREQ_FOR_VG,
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )),
        reads = os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'vg_{}-above{}_or_below{}.reads'.format(
                ALLELE_FREQ_FOR_VG,
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            ))
    run:
        range_bias = '0-{},{}-1'.format(0.5 - BIAS_TAIL_THRDS, 0.5 + BIAS_TAIL_THRDS)
        shell('{PYTHON} {DIR_SCRIPTS}/refbias/find_reads_given_HET.py \
            -s {input.list_path} -v {input.vcf} -f {input.bias} -id {input.list_id} -m {output.reads} -r {range_bias} -vc {input.vcf_complete} -indiv {wildcards.INDIV}')

rule check_find_reads:
    input:
        expand(os.path.join(DIR_RESULTS_BIAS,
            '{INDIV}-' + '{}-above{}_or_below{}.reads.tsv'.format(
            POP_DIRNAME, int(100*(0.5+BIAS_TAIL_THRDS)), int(100*(0.5-BIAS_TAIL_THRDS)))),
            INDIV = INDIV),
        expand(os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'grc-above{}_or_below{}.reads.tsv'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )), INDIV = INDIV),
        expand(os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'major-above{}_or_below{}.reads.tsv'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )), INDIV = INDIV),
        expand(os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per-above{}_or_below{}.reads.tsv'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )), INDIV = INDIV),
        expand(os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per_h2h-above{}_or_below{}.reads.tsv'.format(
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )), INDIV = INDIV),
        expand(os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'vg_{}-above{}_or_below{}.reads.tsv'.format(
                ALLELE_FREQ_FOR_VG,
                int(100*(0.5+BIAS_TAIL_THRDS)),
                int(100*(0.5-BIAS_TAIL_THRDS))
            )), INDIV = INDIV),
        expand(os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per-between_45_55.reads.tsv'
            ), INDIV = INDIV)
    output:
        done = touch(temp(os.path.join(DIR, 'find_biased_reads.done')))
