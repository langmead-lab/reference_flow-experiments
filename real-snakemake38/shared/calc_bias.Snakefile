rule find_reads_overlapping_het_grc:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        bam = os.path.join(DIR_FIRST_PASS, 'wg-GRC-sorted.bam'),
        sam = os.path.join(DIR_FIRST_PASS, 'wg-GRC.sam')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-GRC-{INDIV}_het_no_overlaps.sam')
    shell:
        'module load bedtools;'
        'grep \'^@\' {input.sam} > {output.sam};'
        'bedtools intersect -a {input.bam} -b {input.vcf} | samtools view >> {output.sam}'

rule find_reads_overlapping_het_major:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        bam = os.path.join(DIR_FIRST_PASS, 'wg-major-liftover-sorted.bam'),
        sam = os.path.join(DIR_FIRST_PASS, 'wg-major.sam')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-major-{INDIV}_het_no_overlaps.sam')
    shell:
        'module load bedtools;'
        'grep \'^@\' {input.sam} > {output.sam};'
        'bedtools intersect -a {input.bam} -b {input.vcf} | samtools view >> {output.sam}'

rule calc_major_bias:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        # sam = os.path.join(DIR_FIRST_PASS, 'wg-major-liftover-sorted.sam')
        sam = os.path.join(DIR_FIRST_PASS, 'wg-major-{INDIV}_het_no_overlaps.sam')
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'major-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'major-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'major-refbias.txt')
    shell:
        'echo "major" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
            -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_grc_bias:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        # sam = os.path.join(DIR_FIRST_PASS, 'wg-GRC-sorted.sam')
        sam = os.path.join(DIR_FIRST_PASS, 'wg-GRC-{INDIV}_het_no_overlaps.sam')
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
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        A = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapA-liftover-sorted.sam'),
        B = os.path.join(DIR_FIRST_PASS, 'wg-per-merged-hapB-liftover-sorted.sam')
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'per-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'per-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'per-refbias.txt')
    shell:
        'echo "perA" > {output.list_id};'
        'echo "perB" >> {output.list_id};'
        'ls {input.A} > {output.list_path};'
        'ls {input.B} >> {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
            -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_reffllow_bias:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_SECOND_PASS,
            'wg-refflow-{}-{}-liftover-sorted.sam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        list_path = os.path.join(DIR_SECOND_PASS, 'refflow-{}.paths'.format(POP_DIRNAME)),
        list_id = os.path.join(DIR_SECOND_PASS, 'refflow-{}.ids'.format(POP_DIRNAME)),
        bias = os.path.join(DIR_SECOND_PASS, 'refflow-{}.txt'.format(POP_DIRNAME))
    shell:
        'echo "refflow" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
            -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

# rule calc_refflow_bias:
#     input:
#         vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
#         maj = os.path.join(DIR_FIRST_PASS,
#             'wg-major-mapqgeq{}-liftover-sorted.sam'.format(ALN_MAPQ_THRSD)),
#         second_maj = os.path.join(DIR_SECOND_PASS, '2ndpass-major-liftover-sorted.sam'),
#         second_pop = [os.path.join(DIR_SECOND_PASS,'2ndpass-') + 
#             g + '-liftover-sorted.sam' for g in GROUP]
#     output:
#         list_path = os.path.join(DIR_SECOND_PASS, 'refflow-{}.paths'.format(POP_DIRNAME)),
#         list_id = os.path.join(DIR_SECOND_PASS, 'refflow-{}.ids'.format(POP_DIRNAME)),
#         bias = os.path.join(DIR_SECOND_PASS, 'refflow-{}.txt'.format(POP_DIRNAME))
#     run:
#         #: prepare list_path and list_id
#         second_pass_group = ['maj']
#         for g in GROUP:
#             second_pass_group.append(g)
#         with open(output.list_id, 'w') as f:
#             for g in second_pass_group:
#                 f.write(g + '\n')
#             f.write('maj-mapqgeq{}\n'.format(ALN_MAPQ_THRSD))
#         list_second_pass_lifted_sam = [
#             # os.path.join(DIR_SECOND_PASS, '2ndpass-') + g +
#             os.path.join(DIR, 'experiments/' + wildcards.INDIV +
#             '/' + POP_DIRNAME + '/2ndpass-') + g +
#             '-liftover-sorted.sam' for g in second_pass_group]
#         with open(output.list_path, 'w') as f:
#             for s in list_second_pass_lifted_sam:
#                 # sys.stderr.write(s + '\n')
#                 f.write(s + '\n')
#             f.write(input.maj + '\n')
#         shell('cat {output.list_path};')
#         shell('{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
#                -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}')

rule summarize_grc:
    input:
        os.path.join(DIR_FIRST_PASS,
            'grc-refbias.txt')
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-grc.bias')
    run:
        summarize_allelc_bias(input, output)

rule summarize_major:
    input:
        os.path.join(DIR_FIRST_PASS,
            'major-refbias.txt')
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-major.bias')
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

rule summarize_refflow:
    input:
        os.path.join(DIR_SECOND_PASS, 'refflow-{}.txt'.format(POP_DIRNAME))
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-' + POP_DIRNAME + '.bias')
    run:
        summarize_allelc_bias(input, output)

rule check_refbias_and_write_to_tsv:
    input:
        expand(
            os.path.join(DIR_RESULTS_BIAS, '{INDIV}-grc.bias'),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_RESULTS_BIAS, '{INDIV}-major.bias'),
            INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_RESULTS_BIAS, '{INDIV}-per.bias'),
        #     INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_RESULTS_BIAS, '{INDIV}-' + POP_DIRNAME + '.bias'),
        #     INDIV = INDIV)
        # expand(
        #     os.path.join(DIR_RESULTS_BIAS, '{INDIV}-per.bias'),
        #     INDIV = INDIV),
        # expand(
        #     os.path.join(DIR_RESULTS_BIAS, '{INDIV}-' + POP_DIRNAME + '.bias'),
        #     INDIV = INDIV)
    output:
        tsv = os.path.join(DIR_RESULTS_BIAS, 'bias.tsv'),
        done = touch(temp(os.path.join(DIR, 'allelic_bias.done')))
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
            list_method.append(basename[1][:basename[1].rfind('.bias')])
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
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = os.path.join(DIR, 'wg_filtered.vcf')
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
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = os.path.join(DIR, 'wg_filtered.vcf')
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
            -s {input.list_path} -v {input.vcf} -f {input.bias} -id {input.list_id} -m {output.reads} -r {range_bias} --sample 0.001 -vc {input.vcf_complete} -indiv {wildcards.INDIV}')

rule find_strongly_biased_reads_per:
    input:
        list_path = os.path.join(DIR_FIRST_PASS, 'per-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'per-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'per-refbias.txt'),
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = os.path.join(DIR, 'wg_filtered.vcf')
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

rule find_strongly_biased_reads_grc:
    input:
        list_path = os.path.join(DIR_FIRST_PASS, 'grc-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'grc-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'grc-refbias.txt'),
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = os.path.join(DIR, 'wg_filtered.vcf')
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
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        vcf_complete = os.path.join(DIR, 'wg_filtered.vcf')
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

# rule find_strongly_biased_reads_vg:
#     input:
#         list_path = os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.paths'.format(ALLELE_FREQ_FOR_VG)),
#         list_id = os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.ids'.format(ALLELE_FREQ_FOR_VG)),
#         bias = os.path.join(DIR_FIRST_PASS, 'vg_{}-refbias.txt'.format(ALLELE_FREQ_FOR_VG)),
#         vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
#         vcf_complete = os.path.join(DIR, 'wg_filtered.vcf')
#     output:
#         tsv = os.path.join(
#             DIR_RESULTS_BIAS,
#             '{INDIV}-' + 'vg_{}-above{}_or_below{}.reads.tsv'.format(
#                 ALLELE_FREQ_FOR_VG,
#                 int(100*(0.5+BIAS_TAIL_THRDS)),
#                 int(100*(0.5-BIAS_TAIL_THRDS))
#             )),
#         reads = os.path.join(
#             DIR_RESULTS_BIAS,
#             '{INDIV}-' + 'vg_{}-above{}_or_below{}.reads'.format(
#                 ALLELE_FREQ_FOR_VG,
#                 int(100*(0.5+BIAS_TAIL_THRDS)),
#                 int(100*(0.5-BIAS_TAIL_THRDS))
#             ))
#     run:
#         range_bias = '0-{},{}-1'.format(0.5 - BIAS_TAIL_THRDS, 0.5 + BIAS_TAIL_THRDS)
#         shell('{PYTHON} {DIR_SCRIPTS}/refbias/find_reads_given_HET.py \
#             -s {input.list_path} -v {input.vcf} -f {input.bias} -id {input.list_id} -m {output.reads} -r {range_bias} -vc {input.vcf_complete} -indiv {wildcards.INDIV}')

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
#         expand(os.path.join(
#             DIR_RESULTS_BIAS,
#             '{INDIV}-' + 'vg_{}-above{}_or_below{}.reads.tsv'.format(
#                 ALLELE_FREQ_FOR_VG,
#                 int(100*(0.5+BIAS_TAIL_THRDS)),
#                 int(100*(0.5-BIAS_TAIL_THRDS))
#             )), INDIV = INDIV),
        expand(os.path.join(
            DIR_RESULTS_BIAS,
            '{INDIV}-' + 'per-between_45_55.reads.tsv'
            ), INDIV = INDIV)
    output:
        done = touch(temp(os.path.join(DIR, 'find_biased_reads.done')))
