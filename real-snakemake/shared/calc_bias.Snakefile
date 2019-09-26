rule calc_major_bias:
    input:
        vcf = os.path.join(DIR, 'wg_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_FIRST_PASS, 'wg-h37maj-liftover-sorted.sam')
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
        sam = os.path.join(DIR_FIRST_PASS, 'wg-GRCh37-sorted.sam')
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'grch37-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'grch37-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'grch37-refbias.txt')
    shell:
        'echo "grch37" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
           -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule summarize_grc:
    input:
        os.path.join(DIR_FIRST_PASS,
            'grch37-refbias.txt')
    output:
        os.path.join(DIR_RESULTS_BIAS, '{INDIV}-grch37.bias')
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

rule check_refbias_and_write_to_tsv:
    input:
        expand(
            os.path.join(DIR_RESULTS_BIAS, '{INDIV}-grch37.bias'),
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
        df.to_csv(output.tsv, sep = '\t', index = None)
