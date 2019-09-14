rule liftover_serialize_major:
    input:
        vcf_major = os.path.join(DIR_MAJOR, CHROM + '_h37maj.vcf')
    output:
        lft = os.path.join(DIR_MAJOR, CHROM + '-h37maj.lft')
    params:
        os.path.join(DIR_MAJOR, CHROM + '-h37maj')
    shell:
        'module load gcc/5.5.0;'
        '{LIFTOVER} serialize -v {input.vcf_major} -p {params}'

rule liftover_lift_major:
    input:
        sam = os.path.join(DIR_FIRST_PASS, '{}-h37maj.sam'.format(CHROM)),
        lft = os.path.join(DIR_MAJOR, '{}-h37maj.lft'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover.sam'.format(CHROM))
    params:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover'.format(CHROM))
    shell:
        'module load gcc/5.5.0;'
        '{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}'

rule liftover_serialize_pop_genome:
    input:
        vcf = os.path.join(DIR_POP_GENOME, POP_DIRNAME + '/' +
            CHROM + '_superpop_{GROUP}_' + POP_DIRNAME  + '.vcf')
    output:
        lft = os.path.join(
            DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.lft')
    params:
        os.path.join(
            DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}_' + POP_DIRNAME)
    run:
        shell('module load gcc/5.5.0;')
        shell('{LIFTOVER} serialize -v {input.vcf} -p {params}')
        # dict_indiv_to_pop = build_dict_indiv_to_pop(FAMILY)
        # dict_pop_to_spop = build_dict_pop_to_spop(SPOP)
        # for i, g in enumerate(GROUP):
        #     if dict_pop_to_spop[dict_indiv_to_pop[wildcards.INDIV]] == g:
        #         vcf = input.vcf[i]
        #         sys.stderr.write('{} in {}, to be lifted with {}\n'.format(wildcards.INDIV, g, vcf))
        #         shell('{LIFTOVER} serialize -v {vcf} -p {params}')
        #         break

rule liftover_lift_pop_genome:
    input:
        second_sam_path = os.path.join(DIR_SECOND_PASS, '{0}-h37maj-{1}-{2}.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        lft_pop = expand(os.path.join(
            DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.lft'),
            GROUP = GROUP),
        vcf_pop = expand(os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            CHROM + '_superpop_{GROUP}_' + POP_DIRNAME  + '.vcf'),
            GROUP = GROUP),
        lft_maj = os.path.join(DIR_MAJOR, '{}-h37maj.lft'.format(CHROM))
    output:
        os.path.join(DIR_SECOND_PASS, '2ndpass-h37maj-liftover.sam'),
        [os.path.join(DIR_SECOND_PASS,'2ndpass-') + 
            g + '-liftover.sam' for g in GROUP]
    run:
        shell('module load gcc/5.5.0;')
        list_sam = []
        list_group = []
        #: files should be 
        #: DIR + '/experiments/{INDIV}/{POP_DIRNAME}/2ndpass-{}.sam'
        #: where g should be {GROUP} + 'h37maj'
        with open(input.second_sam_path, 'r') as f:
            for line in f:
                list_sam.append(line.rstrip())
                bn = os.path.basename(line)
                split_bn = os.path.splitext(bn)
                list_group.append(split_bn[0].split('-')[-1])
        for i in range(len(list_sam)):
            sam = list_sam[i]
            prefix = os.path.join(DIR,
                'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME + 
                '/2ndpass-{}-liftover'.format(list_group[i]))
            if list_group[i] == 'h37maj':
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, input.lft_maj))
                shell('{LIFTOVER} lift -a {sam} -l {input.lft_maj} -p {prefix};')
            elif list_group[i] in GROUP:
                for lft in input.lft_pop:
                    pop = os.path.basename(lft)
                    if lft.count(list_group[i]) > 0:
                        break
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, lft))
                try:
                    shell('{LIFTOVER} lift -a {sam} -l {lft} -p {prefix};')
                except:
                    for vcf in input.vcf_pop:
                        vcf_base = os.path.basename(vcf)
                        if vcf_base.count(list_group[i]) > 0:
                            sys.stderr.write('Warning: fail to read from .lft, try building new lft...\n')
                            shell('{LIFTOVER} lift -a {sam} -v {vcf} -p {prefix}')
                            break
        # shell('{liftover} lift -a {input.sam} -l {input.lft} -p {params}')

rule check_liftover:
    input:
        major = expand(
            os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover.sam'.format(CHROM)),
            INDIV = INDIV),
        pop = expand(os.path.join(DIR_SECOND_PASS,
            '2ndpass-{g}-liftover.sam'), g = GROUP, INDIV = INDIV),
        pop_maj = expand(os.path.join(DIR_SECOND_PASS,
            '2ndpass-h37maj-liftover.sam'), INDIV = INDIV)
        # pop = expand(
        #     os.path.join(DIR_SECOND_PASS,
        #     CHROM + '_superpop_{GROUP}_' + '{}-liftover.sam'.format(POP_DIRNAME)),
        #     INDIV = INDIV, GROUP = GROUP)
    output:
        touch(temp(os.path.join(DIR, 'liftover.done')))
