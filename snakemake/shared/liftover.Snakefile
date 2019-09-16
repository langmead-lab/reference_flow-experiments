''' Prepare lft files '''
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

''' Lifting SAMs '''
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

#: Refflow -- first pass
rule liftover_lift_refflow_firstpass:
    input:
        sam = os.path.join(DIR_FIRST_PASS, '{}-h37maj-mapqgeq{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lft = os.path.join(DIR_MAJOR, '{}-h37maj.lft'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-mapqgeq{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD))
    params:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-mapqgeq{}-liftover'.format(CHROM, ALN_MAPQ_THRSD))
    shell:
        'module load gcc/5.5.0;'
        '{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}'

#: Refflow -- second pass
rule liftover_lift_refflow_secondpass:
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

''' Sort using samtools '''
rule sort_lifted_major:
    input:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover-sorted.sam'.format(CHROM))
    shell:
        'samtools sort -@ {THREADS} -o {output} {input}'

rule sort_lifted_refflow_firstpass:
    input:
        os.path.join(DIR_FIRST_PASS,
            '{}-h37maj-mapqgeq{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD))
    output:
        os.path.join(DIR_FIRST_PASS,
            '{}-h37maj-mapqgeq{}-liftover-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD))
    run:
        shell('samtools sort -@ {THREADS} -o {output} {input};')

rule sort_lifted_refflow_secondpass:
    input:
        maj = os.path.join(DIR_SECOND_PASS, '2ndpass-h37maj-liftover.sam'),
        pop = [os.path.join(DIR_SECOND_PASS,'2ndpass-') + 
            g + '-liftover.sam' for g in GROUP]
    output:
        maj = os.path.join(DIR_SECOND_PASS, '2ndpass-h37maj-liftover-sorted.sam'),
        pop = [os.path.join(DIR_SECOND_PASS,'2ndpass-') + 
            g + '-liftover-sorted.sam' for g in GROUP]
    run:
        shell('samtools sort -@ {THREADS} -o {output.maj} {input.maj};')
        for sam in input.pop:
            # fn_out = os.path.basename(sam) + '-sorted.sam'
            fn_out = sam[:sam.find('.sam')] + '-sorted.sam'
            shell('samtools sort -@ {THREADS} -o {fn_out} {sam};')

''' Check points '''
rule check_liftover:
    input:
        major = expand(
            os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover.sam'.format(CHROM)),
            INDIV = INDIV),
        pop = expand(os.path.join(DIR_SECOND_PASS,
            '2ndpass-{g}-liftover.sam'), g = GROUP, INDIV = INDIV),
        pop_maj = expand(os.path.join(DIR_SECOND_PASS,
            '2ndpass-h37maj-liftover.sam'), INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'liftover.done')))

rule check_sort:
    input:
        #: major
        pop_maj = expand(os.path.join(
            DIR_SECOND_PASS, '2ndpass-h37maj-liftover-sorted.sam'),
            INDIV = INDIV),
        #: refflow
        fp = expand(os.path.join(DIR_FIRST_PASS,
            '{}-h37maj-mapqgeq{}-liftover-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD)),
            INDIV = INDIV),
        pop = expand(os.path.join(
            DIR_SECOND_PASS, '2ndpass-{GROUP}-liftover-sorted.sam'),
            INDIV = INDIV, GROUP = GROUP),
        maj = expand(os.path.join(DIR_FIRST_PASS,
            '{}-h37maj-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'sorting.done')))

