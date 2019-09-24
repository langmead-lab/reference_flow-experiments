''' Prepare lft files '''
rule liftover_serialize_major:
    input:
        vcf = os.path.join(DIR_MAJOR, CHROM + '_h37maj.vcf')
    output:
        lft = os.path.join(DIR_MAJOR, CHROM + '-h37maj.lft')
    params:
        os.path.join(DIR_MAJOR, CHROM + '-h37maj')
    shell:
        'module load gcc/5.5.0;'
        '{LIFTOVER} serialize -v {input.vcf} -p {params}'

rule liftover_serialize_per:
    input:
        vcf = PREFIX_VCF_F + '.vcf'
    output:
        lftA = os.path.join(DIR_PER, CHROM + '-perA.lft'),
        lftB = os.path.join(DIR_PER, CHROM + '-perB.lft')
    params:
        A = os.path.join(DIR_PER, CHROM + '-perA'),
        B = os.path.join(DIR_PER, CHROM + '-perB')
    shell:
        'module load gcc/5.5.0;'
        '{LIFTOVER} serialize -v {input.vcf} -p {params.A} -g 0 -s {wildcards.INDIV};'
        '{LIFTOVER} serialize -v {input.vcf} -p {params.B} -g 1 -s {wildcards.INDIV};'

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
        lft = os.path.join(DIR_MAJOR, '{}-h37maj.lft'.format(CHROM)),
        vcf = PREFIX_VCF_F + '.vcf'
    output:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover.sam'.format(CHROM))
    params:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover'.format(CHROM))
    run:
        shell('module load gcc/5.5.0;')
        try:
            shell('{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}')
        except:
            shell('{LIFTOVER} lift -a {input.sam} -v {input.vcf} -p {params} -s {wildcards.INDIV}')

rule liftover_lift_per:
    input:
        samA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA.sam'.format(CHROM)),
        lftA = os.path.join(DIR_PER, CHROM + '-perA.lft'),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB.sam'.format(CHROM)),
        lftB = os.path.join(DIR_PER, CHROM + '-perB.lft'),
        vcf = PREFIX_VCF_F + '.vcf'
    output:
        A = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover.sam'.format(CHROM))
    params:
        A = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover'.format(CHROM))
    run:
        shell('module load gcc/5.5.0;')
        try:
            shell('{LIFTOVER} lift -a {input.samA} -l {input.lftA} -p {params.A};')
        except:
            shell('{LIFTOVER} lift -a {input.samA} -v {input.vcf} -p {params.A} -g 0 -s {wildcards.INDIV};')
        try:
            shell('{LIFTOVER} lift -a {input.samB} -l {input.lftB} -p {params.B};')
        except:
            shell('{LIFTOVER} lift -a {input.samB} -v {input.vcf} -p {params.B} -g 1 -s {wildcards.INDIV};')

rule fix_liftover_per_header:
    input:
        liftA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover.sam'.format(CHROM)),
        liftB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover.sam'.format(CHROM)),
        samA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB.sam'.format(CHROM)),
    output:
        liftA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover_diploid.sam'.format(CHROM)),
        liftB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover_diploid.sam'.format(CHROM))
    run:
        header_A = ''
        with open(input.samA, 'r') as f:
            for line in f:
                if line.startswith('@SQ'):
                    header_A = line
                    break
        header_B = ''
        with open(input.samB, 'r') as f:
            for line in f:
                if line.startswith('@SQ'):
                    header_B = line
                    break
        f_liftA = open(input.liftA, 'r')
        with open(output.liftA, 'w') as f:
            for line in f_liftA:
                if line.startswith('@SQ'):
                    f.write(header_A)
                else:
                    f.write(line)
        f_liftB = open(input.liftB, 'r')
        with open(output.liftB, 'w') as f:
            for line in f_liftB:
                if line.startswith('@SQ'):
                    f.write(header_B)
                else:
                    f.write(line)
        f_liftA.close()
        f_liftB.close()

rule liftover_lift_per_haploid_setting:
    input:
        samA = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid.sam'.format(CHROM)),
        lftA = os.path.join(DIR_PER, CHROM + '-perA.lft'),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid.sam'.format(CHROM)),
        lftB = os.path.join(DIR_PER, CHROM + '-perB.lft'),
        vcf = PREFIX_VCF_F + '.vcf'
    output:
        A = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-liftover.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-liftover.sam'.format(CHROM))
    params:
        A = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-liftover'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-liftover'.format(CHROM))
    run:
        shell('module load gcc/5.5.0;')
        try:
            shell('{LIFTOVER} lift -a {input.samA} -l {input.lftA} -p {params.A};')
        except:
            shell('{LIFTOVER} lift -a {input.samA} -v {input.vcf} -p {params.A} -g 0 -s {wildcards.INDIV};')
        try:
            shell('{LIFTOVER} lift -a {input.samB} -l {input.lftB} -p {params.B};')
        except:
            shell('{LIFTOVER} lift -a {input.samB} -v {input.vcf} -p {params.B} -g 1 -s {wildcards.INDIV};')

rule fix_liftover_per_header_haploid_setting:
    input:
        liftA = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-liftover.sam'.format(CHROM)),
        liftB = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-liftover.sam'.format(CHROM)),
        samA = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid.sam'.format(CHROM)),
    output:
        liftA = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-liftover_diploid.sam'.format(CHROM)),
        liftB = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-liftover_diploid.sam'.format(CHROM))
    run:
        header_A = ''
        with open(input.samA, 'r') as f:
            for line in f:
                if line.startswith('@SQ'):
                    header_A = line
                    break
        header_B = ''
        with open(input.samB, 'r') as f:
            for line in f:
                if line.startswith('@SQ'):
                    header_B = line
                    break
        f_liftA = open(input.liftA, 'r')
        with open(output.liftA, 'w') as f:
            for line in f_liftA:
                if line.startswith('@SQ'):
                    f.write(header_A)
                else:
                    f.write(line)
        f_liftB = open(input.liftB, 'r')
        with open(output.liftB, 'w') as f:
            for line in f_liftB:
                if line.startswith('@SQ'):
                    f.write(header_B)
                else:
                    f.write(line)
        f_liftA.close()
        f_liftB.close()

#: Refflow -- first pass
rule liftover_lift_refflow_firstpass:
    input:
        sam = os.path.join(DIR_FIRST_PASS, '{}-h37maj-mapqgeq{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lft = os.path.join(DIR_MAJOR, '{}-h37maj.lft'.format(CHROM)),
        vcf = PREFIX_VCF_F + '.vcf'
    output:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-mapqgeq{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD))
    params:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-mapqgeq{}-liftover'.format(CHROM, ALN_MAPQ_THRSD))
    run:
        shell('module load gcc/5.5.0;')
        try:
            shell('{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}')
        except:
            shell('{LIFTOVER} lift -a {input.sam} -v {input.vcf} -p {params} -s {wildcards.INDIV};')

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

rule sort_lifted_per:
    input:
        # A = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover.sam'.format(CHROM)),
        # B = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover.sam'.format(CHROM))
        A = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover_diploid.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover_diploid.sam'.format(CHROM))
    output:
        # os.path.join(DIR_FIRST_PASS, '{}-per-liftover-sorted.sam'.format(CHROM))
        A = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover-sorted.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover-sorted.sam'.format(CHROM))
    shell:
        'samtools sort -@ {THREADS} -o {output.A} {input.A};'
        'samtools sort -@ {THREADS} -o {output.B} {input.B}'

rule sort_lifted_per_haploid_setting:
    input:
        A = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-liftover_diploid.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-liftover_diploid.sam'.format(CHROM))
    output:
        A = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-liftover-sorted.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-liftover-sorted.sam'.format(CHROM))
    shell:
        'samtools sort -@ {THREADS} -o {output.A} {input.A};'
        'samtools sort -@ {THREADS} -o {output.B} {input.B}'

rule sort_grc:
    input:
        os.path.join(DIR_FIRST_PASS, '{}-grch37.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-grch37-sorted.sam'.format(CHROM))
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
        perA = expand(
            os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover.sam'.format(CHROM)),
            INDIV = INDIV),
        perB = expand(
            os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover.sam'.format(CHROM)),
            INDIV = INDIV),
        perAh = expand(
            os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-liftover.sam'.format(CHROM)),
            INDIV = INDIV),
        perBh = expand(
            os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-liftover.sam'.format(CHROM)),
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
        maj = expand(os.path.join(DIR_FIRST_PASS,
            '{}-h37maj-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        #: per
        # per = expand(os.path.join(DIR_FIRST_PASS,
        #     '{}-per-liftover-sorted.sam'.format(CHROM)),
        #     INDIV = INDIV),
        perA = expand(os.path.join(DIR_FIRST_PASS,
            '{}-per-merged-hapA-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        perB = expand(os.path.join(DIR_FIRST_PASS,
            '{}-per-merged-hapB-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        perAh = expand(os.path.join(DIR_FIRST_PASS,
            '{}-per_hapA_haploid-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        perBh = expand(os.path.join(DIR_FIRST_PASS,
            '{}-per_hapB_haploid-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        #: grch37
        grc = expand(os.path.join(DIR_FIRST_PASS,
            '{}-grch37-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        #: refflow
        pop_maj = expand(os.path.join(
            DIR_SECOND_PASS, '2ndpass-h37maj-liftover-sorted.sam'),
            INDIV = INDIV),
        fp = expand(os.path.join(DIR_FIRST_PASS,
            '{}-h37maj-mapqgeq{}-liftover-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD)),
            INDIV = INDIV),
        pop = expand(os.path.join(
            DIR_SECOND_PASS, '2ndpass-{GROUP}-liftover-sorted.sam'),
            INDIV = INDIV, GROUP = GROUP)
    output:
        touch(temp(os.path.join(DIR, 'sorting.done')))

