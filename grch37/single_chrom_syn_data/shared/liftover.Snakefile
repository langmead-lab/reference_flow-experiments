''' Local helper functions '''
def convert_chrom(fn_grc, fn_samA, fn_samB, fn_outA, fn_outB):
    with open(fn_grc, 'r') as f:
        for line in f:
            if line.startswith('@SQ'):
                header = line
                break
    fA_out = open(fn_outA, 'w')
    with open(fn_samA, 'r') as fA:
        for line in fA:
            if line.startswith('@SQ'):
                fA_out.write(header)
            elif line[0] == '@':
                fA_out.write(line)
            else:
                line = line.split()
                if line[2] == CHROM + 'A':
                    line[2] = CHROM
                fA_out.write('\t'.join(line) + '\n')
    fB_out = open(fn_outB, 'w')
    with open(fn_samB, 'r') as fB:
        for line in fB:
            if line.startswith('@SQ'):
                fB_out.write(header)
            elif line[0] == '@':
                fB_out.write(line)
            else:
                line = line.split()
                if line[2] == CHROM + 'B':
                    line[2] = CHROM
                fB_out.write('\t'.join(line) + '\n')

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

rule convert_chrom_for_per:
    input:
        sam_grc = os.path.join(DIR_FIRST_PASS, '{}-grch37.sam'.format(CHROM)),
        samA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB.sam'.format(CHROM))
    output:
        samA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-converted.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-converted.sam'.format(CHROM))
    run:
        convert_chrom(input.sam_grc, input.samA, input.samB, output.samA, output.samB)

rule liftover_lift_per:
    input:
        samA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-converted.sam'.format(CHROM)),
        # samA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA.sam'.format(CHROM)),
        lftA = os.path.join(DIR_PER, CHROM + '-perA.lft'),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-converted.sam'.format(CHROM)),
        # samB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB.sam'.format(CHROM)),
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

rule convert_chrom_for_per_haploid_setting:
    input:
        sam_grc = os.path.join(DIR_FIRST_PASS, '{}-grch37.sam'.format(CHROM)),
        samA = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid.sam'.format(CHROM))
    output:
        samA = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-converted.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-converted.sam'.format(CHROM))
    run:
        convert_chrom(input.sam_grc, input.samA, input.samB, output.samA, output.samB)

rule liftover_lift_per_haploid_setting:
    input:
        samA = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-converted.sam'.format(CHROM)),
        lftA = os.path.join(DIR_PER, CHROM + '-perA.lft'),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-converted.sam'.format(CHROM)),
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

rule liftover_refflow_onepass:
    input:
        sam = os.path.join(DIR_FIRST_PASS, CHROM + '-{GROUP}-' + POP_DIRNAME +'.sam'),
        lft = os.path.join(
            DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.lft'),
        vcf = os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            CHROM + '_superpop_{GROUP}_' + POP_DIRNAME  + '.vcf')
    output:
        sam = os.path.join(DIR_FIRST_PASS, CHROM + '-{GROUP}-' + POP_DIRNAME +'-liftover.sam'),
    params:
        os.path.join(DIR_FIRST_PASS, CHROM + '-{GROUP}-' + POP_DIRNAME +'-liftover')
    run:
        shell('module load gcc/5.5.0;')
        try:
            shell('{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}')
        except:
            shell('{LIFTOVER} lift -a {input.sam} -v {input.vcf} -p {params}')

''' If there are multiple files, merge them after liftover ''' 
rule merge_refflow_allinone:
    input:
        maj_fp = os.path.join(DIR_FIRST_PASS, '{}-h37maj-mapqgeq{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        maj_sp = os.path.join(DIR_SECOND_PASS, '2ndpass-h37maj-liftover.sam'),
        pop = [os.path.join(DIR_SECOND_PASS,'2ndpass-') + 
            g + '-liftover.sam' for g in GROUP]
    output:
        os.path.join(DIR_SECOND_PASS, '{}-refflow-{}-{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
    shell:
        'cp {input.maj_fp} {output};'
        'grep -hv "^@" {input.maj_sp} >> {output};'
        'grep -hv "^@" {input.pop} >> {output};'

rule merge_per_allinone:
    input:
        A = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-per-merged-liftover.sam'.format(CHROM))
    shell:
        'cp {input.A} {output};'
        'grep -hv "^@" {input.B} >> {output}'

rule merge_per_haptohap_allinone:
    input:
        A = os.path.join(DIR_FIRST_PASS, '{}-per_hapA_haploid-liftover.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, '{}-per_hapB_haploid-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-per_h2h-merged-liftover.sam'.format(CHROM))
    shell:
        'cp {input.A} {output};'
        'grep -hv "^@" {input.B} >> {output}'

''' Sort using samtools '''
rule sort_lifted_major:
    input:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover-sorted.sam'.format(CHROM))
    threads: 2
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

rule sort_lifted_per:
    input:
        os.path.join(DIR_FIRST_PASS, '{}-per-merged-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-per-merged-liftover-sorted.sam'.format(CHROM))
    threads: 2
    shell:
        'samtools sort -@ {threads} -o {output} {input};'

rule sort_lifted_per_haploid_setting:
    input:
        os.path.join(DIR_FIRST_PASS, '{}-per_h2h-merged-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-per_h2h-merged-liftover-sorted.sam'.format(CHROM))
    threads: 2
    shell:
        'samtools sort -@ {threads} -o {output} {input};'

rule sort_grc:
    input:
        os.path.join(DIR_FIRST_PASS, '{}-grch37.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, '{}-grch37-sorted.sam'.format(CHROM))
    threads: THREADS
    shell:
        'samtools sort -@ {THREADS} -o {output} {input}'

rule sort_refflow_allinone:
    input:
        os.path.join(DIR_SECOND_PASS, '{}-refflow-{}-{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
    output:
        os.path.join(DIR_SECOND_PASS, '{}-refflow-{}-{}-liftover-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
    threads: 2
    run:
        shell('samtools sort -@ {threads} -o {output} {input};')

# rule sort_lifted_refflow_firstpass:
#     input:
#         os.path.join(DIR_FIRST_PASS,
#             '{}-h37maj-mapqgeq{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD))
#     output:
#         os.path.join(DIR_FIRST_PASS,
#             '{}-h37maj-mapqgeq{}-liftover-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD))
#     threads: 2
#     run:
#         shell('samtools sort -@ {threads} -o {output} {input};')
# 
# rule sort_lifted_refflow_secondpass:
#     input:
#         maj = os.path.join(DIR_SECOND_PASS, '2ndpass-h37maj-liftover.sam'),
#         pop = [os.path.join(DIR_SECOND_PASS,'2ndpass-') + 
#             g + '-liftover.sam' for g in GROUP]
#     output:
#         maj = os.path.join(DIR_SECOND_PASS, '2ndpass-h37maj-liftover-sorted.sam'),
#         pop = [os.path.join(DIR_SECOND_PASS,'2ndpass-') + 
#             g + '-liftover-sorted.sam' for g in GROUP]
#     threads: 2
#     run:
#         shell('samtools sort -@ {threads} -o {output.maj} {input.maj};')
#         for sam in input.pop:
#             # fn_out = os.path.basename(sam) + '-sorted.sam'
#             fn_out = sam[:sam.find('.sam')] + '-sorted.sam'
#             shell('samtools sort -@ {threads} -o {fn_out} {sam};')

''' Check points '''
rule check_liftover:
    input:
        major = expand(
            os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover.sam'.format(CHROM)),
            INDIV = INDIV),
        per = expand(
            os.path.join(DIR_FIRST_PASS,
            '{}-per-merged-liftover.sam'.format(CHROM)),
            INDIV = INDIV
            ),
        per_h2h = expand(
            os.path.join(DIR_FIRST_PASS,
            '{}-per_h2h-merged-liftover.sam'.format(CHROM)),
            INDIV = INDIV),
        refflow = expand(
            os.path.join(DIR_SECOND_PASS,
            '{}-refflow-{}-{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV),
        # refflow_onepass = expand(
        #     os.path.join(DIR_FIRST_PASS,
        #     CHROM + '-{GROUP}-' + POP_DIRNAME +'-liftover.sam'),
        #     GROUP = GROUP, INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'liftover.done')))

rule check_sort:
    input:
        #: major
        maj = expand(os.path.join(DIR_FIRST_PASS,
            '{}-h37maj-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        #: per
        per = expand(
            os.path.join(DIR_FIRST_PASS,
            '{}-per-merged-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV
            ),
        per_h2h = expand(
            os.path.join(DIR_FIRST_PASS,
            '{}-per_h2h-merged-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        #: grch37
        grc = expand(os.path.join(DIR_FIRST_PASS,
            '{}-grch37-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        #: refflow
        refflow = expand(
            os.path.join(DIR_SECOND_PASS,
            '{}-refflow-{}-{}-liftover-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'sorting.done')))

