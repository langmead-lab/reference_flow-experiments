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
                elif line[2] == 'chr{}A'.format(CHROM):
                    line[2] = 'chr{}'.format(CHROM)
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
                elif line[2] == 'chr{}B'.format(CHROM):
                    line[2] = 'chr{}'.format(CHROM)
                fB_out.write('\t'.join(line) + '\n')

''' Prepare lft files '''
rule liftover_serialize_major_gnomad:
    input:
        vcf = os.path.join(DIR_MAJOR, 'chr{}-major-gnomad.vcf'.format(CHROM))
    output:
        lft = os.path.join(DIR_MAJOR, 'chr{}-major-gnomad.lft'.format(CHROM))
    params:
        os.path.join(DIR_MAJOR, 'chr{}-major-gnomad'.format(CHROM))
    shell:
        '{LIFTOVER} serialize -v {input.vcf} -p {params}'

rule liftover_serialize_major_onekg:
    input:
        vcf = os.path.join(DIR_MAJOR, 'chr{}-major-1kg.vcf'.format(CHROM)),
        chr_map = LENGTH_MAP
        # chr_map = os.path.join(DIR, 'GRCh38.length_map')
    output:
        lft = os.path.join(DIR_MAJOR, 'chr{}-major-1kg.lft'.format(CHROM))
    params:
        os.path.join(DIR_MAJOR, 'chr{}-major-1kg'.format(CHROM))
    shell:
        '{LIFTOVER} serialize -k {input.chr_map} -v {input.vcf} -p {params}'

rule liftover_serialize_per:
    input:
        vcf = PHASED_VCF_F,
        chr_map = LENGTH_MAP
        # chr_map = os.path.join(DIR, 'GRCh38.length_map')
    output:
        lftA = os.path.join(DIR_PER, 'chr{}-perA.lft'.format(CHROM)),
        lftB = os.path.join(DIR_PER, 'chr{}-perB.lft'.format(CHROM))
    params:
        A = os.path.join(DIR_PER, 'chr{}-perA'.format(CHROM)),
        B = os.path.join(DIR_PER, 'chr{}-perB'.format(CHROM))
    shell:
        '{LIFTOVER} serialize -k {input.chr_map} -v {input.vcf} -p {params.A} -g 0 -s {wildcards.INDIV};'
        '{LIFTOVER} serialize -k {input.chr_map} -v {input.vcf} -p {params.B} -g 1 -s {wildcards.INDIV};'

rule liftover_serialize_pop_genome:
    input:
        vcf = os.path.join(DIR_POP_GENOME, POP_DIRNAME + '/' +
            POP_GENOME_SUFFIX + '.vcf'),
            # CHROM + '_superpop_{GROUP}_' + POP_DIRNAME  + '.vcf')
        chr_map = LENGTH_MAP
        # chr_map = os.path.join(DIR, 'GRCh38.length_map')
    output:
        lft = os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            POP_GENOME_SUFFIX + '.lft')
            # CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.lft')
    params:
        os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            POP_GENOME_SUFFIX)
            # CHROM + '_superpop_{GROUP}_' + POP_DIRNAME)
    run:
        shell('{LIFTOVER} serialize -k {input.chr_map} -v {input.vcf} -p {params}')

''' Lifting SAMs '''
rule liftover_lift_major_gnomad:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad.sam'.format(CHROM)),
        lft = os.path.join(DIR_MAJOR, 'chr{}-major-gnomad.lft'.format(CHROM)),
        vcf = PREFIX_MAJOR + '-gnomad.vcf'
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-liftover.sam'.format(CHROM))
    params:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-liftover'.format(CHROM))
    run:
        shell('{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}')

rule liftover_lift_major_onekg:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg.sam'.format(CHROM)),
        lft = os.path.join(DIR_MAJOR, 'chr{}-major-1kg.lft'.format(CHROM)),
        vcf = PREFIX_MAJOR + '-1kg.vcf'
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-liftover.sam'.format(CHROM))
    params:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-liftover'.format(CHROM))
    run:
        shell('{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}')

rule convert_chrom_for_per:
    input:
        sam_grc = os.path.join(DIR_FIRST_PASS, 'chr{}-grc.sam'.format(CHROM)),
        samA = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB.sam'.format(CHROM))
    output:
        samA = temp(os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA-converted.sam'.format(CHROM))),
        samB = temp(os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB-converted.sam'.format(CHROM)))
    run:
        convert_chrom(input.sam_grc, input.samA, input.samB, output.samA, output.samB)

rule liftover_lift_per:
    input:
        samA = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA-converted.sam'.format(CHROM)),
        lftA = os.path.join(DIR_PER, 'chr{}-perA.lft'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB-converted.sam'.format(CHROM)),
        lftB = os.path.join(DIR_PER, 'chr{}-perB.lft'.format(CHROM)),
        vcf = PHASED_VCF_F
    output:
        A = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA-liftover.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB-liftover.sam'.format(CHROM))
    params:
        A = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA-liftover'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB-liftover'.format(CHROM))
    run:
        shell('{LIFTOVER} lift -a {input.samA} -l {input.lftA} -p {params.A};')
        shell('{LIFTOVER} lift -a {input.samB} -l {input.lftB} -p {params.B};')

rule convert_chrom_for_per_haploid_setting:
    input:
        sam_grc = os.path.join(DIR_FIRST_PASS, 'chr{}-grc.sam'.format(CHROM)),
        samA = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapA_haploid.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapB_haploid.sam'.format(CHROM))
    output:
        samA = temp(os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapA_haploid-converted.sam'.format(CHROM))),
        samB = temp(os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapB_haploid-converted.sam'.format(CHROM)))
    run:
        convert_chrom(input.sam_grc, input.samA, input.samB, output.samA, output.samB)

rule liftover_lift_per_haploid_setting:
    input:
        samA = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapA_haploid-converted.sam'.format(CHROM)),
        lftA = os.path.join(DIR_PER, 'chr{}-perA.lft'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapB_haploid-converted.sam'.format(CHROM)),
        lftB = os.path.join(DIR_PER, 'chr{}-perB.lft'.format(CHROM)),
        vcf = PHASED_VCF_F
    output:
        A = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapA_haploid-liftover.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapB_haploid-liftover.sam'.format(CHROM))
    params:
        A = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapA_haploid-liftover'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapB_haploid-liftover'.format(CHROM))
    run:
        shell('{LIFTOVER} lift -a {input.samA} -l {input.lftA} -p {params.A};')
        shell('{LIFTOVER} lift -a {input.samB} -l {input.lftB} -p {params.B};')

#: Refflow -- first pass
rule liftover_lift_refflow_firstpass_gnomad:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-mapqgeq{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lft = os.path.join(DIR_MAJOR, 'chr{}-major-gnomad.lft'.format(CHROM)),
        vcf = PREFIX_MAJOR + '-gnomad.vcf'
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-mapqgeq{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD))
    params:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-mapqgeq{}-liftover'.format(CHROM, ALN_MAPQ_THRSD))
    run:
        shell('{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}')

#: Refflow -- first pass
rule liftover_lift_refflow_firstpass_onekg:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-mapqgeq{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lft = os.path.join(DIR_MAJOR, 'chr{}-major-1kg.lft'.format(CHROM)),
        vcf = PREFIX_MAJOR + '-1kg.vcf'
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-mapqgeq{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD))
    params:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-mapqgeq{}-liftover'.format(CHROM, ALN_MAPQ_THRSD))
    run:
        shell('{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}')

#: Refflow -- second pass
rule liftover_lift_refflow_secondpass_and_merge_gnomad:
    input:
        maj_fp = os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-mapqgeq{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        second_sam_path = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-gnomad.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        lft_pop = expand(os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            POP_GENOME_SUFFIX + '.lft'),
            # 'chr' + CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.lft'),
            GROUP = GROUP),
        lft_maj = os.path.join(DIR_MAJOR, 'chr{}-major-gnomad.lft'.format(CHROM))
    output:
        lfted_refflow_sam = os.path.join(DIR_SECOND_PASS,
            'chr{}-refflow-{}-{}-liftover-gnomad.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        lfted_major_second_sam = temp(os.path.join(DIR_SECOND_PASS, '2ndpass-gnomad-major-liftover.sam')),
        lfted_group_second_sam = temp([os.path.join(DIR_SECOND_PASS,'2ndpass-gnomad-') + 
            g + '-liftover.sam' for g in GROUP])
    run:
        list_sam = []
        list_group = []
        #: copy lifted first pass sam 
        shell('cp {input.maj_fp} {output.lfted_refflow_sam};')
        #: files should be 
        #: DIR + '/experiments/{INDIV}/{POP_DIRNAME}/2ndpass-gnomad-{}.sam'
        #: where g should be {GROUP} + 'major'
        with open(input.second_sam_path, 'r') as f:
            for line in f:
                list_sam.append(line.rstrip())
                bn = os.path.basename(line)
                split_bn = os.path.splitext(bn)
                list_group.append(split_bn[0].split('-')[-1])
        for i in range(len(list_sam)):
            sam = list_sam[i]
            #: need to specify indiv here manually since snakemake 
            #: does not recognize it if using DIR_SECOND_PASS
            prefix = os.path.join(DIR,
                'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME + 
                '/2ndpass-gnomad-{}-liftover'.format(list_group[i]))
            if list_group[i] == 'major':
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, input.lft_maj))
                shell('{LIFTOVER} lift -a {sam} -l {input.lft_maj} -p {prefix};')
                #: append reads to all-in-one lifted SAM
                shell('grep -hv "^@" {prefix}.sam >> {output.lfted_refflow_sam};')
            elif list_group[i] in GROUP:
                for lft in input.lft_pop:
                    pop = os.path.basename(lft)
                    if lft.count(list_group[i]) > 0:
                        break
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, lft))
                shell('{LIFTOVER} lift -a {sam} -l {lft} -p {prefix};')
                #: append reads to all-in-one lifted SAM
                shell('grep -hv "^@" {prefix}.sam >> {output.lfted_refflow_sam};')

#: Refflow -- second pass
rule liftover_lift_refflow_secondpass_and_merge_onekg:
    input:
        maj_fp = os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-mapqgeq{}-liftover.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        second_sam_path = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-1kg.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        lft_pop = expand(os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            POP_GENOME_SUFFIX + '.lft'),
            # 'chr' + CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.lft'),
            GROUP = GROUP),
        lft_maj = os.path.join(DIR_MAJOR, 'chr{}-major-1kg.lft'.format(CHROM))
    output:
        lfted_refflow_sam = os.path.join(DIR_SECOND_PASS,
            'chr{}-refflow-{}-{}-liftover-1kg.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        lfted_major_second_sam = temp(os.path.join(DIR_SECOND_PASS, '2ndpass-1kg-major-liftover.sam')),
        lfted_group_second_sam = temp([os.path.join(DIR_SECOND_PASS,'2ndpass-1kg-') + 
            g + '-liftover.sam' for g in GROUP])
    run:
        list_sam = []
        list_group = []
        #: copy lifted first pass sam 
        shell('cp {input.maj_fp} {output.lfted_refflow_sam};')
        #: files should be 
        #: DIR + '/experiments/{INDIV}/{POP_DIRNAME}/2ndpass-1kg-{}.sam'
        #: where g should be {GROUP} + 'major'
        with open(input.second_sam_path, 'r') as f:
            for line in f:
                list_sam.append(line.rstrip())
                bn = os.path.basename(line)
                split_bn = os.path.splitext(bn)
                list_group.append(split_bn[0].split('-')[-1])
        for i in range(len(list_sam)):
            sam = list_sam[i]
            #: need to specify indiv here manually since snakemake 
            #: does not recognize it if using DIR_SECOND_PASS
            prefix = os.path.join(DIR,
                'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME + 
                '/2ndpass-1kg-{}-liftover'.format(list_group[i]))
            if list_group[i] == 'major':
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, input.lft_maj))
                shell('{LIFTOVER} lift -a {sam} -l {input.lft_maj} -p {prefix};')
                #: append reads to all-in-one lifted SAM
                shell('grep -hv "^@" {prefix}.sam >> {output.lfted_refflow_sam};')
            elif list_group[i] in GROUP:
                for lft in input.lft_pop:
                    pop = os.path.basename(lft)
                    if lft.count(list_group[i]) > 0:
                        break
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, lft))
                shell('{LIFTOVER} lift -a {sam} -l {lft} -p {prefix};')
                #: append reads to all-in-one lifted SAM
                shell('grep -hv "^@" {prefix}.sam >> {output.lfted_refflow_sam};')

rule liftover_popspecific_onepass:
    input:
        sam = os.path.join(DIR_FIRST_PASS,
            'chr{}'.format(CHROM) + '-{GROUP}-' + POP_DIRNAME +'.sam'),
        lft = os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' + 
            'chr{}'.format(CHROM) + '_superpop_{GROUP}_' + POP_DIRNAME + '.lft'),
        vcf = os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            'chr{}'.format(CHROM) + '_superpop_{GROUP}_' + POP_DIRNAME  + '.vcf')
    output:
        sam = temp(os.path.join(DIR_FIRST_PASS,
            'chr{}'.format(CHROM) + '-{GROUP}-' + POP_DIRNAME +'-liftover.sam'))
    params:
        os.path.join(DIR_FIRST_PASS,
            'chr{}'.format(CHROM) + '-{GROUP}-' + POP_DIRNAME +'-liftover')
    run:
        shell('{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}')

''' If there are multiple files, merge them after liftover ''' 
rule merge_per_allinone:
    input:
        A = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA-liftover.sam'.format(CHROM)),
        B = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-liftover.sam'.format(CHROM))
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
rule sort_lifted_major_gnomad:
    input:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-liftover-sorted.sam'.format(CHROM))
    threads: 2
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

rule sort_lifted_major_onekg:
    input:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-liftover-sorted.sam'.format(CHROM))
    threads: 2
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

rule sort_lifted_per:
    input:
        os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-liftover.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-liftover-sorted.sam'.format(CHROM))
    threads: 2
    shell:
        'samtools sort -@ {threads} -o {output} {input};'

rule sort_lifted_onepass_popspecific:
    input:
        os.path.join(DIR_FIRST_PASS, 'chr{}'.format(CHROM) + '-{GROUP}-' + POP_DIRNAME +'-liftover.sam'),
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}'.format(CHROM) + '-{GROUP}-' + POP_DIRNAME +'-liftover-sorted.sam'),
    threads: 2
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

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
        os.path.join(DIR_FIRST_PASS, 'chr{}-grc.sam'.format(CHROM))
    output:
        os.path.join(DIR_FIRST_PASS, 'chr{}-grc-sorted.sam'.format(CHROM))
    threads: THREADS
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

rule sort_refflow_allinone_gnomad:
    input:
        os.path.join(DIR_SECOND_PASS, 'chr{}-refflow-{}-{}-liftover-gnomad.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
    output:
        os.path.join(DIR_SECOND_PASS, 'chr{}-refflow-{}-{}-liftover-gnomad-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
    threads: 2
    run:
        shell('samtools sort -@ {threads} -o {output} {input};')

rule sort_refflow_allinone_onekg:
    input:
        os.path.join(DIR_SECOND_PASS, 'chr{}-refflow-{}-{}-liftover-1kg.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
    output:
        os.path.join(DIR_SECOND_PASS, 'chr{}-refflow-{}-{}-liftover-1kg-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
    threads: 2
    run:
        shell('samtools sort -@ {threads} -o {output} {input};')

''' Check points '''
rule check_liftover:
    input:
        # major_genomad = expand(
        #     os.path.join(DIR_FIRST_PASS, 'chr{}-major-gnomad-liftover.sam'.format(CHROM)),
        #     INDIV = INDIV),
        major_1kg = expand(
            os.path.join(DIR_FIRST_PASS, 'chr{}-major-1kg-liftover.sam'.format(CHROM)),
            INDIV = INDIV),
        per = expand(
            os.path.join(DIR_FIRST_PASS,
            'chr{}-per-merged-liftover.sam'.format(CHROM)),
            INDIV = INDIV
            ),
        # per_h2h = expand(
        #     os.path.join(DIR_FIRST_PASS,
        #     '{}-per_h2h-merged-liftover.sam'.format(CHROM)),
        #     INDIV = INDIV),
        # refflow_gnomad = expand(
        #     os.path.join(DIR_SECOND_PASS,
        #     'chr{}-refflow-{}-{}-liftover-gnomad.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        #     INDIV = INDIV),
        refflow_onekg = expand(
            os.path.join(DIR_SECOND_PASS,
            'chr{}-refflow-{}-{}-liftover-1kg.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
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
        # maj_gnomad = expand(os.path.join(DIR_FIRST_PASS,
        #     'chr{}-major-gnomad-liftover-sorted.sam'.format(CHROM)),
        #     INDIV = INDIV),
        maj_1kg = expand(os.path.join(DIR_FIRST_PASS,
            'chr{}-major-1kg-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        #: per
        per = expand(
            os.path.join(DIR_FIRST_PASS,
            'chr{}-per-merged-liftover-sorted.sam'.format(CHROM)),
            INDIV = INDIV
            ),
        # per_h2h = expand(
        #     os.path.join(DIR_FIRST_PASS,
        #     '{}-per_h2h-merged-liftover-sorted.sam'.format(CHROM)),
        #     INDIV = INDIV),
        #: grc
        grc = expand(os.path.join(DIR_FIRST_PASS,
            'chr{}-grc-sorted.sam'.format(CHROM)),
            INDIV = INDIV),
        #: refflow
        # refflow_gnomad = expand(
        #     os.path.join(DIR_SECOND_PASS,
        #     'chr{}-refflow-{}-{}-liftover-gnomad-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        #     INDIV = INDIV),
        onepass_popgenome = expand(
            os.path.join(DIR_FIRST_PASS,
            'chr{}'.format(CHROM) + '-{GROUP}-' + POP_DIRNAME +'-liftover-sorted.sam'),
        GROUP = GROUP, INDIV = INDIV),
        refflow_1kg = expand(
            os.path.join(DIR_SECOND_PASS,
            'chr{}-refflow-{}-{}-liftover-1kg-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'sorting.done')))

