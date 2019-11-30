'''
Rules for building major-allele genome
'''
rule build_gnomad_major:
    input:
        genome = GENOME,
        vcf = UNPHASED_VCF_F
    output:
        out_genome = PREFIX_MAJOR + '-gnomad.fa',
        out_var = PREFIX_MAJOR + '-gnomad.var',
        out_vcf = PREFIX_MAJOR + '-gnomad.vcf'
    params:
        AN = DICT_GNOMAD_POP_SIZE['all'],
        OP = PREFIX_MAJOR + '-gnomad'
    shell:
        'echo "Size of population: {params.AN}";'
        '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
        '    --ref {input.genome} --vcf {input.vcf} '
        '    --chrom {CHROM} --out-prefix {params.OP} '
        '    -d gnomad '
        '    --gnomad-ac-field AC --gnomad-pop-count {params.AN} '
        '    --gnomad-af-th 0.5 '
        '    --include-indels'

rule build_gnomad_major_index:
    input:
        major = PREFIX_MAJOR + '-gnomad.fa'
    output:
        PREFIX_MAJOR_IDX + '-gnomad.1.bt2',
        PREFIX_MAJOR_IDX + '-gnomad.2.bt2',
        PREFIX_MAJOR_IDX + '-gnomad.3.bt2',
        PREFIX_MAJOR_IDX + '-gnomad.4.bt2',
        PREFIX_MAJOR_IDX + '-gnomad.rev.1.bt2',
        PREFIX_MAJOR_IDX + '-gnomad.rev.2.bt2'
    threads: THREADS
    params:
        PREFIX_MAJOR_IDX + '-gnomad'
    shell:
        'bowtie2-build --threads {threads} {input.major} {params}'

rule build_onekg_major:
    input:
        genome = GENOME,
        vcf = PHASED_VCF_F,
        chrom_map = CHROM_MAP
        # chrom_map = os.path.join(DIR, 'GRCh38.chrom_map')
    output:
        # chrom_map = os.path.join(DIR_MAJOR, 'chr{}.chrom_map'.format(CHROM)),
        vcf_major = PREFIX_MAJOR_F + '-1kg.vcf',
        out_genome = PREFIX_MAJOR + '-1kg.fa',
        out_var = PREFIX_MAJOR + '-1kg.var',
        out_vcf = PREFIX_MAJOR + '-1kg.vcf'
    params:
        AN = DICT_GNOMAD_POP_SIZE['all'],
        OP = PREFIX_MAJOR + '-1kg'
    shell:
        #: The 1KG hg38 call set uses chrN instead of N for chrom names
        #: so take additional steps to convert it
        # 'echo "{CHROM} chr{CHROM}" > {output.chrom_map};'
        '{BCFTOOLS} view -O v -q 0.5 {input.vcf} -e \'AF = 0.5\' '
        '-v snps,indels -m2 -M2 | {BCFTOOLS} annotate --rename-chrs {input.chrom_map}'
        '> {output.vcf_major};'
        #: if no need to convert chrom name
        # '{BCFTOOLS} view -O v -q 0.5 {input.vcf} -e \'AF = 0.5\''
        # ' -v snps,indels -m2 -M2 '
        # '> {output.vcf_major};'
        '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
        '    --ref {input.genome} --vcf {output.vcf_major} '
        '    --chrom {CHROM} --out-prefix {params.OP} '
        '    --include-indels'

rule build_onekg_major_index:
    input:
        major = PREFIX_MAJOR + '-1kg.fa'
    output:
        PREFIX_MAJOR_IDX + '-1kg.1.bt2',
        PREFIX_MAJOR_IDX + '-1kg.2.bt2',
        PREFIX_MAJOR_IDX + '-1kg.3.bt2',
        PREFIX_MAJOR_IDX + '-1kg.4.bt2',
        PREFIX_MAJOR_IDX + '-1kg.rev.1.bt2',
        PREFIX_MAJOR_IDX + '-1kg.rev.2.bt2'
    threads: THREADS
    params:
        PREFIX_MAJOR_IDX + '-1kg'
    shell:
        'bowtie2-build --threads {threads} {input.major} {params}'

# rule check_major:
#     input:
#         expand(PREFIX_MAJOR_IDX + '-1kg.{idx_item}.bt2', idx_item = IDX_ITEMS),
#         expand(PREFIX_MAJOR_IDX + '-gnomad.{idx_item}.bt2', idx_item = IDX_ITEMS)
#     output:
#         touch(temp(os.path.join(DIR, 'major.done')))

'''
Rules for indexing GRC genome
'''
rule build_grc_index:
    input:
        GENOME
    output:
        expand(os.path.join(
            DIR_GRC_IDX, 'chr{}'.format(CHROM) + '_grc.{i}.bt2'),
            i = IDX_ITEMS)
    params:
        os.path.join(DIR_GRC_IDX, 'chr{}_grc'.format(CHROM))
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params}'

# rule check_grc:
#     input:
#         expand(os.path.join(DIR_GRC_IDX, CHROM + '_grc.{i}.bt2'),
#             i = IDX_ITEMS)
#     output:
#         touch(temp(os.path.join(DIR, 'grc.done')))

'''
Rules for building personalized genome
'''
rule build_per:
    input:
        genome = GENOME,
        vcf = PHASED_VCF_F
    output:
        hapA = PREFIX_PER + '_hapA.fa',
        hapB = PREFIX_PER + '_hapB.fa',
        var = PREFIX_PER + '.var',
        vcf = PREFIX_PER + '.vcf'
    params:
        out_prefix = PREFIX_PER
    shell:
        '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
        '    --ref {input.genome} --vcf {input.vcf} --name {wildcards.INDIV}'
        '    --chrom {CHROM} --out-prefix {params.out_prefix} '
        '    --include-indels'

rule build_per_index:
    input:
        perA = PREFIX_PER + '_hapA.fa',
        perB = PREFIX_PER + '_hapB.fa'
    output:
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapA.1.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapA.2.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapA.3.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapA.4.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapA.rev.1.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapA.rev.2.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapB.1.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapB.2.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapB.3.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapB.4.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapB.rev.1.bt2',
        DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapB.rev.2.bt2'
    params:
        prefix_idxA = DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapA',
        prefix_idxB = DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapB'
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input.perA} {params.prefix_idxA};'
        'bowtie2-build --threads {threads} {input.perB} {params.prefix_idxB}'

rule check_prepare:
    input:
        expand(os.path.join(DIR_GRC_IDX, 'chr{}'.format(CHROM) + '_grc.{idx_item}.bt2'),
            idx_item = IDX_ITEMS),
        expand(PREFIX_MAJOR_IDX + '-1kg.{idx_item}.bt2',
            idx_item = IDX_ITEMS),
        # expand(PREFIX_MAJOR_IDX + '-gnomad.{idx_item}.bt2',
        #     idx_item = IDX_ITEMS),
        expand(
            DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapA.{idx_item}.bt2',
            idx_item = IDX_ITEMS,
            INDIV = INDIV
        ),
        expand(
            DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapB.{idx_item}.bt2',
            idx_item = IDX_ITEMS,
            INDIV = INDIV
        )
    output:
        touch(temp(os.path.join(DIR, 'prepare.done')))
