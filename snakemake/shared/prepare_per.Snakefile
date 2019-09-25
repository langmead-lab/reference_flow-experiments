'''
Rules for building personalized genome
'''
rule build_per:
    input:
        genome = GENOME,
        vcf = PREFIX_VCF_F + '.vcf'
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
        DIR_PER_IDX + CHROM + '-per_hapA.1.bt2',
        DIR_PER_IDX + CHROM + '-per_hapA.2.bt2',
        DIR_PER_IDX + CHROM + '-per_hapA.3.bt2',
        DIR_PER_IDX + CHROM + '-per_hapA.4.bt2',
        DIR_PER_IDX + CHROM + '-per_hapA.rev.1.bt2',
        DIR_PER_IDX + CHROM + '-per_hapA.rev.2.bt2',
        DIR_PER_IDX + CHROM + '-per_hapB.1.bt2',
        DIR_PER_IDX + CHROM + '-per_hapB.2.bt2',
        DIR_PER_IDX + CHROM + '-per_hapB.3.bt2',
        DIR_PER_IDX + CHROM + '-per_hapB.4.bt2',
        DIR_PER_IDX + CHROM + '-per_hapB.rev.1.bt2',
        DIR_PER_IDX + CHROM + '-per_hapB.rev.2.bt2'
    params:
        prefix_idxA = DIR_PER_IDX + CHROM + '-per_hapA',
        prefix_idxB = DIR_PER_IDX + CHROM + '-per_hapB'
    threads: THREADS
    shell:
        'bowtie2-build --threads {THREADS} {input.perA} {params.prefix_idxA};'
        'bowtie2-build --threads {THREADS} {input.perB} {params.prefix_idxB}'

rule check_per:
    input:
        expand(
            DIR_PER_IDX + CHROM + '-per_hapA.{idx_item}.bt2',
            idx_item = IDX_ITEMS,
            INDIV = INDIV
        ),
        expand(
            DIR_PER_IDX + CHROM + '-per_hapB.{idx_item}.bt2',
            idx_item = IDX_ITEMS,
            INDIV = INDIV
        )
    output:
        touch(temp(os.path.join(DIR, 'personalization.done')))

