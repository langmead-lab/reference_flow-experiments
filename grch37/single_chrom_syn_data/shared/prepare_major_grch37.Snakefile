'''
Rules for building major-allele genome
'''
rule build_major:
    input:
        genome = GENOME,
        vcf = PREFIX_VCF_F + '.vcf'
    output:
        vcf_major = PREFIX_MAJOR_F + '.vcf',
        vcf_major_gz = PREFIX_MAJOR_F + '.vcf.gz',
        vcf_major_gz_csi = PREFIX_MAJOR_F + '.vcf.gz.csi',
        out_genome = PREFIX_MAJOR + '.fa',
        out_var = PREFIX_MAJOR + '.var',
        out_vcf = PREFIX_MAJOR + '.vcf'
    threads: THREADS
    shell:
        '{BCFTOOLS} view -O z --threads {THREADS} -q 0.5 {input.vcf} -e \'AF = 0.5\' -v snps,indels -m2 -M2 > '
        '{output.vcf_major_gz};'
        '{BCFTOOLS} index --threads {THREADS} {output.vcf_major_gz};'
        'bgzip -cd {output.vcf_major_gz} > {output.vcf_major};'
        '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
        '    --ref {input.genome} --vcf {output.vcf_major} '
        '    --chrom {CHROM} --out-prefix {PREFIX_MAJOR} '
        '    --include-indels'

rule build_major_index:
    input:
        major = PREFIX_MAJOR + '.fa'
    output:
        PREFIX_MAJOR_IDX + '.1.bt2',
        PREFIX_MAJOR_IDX + '.2.bt2',
        PREFIX_MAJOR_IDX + '.3.bt2',
        PREFIX_MAJOR_IDX + '.4.bt2',
        PREFIX_MAJOR_IDX + '.rev.1.bt2',
        PREFIX_MAJOR_IDX + '.rev.2.bt2'
    threads: THREADS
    shell:
        'bowtie2-build --threads {THREADS} {input.major} {PREFIX_MAJOR_IDX}'

rule check_major:
    input:
        expand(PREFIX_MAJOR_IDX + '.{idx_item}.bt2', idx_item = IDX_ITEMS)
    output:
        touch(temp(os.path.join(DIR, 'major.done')))

'''
Rules for indexing GRCh37
'''
rule build_grc_index:
    input:
        GENOME
    output:
        expand(os.path.join(DIR_GRCH37_IDX, CHROM + '_grch37.{i}.bt2'),
            i = IDX_ITEMS)
    params:
        DIR_GRCH37_IDX + CHROM + '_grch37'
    threads: THREADS
    shell:
        'bowtie2-build --threads {THREADS} {input} {params}'

rule check_grc:
    input:
        expand(os.path.join(DIR_GRCH37_IDX, CHROM + '_grch37.{i}.bt2'),
            i = IDX_ITEMS)
    output:
        touch(temp(os.path.join(DIR, 'grch37.done')))

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
