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
    shell:
        'bowtie2-build --threads {THREADS} {input.major} {PREFIX_MAJOR_IDX}'

rule check_major:
    input:
        expand(PREFIX_MAJOR_IDX + '.{idx_item}.bt2', idx_item = IDX_ITEMS)
    output:
        touch(temp(os.path.join(DIR, 'major.done')))
