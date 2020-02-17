rule build_major:
    input:
        vcf = PREFIX_VCF_F + '.vcf',
        genome_chrom = os.path.join(DIR, 'chr{CHROM}.fa')
    output:
        vcf_major = PREFIX_MAJOR_F + '.vcf',
        # vcf_major_gz = PREFIX_MAJOR_F + '.vcf.gz',
        # vcf_major_gz_csi = PREFIX_MAJOR_F + '.vcf.gz.csi',
        out_genome = PREFIX_MAJOR + '.fa',
        out_var = PREFIX_MAJOR + '.var',
        out_vcf = PREFIX_MAJOR + '.vcf'
    params:
        # out_prefix = os.path.join(DIR, 'major/{CHROM}_maj')
        out_prefix = os.path.join(DIR, 'major/chr{CHROM}_maj')
    shell:
        # '{BCFTOOLS} view -O z -q 0.5 {input.vcf} -e \'AF = 0.5\' -v snps,indels -m2 -M2 > '
        # '{output.vcf_major_gz};'
        # '{BCFTOOLS} index {output.vcf_major_gz};'
        # 'bgzip -cd {output.vcf_major_gz} > {output.vcf_major};'
        '{BCFTOOLS} view -O v -q 0.5 {input.vcf} -e \'AF = 0.5\' -v snps,indels -m2 -M2 > '
        '{output.vcf_major};'
        '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
        '    --ref {input.genome_chrom} --vcf {output.vcf_major} '
        '    --chrom {wildcards.CHROM} --out-prefix {params.out_prefix} '
        '    --include-indels'

rule merge_major_fasta:
    input:
        expand(PREFIX_MAJOR + '.fa', CHROM = CHROM)
    output:
        os.path.join(DIR, 'major/wg-maj.fa')
    shell:
        'cat {input} >> {output}'

rule aggregate_major_vcf:
    input:
        vcf = expand(PREFIX_MAJOR + '.vcf', CHROM = CHROM)
    output:
        vcf = os.path.join(DIR, 'major/wg-maj.vcf')
    shell:
        '{BCFTOOLS} concat -o {output.vcf} {input.vcf}'

rule aggregate_per_vcf:
    input:
        vcf = expand(
            os.path.join(DIR_PER, '{INDIV}-chr{CHROM}-per.vcf'),
            CHROM = CHROM, INDIV = INDIV
        )
    output:
        vcf = os.path.join(DIR_PER, '{INDIV}-wg-per.vcf')
    shell:
        '{BCFTOOLS} concat -o {output.vcf} {input.vcf}'

rule build_per:
    input:
        genome = GENOME,
        vcf = PREFIX_VCF_F + '.vcf'
    output:
        hapA = os.path.join(DIR_PER, '{INDIV}-chr{CHROM}-per_hapA.fa'),
        hapB = os.path.join(DIR_PER, '{INDIV}-chr{CHROM}-per_hapB.fa'),
        var = os.path.join(DIR_PER, '{INDIV}-chr{CHROM}-per.var'),
        vcf = os.path.join(DIR_PER, '{INDIV}-chr{CHROM}-per.vcf')
    params:
        out_prefix = os.path.join(DIR_PER, '{INDIV}-chr{CHROM}-per')
    shell:
        '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
        '    --ref {input.genome} --vcf {input.vcf} --name {wildcards.INDIV}'
        '    --chrom {wildcards.CHROM} --out-prefix {params.out_prefix} '
        '    --include-indels'

'''
Merge fasta for all chromosomes into one whole genome fasta.
The A, B suffixes of chromsomes are removed before merging.
'''
rule merge_per_fasta:
    input:
        hapA = expand(os.path.join(DIR_PER, '{INDIV}-chr{CHROM}-per_hapA.fa'),
            INDIV = INDIV, CHROM = CHROM),
        hapB = expand(os.path.join(DIR_PER, '{INDIV}-chr{CHROM}-per_hapB.fa'),
            INDIV = INDIV, CHROM = CHROM)
    output:
        outA = os.path.join(DIR_PER, 'wg-perA.fa'),
        outB = os.path.join(DIR_PER, 'wg-perB.fa')
    run:
        list_hapA = []
        for h in input.hapA:
            if h.count(wildcards.INDIV) > 0:
                f_outA = open(output.outA, 'a')
                with open(h, 'r') as f:
                    #: remove the suffix of chromosome if it is like 1A, 10A
                    header = f.readline().split()
                    if header[0].endswith('A'):
                        header[0] = header[0][0:len(header[0])-1]
                    fixed = ' '.join(header) + '\n'
                    f_outA.write(fixed)
                    #: write the rest of lines
                    for line in f:
                        f_outA.write(line)
        f_outA.close()
        list_hapB = []
        for h in input.hapB:
            if h.count(wildcards.INDIV) > 0:
                f_outB = open(output.outB, 'a')
                with open(h, 'r') as f:
                    #: remove the suffix of chromosome if it is like 1B, 10B
                    header = f.readline().split()
                    if header[0].endswith('B'):
                        header[0] = header[0][0:len(header[0])-1]
                    fixed = ' '.join(header) + '\n'
                    f_outB.write(fixed)
                    #: write the rest of lines
                    for line in f:
                        f_outB.write(line)
        f_outB.close()
        # list_hapA = []
        # for h in input.hapA:
        #     if h.count(wildcards.INDIV) > 0:
        #         list_hapA.append(h)
        # shell('cat {list_hapA} >> {output.outA};')
        # list_hapB = []
        # for h in input.hapB:
        #     if h.count(wildcards.INDIV) > 0:
        #         list_hapB.append(h)
        # shell('cat {list_hapB} >> {output.outB}')

rule build_major_index:
    input:
        os.path.join(DIR, 'major/wg-maj.fa')
    output:
        expand(
            os.path.join(DIR, 'major/indexes/wg-maj.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        os.path.join(DIR, 'major/indexes/wg-maj')
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params}'
 
rule build_grc_index:
    input:
        genome = GENOME
    output:
        expand(
            os.path.join(DIR, 'grc/wg.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        os.path.join(DIR, 'grc/wg')
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input.genome} {params}'

rule build_perA_index:
    input:
        os.path.join(DIR_PER, 'wg-perA.fa')
    output:
        os.path.join(DIR_PER, 'indexes/wg-perA.1.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perA.2.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perA.3.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perA.4.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perA.rev.1.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perA.rev.2.bt2')
    params:
        os.path.join(DIR_PER, 'indexes/wg-perA')
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params}'

rule build_perB_index:
    input:
        os.path.join(DIR_PER, 'wg-perB.fa')
    output:
        os.path.join(DIR_PER, 'indexes/wg-perB.1.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perB.2.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perB.3.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perB.4.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perB.rev.1.bt2'),
        os.path.join(DIR_PER, 'indexes/wg-perB.rev.2.bt2')
    params:
        os.path.join(DIR_PER, 'indexes/wg-perB')
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params}'

rule check_standard_genomes:
    input:
        expand(
            os.path.join(DIR, 'major/indexes/wg-maj.{idx}.bt2'),
            idx = IDX_ITEMS),
        expand(
            os.path.join(DIR, 'grc/wg.{idx}.bt2'),
            idx = IDX_ITEMS),
        expand(
            os.path.join(DIR_PER, 'indexes/wg-perA.{idx}.bt2'),
            idx = IDX_ITEMS, INDIV = INDIV),
        expand(
            os.path.join(DIR_PER, 'indexes/wg-perB.{idx}.bt2'),
            idx = IDX_ITEMS, INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'prepare_standard_genome.done')))
