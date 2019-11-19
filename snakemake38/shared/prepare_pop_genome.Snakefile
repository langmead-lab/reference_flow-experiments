'''
Rules for building population genomes
'''
rule prepare_pop_indiv:
    output:
        expand(
            os.path.join(DIR, '1KG_indivs/sample_superpop_{GROUP}.txt'),
            GROUP = GROUP
        )
    params:
        prefix = os.path.join(DIR, '1KG_indivs/sample')
    shell:
        '{PYTHON} {DIR_SCRIPTS}/list_indiv_from_pop.py '
        '-p {FAMILY} -sp {SPOP} -op {params.prefix}'

rule build_pop_vcf:
    input:
        vcf = UNPHASED_VCF_F,
        indiv_group = os.path.join(
            DIR,
            '1KG_indivs/sample_superpop_{GROUP}.txt'
        )
    output:
        vcf_gz = os.path.join(DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}.vcf.gz')
    shell:
        '{BCFTOOLS} view -S {input.indiv_group} '
        '--force-samples {input.vcf} -V mnps,other -m2 -M2 | '
        'bgzip > {output.vcf_gz}'

rule get_pop_sample:
    input:
        vcf_gz = os.path.join(DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}.vcf.gz')
    output:
        vcf_header = os.path.join(DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}.samples')
    shell:
        '{BCFTOOLS} view -h {input.vcf_gz} | tail -1 '
        '> {output.vcf_header}'

rule filter_pop_vcf:
    input:
        vcf_gz = os.path.join(DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}.vcf.gz'),
        vcf_header = os.path.join(DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}.samples')
    output:
        vcf = os.path.join(
            DIR_POP_GENOME,
            CHROM + '_superpop_{GROUP}_t' + str(POP_THRSD) + '.vcf'
        )
    run:
        fn = list({input.vcf_header})[0]
        with open(fn, 'r') as f:
            for line in f:
                n = len(line.split()) - 9
                thrsd = int(n * 2 * float(POP_THRSD))
                filt = 'AC > {}'.format(thrsd)
                break
        shell('{BCFTOOLS} view -i "{filt}" \
            -v snps,indels {input.vcf_gz} > {output.vcf};')

# rule build_pop_genome:
#     input:
#         vcf = os.path.join(
#             DIR_POP_GENOME,
#             CHROM + '_superpop_{GROUP}_t' + str(POP_THRSD) + '.vcf'
#         )
#     output:
#         os.path.join(
#             DIR_POP_GENOME_BLOCK,
#             CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.fa'
#         ),
#         os.path.join(
#             DIR_POP_GENOME_BLOCK,
#             CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.var'
#         ),
#         os.path.join(
#             DIR_POP_GENOME_BLOCK,
#             CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.vcf'
#         )
#     params:
#         prefix = os.path.join(
#             DIR_POP_GENOME_BLOCK,
#             CHROM + '_superpop_{GROUP}_' + POP_DIRNAME
#         )
#     run:
#         if POP_STOCHASTIC == 1 and POP_USE_LD == 1:
#             shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
#                 --ref {GENOME} --chrom {CHROM} --vcf {input.vcf} \
#                 --out-prefix {params.prefix} \
#                 --include-indels --stochastic -rs {RAND_SEED} \
#                 --block-size {POP_BLOCK_SIZE} --ld')
#         elif POP_STOCHASTIC == 1:
#             shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
#                 --ref {GENOME} --chrom {CHROM} --vcf {input.vcf} \
#                 --out-prefix {params.prefix} \
#                 --include-indels --stochastic -rs {RAND_SEED} \
#                 --block-size {POP_BLOCK_SIZE}')
#         else:
#             shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
#                 --ref {GENOME} --chrom {CHROM} --vcf {input.vcf} \
#                 --out-prefix {params.prefix} \
#                 --include-indels')

rule build_pop_genome_index:
    input:
        genome = DIR_POP_GENOME_BLOCK + POP_GENOME_SUFFIX + '.fa'
    output:
        DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.1.bt2',
        DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.2.bt2',
        DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.3.bt2',
        DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.4.bt2',
        DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.1.bt2',
        DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.2.bt2'
    params:
        prefix = DIR_POP_GENOME_BLOCK_IDX +
        POP_GENOME_SUFFIX
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input.genome} {params.prefix};'

# rule check_pop_genome:
#     input:
#         expand(
#             DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.{IDX_ITEMS}.bt2',
#             GROUP = GROUP, IDX_ITEMS = IDX_ITEMS
#         )
#     output:
#         touch(temp(os.path.join(DIR, 'prepare_pop_genome.done')))


''' Build pop genomes using gnomad data '''
rule build_pop_genome_gnomad:
    input:
        vcf = UNPHASED_VCF_F
    output:
        os.path.join(
            DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX + '.fa'
            # 'chr' + CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.fa'
        ),
        os.path.join(
            DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX + '.var'
            # 'chr' + CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.var'
        ),
        os.path.join(
            DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX + '.vcf'
            # 'chr' + CHROM + '_superpop_{GROUP}_' + POP_DIRNAME + '.vcf'
        )
    params:
        prefix = os.path.join(
            DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX
            # 'chr' + CHROM + '_superpop_{GROUP}_' + POP_DIRNAME
        )
    run:
        if POP_STOCHASTIC == 1:
            pop_count = DICT_GNOMAD_POP_SIZE[wildcards.GROUP]
            shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
                --ref {GENOME} --chrom {CHROM} --vcf {input.vcf} \
                --out-prefix {params.prefix} \
                --include-indels --stochastic -rs {RAND_SEED} \
                --block-size {POP_BLOCK_SIZE} \
                -d gnomad --gnomad-ac-field AC_{wildcards.GROUP} \
                --gnomad-pop-count {pop_count}')
        # else:
        #     shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
        #         --ref {GENOME} --chrom {CHROM} --vcf {input.vcf} \
        #         --out-prefix {params.prefix} \
        #         --include-indels')

rule check_pop_genome:
    input:
        expand(
            DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.{IDX_ITEMS}.bt2',
            GROUP = GROUP, IDX_ITEMS = IDX_ITEMS
        )
    output:
        touch(temp(os.path.join(DIR, 'prepare_pop_genome.done')))
