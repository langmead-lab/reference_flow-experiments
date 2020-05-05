'''
Prepare genomes and indexes for experimental purposes.
1. Construct diploid personalized reference genomes and their indexes
2. Build GRC indexes
'''

rule build_grc_index:
    input:
        GENOME
    output:
        os.path.join(DIR, 'grc/indexes/' + EXP_LABEL + '.1.bt2'),
        os.path.join(DIR, 'grc/indexes/' + EXP_LABEL + '.2.bt2'),
        os.path.join(DIR, 'grc/indexes/' + EXP_LABEL + '.3.bt2'),
        os.path.join(DIR, 'grc/indexes/' + EXP_LABEL + '.4.bt2'),
        os.path.join(DIR, 'grc/indexes/' + EXP_LABEL + '.rev.1.bt2'),
        os.path.join(DIR, 'grc/indexes/' + EXP_LABEL + '.rev.2.bt2')
    params:
        os.path.join(DIR, 'grc/indexes/' + EXP_LABEL)
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params}'

rule build_per:
    input:
        genome = GENOME,
        vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf')
    output:
        hapA = os.path.join(DIR_PER, EXP_LABEL + '-per_hapA.fa'),
        hapB = os.path.join(DIR_PER, EXP_LABEL + '-per_hapB.fa'),
        var = os.path.join(DIR_PER, EXP_LABEL + '-per.var'),
        vcf = os.path.join(DIR_PER, EXP_LABEL + '-per.vcf')
    params:
        out_prefix = os.path.join(DIR_PER, EXP_LABEL + '-per')
    shell:
        '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
        '    --ref {input.genome} --vcf {input.vcf} --name {wildcards.INDIV}'
        '    --out-prefix {params.out_prefix} '
        '    --include-indels'

rule build_perA_index:
    input:
        os.path.join(DIR_PER, EXP_LABEL + '-per_hapA.fa')
    output:
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perA.1.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perA.2.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perA.3.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perA.4.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perA.rev.1.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perA.rev.2.bt2')
    params:
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perA')
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params}'

rule build_perB_index:
    input:
        os.path.join(DIR_PER, EXP_LABEL + '-per_hapB.fa')
    output:
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perB.1.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perB.2.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perB.3.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perB.4.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perB.rev.1.bt2'),
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perB.rev.2.bt2')
    params:
        os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-perB')
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params}'

rule check_prepare_per:
    input:
        expand(
            os.path.join(DIR_PER, 'indexes/' + EXP_LABEL + '-per{hap}.{IDX_ITEMS}.bt2'),
            IDX_ITEMS = IDX_ITEMS, INDIV = INDIV, hap = ['A', 'B']
        )
    output:
        touch(temp(os.path.join(DIR, 'prepare_per.done')))
