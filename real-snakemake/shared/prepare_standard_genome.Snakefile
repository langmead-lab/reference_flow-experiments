rule merge_major:
    input:
        expand(PREFIX_MAJOR + '.fa', CHROM = CHROM)
    output:
        os.path.join(DIR, 'major/wg_h37maj.fa')
    shell:
        'cat {input} >> {output}'

rule build_major_index:
    input:
        os.path.join(DIR, 'major/wg_h37maj.fa')
    output:
        expand(
            os.path.join(DIR, 'major/indexes/wg_h37maj.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        os.path.join(DIR, 'major/wg_h37maj/indexes')
    shell:
        'bowtie2-build --threads {THREADS} {input} {params}'
 
rule build_grc_index:
    input:
        genome = GENOME
    output:
        expand(
            os.path.join(DIR, 'grch37/wg.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        os.path.join(DIR, 'grch37/wg')
    shell:
        'bowtie2-build --threads {THREADS} {input.genome} {params}'

rule check_standard_genomes:
    input:
        expand(
            os.path.join(DIR, 'major/indexes/wg_h37maj.{idx}.bt2'),
            idx = IDX_ITEMS),
        expand(
            os.path.join(DIR, 'grch37/wg.{idx}.bt2'),
            idx = IDX_ITEMS)
    output:
        touch(temp(os.path.join(DIR, 'prepare_standard_genome.done')))
