rule get_vcf_stats_for_each_individual:
    input:
        os.path.join(DIR_PER, '21-per.vcf')
    output:
        os.path.join(DIR_PER, '21-per.stats')
    shell:
        '{BCFTOOLS} stats {input} > {output}'

rule summarize_indiv_variants:
    input:
        expand(os.path.join(DIR_PER, '21-per.stats'), INDIV = INDIV)
    output:
        os.path.join(DIR_RESULTS, 'variants/var_analysis.tsv')
    run:
        # shell("mkdir -p {DIR_RESULTS}/variants;")
        retrieve_indiv_variants(INDIV, input, output)

rule check_analysis:
    input:
        # expand(os.path.join(DIR_PER, '21-per.stats'), INDIV = INDIV)
        # os.path.join(DIR, 'analysis.tsv.done')
        os.path.join(DIR_RESULTS, 'variants/var_analysis.tsv')
    output:
        touch(temp(os.path.join(DIR, 'analysis.done')))
