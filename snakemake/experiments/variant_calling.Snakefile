""" This Snakefile performs variant calling using GATK.

The pipeline we used here includes (adding read-group),
MarkDuplicates, BQSR and HaplotypeCaller. VQSR is not included.

"""
# picard AddOrReplaceReadGroups I= $SORTED_BAM  O= $SORTED_RG_BAM   RGID=HiSeq2000 RGLB= $METHOD RGPL=illumina RGSM= $SAMPLE RGPU=hg38 TMP_DIR=tmp
# samtools index $SORTED_RG_BAM
# ~/miniconda3/bin/gatk --java-options "-XX:ParallelGCThreads=32" HaplotypeCaller -R $GENOME -I $SORTED_RG_BAM -O $VCF_GZ --tmp-dir tmp --native-pair-hmm-threads 8

PICARD = 'picard'
GATK = '~/miniconda3/bin/gatk' 
'''
GATK requires reads labelled with read groups and then indexed.
Here we use `picard AddOrReplaceReadGroups` to add read groups and `samtools index` to index the new results.
'''
RGID = INDIV
RGPL = 'illumina'
RGPU = 'HFLT3DSXX' # for NYGC # 'readgroup.pu'

# NYGC variant calling pipeline
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20190405_NYGC_b38_pipeline_description.pdf

VARCALL_OBJECTS = [
    EXP_LABEL + '-GRC-sorted',
    EXP_LABEL + '-major-liftover-sorted',
    POP_DIRNAME + '/' + EXP_LABEL + '-refflow-{}-{}-liftover-sorted'.format(ALN_MAPQ_THRSD, POP_DIRNAME),
    EXP_LABEL + '-per-merged-liftover-sorted',
#     EXP_LABEL + '-vg_0.1-pe'
]

rule add_rg:
    input:
        os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '.bam')
    output:
        os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG.bam')
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp/' + INDIV),
        label = INDIV
    shell:
        '{PICARD} AddOrReplaceReadGroups I={input}  O={output}  RGID={RGID} RGLB={params.label} RGPL={RGPL} RGSM={wildcards.INDIV} RGPU={RGPU} TMP_DIR={params.tmp}'

rule mark_dup:
    input:
        os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG.bam')
    output:
        bam = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG-dedup.bam'),
        metric = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '.dedup.metrics')
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp')
    threads: 112
    shell:
        '{PICARD} MarkDuplicates INPUT= {input} OUTPUT= {output.bam} METRICS_FILE= {output.metric} TMP_DIR={params.tmp}'
    
rule build_bqsr_table:
    input:
        bam = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG-dedup.bam'),
        genome = GENOME,
        known_vcf = DBSNP_COMMON
    output:
        table = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '.bqsr.table'),
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp')
    threads: THREADS
    shell:
        '{GATK} BaseRecalibrator -I {input.bam} -R {input.genome} --known-sites {input.known_vcf} -O {output.table} --java-options "-XX:ConcGCThreads={threads}"'

rule apply_bqsr:
    input:
        bam = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG-dedup.bam'),
        genome = GENOME,
        table = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '.bqsr.table'),
    output:
        bam = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG-dedup-bqsr.bam'),
    threads: THREADS
    shell:
        '{GATK} ApplyBQSR -R {input.genome} -I {input.bam} --bqsr-recal-file {input.table} -O {output.bam} --java-options "-XX:ConcGCThreads={threads}"'

# rule index_rg:
#     input:
#         os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG.bam')
#     output:
#         os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG.bam.bai')
#     shell:
#         '{SAMTOOLS} index {input}'

rule index_post_bqsr:
    input:
        os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG-dedup-bqsr.bam'),
    output:
        os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG-dedup-bqsr.bam.bai'),
    shell:
        '{SAMTOOLS} index {input}'

'''
Variant calling using GATK
'''
GENOME_DICT = '.'.join(GENOME.split('.')[:-1]) + '.dict'
rule build_genome_dict:
    '''
    Build `.dict` file for reference genome
    '''
    input:
        GENOME
    output:
        GENOME_DICT
    shell:
        '{PICARD} CreateSequenceDictionary R={input} O={output}'

rule GATK_haplotypecaller:
    input:
        bam = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG-dedup-bqsr.bam'),
        bai = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG-dedup-bqsr.bam.bai'),
        genome_dict = GENOME_DICT
    output:
        vcf_gz = os.path.join(DIR_FIRST_PASS, '{VARCALL_OBJECTS}' + '-RG-dedup-bqsr-hapcal.vcf.gz')
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp')
    threads: THREADS
    shell:
        '{GATK} --java-options "-XX:ParallelGCThreads={threads}" HaplotypeCaller -R {GENOME} \
            -I {input.bam} -O {output.vcf_gz} --tmp-dir {params.tmp} --native-pair-hmm-threads {threads}'

rule check_variant_calling:
    input:
        expand(
            os.path.join(
                DIR_FIRST_PASS,
                '{VARCALL_OBJECTS}' + '-RG-dedup-bqsr-hapcal.vcf.gz'),
                # '{VARCALL_OBJECTS}' + '-RG-dedup.bam'),
                # '{VARCALL_OBJECTS}' + '-RG-dedup-bqsr.bam.bai'),
                # '{VARCALL_OBJECTS}' + '-RG.bam.bai'),
            INDIV = INDIV, VARCALL_OBJECTS = VARCALL_OBJECTS)
    output:
        touch(temp(os.path.join(DIR, 'var_calling.done')))
