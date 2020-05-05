# picard AddOrReplaceReadGroups I= $SORTED_BAM  O= $SORTED_RG_BAM   RGID=HiSeq2000 RGLB= $METHOD RGPL=illumina RGSM= $SAMPLE RGPU=hg38 TMP_DIR=tmp
# samtools index $SORTED_RG_BAM
# ~/miniconda3/bin/gatk --java-options "-XX:ParallelGCThreads=32" HaplotypeCaller -R $GENOME -I $SORTED_RG_BAM -O $VCF_GZ --tmp-dir tmp --native-pair-hmm-threads 8

PICARD = 'picard'
GATK = '~/miniconda3/bin/gatk' 
'''
GATK requires reads labelled with read groups and then indexed.
Here we use `picard AddOrReplaceReadGroups` to add read groups and `samtools index` to index the new results.
'''
RGID = 'readgroup.id'
RGPL = 'illumina.HiSeq2000'
RGPU = 'readgroup.pu'

rule add_rg_grc:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-sorted.bam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-sorted-RG.bam')
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp'),
        method = 'GRC'
    shell:
        '{PICARD} AddOrReplaceReadGroups I={input}  O={output}  RGID={RGID} RGLB={params.method} RGPL={RGPL} RGSM={wildcards.INDIV} RGPU={RGPU} TMP_DIR={params.tmp}'

rule index_rg_grc:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-sorted-RG.bam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-sorted-RG.bam.bai')
    shell:
        '{SAMTOOLS} index {input}'

rule add_rg_major:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted.bam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted-RG.bam')
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp'),
        method = 'Major'
    shell:
        '{PICARD} AddOrReplaceReadGroups I={input}  O={output}  RGID={RGID} RGLB={params.method} RGPL={RGPL} RGSM={wildcards.INDIV} RGPU={RGPU} TMP_DIR={params.tmp}'

rule index_rg_major:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted-RG.bam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted-RG.bam.bai')
    shell:
        '{SAMTOOLS} index {input}'

rule add_rg_refflow:
    input:
        os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted-RG.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    params:
        tmp = os.path.join(DIR_SECOND_PASS, 'tmp'),
        method = 'refflow-{}-{}'.format(ALN_MAPQ_THRSD, POP_DIRNAME)
    shell:
        '{PICARD} AddOrReplaceReadGroups I={input}  O={output}  RGID={RGID} RGLB={params.method} RGPL={RGPL} RGSM={wildcards.INDIV} RGPU={RGPU} TMP_DIR={params.tmp}'

rule index_rg_refflow:
    input:
        os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted-RG.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted-RG.bam.bai'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    shell:
        '{SAMTOOLS} index {input}'

rule add_rg_per:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted.bam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted-RG.bam')
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp'),
        method = 'personalized'
    shell:
        '{PICARD} AddOrReplaceReadGroups I={input}  O={output}  RGID={RGID} RGLB={params.method} RGPL={RGPL} RGSM={wildcards.INDIV} RGPU={RGPU} TMP_DIR={params.tmp}'

rule index_rg_per:
    input:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted-RG.bam')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted-RG.bam.bai')
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

rule GATK_call_grc:
    input:
        bam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-sorted-RG.bam'),
        bai = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC-sorted-RG.bam.bai'),
        genome_dict = GENOME_DICT
    output:
        vcf_gz = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC.vcf.gz')
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp'),
        thread_hmm = 4
    threads: THREADS
    shell:
        '{GATK} --java-options "-XX:ParallelGCThreads={threads}" HaplotypeCaller -R {GENOME} \
            -I {input.bam} -O {output.vcf_gz} --tmp-dir {params.tmp} --native-pair-hmm-threads {params.thread_hmm}'

rule GATK_call_major:
    input:
        bam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted-RG.bam'),
        bai = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-liftover-sorted-RG.bam.bai'),
        genome_dict = GENOME_DICT
    output:
        vcf_gz = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major.vcf.gz')
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp'),
        thread_hmm = 4
    threads: THREADS
    shell:
        '{GATK} --java-options "-XX:ParallelGCThreads={threads}" HaplotypeCaller -R {GENOME} \
            -I {input.bam} -O {output.vcf_gz} --tmp-dir {params.tmp} --native-pair-hmm-threads {params.thread_hmm}'

rule GATK_call_refflow:
    input:
        bam = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted-RG.bam.bai'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        bai = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted-RG.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        genome_dict = GENOME_DICT
    output:
        vcf_gz = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}.vcf.gz'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    params:
        tmp = os.path.join(DIR_SECOND_PASS, 'tmp'),
        thread_hmm = 4
    threads: THREADS
    shell:
        '{GATK} --java-options "-XX:ParallelGCThreads={threads}" HaplotypeCaller -R {GENOME} \
            -I {input.bam} -O {output.vcf_gz} --tmp-dir {params.tmp} --native-pair-hmm-threads {params.thread_hmm}'

rule GATK_call_per:
    input:
        bam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted-RG.bam.bai'),
        bai = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per-merged-liftover-sorted-RG.bam.bai'),
        genome_dict = GENOME_DICT
    output:
        vcf_gz = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per.vcf.gz')
    params:
        tmp = os.path.join(DIR_FIRST_PASS, 'tmp'),
        thread_hmm = 4
    threads: THREADS
    shell:
        '{GATK} --java-options "-XX:ParallelGCThreads={threads}" HaplotypeCaller -R {GENOME} \
            -I {input.bam} -O {output.vcf_gz} --tmp-dir {params.tmp} --native-pair-hmm-threads {params.thread_hmm}'

rule check_variant_calling:
    input:
        expand(
            os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-GRC.vcf.gz'),
            INDIV = INDIV),
#         expand(
#             os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major.vcf.gz'),
#             INDIV = INDIV),
#         expand(
#             os.path.join(DIR_SECOND_PASS,
#             EXP_LABEL + '-refflow-{}-{}.vcf.gz'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
#             INDIV = INDIV),
#         expand(
#             os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-per.vcf.gz'),
#             INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'var_calling.done')))
