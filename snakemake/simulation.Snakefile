'''
The snakemake pipeline for simulation is seperated 
to simplify the main pipeline.
'''

import os
import pandas as pd

configfile: 'simulation_local.yaml'

''' Load from config '''
CHROM = config['CHROM']
INDIV = config['INDIV']
EXP_LABEL = config['EXP_LABEL']
READS1 = config['READS1']
READS2 = config['READS2']
GROUP = config['GROUP']
ALN_MAPQ_THRSD = config['ALN_MAPQ_THRSD']
POP_THRSD = config['POP_THRSD']
POP_STOCHASTIC = config['POP_STOCHASTIC']
POP_BLOCK_SIZE = config['POP_BLOCK_SIZE']
POP_USE_LD = config['POP_USE_LD']

USE_PREBUILT = config['USE_PREBUILT']
SORT_SAM = config['SORT_SAM']

DIR = config['DIR']
GENOME = config['GENOME']
DIR_VCF = config['DIR_VCF']
VCF_PREFIX = config['VCF_PREFIX']
VCF_SUFFIX = config['VCF_SUFFIX']
CHR_PREFIX = config['CHR_PREFIX']
LENGTH_MAP = config['LENGTH_MAP']
CHROM_MAP = config['CHROM_MAP']

FAMILY = config['FAMILY']
SPOP = config['SPOP']
BCFTOOLS = config['BCFTOOLS']
SAMTOOLS = config['SAMTOOLS']
LIFTOVER = config['LIFTOVER']
PYTHON = config['PYTHON']
DIR_SCRIPTS = config['DIR_SCRIPTS']

THREADS = config['THREADS']
RAND_SEED = config['RAND_SEED']

#### Below args are only used in experiments
BEDTOOLS = config['BEDTOOLS']
MASON2 = config['MASON2']
NUM_SIM_READS = config['NUM_SIM_READS']
MAX_SYSTEM_THREADS = config['MAX_SYSTEM_THREADS']

''''''

IDX_ITEMS = ['1', '2', '3', '4', 'rev.1', 'rev.2']
DIR_PER = os.path.join(DIR, 'personalized/{INDIV}/')
PREFIX_PER = os.path.join(DIR_PER, EXP_LABEL)

include: 'experiments/prepare_per_and_grc.Snakefile'

rule all:
    input:
        os.path.join(DIR, 'simulation.done')

rule filter_vcf:
    input:
        vcf = os.path.join(DIR_VCF, VCF_PREFIX + '{CHROM}' + VCF_SUFFIX),
        chrom_map = CHROM_MAP
    output:
        vcf = temp(os.path.join(DIR, '{CHROM}_filtered.vcf'))
    shell:
        # Index VCF
        '{BCFTOOLS} index {input.vcf};'
        # Take PASS variants
        # Does not remove mnps, since they will be needed for constructing personalized reference genome, 
        # and will be removed when building major and refflow references.
        '{BCFTOOLS} view -r {wildcards.CHROM} -c 1 -f PASS {input.vcf} | {BCFTOOLS} annotate --rename-chrs {input.chrom_map} -o {output.vcf}'
        # '{BCFTOOLS} view -r {wildcards.CHROM} -c 1 -f PASS -V mnps,other {input.vcf} | {BCFTOOLS} annotate --rename-chrs {input.chrom_map} -o {output.vcf}'

rule aggregate_vcf:
    input:
        vcf = expand(os.path.join(DIR, '{CHROM}_filtered.vcf'), CHROM = CHROM)
    output:
        vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf')
    shell:
        '{BCFTOOLS} concat -o {output.vcf} {input.vcf}'

'''
Simulate reads from personalized genomes
'''
rule simulate_reads:
    input:
        hapA = PREFIX_PER + '-per_hapA.fa',
        hapB = PREFIX_PER + '-per_hapB.fa'
    output:
        readsA1 = temp(PREFIX_PER + '_hapA_1.fq'),
        readsA2 = temp(PREFIX_PER + '_hapA_2.fq'),
        readsB1 = temp(PREFIX_PER + '_hapB_1.fq'),
        readsB2 = temp(PREFIX_PER + '_hapB_2.fq'),
        samA = PREFIX_PER + '_hapA.sam',
        samB = PREFIX_PER + '_hapB.sam'
    params:
        num = NUM_SIM_READS,
        prefix = PREFIX_PER,
    #: set to MAX_SYSTEM_THREADS to avoid errors due to shared temp files
    threads: MAX_SYSTEM_THREADS
    shell:
        '{MASON2} --num-threads {threads} -ir {input.hapA} -n {params.num} '
        '-o {params.prefix}_hapA_1.fq -or {params.prefix}_hapA_2.fq '
        '-oa {params.prefix}_hapA.sam --read-name-prefix "{params.prefix}_hapA_simulated.";'
        '{MASON2} --num-threads {threads} -ir {input.hapB} -n {params.num} '
        '-o {params.prefix}_hapB_1.fq -or {params.prefix}_hapB_2.fq '
        '-oa {params.prefix}_hapB.sam --read-name-prefix "{params.prefix}_hapB_simulated.";'

rule merge_simulated_reads:
    input:
        readsA1 = PREFIX_PER + '_hapA_1.fq',
        readsA2 = PREFIX_PER + '_hapA_2.fq',
        readsB1 = PREFIX_PER + '_hapB_1.fq',
        readsB2 = PREFIX_PER + '_hapB_2.fq',
    output:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq'
    shell:
        'cat {input.readsA1} > {output.reads1};'
        'cat {input.readsA2} > {output.reads2};'
        'cat {input.readsB1} >> {output.reads1};'
        'cat {input.readsB2} >> {output.reads2};'

rule merge_simulated_sam:
    input:
        samA = PREFIX_PER + '_hapA.sam',
        samB = PREFIX_PER + '_hapB.sam'
    output:
        sam = temp(PREFIX_PER + '.sam'),
        sam1 = PREFIX_PER + '_1.sam',
        sam2 = PREFIX_PER + '_2.sam'
    shell:
        'grep ^@ {input.samA} > {output.sam};'
        #: get headers from B
        'grep ^@ {input.samB} | tail -n +2 >> {output.sam};'
        'grep -v ^@ {input.samA} >> {output.sam};'
        'grep -v ^@ {input.samB} >> {output.sam};'
        'samtools view -h -f 64 {output.sam} -o {output.sam1};'
        'samtools view -h -f 128 {output.sam} -o {output.sam2}'

rule check_simulation:
    input:
        expand(
            PREFIX_PER + '_{seg}.{type}',
            INDIV = INDIV,
            seg = ['1', '2'],
            type = ['fq', 'sam']
        )
    output:
        touch(temp(os.path.join(DIR, 'simulation.done')))


