dir=/home-1/cnaechy1@jhu.edu/scratch/relaxing
fq=$dir/syn_reads/na12878-chr9-phase3_hapA-1M_1.fq
fq2=$dir/syn_reads/na12878-chr9-phase3_hapA-1M_2.fq
g_sam=$dir/syn_reads/na12878-chr9-phase3_hapA-1M.fq.sam
bt_sam=$dir/alignments/na12878-chr9-phase3_hapA-1M-bt2_default.sam
ht_sam=$dir/alignments/na12878-chr9-phase3_hapA-1M-ht2_default.sam
bt_ref=$ARRAY/indexes/hs37d5.fa
ht_ref=$ARRAY/indexes/hisat2/grch37_snp/genome_snp

ht_gref=/scratch/groups/blangme2/jacob/vis_exp_chr9_NA19238/hisat_indexes/popcov10
ht_g_sam=$dir/alignments/na12878-chr9-phase3_popcov10_hapA-1M-ht2_default.sam

bt_per_ref=$dir/na12878/indexes/na12878-chr9-phase3_hapA
bt_per_sam=$dir/alignments/na12878-chr9-phase3_hapA-1M-bt2-per.sam
ht_g_per_ref=/scratch/groups/blangme2/jacob/vis_exp_chr9_NA19238/hisat_indexes/NA12878

# Generate target haplotype
# cd $dir/na12878
# python /home-1/cnaechy1@jhu.edu/scratch/FORGe-experiments/scripts/update_genome.py \
#     --ref $ARRAY/indexes/hs37d5.fa \
#     --vcf $dir/na12878/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf \
#     --chrom 9 \
#     --out-prefix na12878-chr9-phase3 \
#     --name NA12878

# Simulate reads
# /scratch/groups/blangme2/naechyun/software/mason-0.1.2-Linux-x86_64/bin/mason \
#     illumina $dir/na12878/na12878-chr9-phase3_hapA.fa \
#     -N 1000000 -sq -n 100 -hs 0 -hi 0 -hn 1 -mp \
#     -o na12878-chr9-phase3_hapA-1M.fq

# Personalized Genome
#time bowtie2 -p 4 -x $bt_per_ref -U $fq -S $bt_per_sam 
python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $bt_per_sam

#time bowtie2 -p 4 -x $bt_ref -U $fq -S $bt_sam 
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $bt_sam
#time hisat2 -p 4 -x $ht_ref -U $fq -S $ht_sam
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $ht_sam

# Graph
#time hisat2 -p 4 -x $ht_g_per_ref -U $fq -S $ht_g_sam
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $ht_g_sam
#time hisat2 -p 4 -x $ht_gref -U $fq -S $ht_g_sam
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $ht_g_sam
