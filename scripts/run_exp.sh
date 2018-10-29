dir=/home-1/cnaechy1@jhu.edu/scratch/relaxing
fq=$dir/syn_reads/na12878-chr9-phase3_hapA-1M_1.fq
fq2=$dir/syn_reads/na12878-chr9-phase3_hapA-1M_2.fq
g_sam=$dir/syn_reads/na12878-chr9-phase3_hapA-1M.fq.sam
#bt_sam=$dir/alignments/na12878-chr9-phase3_hapA-1M-bt2.sam
#ht_sam=$dir/alignments/na12878-chr9-phase3_hapA-1M-ht2.sam
bt_sam=$dir/alignments/na12878-hapA-wg-bt2.sam
bt_sam_9=$dir/alignments/na12878-hapA-chr9-bt2.sam
bt_ref=$ARRAY/indexes/hs37d5.fa
bt_ref_9=$dir/chr9/chr9

ht_ref=$ARRAY/indexes/hisat2/grch37/genome
ht_ref_snp=$ARRAY/indexes/hisat2/grch37_snp/genome_snp
ht_ref_9=$dir/chr9/chr9
ht_sam=$dir/alignments/na12878-hapA-wg-ht2.sam
ht_sam_9=$dir/alignments/na12878-hapA-chr9-ht2.sam

ht_g_ref=/scratch/groups/blangme2/jacob/vis_exp_chr9_NA19238/hisat_indexes/popcov70
#ht_g_sam=$dir/alignments/na12878-chr9-phase3_popcov10_hapA-1M-ht2.sam
ht_g_sam=$dir/alignments/na12878-chr9-phase3_popcov70_hapA-1M-ht2.sam

bt_per_ref=$dir/na12878/indexes/na12878-chr9-phase3_hapA
bt_per_sam=$dir/alignments/na12878-hapA-chr9-per-bt2.sam
ht_per_ref=$dir/na12878/indexes/na12878-chr9-phase3_hapA
ht_per_sam=$dir/alignments/na12878-hapA-chr9-per-ht2.sam

bwa_per_sam=$dir/alignments/na12878-hapA-chr9-bwamem-per.sam
bwa_sam_9=$dir/alignments/na12878-hapA-chr9-bwamem.sam
bwa_sam=$dir/alignments/na12878-hapA-wg-bwamem.sam
bwa_per_ref=$dir/na12878/na12878-chr9-phase3_hapA.fa
bwa_ref=$ARRAY/indexes/hs37d5.fa
bwa_ref_9=$dir/chr9/chr9.fa

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
echo "bt2, personalized chr9"
time bowtie2 -p 4 -x $bt_per_ref -U $fq -S $bt_per_sam 
python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $bt_per_sam
#echo "ht2, C9_PG"
#time hisat2 -p 4 -x $ht_per_ref -U $fq -S $ht_per_sam --no-spliced-alignment
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $ht_per_sam

#echo "bt2, whole genome"
#time bowtie2 -p 4 -x $bt_ref -U $fq -S $bt_sam 
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $bt_sam
#echo "bt2, chr9"
#time bowtie2 -p 4 -x $bt_ref_9 -U $fq -S $bt_sam_9
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $bt_sam_9

#time bwa mem -t 4 $bwa_per_ref $fq > $bwa_per_sam
#time bwa mem -t 4 $bwa_ref $fq > $bwa_sam
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $bwa_sam
#time bwa mem -t 4 $bwa_ref_9 $fq > $bwa_sam_9
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $bwa_sam_9

#echo "ht2, whole genome (linear genome)"
#time hisat2 -p 4 -x $ht_ref -U $fq -S $ht_sam --no-spliced-alignment
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $ht_sam
#echo "ht2, chr 9"
#time hisat2 -p 4 -x $ht_ref_9 -U $fq -S $ht_sam_9 --no-spliced-alignment
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $ht_sam_9
#echo "ht2, whole genome with snps (provided by ht2)"
#time hisat2 -p 4 -x $ht_ref_snp -U $fq -S $ht_sam --no-spliced-alignment
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $ht_sam

# Graph
#echo "ht2, popcov chr9 genome: 70"
#time hisat2 -p 4 -x $ht_g_ref -U $fq -S $ht_g_sam --no-spliced-alignment
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $ht_g_sam

# Major allele
#ht_ref_9_major=$dir/chr9/major_allele/chr9_major
#ht_maj_sam=$dir/alignments/na12878-hapA-chr9_major-ht2.sam
#bt_ref_9_major=$dir/chr9/major_allele/chr9_major
#bt_maj_sam=$dir/alignments/na12878-hapA-chr9_major-bt2.sam
#bwa_ref_9_major=$dir/chr9/major_allele/chr9_major
#echo "ht2, chr9 major"
#time hisat2 -p 4 -x $ht_ref_9_major -U $fq -S $ht_maj_sam --no-spliced-alignment
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $ht_maj_sam
#echo "bt2, chr9 major"
#time bowtie2 -p 4 -x $bt_ref_9_major -U $fq -S $bt_maj_sam 
#python -O $dir/scripts/analyze_sam.py -s 1 -g $g_sam -n $bt_maj_sam
