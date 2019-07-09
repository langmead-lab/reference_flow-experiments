SCRIPT="/scratch/groups/blangme2/naechyun/relaxing/scripts"
VCF="../21.vcf"
GENOME="../chr21.fa"
CHROM="21"
SAMPLE="NA12878"
NUM_SIM_READS=1000000

set -x

#: marcc modules
module load gcc/5.5.0
module load vcftools

#: keep SNPs and INDELs, remove others
~/bin/bcftools-1.9 view -V mnps,other $VCF > ${CHROM}_remove_mnps_other.vcf
#: TODO
#: see if using bcftools for building consensus genome
#: build personalized genome
python ${SCRIPT}/update_genome.py \
    --ref $GENOME --vcf ${CHROM}_remove_mnps_other.vcf \
    --chrom $CHROM --out-prefix $SAMPLE --name $SAMPLE --include-indels 1

#: simulate reads
mkdir reads_mason2
cd reads_mason2
sh ${SCRIPT}/simulate_reads_mason2.sh ../${SAMPLE}_hapA.fa ../${SAMPLE}_hapB.fa $NUM_SIM_READS
cd ..

#: build major allele reference
~/bin/bcftools-1.9 view -O v --threads 8 -q 0.5 $VCF -e 'AF = 0.5' -v snps,indels | \
    vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | \
    bgzip -@ 8 > ${CHROM}_maj_with_indels.vcf.gz
~/bin/bcftools-1.9 index ${CHROM}_maj_with_indels.vcf.gz
~/bin/bcftools-1.9 consensus -f $GENOME ${CHROM}_maj_with_indels.vcf.gz > ${CHROM}_h37maj_with_indels.fa
bgzip -cd ${CHROM}_maj_with_indels.vcf.gz | python ${SCRIPT}/update_genome.py \
    --ref $GENOME \
    --chrom $CHROM --out-prefix ${CHROM}_h37maj_with_indels \
    --include-indels 1 --var-only 1

#: run test
sh $SCRIPT/test_prepare_genomes.sh $SAMPLE $CHROM
