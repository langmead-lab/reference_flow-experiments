SCRIPT="/scratch/groups/blangme2/naechyun/relaxing/scripts"
VCF="../21.vcf"
GENOME="../chr21.fa"
CHROM="21"
SAMPLE="NA12878"
NUM_SIM_READS=1000000

#: marcc modules
module load gcc/5.5.0
module load vcftools

set -x

#: keep SNPs and INDELs, remove others
~/bin/bcftools view --threads 8 -V mnps,other $VCF > ${CHROM}_remove_mnps_other.vcf
#: build var indexes
python ${SCRIPT}/update_genome.py \
    --ref $GENOME --vcf ${CHROM}_remove_mnps_other.vcf \
    --chrom $CHROM --out-prefix $SAMPLE --name $SAMPLE --include-indels

#: simulate reads
mkdir reads_mason2
cd reads_mason2
sh ${SCRIPT}/simulate_reads_mason2.sh ../${SAMPLE}_hapA.fa ../${SAMPLE}_hapB.fa $NUM_SIM_READS
cd ..

#: build major allele reference
~/bin/bcftools view -O v --threads 8 -q 0.5 $VCF -e 'AF = 0.5' -v snps,indels | \
    vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | \
    bgzip -@ 8 > ${CHROM}_maj.vcf.gz
~/bin/bcftools index ${CHROM}_maj.vcf.gz
~/bin/bcftools consensus -f $GENOME ${CHROM}_maj.vcf.gz > ${CHROM}_h37maj.fa
bgzip -cd ${CHROM}_maj.vcf.gz | python ${SCRIPT}/update_genome.py \
    --ref $GENOME \
    --chrom $CHROM --out-prefix ${CHROM}_h37maj \
    --include-indels --var-only 1

#: run test
sh $SCRIPT/test_prepare_genomes.sh $SAMPLE $CHROM
