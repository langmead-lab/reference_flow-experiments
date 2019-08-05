usage(){
    echo "Usage: $(basename $0) [-cnpS] -f fasta -s id -v vcf"
    echo "------ Requirements -----"
    echo "  -c  name of chromosome"
    echo "  -f  reference FASTA file"
    echo "  -s  target individual in 1000 Genomes Project"
    echo "  -v  1000 Genomes Project VCF file for target chromosome"
    echo "------ Options -----"
    echo "  -n  number of reads simulated [1000000]"
    echo "  -p  number of threads [8]"
    echo "  -S  path of scripts ['$'REL/scripts]"
    exit
}

SCRIPTS="/scratch/groups/blangme2/naechyun/relaxing/scripts"
THREADS="8"
NUM_SIM_READS=1000000
#VCF="../21.vcf"
#GENOME="../chr21.fa"
#CHROM="21"
#SAMPLE="NA12878"

while getopts n:p:S:c:f:s:v:h: option
do
    case "${option}"
    in
    #: options
    n) 
        NUM_SIM_READS=${OPTARG}
        echo "Set number of simulated reads -> $NUM_SIM_READS"
        ;;
    p) 
        THREADS=${OPTARG}
        echo "Set number of threads -> $THREADS"
        ;;
    S) 
        SCRIPTS=${OPTARG}
        echo "Set script directory -> $SCRIPTS"
        ;;

    #: requirements
    c) CHROM=${OPTARG};;
    f) GENOME=${OPTARG};;
    s) SAMPLE=${OPTARG};;
    v) VCF=${OPTARG};;

    #: help
    h)
        usage;;
    *)
        usage
    esac
done

if [ -z ${GENOME+x} ]
then
    echo "Error: required input GENOME (-f) is not set"
    usage
fi
if [ -z ${SAMPLE+x} ]
then
    echo "Error: required input SAMPLE (-s) is not set"
    usage
fi
if [ -z ${VCF+x} ]
then
    echo "Error: required input VCF (-v) is not set"
    usage
fi

#: exit when any command fails
set -e

#: marcc modules
module load gcc/5.5.0
module load vcftools

# set -x

#: keep SNPs and INDELs, remove others
~/bin/bcftools view --threads $THREADS -V mnps,other $VCF > ${CHROM}_remove_mnps_other.vcf
#: build var indexes
python ${SCRIPTS}/update_genome.py \
    --ref $GENOME --vcf ${CHROM}_remove_mnps_other.vcf \
    --chrom $CHROM --out-prefix $SAMPLE --name $SAMPLE --include-indels

#: simulate reads
mkdir -p reads_mason2
cd reads_mason2
sh ${SCRIPTS}/simulate_reads_mason2.sh ../${SAMPLE}_hapA.fa ../${SAMPLE}_hapB.fa $NUM_SIM_READS
cd ..

#: build major allele reference
~/bin/bcftools view -O v --threads $THREADS -q 0.5 $VCF -e 'AF = 0.5' -v snps,indels | \
    vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | \
    bgzip -@ $THREADS > ${CHROM}_maj.vcf.gz
~/bin/bcftools index ${CHROM}_maj.vcf.gz
~/bin/bcftools consensus -f $GENOME ${CHROM}_maj.vcf.gz > ${CHROM}_h37maj.fa
bgzip -cd ${CHROM}_maj.vcf.gz | python ${SCRIPTS}/update_genome.py \
    --ref $GENOME \
    --chrom $CHROM --out-prefix ${CHROM}_h37maj \
    --include-indels --var-only

#: run test
# sh $SCRIPTS/test_prepare_genomes.sh $SAMPLE $CHROM

#: build bowtie2 indexes
mkdir -p indexes
bowtie2-build --threads $THREADS ${CHROM}_h37maj.fa ${CHROM}_h37maj
