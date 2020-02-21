usage(){
    #echo "Usage: $(basename $0) [-blrtpS] -c chrom -C cat -f fasta -v vcf"
    echo "Usage: $(basename $0) [-lr] [-b INT] [-t INT] [-p INT] [-S path] [-x STR] -c chrom -C cat -f fasta -v vcf"
    echo "------ Requirements -----"
    echo "  -c  name of chromosome"
    echo "  -C  category {superpop/pop}"
    echo "  -f  reference FASTA file"
    echo "  -v  1000 Genomes Project VCF file for target chromosome"
    echo "------ Options -----"
    echo "  -b  size of block for stochastic genomes [1]"
    echo "  -l  blocking with pseudo-LD, set to enable [False]"
    echo "  -r  stochastic population genomes, set to enable [False]"
    echo "  -t  population frequency cutoff level [1]"
    echo "  -p  number of threads [8]"
    echo "  -P  directory for 1KG population individuals []"
    echo "  -S  path of scripts []"
    echo "  -x  indivs to exclude, separate by comma [None]"
    exit
}

DIR_SAMPLE="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/samples_in_each_pop_superpop"
SCRIPTS="/scratch/groups/blangme2/naechyun/relaxing/scripts"
THREADS=8
FRAC=0
BLOCK_SIZE=1
STOCHASTIC=0
LD=0
EXCLUDE_ID=""

while getopts c:C:f:s:v:b:t:p:P:S:x:rlh option
do
    case "${option}"
    in
    #: options
    b) 
        BLOCK_SIZE=${OPTARG}
        echo "Set stochastic block size -> $BLOCK_SIZE"
        ;;
    l)
        LD=1
        echo "Use pseudo-LD blocking"
        ;;
    r)
        STOCHASTIC=1
        echo "Use stochastic populations genomes"
        ;;
    t)
        FRAC=${OPTARG}
        echo "Set frequency cutoff level -> $FRAC"
        ;;
    p) 
        THREADS=${OPTARG}
        echo "Set number of threads -> $THREADS"
        ;;
    S) 
        SCRIPTS=${OPTARG}
        echo "Set script directory -> $SCRIPTS"
        ;;
    P) 
        DIR_SAMPLE=${OPTARG}
        echo "Set 1KG sample directory -> $DIR_SAMPLE"
        ;;
    x) 
        EXCLUDE_ID=${OPTARG}
        echo "Set indivs to exclude -> $EXCLUDE_ID"
        ;;

    #: requirements
    c) CHROM=${OPTARG};;
    C) CAT=${OPTARG};;
    f) GENOME=${OPTARG};;
    v) VCF=${OPTARG};;

    #: help
    h)
        usage;;
    *)
        usage
    esac
done
shift "$(($OPTIND -1))"

#: check if required fields are given
if [ -z ${CHROM+x} ]
then
    echo "error: required input chrom (-c) is not set"
    usage
fi
if [ -z ${CAT+x} ]
then
    echo "error: required input cat (-C) is not set"
    usage
fi
if [ -z ${GENOME+x} ]
then
    echo "error: required input genome (-f) is not set"
    usage
fi
if [ -z ${VCF+x} ]
then
    echo "Error: required input VCF (-v) is not set"
    usage
fi

#: set population
if [ $CAT == "pop" ]
then
    array=( "ACB" "ASW" "BEB" "CDX" "CEU" "CHB" "CHS" "CLM" "ESN" "FIN" "GBR" "GIH" "GWD" "IBS" "ITU" "JPT" "KHV" "LWK" "MSL" "MXL" "PEL" "PJL" "PUR" "STU" "TSI" "YRI" )
elif [ $CAT == "superpop" ]
then
    #: for testing purpose
    #array=( "EUR" )
    array=( "EUR" "AMR" "SAS" "AFR" "EAS" )
else
    echo "ERROR: invalid input: $CAT"
    usage
fi

#: marcc
module load gcc/5.5.0
module load vcftools

#: exit when any command fails
set -e
# set -x 

for s in "${array[@]}"
do
    if [ -f "${CHROM}_${CAT}_${s}.vcf.gz" ]; then
        echo "Use existing ${CHROM}_${CAT}_${s}.vcf.gz"
    else
        ~/bin/bcftools view -S $DIR_SAMPLE/samples_${CAT}_${s}.txt --force-samples $VCF -V mnps,other |\
            vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | \
            bgzip -@ $THREADS > ${CHROM}_${CAT}_${s}.vcf.gz
    fi
    if [ -f "vcf_${CAT}_${s}_samples.txt" ]; then
        echo "Use existing vcf_${CAT}_${s}_samples.txt"
    else
        ~/bin/bcftools view -h ${CHROM}_${CAT}_${s}.vcf.gz | tail -1 > vcf_${CAT}_${s}_samples.txt
    fi

    THRSD=0
    if [ -f "${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz" ]; then
        echo "Use existing ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz"
    else
        if [ "$FRAC" -ne "0" ]; then
            THRSD=$(( ($(cat vcf_${CAT}_${s}_samples.txt | wc -w) - 9) * 2 / $FRAC ))
        fi
        FILTER="AC > $THRSD"
        ~/bin/bcftools view -i "$FILTER" -v snps,indels ${CHROM}_${CAT}_${s}.vcf.gz |\
        bgzip -@ $THREADS > ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz
        ~/bin/bcftools index ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz
    fi
    
    echo "Threshold = $THRSD"

    #: non-stochastic update
    if [ $STOCHASTIC == "0" ]; then
        bgzip -@ $THREADS -cd ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz |\
        python3.7 $SCRIPTS/update_genome.py \
            --ref $GENOME --chrom $CHROM --out-prefix ${CHROM}_${CAT}_${s}_thrsd${FRAC} \
            --include-indels
    elif [ $STOCHASTIC == "1" ]; then
        if [ $LD == "0" ]; then
            bgzip -@ $THREADS -cd ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz |\
            python3.7 $SCRIPTS/update_genome.py \
                --ref $GENOME --chrom $CHROM --out-prefix ${CHROM}_${CAT}_${s}_thrsd${FRAC}_stochastic_b${BLOCK_SIZE} \
                --include-indels --stochastic -rs 0 --block-size $BLOCK_SIZE
        else
            bgzip -@ $THREADS -cd ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz |\
            python3.7 $SCRIPTS/update_genome.py \
                --ref $GENOME --chrom $CHROM --out-prefix ${CHROM}_${CAT}_${s}_thrsd${FRAC}_stochastic_b${BLOCK_SIZE} \
                --include-indels --stochastic -rs 0 --block-size $BLOCK_SIZE --ld --exclude-name $EXCLUDE_ID
        fi
    fi
done
