usage(){
    echo "Usage: $(basename $0) [-brtpS] -f fasta -c chrome -v vcf -C cat"
    echo "------ Requirements -----"
    echo "  -c  name of chromosome"
    echo "  -C  category {superpop/pop} [superpop]"
    echo "  -f  reference FASTA file"
    echo "  -s  target individual in 1000 Genomes Project"
    echo "  -v  1000 Genomes Project VCF file for target chromosome"
    echo "------ Options -----"
    echo "  -b  size of block for stochastic genomes [1]"
    echo "  -r  stochastic population genomes, set to enable [False]"
    echo "  -t  population frequency cutoff level [1]"
    echo "  -p  number of threads [8]"
    echo "  -P  directory for 1KG population individuals"
    echo "  -S  path of scripts ['$'REL/scripts]"
    exit
}

#VCF="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/21.vcf"
#CHROM="21"
DIR_SAMPLE="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/samples_in_each_pop_superpop"
SCRIPTS="$REL/scripts"
THREADS=8
#CAT=$1
FRAC=0
BLOCK_SIZE=1
STOCHASTIC=0
#GENOME="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/chr21.fa"

while getopts c:C:f:s:v:b:r:t:p:S:h: option
do
    case "${option}"
    in
    #: options
    b) 
        BLOCK_SIZE=${OPTARG}
        echo "Set stochastic block size -> $BLOCK_SIZE"
        ;;
    r)
        echo "Use stochastic populations genomes"
        STOCHASTIC=1
        ;;
    t)
        FRAC=${OPTARG}
        echo "Set frequency cutoff level -> $FRAC"
        ;;
    p) 
        echo "Set number of threads -> $THREADS"
        THREADS=${OPTARG}
        ;;
    S) 
        echo "Set script directory -> $SCRIPTS"
        SCRIPTS=${OPTARG}
        ;;
    P) 
        echo "Set 1KG sample directory -> $DIR_SAMPLE"
        DIR_SAMPLE=${OPTARG}
        ;;

    #: requirements
    c) CHROM=${OPTARG};;
    C) CAT=${OPTARG};;
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
# usage(){
#     echo "Usage $0 <mode> <frac> <stochastic> <block_size>"
#     echo "    mode: operating mode [\"pop/superpop\"]"
#     echo "    frac: frac := #haps/min_frequency, e.g. 2: major allele, 4: first quantile and above [INT]"
#     echo "    stochastic: 1 to enable stochastic genome update; 0 to disable [INT]"
#     echo "    block_size: size of stochastic update blocks; set anything for deterministic updates [INT]"
#     exit
# }
# 
# if [ "$#" -ne 4 ]
# then
#     echo "ERROR: incorrect number of arguments"
#     usage
# fi

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
        # if [ "$FRAC" == "0" ]; then
        #     THRSD=0
        # else
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
        # ~/bin/bcftools consensus -f $GENOME ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz > ${CHROM}_${CAT}_${s}_thrsd${FRAC}.fa
        bgzip -@ $THREADS -cd ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz |\
            python3.7 $REL/scripts/update_genome.py \
                --ref $GENOME --chrom $CHROM --out-prefix ${CHROM}_${CAT}_${s}_thrsd${FRAC} \
                --include-indels
    elif [ $STOCHASTIC == "1" ]; then
        bgzip -@ $THREADS -cd ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz |\
            python3.7 $REL/scripts/update_genome.py \
                --ref $GENOME --chrom $CHROM --out-prefix ${CHROM}_${CAT}_${s}_thrsd${FRAC}_stochastic_b${BLOCK_SIZE} \
                --include-indels --stochastic -rs 0 --block-size $BLOCK_SIZE
    fi
done
