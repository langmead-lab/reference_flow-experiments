VCF="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/21.vcf"
CHROM="21"
DIR_SAMPLE="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/samples_in_each_pop_superpop"
CAT=$1
FRAC=$2
GENOME="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/chr21.fa"

usage(){
    echo "Usage $0 <mode> <frac> <stochastic> <block_size>"
    echo "    mode: operating mode [\"pop/superpop\"]"
    echo "    frac: frac := #haps/min_frequency, e.g. 2: major allele, 4: first quantile and above [INT]"
    echo "    stochastic: 1 to enable stochastic genome update; 0 to disable [INT]"
    echo "    block_size: size of stochastic update blocks; set anything for deterministic updates [INT]"
    exit
}

if [ "$#" -ne 4 ]
then
    echo "ERROR: incorrect number of arguments"
    usage
fi

if [ $CAT == "pop" ]
then
    array=( "ACB" "ASW" "BEB" "CDX" "CEU" "CHB" "CHS" "CLM" "ESN" "FIN" "GBR" "GIH" "GWD" "IBS" "ITU" "JPT" "KHV" "LWK" "MSL" "MXL" "PEL" "PJL" "PUR" "STU" "TSI" "YRI" )
elif [ $CAT == "superpop" ]
then
    #: for testing purpose
    #array=( "EUR" )
    array=( "EUR" "AMR" "SAS" "AFR" "EAS" )
else
    echo "ERROR: invalid input: $1"
    usage
fi

module load gcc/5.5.0
module load vcftools

# set -x 

for s in "${array[@]}"
do
    if [ -f "${CHROM}_${CAT}_${s}.vcf.gz" ]; then
        echo "Use existing ${CHROM}_${CAT}_${s}.vcf.gz"
    else
        ~/bin/bcftools view -S $DIR_SAMPLE/samples_${CAT}_${s}.txt --force-samples $VCF -V mnps,other |\
            vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | \
            bgzip -@ 8 > ${CHROM}_${CAT}_${s}.vcf.gz
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
        bgzip -@ 8 > ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz
        ~/bin/bcftools index ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz
    fi
    
    echo "Threshold = $THRSD"

    #: non-stochastic update
    if [ $3 == "0" ]; then
        # ~/bin/bcftools consensus -f $GENOME ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz > ${CHROM}_${CAT}_${s}_thrsd${FRAC}.fa
        bgzip -cd ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz |\
            python3.7 $REL/scripts/update_genome.py \
                --ref $GENOME --chrom $CHROM --out-prefix ${CHROM}_${CAT}_${s}_thrsd${FRAC} \
                --include-indels 1
    elif [ $3 == "1" ]; then
        bgzip -cd ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf.gz |\
            python3.7 $REL/scripts/update_genome.py \
                --ref $GENOME --chrom $CHROM --out-prefix ${CHROM}_${CAT}_${s}_thrsd${FRAC}_stochastic_b${4} \
                --include-indels 1 --stochastic -rs 0 --block-size $4
    fi
done
