VCF="../21.vcf"
CHROM="21"
DIR_SAMPLE="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/samples_in_each_pop_superpop"
CAT=$1
FRAC=$2
if [ $CAT == "pop" ]
then
    array=( "ACB" "ASW" "BEB" "CDX" "CEU" "CHB" "CHS" "CLM" "ESN" "FIN" "GBR" "GIH" "GWD" "IBS" "ITU" "JPT" "KHV" "LWK" "MSL" "MXL" "PEL" "PJL" "PUR" "STU" "TSI" "YRI" )
elif [ $CAT == "superpop" ]
then
    array=( "EUR" "AMR" "SAS" "AFR" "EAS" )
else
    echo "ERROR: invalid input: $1"
    echo "Usage $0 mode frac"
    echo "    mode: operating mode [\"pop/superpop\"]"
    echo "    frac: frequency cutoff, e.g. 2: major allele [INT]"
fi

module load gcc/5.5.0
module load vcftools

for s in "${array[@]}"
do
    ~/bin/bcftools-1.9 view -S $DIR_SAMPLE/samples_${CAT}_${s}.txt --force-samples $VCF -V mnps,other |\
        vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | \
        bgzip -@ 8 > ${CHROM}_${CAT}_${s}.vcf.gz
    ~/bin/bcftools-1.9 view -h ${CHROM}_${CAT}_${s}.vcf.gz | tail -1 > vcf_${CAT}_${s}_samples.txt
    #THRSD=$(( $(cat vcf_${CAT}_${s}_samples.txt | wc -w) - 9 ))
    THRSD=$(( ($(cat vcf_${CAT}_${s}_samples.txt | wc -w) - 9) * 2 / $FRAC ))
    FILTER="AC > $THRSD"
    #echo $FILTER
    ~/bin/bcftools-1.9 view -i "$FILTER" -v snps,indels ${CHROM}_${CAT}_${s}.vcf.gz > ${CHROM}_${CAT}_${s}_thrsd${FRAC}.vcf
done
