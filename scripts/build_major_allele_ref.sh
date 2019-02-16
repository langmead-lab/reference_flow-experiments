if [ "$#" -ne "3" ]; then
    echo "Usage: $0 prefix vcf_path index_path" >&3
else
    msg=$(bcftools -h | wc -l)
    if [ "$msg" -gt "1" ]; then
        echo "bcftools has already been loaded";
    else
        echo "loading bcftools...";
        module load bcftools;
    fi
    PREFIX=$1
    VCF="${PREFIX}.vcf.gz"
    MA_REF="${PREFIX}_major_allele.fa"
    bcftools view -O z -v snps,indels --threads 8 -q 0.5 $2 > $VCF
    bcftools index $VCF
    bcftools consensus -f $3 $VCF > $MA_REF
fi
