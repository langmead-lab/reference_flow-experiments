if [ "$#" -ne "3" ]; then
    echo "Usage: $0 prefix vcf_path fasta_path" >&3
else
    msg=$(bcftools -h | wc -l)
    if [ "$msg" -gt "1" ]; then
        echo "bcftools has already been loaded";
    else
        echo "loading bcftools...";
        module load bcftools;
    fi
    msg=$(vcftools -h | wc -l)
    if [ "$msg" -gt "1" ]; then
        echo "vcftools has already been loaded";
    else
        echo "loading vcftools...";
        module load vcftools;
    fi
    msg=$(samtools -h | wc -l)
    if [ "$msg" -gt "1" ]; then
        echo "samtools has already been loaded";
    else
        echo "loading samtools...";
        module load samtools;
    fi
    PREFIX=$1
    VCF="${PREFIX}.vcf.gz"
    VCF_R="${PREFIX}.recode.vcf"
    VCF_R_GZ="${PREFIX}.recode.vcf.gz"
    MA_REF="${PREFIX}_major_allele.fa"
    NUM_TH=16
    #: bcftools dont solve mnps
    bcftools view -O z -v snps --threads $NUM_TH -q 0.5 $2 > $VCF
    #: vcftools need sample to calculate frequency (cannot only look at INFO:AF)
    # vcftools --gzvcf $2 --non-ref-af 0.5 --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out $PREFIX --remove-indels
    # bgzip $VCF -@ $NUM_TH

    #: [start] this works
    #vcftools --gzvcf $VCF --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out $PREFIX --remove-indels
    #bgzip $VCF_R -@ $NUM_TH
    #bcftools index $VCF_R_GZ
    #bcftools consensus -f $3 $VCF_R_GZ > $MA_REF
    #: [end]

    vcftools --gzvcf $VCF --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout --remove-indels | bgzip -@ $NUM_TH > $VCF_R_GZ
    bcftools index $VCF_R_GZ
    #bcftools consensus -f $3 $VCF_R_GZ > $MA_REF
    samtools faidx $3 $1 | bcftools consensus $VCF_GZ > $MA_REF
fi
