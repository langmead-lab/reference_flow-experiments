if [ "$#" -ne "4" ]; then
    echo "Usage: $0 prefix vcf fasta_prefix gatk_path" >&3
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
    VCF="${PREFIX}.vcf"
    VCF_GZ="${PREFIX}.vcf.gz"
    VCF_R="${PREFIX}.recode.vcf"
    VCF_R_GZ="${PREFIX}.recode.vcf.gz"
    MA_REF="${PREFIX}_major_allele.fa"
    NUM_TH=16
    #: bcftools dont solve mnps
    # bcftools view -O z -v snps --threads $NUM_TH -q 0.5 $2 > $VCF_GZ

    #: vcftools need sample to calculate frequency (cannot only look at INFO:AF)
    # vcftools --gzvcf $2 --non-ref-af 0.5 --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out $PREFIX --remove-indels
    # bgzip $VCF -@ $NUM_TH

    #: [start]
    #vcftools --gzvcf $VCF_GZ --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out $PREFIX --remove-indels
    #bgzip $VCF_R -@ $NUM_TH
    #bcftools index $VCF_R_GZ
    #bcftools consensus -f ${3}.fa $VCF_R_GZ > $MA_REF
    #: [end]

    bcftools view -O v -v snps --threads $NUM_TH -q 0.5 $2 > $VCF_GZ
    vcftools --gzvcf $VCF_GZ --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout --remove-indels | bgzip -@ $NUM_TH > $VCF_R_GZ
    bcftools index $VCF_R_GZ
    $4 CreateSequenceDictionary -R ${3}.fa -O ${3}.dict
    $4 IndexFeatureFile -F $VCF_R_GZ
    $4 FastaAlternateReferenceMaker -R ${3}.fa -O $MA_REF -V $VCF_R_GZ
    sed -iE 's/>[0-9+]* />/; s/:/ /' $MA_REF

    #: test major allele reference using chromosome vcfs
    #sh $REL/scripts/test_ma.sh $MA_REF
    for i in $(seq 1 22); do
        bgzip -cd $GS/1kg_vcf/${i}.vcf.gz | \
            python $REL/scripts/vcf_processing.py --out_var_loc 1 --min_af 0.5 --rand_th 1 | \
            python $REL/scripts/test_major_allele_ref.py -m $MA_REF -r ${3}.fa  >> test_major_allele_ref.txt
    done
fi
