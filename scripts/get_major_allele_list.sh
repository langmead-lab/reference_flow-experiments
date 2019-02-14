if [ "$#" -ne "2" ]; then
    echo "Usage: $0 vcf out_prefix" >&2
else
    case $1 in
        *gz) 
            vcftools --freq --gzvcf $1 --non-ref-af 0.5 --out $2
            ;;
        *vcf)
            vcftools --freq --vcf $1 --non-ref-af 0.5 --out $2
            ;;
    esac
fi
