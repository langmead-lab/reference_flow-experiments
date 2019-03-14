if [ "$#" -ne "1" ]; then
    echo "Usage: $0 fasta" >&3
else
    for i in $(seq 1 22) X Y; do
        echo ${i}ing
        # bgzip -cd $GS/1kg_vcf/${i}.vcf.gz | python ../scripts/vcf_processing.py --out_var_loc 1 --min_af 0.5 --rand_th 0.01 > chr${i}_var_001.txt
        # python ../scripts/test_major_allele_ref.py -m wgs_gatk.fa -r $ARRAY/indexes/hg19.fa -c chr${i}_var_001.txt >> test_major_allele_ref.txt
        echo ${i} >> test_major_allele_ref.txt 
        bgzip -cd $GS/1kg_vcf/${i}.vcf.gz | \
            python $REL/scripts/vcf_processing.py --out_var_loc 1 --min_af 0.5 --rand_th 1 --max_allele_len 1| \
            python $REL/scripts/test_major_allele_ref.py -m $1 -r $ARRAY/indexes/hg19.fa  >> test_major_allele_ref.txt
    done
fi
