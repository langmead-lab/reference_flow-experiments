# METHOD="vg_from_randflow_ld_six_p13"
# parallel -j 5 -k "$NCA/software/vg/vg map -d vg_from_randflow_ld/six_genomes/chr21-randflow_ld -f simulation/{1}/chr21-per_1.fq -t 16 --surject-to bam | samtools sort -@ 4 -O sam > experiments/{1}/chr21-${METHOD}-sorted.sam" ::: `ls experiments`
# parallel -j 20 -k "python -O $NC2/refflow-exp/scripts/analyze_diploid_indels.py -c chr21 -g simulation/{1}/chr21-per_1.sam -vr simulation/{1}/chr21-per.var -n experiments/{1}/chr21-${METHOD}-sorted.sam > experiments/{1}/chr21-${METHOD}.acc_log" ::: `ls experiments`
# parallel -j 30 -k "python summarize_mapping_acc.py -i experiments/{1}/chr21-${METHOD}.acc_log -o results/{1}-chr21-${METHOD}.acc" ::: `ls experiments`
# METHOD="vg_from_randflow_ld_p13"
# # parallel -j 5 -k "$NCA/software/vg/vg map -d vg_from_randflow_ld/p13/chr21-thrds0_S1_b1000_ld1-merged-p13 -f simulation/{1}/chr21-per_1.fq -t 16 --surject-to bam | samtools sort -@ 4 -O sam > experiments/{1}/${METHOD}-sorted.sam" ::: `ls experiments`
# # parallel -j 16 -k "samtools sort -@ 4 -O sam -o experiments/{1}/chr21-vg_from_randflow_ld-sorted.sam experiments/{1}/chr21-vg_from_randflow_ld.bam" ::: `ls experiments`
# parallel -k "mv experiments/{1}/${METHOD}-sorted.sam experiments/{1}/chr21-${METHOD}-sorted.sam" ::: `ls experiments`
# parallel -j 35 -k "python -O $NC2/refflow-exp/scripts/analyze_diploid_indels.py -c chr21 -g simulation/{1}/chr21-per_1.sam -vr simulation/{1}/chr21-per.var -n experiments/{1}/chr21-${METHOD}-sorted.sam > experiments/{1}/chr21-${METHOD}.acc_log" ::: `ls experiments`
# parallel -j 30 -k "python summarize_mapping_acc.py -i experiments/{1}/chr21-${METHOD}.acc_log -o results/{1}-chr21-${METHOD}.acc" ::: `ls experiments`
# METHOD="vg_0.1"
# parallel -j 5 -k "$NCA/software/vg/vg map -d vg/21-0.1 -f simulation/{1}/chr21-per_1.fq -t 16 --surject-to bam > experiments/{1}/chr21-${METHOD}.bam" ::: `ls experiments`
# parallel -j 16 -k "samtools sort -@ 4 -O sam -o experiments/{1}/chr21-${METHOD}-sorted.sam experiments/{1}/chr21-${METHOD}.bam" ::: `ls experiments`
# parallel -j 16 -k "python -O $NC2/refflow-exp/scripts/analyze_diploid_indels.py -c chr21 -g simulation/{1}/chr21-per_1.sam -vr simulation/{1}/chr21-per.var -n experiments/{1}/chr21-${METHOD}-sorted.sam > experiments/{1}/chr21-${METHOD}.acc_log" ::: `ls experiments`
# parallel -j 30 -k "python summarize_mapping_acc.py -i experiments/{1}/chr21-${METHOD}.acc_log -o results/{1}-chr21-${METHOD}.acc" ::: `ls experiments`
# METHOD="vg_linear_p13"
# parallel -j 5 -k "$NCA/software/vg/vg map -d vg_linear/p13/linear_p13 -f simulation/{1}/chr21-per_1.fq -t 16 --surject-to bam | samtools sort -@ 4 -O sam > experiments/{1}/chr21-${METHOD}-sorted.sam" ::: `ls experiments`
# parallel -j 16 -k "python -O $NC2/refflow-exp/scripts/analyze_diploid_indels.py -c chr21 -g simulation/{1}/chr21-per_1.sam -vr simulation/{1}/chr21-per.var -n experiments/{1}/chr21-${METHOD}-sorted.sam > experiments/{1}/chr21-${METHOD}.acc_log" ::: `ls experiments`
# parallel -j 30 -k "python summarize_mapping_acc.py -i experiments/{1}/chr21-${METHOD}.acc_log -o results/{1}-chr21-${METHOD}.acc" ::: `ls experiments`
# python summarize_as_tsv.py --fn_family $NC2/refflow-exp/reference_flow/resources/20130606_g1k.ped --fn_spop $NC2/refflow-exp/reference_flow/resources/1kg.superpopulation --dir_acc results --output_tsv acc.tsv
# METHOD="hisat_0.1"
# parallel -j 5 -k "hisat2 -p 16 -x hisat2/indexes/chr21-0.1 -k 10 --no-spliced-alignment --no-temp-splicesite -U simulation/{1}/chr21-per_1.fq | samtools view -h -F 256 | samtools sort -@ 4 -O sam -o experiments/{1}/chr21-${METHOD}-sorted.sam" ::: `ls experiments`
# parallel -j 16 -k "python -O $NC2/refflow-exp/scripts/analyze_diploid_indels.py -c chr21 -g simulation/{1}/chr21-per_1.sam -vr simulation/{1}/chr21-per.var -n experiments/{1}/chr21-${METHOD}-sorted.sam > experiments/{1}/chr21-${METHOD}.acc_log" ::: `ls experiments`
# parallel -j 30 -k "python summarize_mapping_acc.py -i experiments/{1}/chr21-${METHOD}.acc_log -o results/{1}-chr21-${METHOD}.acc" ::: `ls experiments`

# Different random seeds
# for i in $(seq 0 14); do
#     METHOD="randflow_ld-rand${i}"
#     parallel -j 30 -k "python summarize_mapping_acc.py -i experiments/{1}/chr21-${METHOD}.acc_log -o results/{1}-chr21-${METHOD}.acc" ::: `ls experiments`
# done

METHOD="$1" #"per10" # "randflow_26" # "randflow_ld_26"
# parallel -j 25 -k "python summarize_mapping_acc.py -i experiments/{1}/chr21-${METHOD}.acc_log -o results/num_incorrect/{1}-chr21-${METHOD}.num_incorrect -m num_incorrect" ::: `ls experiments`
parallel -j 25 -k "python summarize_mapping_acc.py -i experiments/{1}/chr21-${METHOD}.acc_log -o results/num_unaligned/{1}-chr21-${METHOD}.unaligned -m unaligned" ::: `ls experiments`

# Using population frequencies instead of super population
# METHOD="randflow" #"per10" # "randflow_26" # "randflow_ld_26"
# parallel -j 25 -k "python summarize_mapping_acc.py -i experiments/{1}/chr21-${METHOD}.acc_log -o results/{1}-chr21-${METHOD}.acc" ::: `ls experiments`
# python summarize_as_tsv.py --fn_family $NC2/refflow-exp/reference_flow/resources/20130606_g1k.ped --fn_spop $NC2/refflow-exp/reference_flow/resources/1kg.superpopulation --dir_acc results --output_tsv acc.tsv
