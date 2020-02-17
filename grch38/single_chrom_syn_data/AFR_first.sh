indiv='HG01974'

# ['NA19312', 'HG01974', 'NA20541', 'HG01863', 'HG03756']

bowtie2 --threads 20 -x ../../../../hg38_100/pop_genome/thrds0.5_S0_b1_ld0/indexes/chr21_superpop_AFR_thrds0.5_S0_b1_ld0 -U ../../../simulation/${indiv}/chr21-per_1.fq -S chr21-AFR.sam
awk -v var=10 '{{ if ($5 >= var || $1 ~ /^@/) {{ print }} }}' chr21-AFR.sam > chr21-AFR-mapqgeq10.sam
awk -v var=10 '{{ if ($5 < var || $1 ~ /^@/) {{ print }} }}' chr21-AFR.sam > chr21-AFR-mapqlt10.sam
samtools fastq chr21-AFR-mapqlt10.sam > chr21-AFR-mapqlt10.fq

for i in EUR AMR SAS EAS; do bowtie2 --reorder --threads 20 -x ../../../../hg38_100/pop_genome/thrds0.5_S0_b1_ld0/indexes/chr21_superpop_${i}_thrds0.5_S0_b1_ld0 -U chr21-AFR-mapqlt10.fq -S chr21-AFR-10-${i}.sam; done
bowtie2 --threads 20 --reorder -x ../../../major/indexes/chr21-major-1kg -U chr21-AFR-mapqlt10.fq -S chr21-AFR-10-maj.sam

for i in AFR EUR AMR SAS EAS maj; do 
    echo ${i} >> chr21-AFR_first.ids
done

ls chr21-AFR-mapqlt10.sam >> chr21-AFR_first.paths
for i in EUR AMR SAS EAS maj; do 
    ls chr21-AFR-10-${i}.sam >> chr21-AFR_first.paths
done

python $REL2/scripts/merge_incremental.py -ns chr21-AFR_first.paths -ids chr21-AFR_first.ids -rs 0 -p 2ndpass-AFR_first -l chr21-AFR_first.merge_paths

$NC2/liftover/liftover lift -a chr21-AFR-mapqgeq10.sam -l ../../../../hg38_100/pop_genome/thrds0.5_S0_b1_ld0/chr21_superpop_AFR_thrds0.5_S0_b1_ld0.lft -p chr21-AFR-mapqgeq10-liftover
for i in AFR EUR AMR SAS EAS; do 
    $NC2/liftover/liftover lift -a 2ndpass-AFR_first-${i}.sam -l ../../../../hg38_100/pop_genome/thrds0.5_S0_b1_ld0/chr21_superpop_${i}_thrds0.5_S0_b1_ld0.lft -p 2ndpass-AFR_first-${i}-liftover
done
$NC2/liftover/liftover lift -a 2ndpass-AFR_first-maj.sam -l ../../../major/chr21-major-1kg.lft -p 2ndpass-AFR_first-maj-liftover

cp chr21-AFR-mapqgeq10-liftover.sam chr21-AFR_first-refflow.sam
for i in AFR EUR AMR SAS EAS maj; do 
    cat 2ndpass-AFR_first-${i}-liftover.sam >> chr21-AFR_first-refflow.sam
done

python -O $REL2/scripts/analyze_diploid_indels.py -c 21 -g ../../../simulation/${indiv}/chr21-per_1.sam -vr ../../../simulation/${indiv}/chr21-per.var -p 0 -n chr21-AFR_first-refflow.sam
