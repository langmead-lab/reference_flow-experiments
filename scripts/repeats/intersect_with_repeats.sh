# Get fraction of biased sites overlapped with the entire rmsk annotated regions.
for i in `ls *-NA12878-biased_sites.bed`; do
    wc -l $i; bedtools intersect -a $i -b rmsk.bed | wc -l
done

# Get fraction of biased sites overlapping different repeats.
for j in GRC major randflow_ld randflow_ld_26 vg per; do
    echo $j
    ALL=`bedtools intersect -a $j-NA12878-biased_sites.bed -b rmsk.bed | wc -l`
    for i in `ls rmsk-*.bed`; do
        bedtools intersect -a $j-NA12878-biased_sites.bed -b $i > ${j}-${i}
        CURRENT=`cat $j-$i | wc -l`
        echo $i
        echo $(($CURRENT / ${ALL}.))
    done
done

# ALL=`bedtools intersect -a per-NA12878-biased_sites.bed -b rmsk.bed| wc -l`
# for i in `ls rmsk-*.bed`; do 
#     echo $i
#     echo $((`bedtools intersect -a per-NA12878-biased_sites.bed -b $i | wc -l` / ${ALL}.)) 
# done
