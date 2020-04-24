# By Nae-Chyun Chen, 2020
# 
# This script compares two FASTA files and report if there are differences.
# By default we compare chr1-22, chrX and chrY.
# Temp files are deleted if contigs from two files are identical.
# Final output is stored as a log file with a timestamp.
#
# It is suggested that the FASTA files have no other chars than ACGTN.
#

NEW=$1 #"wg-major.fa"
OLD=$2 #"wg-maj-old.fa"
LOGF=`date +%Y_%m_%d_%H_%m_%S`
echo "$0 $1 $2" > ${LOGF}.diff.log

for ii in  $(seq 1 22) X Y; do
i=chr${ii};
samtools faidx $NEW $i > new_${i}.fa
samtools faidx $OLD $i > old_${i}.fa
diff new_${i}.fa old_${i}.fa > ${i}.diff.log
if [[ `wc -l <${i}.diff.log` -eq 0 ]]; # if diff file is empty
then
    echo "${i} passed" >> ${LOGF}.diff.log
    rm new_${i}.fa old_${i}.fa ${i}.diff.log
else
    echo "============"
    echo "${i} failed" >> ${LOGF}.diff.log
    echo "============"
fi
done
