# set -x
# Usage: sh run_refflow_rand.sh HG00103 1
# or parallel -j 5 -k "sh run_refflow_rand.sh {1} {2}" ::: `ls experiments` ::: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
# $1: sample name
# $2: random seed

THREADS=4
# Bowtie2 index: BT_INDEX_PREFIX + $POP + BT_INDEX_SUFFIX
BT_INDEX_PREFIX="pop_genome/thrds0_S1_b1000_ld1-rand${2}/indexes/chr21-"
BT_INDEX_SUFFIX="-randflow_ld"

READS="experiments/${1}/chr21-major-1kg-mapqlt10.fq" # experiments/HG00103/chr21-major-1kg-mapqlt10.fq
DIR="experiments/${1}"
RF_DIR="experiments/${1}/randflow_ld-${2}"
mkdir -p $RF_DIR

if [[ -s experiments/${1}/chr21-randflow_ld-rand${2}.acc_log ]]; then
    echo "${1}/${2} has already been processed. Skipped"
    exit
fi

# Second pass
ls ${DIR}/chr21-major-1kg-mapqlt10.sam > ${RF_DIR}/randflow_ld.path
echo "major" > ${RF_DIR}/randflow_ld.ids

for pop in AFR AMR EAS EUR SAS; do
    bowtie2 --reorder -p $THREADS -x ${BT_INDEX_PREFIX}${pop}${BT_INDEX_SUFFIX} -U $READS -S ${RF_DIR}/raw-${pop}.sam
    ls ${RF_DIR}/raw-${pop}.sam >> ${RF_DIR}/randflow_ld.path
    echo "${pop}" >> ${RF_DIR}/randflow_ld.ids
done
python3.6 $NC2/refflow-exp/reference_flow/src/merge_incremental.py \
    -ns ${RF_DIR}/randflow_ld.path \
    -ids ${RF_DIR}/randflow_ld.ids \
    -rs 0 -p ${RF_DIR}/second_pass -l ${RF_DIR}/randflow_ld.output_path

cp ${DIR}/chr21-major-1kg-mapqgeq10-liftover.sam ${DIR}/chr21-randflow_ld-rand${2}.sam
for pop in major ; do
    leviosam lift -a ${RF_DIR}/second_pass-${pop}.sam \
        -l major/chr21-major-1kg.lft -p ${RF_DIR}/second_pass-${pop}-lifted -t $THREADS
    grep -hv "^@" ${RF_DIR}/second_pass-${pop}-lifted.sam >> ${DIR}/chr21-randflow_ld-rand${2}.sam
done

for pop in AFR AMR EAS EUR SAS; do
    leviosam lift -a ${RF_DIR}/second_pass-${pop}.sam \
        -l pop_genome/thrds0_S1_b1000_ld1-rand${2}/chr21-${pop}-randflow_ld.lft \
        -p ${RF_DIR}/second_pass-${pop}-lifted -t $THREADS
    grep -hv "^@" ${RF_DIR}/second_pass-${pop}-lifted.sam >> ${DIR}/chr21-randflow_ld-rand${2}.sam
done

python3.6 -O $NC2/refflow-exp/scripts/analyze_diploid_indels.py -c chr21 \
    -g simulation/${1}/chr21-per_1.sam -vr simulation/${1}/chr21-per.var \
    -n experiments/${1}/chr21-randflow_ld-rand${2}.sam > experiments/${1}/chr21-randflow_ld-rand${2}.acc_log
