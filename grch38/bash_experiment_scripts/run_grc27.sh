set -x
# Usage: sh run_grc27.sh HG00103
# or parallel -j 5 -k "sh run_grc27.sh {1}" ::: `ls experiments`
# $1: sample name

THREADS=8
BT_INDEX_A="grc/indexes/chr21_grc"

ABBREV="grc27"
READS="simulation/${1}/chr21-per_1.fq" # "experiments/${1}/chr21-major-1kg-mapqlt10.fq" # experiments/HG00103/chr21-major-1kg-mapqlt10.fq
DIR="experiments/${1}"
RF_DIR="experiments/${1}/${ABBREV}"
mkdir -p $RF_DIR

if [[ -s experiments/${1}/chr21-${ABBREV}.acc_log ]]; then
    echo "${1} has already been processed. Skipped"
    exit
fi

# Align to personalized reference genomes using different random seeds
rm ${DIR}/${ABBREV}.path
rm ${DIR}/${ABBREV}.ids
parallel -j 3 -k "bowtie2 --reorder -p $THREADS -x ${BT_INDEX_A} -U $READS -S ${DIR}/chr21-${ABBREV}-seed{}.sam --seed {}" ::: `echo $(seq 0 26)`
for s in $(seq 0 26); do
    # bowtie2 --reorder -p $THREADS -x ${BT_INDEX_A} -U $READS -S ${DIR}/chr21-${ABBREV}-seed${s}.sam --seed ${s}
    ls ${DIR}/chr21-${ABBREV}-seed${s}.sam >> ${DIR}/${ABBREV}.path
    echo "seed${s}" >> ${DIR}/${ABBREV}.ids
done

# Merge pre-lifted results
python3.6 $NC2/refflow-exp/reference_flow/src/merge_incremental.py \
    -ns ${DIR}/${ABBREV}.path \
    -ids ${DIR}/${ABBREV}.ids \
    -rs 0 -p ${DIR}/chr21-${ABBREV}-merged -l ${DIR}/${ABBREV}.output_path
MERGE_LIST=""
for s in $(seq 0 26); do
    rm ${DIR}/chr21-${ABBREV}-seed${s}.sam
    MERGE_LIST="${MERGE_LIST} ${DIR}/chr21-${ABBREV}-merged-seed${s}.sam"
done

# Merge alignments to the same haplotype and perform liftover using levioSAM
samtools merge ${DIR}/chr21-${ABBREV}.sam $MERGE_LIST
for s in $(seq 0 26); do
    rm ${DIR}/chr21-${ABBREV}-merged-seed${s}.sam
done

python3.6 -O $NC2/refflow-exp/scripts/analyze_diploid_indels.py -c chr21 \
    -g simulation/${1}/chr21-per_1.sam -vr simulation/${1}/chr21-per.var \
    -n experiments/${1}/chr21-${ABBREV}.sam > experiments/${1}/chr21-${ABBREV}.acc_log

samtools view -hb -o ${DIR}/chr21-${ABBREV}.bam ${DIR}/chr21-${ABBREV}.sam
rm ${DIR}/chr21-${ABBREV}.sam
