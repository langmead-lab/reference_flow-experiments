# set -x
# Usage: sh run_refflow_26.sh HG00103
# or parallel -j 5 -k "sh run_refflow_26.sh {1}" ::: `ls experiments`
# $1: sample name

THREADS=16
# RandFlow
METHOD="thrds0_S1_b1_ld0"
# Bowtie2 index: {BT_INDEX_PREFIX}-${POP}-{BT_INDEX_SUFFIX}
BT_INDEX_PREFIX="pop_genome26/${METHOD}/indexes/chr21"
BT_INDEX_SUFFIX="randflow"

# RandFlow-LD
# METHOD="thrds0_S1_b1000_ld1"
# # Bowtie2 index: {BT_INDEX_PREFIX}-${POP}-{BT_INDEX_SUFFIX}
# BT_INDEX_PREFIX="pop_genome26/${METHOD}/indexes/chr21"
# BT_INDEX_SUFFIX="randflow_ld"

ABBREV="${BT_INDEX_SUFFIX}_26"
READS="experiments/${1}/chr21-major-1kg-mapqlt10.fq" # experiments/HG00103/chr21-major-1kg-mapqlt10.fq
DIR="experiments/${1}"
RF_DIR="experiments/${1}/${ABBREV}"
mkdir -p $RF_DIR

if [[ -s experiments/${1}/chr21-${ABBREV}.acc_log ]]; then
    echo "${1} has already been processed. Skipped"
    exit
fi

# Second pass
ls ${DIR}/chr21-major-1kg-mapqlt10.sam > ${RF_DIR}/${ABBREV}.path
echo "major" > ${RF_DIR}/${ABBREV}.ids

for pop in CHB JPT CHS CDX KHV CEU TSI FIN GBR IBS YRI LWK GWD MSL ESN ASW ACB MXL PUR CLM PEL GIH PJL BEB STU ITU; do
    bowtie2 --reorder -p $THREADS -x ${BT_INDEX_PREFIX}-${pop}-${BT_INDEX_SUFFIX} -U $READS -S ${RF_DIR}/raw-${pop}.sam
    ls ${RF_DIR}/raw-${pop}.sam >> ${RF_DIR}/${ABBREV}.path
    echo "${pop}" >> ${RF_DIR}/${ABBREV}.ids
done
python3.6 $NC2/refflow-exp/reference_flow/src/merge_incremental.py \
    -ns ${RF_DIR}/${ABBREV}.path \
    -ids ${RF_DIR}/${ABBREV}.ids \
    -rs 0 -p ${RF_DIR}/second_pass -l ${RF_DIR}/${ABBREV}.output_path

cp ${DIR}/chr21-major-1kg-mapqgeq10-liftover.sam ${DIR}/chr21-${ABBREV}.sam
for pop in major ; do
    leviosam lift -a ${RF_DIR}/second_pass-${pop}.sam \
        -l major/chr21-major-1kg.lft -p ${RF_DIR}/second_pass-${pop}-lifted -t $THREADS
    grep -hv "^@" ${RF_DIR}/second_pass-${pop}-lifted.sam >> ${DIR}/chr21-${ABBREV}.sam
done

for pop in CHB JPT CHS CDX KHV CEU TSI FIN GBR IBS YRI LWK GWD MSL ESN ASW ACB MXL PUR CLM PEL GIH PJL BEB STU ITU; do
    leviosam lift -a ${RF_DIR}/second_pass-${pop}.sam \
        -l pop_genome26/${METHOD}/chr21-${pop}-${BT_INDEX_SUFFIX}.lft \
        -p ${RF_DIR}/second_pass-${pop}-lifted -t $THREADS
    grep -hv "^@" ${RF_DIR}/second_pass-${pop}-lifted.sam >> ${DIR}/chr21-${ABBREV}.sam
done

python3.6 -O $NC2/refflow-exp/scripts/analyze_diploid_indels.py -c chr21 \
    -g simulation/${1}/chr21-per_1.sam -vr simulation/${1}/chr21-per.var \
    -n experiments/${1}/chr21-${ABBREV}.sam > experiments/${1}/chr21-${ABBREV}.acc_log

samtools view -hb -o ${DIR}/chr21-${ABBREV}.bam ${DIR}/chr21-${ABBREV}.sam
rm ${DIR}/chr21-${ABBREV}.sam
