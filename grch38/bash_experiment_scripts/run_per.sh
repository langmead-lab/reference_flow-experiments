# set -x
# Usage: sh run_per.sh HG00103
# or parallel -j 5 -k "sh run_per.sh {1}" ::: `ls experiments`
# $1: sample name

THREADS=8
BT_INDEX_A="simulation/${1}/indexes/chr21-per_hapA-no_hap"
BT_INDEX_B="simulation/${1}/indexes/chr21-per_hapB-no_hap"
#BT_INDEX_A="/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/refflow-exp/snakemake/sim_chr21/personalized/${1}/indexes/chr21-perA"
#BT_INDEX_B="/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/refflow-exp/snakemake/sim_chr21/personalized/${1}/indexes/chr21-perB"

ABBREV="per10"
READS="simulation/${1}/chr21-per_1.fq" # "experiments/${1}/chr21-major-1kg-mapqlt10.fq" # experiments/HG00103/chr21-major-1kg-mapqlt10.fq
DIR="experiments/${1}"
RF_DIR="experiments/${1}/${ABBREV}"
mkdir -p $RF_DIR

if [[ -s experiments/${1}/chr21-${ABBREV}.acc_log ]]; then
    echo "${1} has already been processed. Skipped"
    exit
fi

# Pre-processing; build indexes without haplotype labels
sed -e 's/^>chr21A/>chr21/g' simulation/${1}/chr21-per_hapA.fa > simulation/${1}/chr21-per_hapA-no_hap.fa
bowtie2-build --threads 16 simulation/${1}/chr21-per_hapA-no_hap.fa simulation/${1}/indexes/chr21-per_hapA-no_hap
sed -e 's/^>chr21B/>chr21/g' simulation/${1}/chr21-per_hapB.fa > simulation/${1}/chr21-per_hapB-no_hap.fa
bowtie2-build --threads 16 simulation/${1}/chr21-per_hapB-no_hap.fa simulation/${1}/indexes/chr21-per_hapB-no_hap

# Align to personalized reference genomes using different random seeds
rm ${DIR}/${ABBREV}.path
rm ${DIR}/${ABBREV}.ids
for s in 0 1 2 3 4; do
    bowtie2 --reorder -p $THREADS -x ${BT_INDEX_A} -U $READS -S ${DIR}/chr21-${ABBREV}-hapA-seed${s}.sam --seed ${s}
    ls ${DIR}/chr21-${ABBREV}-hapA-seed${s}.sam >> ${DIR}/${ABBREV}.path
    echo "hapA-seed${s}" >> ${DIR}/${ABBREV}.ids
    bowtie2 --reorder -p $THREADS -x ${BT_INDEX_B} -U $READS -S ${DIR}/chr21-${ABBREV}-hapB-seed${s}.sam --seed ${s}
    ls ${DIR}/chr21-${ABBREV}-hapB-seed${s}.sam >> ${DIR}/${ABBREV}.path
    echo "hapB-seed${s}" >> ${DIR}/${ABBREV}.ids
done

# Merge pre-lifted results
python3.6 $NC2/refflow-exp/reference_flow/src/merge_incremental.py \
    -ns ${DIR}/${ABBREV}.path \
    -ids ${DIR}/${ABBREV}.ids \
    -rs 0 -p ${DIR}/chr21-${ABBREV}-merged -l ${DIR}/${ABBREV}.output_path
for s in 0 1 2 3 4; do
    rm ${DIR}/chr21-${ABBREV}-hapA-seed${s}.sam
    rm ${DIR}/chr21-${ABBREV}-hapB-seed${s}.sam
done

# Merge alignments to the same haplotype and perform liftover using levioSAM
for h in A B; do
    samtools merge ${DIR}/chr21-${ABBREV}-hap${h}.sam ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed0.sam \
        ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed1.sam ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed2.sam \
        ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed3.sam ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed4.sam
    rm ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed0.sam ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed1.sam ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed2.sam ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed3.sam ${DIR}/chr21-${ABBREV}-merged-hap${h}-seed4.sam
    # leviosam lift -a ${DIR}/chr21-${ABBREV}-hap${h}.sam -l /net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/refflow-exp/snakemake/sim_chr21/personalized/${1}/chr21-per${h}.lft \
    leviosam lift -a ${DIR}/chr21-${ABBREV}-hap${h}.sam -l simulation/${1}/chr21-per${h}.lft \
        -p ${DIR}/chr21-${ABBREV}-hap${h}-lifted -t $THREADS
    rm ${DIR}/chr21-${ABBREV}-hap${h}.sam
done

samtools merge ${DIR}/chr21-${ABBREV}-lifted.sam ${DIR}/chr21-${ABBREV}-hapA-lifted.sam ${DIR}/chr21-${ABBREV}-hapB-lifted.sam
rm ${DIR}/chr21-${ABBREV}-hapA-lifted.sam ${DIR}/chr21-${ABBREV}-hapB-lifted.sam

python3.6 -O $NC2/refflow-exp/scripts/analyze_diploid_indels.py -c chr21 \
    -g simulation/${1}/chr21-per_1.sam -vr simulation/${1}/chr21-per.var \
    -n experiments/${1}/chr21-${ABBREV}-lifted.sam > experiments/${1}/chr21-${ABBREV}.acc_log

samtools view -hb -o ${DIR}/chr21-${ABBREV}-lifted.bam ${DIR}/chr21-${ABBREV}-lifted.sam
rm ${DIR}/chr21-${ABBREV}-lifted.sam
