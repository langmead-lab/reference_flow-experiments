# This script merges multiple VCF.gz files to create a joint variant set.
# It is useful when performing experiments that we'd like to use multiple
# pan-genomic methods based on the same set of variants.

set -x

THREADS=16
# Adjust this variable for different VCF sets.
OUTPUT_PREFIX=wg-thrds0_S1_b1000_ld1

if [ -f ${OUTPUT_PREFIX}-merged.vcf ]; then
    echo "${OUTPUT_PREFIX}-merged.vcf exists. Remove the file and re-run this"
    echo "script."
    exit
fi

for FILE in "$@"; do
    if [ ! -f ${FILE}.gz ]; then
        bgzip -c $FILE > ${FILE}.gz
    fi
    if [ ! -f ${FILE}.gz.csi ]; then
        bcftools index ${FILE}.gz
    fi
    echo ${FILE}*
done

# Add ".gz" prefixes to VCF files.
set -- "${@/%/.gz}"
bcftools merge "${@}" --force-samples --threads $THREADS > \
    ${OUTPUT_PREFIX}-merged.vcf
