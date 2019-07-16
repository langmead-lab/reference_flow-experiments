display_usage() {
    echo "This script simulates reads from a diploid genome using mason_simulator"
    echo ""
    echo "Usage: $0 genome_hapA genome_hapB num_sim_read"
    echo "    genome_hapA: fasta file for the first haplotype"
    echo "    genome_hapB: fasta file for the second haplotype"
    echo "    num_sim_read: number of reads to simulate for each haplotype"
}
if [ "$#" -ne "3" ]; then
    display_usage
    exit
else
    set -x
    FILENAMEA=`echo $1 | rev | cut -d '/' -f1 | rev`
    PREFIXA=`echo ${FILENAMEA} | cut -d. -f1`
    mason_simulator -ir $1 -n $3 -o ${PREFIXA}_1.fq -or ${PREFIXA}_2.fq -oa ${PREFIXA}.sam --read-name-prefix "${PREFIXA}_simulated."

    FILENAMEB=`echo $2 | rev | cut -d '/' -f1 | rev`
    PREFIXB=`echo ${FILENAMEB} | cut -d. -f1`
    mason_simulator -ir $2 -n $3 -o ${PREFIXB}_1.fq -or ${PREFIXB}_2.fq -oa ${PREFIXB}.sam --read-name-prefix "${PREFIXB}_simulated."

    SAMPLEA=`cut -d_ -f1 <<< ${PREFIXA}`
    SAMPLEB=`cut -d_ -f1 <<< ${PREFIXB}`
    if [ $SAMPLEA != $SAMPLEB ]; then
        echo "Error: sample A (${SAMPLEA}) doesn't match sample B (${SAMPLEB})"
        exit
    else
        SAMPLE=$SAMPLEA
        cat ${PREFIXA}_1.fq > ${SAMPLE}_1.fq
        cat ${PREFIXB}_1.fq >> ${SAMPLE}_1.fq
        cat ${PREFIXA}_2.fq > ${SAMPLE}_2.fq
        cat ${PREFIXB}_2.fq >> ${SAMPLE}_2.fq
        #cat ${PREFIXA}.sam > ${SAMPLE}.sam
        #cat ${PREFIXB}.sam >> ${SAMPLE}.sam
        grep ^@ ${PREFIXA}.sam > ${SAMPLE}.sam
        grep ^@ ${PREFIXB}.sam | tail -n +2 >> ${SAMPLE}.sam
        grep -v ^@ ${PREFIXA}.sam >> ${SAMPLE}.sam
        grep -v ^@ ${PREFIXB}.sam >> ${SAMPLE}.sam
    fi
fi
