display_usage() {
    echo "This script simulates reads from a diploid genome using mason_simulator"
    echo ""
    echo "Usage: $0 genome_hapA genome_hapB num_sim_read out_dir"
    echo "    genome_hapA: fasta file for the first haplotype"
    echo "    genome_hapB: fasta file for the second haplotype"
    echo "    num_sim_read: number of reads to simulate for each haplotype"
	echo "    out_dir: directory to write output files to"
}

MASON=/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/software/mason_simulator

#Editted to write files to specified output directory

mkdir -p "$4"

if [ "$#" -ne "4" ]; then
    display_usage
    exit
else
    set -x
    FILENAMEA=`echo $1 | rev | cut -d '/' -f1 | rev`
    PREFIXA=`echo ${FILENAMEA} | cut -d. -f1`
	OUTFIXA=$4${PREFIXA}
    $MASON -ir $1 -n $3 -o ${OUTFIXA}_1.fq -or ${OUTFIXA}_2.fq -oa ${OUTFIXA}.sam --read-name-prefix "${OUTFIXA}_simulated."

    FILENAMEB=`echo $2 | rev | cut -d '/' -f1 | rev`
    PREFIXB=`echo ${FILENAMEB} | cut -d. -f1`
    $MASON -ir $2 -n $3 -o ${OUTFIXB}_1.fq -or ${OUTFIXB}_2.fq -oa ${OUTFIXB}.sam --read-name-prefix "${OUTFIXB}_simulated."

    SAMPLEA=`cut -d_ -f1 <<< ${PREFIXA}`
    SAMPLEB=`cut -d_ -f1 <<< ${PREFIXB}`
    if [ $SAMPLEA != $SAMPLEB ]; then
        echo "Error: sample A (${SAMPLEA}) doesn't match sample B (${SAMPLEB})"
        exit
    else
        SAMPLE=$4${SAMPLEA}
        cat ${OUTFIXA}_1.fq > ${SAMPLE}_1.fq
        cat ${OUTFIXB}_1.fq >> ${SAMPLE}_1.fq
        cat ${OUTFIXA}_2.fq > ${SAMPLE}_2.fq
        cat ${OUTFIXB}_2.fq >> ${SAMPLE}_2.fq
        #cat ${PREFIXA}.sam > ${SAMPLE}.sam
        #cat ${PREFIXB}.sam >> ${SAMPLE}.sam
        grep ^@ ${OUTFIXA}.sam > ${SAMPLE}.sam
        grep ^@ ${OUTFIXB}.sam | tail -n +2 >> ${SAMPLE}.sam
        grep -v ^@ ${OUTFIXA}.sam >> ${SAMPLE}.sam
        grep -v ^@ ${OUTFIXB}.sam >> ${SAMPLE}.sam
    fi
fi
