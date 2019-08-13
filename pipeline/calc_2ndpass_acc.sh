usage(){
    echo "Usage: $(basename $0) [-CpS] -c chrom -g gold_sam -V path_var_reads -P path_var_genome"
    echo "------ Requirements -----"
    echo "  -c  name of chromosome"
    echo "  -g  gold SAM file provided by simulator"
    echo "  -V  path of VAR for simulated individual (reads)"
    echo "  -P  path of VAR for genomes to be aligned"
    echo "------ Options -----"
    echo "  -C  category {superpop/pop} [superpop]"
    echo "  -p  number of threads [8]"
    echo "  -S  path of scripts ['$'R]"
    exit
}

SCRIPTS="/scratch/groups/blangme2/naechyun/relaxing/scripts"
THREADS=8
CAT="superpop"

# VAR="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline"
# VAR="/storage2/naechyun/refflow/chr9/pipeline"
#POP="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/experiments/stochastic_seed0"
#GOLD="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/reads_mason2/NA12878_1-h37maj-mapql10.sam"

while getopts c:g:V:P:C:p:S:h: option
do
    case "${option}"
    in
    #: options
    C) 
        CAT=${OPTARG}
        echo "Set category -> $CAT"
        ;;
    p) 
        THREADS=${OPTARG}
        echo "Set number of threads -> $THREADS"
        ;;
    S) 
        SCRIPTS=${OPTARG}
        echo "Set script directory -> $SCRIPTS"
        ;;

    #: requirements
    c) CHROM=${OPTARG};;
    g) GOLD=${OPTARG};;
    V) VAR=${OPTARG};;
    P) POP=${OPTARG};;

    #: help
    h)
        usage;;
    *)
        usage
    esac
done

#: check if required fields are given
if [ -z ${CHROM+x} ]
then
    echo "error: required input chrom (-c) is not set"
    usage
fi
if [ -z ${GOLD+x} ]
then
    echo "error: required input gold (-g) is not set"
    usage
fi
if [ -z ${VAR+x} ]
then
    echo "error: required input var (-V) is not set"
    usage
fi
if [ -z ${POP+x} ]
then
    echo "error: required input population var (-P) is not set"
    usage
fi

#: set population
if [ $CAT == "pop" ]
then
    array=( "ACB" "ASW" "BEB" "CDX" "CEU" "CHB" "CHS" "CLM" "ESN" "FIN" "GBR" "GIH" "GWD" "IBS" "ITU" "JPT" "KHV" "LWK" "MSL" "MXL" "PEL" "PJL" "PUR" "STU" "TSI" "YRI" )
elif [ $CAT == "superpop" ]
then
    #: for testing purpose
    #array=( "EUR" )
    array=( "EUR" "AMR" "SAS" "AFR" "EAS" )
else
    echo "ERROR: invalid input: $CAT"
    usage
fi

TP=0
SUM=0

update(){
    #grep -Eo '[0-9]+'
    TP=$(( $TP + `echo $TMP | grep -Eo '[0-9]+' | tail -2 | head -1`))
    SUM=$(( $SUM + `echo $TMP | grep -Eo '[0-9]+' | tail -1`))
    echo "TP=$TP"
    echo "SUM=$SUM"
}

WC=`wc -l merging.ids | cut -d ' ' -f1`
ARR=(`tail -$(($WC+1)) tmp_merge.log | head -$(($WC))`)
for i in $(seq 0 4)
do
    TMP=`python3.7 -O $SCRIPTS/analyze_diploid_indels.py \
        -c $CHROM -g $GOLD -p 0 -vr $VAR \
        -vs $POP/${CHROM}_${CAT}_${array[$i]}*.var -n ${ARR[$i]} \
        | grep 'sensitivity_all'`
    # TMP=`python3.7 -O $REL/scripts/analyze_diploid_indels.py -c $CHROM -g $GOLD -p 0 -vr $VAR/${SAMPLE}.var -vs $POP/*${array[$i]}*.var -n ${ARR[$i]} | grep 'sensitivity_all'`
    # TMP=`python3.7 -O $REL/scripts/analyze_diploid_indels.py -c $CHROM -g $GOLD -p 0 -vr $VAR/${SAMPLE}.var -vs $POP/${CHROM}_superpop_"${array[$i]}"_thrsd0_stochastic.var -n ${ARR[$i]} | grep 'sensitivity_all'`
    update
done

echo "RECALL = `bc -l <<< "$TP / $SUM"`"
