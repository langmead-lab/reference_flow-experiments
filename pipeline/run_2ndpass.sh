usage(){
    echo "Usage: $(basename $0) [-cCth] -i index_dir -f fastq"
    echo "------ Requirements -----"
    echo "  -i  directory of indexes"
    echo "  -f  input FASTQ file"
    echo "  -p  prefix for output SAM file"
    echo "------ Options -----"
    echo "  -b  size of block for stochastic genomes [1]"
    echo "  -c  name of chromosome [21]"
    echo "  -C  category {superpop} [superpop]"
    echo "  -s  path of scripts ['$'REL/scripts]"
    echo "  -t  population frequency threshold [0]"
}

#: default
CHROM="21"
CAT="superpop"
FRAC="0"
BLOCK_SIZE="1"
SCRIPTS="$REL/scripts"

while getopts b:c:C:s:t:i:f:p:h: option
do
    case "${option}"
    in
    #: options
    b) 
        BLOCK_SIZE=${OPTARG}
        echo "Set stochastic block size -> $BLOCK_SIZE"
        ;;
    c) 
        echo "Set chromosome name -> $CHROM"
        CHROM=${OPTARG}
        ;;
    C) 
        echo "Set population category -> $CAT"
        CAT=${OPTARG}
        ;;
    t) 
        echo "Set frequency cutoff level -> $FRAC"
        FRAC=${OPTARG}
        ;;
    s) 
        echo "Set script directory -> $SCRIPTS"
        SCRIPTS=${OPTARG}
        ;;

    #: requirements
    i) INDEX_DIR=${OPTARG};;
    f) FASTQ=${OPTARG};;
    p) OUT_PREFIX=${OPTARG};;

    h) usage;;
    *) 
        usage
        exit 1
    esac
done

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

if [[ -f merging.paths ]]
then
    echo "Warning: there's an existing merging.paths file"
    rm -i merging.paths
fi
if [[ -f merging.ids ]]
then
    echo "Warning: there's an existing merging.ids file"
    rm -i merging.ids
fi

#: alignment
for s in "${array[@]}"
do
    SAM="$OUT_PREFIX-${s}s${FRAC}_b${BLOCK_SIZE}.sam"
    #bowtie2 -p 8 -x $INDEX_DIR/${CHROM}_${CAT}_${s}_thrsd${FRAC}_stochastic_b${BLOCK_SIZE} -U $FASTQ -S $SAM --seed 0
    #: use 1 thread for experiments
    bowtie2 -x $INDEX_DIR/${CHROM}_${CAT}_${s}_thrsd${FRAC}_stochastic_b${BLOCK_SIZE} -U $FASTQ -S $SAM --seed 0
    echo `pwd`/$SAM >> merging.paths
    echo ${s}s${FRAC} >> merging.ids
done


if [[ -f tmp_merge.log ]]
then
    echo "Warning: there's an existing rmp_merge.log file"
    rm -i tmp_merge.log
fi

#: merge alignments
python $SCRIPTS/merge_sam.py -ns merging.paths -ids merging.ids -rs 0 >> tmp_merge.log
WC=`wc -l merging.ids | cut -d ' ' -f1`
ARR=(`tail -$(($WC+1)) tmp_merge.log | head -$(($WC))`)
for s in "${ARR[@]}"
do
    echo $s
done