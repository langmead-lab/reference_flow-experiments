usage(){
    echo "Usage: $(basename $0) [-cCth] -i fastq -I index_dir"
    echo "------ Requirements -----"
    echo "  -i  input FASTQ file"
    echo "  -I  directory of indexes"
    echo "  -p  prefix for output SAM file"
    echo "------ Options -----"
    echo "  -b  size of block for stochastic genomes [1]"
    echo "  -c  name of chromosome [21]"
    echo "  -C  category {superpop} [superpop]"
    echo "  -S  path of scripts ['$'REL/scripts]"
    echo "  -t  population frequency threshold [0]"
}

#: default
CHROM="21"
CAT="superpop"
FRAC="0"
BLOCK_SIZE="1"
SCRIPTS="$REL/scripts"

while getopts b:c:C:S:t:i:I:p:h: option
do
    case "${option}"
    in
    #: options
    b) 
        BLOCK_SIZE=${OPTARG}
        echo "Set stochastic block size -> $BLOCK_SIZE"
        ;;
    c) 
        CHROM=${OPTARG}
        echo "Set chromosome name -> $CHROM"
        ;;
    C) 
        CAT=${OPTARG}
        echo "Set population category -> $CAT"
        ;;
    S) 
        SCRIPTS=${OPTARG}
        echo "Set script directory -> $SCRIPTS"
        ;;
    t)
        FRAC=${OPTARG}
        echo "Set frequency cutoff level -> $FRAC"
        ;;

    #: requirements
    I) INDEX_DIR=${OPTARG};;
    i) FASTQ=${OPTARG};;
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

#: exit when any command fails
set -e

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
python3.7 $SCRIPTS/merge_sam.py -ns merging.paths -ids merging.ids -rs 0 >> tmp_merge.log
WC=`wc -l merging.ids | cut -d ' ' -f1`
ARR=(`tail -$(($WC+1)) tmp_merge.log | head -$(($WC))`)
for s in "${ARR[@]}"
do
    echo $s
done
