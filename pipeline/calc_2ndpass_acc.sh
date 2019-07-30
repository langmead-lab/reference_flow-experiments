REL="$REL"
VAR="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline"
#POP="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/experiments/stochastic_seed0"
#GOLD="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/reads_mason2/NA12878_1-h37maj-mapql10.sam"
#EXP="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/experiments/exp_seed0"
GOLD=$1
POP=$2

if [ "$#" -ne 2 ]
then
    echo "ERROR: incorrect number of arguments (should be 2)"
    echo "Usage:"
    echo "$(basename $0) GOLD POP_VAR_DIR"
    exit
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

array=( "EUR" "AMR" "SAS" "AFR" "EAS" )
WC=`wc -l merging.ids | cut -d ' ' -f1`
ARR=(`tail -$(($WC+1)) tmp_merge.log | head -$(($WC))`)
for i in $(seq 0 4)
do
    TMP=`python -O $REL/scripts/analyze_diploid_indels.py -c 21 -g $GOLD -p 0 -vr $VAR/NA12878.var -vs $POP/*${array[$i]}*.var -n ${ARR[$i]} | grep 'sensitivity_all'`
    # TMP=`python -O $REL/scripts/analyze_diploid_indels.py -c 21 -g $GOLD -p 0 -vr $VAR/NA12878.var -vs $POP/21_superpop_"${array[$i]}"_thrsd0_stochastic.var -n ${ARR[$i]} | grep 'sensitivity_all'`
    update
done

echo "RECALL = `bc -l <<< "$TP / $SUM"`"
