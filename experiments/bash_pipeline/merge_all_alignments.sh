STORAGE=/net/langmead-bigmem-ib.bluecrab.cluster/storage
SCRATCH=/home-1/bsolomo9\@jhu.edu/scratch/bsolomo9/genome_relaxation/
#RELAX=/scratch/groups/blangme2/naechyun/relaxing/
RELAX=/home-1/bsolomo9\@jhu.edu/genome_relaxation/

INDEX=/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/indexes/

#module load samtools
#module load bowtie2
#module load python/3.7-anaconda

# NOT INCLUDED IN THIS SCRIPT: GENERATING THE hapA.fa / hapB.fa / .var files using update_genome.py
# the fast way to do that is through the mass_buildhap_from_vcf.sh file

#ALSO FORGOT THE MASON SIM

# If want dual fastqs cat them or write a separate script. This is the 'fast' version.
# Mostly taken from /net/langmead-bigmem-ib.bluecrab.cluster/storage/bsolomo9/1000G_R100_03-14-2019
FASTQ="$1"

# /net/langmead-bigmem-ib.bluecrab.cluster/storage/bsolomo9/relaxation/1000G_R100_03-14-2019_out/
OUTDIR="$2"
THRESH="$3"
CHROM=21
CHROM_VAR="/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/refflow/chr21/pipeline_100/major/21_h37maj.var"

FNAME=$(basename $FASTQ)
OUTPREFIX="${OUTDIR}/${FNAME%.*}"

OUTSAM=${OUTPREFIX}-h37maj.sam
PASS1_LOG=${OUTPREFIX}-h37maj.log
# Run the first pass alignment
if [ -f $OUTSAM ]; then
	echo "$OUTSAM already exists"
else
	echo "First pass alignment of $FNAME"
	$TIME -v bowtie2 -x ${INDEX}/21_h37maj -S $OUTSAM -U $FASTQ 2> $PASS1_LOG 
	echo "Alignment written to ${OUTSAM}"
fi

#Get second pass reads
FILTER_FQ="${OUTSAM%.*}-mapql${THRESH}.fq"

if [ -f $FILTER_FQ ]; then
	echo "$FILTER_FQ already exists"
else
	echo "Splitting $(basename $OUTSAM)"
	./get_lowq_by_mapq.sh $OUTSAM $THRESH
	echo "2nd Pass reads written to $FILTER_FQ"
fi


# TEMP FILES MEAN MOVE OUT OF HOME DIRECTORY (OR WHEREVER THIS SCRIPT IS) AND INTO OUTPUT DIRECTORY
echo "Moving to $OUTDIR"
cd $OUTDIR

# Run second pass alignment
BSIZE=1
s=1
FRAC=0

#THIS NEEDS TO BE ADJUSTED TO LOOK FOR THE RIGHT THING
PASS2_OUTSAM="${OUTSAM%.*}_pass2-t${THRESH}-${s}s${FRAC}_b${BSIZE}.sam"

#PASS2_OUTSAM=${OUTSAM%.*}-2pass_BSIZE${BSIZE}.sam
LOG_FILE=${PASS2_OUTSAM%.*}.log
BLOCK_INDEX=/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline


# THIS WILL ALWAYS FAIL BECAUSE THATS NOT THE RIGHT OUTPUT SAM
if [ -f $PASS2_OUTSAM ]; then
	echo "$PASS2_OUTSAM already exists"
else
	echo "Running relaxation alignment... ($PASS2_OUTSAM)"
	echo "TEMP DISABLED FOR TESTING SHOULD NOT SEE THIS"
	#rm merging.*
	#rm tmp_merge.log
	#$TIME -v $RELAX/pipeline/global_run_2ndpass.sh -b ${BSIZE} -I $BLOCK_INDEX/stochastic_pop_genomes/block_${BSIZE}/indexes -i $FILTER_FQ -p "${OUTSAM%.*}_pass2-t${THRESH}" -S $RELAX/scripts/ 2> $LOG_FILE 
fi
# Calculate accuracy
GOLD_SAM=${FASTQ%_*}.sam
VAR_FILE=${FASTQ%_*}.var
VAR_GENOMES=$BLOCK_INDEX/stochastic_pop_genomes/block_${BSIZE}
ACC_OUT=${LOG_FILE%.*}_accuracy.txt

echo $GOLD_SAM
echo $VAR_FILE

#VAR_GENOMES=$(dirname $PASS2_OUTSAM)
echo "Calculating 2nd-pass accuracy"
echo "TEMP DISABLED FOR TESTING"
#$RELAX/pipeline/calc_2ndpass_acc.sh -c $CHROM -g $GOLD_SAM -V $VAR_FILE -P $VAR_GENOMES -S $RELAX/scripts/ > $ACC_OUT 


echo "Calculating pre-merge accuracies"
ACC_1P=${OUTSAM%.*}_premerge_acc.txt
echo "Calculating 1st-pass accuracy"
echo "TEMP DISABLED FOR TESTING"
#python3.7 $RELAX/scripts/analyze_diploid_indels.py -c $CHROM -n $OUTSAM -vr $VAR_FILE -vs $CHROM_VAR -p 0 -g $GOLD_SAM > $ACC_1P

echo "Calculating merged 2nd-pass accuracy"
#LIST_2PASS=$(find ${OUTDIR} -name "${FNAME%.*}-h37maj_pass2-t${THRESH}*\+*.sam")

AFR_SAM=$(find ${OUTDIR} -name "${FNAME%.*}-h37maj_pass2-t${THRESH}*AFR*\+*.sam")
AMR_SAM=$(find ${OUTDIR} -name "${FNAME%.*}-h37maj_pass2-t${THRESH}*AMR*\+*.sam")
EAS_SAM=$(find ${OUTDIR} -name "${FNAME%.*}-h37maj_pass2-t${THRESH}*EAS*\+*.sam")
EUR_SAM=$(find ${OUTDIR} -name "${FNAME%.*}-h37maj_pass2-t${THRESH}*EUR*\+*.sam")
SAS_SAM=$(find ${OUTDIR} -name "${FNAME%.*}-h37maj_pass2-t${THRESH}*SAS*\+*.sam")

AFR_VAR="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/stochastic_pop_genomes/block_1/21_superpop_AFR_thrsd0_stochastic_b1.var"
AMR_VAR="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/stochastic_pop_genomes/block_1/21_superpop_AMR_thrsd0_stochastic_b1.var"
EAS_VAR="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/stochastic_pop_genomes/block_1/21_superpop_EAS_thrsd0_stochastic_b1.var"
EUR_VAR="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/stochastic_pop_genomes/block_1/21_superpop_EUR_thrsd0_stochastic_b1.var"
SAS_VAR="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline/stochastic_pop_genomes/block_1/21_superpop_SAS_thrsd0_stochastic_b1.var"


for TSAM in $AFR_SAM; do
ACC_2P=${TSAM%.*}_premerge_acc.txt
python3.7 $RELAX/scripts/analyze_diploid_indels.py -c $CHROM -n $TSAM -vr $VAR_FILE -vs $AFR_VAR -p 0 -g $GOLD_SAM > $ACC_2P
done

for TSAM in $AMR_SAM; do
ACC_2P=${TSAM%.*}_premerge_acc.txt
python3.7 $RELAX/scripts/analyze_diploid_indels.py -c $CHROM -n $TSAM -vr $VAR_FILE -vs $AMR_VAR -p 0 -g $GOLD_SAM > $ACC_2P
done

for TSAM in $EAS_SAM; do
ACC_2P=${TSAM%.*}_premerge_acc.txt
python3.7 $RELAX/scripts/analyze_diploid_indels.py -c $CHROM -n $TSAM -vr $VAR_FILE -vs $EAS_VAR -p 0 -g $GOLD_SAM > $ACC_2P
done

for TSAM in $EUR_SAM; do
ACC_2P=${TSAM%.*}_premerge_acc.txt
python3.7 $RELAX/scripts/analyze_diploid_indels.py -c $CHROM -n $TSAM -vr $VAR_FILE -vs $EUR_VAR -p 0 -g $GOLD_SAM > $ACC_2P
done

for TSAM in $SAS_SAM; do
ACC_2P=${TSAM%.*}_premerge_acc.txt
python3.7 $RELAX/scripts/analyze_diploid_indels.py -c $CHROM -n $TSAM -vr $VAR_FILE -vs $SAS_VAR -p 0 -g $GOLD_SAM > $ACC_2P
done


#for TSAM in $(find ${OUTDIR} -name "${FNAME%.*}-h37maj_pass2-t${THRESH}*\+*.sam"); do
#	ACC_2P=${TSAM%.*}_premerge_acc.txt
#	python3.7 $RELAX/scripts/analyze_diploid_indels.py -c $CHROM -n $TSAM -vr $VAR_FILE -vs $CHROM_VAR -p 0 -g $GOLD_SAM > $ACC_2P
#done
