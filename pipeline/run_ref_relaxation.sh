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

# If want dual fastqs cat them or write a separate script. This is the 'fast' version.
# Mostly taken from /net/langmead-bigmem-ib.bluecrab.cluster/storage/bsolomo9/1000G_R100_03-14-2019
FASTQ="$1"

# /net/langmead-bigmem-ib.bluecrab.cluster/storage/bsolomo9/relaxation/1000G_R100_03-14-2019_out/
OUTDIR="$2"
THRESH="$3"

FNAME=$(basename $FASTQ)
OUTPREFIX="${OUTDIR}/${FNAME%.*}"

OUTSAM=${OUTPREFIX}-h37maj.sam
# Run the first pass alignment
if [ -f $OUTSAM ]; then
	echo "$OUTSAM already exists"
else
	echo "First pass alignment of $FNAME"
	bowtie2 -x ${INDEX}/21_h37maj -S $OUTSAM -U $FASTQ
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
PASS2_OUTSAM="${OUTSAM%.*}_pass2-t${THRESH}-${s}s${FRAC}_b${BSIZE}.sam"

#PASS2_OUTSAM=${OUTSAM%.*}-2pass_BSIZE${BSIZE}.sam
LOG_FILE=${PASS2_OUTSAM%.*}.log
BLOCK_INDEX=/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/pipeline


if [ -f $PASS2_OUTSAM ]; then
	echo "$PASS2_OUTSAM already exists"
else
	echo "Running relaxation alignment..."
	rm merging.*
	rm tmp_merge.log
	$TIME -v $RELAX/pipeline/global_run_2ndpass.sh -b ${BSIZE} -I $BLOCK_INDEX/stochastic_pop_genomes/block_${BSIZE}/indexes -i $FILTER_FQ -p "${OUTSAM%.*}_pass2-t${THRESH}" -S $RELAX/scripts/ 2> $LOG_FILE 
fi
# Calculate accuracy
GOLD_SAM=${FASTQ%_*}.sam
VAR_FILE=${FASTQ%_*}.var
VAR_GENOMES=$BLOCK_INDEX/stochastic_pop_genomes/block_${BSIZE}
ACC_OUT=${LOG_FILE%.*}_accuracy.txt

echo $GOLD_SAM
echo $VAR_FILE

#VAR_GENOMES=$(dirname $PASS2_OUTSAM)
echo "Calculating final accuracy"
$RELAX/pipeline/calc_2ndpass_acc.sh -c 21 -g $GOLD_SAM -V $VAR_FILE -P $VAR_GENOMES -S $RELAX/scripts/ > $ACC_OUT 
