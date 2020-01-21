STORAGE=/net/langmead-bigmem-ib.bluecrab.cluster/storage/bsolomo9
SCRATCH=/home-1/bsolomo9\@jhu.edu/scratch/bsolomo9/genome_relaxation/
#RELAX=/scratch/groups/blangme2/naechyun/relaxing/
RELAX=/home-1/bsolomo9\@jhu.edu/genome_relaxation/
DATE=$(date +"%m-%d-%Y")
#RANDSIZE="$2"
#RANDSET=${SCRATCH}/1000G_R${RANDSIZE}_${DATE}.txt
#RANDSTORAGE=${STORAGE}/bsolomo9/1000G_R${RANDSIZE}_${DATE}
DASHING=/home-1/bsolomo9@jhu.edu/software/dashing/dashing
#mkdir -p $RANDSTORAGE
DIR="$1"
DASHLIST=${STORAGE}/"dashing_list_${DATE}.txt"

K=31 #kmer size
P=10 #threads
S=10 #sketch size

# For 1k Genomes from update_genomes.py / VCF files:
#SUFFIX=".fa"

# Build list of all fq in directory
find $DIR -name "*.fq" > $DASHLIST

# Build a dash sketch for each fq in directory
$DASHING sketch -k${K} -p${P} -S${S} -F $DASHLIST

# Union paired end sketches to produce joint sketch
SUFFIX="_indel-1M_1.fq.w.${K}.spacing.${S}.hll"

for FILE1 in $(find $DIR -name "*${SUFFIX}"); do
	FNAME=${FILE1%_1.fq.w.${K}.spacing.${S}.hll}
	FILE2="${FILE1%_1.fq.w.${K}.spacing.${S}.hll}_2.fq.w.${K}.spacing.${S}.hll"
	$DASHING union $FILE1 $FILE2 -o "${FNAME}.hll"
done

