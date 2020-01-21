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
OUTLIST="$2"
DASHLIST=${STORAGE}"/diploid_dashing/dashing_list_${DATE}.txt"
D_DASH_OUT=${STORAGE}"/diploid_dashing/${DATE}/"

mkdir -p $D_DASH_OUT

K=31 #kmer size
P=10 #threads
S=10 #sketch size

# For 1k Genomes from update_genomes.py / VCF files:
#SUFFIX=".fa"

# Union paired end sketches to produce joint sketch
SUFFIX="hapA_indel-1M.hll"
for FILE1 in $(find ${DIR} -name "*${SUFFIX}"); do
	FNAME=${FILE1%${SUFFIX}}
	FILE2="${FNAME}hapB_indel-1M.hll"
	$DASHING union $FILE1 $FILE2 -o "${D_DASH_OUT}/${FNAME##*/}.hll" 
done

SUFFIX=".hll"
#Make outlist empty file
echo "" > $OUTLIST

for FILE1 in $(find $D_DASH_OUT -name "*${SUFFIX}"); do
	FNAME=${FILE1%${SUFFIX}}
	echo "$FNAME $($DASHING view $FILE1)" >> $OUTLIST
done
