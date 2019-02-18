STORAGE=/net/langmead-bigmem-ib.bluecrab.cluster/storage
SCRATCH=/home-1/bsolomo9\@jhu.edu/scratch/bsolomo9/genome_relaxation/
#RELAX=/scratch/groups/blangme2/naechyun/relaxing/
RELAX=/home-1/bsolomo9\@jhu.edu/genome_relaxation/
#DATE="$1"
#$(date +"%m-%d-%Y")
#RANDSIZE="$2"
#RANDSET=${SCRATCH}/1000G_R${RANDSIZE}_${DATE}.txt
#RANDSTORAGE=${STORAGE}/bsolomo9/1000G_R${RANDSIZE}_${DATE}
MASH=/home-1/bsolomo9@jhu.edu/software/mash-Linux64-v2.1/mash
#mkdir -p $RANDSTORAGE
DIR="$1"
SUFFIX="_indel-1M_1.fq"

# Run script to select a random set
#touch $RANDSET 
#python getRand.py $SCRATCH/phase3_names.txt ${RANDSIZE} > $RANDSET

# Run script to generate updated reference genome with VCF and produce -hapA, -hapB, and .var files prefixed by out-prefix
#	python $RELAX/scripts/update_genome.py --ref $STORAGE/indexes/hs37d5.fa --vcf $STORAGE/naechyun/1000Genomes/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --chrom 9 --out-prefix ${RANDSTORAGE}/$NAME --name $NAME --include-indels 1


# Run script to take hap files and simulate 1M reads in fastq format
for FILE1 in $(find $DIR -name "*${SUFFIX}"); do
	FILE2="${FILE1%_1.fq}_2.fq"
	#echo $FILE1
	#echo $FILE2
	cat $FILE1 $FILE2 | $MASH sketch - -o "${FILE1%_1.fq}.msh"
done

