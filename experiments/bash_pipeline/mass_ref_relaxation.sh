STORAGE=/net/langmead-bigmem-ib.bluecrab.cluster/storage
SCRATCH=/home-1/bsolomo9\@jhu.edu/scratch/bsolomo9/genome_relaxation/
#RELAX=/scratch/groups/blangme2/naechyun/relaxing/
RELAX=/home-1/bsolomo9\@jhu.edu/genome_relaxation/
#DATE="$1"
#$(date +"%m-%d-%Y")

DIR="$1"
OUTDIR="$2"
THRESH="$3"
#RANDSET=${SCRATCH}/1000G_R${RANDSIZE}_${DATE}.txt
#RANDSTORAGE=${STORAGE}/bsolomo9/1000G_R${RANDSIZE}_${DATE}

#mkdir -p $RANDSTORAGE

# Run script to select a random set
#touch $RANDSET 
#python getRand.py $SCRATCH/phase3_names.txt ${RANDSIZE} > $RANDSET

# Run script to generate updated reference genome with VCF and produce -hapA, -hapB, and .var files prefixed by out-prefix
#	python $RELAX/scripts/update_genome.py --ref $STORAGE/indexes/hs37d5.fa --vcf $STORAGE/naechyun/1000Genomes/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --chrom 9 --out-prefix ${RANDSTORAGE}/$NAME --name $NAME --include-indels

# Run script to take hap files and simulate 1M reads in fastq format


# THIS ONLY WORKS IF OYU MOVE ALL THE HAPA/HAPB _1/_2.fq TO A DIFFERENT DIRECTORY
# (run mv *_hapA*.fq temp_hap_stuff/ )
# mv *_hapB*.fq temp_hap_stuff/

for FILE in $(find $DIR -maxdepth 1 -name "*_1.fq")
do
	echo "Running ref relaxation for $FILE"
	${RELAX}/pipeline/run_ref_relaxation.sh $FILE ${OUTDIR} $THRESH
done
