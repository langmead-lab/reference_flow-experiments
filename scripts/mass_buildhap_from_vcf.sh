STORAGE=/net/langmead-bigmem-ib.bluecrab.cluster/storage
SCRATCH=/home-1/bsolomo9\@jhu.edu/scratch/bsolomo9/genome_relaxation/
#RELAX=/scratch/groups/blangme2/naechyun/relaxing/
RELAX=/home-1/bsolomo9\@jhu.edu/genome_relaxation/
DATE=$(date +"%m-%d-%Y")


# Triple user input instead of hard-coded 
VCF_FILE="$1" #$STORAGE/naechyun/1000Genomes/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
CHROM_NUM="$2" #9
RANDSIZE="$3" #100
RANDSET=${SCRATCH}/1000G_R${RANDSIZE}_${DATE}.txt
RANDSTORAGE=${STORAGE}/bsolomo9/1000G_R${RANDSIZE}_${DATE}

# RANDSTORAGE automatically pulls down dates. Don't run the same randsize on the same day.
mkdir -p $RANDSTORAGE

# Run script to select a random set
touch $RANDSET 
python getRand.py $SCRATCH/phase3_names.txt ${RANDSIZE} > $RANDSET

# Run script to generate updated reference genome with VCF and produce -hapA, -hapB, and .var files prefixed by out-prefix

while read -r NAME
do
	python $RELAX/scripts/update_genome.py --ref $STORAGE/indexes/hs37d5.fa --vcf $VCF_FILE --chrom $CHROM_NUM --out-prefix ${RANDSTORAGE}/$NAME --name $NAME --include-indels 0
done < $RANDSET

#/scratch/groups/blangme2/naechyun/software/mason-0.1.2-Linux-x86_64/bin/mason \
#    illumina $RELAX/na12878/indels/na12878-chr9-indel-hapB.fa \
#    -N 1000000 -sq -n 100 -hs 0 -hi 0 -hn 1 -mp \
#    -o na12878-chr9-phase3_hapB_indel-1M.fq


