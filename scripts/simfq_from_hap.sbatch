#!/usr/bin/bash

#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node 2

STORAGE=/net/langmead-bigmem-ib.bluecrab.cluster/storage
SCRATCH=/home-1/bsolomo9\@jhu.edu/scratch/bsolomo9/genome_relaxation/
#RELAX=/scratch/groups/blangme2/naechyun/relaxing/
RELAX=/home-1/bsolomo9\@jhu.edu/genome_relaxation/
DATE="$1"
#$(date +"%m-%d-%Y")
RANDSIZE="$2"
NAME="$3"
#RANDSET=${SCRATCH}/1000G_R${RANDSIZE}_${DATE}.txt
RANDSTORAGE=${STORAGE}/bsolomo9/1000G_R${RANDSIZE}_${DATE}
TIME=/home-1/bsolomo9@jhu.edu/bin/time
# Run script to take hap files and simulate 1M reads in fastq format
	echo "Constructing Hap A for $NAME"
	($TIME -v /scratch/groups/blangme2/naechyun/software/mason-0.1.2-Linux-x86_64/bin/mason \
    	illumina ${RANDSTORAGE}/${NAME}_hapA.fa \
   		-N 1000000 -sq -n 100 -hs 0 -hi 0 -hn 1 -mp \
    	-o ${RANDSTORAGE}/${NAME}_hapA_indel-1M.fq) > ${RANDSTORAGE}/${NAME}_hapA.log 2> ${RANDSTORAGE}/${NAME}_hapA_time.log

	gzip ${RANDSTORAGE}/${NAME}_hapA_indel-1M.fq

	echo "Constructing Hap B for $NAME"
	($TIME -v /scratch/groups/blangme2/naechyun/software/mason-0.1.2-Linux-x86_64/bin/mason \
    	illumina ${RANDSTORAGE}/${NAME}_hapB.fa \
    	-N 1000000 -sq -n 100 -hs 0 -hi 0 -hn 1 -mp \
    	-o ${RANDSTORAGE}/${NAME}_hapB_indel-1M.fq) > ${RANDSTORAGE}/${NAME}_hapB.log 2> ${RANDSTORAGE}/${NAME}_hapB_time.log

	gzip ${RANDSTORAGE}/${NAME}_hapB_indel-1M.fq
