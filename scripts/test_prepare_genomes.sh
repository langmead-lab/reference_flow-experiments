### This script is to test if the results of prepare_genomes.sh match that using bcftools.
### Currently only support chr21 for NA12878
PREP_PATH="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/test_prep"
INDIV="NA12878"
CHROM="21"
# INDIV=$1
# CHROM=$2

if [ "0" -eq "$(diff ${INDIV}_hapA.fa ${PREP_PATH}/bcftools_1.9_206-${INDIV}_1.fa | wc -l)" ];
then
    echo "test -hapA.fa: pass"
else
    echo "test -hapA.fa: failed"
    exit
fi

if [ 0 -eq "$(diff ${INDIV}_hapB.fa ${PREP_PATH}/bcftools_1.9_206-${INDIV}_2.fa | wc -l)" ];
then
    echo "test -hapB.fa: pass"
else
    echo "test -hapB.fa: failed"
    exit
fi

if [ 0 -eq "$(diff ${INDIV}.var ${PREP_PATH}/${INDIV}.var | wc -l)" ];
then
    echo "test personalized.var: pass"
else
    echo "test personalized.var: failed"
    exit
fi

if [ 0 -eq "$(diff ${CHROM}_h37maj.fa ${PREP_PATH}/${CHROM}_h37maj.fa | wc -l)" ];
then
    echo "test maj_allele_ref.fa: pass"
else
    echo "test maj_allele_ref.fa: failed"
    exit
fi

if [ 0 -eq "$(diff ${CHROM}_h37maj.var ${PREP_PATH}/${CHROM}_h37maj.var | wc -l)" ];
then
    echo "test maj_allele_ref.var: pass"
else
    echo "test maj_allele_ref.var: failed"
    exit
fi
