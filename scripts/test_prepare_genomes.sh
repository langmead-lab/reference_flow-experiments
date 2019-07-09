PREP_PATH="/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/relaxation/chr21/test_prep"
INDIV=$1
CHROM=$2

if [ "0" -eq "$(diff ${INDIV}_hapA.fa ${PREP_PATH}/${INDIV}_hapA.fa | wc -l)" ];
then
    echo "test -hapA.fa: pass"
else
    echo "test -hapA.fa: failed"
    exit
fi

if [ 0 -eq "$(diff ${INDIV}_hapB.fa ${PREP_PATH}/${INDIV}_hapB.fa | wc -l)" ];
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

if [ 0 -eq "$(diff ${CHROM}_h37maj_with_indels.fa ${PREP_PATH}/${CHROM}_h37maj_with_indels.fa | wc -l)" ];
then
    echo "test maj_allele_ref.fa: pass"
else
    echo "test maj_allele_ref.fa: failed"
    exit
fi

if [ 0 -eq "$(diff ${CHROM}_h37maj_with_indels.var ${PREP_PATH}/${CHROM}_h37maj_with_indels.var | wc -l)" ];
then
    echo "test maj_allele_ref.var: pass"
else
    echo "test maj_allele_ref.var: failed"
    exit
fi
