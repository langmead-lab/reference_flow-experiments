python ../scripts/build_erg.py -r genome.fa -m gold_hapA.fa -a gold_hapB.fa -v gold.var --print-mg 1 > test_NA12878_hapAB.fa

if cmp -s "gold_hapAB.fa" "test_NA12878_hapAB.fa"; then
    echo "hapAB.fa -- PASS"
else
    echo "hapAB.fa -- FAILED"
    exit
fi
# grep -v '##' test_NA12878.vcf | cut -f 1-7 > test_NA12878_v.vcf
# if cmp -s "gold.vcf" "test_NA12878_v.vcf"; then
#     echo ".vcf -- PASS"
# else
#     echo ".vcf -- FAILED"
#     exit
# fi

echo "------------------------------"
echo "$0 has PASSED!"
#ls test_NA12878*
read -p "Delete generated test files?
--
`ls test_NA12878*`
--
[y/N] " if_del
case $if_del in
    [Yy]* ) rm test_NA12878*; break;;
    * ) exit;;
esac
