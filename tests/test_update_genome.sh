python ../scripts/update_genome.py -r genome.fa -v variants.vcf -c 21 -op test_NA12878 -i -s NA12878

if cmp -s "gold_hapA.fa" "test_NA12878_hapA.fa"; then
    echo "hapA.fa -- PASS"
else
    echo "hapA.fa -- FAILED"
    exit
fi
if cmp -s "gold_hapB.fa" "test_NA12878_hapB.fa"; then
    echo "hapB.fa -- PASS"
else
    echo "hapB.fa -- FAILED"
    exit
fi
if cmp -s "gold.var" "test_NA12878.var"; then
    echo ".var -- PASS"
else
    echo ".var -- FAILED"
    exit
fi
grep -v '##' test_NA12878.vcf | cut -f 1-7 > test_NA12878_v.vcf
if cmp -s "gold.vcf" "test_NA12878_v.vcf"; then
    echo ".vcf -- PASS"
else
    echo ".vcf -- FAILED"
    exit
fi

echo "------------------------------"
echo "test_update_genome has PASSED!"
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
