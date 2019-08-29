## Reference Bias Related Scripts
### VCF Preparation

Extracts HET variants carried by an 1000 Genomes Project individual.
Only bi-allelic heterozygous sites (HETs) are considered.

`bcftools view -s NA12878 21.vcf | bcftools view -i "AC>0" -v snps -g het -m2 -M2 > 21_NA12878.vcf`

Here we use bcftools 1.9-206-g4694164 and htslib 1.9-258-ga428aa2

### Calculated Reference Bias

`python lift_ref_flow.py -v 21_NA12878.vcf -s sam.list -n name.list -f chr21.fa -o bias.txt`

* `sam.list`: paths to the SAM files of interest, each line has one path.

* `name.list`: names for the SAM files. `lift_ref_flow.py` also reports separate reference bias logs for each SAM file, named after the prefixes specified in `name.list`

### Extract Reads Aligned to Biased Sites

`python find_reads_given_HET.py -s sam.list -v 21_NA12878.vcf -f bias.txt -o name.list -m above080_or_below020.reads -r 0-0.2,0.8-1`

List all the reads covering HET sites with base bias <= 0.2 or base bias >= 0.8.

`python find_reads_given_HET.py -s sam.list -v 21_NA12878.vcf -f bias.txt -o name.list -m above045_and_below055.reads -r 0.45-0.55 --sample 0.01`

List all the reads covering HET sites with base bias <= 0.55 or base bias >= 0.45. Since there can be many reads passing this filter, setting `--sample` to 0.01 randomly selects 1% of the reads.
