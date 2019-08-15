## Reference Bias Related Scripts
### VCF Preparation

Extracts HET variants carried by an 1000 Genomes Project individual.
Only bi-allelic heterozygous sites (HETs) are considered.

`bcftools view -s NA12878 21.vcf | bcftools view -i "AC>0" -v snps -g het -m2 -M2 > 21_NA12878.vcf`

Here we use bcftools 1.9-206-g4694164 and htslib 1.9-258-ga428aa2

### Calculated Reference Bias

`python lift_ref_flow.py -v 21.vcf -s list.sams -n list.names -f chr21.fa -o output.txt`

* `list.sams`: paths to the SAM files of interest, each line has one path.

* `list.names`: names for the SAM files. `lift_ref_flow.py` also reports separate reference bias logs for each SAM file, named after the prefixes specified in `list.names`
