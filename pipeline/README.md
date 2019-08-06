# Examples of how to run the pipeline #

## Build personalized genomes and simulate reads from them, build major allele reference ##

Take chromosome 21 and individual NA12878 as an example:

`prepare_genomes.sh -c 21 -f chr21.fa -s NA12878 -v 21.vcf`

Outputs:

- `NA12878_hapA.fa`, `NA12878_hapB.fa`: personalized genomes for each NA12878 haplotype. Both SNPs and INDELs are considered.
- `NA12878.var`: recording variants and corresponding offsets (wrt to chr21.fa) used to build each haplotype. This is used later for calculating alignment accuracy.
- `NA12878.vcf`: VCF for included variants.
- `21_h37maj.fa`: major allele reference.
- `21_h37maj.vcf`, `21_h37maj.var`: VAR and VCF for major allele reference.
- `indexes/*`: bowtie2 indexes for 21_h37maj

## Build population genomes ##

`build_pop_vcf.sh -r -c 21 -C superpop -f chr21.fa -v 21.vcf`

Build stochastic population genomes for each superpopulation (AMR, EUR, EAS, SAS, AFR) in 1000 Genomes Project.
Remove `-r` to build major population genomes. Use `buiild_pop_vcf.sh -h` to see usage.

Outputs:

- `21_superpop_*_thrsd0_stochastic_b1.fa/var/vcf`: FASTA, VAR and VCF for stochastic population genomes.
The values after _thrds__ and _b__ depend on `-t` and `-b` options.
