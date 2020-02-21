# Reference Flow
The reference flow method utilizes population genomic data to enhance alignment accuracy and reduce reference bias with high computational efficiency.
We first aligned reads to a major-allele reference genome based on the 1000 Genomes Project GRCh38 call set.
We assigned unmapped and ambiguous reads determined by mapping quality threshold to the "deferred" group.
These reads are re-aligned using a set of population genomes based on "superpopulation" labels in the 1000 Genomes Project.
We finally merged all the reads into an unified SAM output, which is based on the coordinate system of GRCh38.


## Snakemake
The reference flow method is built based on Snakemake for efficient and scalable computing.
The GRCh37 and GRCh38 pipelines are put under `grch37` and `grch38` directories, respectively.
Users may modify the `*.yaml` files based on the environment and run `snakemake -np` to see if configurations are set correctly.

Finally, run

```
snakemake -j 32
```

to start the pipeline. Option `-j` specifies the number of threads used. In this example, `32` threads are used.

