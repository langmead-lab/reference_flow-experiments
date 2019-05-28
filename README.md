# Reference Flow
The reference flow method utilizes population genomic data to enhance alignment accuracy and reduce reference bias with high computational efficiency.
Based on selection of ambiguous alignments and re-alignment to a set of genomes covering a wider variant space, reference flow aligns reads as accurately as the "personalized" approach with high computational efficiency.


## Usage of utility scripts for simulated data
### Analyze the accuracy of alignments
#### analyze_diploid_indels.py 
Reads a sam file and analyzes its accuracy (sensitivity and precision).
Diploid genomes including indels are supported.
Set `-O` to avoid printing debug messages.

Important arguments:

Command                       | Content
----------------------------- | ----------------
-n SAM, --sam SAM             | target .sam file
-g GOLDEN, --golden GOLDEN    | golden .sam file
-vs VAR_SAMPLE, --var_sample VAR_SAMPLE | .var file for the personlized genome
-p PERSONALIZED, --personalized PERSONALIZED (0) | set 2 for diploid genomes, set 0 (default) for haploid genomes
-c CHRM, --chrm CHRM          | target chromosome name

Example:

`python -O analyze_diploid_indels.py 
--chrm 21 
--sam chr21-10M-h37maj-mapql10.sam
--golden chr21-h37maj-lowq.fq.sam 
--var_reads chr21_na12878.var 
--var_sample chr21_h37maj.var 
--personalized 0`

### Look into multi-mapped alignments
#### analyze_multimappers.py
Reads a top-q target alignment file (can be generated with `-k q` input of Bowtie2) and reports properties of interest.

Example:

`python analyze_multimappers.py 
-stats HG01047-should-be-correct-k3.sam-stats.pkl 
-sam HG01047-should-be-correct-k3.sam 
-vr chr21_na12878.var 
-vs HG01047.var 
-g chr21-per-incorrect.fq.sam 
-o analyze_multimappers_HG01047.log`

### Look into alignments and group by different characteristics
#### extract_unique_var.py
Reads two .var files and outputs unique variants w.r.t each genome.
This step is important for following `calculate_as_from_simulation.py` script since a variant existing on both simulation and target genome does not make a difference.

Example:

`python extract_unique_var.py
--fn_var1 chr21_na12878.var 
--fn_var2 chr21_h37maj.var 
--unique 1
`

#### calculate_as_from_simulation.py
Supported characteristics:
- Whether aligned score (_AS:i_) is higher than simulation.
- Number of sequencing errors added during simulation
- Number of unique overlapping variants on simulated genome and target genome
- Alignment correctness

`python calculate_as_from_simulation.py 
--fn_sam chr21-10M-h37maj-mapql10.sam 
--fn_sam_stats chr21-10M-h37maj-mapql10.sam-stats.pkl 
--fn_gold chr21-h37maj-lowq.fq.sam`
