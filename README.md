# Reference Flow
The reference flow method utilizes population genomic data to enhance alignment accuracy and reduce reference bias with high computational efficiency.
Based on selection of ambiguous alignments and re-alignment to a set of genomes covering a wider variant space, reference flow aligns reads as accurately as the "personalized" approach with high computational efficiency.


## Usage of utility scripts for simulated data
### analyze_diploid_indels.py
Reads a sam file and analyzes its accuracy.

Example:

`python -O analyze_diploid_indels.py 
--chrm 21 
--sam chr21-10M-h37maj-mapql10.sam
--golden chr21-h37maj-lowq.fq.sam 
--var_reads chr21_na12878.var 
--var_sample chr21_h37maj.var 
--personalized 0`

### calculate_as_from_simulation.py
Compares the alignment score (_AS:i_) with simulation. Also reports number of overlapping simulated sequencing errors and alignment correctness.

`python calculate_as_from_simulation.py 
--fn_sam chr21-10M-h37maj-mapql10.sam 
--fn_sam_stats chr21-10M-h37maj-mapql10.sam-stats.pkl 
--fn_gold chr21-h37maj-lowq.fq.sam`
