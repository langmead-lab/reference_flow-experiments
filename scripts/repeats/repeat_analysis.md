#### Download RepeatMasker annotation
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
gunzip rmsk.txt.gz
```

The annotation is provided by UCSC, from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/.

#### Convert annotation to BED format
This can be done with script `extract_repeats_from_rmsk.py`. 
This script supports querying different repeat labels in `rmsk.txt`, including repName, repClass and repFamily.

An example to extract the L1 repeat family in BED file `rmsk-L1.bed`:
```
python extract_repeats_from_rmsk.py -f rmsk.txt -l repFamily -r L1 -o rmsk-L1.bed
```

#### Convert an allelic bias file to BED format
This can be done with `extract_biased_sites.py`
`extract_biased_sites.py` and `extract_repeats_from_rmsk.py` are wrapped in `run_extract_biased_sites.sh`.
Make sure to change paths to your local directory in the shell script beforing executing.

#### Intersecting bias BED files with repeat annotations
This is wrapped in shell script `intersect_with_repeats.sh`
