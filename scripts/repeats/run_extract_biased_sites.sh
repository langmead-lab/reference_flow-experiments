DIR="../../snakemake/SRR622457_se/experiments/NA12878"
DIR_RPT="../../snakemake/SRR622457_se/experiments/NA12878/repeat_analysis"
python extract_biased_sites.py -f ${DIR}/wg-GRC-NA12878.allele.bias -o ${DIR_RPT}/GRC-NA12878-biased_sites.bed
python extract_biased_sites.py -f ${DIR}/wg-major-liftover-NA12878.allele.bias -o ${DIR_RPT}/major-NA12878-biased_sites.bed
python extract_biased_sites.py -f ${DIR}/wg-refflow-10-thrds0_S1_b1000_ld1-NA12878.allele.bias -o ${DIR_RPT}/randflow_ld-NA12878-biased_sites.bed
python extract_biased_sites.py -f ${DIR}/wg-per-merged-liftover-NA12878.allele.bias -o ${DIR_RPT}/per-NA12878-biased_sites.bed

python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repFamily -r L1 -o ${DIR_RPT}/rmsk-L1.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repFamily -r L2 -o ${DIR_RPT}/rmsk-L2.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repFamily -r Alu -o ${DIR_RPT}/rmsk-Alu.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repFamily -r ERVL -o ${DIR_RPT}/rmsk-ERVL.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repFamily -r ERV1 -o ${DIR_RPT}/rmsk-ERV1.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repFamily -r MIR -o ${DIR_RPT}/rmsk-MIR.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repClass -r Satellite -o ${DIR_RPT}/rmsk-Satellite.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repClass -r SINE -o ${DIR_RPT}/rmsk-SINE.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repClass -r LINE -o ${DIR_RPT}/rmsk-LINE.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repClass -r LTR -o ${DIR_RPT}/rmsk-LTR.bed
python extract_repeats_from_rmsk.py -f ${DIR_RPT}/rmsk.txt -l repClass -r Simple_repeat -o ${DIR_RPT}/rmsk-Simple_repeat.bed
