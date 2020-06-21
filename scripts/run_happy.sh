# Compare the accuracy of a VCF file with the truth variant set using hap.py.
# This needs to be run in a python2 environment.

# GIAB-GRCh38-HG001/NA12878 call set.
TRUE_VCF=/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/bias_inspector/variant_analysis/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
CONF_REGION=/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/bias_inspector/variant_analysis/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed
# NA12878 variants from the 1KG GRCh38 calls set.
# TRUE_VCF=/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/refflow-exp/snakemake/real_2/personalized/NA12878/wg-per.vcf

QUERY_VCF=$1
OUTPUT_PREFIX=$2
REF=/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/refflow-exp/reference_flow/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
python $NC2/software/hap.py/bin/hap.py $TRUE_VCF $1 -o $OUTPUT_PREFIX --threads 24 -r $REF -f $CONF_REGION
