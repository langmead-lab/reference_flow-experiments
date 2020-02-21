usage(){
echo "Usage: $0 <sam>"
echo "      sam: FILENAME"
echo "      please modify this script if you'd like to use another MAPQ cutoff"
#echo "Usage: $0 <threshold> <sam>"
#echo "      threshold: INT, outputs will be UNMAPPED+MAPQ_0:<threshold-1> and MAPQ_<threshold>:"
#echo "      sam: FILENAME"
}
if [ "$#" -ne 1 ]
then
    echo "ERROR: incorrect number of arguments"
    usage
    exit
fi

filename=$(basename -- "$1")
filename="${filename%.*}"

#: changes above three lines manually for different MAPQ threshold
THRSD="10"
awk '{ if ($5 < 10 || $1 ~ /^@/) { print } }' $1 > $filename-mapql${THRSD}.sam
awk '{ if ($5 >= 10 || $1 ~ /^@/) { print } }' $1 > $filename-mapqgeq${THRSD}.sam

module load samtools
samtools fastq $filename-mapql${THRSD}.sam > $filename-mapql${THRSD}.fq
