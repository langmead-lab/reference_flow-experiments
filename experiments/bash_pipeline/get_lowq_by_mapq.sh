usage(){
echo "Usage: $0 <sam> <thresh>"
echo "      sam: FILENAME"
echo "      thresh: INT, outputs will be in mapql<THRESH> and mapgeq<THRESH>"
#echo "Usage: $0 <threshold> <sam>"
#echo "      threshold: INT, outputs will be UNMAPPED+MAPQ_0:<threshold-1> and MAPQ_<threshold>:"
#echo "      sam: FILENAME"
}
if [ "$#" -ne 2 ]
then
    echo "ERROR: incorrect number of arguments"
    usage
    exit
fi

#filename=$(basename -- "$1")
filename="${1%.*}"


#: changes above three lines manually for different MAPQ threshold
THRSD="$2"

awk -v var="$THRSD" '{ if ($5 < var || $1 ~ /^@/) { print } }' $1 > $filename-mapql${THRSD}.sam
awk -v var="$THRSD" '{ if ($5 >= var || $1 ~ /^@/) { print } }' $1 > $filename-mapqgeq${THRSD}.sam

#Just make it mandatory to have samtools
#module load samtools
samtools fastq $filename-mapql${THRSD}.sam > $filename-mapql${THRSD}.fq
