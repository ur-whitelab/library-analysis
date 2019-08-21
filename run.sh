#!/bin/bash
if [ "$1" = "" ] || [ "$2" = "" ]
then
    echo "Missing positional parameter(s)! Usage: run.sh [fastq_filename] [barcode_filename]"
    exit 1
fi
BASENAME="$1"
DIRNAME="${BASENAME%.*}"
mkdir $DIRNAME
cp $1 "$DIRNAME"
cp $2 "$DIRNAME"
cd $DIRNAME
python3 ../parse_fastq.py $1 $2
for a in $(ls *SEQS.txt)
do
    python3 ../calc_summary_composition.py $a $2
done
