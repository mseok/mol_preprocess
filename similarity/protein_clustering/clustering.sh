#!/bin/bash
if ! command -v curl &> /dev/null
then
    echo "'cd-hit' could not be found. You need to install 'cd-hit' to cluster the proteins."
    exit 1
fi

CUTOFF=$1
PDBBIND_DIR=$2
declare -A CUTOFF2N
CUTOFF2N["0.4"]="2"
CUTOFF2N["0.5"]="3"
CUTOFF2N["0.6"]="4"
CUTOFF2N["0.7"]="5"
CUTOFF2N["0.8"]="5"
CUTOFF2N["0.9"]="5"
CUTOFF2N["1.0"]="5"
N="${CUTOFF2N[$CUTOFF]}"
if [ -z $N ]
then
  echo "CUTOFF should be included in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]"
  exit 1
fi

./pdb2fasta.py -i ./target_list.txt -o ./target.fasta -d $PDBBIND_DIR
cd-hit -i ./target.fasta -o target -c $CUTOFF -n $N
# Valid keys
./filter.py -i ./target.clstr -o valid_target.txt -d $PDBBIND_DIR --filter_valid
# Invalid keys
# ./filter.py -i ./target.clstr -o valid_target.txt -d $PDBBIND_DIR
