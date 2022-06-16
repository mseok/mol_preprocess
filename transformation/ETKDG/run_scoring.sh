#!/bin/bash
PDBBIND_V2020=/home/mseok/data/PDBbind/pdbbind_v2020/refined-set.processed

scoring() {
  echo -en "$1\t"
  pdb=$(basename $(dirname $1))
  smina \
    --local_only \
    --score_only \
    -r ${PDBBIND_V2020}/${pdb}/${pdb}_protein.pdb \
    -l $1 \
    --autobox_ligand ${PDBBIND_V2020}/${pdb}/${pdb}_ligand.sdf \
    | grep "Affinity" \
    | awk '{print $2}'
}


export -f scoring

parallel -j4 -a <(echo ./etkdg/*/*.mol2 | xargs ls | shuf) scoring > ./docking_score.txt
