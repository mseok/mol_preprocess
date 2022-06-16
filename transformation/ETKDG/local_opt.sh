#!/bin/bash
PDBBIND_V2020=/home/mseok/data/PDBbind/pdbbind_v2020/refined-set.processed

smina_local_opt() {
  target=${PDBBIND_V2020}/${1}/${1}_protein.pdb
  ligand=${PDBBIND_V2020}/${1}/${1}_ligand.sdf
  new_file=$(dirname $1)/optimized.sdf
  smina \
    -r $target \
    -l $1 \
    --autobox_ligand $ligand \
    --autobox_add 4 \
    -o $new_file \
    --local_only \
    --minimize
  done
}


export -f smina_local_opt

# Single execution
# smina_local_opt ./etkdg/10gs
# Multiprocessing
# parallel -j4 -a <(ls ./etkdg/* -d) smina_local_opt
