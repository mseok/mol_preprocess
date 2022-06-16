#!/bin/bash

# Directory
DATA_HOME="/home/mseok/data/PDBbind/pdbbind_v2020"
REFINED_DIR="${DATA_HOME}/refined-set.processed"
RESULT_DIR="${DATA_HOME}/analyze/docking"

# SMINA Values
TIMEOUT=600
NUM_MODES=9
SEED=0

_smina() {
  pdb=$(echo $1 | cut -d "_" -f1)
  idx=$(echo $1 | cut -d "_" -f2)
  if [ ! -d "${RESULT_DIR}/${pdb}" ]; then
    mkdir "${RESULT_DIR}/${pdb}" 
  fi

  timeout $TIMEOUT \
    smina \
    -r "${REFINED_DIR}/${pdb}/${pdb}_protein.pdb" \
    -l "${REFINED_DIR}/${pdb}/${pdb}_ligand.sdf" \
    -o "${RESULT_DIR}/${pdb}/optimized_${idx}.sdf" \
    --autobox_ligand "${REFINED_DIR}/${pdb}/${pdb}_ligand.sdf" \
    --num_modes $NUM_MODES \
    --seed $SEED
}


export -f _smina

# Single execution
# _smina $1
# Multiprocessing
# parallel -j4 -a <(ls ${REFINED_DIR}) _smina
