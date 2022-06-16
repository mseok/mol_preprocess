#!/bin/bash
PDBBIND_V2020=/home/mseok/data/PDBbind/pdbbind_v2020/refined-set.processed

dockrmsd() {
  dir=$(dirname $1)
  pdb=$(basename $dir)
  filename=${dir}/$(basename $1 .sdf)_.mol2
  filename_all=${dir}/$(basename $1 .sdf)_*.mol2
  obabel $1 -O $filename -m -xu 2>/dev/null
  for file in $(ls $filename_all);
  do
    echo -en "$file\t"
    DockRMSD $file ${PDBBIND_V2020}/${pdb}/${pdb}_ligand.mol2 \
      | grep "Calculated Docking RMSD" \
      | awk '{print $4}'
  done
}


export -f dockrmsd

# Single execution
# dockrmsd ./etkdg/10gs/optimized.sdf
# Multiprocessing
# parallel -j4 -a <(ls ./etkdg/*/optimized.sdf -d) dockrmsd > docking_rmsd.txt
