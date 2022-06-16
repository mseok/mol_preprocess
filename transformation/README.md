# Transformation
## Random transformation with RDKit

1. `translation.py`  
  Translation in certain RMSD range.
2. `rotation.py`  
  Rotation in certain RMSD range.
3. `dihedral.py`  
  Dihedral transformation based on the rotable bonds.
4. `check_rmsd.py`  
  Compare the RMSD between original molecule and generated molecule.

## ETKDG conformer generation with RDKit
1. `ETKDG.py`  
  ETKDG conformer generation with RDKit.
2. `local_opt.sh`  
  Local optimization of generated conformers with [SMINA](https://www.google.com/search?client=safari&rls=en&q=smina&ie=UTF-8&oe=UTF-8).
3. `run_dockrmsd.sh`  
  RMSD based filteration with [dockrmsd](https://zhanggroup.org/DockRMSD/).
4. `run_scoring.sh`  
  SMINA scoring of generated conformers.
