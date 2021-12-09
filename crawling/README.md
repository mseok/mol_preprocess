# Crawling

ligand expo: http://ligand-expo.rcsb.org/ld-download.html  
ligand stereo SMILES data: http://ligand-expo.rcsb.org/dictionaries/Components-smiles-stereo-cactvs.smi  
ligand SMILES data: http://ligand-expo.rcsb.org/dictionaries/Components-smiles-cactvs.smi  

1. `download_protein.py`  
  Download cif or pdb files from RCSB. Each pdb file contain corresponding ligand id at `HET`.
2. `download_ligand_sdf.py`  
  Download sdf for corresponding existing cif file from RCSB.
3. `find_ligand_of_pdb.py`  
  Get ligand id from each pdb files.
4. `get_ligand_set.py`  
  Gather the ligand information to single file with comma-separated.
5. `download_ligand_smi.py`  
  Download canonical SMILES from the pdbbind homepage.
