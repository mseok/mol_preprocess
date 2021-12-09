# Crawling

ligand expo: http://ligand-expo.rcsb.org/ld-download.html  
ligand stereo SMILES data: http://ligand-expo.rcsb.org/dictionaries/Components-smiles-stereo-cactvs.smi  
ligand SMILES data: http://ligand-expo.rcsb.org/dictionaries/Components-smiles-cactvs.smi  

1. `download_protein_pdb.py`  
  Download pdb files from RCSB. Each pdb file contain corresponding ligand id at `HET`.
2. `find_ligand_of_pdb.py`  
  Get ligand id from each pdb files.
3. `get_ligand_set.py`  
  Gather the ligand information to single file with comma-separated.
4. `download_ligand_smi.py`  
  Download canonical SMILES from the pdbbind homepage.
