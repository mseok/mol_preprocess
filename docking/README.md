# Docking

1. `docking.py`  
  Execute the docking with `SMINA`(needs `open-babel` for preprocessing).  
  Note that you can use almost any file type for docking input and outputs.  
  In this script, we used `pdbqt` for both ligand and protein(receptor).  
  But you can also use `pdb` or `sdf` for input and output.
2. `split_sdf.py`  
  Use `open-babel` to split the docking result `pdbqt` files into `sdf`.
3. `run_smina.sh`
  `SMINA` execution code in bash language. Concise multiprocessing with GNU parallel.
