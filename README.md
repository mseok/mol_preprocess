# mol_preprocess

```bash
.
├── crawling
│   ├── download_ligand_sdf.py
│   ├── download_ligand_smi.py
│   ├── download_protein.py
│   ├── find_ligand_of_pdb.py
│   ├── get_ligand_set.py
│   └── README.md
├── docking
│   ├── docking.py
│   └── split_sdf.py
├── LICENSE
├── preprocess
│   └── preprocess.py
├── README.md
├── similarity
│   ├── ligand
│   │   ├── cat_fps.py
│   │   ├── extract_valids.py
│   │   ├── generate_fps.py
│   │   ├── gpu_similarity.py
│   │   └── README.md
│   └── protein_clustering
│       ├── clustering.sh
│       ├── filter.py
│       ├── parser.py
│       ├── pdb2fasta.py
│       ├── README.md
│       └── target_list.py
└── transformation
    ├── check_rmsd.py
    ├── dihedral.py
    ├── rotation.py
    └── translation.py
```

1. crawling  
  Crawling SMILES, PDB, SDF and etc from the [rcsb pdb](https://www.rcsb.org) or [pdbbind](http://www.pdbbind.org.cn).  
2. docking  
  Scripts for docking with SMINA.
3. preprocess  
  Preprocess code for docked data.
4. similarity  
  - Ligand  
    Small molecules filtering based on the Tanimoto Similarity.
  - Protein  
    Protein filtering based on the cd-hit.
5. transformation  
  Code for tranform the molecule conformation.
