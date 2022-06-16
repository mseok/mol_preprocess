#!/usr/bin/env python
import os
from copy import deepcopy
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Mol, rdDistGeom, rdMolAlign


def read_mol(file: str, add_hs: bool = True) -> Mol:
    mol = Chem.SDMolSupplier(file)[0]
    if add_hs:
        mol = Chem.AddHs(mol, addCoords=True)
    return mol


def get_parameters(
    seed: int = 2022,
    use_random_coords: bool = True,
    num_theads: int = 1,
    prune_rms_thresh: float = 1,
    box_size_mult: int = 2,
    clear_confs: bool = False,
):
    ps = AllChem.srETKDGv3()
    ps.randomSeed = seed
    ps.useRandomCoords = use_random_coords
    ps.numThreads = num_theads
    ps.pruneRmsThresh = prune_rms_thresh
    ps.boxSizeMult = box_size_mult
    ps.clearConfs = clear_confs
    return ps


def main():
    path = DATA_DIR / KEY
    path.mkdir(exist_ok=True)

    ligand = read_mol(f"../refined-set.processed/{KEY}/{KEY}_ligand.sdf")
    ps = get_parameters(clear_confs=False)

    confids = AllChem.EmbedMultipleConfs(ligand, 500, ps)

    rdMolAlign.AlignMolConformers(ligand)
    file_path = path / f"{KEY}.sdf"
    writer = Chem.SDWriter(str(file_path))
    for confid in confids[1:]:
        writer.write(ligand, confId=confid)
    writer.close()
    return


if __name__ == "__main__":
    import sys

    DATA_DIR = Path("./etkdg")
    KEY = sys.argv[1]
    main()
