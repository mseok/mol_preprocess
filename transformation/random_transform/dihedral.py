import argparse
import glob
from copy import deepcopy
from multiprocessing import Pool
from os import mkdir
from os.path import exists, join
from random import shuffle, uniform
from typing import List, Optional, Tuple

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Atom, Bond, Mol
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdMolAlign import CalcRMS
from rdkit.Chem.rdMolTransforms import GetDihedralRad, SetDihedralRad

RDLogger.DisableLog("rdApp.*")


def get_rotable_bonds(mol: Mol) -> List[Tuple[int]]:
    rotable_bonds = mol.GetSubstructMatches(RotatableBondSmarts)
    rotable_bonds = list(sorted(rotable_bonds))
    return rotable_bonds


def get_adjacent_atoms(
    atom: Atom, exclude_atom_idx: Optional[int] = None
) -> List[Atom]:
    bonds: List[Bond] = list(atom.GetBonds())
    bonded_atoms = [bond.GetOtherAtom(atom) for bond in bonds]
    bonded_atoms = [atom for atom in bonded_atoms if atom.GetIdx() != exclude_atom_idx]
    return bonded_atoms


def dihedral(mol: Mol, _max: float = 2.0, _min: float = 0.5) -> Mol:
    rotable_bonds: List[Tuple[int]] = get_rotable_bonds(mol)
    if not rotable_bonds:
        return
    _mol = deepcopy(mol)
    atoms = _mol.GetAtoms()
    conf: Chem.Conformer = _mol.GetConformer()
    pos = conf.GetPositions()
    _new_pos = pos.copy()
    n = len(rotable_bonds)
    variance = np.pi / n
    rmsd = 0
    count = 0
    while rmsd <= _min:
        shuffle(rotable_bonds)
        for indices in rotable_bonds:
            begin_idx, end_idx = indices
            begin_atom: Atom = atoms[begin_idx]
            end_atom: Atom = atoms[end_idx]
            begin_neighbors = get_adjacent_atoms(begin_atom, end_idx)
            end_neighbors = get_adjacent_atoms(end_atom, begin_idx)
            shuffle(begin_neighbors)
            shuffle(end_neighbors)
            atom_indices = [
                begin_neighbors[0].GetIdx(),
                begin_idx,
                end_idx,
                end_neighbors[0].GetIdx(),
            ]
            dihedral = GetDihedralRad(conf, *atom_indices)
            diff = uniform(0, variance)
            SetDihedralRad(conf, *atom_indices, dihedral + diff)
            new_pos = conf.GetPositions()
            rmsd = np.sqrt(((new_pos - pos) ** 2).sum() / pos.shape[0])
            if rmsd > _max:
                rmsd = np.sqrt(((_new_pos - pos) ** 2).sum() / pos.shape[0])
                break
            _new_pos = new_pos
        count += 1
        if count > 500:
            return None
    for idx, pos in enumerate(_new_pos):
        conf.SetAtomPosition(idx, pos)
    return _mol, rmsd


def mol2sdf(mol: Mol, file: str) -> None:
    w = Chem.SDWriter(file)
    w.write(mol)
    w.close()
    return


def worker(file: str) -> None:
    count = 0
    success_count = 0
    try:
        mol = Chem.SDMolSupplier(file)[0]
    except Exception as e:
        print(file, e)
    finally:
        if mol is None:
            return
        while count < args.num_gen ** 2:
            result = dihedral(mol, args.max, args.min)
            if result:
                _mol, rmsd = result
                key = file.split("/")[-1].split("_")[0]
                new_file = join(args.out_dir, f"{key}_dihederal_{success_count}.sdf")
                mol2sdf(_mol, new_file)
                success_count += 1
                if success_count == args.num_gen:
                    break
            count += 1
        return


def mp(ligands: List[str], ncpu: int = 4) -> None:
    pool = Pool(ncpu)
    results = pool.map_async(worker, ligands)
    results.wait()
    pool.close()
    pool.join()
    return


def main(args: argparse.Namespace) -> None:
    ligands = glob.glob(join(args.pdbbind_dir, "????/????_ligand.sdf"))
    if not exists(args.out_dir):
        mkdir(args.out_dir)
    mp(ligands, args.ncpu)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--pdbbind_dir", type=str, default="/home/wykgroup/mseok/data/refined-set"
    )
    parser.add_argument("-o", "--out_dir", type=str, default="./sdfs")
    parser.add_argument("--max", type=float, default=2)
    parser.add_argument("--min", type=float, default=0.5)
    parser.add_argument("--ncpu", type=int, default=4)
    parser.add_argument("-n", "--num_gen", type=int, default=5)
    args = parser.parse_args()

    main(args)
