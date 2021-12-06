from random import shuffle, uniform
from typing import List, Optional, Tuple

from rdkit.Chem import Atom, Bond, Mol
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdMolAlign import CalcRMS
from rdkit.Chem.rdMolTransforms import GetDihedralDeg, SetDihedralDeg


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


def dihedral(mol: Mol, variance: float = 1.0, threshold: float = 0.5) -> Mol:
    rotable_bonds: List[Tuple[int]] = get_rotable_bonds(mol)
    _mol = deepcopy(mol)
    atoms = _mol.GetAtoms()
    _conf: Chem.Conformer = _mol.GetConformer()
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
        dihedral = GetDihedralDeg(_conf, *atom_indices)
        diff = uniform(-variance, variance)
        positions = _conf.GetPositions()
        SetDihedralDeg(_conf, *atom_indices, dihedral + diff * 5)
        rmsd = CalcRMS(_mol, mol)
        if rmsd >= threshold:
            for idx, pos in enumerate(positions):
                _conf.SetAtomPosition(idx, pos)
            break
    return _mol, rmsd


if __name__ == "__main__":
    import os
    import sys

    file = sys.argv[1]
    mol: Mol = Chem.MolFromMol2File(file)
    if not os.path.exists("output"):
        os.mkdir("output")

    print("dihedral")
    for i in range(N):
        _mol, rmsd = dihedral(mol, 5)
        # print(rmsd)
        # w = Chem.SDWriter(f"output/{i}/dihedral.sdf")
        # w.write(_mol)
        # w.close()
