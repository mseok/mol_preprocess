from copy import deepcopy
from random import shuffle

import numpy as np
from rdkit.Chem import Mol

ROTATIONS = {0: rot_ar_x, 1: rot_ar_y, 2: rot_ar_z}


def rot_ar_x(radi: float) -> np.ndarray:
    return np.array(
        [
            [1, 0, 0],
            [0, np.cos(radi), -np.sin(radi)],
            [0, np.sin(radi), np.cos(radi)],
        ],
        dtype=np.double,
    )


def rot_ar_y(radi: float) -> np.ndarray:
    return np.array(
        [
            [np.cos(radi), 0, np.sin(radi)],
            [0, 1, 0],
            [-np.sin(radi), 0, np.cos(radi)],
        ],
        dtype=np.double,
    )


def rot_ar_z(radi: float) -> np.ndarray:
    return np.array(
        [
            [np.cos(radi), -np.sin(radi), 0],
            [np.sin(radi), np.cos(radi), 0],
            [0, 0, 1],
        ],
        dtype=np.double,
    )


def get_center_of_mass(mol: Mol) -> np.ndarray:
    weights = []
    for atom in mol.GetAtoms():
        weight = Chem.Atom(atom.GetSymbol()).GetMass()
        weights.append(weight)
    conf = mol.GetConformer()
    pos = conf.GetPositions()
    pos = np.array(pos)
    weights = np.array(weights)
    weights = weights.reshape(-1, 1)
    com = (weights * pos).sum(0) / weights.sum()
    return com


def rotate(mol: Mol) -> Mol:
    _mol = deepcopy(mol)
    com = get_center_of_mass(_mol)
    conf = _mol.GetConformer()
    pos = conf.GetPositions()
    _pos = pos - com
    new_pos = _pos.copy()
    indices = list(range(3))
    shuffle(indices)
    for idx in indices:
        _new_pos = new_pos @ ROTATIONS[idx](uniform(0, np.pi / 4))
        rmsd = np.sqrt(((_new_pos - _pos) ** 2).sum() / pos.shape[0])
        if rmsd >= 2:
            break
        new_pos = _new_pos
    new_pos += com
    for idx, pos in enumerate(new_pos):
        conf.SetAtomPosition(idx, pos)
    return _mol


if __name__ == "__main__":
    import os
    import sys

    file = sys.argv[1]
    mol: Mol = Chem.MolFromMol2File(file)
    if not os.path.exists("output"):
        os.mkdir("output")

    print("rotate")
    for i in range(N):
        print("rotation: ", i)
        _mol = rotate(mol)
        # print(CalcRMS(_mol, mol))
        # w = Chem.SDWriter(f"output/{i}/rotation.sdf")
        # w.write(_mol)
        # w.close()
