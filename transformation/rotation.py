import argparse
import glob
from copy import deepcopy
from multiprocessing import Pool
from os import mkdir
from os.path import exists, join
from random import shuffle, uniform
from typing import List

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Mol
from rdkit.Chem.rdMolAlign import CalcRMS

RDLogger.DisableLog("rdApp.*")


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


def rotate(mol: Mol, _max: float = 2.0, _min: float = 0.5) -> Mol:
    _mol = deepcopy(mol)
    com = get_center_of_mass(_mol)
    conf = _mol.GetConformer()
    pos = conf.GetPositions()
    _pos = pos - com
    new_pos = _pos.copy()
    indices = list(range(3))
    rmsd = 0
    variance = np.pi
    while rmsd <= _min:
        shuffle(indices)
        for idx in indices:
            _new_pos = new_pos @ ROTATIONS[idx](uniform(0, variance))
            rmsd = np.sqrt(((_new_pos - _pos) ** 2).sum() / pos.shape[0])
            if rmsd > _max:
                rmsd = np.sqrt(((new_pos - _pos) ** 2).sum() / pos.shape[0])
                break
            new_pos = _new_pos
        variance = variance * 0.95
    new_pos += com
    for idx, pos in enumerate(new_pos):
        conf.SetAtomPosition(idx, pos)
    return _mol, rmsd


def mol2sdf(mol: Mol, file: str) -> None:
    w = Chem.SDWriter(file)
    w.write(mol)
    w.close()
    return


def worker(file: str) -> None:
    count = 0
    try:
        mol = Chem.SDMolSupplier(file)[0]
    except Exception as e:
        print(file, e)
    finally:
        if mol is None:
            return
        while count < args.num_gen:
            _mol, rmsd = rotate(mol, args.max, args.min)
            key = file.split("/")[-1].split("_")[0]
            new_file = join(args.out_dir, f"{key}_rotation_{count}.sdf")
            mol2sdf(_mol, new_file)
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

    ROTATIONS = {0: rot_ar_x, 1: rot_ar_y, 2: rot_ar_z}

    main(args)
