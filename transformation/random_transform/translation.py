import argparse
import glob
from copy import deepcopy
from multiprocessing import Pool
from os import mkdir
from os.path import exists, join
from random import uniform
from typing import List

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Mol

RDLogger.DisableLog("rdApp.*")


def translate(mol: Mol, _max: float = 2.0, _min: float = 0.5) -> Mol:
    _mol = deepcopy(mol)
    conf: Chem.Conformer = _mol.GetConformer()
    dist = 0
    while dist <= _min:
        variance = np.random.uniform(-1, 1, (3,))
        dist = np.sqrt((variance ** 2).sum())
        if dist > _max:
            variance = uniform(_min, _max) * variance / dist
    pos = conf.GetPositions()
    for idx in range(conf.GetNumAtoms()):
        conf.SetAtomPosition(idx, (pos[idx] + variance))
    return _mol, dist


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
            _mol, rmsd = translate(mol, args.max, args.min)
            key = file.split("/")[-1].split("_")[0]
            new_file = join(args.out_dir, f"{key}_translation_{count}.sdf")
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

    main(args)
