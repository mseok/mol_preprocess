import argparse
import glob
import os
from multiprocessing import Pool
from typing import List

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def get_rmsd(mol1: Chem.Mol, mol2: Chem.Mol) -> float:
    pos1 = mol1.GetConformer().GetPositions()
    pos2 = mol2.GetConformer().GetPositions()
    rmsd = np.sqrt(((pos1 - pos2) ** 2).sum() / pos1.shape[0])
    return rmsd


def worker(key: str) -> None:
    file1 = f"/home/wykgroup/mseok/data/refined-set/{key}/{key}_ligand.sdf"
    try:
        mol1 = Chem.SDMolSupplier(file1)[0]
    except Excpetion:
        pass
    finally:
        if mol1 is None:
            return
        mol1 = Chem.RemoveHs(mol1)
        files2 = glob.glob(f"./{args.target_dir}/{key}_*_*.sdf")
        for file2 in files2:
            mol2 = Chem.SDMolSupplier(file2)[0]
            if mol2 is None:
                os.remove(file2)
                # print(file2, "None")
                return
            rmsd = get_rmsd(mol1, mol2)
            if rmsd < args.max and rmsd >= args.min:
                # print(file2, rmsd)
                continue
            else:
                print(file2, "RSMD: ", rmsd)
                # os.remove(file2, "RMSD")


def mp(ligands: List[str], ncpu: int = 4) -> None:
    pool = Pool(ncpu)
    results = pool.map_async(worker, ligands)
    results.wait()
    pool.close()
    pool.join()
    return


def main(args: argparse.Namespace) -> None:
    ligands = glob.glob(os.path.join(args.pdbbind_dir, "????"))
    ligands = [ligand.split("/")[-1] for ligand in ligands]
    print(ligands)
    mp(ligands, 12)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdbbind_dir", type=str, default="/home/wykgroup/mseok/data/refined-set"
    )
    parser.add_argument("-d", "--target_dir", type=str, default="./sdfs")
    parser.add_argument("--max", type=float, default=2)
    parser.add_argument("--min", type=float, default=0.5)
    parser.add_argument("--ncpu", type=int, default=4)
    args = parser.parse_args()

    main(args)
