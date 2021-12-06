import argparse
import os
import pickle
import sys
import time
from multiprocessing import Pool
from pathlib import PosixPath
from subprocess import DEVNULL, PIPE, Popen
from typing import AnyStr, List, Tuple, Union

from rdkit import Chem
from rdkit.Chem import AllChem, SDMolSupplier, SDWriter
from rdkit.Chem.rdmolfiles import PDBWriter


def run_proc(cmd: List[str], timeout: int = 1200) -> Tuple[AnyStr, AnyStr]:
    proc = Popen(cmd, stdout=DEVNULL, stderr=PIPE, encoding="utf-8")
    out, err = proc.communicate(timeout=timeout)
    return out, err


def save_as_pdb(mol: Chem.Mol, pdb_name: str, min_idx: Union[int, None] = None) -> None:
    w = PDBWriter(pdb_name)
    if min_idx is not None:
        w.write(mol, min_idx)
    else:
        w.write(mol)
    w.close()
    return


def uff_calc(mol: Chem.Mol) -> Tuple[Chem.Mol, int]:
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=20)
    cenergy = []
    for conf in cids:
        converged = not AllChem.UFFOptimizeMolecule(mol, confId=conf)
        cenergy.append(AllChem.UFFGetMoleculeForceField(mol, confId=conf).CalcEnergy())
    min_idx = cenergy.index(min(cenergy))
    mol = Chem.RemoveHs(mol)
    return mol, min_idx


def docking(keys: List[str]) -> None:
    pdb_id = keys[0]
    protein = os.path.join(FLAGS.pdbbind_dir, pdb_id, f"{pdb_id}_protein_nowater.pdb")
    protein_pdbqt = os.path.join(
        FLAGS.pdbbind_dir, pdb_id, f"{pdb_id}_protein_nowater.pdbqt"
    )
    ligand = os.path.join(FLAGS.pdbbind_dir, pdb_id, f"{pdb_id}_ligand.sdf")
    if not os.path.exists(ligand):
        protein = os.path.join(
            FLAGS.pdbbind_general_dir, pdb_id, f"{pdb_id}_protein_nowater.pdb"
        )
        protein_pdbqt = os.path.join(
            FLAGS.pdbbind_general_dir, pdb_id, f"{pdb_id}_protein_nowater.pdbqt"
        )
        ligand = os.path.join(FLAGS.pdbbind_general_dir, pdb_id, f"{pdb_id}_ligand.sdf")
    log = os.path.join(FLAGS.log_dir, f"{pdb_id}.log")
    pdb = os.path.join(FLAGS.pdb_dir, f"{pdb_id}.pdb")
    pdbqt = os.path.join(FLAGS.pdbqt_dir, f"{pdb_id}.pdbqt")
    result_pdbqt = os.path.join(FLAGS.result_pdbqt_dir, f"{pdb_id}.pdbqt")

    # Pass condition
    if os.path.exists(result_pdbqt) and os.path.getsize(result_pdbqt):
        return

    # Generate 3D structure of ligand
    m = SDMolSupplier(ligand)[0]
    if m is None:
        return
    save_as_pdb(m, pdb)

    # Pdb to pdbqt (both of ligand and protein)
    if not os.path.exists(pdbqt):
        command = f"obabel {pdb} -O {pdbqt}".split()
        run_proc(command)
    if not os.path.exists(protein_pdbqt):
        command = f"obabel {protein} -O {protein_pdbqt}".split()
        run_proc(command)

    command = f"smina \
            -r {protein_pdbqt} \
            -l {pdbqt} \
            --autobox_ligand {ligand} \
            --autobox_add 8 \
            --exhaustiveness 8 \
            --log {log} \
            -o {result_pdbqt} \
            --cpu 1 \
            --num_modes 100 \
            --seed 0"
    _, err = run_proc(command.split())
    if len(err.split()) > 1:
        raise Exception("docking failed")
    return


def run_docking(keys: List[str]) -> None:
    try:
        docking(keys)
        print(keys, "done")
    except Exception as e:
        print(keys, e)
    finally:
        return


def mp_pool(ncpu: int, keys: List[List[str]]) -> None:
    pool = Pool(ncpu)
    results = pool.map_async(run_docking, keys)
    results.wait()
    pool.close()
    pool.join()
    results.get()
    return


def parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdbbind_dir",
        type=PosixPath,
        default="/home/wykgroup/mseok/data/PDBbind/2020/refined-set/",
    )
    parser.add_argument(
        "--pdbbind_general_dir",
        type=PosixPath,
        default="/home/wykgroup/mseok/data/PDBbind/2020/v2020-other-PL/",
    )
    parser.add_argument("--pdb_dir", type=PosixPath, default="./pdb")
    parser.add_argument("--pdbqt_dir", type=PosixPath, default="./pdbqt")
    parser.add_argument("--log_dir", type=PosixPath, default="./log")
    parser.add_argument("--result_pdbqt_dir", type=PosixPath, default="./result_pdbqt")
    parser.add_argument("--ncpu", type=int, default=4)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument(
        "--candidates_fn", type=PosixPath, default="./candidates/docking_candidates.txt"
    )
    FLAGS, _ = parser.parse_known_args()
    return FLAGS


def parse_candidate_fn(fn: PosixPath) -> List[List[str]]:
    with open(fn, "r") as f:
        lines = f.readlines()
        lines = [line.split() for line in lines]
    return lines


def main(FLAGS: argparse.Namespace) -> None:
    for key, value in vars(FLAGS).items():
        if "dir" in key and not os.path.exists(value):
            os.mkdir(value)
    keys = parse_candidate_fn(FLAGS.candidates_fn)
    if FLAGS.start is not None and FLAGS.end is not None:
        keys = keys[FLAGS.start : FLAGS.end]
    mp_pool(FLAGS.ncpu, keys)
    return


if __name__ == "__main__":
    FLAGS = parser()
    main(FLAGS)
