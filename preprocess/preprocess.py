import argparse
import glob
import os
import pickle
import sys
import tempfile
import time
from multiprocessing import Pool
from pathlib import PosixPath
from typing import Dict, List, Optional

import numpy as np
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBIO import Select
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.SaltRemover import SaltRemover

RDLogger.DisableLog("rdApp.*")
from scipy.spatial import distance_matrix


def remove_water(m: Chem.Mol) -> Chem.Mol:
    remover = SaltRemover(defnData="[O]")
    return remover.StripMol(m)


def extract(ligand: Chem.Mol, pdb: str) -> Chem.Mol:
    parser = PDBParser()
    if not os.path.exists(pdb):
        return None
    structure = parser.get_structure("protein", pdb)
    ligand_positions = ligand.GetConformer().GetPositions()

    class ResidueSelect(Select):
        def accept_residue(self, residue):
            residue_positions = np.array(
                [
                    np.array(list(atom.get_vector()))
                    for atom in residue.get_atoms()
                    if "H" not in atom.get_id()
                ]
            )
            if len(residue_positions.shape) < 2:
                return 0
            n_residues = residue_positions.shape[0]
            n_ligands = ligand_positions.shape[0]
            rp = np.repeat(np.expand_dims(residue_positions, 1), n_ligands, 1)
            lp = np.repeat(np.expand_dims(ligand_positions, 0), n_residues, 0)
            dm = np.sqrt(((rp - lp) ** 2).sum(-1))
            min_dis = np.min(dm)
            if min_dis < 5.0:
                return 1
            else:
                return 0

    io = PDBIO()
    io.set_structure(structure)
    _, path = tempfile.mkstemp(suffix=".pdb", prefix="BS_tmp_", dir="/tmp/")
    io.save(path, ResidueSelect())
    m2 = Chem.MolFromPDBFile(path)
    os.unlink(path)
    return m2


def preprocessor(l: str) -> None:
    if not os.path.exists(l):
        return
    pdb_fn = l.split("/")[-1].split(".")[0]
    ligand_key = "_".join(pdb_fn.split("_")[:-1])
    target = ligand_key.split("-")[0]

    save_dir = FLAGS.save_dir
    if os.path.exists(f"{save_dir}/{pdb_fn}"):
        return
    bs_pdb_fn = os.path.join(FLAGS.pdbbind_dir, target, f"{target}_protein_nowater.pdb")

    temp = Chem.MolFromSmiles(key_to_smi[ligand_key])
    m1 = Chem.SDMolSupplier(l, sanitize=False)[0]
    m1 = Chem.RemoveHs(m1, sanitize=False)
    m1 = AssignBondOrdersFromTemplate(temp, m1)  # changed part
    if m1 is None:
        print(f"{key} no mol from sdf!")
        return

    # extract binding pocket
    m2 = extract(m1, bs_pdb_fn)
    if m2 is None:
        return
    m2 = remove_water(m2)

    if len(m1.GetConformers()) == 0:
        return
    if len(m2.GetConformers()) == 0:
        return
    with open(os.path.join(FLAGS.save_dir, pdb_fn), "wb") as fp:
        pickle.dump((m1, m2), fp)
    return


def worker(l: str) -> None:
    try:
        preprocessor(l)
        print(l, "done")
    except Exception as e:
        print(l, "failed", e)
    finally:
        return


def mp_pool(ncpu: int, keys: List[List[str]]) -> None:
    pool = Pool(ncpu)
    results = pool.map_async(worker, keys)
    results.wait()
    pool.close()
    pool.join()
    results.get()
    return


def parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--save_dir", type=PosixPath, default="./preprocessed")
    parser.add_argument(
        "--pdbbind_dir",
        type=PosixPath,
        default="/home/mseok/data/PDBbind/2020/refined-set/",
    )
    parser.add_argument(
        "--candidates_fn",
        type=PosixPath,
        default="candidates/preprocess_candidates.txt",
    )
    parser.add_argument(
        "--docking_candidates_fn",
        type=PosixPath,
        default="candidates/docking_candidates.txt",
    )
    parser.add_argument("--ncpu", type=int, default=4)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    FLAGS, _ = parser.parse_known_args()
    return FLAGS


def parse_candidate_fn(fn: PosixPath) -> List[List[str]]:
    with open(fn, "r") as f:
        lines = f.readlines()
        lines = [line.split()[0] for line in lines]
    return lines


def get_smiles_dic(fn: PosixPath) -> Dict[str, str]:
    with open(fn, "r") as f:
        lines = f.readlines()
        lines = [t.split() for t in lines]
        lines = [["-".join(t[:-1]), t[-1]] for t in lines]
    return dict(lines)


def main(FLAGS: argparse.Namespace) -> None:
    if not os.path.exists(FLAGS.save_dir):
        os.makedirs(FLAGS.save_dir, exist_ok=True)
    keys = parse_candidate_fn(FLAGS.candidates_fn)
    if FLAGS.start is not None and FLAGS.end is not None:
        keys = keys[FLAGS.start : FLAGS.end]
    global key_to_smi
    key_to_smi = get_smiles_dic(FLAGS.docking_candidates_fn)
    mp_pool(FLAGS.ncpu, keys)
    return


if __name__ == "__main__":
    global FLAGS
    FLAGS = parser()
    main(FLAGS)
