import argparse
import glob
from multiprocessing import Process, Queue, Manager
import os
import pickle
import sys
import time

import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import RDLogger
from rdkit.Chem.Descriptors import ExactMolWt

RDLogger.DisableLog("rdApp.*")


def generate_fp(key, smiles):
    try:
        fn = os.path.join(args.dirname, key)
        if os.path.exists(fn):
            return
        if "." in smiles:
            return
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return
        fs = AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024)
        np.save(fn, np.array(fs))
    except Exception as _:
        return


def worker(queue):
    done = queue.empty()
    while not done:
        args = queue.get()
        lines, idx = args
        for line in lines:
            key, smiles = line
            generate_fp(key, smiles)
        done = queue.empty()
    return


def initialize_queue(data, NTASKS=20, shared_dict=None):
    queue = Queue()
    length = len(data)
    for i in range(NTASKS):
        length_per_proc = length // NTASKS
        st = i * length_per_proc
        end = (i + 1) * length_per_proc if i != NTASKS - 1 else length
        data_per_proc = data[st:end]
        if shared_dict is not None:
            args = (data_per_proc, shared_dict, i)
        else:
            args = (data_per_proc, i)
        queue.put(args)
    return queue


def initialize_proc(queue, fn, NCPU=4):
    procs = []
    for _ in range(NCPU):
        proc = Process(target=fn, args=(queue,))
        procs.append(proc)
        proc.start()
        time.sleep(0.5)
    return procs


def write_result(keys, smiles_list, clusters, fn):
    with open(fn, "w") as w:
        for c in clusters:
            ligand = keys[c[0]]
            smiles = smiles_list[c[0]]
            w.write(f"{ligand}\t{smiles}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn", type=str, help="filename")
    parser.add_argument("--dirname", type=str, help="directory name")
    parser.add_argument("--ncpu", type=int, help="ncpu")
    args, _ = parser.parse_known_args()

    NCPU = args.ncpu
    NTASKS = args.ncpu * 3

    with open(args.fn, "r") as f:
        lines = f.readlines()
        lines = [line.split() for line in lines]

    if not os.path.exists(args.dirname):
        os.mkdir(args.dirname)
    queue = initialize_queue(lines, NTASKS)
    procs = initialize_proc(queue, worker, NCPU)
    for proc in procs:
        proc.join()
