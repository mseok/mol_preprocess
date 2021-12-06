import argparse
import glob
import os
import sys
from multiprocessing import Pool
from pathlib import PosixPath
from typing import List


def run(fn: str) -> None:
    key = "_".join(fn.split("/")[-1].split(".")[0].split("_"))
    os.system(f"obabel {fn} -O ./result_sdf/{key}_.sdf -xrp -m")
    return


def parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result_pdbqt_dir", type=PosixPath, default="./result_pdbqt")
    parser.add_argument("--result_sdf_dir", type=PosixPath, default="./result_sdf")
    parser.add_argument("--ncpu", type=int, default=4)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    FLAGS, _ = parser.parse_known_args()
    return FLAGS


def mp_pool(ncpu: int, fns: List[str]) -> None:
    pool = Pool(ncpu)
    r = pool.map_async(run, fns)
    r.wait()
    pool.close()
    pool.join()
    return


def main(FLAGS: argparse.Namespace) -> None:
    fns = glob.glob(os.path.join(FLAGS.result_pdbqt_dir, "*.pdbqt"))
    fns = fns[FLAGS.start : FLAGS.end]
    mp_pool(FLAGS.ncpu, fns)
    return


if __name__ == "__main__":
    FLAGS = parser()
    if not os.path.exists(FLAGS.result_sdf_dir):
        os.makedirs(FLAGS.result_sdf_dir, exist_ok=True)
    main(FLAGS)
