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
    args, _ = parser.parse_known_args()
    return args


def mp_pool(ncpu: int, fns: List[str]) -> None:
    pool = Pool(ncpu)
    r = pool.map_async(run, fns)
    r.wait()
    pool.close()
    pool.join()
    return


def main(args: argparse.Namespace) -> None:
    fns = glob.glob(os.path.join(args.result_pdbqt_dir, "*.pdbqt"))
    fns = fns[args.start : args.end]
    mp_pool(args.ncpu, fns)
    return


if __name__ == "__main__":
    args = parser()
    if not os.path.exists(args.result_sdf_dir):
        os.makedirs(args.result_sdf_dir, exist_ok=True)
    main(args)
