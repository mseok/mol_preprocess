import argparse
import glob
from multiprocessing import Pool
import os
import pickle

import numpy as np


def worker(fn):
    return fn, np.load(fn)


def main(result_fn, fns, ncpu):
    pool = Pool(ncpu)
    results = pool.map_async(worker, fns)
    pool.close()
    pool.join()
    data = results.get()
    used_fns, bits = zip(*data)
    np.save(result_fn + ".npy", bits)
    with open(result_fn + ".pkl", "wb") as w:
        pickle.dump(list(used_fns), w)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fps_dir", type=str, help="fp directories")
    parser.add_argument("--result_fn", type=str, help="result filename")
    parser.add_argument("--ncpu", type=int, help="ncpu", default=1)
    args, _ = parser.parse_known_args()

    PATH = "samples"
    if not os.path.exists(PATH):
        os.mkdir(PATH)
    if "hits" in args.fps_dir:
        # hits
        fns = glob.glob(os.path.join(args.fps_dir, "*.npy"))
        targets = list(set([fn.split("/")[-1].split("-")[0] for fn in fns]))
        for target in targets:
            fns = glob.glob(os.path.join(args.fps_dir, f"{target}*.npy"))
            fns = sorted(fns)
            result_fn = os.path.join(PATH, target)
            main(result_fn, fns, args.ncpu)
    else:
        fns = glob.glob(os.path.join(args.fps_dir, "*.npy"))
        fns = sorted(fns)
        result_fn = os.path.join(PATH, args.result_fn)
        main(result_fn, fns, args.ncpu)
