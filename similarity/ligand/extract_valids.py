import glob
import os
import pickle

import numpy as np

fns = glob.glob("./similarity/*.npy")
hits = []
rands = []
for fn in fns:
    target = fn.split("/")[-1].split(".")[0]
    with open(os.path.join("samples", target + ".pkl"), "rb") as f:
        keys = pickle.load(f)
    indices = np.load(fn)
    # fps_hits_candidates/adrb2-hits_187.npy
    for idx in indices:
        fp_fn = keys[idx]
        new_key = fp_fn.split("/")[-1].split(".")[0]
        if "STOCK" in new_key:
            rands.append(new_key)
        else:
            hits.append(new_key)

with open("total_hits.txt", "w") as w:
    for hit in hits:
        w.write(hit + "\n")
with open("total_rands.txt", "w") as w:
    for rand in rands:
        w.write(rand + "\n")
