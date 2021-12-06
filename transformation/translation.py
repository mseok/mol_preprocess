from copy import deepcopy

import numpy as np
from rdkit import Chem
from rdkit.Chem import Mol


def translate(mol: Mol) -> Mol:
    _mol = deepcopy(mol)
    conf: Chem.Conformer = _mol.GetConformer()
    variance = np.random.uniform(-1, 1, (3,))
    dist = np.sqrt((variance ** 2).sum())
    if dist >= 2:
        variance = 2 * variance / dist
    pos = conf.GetPositions()
    for idx in range(conf.GetNumAtoms()):
        conf.SetAtomPosition(idx, (pos[idx] + variance))
    return _mol


if __name__ == "__main__":
    import os
    import sys

    file = sys.argv[1]
    mol: Mol = Chem.MolFromMol2File(file)
    if not os.path.exists("output"):
        os.mkdir("output")

    # Translation
    print("translate")
    for i in range(N):
        _mol = translate(mol)
        # print(CalcRMS(_mol, mol))
        # w = Chem.SDWriter(f"output/{i}/translate.sdf")
        # w.write(_mol)
        # w.close()
