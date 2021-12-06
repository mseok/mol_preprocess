import argparse
import glob
import os
from subprocess import run


def main(args: argparse.Namespace) -> None:
    files = os.path.join(args.pdbbind_dir, "????/????_protein.pdb")
    pdbs = glob.glob(files)
    for pdb in pdbs:
        if not os.path.exists(f"{pdb.split('.')[0]}_nowater.pdb"):
            cmd = f"grep -v HOH {pdb}".split()
            with open(f"{pdb.split('.')[0]}_nowater.pdb", "w"):
                run(cmd, stdout=w)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--pdbbind_dir", type=str)
    args = parser.parse_args()

    main(args)
