import argparse
import os


def main(args: argparse.Namespace):
    if not os.path.exists(args.smiles_file):
        URL = "http://ligand-expo.rcsb.org/ld-download.html"
        raise FileNotFoundError(f"Download smiles file from {URL}.")

    with open(args.in_file, "r") as f:
        lines = f.readlines()
        ligands = []
        lines = [line.split()[1:] for line in lines]
        for line in lines:
            ligands += line
        ligands = list(set(ligands))

    with open(args.smiles_file, "r") as f:
        lines = f.readlines()
        lines = [line.split()[:2] for line in lines]
        dic = dict([[line[1], line[0]] for line in lines])

    with open(args.in_file, "w") as w:
        print(len(ligands))
        ligands = [ligand for ligand in ligands if ligand in dic]
        print(len(ligands))
        w.write(",".join(ligands))
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_file", type=str, help="input file")
    parser.add_argument("-o", "--out_file", type=str, help="output file")
    parser.add_argument(
        "-s",
        "--smiles_file",
        type=str,
        help="smiles file",
        default="./Components-smiles-stereo-cactvs.smi",
    )
    args = parser.parse_args()

    main(args)
