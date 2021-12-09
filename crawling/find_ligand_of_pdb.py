import argparse
import glob
import os


def get_ligand_ids(file: str):
    with open(file, "r") as f:
        lines = f.readlines()
        lines = [line.split() for line in lines]
        lines = [line[1] for line in lines if "HET" == line[0]]
    ligand_ids = list(set(lines))
    return ligand_ids


def get_all_ligand_ids(_dir: str):
    dic = {}
    files = glob.glob(os.path.join(_dir, "*.pdb"))
    for file in files:
        ligand_ids = get_ligand_ids(file)
        dic[file.split("/")[-1].split(".")[0]] = ligand_ids
    return dic


def write_result(dic: dict, file: str):
    with open(file, "w") as w:
        sorted_dic = sorted(list(dic.items()))
        for key, values in sorted_dic:
            values = "\t".join(values)
            w.write(f"{key}\t{values}\n")


def main(args: argparse.Namespace):
    ligands = get_ligand_ids(args.dir)
    write_result(ligands, args.out_file)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", type=str, help="input data dir")
    parser.add_argument("-o", "--out_file", type=str, help="output file")
    args = parser.parse_args()

    main(args)
