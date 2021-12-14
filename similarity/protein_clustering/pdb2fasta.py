#!/usr/bin/env python
import argparse
import glob
import os
import re
import sys
from typing import List


def main(filename: str) -> List[str]:
    pattern = "^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])"
    ca_pattern = re.compile(pattern)
    chain_dict = dict()
    chain_list = []

    fp = open(filename, "r")
    i = 0
    for line in fp.read().splitlines():
        if line.startswith("ENDMDL"):
            break
        match_list = ca_pattern.findall(line)
        if match_list:
            resn = match_list[0][0] + match_list[0][2]
            # there can be multiple chains
            chain = match_list[0][1] + match_list[0][3]
            if chain in chain_dict:
                if resn not in AA3TO1:
                    resn = "ANY"
                chain_dict[chain] += AA3TO1[resn]
            else:
                if resn not in AA3TO1:
                    resn = "ANY"
                chain_dict[chain] = AA3TO1[resn]
                chain_list.append(chain)
    fp.close()

    lines = []
    fn = "/".join(filename.split("/")[-2:])
    for chain in chain_list:
        lines.append("> %s:%s\n%s\n" % (fn, chain, chain_dict[chain]))
    return lines


if __name__ == "__main__":
    AA3TO1 = {
        "ALA": "A",
        "VAL": "V",
        "PHE": "F",
        "PRO": "P",
        "MET": "M",
        "ILE": "I",
        "LEU": "L",
        "ASP": "D",
        "GLU": "E",
        "LYS": "K",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "TYR": "Y",
        "HIS": "H",
        "CYS": "C",
        "ASN": "N",
        "GLN": "Q",
        "TRP": "W",
        "GLY": "G",
        "MSE": "M",
        "ANY": "X",
    }

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str)
    parser.add_argument("-o", "--output_file", type=str)
    parser.add_argument("-d", "--pdbbind_dir", type=str)
    args = parser.parse_args()

    if args.output_file:
        with open(args.input_file, "r") as f:
            keys = f.readlines()
            keys = [key.split()[0] for key in keys]
        pdbbind_dir = os.path.abspath(args.pdbbind_dir)
        targets = [
            f"{pdbbind_dir}/{key}/{key}_protein_nowater.pdb" for key in keys
        ]
        fasta_list = []
        for target in targets:
            fasta = main(target)
            fasta_list.append(fasta)
        with open(args.output_file, "w") as w:
            for fasta in fasta_list:
                w.write("".join(fasta))
                if fasta_list.index(fasta) != len(fasta_list) - 1:
                    w.write("\n")
