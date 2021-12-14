#!/usr/bin/env python
import argparse
import glob
from collections import defaultdict
from itertools import product


def main(args: argparse.Namespace) -> None:
    if args.pdbbind_dir:
        pdbbind_dir = os.path.abspath(args.pdbbind_dir)
        test_keys = glob.glob(f"{pdbbind_dir}/????")
        test_keys = [test_key.split("/")[-1] for test_key in test_keys]

    with open(args.input_file, "r") as f:
        lines = f.readlines()

    cluster_dict = defaultdict(list)
    for line in lines:
        if ">" == line[0]:
            key = line.split()[-1]
        else:
            cluster_dict[key].append(line.split()[3].split("/")[0])

    if args.filter_valid:
        valid_train_keys = []
        valid_test_keys = []
        for i in cluster_dict:
            flag = 0
            clusters = list(set(cluster_dict[i]))
            test_clusters = list(set(clusters) & set(test_keys))
            if test_clusters:
                valid_test_keys.append(test_clusters[0])
            else:
                valid_train_keys.append(clusters[0])
        with open("train_" + args.output_file, "w") as w:
            for key in valid_train_keys:
                w.write(key + "\n")
        with open("test_" + args.output_file, "w") as w:
            for key in valid_test_keys:
                w.write(key + "\n")
    else:
        invalid_keys = []
        for i in cluster_dict:
            flag = 0
            clusters = list(set(cluster_dict[i]))
            ith_invalid_keys = list(product(clusters, clusters))
            ith_invalid_keys = [
                ith_invalid_key
                for ith_invalid_key in ith_invalid_keys
                if ith_invalid_key[0] != ith_invalid_key[1]
            ]
            ith_invalid_keys = [
                "_".join(ith_invalid_key) for ith_invalid_key in ith_invalid_keys
            ]
            invalid_keys += ith_invalid_keys
        with open(args.output_file, "w") as w:
            for key in invalid_keys:
                w.write(key + "\n")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str)
    parser.add_argument("-o", "--output_file", type=str)
    parser.add_argument("-d", "--pdbbind_dir", type=str)
    parser.add_argument("--filter_valid", action="store_true")
    args = parser.parse_args()

    main(args)
