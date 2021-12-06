import argparse
import glob


def main(args: argparse.Namespace) -> None:
    keys = glob.glob(f"{args.pdbbind_dir}/????")
    keys = [key.split("/")[-1] for key in keys]
    with open(args.output_file, "w") as w:
        for key in keys:
            w.write(key + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_file", type=str, default="target_list.txt")
    parser.add_argument("-d", "--pdbbind_dir", type=str)
    args = parser.parse_args()
