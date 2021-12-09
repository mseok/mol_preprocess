import argparse
import glob
import os
from subprocess import run, DEVNULL
from typing import Dict

PDBX_NONPOLY_SCHEME = [
    "asym_id",
    "entity_id",
    "mon_id",
    "ndb_seq_num",
    "pdb_seq_num",
    "auth_seq_num",
    "pdb_mon_id",
    "auth_mon_id",
    "pdb_strand_id",
    "pdb_ins_code",
]
KEY_DICT = {
    "pdb_strand_id": "auth_asym_id",
    "auth_seq_num": "auth_seq_id",
    "pdb_mon_id": "ligand_key",
}
KEYWORD = "_pdbx_nonpoly_scheme.pdb_ins_code"


def get_base_url() -> str:
    base_url = "https://models.rcsb.org/v1/{}/ligand?"
    data = [
        "auth_asym_id={}",
        "auth_seq_id={}",
        "encoding=sdf",
        "copy_all_categories=true",
        "filename={}",
    ]
    data = "&".join(data)
    base_url = base_url + data
    return base_url


def read_cif(file: str):
    with open(file, "r") as f:
        lines = f.readlines()
        lines = [line for line in lines]
        start = [line for line in lines if KEYWORD in line][0]
        idx = lines.index(start) + 1
        data = dict(zip(PDBX_NONPOLY_SCHEME, lines[idx].split()))
        required_data = dict()
        for key, new_key in KEY_DICT.items():
            required_data[new_key] = data[key]
    return required_data


def download(data: Dict[str, str], pdbid: str):
    auth_asym_id = data["auth_asym_id"]
    auth_seq_id = data["auth_seq_id"]
    file = f"{pdbid}_{auth_asym_id}_{data['ligand_key']}.sdf"
    url = BASE_URL.format(pdbid, auth_asym_id, auth_seq_id, file)
    file = os.path.join(args.out_dir, file)
    cmd = ["curl", "-L", url, "-o", file]
    run(cmd, stdout=DEVNULL, stderr=DEVNULL)
    return


def main():
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    cifs = glob.glob(os.path.join(args.cif_dir, "*.cif"))
    for cif in cifs:
        pdbid = cif.split("/")[-1].split(".")[0]
        data = read_cif(cif)
        download(data, pdbid)
    return


if __name__ == "__main__":
    BASE_URL = get_base_url()
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cif_dir", type=str)
    parser.add_argument("-o", "--out_dir", type=str)
    args = parser.parse_args()

    main()
