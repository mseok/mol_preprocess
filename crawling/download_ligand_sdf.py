import argparse
import glob
import os
import time
import typing as tp
import traceback
from io import BytesIO
from asyncio import (
    ensure_future,
    gather,
    get_event_loop,
    Semaphore,
    sleep,
    TimeoutError,
)

import aiohttp
from lxml import etree, html

BASE_URL = "http://www.pdbbind.org.cn/quickpdb.php?quickpdb={}"
SEMA_LIMIT = 50
TCPCONNECTOR_LIMIT = 50

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
    "asym_id": "label_asym_id",
    "pdb_seq_num": "auth_seq_id",
    "pdb_mon_id": "ligand_key",
}
KEYWORD = "_pdbx_nonpoly_scheme.pdb_ins_code"


def get_base_url() -> str:
    base_url = "https://models.rcsb.org/v1/{}/ligand?"
    data = [
        "label_asym_id={}",
        "auth_seq_id={}",
        "encoding=sdf",
        "copy_all_categories=true",
        "filename={}",
    ]
    data = "&".join(data)
    # base_url = base_url + data
    base_url = base_url + "{}"
    return base_url


async def fetch(session, semaphore, data):
    async with semaphore:
        data: dict
        pdbid = data.pop("pdbid")
        filename = data["filename"]
        data = "&".join(["=".join(item) for item in data.items()])
        url = BASE_URL.format(pdbid, data)
        try:
            async with session.get(url) as response:
                await sleep(0)
                content = await response.text(errors="ignore")
                with open(os.path.join(args.out_dir, filename), "w") as w:
                    w.write(content)
        except TimeoutError:
            print(filename, "timeout")
            return filename, "None"
        except Exception:
            traceback.print_exc()
            print(filename, "exception")
            return filename, "exception"


async def main():
    semaphore = Semaphore(SEMA_LIMIT)
    connector = aiohttp.TCPConnector(limit=TCPCONNECTOR_LIMIT)
    client = aiohttp.ClientSession(connector=connector)
    futures = []
    async with client as session:
        for file in files:
            data = read_cif(file)
            try:
                if isinstance(data, dict):
                    futures.append(ensure_future(fetch(session, semaphore, data)))
                elif isinstance(data, list):
                    for _data in data:
                        futures.append(ensure_future(fetch(session, semaphore, _data)))
            except aiohttp.ClientConnectionError:
                # something went wrong with the exception, decide on what to do next
                client.close()
                semaphore.release()
                client = aiohttp.ClientSession(connector=connector)
            except aiohttp.ClientError:
                # something went wrong in general. Not a connection error, that was handled
                # above.
                client.close()
                semaphore.release()
                client = aiohttp.ClientSession(connector=connector)
        res = await gather(*futures)
        return res


def read_cif(file: str):
    with open(file, "r") as f:
        lines = f.readlines()
        lines = [line for line in lines]
        start = [line for line in lines if KEYWORD in line][0]
        idx = lines.index(start) + 1
        data = lines[idx].split()
        length = len(data)
        if length < len(PDBX_NONPOLY_SCHEME):
            start = idx - 10
            end = idx
            lines = lines[start:end]
            data = [line.split()[1] for line in lines]
            data = dict(zip(PDBX_NONPOLY_SCHEME, data))
            required_data = refine_data(data, file)
            return required_data
        elif length == len(PDBX_NONPOLY_SCHEME):
            lines = lines[idx:]
            required_datas = []
            for line in lines:
                if "HOH" in line:
                    continue
                data = line.split()
                if data[0] == "#":
                    break
                data = dict(zip(PDBX_NONPOLY_SCHEME, data))
                required_data = refine_data(data, file)
                required_datas.append(required_data)
            return required_datas
        else:
            raise Exception("file: ", file)


def refine_data(data: tp.Dict, file: str):
    required_data = dict()
    for key, new_key in KEY_DICT.items():
        required_data[new_key] = data[key]
    pdbid = file.split("/")[-1].split(".")[0]
    required_data["pdbid"] = pdbid
    required_data["encoding"] = "sdf"
    label_asym_id = required_data["label_asym_id"]
    ligand_key = required_data["ligand_key"]
    required_data["filename"] = f"{pdbid}_{label_asym_id}_{ligand_key}.sdf"
    return required_data


if __name__ == "__main__":
    BASE_URL = get_base_url()
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cif_dir", type=str, default="./rcsb_cifs/")
    parser.add_argument("-o", "--out_dir", type=str, default="./rcsb_sdfs/")
    args = parser.parse_args()

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    files = glob.glob(os.path.join(args.cif_dir, "*.cif"))
    loop = get_event_loop()
    start = time.time()
    results = loop.run_until_complete(main())
    end = time.time()
    print(f"elapsed time = {end - start}s")
    loop.close()
