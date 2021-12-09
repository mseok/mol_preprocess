#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from asyncio import Semaphore, ensure_future, gather, run
from io import BytesIO
from lzma import open as lzma_open
from struct import calcsize, unpack
import os

from aiofile import AIOFile
from aiohttp import ClientSession, TCPConnector

http_ok = [200]
SEMA_LIMIT = 200
CONNECTOR_LIMIT = 50
BASE_URL = "https://files.rcsb.org/download/{}.{}"


async def download():
    tasks = list()
    sem = Semaphore(SEMA_LIMIT)
    connector = TCPConnector(limit=CONNECTOR_LIMIT)
    async with ClientSession(connector=connector) as session:
        for file in files:
            tasks.append(
                ensure_future(
                    get_file(
                        file,
                        session=session,
                        sem=sem,
                    )
                )
            )
        return await gather(*tasks)


async def get_file(_id, session, sem):
    file = os.path.join(args.out_dir, f"{_id}.{args.filetype}")
    if os.path.exists(file):
        return
    url = BASE_URL.format(_id, args.filetype)
    async with sem:
        async with session.get(url) as response:
            content = await response.read()
        if response.status not in http_ok:
            print(f"Scraping {url} failed due to the return code {response.status}")
            return
        with open(file, "w") as w:
            w.write(content.decode())
        return


if __name__ == "__main__":
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--file",
        help="the input file containing a comma-separated list of PDB ids",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        help="the output dir, default: current dir",
        type=str,
        default="./",
    )
    parser.add_argument(
        "-t", "--filetype", help="file type", type=str, choices=["cif", "pdb"]
    )
    args = parser.parse_args()

    with open(args.file, "r") as f:
        files = f.readline()
        files = files.split(",")

    run(download())
