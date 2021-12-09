import argparse
import time
import typing as tp
import traceback
import os
from io import BytesIO
from asyncio import ensure_future, gather, get_event_loop, Semaphore, sleep, TimeoutError

import aiohttp
from lxml import etree, html

BASE_URL = "http://www.pdbbind.org.cn/quickpdb.php?quickpdb={}"
SEMA_LIMIT = 20
TCPCONNECTOR_LIMIT = 50


async def fetch(session, pdb, semaphore):
    async with semaphore:
        try:
            async with session.get(BASE_URL.format(pdb)) as response:
                await sleep(0)
                content = await response.text(errors="ignore")
                tree = html.fromstring(content)
                try:
                    textareas = tree.xpath(".//tr/*/textarea")
                    smiles = textareas[0].text
                except ValueError or Exception:
                    smiles = "None"
                finally:
                    print(pdb, smiles)
                    return pdb, smiles
        except TimeoutError:
            print(pdb, "timeout")
            return pdb, "None"
        except Exception:
            traceback.print_exc()
            print(pdb, "exception")
            return pdb, "exception"


async def main():
    semaphore = Semaphore(SEMA_LIMIT)
    connector = aiohttp.TCPConnector(limit=TCPCONNECTOR_LIMIT)
    client = aiohttp.ClientSession(connector=connector)
    futures = []
    async with client as session:
        for pdb in pdbs:
            try:
                futures.append(ensure_future(fetch(session, pdb, semaphore)))
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
        # futures = [ensure_future(fetch(session, pdb, semaphore)) for pdb in pdbs]
        res = await gather(*futures)
        return res


def read_file(file):
    with open(file, "r") as f:
        data = f.readline()
        data = data.split(",")
    return data


def write_file(data, file, mode="w"):
    failed = []
    with open(file, mode) as w:
        for result in data:
            if "None" == result[1] or "exception" == result[1]:
                failed.append(result[0])
                continue
            w.write("\t".join(result) + "\n")
    return failed


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--file",
        help="the input file containing a comma-separated list of PDB ids",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--out_file",
        help="the output file",
        type=str,
    )
    args = parser.parse_args()

    pdbs = read_file(args.file)
    loop = get_event_loop()
    start = time.time()
    results = loop.run_until_complete(main())
    end = time.time()
    print(f"elapsed time = {end - start}s")
    loop.close()
    failed = write_file(results, args.out_file)
