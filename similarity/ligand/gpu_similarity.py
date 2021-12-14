import argparse
import glob
import os
import sys

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader


def cal_similarity(vec1, vec2):
    # [N1, 1024], [N2, 1024] -> [N2, N1]
    or1 = torch.logical_not(vec1).type(torch.float16).to(device)
    or2 = torch.logical_not(vec2).type(torch.float16).to(device)
    and_out = torch.matmul(vec1, vec2.T)
    or_out = torch.matmul(or1, or2.T)
    or_out = N_FEATURES - or_out
    similarity = and_out / or_out
    return similarity


def get_dataloader(dataset, batch_size):
    return DataLoader(dataset, batch_size=batch_size, shuffle=False)


class BitData(Dataset):
    def __init__(self, fn):
        self.bits = self._load(fn)

    def __len__(self):
        return len(self.bits)

    def __getitem__(self, idx):
        return self.bits[idx]

    def _load(self, fn):
        bits = np.load(fn)
        return bits


def main(dataset1, dataset2, batch1, batch2, self_similarity=False):
    indices = []
    accum_idx = 0
    dataloader1 = iter(get_dataloader(dataset1, batch1))
    for idx1, data1 in enumerate(dataloader1):
        dataloader2 = iter(get_dataloader(dataset2, batch2))
        data1 = data1.type(torch.float16).to(device)
        sim_tmp = torch.zeros(data1.shape[0], batch2).to(device)
        for idx2, data2 in enumerate(dataloader2):
            data2 = data2.type(torch.float16).to(device)
            sim = cal_similarity(data1, data2)
            if self_similarity and idx1 == idx2:
                diagonal = torch.eye(data1.shape[0])
                diagonal = diagonal.type(torch.float16).to(device)
                sim -= diagonal
            val_sim_tmp = torch.where(sim > THRESHOLD, 1, 0)
            # print(val_sim_tmp)
            # print(idx1, idx2, data2.shape, sim.shape, val_sim_tmp.shape)
            if sim_tmp.shape != val_sim_tmp.shape:
                sim_tmp = torch.cat((sim_tmp, val_sim_tmp), dim=1)
            else:
                sim_tmp = torch.logical_or(sim_tmp, val_sim_tmp)
            # print(sim_tmp.shape)
        sim_idx = torch.sum(sim_tmp, dim=1).cpu()
        sim_idx = torch.where(sim_idx == 0)[0]
        sim_idx += accum_idx
        accum_idx += len(data1)
        indices.append(sim_idx)
        # print(sim_idx.shape)
    indices = torch.cat(indices, dim=0)
    # print(len(indices))
    return indices


def check_current_running_device(ngpu):
    from subprocess import Popen
    from subprocess import PIPE

    for i in range(ngpu):
        nvidia_smi = Popen(["nvidia-smi", "-i", str(i)], stdout=PIPE)
        grep = Popen(["grep", "No running"], stdin=nvidia_smi.stdout, stdout=PIPE)
        count = Popen(["wc", "-l"], stdin=grep.stdout, stdout=PIPE, encoding="utf-8")
        output = int(count.communicate()[0].split("\n")[0])
        if output:
            break
    return i


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=str, help="compare filename")
    parser.add_argument("--ref", type=str, help="refence filename")
    parser.add_argument("--ngpus", type=int, help="number of gpus in node")
    args, _ = parser.parse_known_args()

    PATH = "similarity"
    if not os.path.exists(PATH):
        os.mkdir(PATH)
    device_id = check_current_running_device(args.ngpus)
    os.environ["CUDA_VISIBLE_DEVICES"] = str(device_id)
    device = torch.device("cuda")
    THRESHOLD = 0.7
    N_FEATURES = 1024
    save_fn = os.path.join(PATH, args.target.split("/")[-1])
    dataset1 = BitData(args.target)
    dataset2 = BitData(args.ref)
    print("similsrity between train dataset")
    print(args.target, "length: ", len(dataset1))
    print(args.ref, "length: ", len(dataset2))

    total1 = torch.zeros(len(dataset1))
    indices1 = main(dataset1, dataset2, 10000, 20000).long()
    total1[indices1] = 1

    print("self similarity")
    total2 = torch.zeros(len(dataset1))
    indices2 = main(dataset1, dataset1, 10000, 10000, True).long()
    total2[indices2] = 1

    total = torch.logical_and(total1, total2)
    total = torch.where(total == 1)[0].numpy()
    np.save(save_fn, total)
