import textdistance as td
import funcDictionary as dic
import numpy as np

import sys
import argparse

import plotData

def ref_levenshtein(parser_N, src_name, ld_max):
    """
    Info: calculate the adjacency matrix with the normal textdistance
        Levenshtein distance to verify if our algorithm is correct included
    Args:  ii
    Returns: -
    """

    plotData.clusterMinLen = 1
    plotData.N = parser_N#10**6#4*10**3
    plotData.min_ldVal = -1
    plotData.maxVal = 3
    plotData.max_ldVal = plotData.maxVal
    gpu_l = 8000 #5000
    step = 0
    _, _, seq, _, _ = dic.loadSequence(step, plotData, isExtractNum=False, src=src_name)

    a = np.zeros([len(seq), len(seq)])

    for i in range(len(seq)):
        for j in range(len(seq)):
            dist = td.levenshtein(seq[i], seq[j])
            a[i][j] = 1 if dist < ld_max else 0

    name = "two_number_indices_ref.txt"
    with open(name, 'w') as f:
        np.savetxt(f, a, fmt='%i')
    print("\n saved under: ", name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--N', help="")
    parser.add_argument('--src', help="")
    parser.add_argument('--ld_max', help="<ld_max = 1 else 0")
    args = parser.parse_args()

    parser_N = int(args.N)
    src_name = args.src
    ld_max = int(args.ld_max) if not args.ld_max == None else 2
    ref_levenshtein(parser_N, src_name, ld_max)
