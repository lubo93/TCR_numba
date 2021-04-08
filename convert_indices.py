# INFO: convert indices: 

import numpy as np
import make_net 

import argparse
import sys

if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', help="name of the single-number-index file, default data/sparse.txt")
    parser.add_argument('--N_part', help="number of sequences per single block, \
            default 4*10**3, better don't change that")
    parser.add_argument('--len_xy', help="number of single blocks in x and y direction respectively. \
            number of processed (not loaded) sequences: len_xy * N_part")
    parser.add_argument('--target_name', help="target name, under which the new list is saved")
    args = parser.parse_args()
    
    name = args.name if not args.name == None else "data/sparse.txt"
    len_x = args.len_xy if not args.name == None else 10
    single_sidelength = args.len_xy if not args.len_xy == None else 4*10**3
    target_name = args.target_name if not args.target_name == None else "two_number_indices.txt"

    data = make_net.convert_edges(name=name, len_x=len_x, single_sidelength=single_sidelength)
    with open(target_name, "w") as f: 
        np.savetxt(f, data)
    print("\n saved under \"", target_name, "\".")

