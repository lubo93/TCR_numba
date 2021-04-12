# INFO: convert indices: 

import numpy as np
import make_net 

import argparse
import sys
import scipy as sp

def dense_matrix(data, len_xy):
    """
    Info: 
    Args:
    Returns: 
    """
    
    print("\n data_in: ", data)
    rows = [el[0] for el in data]
    cols = [el[1] for el in data]
    vals = [1 for el in data]
    
    dense_matrix = sp.sparse.coo_matrix((vals, (rows, cols)), shape=(len_xy, len_xy)).toarray()
    return dense_matrix

if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', help="name of the single-number-index file, default data/sparse.txt")
    parser.add_argument('--single_sidelength', help="number of sequences per single block, \
            default 4*10**3, better don't change that")
    parser.add_argument('--len_xy', help="number of single blocks in x and y direction respectively. \
            number of processed (not loaded) sequences: len_xy * N_part")
    parser.add_argument('--target_name', help="target name, under which the new list is saved")
    parser.add_argument('--format', help="sparse: makes the two-number-indices list; \
            dense: makes the entire matrix in dense format - don't use that for very \
            large data, since it consumes high amounts of memory")
    args = parser.parse_args()
    
    name = args.name if not args.name == None else "data/sparse.txt"
    len_x = int(args.len_xy) if not args.len_xy == None else 10
    single_sidelength = int(args.single_sidelength) if not args.single_sidelength == None else 4*10**3
    target_name = args.target_name if not args.target_name == None else "two_number_indices.txt"
    is_dense = True if args.format == "dense" else False
    
    data = make_net.convert_edges(name=name, len_x=len_x, single_sidelength=single_sidelength)
    print("\n data A: ", data)
    if is_dense:
        data = dense_matrix(data, len_x*single_sidelength)

    with open(target_name, "w") as f: 
        np.savetxt(f, data, fmt='%i')
    print("\n saved under \"", target_name, "\".")

