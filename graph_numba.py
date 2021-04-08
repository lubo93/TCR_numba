from numba import jit, njit, vectorize
import time
import matplotlib.pyplot as plt
import numpy as np

import sys
import argparse 

from concurrent.futures import ThreadPoolExecutor

from numba import cuda

import funcDictionary as dic
#import funcDictionaryLevenshtein as dicLeven
import plotData
import textdistance as td

import psutil
from scipy import sparse
import pandas as pd
import math

from graph_support import *

def gpu_sparse(strings1, strings2, num_strings, blockspergrid, threadsperblock):
    """
    Info: calculate the sparse Levenshtein distance matrix of 
        two arrays "strings1" and "strings2"  on the GPU
    Args: strings1: first array of strings, encoded in ASCII numbers, 
            length_of_array = num_strings * length_per_string, with 
            length_per_string = 30
        strings2: first array of strings, encoded in ASCII numbers
        num_strings: number of strings in each array, usually 4*10**3, 
        blockspergrid: a parameter for the GPU calculation
        threadsperblock: a parameter for the GPU calculation
    Returns: h_sparse: list of indices, where the matrix entry, is 1. 
            The entry is 1, if and only if the Levenshtein distance between 
            the two corresponding strings is smaller than or equal to the 
            threshold distance max_ld = 2.
    """
    # TODO: check, if our condition is <2 or <=2
    strings1 = strings1.astype(dtype=np.int8)
    strings2 = strings2.astype(dtype=np.int8)

    dist = np.zeros(num_strings**2)#**2
    dist = dist.astype(dtype=np.int8)

    d_strings1 = cuda.to_device(strings1)
    d_strings2 = cuda.to_device(strings2)

    d_num_strings = cuda.to_device(num_strings)
    d_dist = cuda.to_device(dist)

    max_ld = 2
    d_max_ld = cuda.to_device(max_ld)
    
    cuda.synchronize()
    levenshtein[blockspergrid, threadsperblock](d_strings1, d_strings2, num_strings, d_dist, max_ld)
    cuda.synchronize()
    h_dist = d_dist.copy_to_host()
    h_sparse = np.asarray(np.where(h_dist!=0), dtype=np.int32)
    h_sparse = h_sparse[0]

    return h_sparse

def total_idx(el, N_part, i, j, len_x, len_y):
    """
    Info: calculate the total one-number-index from the one-number-index
        in the small block matrix with sidelength N_part
    Args: el: one-number-index in the long block
        N_part: sidelength of single block (usually 4*10**3) 
        i: row index of the block
        j: column index of the block
        len_x: number of blocks in one row
        len_y: number of blocks in one column
    Returns: idx: Total index, one number
    """
    offset_index = len_x * i * N_part**2
    el_x = int(el/N_part)
    el_y = int(el%N_part)
    el_A = el_x * N_part * len_y
    el_B = el_y + j*N_part
    idx = offset_index + el_A + el_B
    return idx

def adjacency_matrix(seq, name="sparse.txt", idx_max=0, name_params="data/sparse_params.txt", N_part=4*10**3, len_xy=10):
    """
    Info: calculating the entire adjacency matrix and save it
    Args: seq: the array of strings
        name: name of the txt file to save the one-number-index
        idx_max: number of sequences, which are replaced by virus 
            recognizing sequences (0, except if you do virus analysis)
        name_params: name of the txt file to save idx_max.
        N_part: sidelength of a single block, default: 4*10**3
        len_xy: number of single blocks in x and y direction respectively, default: 10
    Returns: -
    """

    if not idx_max == 0:
        with open(name_params, 'w') as f:
            np.savetxt(f, np.array([idx_max]), fmt='%i')

    Ns = np.multiply([0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0], 10**3)
    maxVals = [1, 2, 3, 4, 5]
    #N_part = 4*10**3
    num_strings = N_part
    threadsperblock = 256# 24
    blockspergrid = ((num_strings)**2 + (threadsperblock - 1)) // threadsperblock

    print("\n starting ascii ...")
    strings_1, strings_2 = convert_to_ascii(seq, seq, 30)
    print("\n seq: ", len(seq), "; len(strings_1): ", len(strings_1), "; len(strings_2): ", len(strings_2))
    print("\n ")
    print("\n finished ascii")

    #@cuda.jit('void(int8[:], int8, int8[:])')
    tA = time.time()
    
    # TODO: make sure, you always iterate over the entire length of seq
    with open(name, 'w') as f:
        len_x, len_y = len_xy, len_xy
        for i in range(len_x):#250
            print("\n now setting to zero ...")
            h_total = np.array([])
            print("\n ... finished setting to zero")
            for j in range(i, len_y):#250
                strings1 = strings_1[(i*N_part*30):((i+1)*N_part*30)]
                strings2 = strings_2[(j*N_part*30):((j+1)*N_part*30)]
                
                h_sparse = gpu_sparse(strings1, strings2, num_strings, blockspergrid, threadsperblock)
                print("\n np.shape(h_sparse): ", np.shape(h_sparse), "; np.shape(h_total): ", np.shape(h_total))

                # info: addding the offset index of the block
                #offset_index = (len_y * i + j) * N_part**2
                #h_sparse = [el + offset_index for el in h_sparse]
                h_sparse = [total_idx(el, N_part, i, j, len_x, len_y) for el in h_sparse]
                
                # info: removing self edges
                h_sparse = contract_array(h_sparse, idx_max, N_part*len_y)

                h_total = np.concatenate([h_total, h_sparse])
                print("\n i: ", i, "j: ", j)
                
            s = "\n " + str(h_total)
            print("\n h_total[0]: ", h_total[0], "; h_total[1]: ", \
                    h_total[1], "; el%N_part: ", h_total[0]%N_part, \
                    "; int(el/N_part): ", int(h_total[1]/N_part))
            #-h_total = [[el%(N_part*len_x), int(el/(N_part*len_x))] for el in h_total]
            #f.write(s)
            np.savetxt(f, h_total, fmt='%i')
    tB = time.time()

    print("\n dt: ", tB - tA) 
    print("\n blockspergrid: ", blockspergrid)

if __name__ == "__main__":
    fig1, ax1 = plt.subplots()

    #Ns = [1*10**1, 2*10**1]#, 10**4, 10**5]
    Ns = [10**1]
    t_py = []
    t_gpu = []
    
    # TODO: establish command line argument taking
    parser = argparse.ArgumentParser()
    parser.add_argument('--N', help="number of sequences to be loaded, default 10**6")
    parser.add_argument('--N_part', help="number of sequences per single block, \
            default 4*10**3, better don't change that")
    parser.add_argument('--len_xy', help="number of single blocks in x and y direction respectively. \
            number of processed (not loaded) sequences: len_xy * N_part")
    parser.add_argument('--name', help="name of the single index edge list file")
    parser.add_argument('--name_params', help="name of the name_params file")
    parser.add_argument('--with_pathogens', help="1: with pathogens, 0 (default): without pathogens")
    # TODO: make sure that the path is also set correctly when 
    #     with_pathogens=1, currently he just saves under the 
    #     default name and path
    args = parser.parse_args()

    parser_N = int(args.N) if not args.N == None else 10**6
    parser_N_part = int(args.N_part) if not args.N_part == None else 4*10**3
    parser_len_xy = int(args.len_xy) if not args.len_xy == None else 10
    parser_name = args.name if not args.name == None else "data/sparse.txt"
    parser_name_params = args.name_params if not args.name_params == None else "data/sparse_params.txt"
    with_pathogens = True if args.name_params == 1 else 0
    
    for N in Ns:
        plotData.clusterMinLen = 1
        plotData.N = parser_N#10**6#4*10**3
        plotData.min_ldVal = -1
        plotData.maxVal = 3
        plotData.max_ldVal = plotData.maxVal
        gpu_l = 8000 #5000
        step = 0
        a, a4, seq, filename, ls = dic.loadSequence(step, plotData, isExtractNum=False)
        
        #with_pathogens = True
        if with_pathogens:
            pathogens_savenames = ["Influenza", "HIV", "Alzheimer", \
                    "Parkinson", "Tuberculosis", "COVID-19"]
            pathogens = ["Influenza", "Human immunodeficiency virus (HIV)", \
                    "Alzheimer's disease", "Parkinson disease", \
                    "M.Tuberculosis", "COVID-19"]
            for pathogen, save_name in zip(pathogens, pathogens_savenames):
                index = 0
                print("\n ------------------- pathogen index: ", pathogens.index(pathogen))
                seq, idx_max = insert_pathogens(seq, pathogen)# if with_pathogens == True else seq   

                name = "data/" + save_name + "_small.txt"
                name_params = "data/" + save_name + "_params.txt"
                #-print("\n seq[:500]: ", seq[:500])
                adjacency_matrix(seq, name=name, idx_max=idx_max, name_params=name_params, N_part=parser_N_part, len_xy=parser_len_xy)
        else:
            name = parser_name#"data/sparse.txt"
            adjacency_matrix(seq, name=name, N_part=parser_N_part, len_xy=parser_len_xy)
            print("\n saved under name: ", name)
