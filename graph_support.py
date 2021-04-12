# info: contains the functions for graph_numba.py

from numba import jit, njit, vectorize
import time
import matplotlib.pyplot as plt
import numpy as np

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

def isNaN(string):
    return string!=string

def convert_to_ascii(strings1, strings2, max_len):
    """
    info: take two arrays of strings and convert both into 
        the special ascii format, where each character
        is saved as an ubyte (unsigned byte), which is 8 
        bits.
        Actually ASCII would only need 7 bits, so there would
        still be little space for improvement.
        It is assumed, that no string is longer than max_len, which
        is usually 30. If a string is shorter, then "0"-characters
        will be appended until the length is max_len.
        Thus, all strings have the same length, which is 
        necessary in order to define equidistant buffers later
        for the GPU.
        (I think also there is a way to improve that.)
    inputs: strings1: array of strings
        strings2: array of strings
        max_len: maximum length, that each string can have 
        (usually 30)
    outputs: ascii_total1, ascii_total2:
            two string arrays converted to the special ascii_type, where
            each element has the same length max_len
    """

    #ascii_total = np.array([])
    ascii_total = list()
    #strings = strings.astype(np.float32)

    for i in range(len(strings1)):
        s = strings1[i]
        if len(s) > max_len:
            print("\n info: an s is longer then max_len, we cut it")
            s = s[:max_len]
        for j in range(len(s), max_len):
            s = s + '0'
        ascii_code = [ord(c) for c in s]
        for j in range(len(ascii_code)):
            ascii_total.append(ascii_code[j])
        #ascii_total = ascii_total + ascii_code
    ascii_total = np.array(ascii_total)
    ascii_total1 = ascii_total.astype(np.ubyte)

    ascii_total = list()
    #strings = strings.astype(np.float32)

    for i in range(len(strings2)):
        s = strings2[i]
        if len(s) > max_len:
            print("\n info: an s is longer then max_len, we cut it")
            s = s[:max_len]
        for j in range(len(s), max_len):
            s = s + '0'
        ascii_code = [ord(c) for c in s]
        for j in range(len(ascii_code)):
            ascii_total.append(ascii_code[j])
        #ascii_total = ascii_total + ascii_code
    ascii_total = np.array(ascii_total)
    ascii_total2 = ascii_total.astype(np.ubyte)

    return ascii_total1, ascii_total2

@cuda.jit('void(int8[:], int8[:], int32, int8[:], int8)')
def levenshtein(strings1, strings2, num_strings, dist, max_ld):
    """
    Info: GPU kernel - computes the Levenshtein distance based adjacency 
        matrix for the string lists strings1 and strings2
    Args: strings1: list of strings, where each string has 30 letters 
            in ascii format
        strings2: list of strings, where each string has 30 letters 
            in ascii format
        num_strings: number of strings in each of the arrays strings1 and strings2
            (hence we assume, that strings1 and strings2 have the same length)
        dist: empty buffer for the adjacency matrix
        max_ld: threshold distance
    Returns: -
    """ 

    test = 0

    tx = cuda.threadIdx.x
    ty = cuda.blockIdx.x
    bw = cuda.blockDim.x
    gid = ty*bw + tx
    
    lx = int(gid%num_strings)
    ly = int(gid/num_strings)

    m = 30
    n = 30
    m_done = False
    n_done = False
    
    dprev = cuda.local.array(shape=(30), dtype=np.int8)
    dnew = cuda.local.array(shape=(30), dtype=np.int8)
    els = cuda.local.array(shape=(3), dtype=np.int8)

    dp0, dp1, dp2, dp3, dp4, dp5, dp6, dp7, dp8, dp9, dp10, dp11, dp12, dp13, dp14, dp15, dp16, dp17, dp18, \
            dp19, dp20, dp21, dp22, dp23, dp24, dp25, dp26, dp27, dp28, dp29 = 0, 1, 2, 3, 4, 5, 6, 7, 8, \
            9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29

    dn0, dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8, dn9, dn10, dn11, dn12, dn13, dn14, dn15, dn16, dn17, dn18, \
            dn19, dn20, dn21, dn22, dn23, dn24, dn25, dn26, dn27, dn28, dn29 = 0, 1, 2, 3, 4, 5, 6, 7, 8, \
            9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29
    
    dif = cuda.local.array(shape=1, dtype=np.int8)
    strings1_local = cuda.local.array(shape=30, dtype=np.int8)
    strings2_local = cuda.local.array(shape=30, dtype=np.int8)

    for i in range(30):
        strings1_local[i] = strings1[lx*30 + i]
        strings2_local[i] = strings2[ly*30 + i]
        if strings1[lx*30 + i] == 48:
            if m_done == False:
                m = i
                m_done = True
            else:
                pass
        if strings2[ly*30 + i] == 48:
            if n_done == False:
                n = i
                n_done = True
            else:
                pass
    
    for j in range(m+1):
        dprev[j] = j
    d_val = -1
    
    min_val = 0

    for i in range(1, n+1):#n+1):
        dnew[0] = i
        dn0 = i

        min_val = 30
        for j in range(1, m+1):
            if strings2_local[i-1] == strings1_local[j-1]:
                substitutionCost = 0
            else:
                substitutionCost = 1
            """
            el1 = dprev[j]# +1
            el2 = dnew[j-1]# +1
            el3 = dprev[j-1] + substitutionCost
            
            if el2 < el1:
                el1 = el2
            if el3 < el1 + 1:
                el1 = el3
            else: 
                el1 += 1
            dnew[j] = el1
            """

            els[0] = dprev[j]# +1
            els[1] = dnew[j-1]# +1
            els[2] = dprev[j-1] + substitutionCost

            if els[1] < els[0]:
                els[0] = els[1]
            if els[2] < els[0] + 1:
                els[0] = els[2]
            else:
                els[0] += 1
            dnew[j] = els[0]
            
            
            if dnew[j] < min_val:
                min_val = dnew[j]
            
            #d_val = el1
        
        dprev[0] = i 
        
        for j in range(1, m+1):
            dprev[j] = dnew[j]
        
        if min_val > max_ld:
            break
        #print("\n A gid: ", gid, "; i: ", i, "; j: ", j, "; in min_val")

    #dist[gid] = d_val

    dist[gid] = 1 if dnew[m]<max_ld else 0

def insert_pathogens(seq, pathogen_name):
    """
    Info: replace the first strings in "seq" by the pathogen cognate CDR3 
        of the pathogen with name "pathogen_name", where the corresponding 
        cognate CDR3 are read from the "virus_databank.csv" file
    Args: seq: normal array of strings of CDR3 sequences (normal format)
        pathogen_name: name of the pathogen
    Returns: seq: new sequence with the replaced CDR3
        idx: number of replaced sequences
    """
   
    data = pd.read_csv('virus_databank.csv', low_memory=False)
    CDR3_beta_total = list(data["CDR3.beta.aa"])
    CDR3_beta_names = list(data["Pathology"])
    
    idx = 0
    for i in range(len(CDR3_beta_total)):
        #print("\n i: ", i, "; CDR3_beta_names[i]: ", CDR3_beta_names[i], "; pathogen_name: ", pathogen_name)
        if CDR3_beta_names[i] == pathogen_name:
            #print("\n has a chance at i=", i)
            #print("\n i: ", i, "; CDR3_beta_names[i]: ", CDR3_beta_names[i], "; pathogen_name: ", pathogen_name)
            #-try:
            if True: 
                #print("\n CDR3_beta_total[i]: ", CDR3_beta_total[i], "; isnan: ", isNaN(CDR3_beta_total[i]))
                #if not math.isnan(CDR3_beta_total[i]):
                if not isNaN(CDR3_beta_total[i]):
                    seq[idx] = CDR3_beta_total[i]
                    idx += 1
            #-except: 
            #-    pass
    # info: increasing by one in order to return a length value, 
    #    which is usually 1 higher than the highest index 
    idx += 1
    return seq, idx

def contract_array(h_total, idx_max, side_len):
    """
    Info: remove the edges, which just connect two antigen cognate CDR3
    Args: h_total: sparse adjacency matrix
        idx_max: number of sequences which are replaced by antigen cognate 
            CDR3 sequence
        side_len: length of string arrays, whose distance matrix is calculated
            (i.e., side length of distance matrix)
    Returns: new sparse adjacency matrix with removed strings
    """
    # perh. TODO: change the name to something more fitting, e.g. "remove self_edges"
    
    i = 0
    while i < len(h_total):
        x = h_total[i]%side_len
        y = int(h_total[i]/side_len)
        if (x < idx_max and y < idx_max):
            h_total = np.delete(h_total, i)
            i-=1
            #-print("\n the was a self edge among the pathogens")
        i+=1
    return h_total
