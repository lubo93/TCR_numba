# try to look, when a common total cluster is created

import funcDictionary as dic
import funcDictionaryLevenshtein as dicLeven
import plotData

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

import time

import networkx as nx
import sys

import imnet

import matplotlib.patches as mpatches
#from matplotlib.pypot import gca, show

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from collections import OrderedDict

sys.path.append("./phasetransition/")
import plot1 as p1
import plot2 as p2
import plot3 as p3
import plot4 as p4

import sys
import argparse

import graphics 
# todo: make an overall class plot

from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
# info: example command from the command line:
#     python performance.py --name="data/spar"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", help="name under which the edge list is saved, usually: data/sparse_performance.txt")
    parser.add_argument("--idx_max", help="maximum index of a pathogen sequence, usually 0")
    parser.add_argument("--name_params", help="name under which to save the parameters, usually: data/sparse_params_performance.txt")
    parser.add_argument("--parser_N_part", help="number of sequences per single block, usually 4000 ")
    parser.add_argument("--parser_len_xy", help="number of single blocks in x and y direction respectively. \
            number of processed (not loaded) sequences: len_xy * N_part")
    
    args = parser.parse_args()
    name = args.name if not args.name == None else "data/sparse_performance.txt"
    idx_max = int(args.idx_max) if not args.idx_max == None else 0
    name_params = args.name_params if not args.name_params == None else "data/sparse_params_performance.txt"
    parser_N_part = int(args.parser_N_part) if not args.parser_N_part == None else 10**3
    parser_len_xy = int(args.parser_len_xy) if not args.parser_len_xy == None else 4
    
    plotData = plotData.PlotDataLevenshtein()
    
    graphics.design(2, 0.5)
    
    # info: set parameters
    plotData.clusterMinLen = 1
    plotData.N = 1*10**3
    plotData.min_ldVal = -1
    plotData.maxVal = 3
    plotData.max_ldVal = plotData.maxVal
    gpu_l = 8000 #5000
    step = 0
    
    #Ns = np.multiply([8, 16], 10**3)
    Ns = np.multiply([0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0], 10**3)
    #Ns = np.multiply([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 4.0, 5.0], 10**3)
    
    #maxVals = [1, 2]
    maxVals = [1, 2, 3, 4, 5]
    
    xlabel = 0.02
    ylabel = 0.92
    
    DxLabel = 0.115
    DyLabel = 0.1
    
    # info: plotting the phase transition -> over N for fixed max_ldVal
    
    #fig, ax = plt.subplots(nrows = 1, ncols = 2) - info: two plots - only for masterthesis
    fig, ax = plt.subplots(nrows = 1, ncols = 2)
    
    #     but not in the paper
    #ax, ax_inset = p4.plot4(ax, ax_inset)
    ax[0], ax[1] = p4.plot4(ax[0], ax[1], name=name, idx_max=idx_max, name_params=name_params, parser_N_part=parser_N_part, parser_len_xy=parser_len_xy)
    
    #-ax[1][1] = p4.plot4(ax[1][1])
    
    print("\n /home/paul/Documents/imnet-master2/phasetransition.py: now going to show")
    plt.show()    
    
    fig.tight_layout()
    fig.savefig("performance.png", dpi=300, bbox_inches="tight")
    print("\n /home/paul/Documents/imnet-master2/phasetransition.py: saved")
