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

import graph_numba as gn

# todo: make an overall class plot

plotData = plotData.PlotDataLevenshtein()

rcParams['axes.labelsize'] = 50
rcParams['xtick.labelsize'] = 46
rcParams['ytick.labelsize'] = 46
rcParams['legend.fontsize'] = 46
#rcParams['font.family'] = 'sans-serif'
rcParams['font.serif'] = ['Helvetica']
rcParams['text.usetex'] = True
rcParams['figure.figsize'] = 12, 8

# info: set parameters
plotData.clusterMinLen = 1
plotData.N = 1*10**3
plotData.min_ldVal = -1
plotData.maxVal = 3
plotData.max_ldVal = plotData.maxVal
gpu_l = 8000 #5000
step = 0

Ns = np.multiply([0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0], 10**3)
#Ns = np.multiply([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 4.0, 5.0], 10**3)
maxVals = [1, 2, 3, 4, 5]

#fig, ax = plt.subplots(nrows = 2, ncols = 2)

xlabel = 0.02
ylabel = 0.92

DxLabel = 0.115
DyLabel = 0.1

# info: plotting the phase transition -> over N for fixed max_ldVal

# info: -----------------------plot[1][1]---------------------------------------------
def plot4(ax_11, ax_2, name=None, idx_max=None, name_params=None, parser_N_part=None, parser_len_xy=None):
    # Let's do computing time analysis:
    #-Ns = [0.5*10**3, 1.0*10**3, 2.0*10**3, 3.0*10**3, 4.0*10**3, 5.0*10**3, 8.0*10**3, 16*10**3, 24*10**3, 32*10**3]
    #Ns = [0.5*10**3, 1.0*10**3, 2.0*10**3, 3*10**3, 4*10**3, 5*10**3, 8*10**3, 16*10**3, 32*10**3]
    Ns = [4*10**3, 8*10**3, 12*10**3]
    
    # info. maxVals should always only contain one array
    maxVals = [3]
    dts = list()
    dts_fast = list()
    dts_imnet = list()
    
    for i in range(len(Ns)):
        plotData.N = Ns[i]
        
        for i in range(len(maxVals)):  # this should always contain only 1 array
            print("\n paul in phasetransition.py: ", maxVals[i])
        
            plotData.max_ldVal = maxVals[i]
            print("\n gpu ...")
            t0_load = time.time()
            a, a4, seq, filename, ls = dic.loadSequence(step, plotData, isExtractNum=False)
            # info: of course t1_load should be always quite the same as t0_graph, because
            #     they are just next to each other
            t1_load = time.time()
            dt_load = t1_load - t0_load
            
            # TODO: What's a good value for parser_N_part and len_xy? Values must equal the 
            #    values of imnet in order to ensure comparability
            #-name, idx_max, name_params, parser_N_part, parser_len_xy = "dummy_file.txt", 0, \
            #-        "dummy_file_2.txt", 1000, 5
            parser_len_xy = int(plotData.N/parser_N_part)

            t0_graph = time.time()
            gn.adjacency_matrix(seq, name=name, idx_max=idx_max, name_params=name_params, N_part=parser_N_part, len_xy=parser_len_xy)
            #-g = dicLeven.make_graph(seq, min_ld=plotData.min_ldVal,
            #                        max_ld=plotData.max_ldVal, gpu_l=gpu_l, mode="slow")
            t1_graph = time.time()
            dt_graph = t1_graph - t0_graph
            dts.append(dt_graph)
            """-
            t0_graph = time.time()
            g = dicLeven.make_graph(seq, min_ld=plotData.min_ldVal,
                                    max_ld=plotData.max_ldVal, gpu_l=gpu_l, mode="fast")
            t1_graph = time.time()
            dt_graph = t1_graph - t0_graph
            dts_fast.append(dt_graph)
            """
            #-print("\n len(g.edges(): ", len(g.edges()))
            print("\n imnet ...")
            #todo: is min_ld here also min_ld< or min_ld<=? same with max_ldVal
            t0_imnet = time.time()
            g = imnet.process_strings.generate_graph(seq, min_ld=plotData.min_ldVal, 
                    max_ld=plotData.max_ldVal)
            t1_imnet = time.time()
            dt_imnet = t1_imnet - t0_imnet
            dts_imnet.append(dt_imnet)
    
    max_gpu = max(dts)
    max_imnet = max(dts_imnet)
    max_total = max([max_gpu, max_imnet])

    dts_ratio = list()
    dts_ratio_fast = list()
    for i in range(len(dts)):
        dts_ratio.append(dts[i]/dts_imnet[i])
        #-dts_ratio_fast.append(dts_fast[i]/dts_imnet[i])
    ax_2.plot(Ns, [1/el for el in dts_ratio])
    #-ax_2.plot(Ns, [1/el for el in dts_ratio_fast])
    ax_2.set_xlabel("$N$")
    #ax_2.set_ylabel("$t_{gpu}/t_{imnet}$")
    ax_2.set_ylabel("speedup")
    ax_2.legend(["gpu"], frameon=False)
    ax_2.set_ylim(bottom=0.0)
    ax_2.set_xlim(left=0.0)
    ax_11.plot(Ns, dts)
    #-ax_11.plot(Ns, dts_fast)
    ax_11.plot(Ns, dts_imnet)

    ax_11.legend(["gpu", "imnet"], frameon=False)
    ax_11.set_xlabel('$N$')
    ax_11.set_ylabel('$t$')
    ax_11.set_xlim(0, Ns[len(Ns)-1])
    ax_11.set_ylim(0, max_total)
    
    ax_11.add_patch(Polygon([[1.0, 1.0], [DxLabel, 1.0], [DxLabel, 1.0 - DyLabel], [0.0, 1.0 - DyLabel]],
        closed=False, fill=True, color="white", alpha=0.9))
    ax_11.text(xlabel, ylabel, '(a)', transform = ax_11.transAxes)
    ax_2.text(xlabel, ylabel, '(b)', transform = ax_2.transAxes)

    return ax_11, ax_2
