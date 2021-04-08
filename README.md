To run the script with default parameters, we type

    python graph_numba.py 

You receive an overview over parameters by typing

    python graph_numba.py -h

The parameter N sets the number of CDR3 beta sequences, which is loaded. N_part refers to the length of a single block. Moreover, len_xy defines the number of blocks in a single row. N should always be chosen, so that

    N >= len_xy * N_part.

The script results in produces a txt file, containing the sparse adjacency matrix, that contains the single-number indices of the edges. Moreover, the script ii
