To run the script with default parameters, we type

    python graph_numba.py 

You receive an overview over parameters by typing

    python graph_numba.py -h

The parameter N sets the number of CDR3 beta sequences, which is loaded. N_part refers to the length of a single block. Moreover, len_xy defines the number of blocks in a single row. N should always be chosen, so that

    N >= len_xy * N_part.

The script results in produces a txt file, containing the sparse adjacency matrix, that contains the single-number indices of the edges. Moreover, the script ii

# Graph_tool
Necessary for running the file make_net.py is the package "graph_tool". 
On Linux it can be installed via 
    
    conda create --name gt -c conda-forge graph-tool
    conda activate gt
    
For detailed instructions see the homepage https://graph-tool.skewed.de/.
Since "graph_tool" is incompatible with numba, "graph_tool" has to be deactivated to run the 
files "graph_numba.py" or "convert_indices.py". You can do that via

    conda deactivate
    conda activate
    
Make sure that at the beginning of each line in the terminal we see "(base)", not "(gt)". To run 
"make_net.py" you can activate "graph_tool" again using 

    conda activate gt

# Use the command line to execute files. 
The scripts can be run from the Command line. Examples are

python graph_numba.py --N=10 --N_part=2 --len_xy=5 --src=test_data.txt

python convert_indices.py --single_sidelength=2 --len_xy=5 --format=dense

python graph_ref.py --N=10 --src="test_data.txt"

Don't write input numbers as exponentials, e.g. write 5000 instead of 5*10**3.

// more description ...
