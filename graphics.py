from matplotlib.colors import LinearSegmentedColormap
import numpy as np

#import matplotlib.pyplot as plt
from matplotlib import rcParams


def background(ax_00):
    """
    info: make a gray color gradient for the background
    input: ax_00
    output: ax_00
    """

    cmap = LinearSegmentedColormap.from_list('mycmap', ['gray', 'white'])
    
    lx = 101
    ly = 101
    z = np.zeros([lx, ly])
    for i in range(lx):
        for j in range(ly):
            z[i][j] = i + j

    x = list(range(lx))
    y = list(range(ly))
    x = [el/100 for el in x]
    y = [el/100 for el in y]
    x, y = np.meshgrid(x, y)

    ax_00.pcolor(x, y, z, cmap=cmap, transform=ax_00.transAxes)

    return ax_00

def design(factor=2, factor_ratio=1):
    """
    info: 
    input:
    output: 
    """
    
    # set nice figure sizes
    fig_width_pt = 245*factor    # Get this from LaTeX using \showthe\columnwidth
    golden_mean = (np.sqrt(5.) - 1.) / 2.  # Aesthetic ratio
    ratio = golden_mean
    inches_per_pt = 1. / 72.27  # Convert pt to inches
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width*ratio*factor_ratio # height in inches
    fig_size = [fig_width, fig_height]
    rcParams.update({'figure.figsize': fig_size})
    
    params = { #'backend': 'ps',
            'font.family': 'serif',
            #'font.serif': 'Latin Modern Roman',
            #'font.family': 'sans serif',
            'font.serif': 'Helvetica',
            'font.size': 10,
            'axes.labelsize': 'medium',
            'axes.titlesize': 'medium',
            'legend.fontsize': 'medium',
            'xtick.labelsize': 'small',
            'ytick.labelsize': 'small',
            'savefig.dpi': 15,
            'text.usetex': True
    }
    
    # tell matplotlib about your params
    rcParams.update(params)
    
    # sphinx_gallery_thumbnail_number = 7
    
    #plt.rcParams['savefig.facecolor'] = "0.8"


    """
    params = { #'backend': 'ps',
        'font.family': 'serif',
        #'font.serif': 'Latin Modern Roman',
        #'font.family': 'sans serif',
        'font.serif': 'Helvetica',
        'font.size': 10,
        'axes.labelsize': 'medium',
        'axes.titlesize': 'medium',
        'legend.fontsize': 'medium',
        'xtick.labelsize': 'small',
        'ytick.labelsize': 'small',
        'savefig.dpi': 15,
        'text.usetex': True
    }
    ""
    rcParams['axes.labelsize'] = 50
    rcParams['xtick.labelsize'] = 46
    rcParams['ytick.labelsize'] = 46
    rcParams['legend.fontsize'] = 46
    rcParams['font.family'] = 'sans serif'
    rcParams['font.serif'] = ['Helvetica']
    rcParams['text.usetex'] = True
    rcParams['figure.figsize'] = 12, 8
    ""

    # set nice figure sizes
    fig_width_pt = 245    # Get this from LaTeX using \showthe\columnwidth
    golden_mean = (np.sqrt(5.) - 1.) / 2.  # Aesthetic ratio
    ratio = golden_mean
    inches_per_pt = 1. / 72.27  # Convert pt to inches
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width*ratio  # height in inches
    fig_size = [fig_width, fig_height]
    """

    return fig_size, params 

if __name__ == "__main__": 
    print("\n doing main")
    design()
