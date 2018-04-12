from matplotlib import pyplot as plt
from cycler import cycler


def update_rcparams():
    plt.rcParams['xtick.labelsize'] = 15
    plt.rcParams['ytick.labelsize'] = 15
    plt.rcParams['xtick.color'] = 'black'
    plt.rcParams['ytick.color'] = 'black'
    plt.rcParams['axes.labelsize'] = 15
    plt.rcParams['axes.labelcolor'] = 'black'
    plt.rcParams['axes.titlesize'] = 20
    plt.rcParams['axes.facecolor'] = '#eeeeee'
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['text.color'] = 'black'
    plt.rcParams['grid.color'] = 'white'
    plt.rcParams['figure.edgecolor'] = 'black'
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 2
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.major.width'] = 2
    plt.rcParams['axes.grid'] = False
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['legend.edgecolor'] = 'k'
    plt.rcParams['legend.framealpha'] = .8
    plt.rcParams['legend.fancybox'] = True
