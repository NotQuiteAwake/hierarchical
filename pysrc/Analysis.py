"""!
@file
@brief Python-side processing of C++ produced data
"""

import statistics
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import scipy as sp
import subprocess

import numpy.typing as npt
import typing as tp
from numbers import Number

import IO
from Grid import Grid

# int_type name mapping
intMapping:dict[str, str] = {'fmm':'FMM',
                              'brute':'Brute-force',
                              'bh':'Barnes-Hut'
                              }

combScript:str = "scripts/combine-vid.sh"

def approx(a:float, b:float) -> bool:
    """!
    @brief Check whether two numbers are approximately equal

    Can be used to test if a number is very nearly an integer.
    """

    eps:float = 1e-5
    return abs(a - b) < eps


def GetResStats(paramList:list, res:dict, int_type:str) -> tuple[list, list]:
    """!
    @brief Get basic statistics for timing results

    @param paramList Order in which to access all param keys in timing results
    @param res Dictionary containing all timing results
    @param int_type Interaction type to consider

    @return lists of average and stdev at each param, ordered as in paramList
    """

    avg_list:list = []
    std_list:list = []
    int_name = intMapping[int_type]

    for p in paramList:
        avg_list.append(np.average(res[p][int_type]))
        std_list.append(np.std(res[p][int_type]))

    return avg_list, std_list


def AnalyseParam(fileName:str,
                 figDir:str,
                 loopPlt:tp.Callable,
                 finalPlt:tp.Callable):
    """!
    @brief The master function for time complexity analysis.

    @param fileName Path to timing results file.
    @param figDir Directory to save figures to.
    @param loopPlt Plotting function to call once for each interaction type.
    @param finalPlt Plotting function finish up and save the plot.
    """

    res = IO.LoadParamResults(fileName)
    IO.MakeDir(figDir)

    param_list:list = sorted(list(res.keys()))
    param_min:int = min(param_list)
    param_max:int = max(param_list)

    # types of interaction up for plotting
    int_types:list = ['brute', 'bh', 'fmm']
    
    if 'brute' not in res[param_max].keys():
        int_types.remove('brute')

        avg_time_brute:float = np.average(res[param_min]['brute'])
        plt.plot(param_list,
                 [avg_time_brute for t in param_list],
                 label = intMapping['brute'])

    for int_type in int_types:
        int_name = intMapping[int_type]
        avg_list, std_list = GetResStats(param_list, res, int_type)
        loopPlt(param_list, avg_list, std_list, int_name)
    
    finalPlt(param_min, param_max)
    plt.clf()


def AnalyseN(fileName:str, figDir:str):
    """!
    @brief Analyse time complexity dependence on number of particles

    @param fileName Path to timing results file
    @param figDir Directory to save figure to
    """

    def loopPlt(n_list, avg_list, std_list, int_name):
        plt.figure(0)
        power_law = lambda x, a, b: a * (x**b)
        a, b = sp.optimize.curve_fit(power_law, n_list, avg_list)[0]
        color = plt.errorbar(n_list,
                             avg_list,
                             yerr = std_list,
                             label = int_name,
                             linestyle = '',
                             marker = 'x'
                             )[0].get_color()
        n_dense = np.linspace(min(n_list), max(n_list), 1000)
        plt.plot(n_dense,
                 [power_law(n, a, b) for n in n_dense],
                 label = f'{int_name} fit, pow = {b:.2f}',
                 color = color)

        plt.figure(1)

        linear = lambda x, a, b: a * x + b
        ln_n = np.log(n_list)
        ln_t = np.log(avg_list)
        a, b = sp.optimize.curve_fit(linear, ln_n, ln_t)[0]

        color = plt.plot(ln_n,
                         ln_t,
                         label = int_name,
                         linestyle = '',
                         marker = 'x')[0].get_color()

        ln_n_dense = np.linspace(ln_n[0], ln_n[-1], 1000)
        plt.plot(ln_n_dense,
                 [linear(lnn, a, b) for lnn in ln_n_dense],
                 label = f'{int_name} ln fit, pow = {a:.2f}',
                 color = color)

    def finalPlt(n_min, n_max):
        plt.figure(0)
        plt.title('Computation time against number of particles $n$')
        plt.legend()
        plt.xlim(0, n_max)
        plt.ylim(0)
        # plt.ylim(0, y_max * 1.1)
        plt.xlabel('Number of particles $n$')
        plt.ylabel('Average computation time $t / \\mu s$')
        plt.savefig(figDir + 'time_n.pdf')
        plt.clf()

        plt.figure(1)
        plt.title('$\\ln{t}$ against $\\ln{N}$')
        plt.legend()
        plt.xlim(np.log(n_min), np.log(n_max))
        # plt.ylim(np.log(y_min), np.log(y_max) * 1.1)
        plt.xlabel('$\\ln{N}$')
        plt.ylabel('$\\ln{t}$')
        # plt.show()
        plt.savefig(figDir + 'time_n_lnln.pdf')

    AnalyseParam(fileName, figDir, loopPlt, finalPlt)


def AnalyseP(fileName:str, figDir:str):
    """!
    @brief Analyse time complexity as function of multipole expansion order

    @param fileName Path to timing results file
    @param figDir Directory to save the plots to
    """

    def loopPlt(p_list, avg_list, std_list, int_name):
        power_law = lambda x, a, b: a * (x**b)
        a, b = sp.optimize.curve_fit(power_law,
                                     p_list,
                                     avg_list)[0]
        p_dense = np.linspace(p_list[0], p_list[-1], 1000)
        color = plt.plot(p_dense,
                         [power_law(p, a, b) for p in p_dense],
                         label = f'{int_name}, pow = {b:.2f}')[0].get_color()
        plt.errorbar(p_list,
                     avg_list,
                     std_list,
                     label = int_name,
                     linestyle = '',
                     marker = 'x',
                     color = color)

    def finalPlt(p_min, p_max):
        plt.legend()
        plt.title('Computation time against order of expansion $p$')
        plt.xlabel('Order of expansion $p$')
        plt.ylabel('Average computation time $t / \\mu s$')
        plt.xlim(p_min, p_max)
        plt.ylim(0)
        plt.savefig(figDir + 'time_p.pdf')

    AnalyseParam(fileName, figDir, loopPlt, finalPlt)


def AnalyseTheta(fileName:str, figDir:str):
    """!
    @brief Analyse time complexity as function of opening angle theta

    @param fileName Path to timing results file
    @param figDir Directory to save the plots to
    """
    def loopPlt(theta_list, avg_list, std_list, int_name):
        plt.errorbar(theta_list,
                     avg_list,
                     std_list,
                     label = int_name,
                     linestyle = '--',
                     marker = 'x')

    def finalPlt(theta_min, theta_max):
        plt.title('Computation time against opening angle $\\theta$')
        plt.xlabel('Opening angle $\\theta$')
        plt.ylabel('Average computation time / $\\mu s$')
        plt.tight_layout()
        plt.legend()
        plt.gca().set_xscale("log")
        plt.xlim(theta_min, theta_max)
        plt.ylim(0)
        plt.savefig(f'{figDir}time_theta.pdf')

    AnalyseParam(fileName, figDir, loopPlt, finalPlt)


def GetErrors(fn_brute:str, fn_inter:str) -> list[float]:
    #gl: grid list
    """!
    @brief Calculate relative errors in calculated acceleration

    @param fn_brute Path to brute force file taken as "ground truth".
    @param fn_inter Path to other interaction file to be compared to brute force

    @return List of errors (all three components, all particles, all repeats)
    """
    gl_brute:list = IO.LoadGrids(fn_brute)
    gl_inter:list = IO.LoadGrids(fn_inter)

    error = []

    for g_brute, g_inter in zip(gl_brute, gl_inter):
        for p_brute, p_inter in zip(g_brute.mParticles,
                                        g_inter.mParticles):
            
            ac_brute = p_brute.accel
            ac_inter = p_inter.accel
            err = np.abs((ac_brute - ac_inter) / (ac_brute))
                  
            error.append(np.mean(err))

    return error


def AnalyseParamError(folderName:str,
                      figDir:str,
                      paramName:str,
                      xlabel:str,
                      logScale:bool = True
                      ):
    """!
    @brief Analyse error distribution and change as a parameter is varied

    Plots both one graph with mean errors of all non-brute interaction types
    compared, and one for each interaction containing change of mean, 5th and
    95th percentile errors against parameter. Colours in the first graph is made
    consistent with timing plots.

    If the folder contains fewer brute calculations than the other interactions,
    it is assumed brute is not affected by parameter varied and everything is
    compared against brute results of the smallest parameter (for same repeats)

    If the recorded param is very close to an integer, it is converted to an
    integer for better plotting results.

    @param folderName Folder containing dump files with one param varying
    @param figDir Directory to save figures (created if doesn't exist.)
    @param paramName Name of parameter varied (for setting text on plots)
    @param xlabel xlabel to use for plots
    @param logScale Whether to use log scale for x-axis of plots.
    """

    IO.MakeDir(figDir)

    fn_by_inter = IO.LoadDumpFolder(folderName)

    param_list:list[float] = list(fn_by_inter['bh'].keys())
    param_min:float = min(param_list)
    param_max:float = max(param_list)

    # process error
    for int_type, param_dict in fn_by_inter.items():
        if int_type == 'brute':
            continue

        int_name:str = intMapping[int_type]
        print(int_type)

        # all errors at each p
        allerr:dict[float, list[float]] = {}

        for param in param_list:
            # displayed param
            param_disp = param
            if approx(int(param_disp), param_disp):
                # if it looks like an interger make it so
                param_disp = int(param_disp)
            else:
                # 3 sig fig.
                # I was just about to praise dynamic typing, when I typed the
                # variable wrong and couldn't see why param_disp isn't updated.
                param_disp = f'{param:.3}'
            print(param_disp)

            # fn: file_name
            fn_brute = fn_by_inter['brute'][param_min]
            # some may have brute calculated only once
            if len(fn_by_inter['brute'].keys()) == len(param_list):
                fn_brute = fn_by_inter['brute'][param]
            fn_inter = fn_by_inter[int_type][param]
            allerr[param] = GetErrors(fn_brute, fn_inter) 

            l = np.percentile(allerr[param], 5) / 10
            r = np.percentile(allerr[param], 95) * 10
            if r == 0:
                print('Skipped (r = 0)')
                continue
            if l == 0:
                l = 1e-16

            divs:int = 30

            plt.title(f"Error distribution for {xlabel} = {param_disp}, " + 
                      f"{int_name}")
            plt.xlabel('Error')
            plt.ylabel('Number of particles')
            plt.hist(allerr[param],
                     bins=np.logspace(np.log10(l), np.log10(r), divs))
            plt.xlim(l, r)
            plt.gca().set_xscale("log")
            plt.savefig(f'{figDir}{int_type}_hist_{param_disp}.pdf')
            plt.clf()

        perc5_list = [np.percentile(allerr[t], 5) for t in param_list]
        perc95_list = [np.percentile(allerr[t], 95) for t in param_list]
        mean_err_list = [np.mean(allerr[t]) for t in param_list]

        def plotter(y:list, err_type:str, **kwargs):
            plt.plot(param_list,
                     y,
                     label = err_type,
                     linestyle = '--',
                     marker = 'x',
                     **kwargs)

        plt.figure(1)
        # matplotlib order for AnaylseX plots
        color_map = {'brute': 'blue', 'bh': "orange", "fmm": 'green'}
        plotter(mean_err_list, int_name, color = color_map[int_type])

        plt.figure(0)
        plotter(mean_err_list, 'mean error')
        plotter(perc5_list, '5th percentile error')
        plotter(perc95_list, '95th percentile error')
        plt.gca().set_yscale('log')
        if logScale:
            plt.gca().set_xscale("log")
        plt.title(f"Error against {xlabel} ({int_name})")
        plt.xlabel(f'{xlabel.capitalize()}')
        plt.ylabel('Relative error')
        plt.legend()
        plt.xlim(param_min, param_max)

        plt.savefig(f'{figDir}err_{paramName}_{int_type}.pdf')
        plt.clf()

    plt.figure(1)
    plt.gca().set_yscale('log')
    if logScale:
        plt.gca().set_xscale("log")
    plt.title(f"Error against {xlabel}")
    plt.xlabel(f'{xlabel.capitalize()}')
    plt.ylabel('Relative error')
    plt.legend()
    plt.xlim(param_min, param_max)

    plt.savefig(f'{figDir}err_{paramName}.pdf')
    plt.clf()


def VisualiseGrid(grid:Grid,
                  scale: float,
                  title:str = None,
                  fig = plt.figure()):
    """!
    @brief Plot particles in a grid on a 3D scatter plot

    Scale determines the limits of the plot. It is assumed that the system is
    centred at the origin.

    This is currently used mainly for animations. To plot a single Grid using
    this function, one might have to do:

    fig = VisualiseGrid(...)
    fig.show()
    input()

    To block the main process.

    @param grid Grid containing all particles of interest
    @param scale Scale to do the plot at (Some particles might lie outside)
    @param title Title of plot
    @param fig Figure object to operate on. Will create one if it's not supplied
    
    @return Figure object with a 3D scatter plot
    """
    mass_list:list = [p.mass for p in grid.mParticles]
    pos_list:list = [p.pos for p in grid.mParticles] 
    
    if len(fig.axes) > 0:
        ax = fig.gca()
        ax.clear()
    else:
        ax = fig.add_subplot(projection = '3d')

    ax.scatter([p[0] for p in pos_list],
               [p[1] for p in pos_list],
               [p[2] for p in pos_list],
               marker = 'o',
               s = np.log(mass_list),
               alpha = 0.8)
    
    ax.set_title(title)
    ax.set_xlim(-scale, scale)
    ax.set_ylim(-scale, scale)
    ax.set_zlim(-scale, scale)

    return fig


def AnimateGrid(grids:dict[float, Grid],
                scale:float,
                figName:str,
                int_type:str = None):
    """!
    @brief Animate the time evolution of a system, saving to disk.

    @param grids Dictionary of system snapshots at each timestamp.
    @param scale Characteristic scale of system
    @param figName Name of the exported animation
    @param int_type Type of interaction for this run
    """
    grid_list:list[Grid] = list(grids.values())
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')

    _title:str = f'{intMapping[int_type]} ' if int_type is not None else ""
    t_list = list(grids.keys())

    def animate(i:int):
        if (i != 0 and i % 10 == 0):
            print()

        print(f'{t_list[i]:.2f}', end = ' ', flush = True)
        VisualiseGrid(grid_list[i],
                      scale,
                      title = _title + f"t = {t_list[i]:.2f}",
                      fig = fig)

    ani = FuncAnimation(fig,
                        animate,
                        frames = len(grids),
                        interval = 30,
                        repeat = False)
    ffwriter = FFMpegWriter(fps=60,
                            bitrate=5000)
    ani.save(figName,
             writer=ffwriter)
    
    # clean after ourselves
    fig.clf()
    print()
    

def GridSnapshots(grids:dict[float, Grid],
                  scale:float,
                  figDir:str,
                  int_type:str = None):
    """!
    @brief Save snapshots of time evolution at interval timestamps.

    @param grids Grids at each timestamp of time evolution
    @param scale Characteristic scale of system
    @param figDir Directory to save figures to
    @param int_type Type of interaction.
    """
    IO.MakeDir(figDir)

    t_list = list(grids.keys())
    t_min:float = min(t_list)
    t_max:float = max(t_list)
    snap_t_list:list[float] = [t for t in t_list if approx(int(t), t)]
    snapshot_list:list[Grid] = [grids[t] for t in snap_t_list]
    
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')

    for t, snapshot in zip(snap_t_list, snapshot_list):
        VisualiseGrid(snapshot,
                      scale,
                      f"{intMapping[int_type]} t = {t:.2f}",
                      fig)
        fig.tight_layout()
        fig.savefig(f'{figDir}{int_type}_{t}.pdf')

    fig.clf()

def AnalyseEvo(folderName:str, figDir:str):
    """!
    @brief For time evolution, generate videos and plots of conserved quantities

    Current plots include energy and components of angular momentum against time 

    @param folderName Folder from which to read time evolution results
    @param figDir Directory to save figures and videos to.

    """
    IO.MakeDir(figDir)

    file_names:list[str] = [folderName + f for f in os.listdir(folderName)]
    file_names = sorted([f for f in file_names if os.path.isfile(f)])

    dim:list[str] = ['x', 'y', 'z']

    for file_name in file_names:
        print(f'Load {file_name}')
        int_type, stats_list, grids, scale = IO.LoadEvo(file_name)
        print(int_type)
        int_name = intMapping[int_type]
        print('Animation')
        AnimateGrid(grids, scale, f'{figDir}{int_type}_ani.mp4', int_type)

        t_list = [s.t for s in stats_list]
        E_list = [s.PE + s.KE for s in stats_list]
        L_list = [s.L for s in stats_list]

        print("Plots")
        plt.figure(0)
        plt.plot(t_list, E_list, label = int_name)
        plt.xlim(min(t_list), max(t_list))
    
        for i in range(3):
            plt.figure(i + 1)
            plt.plot(t_list, [l[i] for l in L_list], label = int_name)
            plt.xlim(min(t_list), max(t_list))

        print('Snapshots')
        GridSnapshots(grids, scale, f'{figDir}snap_{int_type}/', int_type)

    print('Combine videos')
    subprocess.call([combScript, f'{figDir}/', '-y'],
                    stdin = subprocess.DEVNULL)

    print('Conservation')
    plt.figure(0)
    plt.title("Total energy against time")
    plt.xlabel("Time $t$")
    plt.ylabel("Energy $E$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{figDir}E_t.pdf')
    plt.clf()

    for i in range(3):
        plt.figure(i + 1)
        plt.title(f"Total angular momentum in {dim[i]} against time")
        plt.xlabel("Time $t$")
        plt.ylabel(f"Angular momentum $L_{dim[i]}$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{figDir}L_{dim[i]}_t.pdf')
        plt.clf()

