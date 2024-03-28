import statistics
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import scipy as sp

import numpy.typing as npt
import typing as tp
from numbers import Number

import IO
from Grid import Grid

# int_type name mapping
int_mapping:dict[str, str] = {'fmm':'FMM',
                              'brute':'Brute-force',
                              'bh':'Barnes-Hut'
                              }

def GetResStats(paramList:list, res:dict, int_type:str):
    avg_list:list = []
    std_list:list = []
    int_name = int_mapping[int_type]

    for p in paramList:
        avg_list.append(np.average(res[p][int_type]))
        std_list.append(np.std(res[p][int_type]))

    return avg_list, std_list


def AnalyseParam(fileName:str,
                 figDir:str,
                 loopPlt:tp.Callable,
                 finalPlt:tp.Callable):
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
                 label = 'brute')

    for int_type in int_types:
        int_name = int_mapping[int_type]
        avg_list, std_list = GetResStats(param_list, res, int_type)
        loopPlt(param_list, avg_list, std_list, int_name)
    
    finalPlt(param_min, param_max)
    plt.clf()


def AnalyseN(fileName:str, figDir:str):
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
        plt.title('Average time against number of masses')
        plt.legend()
        plt.xlim(0, n_max)
        plt.ylim(0)
        # plt.ylim(0, y_max * 1.1)
        plt.xlabel('Number of masses $N$')
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
        plt.title('Computation time against order of multipole expansion')
        plt.xlabel('Order of multipole expansion, $p$')
        plt.ylabel('Average computation time $t / \\mu s$')
        plt.xlim(p_min, p_max)
        plt.ylim(0)
        plt.savefig(figDir + 'time_p.pdf')

    AnalyseParam(fileName, figDir, loopPlt, finalPlt)


def AnalyseTheta(fileName:str, figDir:str):
    def loopPlt(theta_list, avg_list, std_list, int_name):
        plt.errorbar(theta_list,
                     avg_list,
                     std_list,
                     label = int_name,
                     linestyle = '--',
                     marker = 'x')

    def finalPlt(theta_min, theta_max):
        plt.title('Computation time against opening angle $\\theta$')
        plt.xlabel('$\\theta$')
        plt.ylabel('Average computation time / $\\mu s$')
        plt.tight_layout()
        plt.legend()
        plt.gca().set_xscale("log")
        plt.xlim(theta_min, theta_max)
        plt.ylim(0)
        plt.savefig(f'{figDir}time_theta.pdf')

    AnalyseParam(fileName, figDir, loopPlt, finalPlt)


def AnalyseParamError(folderName:str,
                      figDir:str,
                      paramName:str,
                      xlabel:str,
                      logScale:bool = True
                      ):
    IO.MakeDir(figDir)

    file_names:list[str] = [folderName + f for f in os.listdir(folderName)]
    file_names = sorted([f for f in file_names if os.path.isfile(f)])

    # result by type of interaction
    fn_by_inter:dict[str, dict[float, str]] = {
            'brute': {},
            'bh': {},
            'fmm': {}
            }

    for file_name in file_names:
        params:list[str] = file_name.split('/')[-1].strip('.dump').split('_')
        int_type:str = params[0]
        n:int = int(params[1])
        param:float = float(params[2])

        fn_by_inter[int_type][param] = file_name

    param_list:list[float] = sorted(list(fn_by_inter['bh'].keys()))
    param_min:int = min(param_list)
    param_max:int = max(param_list)

    # process error
    for int_type, param_dict in fn_by_inter.items():
        if int_type == 'brute':
            continue

        int_name:str = int_mapping[int_type]
        print(int_type)

        # all errors at each p
        allerr:dict[float, list[float]] = {}

        for param in param_list:
            print(param)
            # fn: file_name
            fn_brute = fn_by_inter['brute'][param_min]
            # some may have brute calculated only once
            if len(fn_by_inter['brute'].keys()) == len(param_list):
                fn_brute = fn_by_inter['brute'][param]
            fn_inter = fn_by_inter[int_type][param]

            #gl: grid list
            gl_brute:list = IO.LoadGrids(fn_brute)
            gl_inter:list = IO.LoadGrids(fn_inter)

            allerr[param] = []
            
            for g_brute, g_inter in zip(gl_brute, gl_inter):
                err_rep = []
                for p_brute, p_inter in zip(g_brute.mParticles,
                                                g_inter.mParticles):
                    
                    ac_brute = p_brute.accel
                    ac_inter = p_inter.accel
                    err = np.abs((ac_brute - ac_inter) / (ac_brute))
                          
                    allerr[param].append(np.mean(err))
    
            l = np.percentile(allerr[param], 5) / 10
            r = np.percentile(allerr[param], 95) * 10
            if r == 0:
                print('Skipped (r = 0)')
                continue
            if l == 0:
                l = 1e-16

            divs:int = 30
            plt.title(f'Error distribution ({paramName} = {param}, {int_name})')
            plt.xlabel('error')
            plt.ylabel('Number of particles')
            plt.hist(allerr[param],
                     bins=np.logspace(np.log10(l), np.log10(r), divs))
            plt.xlim(l, r)
            plt.gca().set_xscale("log")
            plt.savefig(f'{figDir}{int_type}_hist_{param}.pdf')
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
        plt.ylabel('relative error')
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
    plt.ylabel('relative error')
    plt.legend()
    plt.xlim(param_min, param_max)

    plt.savefig(f'{figDir}err_{paramName}.pdf')
    plt.clf()

# TODO: BROKEN with plt.show()?
def VisualiseGrid(grid:Grid,
                  scale: float,
                  title:str = None,
                  fig = plt.figure()):
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
    grid_list:list[Grid] = list(grids.values())
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')

    _title:str = f'{int_mapping[int_type]} ' if int_type is not None else ""
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
    

# the sc2 player, yes.
class Stats:
    t:float
    timing:int
    PE:float
    KE:float
    L:npt.NDArray
    P:npt.NDArray

def LoadEvo(fileName:str)->tuple[str, list[Stats], dict[float, Grid], float]:
    with open(fileName) as file:
        int_type:str = file.readline().strip()
        step_cnts, step, scale =[IO.LoadFloat(x)
                                 for x in file.readline().split()]
        step_cnts = int(step_cnts)

        stats_list:list[Stats] = []
        grids:dict[float, Grid] = {}
        for i in range(step_cnts):
            t:float = i * step
            grids[t] = IO.LoadGrid(file)

            if (i != 0 and i % 10 == 0):
                print()
            print(f'{t:.2f}', end = ' ', flush = True)

            stats = Stats()
            params:list = file.readline().split()
            stats.t = t
            stats.timing = int(params[0])
            stats.PE = IO.LoadFloat(params[1])
            stats.KE = IO.LoadFloat(params[2])
            stats.L = IO.LoadVec(file)
            stats.P = IO.LoadVec(file)

            stats_list.append(stats)

        print()
        return int_type, stats_list, grids, scale


def GridSnapshots(grids:dict[float, Grid],
                  scale:float,
                  figDir:str,
                  int_type:str = None):
    IO.MakeDir(figDir)

    eps:float = 1e-5
    def approx(a:float, b:float):
        return abs(a - b) < eps

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
                      f"{int_mapping[int_type]} t = {t:.2f}",
                      fig)
        fig.savefig(f'{figDir}{int_type}_{t}.pdf')

    fig.clf()

def AnalyseEvo(folderName:str, figDir:str):
    IO.MakeDir(figDir)

    file_names:list[str] = [folderName + f for f in os.listdir(folderName)]
    file_names = sorted([f for f in file_names if os.path.isfile(f)])

    dim:list[str] = ['x', 'y', 'z']

    for file_name in file_names:
        print(f'Load {file_name}')
        int_type, stats_list, grids, scale = LoadEvo(file_name)
        print(int_type)
        int_name = int_mapping[int_type]
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

    plt.figure(0)
    plt.title("Total energy against time")
    plt.xlabel("Time $t$")
    plt.ylabel("Energy $E$")
    plt.legend()
    plt.savefig(f'{figDir}E_t.pdf')
    plt.clf()

    for i in range(3):
        plt.figure(i + 1)
        plt.title(f"Total angular momentum in {dim[i]} against time")
        plt.xlabel("Time $t$")
        plt.ylabel(f"Angular momentum $L_{dim[i]}$")
        plt.legend()
        plt.savefig(f'{figDir}L_{dim[i]}_t.pdf')
        plt.clf()

