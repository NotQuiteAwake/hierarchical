import statistics
import os
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import scipy as sp

import IO
from Grid import Grid

def AnalyseComplexity(fileName:str, figDir:str):
    res = IO.LoadTimingResults(fileName) 

    x_max:int = 0
    x_min:int = float('inf')
    y_max:float = 0
    y_min:float = float('inf') 
    for int_name, int_results in res.items():
        n_list:list = []
        timing_list:list = []
        timing_stdev:list = []
        for n, timed in int_results.items():
            mean = statistics.mean(timed) 
            stdev = statistics.stdev(timed)

            n_list.append((n))
            timing_list.append((mean))
            timing_stdev.append(stdev)
            
            x_max = max(x_max, n)
            x_min = min(x_min, n)
            y_max = max(y_max, mean)
            y_min = min(y_min, mean)

        plt.figure(0)
        power_law = lambda x, a, b: a * (x**b)
        a, b = sp.optimize.curve_fit(power_law, n_list, timing_list)[0]
        color = plt.errorbar(n_list,
                             timing_list,
                             yerr = timing_stdev,
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
        ln_t = np.log(timing_list)
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
         
    plt.figure(0)
    plt.title('Average time against number of masses')
    plt.legend()
    plt.xlim(0, x_max)
    plt.ylim(0, y_max * 1.1)
    plt.xlabel('Number of masses $N$')
    plt.ylabel('Average calculation time $t / \\mu s$')
    plt.savefig(figDir + 'time_mass.pdf')
    plt.clf()

    plt.figure(1)
    plt.title('$\\ln{t}$ against $\\ln{N}$')
    plt.legend()
    plt.xlim(np.log(x_min), np.log(x_max))
    plt.ylim(np.log(y_min), np.log(y_max) * 1.1)
    plt.xlabel('$\\ln{N}$')
    plt.ylabel('$\\ln{t}$')
    # plt.show()
    plt.savefig(figDir + 'time_mass_lnln.pdf')
    plt.clf()


def AnalyseError(folderName:str, figDir:str):
    file_names:list[str] = [folderName + f for f in os.listdir(folderName)]
    file_names = sorted([f for f in file_names if os.path.isfile(f)])

    # result by type of interaction
    fn_by_inter:dict[str, dict[int, str]] = {
            'brute': {},
            'bh': {},
            'fmm': {}
            }

    # because I forgot to print this in the .dump files. Haha.
    NUM_REPEATS:int = 10

    # file read-in
    for file_name in file_names:
        param:list[str] = file_name.split('/')[-1].split('.')[0].split('_')
        int_type:str = param[0]
        n:int = int(param[1])
            
        fn_by_inter[int_type][n] = file_name

    # process error
    for int_type, n_dict in fn_by_inter.items():
        if int_type == 'brute':
            continue

        print(int_type)

        # mean error by axis
        mean_err_ax:list[dict[int, float]] = [{} for i in range(3)]
        # error at 'perc' percentile by axis
        perc_err_ax:list[dict[int, float]] = [{} for i in range(3)]
        perc = 75
        # all errors at each n
        allerr:dict[int, list[float]] = {}
        n_list:list[int] = []

        for n in n_dict.keys():
            n_list.append(n)

            print(n)

            # fn: file_name
            fn_brute = fn_by_inter['brute'][n]
            fn_inter = fn_by_inter[int_type][n]

            #gl: grid list
            gl_brute:list = IO.LoadGrids(fn_brute, NUM_REPEATS)
            gl_inter:list = IO.LoadGrids(fn_inter, NUM_REPEATS)

            err_n:list[list] = []
            allerr[n] = []
            
            for g_brute, g_inter in zip(gl_brute, gl_inter):
                err_rep = []
                for p_brute, p_inter in zip(g_brute.mParticles,
                                                g_inter.mParticles):
                    
                    ac_brute = p_brute.accel
                    ac_inter = p_inter.accel
                    err = np.abs((ac_brute - ac_inter) / (ac_brute))
                          
                    allerr[n].append(np.mean(err))
                    err_rep.append(err)

                err_n.append(err_rep)

            for i in range(3):
                # error for n and axis i
                err_n_ax:list = [err_n[j][i] for j in range(len(err_n))]
                mean_err_ax[i][n] = np.mean(err_n_ax)
                perc_err_ax[i][n] = np.percentile(err_n_ax, perc)

        for i in range(3):
            mean_err_list = list(mean_err_ax[i].values())
            perc_err_list = list(perc_err_ax[i].values())
            plt.scatter(np.log(n_list),
                        mean_err_list,
                        label = 'mean')
            plt.scatter(np.log(n_list),
                        perc_err_list,
                        label = '75 percentile')
            plt.legend()
            plt.savefig(figDir + int_type + '_err_ax_' + str(i) + '.pdf')
            plt.clf()
    
        for n in sorted(n_list): 
            # these sometimes have absurdly small errors
            if n < 500:
                continue
            l = np.percentile(allerr[n], 5) / 10
            r = np.percentile(allerr[n], 95) * 10
            divs:int = 30
            plt.title(f'Error distribution (N = {n}, {int_type})')
            plt.xlabel('error')
            plt.ylabel('Number of particles')
            plt.hist(allerr[n], bins=np.logspace(np.log10(l), np.log10(r), divs))
            plt.xlim(l, r)
            plt.gca().set_xscale("log")
            plt.savefig(f'{figDir}{int_type}_hist_{n}.pdf')
            plt.clf()
        

def AnalyseExpansionOrder(fileName:str, figDir:str):
    n, res = IO.LoadExpansionOrderResults(fileName)

    p_list:list = list(res.keys())
    min_p:int = min(p_list)
    max_p:int = max(p_list)

    avg_time_brute:float = np.average(res[min_p]['brute'])
    plt.plot(p_list, [avg_time_brute for p in p_list], label = 'brute')

    for int_type in ['bh', 'fmm']:
        avg_list:list[float] = []
        std_list:list[float] = []
        for p in p_list:
            avg_list.append(np.average(res[p][int_type]))
            std_list.append(np.std(res[p][int_type]))


        power_law = lambda x, a, b: a * (x**b)
        a, b = sp.optimize.curve_fit(power_law,
                                     p_list,
                                     avg_list)[0]
        p_dense = np.linspace(p_list[0], p_list[-1], 1000)
        color = plt.plot(p_dense,
                         [power_law(p, a, b) for p in p_dense],
                         label = f'{int_type}, pow = {b:.2f}')[0].get_color()
        plt.errorbar(p_list,
                     avg_list,
                     std_list,
                     label = int_type,
                     linestyle = '',
                     marker = 'x',
                     color = color)

    plt.legend()
    plt.title('Calculation time against order of multipole expansion')
    plt.xlabel('Order of multipole expansion, $p$')
    plt.ylabel('Average run-time $t / \\mu s$')
    plt.xlim(min_p, max_p)
    plt.ylim(0)
    plt.savefig(figDir + 'time_p.pdf')
    plt.clf()


def AnalyseExpansionError(folderName:str, figDir:str):
    file_names:list[str] = [folderName + f for f in os.listdir(folderName)]
    file_names = sorted([f for f in file_names if os.path.isfile(f)])

    # result by type of interaction
    fn_by_inter:dict[str, dict[int, str]] = {
            'brute': {},
            'bh': {},
            'fmm': {}
            }

    # because I forgot to print this in the .dump files. Haha.
    NUM_REPEATS:int = 10

    # file read-in
    for file_name in file_names:
        param:list[str] = file_name.split('/')[-1].split('.')[0].split('_')
        int_type:str = param[0]
        n:int = int(param[1])
        p:int = int(param[2])
            
        fn_by_inter[int_type][p] = file_name

    p_list:list[int] = sorted(list(fn_by_inter['bh'].keys()))
    p_min:int = min(p_list)
    p_max:int = max(p_list)

    # process error
    for int_type, p_dict in fn_by_inter.items():
        if int_type == 'brute':
            continue

        print(int_type)

        # all errors at each p
        allerr:dict[int, list[float]] = {}

        for p in p_list:
            print(p)

            # fn: file_name
            fn_brute = fn_by_inter['brute'][p_min]
            fn_inter = fn_by_inter[int_type][p]

            #gl: grid list
            gl_brute:list = IO.LoadGrids(fn_brute, NUM_REPEATS)
            gl_inter:list = IO.LoadGrids(fn_inter, NUM_REPEATS)

            allerr[p] = []
            
            for g_brute, g_inter in zip(gl_brute, gl_inter):
                err_rep = []
                for p_brute, p_inter in zip(g_brute.mParticles,
                                                g_inter.mParticles):
                    
                    ac_brute = p_brute.accel
                    ac_inter = p_inter.accel
                    err = np.abs((ac_brute - ac_inter) / (ac_brute))
                          
                    allerr[p].append(np.mean(err))
    
        for p in p_list:
            l = np.percentile(allerr[p], 5) / 10
            r = np.percentile(allerr[p], 95) * 10
            divs:int = 30
            plt.title(f'Error distribution (p = {p}, {int_type})')
            plt.xlabel('error')
            plt.ylabel('Number of particles')
            plt.hist(allerr[p],
                     bins=np.logspace(np.log10(l), np.log10(r), divs))
            plt.xlim(l, r)
            plt.gca().set_xscale("log")
            plt.savefig(f'{figDir}{int_type}_hist_{p}.pdf')
            plt.clf()

        perc5_list = [np.percentile(allerr[p], 5) for p in p_list]
        perc95_list = [np.percentile(allerr[p], 95) for p in p_list]
        mean_err_list = [np.mean(allerr[p]) for p in p_list]

        def plotter(y:list, err_type:str):
            lny = np.log(y)

            color = plt.plot(p_list,
                             lny,
                             label = err_type,
                             linestyle = '--',
                             marker = 'x')[0].get_color()

        plotter(mean_err_list, 'mean error')
        plotter(perc5_list, '5th percentile error')
        plotter(perc95_list, '95th percentile error')
        plt.title(f"Error against order of expansion $p$ ({int_type})")
        plt.xlabel('Order $p$')
        plt.ylabel('ln(relative error)')
        plt.legend()
        plt.xlim(min(p_list), max(p_list))

        plt.savefig(f'{figDir}err_p_{int_type}')
        plt.clf()


def AnalyseTheta(fileName:str, figDir:str):
    n, res = IO.LoadThetaResults(fileName)

    theta_list:list = sorted(list(res.keys()))
    theta_min:int = min(theta_list)
    theta_max:int = max(theta_list)

    avg_time_brute:float = np.average(res[theta_min]['brute'])
    plt.plot(theta_list, [avg_time_brute for t in theta_list], label = 'brute')

    for int_type in ['bh', 'fmm']:
        avg_list:list[float] = []
        std_list:list[float] = []
        for t in theta_list:
            avg_list.append(np.average(res[t][int_type]))
            std_list.append(np.std(res[t][int_type]))


       # power_law = lambda x, a, b: a * (x**b)
       # a, b = sp.optimize.curve_fit(power_law,
       #                              theta_list,
       #                              avg_list)[0]
       # p_dense = np.linspace(p_list[0], p_list[-1], 1000)
       # color = plt.plot(p_dense,
       #                  [power_law(p, a, b) for p in p_dense],
       #                  label = f'{int_type}, pow = {b:.2f}')[0].get_color()
        plt.errorbar(theta_list,
                     avg_list,
                     std_list,
                     label = int_type,
                     linestyle = '--',
                     marker = 'x')

    plt.title('Computation time against opening angle')
    plt.xlabel('Opening angle $\\theta$')
    plt.ylabel('Time / $\\mu s$')
    plt.legend()
    plt.xlim(theta_min, theta_max)
    plt.savefig(f'{figDir}time_theta.pdf')
    plt.clf()


def AnalyseThetaError(folderName:str, figDir:str):
    file_names:list[str] = [folderName + f for f in os.listdir(folderName)]
    file_names = sorted([f for f in file_names if os.path.isfile(f)])

    # result by type of interaction
    fn_by_inter:dict[str, dict[float, str]] = {
            'brute': {},
            'bh': {},
            'fmm': {}
            }

    # because I forgot to print this in the .dump files. Haha.
    NUM_REPEATS:int = 10

    # file read-in
    for file_name in file_names:
        param:list[str] = file_name.split('/')[-1].strip('.dump').split('_')
        int_type:str = param[0]
        n:int = int(param[1])
        theta:float = float(param[2])
            
        fn_by_inter[int_type][theta] = file_name

    theta_list:list[float] = sorted(list(fn_by_inter['bh'].keys()))
    theta_min:float = min(theta_list)
    theta_max:float = max(theta_list)

    # process error
    for int_type, p_dict in fn_by_inter.items():
        if int_type == 'brute':
            continue

        print(int_type)

        # all errors at each p
        allerr:dict[float, list[float]] = {}

        for theta in theta_list:
            print(theta)

            # fn: file_name
            fn_brute = fn_by_inter['brute'][theta_min]
            fn_inter = fn_by_inter[int_type][theta]

            #gl: grid list
            gl_brute:list = IO.LoadGrids(fn_brute, NUM_REPEATS)
            gl_inter:list = IO.LoadGrids(fn_inter, NUM_REPEATS)

            allerr[theta] = []
            
            for g_brute, g_inter in zip(gl_brute, gl_inter):
                err_rep = []
                for p_brute, p_inter in zip(g_brute.mParticles,
                                                g_inter.mParticles):
                    
                    ac_brute = p_brute.accel
                    ac_inter = p_inter.accel
                    err = np.abs((ac_brute - ac_inter) / (ac_brute))
                          
                    allerr[theta].append(np.mean(err))
    
        for theta in theta_list:
            l = np.percentile(allerr[theta], 5) / 10
            r = np.percentile(allerr[theta], 95) * 10
            divs:int = 30
            plt.title(f'Error distribution (theta = {theta}, {int_type})')
            plt.xlabel('error')
            plt.ylabel('Number of particles')
            plt.hist(allerr[theta],
                     bins=np.logspace(np.log10(l), np.log10(r), divs))
            plt.xlim(l, r)
            plt.gca().set_xscale("log")
            plt.savefig(f'{figDir}{int_type}_hist_{theta}.pdf')
            plt.clf()

        perc5_list = [np.percentile(allerr[t], 5) for t in theta_list]
        perc95_list = [np.percentile(allerr[t], 95) for t in theta_list]
        mean_err_list = [np.mean(allerr[t]) for t in theta_list]

        def plotter(y:list, err_type:str):
            lny = np.log(y)

            color = plt.plot(theta_list,
                             lny,
                             label = err_type,
                             linestyle = '--',
                             marker = 'x')[0].get_color()

        plotter(mean_err_list, 'mean error')
        plotter(perc5_list, '5th percentile error')
        plotter(perc95_list, '95th percentile error')
        plt.title(f"Error against opening angle $\\theta$")
        plt.xlabel('Opening angle $\\theta$')
        plt.ylabel('ln(relative error)')
        plt.legend()
        plt.xlim(theta_min, theta_max)

        plt.savefig(f'{figDir}err_theta_{int_type}.pdf')
        plt.clf()

def VisualiseGrid(grid:Grid,
                  scale: float,
                  title:str = None,
                  fig = plt.figure()):
    mass_list:list = [p.mass for p in grid.mParticles]
    pos_list:list = [p.pos for p in grid.mParticles] 
    
    x_list:list = [p[0] for p in pos_list]
    y_list:list = [p[1] for p in pos_list]
    z_list:list = [p[2] for p in pos_list]

    ax = fig.gca()
    ax.clear()
    ax.scatter(x_list,
               y_list,
               z_list,
               marker = 'o',
               s = np.log(mass_list),
               alpha = 0.8)
    
    ax.set_title(title)
    # TODO: have c++ end output a scale parameter
    ax.set_xlim(-scale, scale)
    ax.set_ylim(-scale, scale)
    ax.set_zlim(-scale, scale)


def AnimateGrid(grids:dict[float, Grid],
                scale:float,
                figName:str,
                int_type:str = None):
    grid_list:list[Grid] = list(grids.values())
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')

    _title:str = ""
    t_list = list(grids.keys())
    if int_type is not None:
        _title = int_type + " "

    def animate(i:int):
        if (i != 0 and i % 10 == 0):
            print()

        print(f'{t_list[i]:.2f}', "", end = '')
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
        int_name:str = file.readline().strip()
        step_cnts, step, scale =[IO.LoadFloat(x)
                                 for x in file.readline().split()]
        step_cnts = int(step_cnts)

        stats_list:list[Stats] = []
        grids:dict[float, Grid] = {}
        for i in range(step_cnts):
            t:float = i * step
            grids[t] = IO.LoadGrid(file)

            stats = Stats()
            params:list = file.readline().split()
            stats.t = t
            stats.timing = int(params[0])
            stats.PE = IO.LoadFloat(params[1])
            stats.KE = IO.LoadFloat(params[2])
            stats.L = IO.LoadVec(file)
            stats.P = IO.LoadVec(file)

            stats_list.append(stats)

        return int_name, stats_list, grids, scale

def AnalyseEvo(folderName:str, figDir:str):
    file_names:list[str] = [folderName + f for f in os.listdir(folderName)]
    file_names = sorted([f for f in file_names if os.path.isfile(f)])

    dim:list[str] = ['x', 'y', 'z']

    for file_name in file_names:
        int_name, stats_list, grids, scale = LoadEvo(file_name)
        print(int_name)
        print('Animation')
        AnimateGrid(grids, scale, f'{figDir}{int_name}_ani.mp4', int_name)

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
