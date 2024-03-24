import statistics
import os
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
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
            plt.xlabel('% error')
            plt.ylabel('Number of particles')
            plt.hist(allerr[n], bins=np.logspace(np.log10(l), np.log10(r), divs))
            plt.gca().set_xscale("log")
            plt.savefig(figDir + int_type + '_hist_' + str(n) + '.pdf')
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

