import IO
import statistics
import numpy as np
import matplotlib.pyplot as plt

def AnalyseComplexity(fileName):
    res = IO.LoadTimingResults(fileName) 
    
    for int_name, int_results in res.items():
        n_list:list = []
        timing_list:list = []
        timing_stdev:list = []
        for n, timed in int_results.items():
            mean = statistics.mean(timed) 
            stdev = statistics.stdev(timed)

            n_list.append((n))
            timing_list.append((mean))
        
        plt.errorbar(n_list,
                     timing_list,
                     # yerr = timing_stdev,
                     label = int_name
                     )
         
    plt.title('Average time against number of masses')
    plt.legend()
    plt.xlim()
    plt.show()
    
