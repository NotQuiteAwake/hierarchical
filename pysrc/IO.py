"""!
@file
@brief IO Library for loading calculations done in C++
"""

import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
import os
from Grid import Grid
from Particle import Particle

def MakeDir(dir:str):
    try:
        os.mkdir(dir)
    except FileExistsError:
        pass

def LoadFloat(floatStr:str) -> float:
    """!
    @brief Load a float from string, either in decimal or hexfloat format.

    @param floatStr String containing a float.

    @return the float number.
    """
    num = 0
    try:
        num = float(floatStr)
    except ValueError:
        num = float.fromhex(floatStr)
    return num


def LoadOctant(file):
    octant = np.zeros([3, 2])
    line = file.readline()
    if int(line) == 0:
        return octant

    for i in range(3):
        line:str = file.readline()
        octant[i] = [LoadFloat(x) for x in line.split()]

    return octant
        

def LoadVec(file) -> npt.NDArray:
    return np.array([LoadFloat(x) for x in file.readline().split()])


def LoadParticle(file) -> Particle:
    mass, charge = [LoadFloat(x) for x in file.readline().split()]
    par = Particle(mass, charge) 
    par.pos = LoadVec(file)
    par.vel = LoadVec(file)
    par.accel = LoadVec(file)
    return par


def LoadGrid(file) -> Grid:
    grid = Grid()

    line:str = file.readline()
    while (line == ""):
        line = file.readline()

    assert line.strip() == "BEGIN GRID"

    grid.mMaxLim = LoadOctant(file)
    grid.mOctant = LoadOctant(file)

    size:int = int(file.readline())

    for i in range(size):
        grid.mParticles.append(LoadParticle(file))

    assert file.readline().strip() == "END GRID"

    return grid


def LoadGrids(fileName:str) -> list[Grid]:
    """!
    @brief Load a series of grids.

    @return A list of the series of grids.
    """
    grids:list = []
    with open(fileName) as file:
        repeats:int = int(file.readline())
        for i in range(repeats):
            grids.append(LoadGrid(file))

    return grids


def LoadParamResults(fileName:str) -> dict:
    """!
    @brief Load timing results for a certain parameter produced by C++

    The returned dictinoary is nested. Its outermost key is the parameter value,
    and the second is interaction type ('brute', 'bh' and 'fmm').

    The value to the second key is a list of timing repeats, integers in
    microseconds.

    @return A nested dictionary of the result. 
    """
    with open(fileName) as file:
        num_params = int(file.readline())
        res:dict = {}

        for i in range(num_params):
            param, int_types = [LoadFloat(x) for x in file.readline().split()]
            int_types = int(int_types)
            res[param] = {}

            for j in range(int_types):
                int_type, repeat = file.readline().split()
                repeat = int(repeat)

                res[param][int_type] = [int(x) for x in file.readline().split()]
    
    return res


def LoadDumpFolder(folderName:str) -> dict:
    """!
    @brief Load and organise fileNames in the dump folder.

    @param folderName Path to folder containing dump files.

    @return A nested dictionary. First key interaction type, second key param
    value (needs context for which parameter it is), value is file name.
    """

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

    for int_type in list(fn_by_inter.keys()):
        param_dict = fn_by_inter[int_type]
        fn_by_inter[int_type] = dict(sorted(param_dict.items()))
        
    return fn_by_inter


# the sc2 player, yes.
class Stats:
    """!
    @brief Class storing some statistics for a timestep of time evolution

    No calculation is carried out Python-end, C++ does all the heavy lifting.
    """

    t:float
    timing:int
    PE:float
    KE:float
    L:npt.NDArray
    P:npt.NDArray


def LoadEvo(fileName:str)->tuple[str, list[Stats], dict[float, Grid], float]:
    """!
    @brief Load a file containing time evolution of a system.

    @return interaction type, a list of Stats, grid at each time, and scale.
    """

    with open(fileName) as file:
        int_type:str = file.readline().strip()
        step_cnts, step, scale =[LoadFloat(x)
                                 for x in file.readline().split()]
        step_cnts = int(step_cnts)

        stats_list:list[Stats] = []
        grids:dict[float, Grid] = {}
        for i in range(step_cnts):
            t:float = i * step
            grids[t] = LoadGrid(file)

            if (i != 0 and i % 10 == 0):
                print()
            print(f'{t:.2f}', end = ' ', flush = True)

            stats = Stats()
            params:list = file.readline().split()
            stats.t = t
            stats.timing = int(params[0])
            stats.PE = LoadFloat(params[1])
            stats.KE = LoadFloat(params[2])
            stats.L = LoadVec(file)
            stats.P = LoadVec(file)

            stats_list.append(stats)

        print()
        return int_type, stats_list, grids, scale


