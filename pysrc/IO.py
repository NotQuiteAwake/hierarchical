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
    grids:list = []
    with open(fileName) as file:
        repeats:int = int(file.readline())
        for i in range(repeats):
            grids.append(LoadGrid(file))

    return grids


def LoadParamResults(fileName:str) -> dict:
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

