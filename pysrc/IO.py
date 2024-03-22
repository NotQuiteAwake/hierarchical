import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from Grid import Grid
from Particle import Particle

def LoadOctant(file):
    octant = np.zeros([3, 2])
    for i in range(3):
        line:str = file.readline()
        octant[i] = [float(x) for x in line.split()]

    return octant
        
def LoadVec(file) -> npt.ArrayLike:
    return np.array([float(x) for x in file.readline().split()])

def LoadParticle(file) -> Particle:
    mass, charge = [float(x) for x in file.readline().split()]
    par = Particle(mass, charge) 
    par.pos = LoadVec(file)
    par.vel = LoadVec(file)
    par.accel = LoadVec(file)
    return par

def LoadGrid(fileName: str) -> Grid:
    grid = Grid()
    with open(fileName) as file:
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

def LoadTimingResults(fileName:str) -> dict[str, dict[int, list[int]]]:
    with open(fileName) as file:
        line:str = file.readline()
        int_types = int(line)
        int_results:dict = {}

        for i in range(int_types):
            line = file.readline()
            int_name, runs = line.split()
            runs = int(runs) + 1
            
            scaling:dict = {}

            for j in range(runs):
                line = file.readline()
                print(line)
                n, repeats = [int(x) for x in line.split()] 
                
                line = file.readline()
                timed = [int(x) for x in line.split()]
                    
                scaling[n] = timed
            
            int_results[int_name] = scaling

        return int_results

