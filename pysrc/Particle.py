import numpy as np
import numpy.typing as npt

class Particle:
    mDim:int = 3
    pos:npt.ArrayLike = np.zeros(mDim)
    vel:npt.ArrayLike = np.zeros(mDim)
    accel:npt.ArrayLike = np.zeros(mDim)
    mass:float
    charge:float

    def __init__(self, _mass:float, _charge:float):
        mass = _mass
        charge = _charge
