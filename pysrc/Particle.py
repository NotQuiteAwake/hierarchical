"""!
@file
@brief Python side boiled-down representation of C++ Particle class.

Represents a particle
"""

import numpy as np
import numpy.typing as npt

class Particle:
    """!
    @brief Python side boiled-down representation of C++ Particle class.

    Represents a particle
    """

    mDim:int = 3
    pos:npt.ArrayLike = np.zeros(mDim)
    vel:npt.ArrayLike = np.zeros(mDim)
    accel:npt.ArrayLike = np.zeros(mDim)
    mass:float
    charge:float

    def __init__(self, mass:float, charge:float):
        self.mass = mass
        self.charge = charge
