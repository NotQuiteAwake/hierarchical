import numpy as np
import numpy.typing as npt

from Particle import Particle

class Grid:
    # OH WHY DO I HAVE TO BE SCREWED OVER BY SOMETHING LIKE THIS EVERY TIME?
    # I WAS LUCKY TO HAVE USED DATACLASSES IN THE SUMMER :))))))))))
    # mParticles:list[Particle] = []
    
    def __init__(self):
        self.mDim:int = 3
        self.mMaxLim:npt.ArrayLike = np.zeros([self.mDim, 2])
        self.mOctant:npt.ArrayLike = np.zeros([self.mDim, 2])
        self.mParticles:list[Particle] = []

    def GetSize(self):
        return len(self.mParticles)

    def __getitem__(self, ind:int):
        return self.mParticles[ind]
