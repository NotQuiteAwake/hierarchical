import numpy as np
import numpy.typing as npt
from Particle import Particle

class Grid:
    mDim:int = 3
    mMaxLim:npt.ArrayLike = np.zeros([mDim, 2])
    mOctant:npt.ArrayLike = np.zeros([mDim, 2])
    mParticles:list[Particle] = []
    
    def GetSize(self):
        return len(self.mParticles)

    def __getitem__(self, ind:int):
        return self.mParticles[ind]
