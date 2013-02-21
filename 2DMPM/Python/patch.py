import numpy as np

#===============================================================================
class Patch:
    # The computational domain - called Patch to match with Vaango/Uintah
    def __init__(self,X0,X1,Nc,nGhost,t0,dt,th,bc):
        dim = 2
        self.X0 = X0                 # Bottom corner of patch domain
        self.X1 = X1                 # Top corner of patch domain
        self.Nc = Nc+2*nGhost        # Vector of node counts
        self.thick = th              # Thickness
        self.nGhost = nGhost         # Number of Ghost nodes
        self.dX = (X1-X0)/(Nc+1.0)   # Cell size
        self.gridBC = bc             # Grid boundary conditions
        self.t = t0                  # Time
        self.dt = dt                 # Time increment
        self.it = 0                  # Timestep
            
    def initGrid(self, dw):
        for jj in range(self.Nc[1]):
            yy = (jj-self.nGhost)*self.dX[1] + self.X0[1]
            for ii in range(self.Nc[0]):
                xx = (ii-self.nGhost)*self.dX[0] + self.X0[0]
                XX = np.array( [[xx], [yy]] )
                dw.addNode( XX )