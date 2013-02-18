import numpy as np


#===============================================================================
class Patch:
    # Patch that computation takes place on
    def __init__(self,X0,X1,Nc,nGhost):
        dim = 2
        self.X0 = X0                 # Bottom corner of patch domain
        self.X1 = X1                 # Top corner of patch domain
        self.Nc = Nc                 # Vector of node counts
        self.nGhost = nGhost         # Number of Ghost nodes
        self.matList = []            # List of materials (objects)
        self.dX = (X1-X0)/(Nc+1.0)   # Cell size
        self.nodes = []              # Node List
        self.initGrid()              # Create Grid of Nodes
        
    def addMaterial( mat ):
        self.matList.append( mat )
        
    def initGrid(self):
        nX = self.Nc[0] + 2*self.nGhost   # Add Ghost node rows
        nY = self.Nc[1] + 2*self.nGhost
        
        for jj in range(nY):
            yy = (jj-self.nGhost)*self.dX[1] + self.X0[1]
            for ii in range(nX):
                xx = (ii-self.nGhost)*self.dX[0] + self.X0[0]
                XX = np.array( [[xx], [yy]] )
                idx = ii + nX * jj
                self.nodes.append( Node(idx,XX) )
        

#===============================================================================
class Particle:
    # Particle object - container for particle variables
    def __init__(self,pos,mass,vol):
        dim = 2
        
        Z1 = np.zeros([dim,1])
        Z2 = np.zeros([dim,dim])
    
        self.pCon = []         # Nodes and weights to which the particle contrib
        self.pF = Z2.copy()    # Deformation Gradient
        self.pGv = Z2.copy()   # Velocity Gradient
        self.pVS = Z2.copy()   # Volume*Stress
        self.pX = pos          # Initial Position
        self.px = pos          # Position
        self.pv = Z1.copy()    # Velocity
        self.pfe = Z1.copy()   # External Force
        self.pw = Z1.copy()    # Momentum
        self.pVI = Z1.copy()   # Velocity increment
        self.pxI = Z1.copy()   # Position increment
        self.pm = mass         # Mass
        self.pVol = vol        # Initial Volume
        self.pJ = 0.0          # Jacobian

        
#===============================================================================
class Node:
    #  Node object - container for node variables
    def __init__(self,idx,pos):
        dim = 2;
        Z1 = np.zeros([dim,1])
        
        self.gIdx = idx       # Index
        self.gm = 0.0         # Mass        
        self.gx = pos         # Position
        self.gv = Z1.copy()   # Velocity
        self.gw = Z1.copy()   # Momentum        
        self.gfe = Z1.copy()  # External Force
        self.gfi = Z1.copy()  # Internal Force
        self.ga = Z1.copy()   # Grid acceleration    

        
#===============================================================================        
def test_patch_init():
    v0 = np.ones([2,1])
    x0 = 0.0 * v0
    x1 = 1.0 * v0
    Nc = 10.0 * v0
    nGhost = 2
    test_patch = Patch( x0, x1, Nc, nGhost )
    return test_patch 


#===============================================================================
tp = test_patch_init()