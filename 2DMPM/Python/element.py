import numpy as np
import shape2 as sh
import collections

#===============================================================================
class Patch:
    # Patch that computation takes place on
    def __init__(self,X0,X1,Nc,nGhost,th):
        dim = 2
        self.X0 = X0                 # Bottom corner of patch domain
        self.X1 = X1                 # Top corner of patch domain
        self.Nc = Nc+2*nGhost        # Vector of node counts
        self.thick = th              # Thickness
        self.nGhost = nGhost         # Number of Ghost nodes
        self.matList = []            # List of materials (objects)
        self.dX = (X1-X0)/(Nc+1.0)   # Cell size
        self.nodes = []              # Node List
        
    def addMaterial( self, mat ):
        self.matList.append( mat )
        
    def initGrid(self):
        for jj in range(self.Nc[1]):
            yy = (jj-self.nGhost)*self.dX[1] + self.X0[1]
            for ii in range(self.Nc[0]):
                xx = (ii-self.nGhost)*self.dX[0] + self.X0[0]
                XX = np.array( [[xx], [yy]] )
                idx = ii + nX * jj
                self.nodes.append( Node(idx,XX) )
	

#===============================================================================
pContrib = collections.namedtuple('contrib','idx w dx dy')


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

		def getCell( self, patch, shape ):                              # Gets lower left node of 4-cell block
			x_sc = (self.px - patch.X0)/patch.dX + patch.nGhost
			idx  = np.floor(x_sc)
			rem  = x_sc - 1.*idx
			ii   = idx[0] if (rem[0]>=0.5) else idx[0]-1
			jj   = idx[1] if (rem[1]>=0.5) else idx[1]-1
			return jj * patch.Nc[0] + ii

		def updateSG( patch, part, cell ):                       # Update the S and G values for a particle position
			r = part.px - patch.node[cell].gx
			dxdy = np.array([patch.dX[0]/patch.dX[1],
				patch.dX[1]/patch.dX[0]])                        # x/y asymmetry
			dfac = 2.0**self.dim                                 # Scale factor - depends on dimension
			lp = np.sqrt( part.pV / (dfac*patch.thick*dxdy) )    # Particle size to integrate over
			S = getS( r, patch.dX, lp * np.diag( part.dF ) )
			G = getG( r, patch.dX, lp * np.diag( part.dF ) )
			return(S,G)

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
    test_patch = Patch( x0, x1, Nc, nGhost, 1.0)
    return test_patch 


#===============================================================================
tp = test_patch_init()
