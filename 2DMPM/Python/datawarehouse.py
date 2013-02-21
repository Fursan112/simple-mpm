import numpy as np
import collections

class DataWarehouse:
    # Holds all the data particle and node for an individual timestep used for
    # data access
    def __init__(self, t, dim=2):
	self.t = t
	self.dim = dim
	
	# Particle variable lists
	self.nParts = 0    # Number of particles
	self.pX   = []     # Initial Position
	self.px   = []     # Position	
	self.pCon = []     # Nodes and weights to which the particle contributes
	self.pMat = []     # Material ID
	self.pF   = []     # Deformation Gradient
	self.pGv  = []     # Velocity Gradient
	self.pVS  = []     # Volume*Stress
	self.pv   = []     # Velocity
	self.pfe  = []     # External Force
	self.pw   = []     # Momentum
	self.pVI  = []     # Velocity increment
	self.pxI  = []     # Position increment
	self.pm   = []     # Mass
	self.pVol = []     # Initial Volume
	self.pJ   = []     # Jacobian	
        
	# Node variable lists
	self.nNodes = 0      # Number of nodes
        self.gx = []         # Position	
        self.gm = []         # Mass        
        self.gv = []         # Velocity
        self.gw = []         # Momentum        
        self.gfe = []        # External Force
        self.gfi = []        # Internal Force
        self.ga = []         # Grid acceleration    

	
    def getData( self, name ):
	source = getattr(self, name)
	return source

    
    def getMatIndex( self, matid ):
	if matid in self.pMat:	
	    matIdx = findall( self.pMat, matid )
	
	return matIdx	
 
    
    def addParticle( self, mat, pos, mass, vol ):
	self.nParts += 1
	
	Z1 = np.zeros([self.dim,1])
	Z2 = np.zeros([self.dim,self.dim])
	I1 = np.ones((1,self.dim))[0]
	I2 = np.diag(I1)

	self.pX.append( pos )          # Initial Position
	self.px.append( pos )          # Position
	self.pCon.append( [] )         # Nodes and weights
	self.pMat.append( mat )        # Material ID
	self.pF.append( I2.copy() )    # Deformation Gradient
	self.pGv.append( Z2.copy() )   # Velocity Gradient
	self.pVS.append( Z2.copy() )   # Volume*Stress
	self.pv.append( Z1.copy() )    # Velocity
	self.pfe.append( Z1.copy() )   # External Force
	self.pw.append( Z1.copy() )    # Momentum
	self.pVI.append( Z1.copy() )   # Velocity increment
	self.pxI.append( Z1.copy() )   # Position increment
	self.pm.append( mass )         # Mass
	self.pVol.append( vol )        # Initial Volume
	self.pJ.append( 0.0 )          # Jacobian

	
    def addNode( self, pos ):
	self.nNodes += 1 
	
	Z1 = np.zeros([self.dim,1])	
	self.gx.append( pos )          # Position
	self.gm.append( 0.0 )          # Mass        
	self.gv.append( Z1.copy )      # Velocity
	self.gw.append( Z1.copy )      # Momentum        
	self.gfe.append( Z1.copy )     # External Force
	self.gfi.append( Z1.copy )     # Internal Force
	self.ga.append( Z1.copy )      # Grid acceleration    		


    def resetNodes( self ):
	#  Reset nodal variables (other than position) to zero
	for ii in range(self.nNodes):
	    Z1 = np.zeros([self.dim,1])	
	    self.gm[ii] = 0.0           # Mass        
	    self.gv[ii] = Z1.copy       # Velocity
	    self.gw[ii] = Z1.copy       # Momentum        
	    self.gfe[ii] = Z1.copy      # External Force
	    self.gfi[ii] = Z1.copy      # Internal Force
	    self.ga[ii] = Z1.copy       # Grid acceleration    		

    def saveDataAndAdvance( self, dt ):
	self.t += dt
	self.resetNodes()

    @staticmethod
    def findall(L, value, start=0):
	idx = []
	ii = start - 1
	try:
	    ii = L.index(value, i+1)
	    idx.append(ii)
	except ValueError:
	    pass
	return idx

#===============================================================================
pContrib = collections.namedtuple('contrib','idx w grad')