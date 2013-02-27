import numpy as np
import scipy.io as sio
import collections


#===============================================================================	
class DataWarehouse:
    # Holds all the data particle and node for an individual timestep used for
    # data access
    def __init__(self, t, idx, bUseMat=False, dim=2, nzeros=6):
	self.t = t
	self.idx = idx
	self.dim = dim
	self.useVTK = tryImportVTK()
	self.useMat  = bUseMat
	self.nzeros  = nzeros
	
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
	self.pvI  = []     # Velocity increment
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


    def saveDataAndAdvance( self, dt, fName ):
	if( self.useMat or (not self.useVTK) ):
	    saveDataMat(fName)
	else:
	    saveDataVTK(fName)
	self.resetNodes()
	
	self.t += dt
	self.idx += 1	


    def saveDataMat( self, fName ):
	fName = fName + str(self.idx).zfill(self.nzeros) + '.mat'
	fl1 = ('px','pMat','pF','pGv','pVS','pv')
	fl2 = ('pfe','pw','pvI','pxI','pm','pVol','pJ')
	fl3 = ('gx','gm','gv','gw','gfe','gfi','ga')
	flist = fl1+fl2+fl3
	fDict = {'pX':self.pX}
	for ii in flist:
	    fDict[ii] = getAttr(self,ii)
	
	sio.savemat( fName, fDict )
	
    def saveDataVTK( self, fName ):
	fName = fName + str(self.idx).zfill(self.nzeros) + '.vtu'
	
	px,py,pz = [], [], []
	vs11,vs12,vs21,vs22 = [], [], [], []
	for ii in range(len(self.px)):
	    px.append( self.px[ii][0] )
	    py.append( self.px[ii][1] )
	    pz.append( 0.0 )
	    vs11.append( self.pVS[ii][0,0] )
	    vs12.append( self.pVS[ii][0,1] )
	    vs21.append( self.pVS[ii][1,0] )
	    vs22.append( self.pVS[ii][1,1] )
	    
	px = np.array(px)
	py = np.array(py)
	pz = np.array(pz)
	vs11 = np.array(vs11)
	vs12 = np.array(vs12)
	vs21 = np.array(vs21)
	vs22 = np.array(vs22)
	vsdat = {"vs11":vs11, "vs12":vs12, "vs21":vs21, "vs22":vs22}
	pointsToVTK(fName, px, py, pz, data = vsdat)
	
	
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

def tryImportVTK():
    try:
	from evtk import pointsToVTK
	bUseVTK = True
    except ImportError:
	bUseVTK = False
    return bUseVTK    