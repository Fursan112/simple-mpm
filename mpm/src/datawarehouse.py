import numpy as np
import scipy.io as sio
import collections
import os

#===============================================================================	
class DataWarehouse:
    # Holds all the data particle and node for an individual timestep used for
    # data access
    def __init__(self, 
                 t, idx, sidx, tout, ddir='.', bUseMat=False, dim=2, nzeros=4):
	self.t = t
	self.idx = idx
	self.saveidx = sidx
	self.savet = tout
	self.dim = dim
	self.ddir = ddir
	self.useVTK = tryImportVTK()
	self.useMat  = bUseMat
	self.nzeros  = nzeros

	try:
	    os.mkdir( ddir )	    
	except Exception:
	    tmp = 0
	    
	# Particle variable lists
	self.nParts = 0    # Number of particles
	self.pX   = []     # Initial Position
	self.px   = []     # Position	
	self.pCon = []     # Nodes and weights to which the particle contributes
	self.pMat = []     # Material ID
	self.pF   = []     # Deformation Gradient
	self.pGv  = []     # Velocity Gradient
	self.pVS  = []     # Volume*Stress
	self.pfe  = []     # External Force
	self.pw   = []     # Momentum
	self.pvI  = []     # Velocity increment
	self.pxI  = []     # Position increment
	self.pm   = []     # Mass
	self.pVol = []     # Initial Volume
        
	# Node variable lists
	self.nNodes = 0      # Number of nodes
        self.gx = []         # Position	
        self.gm = []         # Mass        
        self.gv = []         # Velocity
        self.gw = []         # Momentum        
        self.gfe = []        # External Force
        self.gfi = []        # Internal Force
        self.ga = []         # Grid acceleration    

    def initArrays( self ):
	pList = ['pX','px','pF','pGv','pVS','pfe','pw','pvI','pxI','pm','pVol']
	nList = ['gx','gm','gv','gw','gfe','gfi','ga']
	
	for lattr in pList+nList:
	    attr = getattr(self,lattr)
	    attr = np.array(attr)
	
    def getData( self, name ):
	source = getattr(self, name)
	return source

    
    def getMatIndex( self, matid ):
	matIdx = []
	if matid in self.pMat:	
	    matIdx = self.findall( self.pMat, matid )
	
	return matIdx	
 
    
    def addParticle( self, mat, pos, mass, vol ):
	self.nParts += 1
	
	Z1 = np.zeros(self.dim)
	Z2 = np.zeros([self.dim,self.dim])
	I1 = np.ones(self.dim)
	I2 = np.diag(I1)

	self.pX.append( pos.copy() )          # Initial Position
	self.px.append( pos.copy() )          # Position
	self.pCon.append( [] )         # Nodes and weights
	self.pMat.append( mat )        # Material ID
	self.pF.append( I2.copy() )    # Deformation Gradient
	self.pGv.append( Z2.copy() )   # Velocity Gradient
	self.pVS.append( Z2.copy() )   # Volume*Stress
	self.pfe.append( Z1.copy() )   # External Force
	self.pw.append( Z1.copy() )    # Momentum
	self.pvI.append( Z1.copy() )   # Velocity increment
	self.pxI.append( Z1.copy() )   # Position increment
	self.pm.append( mass )         # Mass
	self.pVol.append( vol )        # Initial Volume

	
    def addNode( self, pos ):
	self.nNodes += 1 
	
	Z1 = np.zeros(self.dim)	
	self.gx.append( pos )          # Position
	self.gm.append( 0.0 )          # Mass        
	self.gv.append( Z1.copy() )      # Velocity
	self.gw.append( Z1.copy() )      # Momentum        
	self.gfe.append( Z1.copy() )     # External Force
	self.gfi.append( Z1.copy() )     # Internal Force
	self.ga.append( Z1.copy() )      # Grid acceleration    		


    def resetNodes( self ):
	#  Reset nodal variables (other than position) to zero
	for ii in range(self.nNodes):
	    self.gm[ii] *= 0.0                    
	    self.gv[ii] *= 0.0
	    self.gw[ii] *= 0.0
	    self.gfe[ii] *= 0.0
	    self.gfi[ii] *= 0.0
	    self.ga[ii] *= 0.0    		


    def saveDataAndAdvance( self, dt, fName ):
	rem = self.t % self.savet / self.savet
	tol = dt/2./self.savet
	if (rem < tol) or ((1-rem) < tol):
	    self.saveData(fName)
	    self.saveidx += 1
	
	self.t += dt
	self.idx += 1	


    def saveData( self, fOut ):
	fName = self.ddir + '/' + fOut + str(self.saveidx).zfill(self.nzeros)
	fNameNode = self.ddir + '/' + fOut + "_n" + str(self.saveidx).zfill(self.nzeros)
	self.saveDataVTK(fName)
	self.saveNodeDataVTK(fNameNode)	

	
    def saveDataVTK( self, fName ):
	from evtk.hl import pointsToVTK
	
	px,py,pz = [], [], []
	vs11,vs12,vs21,vs22 = [], [], [], []
	vx, vy, vv = [], [], []
	for ii in range(len(self.px)):
	    px.append( self.px[ii][0] )
	    py.append( self.px[ii][1] )
	    pz.append( 0.0 )
	    vs11.append( self.pVS[ii][0,0] )
	    vs12.append( self.pVS[ii][0,1] )
	    vs21.append( self.pVS[ii][1,0] )
	    vs22.append( self.pVS[ii][1,1] )
	    vx.append( self.pxI[ii][0] )
	    vy.append( self.pxI[ii][1] )
	    vv.append( np.sign(self.pxI[ii][0]) * np.dot(self.pxI[ii],self.pxI[ii] ) )
	    
	px = np.array(px)
	py = np.array(py)
	pz = np.array(pz)
	vs11 = np.array(vs11)
	vs12 = np.array(vs12)
	vs21 = np.array(vs21)
	vs22 = np.array(vs22)
	vx = np.array(vx)
	vy = np.array(vy)
	vv = np.array(vv)
	
	vsdat = {"vs11":vs11, "vs12":vs12, "vs21":vs21, "vs22":vs22, "vx":vx, "vy":vy, "v":vv}
	pointsToVTK(fName, px, py, pz, data = vsdat)
	
	
    def saveNodeDataVTK( self, fName ):
	from evtk.hl import pointsToVTK
		
	gx,gy,gz = [], [], []
	vx,vy,ax,ay = [], [], [], []
	for ii in range(len(self.gx)):
	    gx.append( self.gx[ii][0] )
	    gy.append( self.gx[ii][1] )
	    gz.append( 0.0 )
	    vx.append( self.gv[ii][0] )
	    vy.append( self.gfi[ii][0] )
	    ax.append( self.ga[ii][0] )
	    ay.append( self.gw[ii][0] )		
		    
	gx = np.array(gx)
	gy = np.array(gy)
	gz = np.array(gz)
	vx = np.array(vx)
	vy = np.array(vy)
	ax = np.array(ax)
	ay = np.array(ay)	    
		
	vsdat = {"vx":vx, "vy":vy, "ax":ax, "ay":ay}
	pointsToVTK(fName, gx, gy, gz, data = vsdat)	
	
	
    def findall(self, L, value, start=0):
	idx = []
	ii = start - 1
	while True:
	    try:
		ii = L.index(value, ii+1)
		idx.append(ii)
	    except ValueError:
		break
	return idx

#===============================================================================
pContrib = collections.namedtuple('contrib','idx w grad')

def tryImportVTK():
    try:
	from evtk.hl import pointsToVTK
	bUseVTK = True
    except ImportError:
	bUseVTK = False
    return bUseVTK    