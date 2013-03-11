import numpy as np
from datawarehouse import pContrib


#===============================================================================
class Shape:
    #  Shape functions - compute nodal contributions to particle values
    def __init__(self):
	self.dim = 2;
	self.S = np.zeros([self.dim,1])    # Value of Shape function
	self.G = np.zeros([self.dim,1])    # Value of Shape function derivative	


#===============================================================================
class GIMP(Shape):
    def __init__(self):
	self.nSupport = 9
	self.nGhost = 2
	Shape.__init__(self)

    def uS( self, x, h, l ):                    # Smooth spline interpolant
	r = np.abs(x)
	if( r < l ):      S = 1. - (r*r+l*l) / (2.*h*l)
	elif( r < h-l ):  S = 1. - r/h
	elif( r < h+l ):  S = (h+l-r) * (h+l-r) / (4.*h*l)
	else:             S = 0.
	return S

    def uG( self, x, h, l ):           # Derivative of S
	r = np.abs(x)
	sgnx = np.sign(x)
	if( r < l ):      G = -x/(h*l)
	elif( r < h-l ):  G = -sgnx/h
	elif( r < h+l ):  G = (h+l-r) / (-2.*sgnx*h*l)
	else:             G = 0.
	return G
	

    def getCell( self, dw, patch, idx ):    
	# Gets lower left node of 4-cell block
	pos = dw.px[idx]
	x_sc = (pos - patch.X0)/patch.dX + patch.nGhost
	idx  = np.floor(x_sc)
	print idx
	print pos
	print patch.dX
	rem  = x_sc - 1.*idx
	ii   = idx[0] if (rem[0]>=0.5) else idx[0]-1
	jj   = idx[1] if (rem[1]>=0.5) else idx[1]-1
	return int(jj * patch.Nc[0] + ii)
    

    def updateContribList( self, dw, patch ):
	# Update node contribution list
	nx = patch.Nc[0]
	h = patch.dX
	dxdy = h/h[::-1]
	idxs = np.array([0,1,2,nx,nx+1,nx+2,2*nx,2*nx+1,2*nx+2])
	
	for ii in range(len(dw.pCon)):
	    dw.pCon[ii] = []
	    cc = self.getCell( dw, patch, ii )
	    
	    px = dw.px[ii]
	    lp = np.sqrt( dw.pVol[ii] / (4.0*patch.thick*dxdy) ) 
	    l = lp * np.diag( dw.pF[ii] )		

	    for jj in range(len(idxs)):	
		idx = cc + idxs[jj]
		
		r = px - dw.gx[idx]	
		S = np.zeros(r.size)
		G = np.zeros(r.size)
		for kk in range(len(r)):
		    S[kk] = self.uS( r[kk], h[kk], l[kk] )
		    G[kk] = self.uG( r[kk], h[kk], l[kk] )	
		w = S[0]*S[1]
		grad = G * S[::-1]                        # Grad = Gx*Sy, Gy*Sx		
		
		dw.pCon[ii].append( pContrib( idx, w, grad ) ) 