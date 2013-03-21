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

    
    def uSG( self, x, h, l ):
	r = np.abs(x)
	sgnx = np.sign(x)
	if (r<l):      S = 1. - (r*r+l*l) / (2.*h*l);  G = -x/(h*l)
	elif(r<h-l):   S = 1. - r/h;                   G = -sgnx/h
	elif(r<h+l):   S = (h+l-r)*(h+l-r) / (4.*h*l); G = (h+l-r) / (-2.*sgnx*h*l)
	else:          S = G = 0.
	return (S,G)
	

    def getCell( self, dw, patch, idx ):    
	# Gets lower left node of 4-cell block
	pos = dw.px[idx]
	x_sc = (pos - patch.X0)/patch.dX + patch.nGhost
	idx = np.floor(x_sc)
	rem = (x_sc - 1.*idx) >= 0.5
	ii = idx[0] if rem[0] else idx[0]-1
	jj = idx[1] if rem[1] else idx[1]-1
	
	return int(jj * patch.Nc[0] + ii)
    

    def updateContribList( self, dw, patch, mIdx ):
	# Update node contribution list
	nx = patch.Nc[0]
	h = patch.dX
	dxdy = h[::-1]/h
	idxs = np.array([0,1,2,nx,nx+1,nx+2,2*nx,2*nx+1,2*nx+2])
	S = np.zeros(h.size)
	G = np.zeros(h.size)	

	for ii in mIdx:
	    cc = self.getCell( dw, patch, ii )
	    px = dw.px[ii]
	    lp = np.sqrt( dw.pVol[ii] / (4.0*patch.thick*dxdy) ) 
	    l = lp * np.diag( dw.pF[ii] )		

	    for jj in range(9):	
		idx = idxs[jj] + cc 
		r = px - dw.gx[idx]	
		
		for kk in range(len(r)):
		    S[kk],G[kk] = self.uSG( r[kk], h[kk], l[kk] )
		    
		w = S[0]*S[1]
		grad = G * S[::-1]                        # Grad = Gx*Sy, Gy*Sx		
		
		dw.cIdx[ii][jj] = idx
		dw.cW[ii][jj] = w
		dw.cGrad[ii][jj] = grad