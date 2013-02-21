import numpy as np
from datawarehouse import pContrib


#===============================================================================
class Shape:
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

	@staticmethod
	def uS(x,h,l):                    # Smooth spline interpolant
		r = np.abs(x)
		if( r < l ):      S = 1. - (r*r+l*l) / (2.*h*l)
		elif( r < h-l ):  S = 1. - r/h
		elif( r < h+l ):  S = (h+l-r) * (h+l-r) / (4.*h*l)
		else:             S = 0.
		return S

	@staticmethod
	def uG(x,h,l):                    # Derivative of S
		r = np.abs(x)
		sgnx = np.sign(x)
		if( r < l ):      G = -x/(h*l)
		elif( r < h-l ):  G = -sgnx/h
		elif( r < h+l ):  G = (h+l-r) / (-2.*sgnx*h*l)
		else:             G = 0.
		return G
	
	@staticmethod
	def getS( x, h, l ):              
	# Get S for a point vector
		S = np.zeros([length(x),1])
		for ii in range(length(x)):
			S[ii] = uS( x[ii], h[ii], l[ii] )
	
	@staticmethod
	def getG( x, h, l ):              
	# Get G for a point vector
		G = np.zeros([length(x),1])
		for ii in range(length(x)):
			G[ii] = uG( x[ii], h[ii], l[ii] )

	@staticmethod
	def getCell( dw, patch, idx ):    
	# Gets lower left node of 4-cell block
		pos = dw.px[idx]
		x_sc = (pos - patch.X0)/patch.dX + patch.nGhost
		idx  = np.floor(x_sc)
		rem  = x_sc - 1.*idx
		ii   = idx[0] if (rem[0]>=0.5) else idx[0]-1
		jj   = idx[1] if (rem[1]>=0.5) else idx[1]-1
		return jj * patch.Nc[0] + ii

	@staticmethod
	def updateSG( dw, patch, idx, cell ):           
	# Update the S and G values for a particle position
		r = dw.px[idx] - dw.gx[cell]
		dxdy = np.array([patch.dX[0]/patch.dX[1],
			patch.dX[1]/patch.dX[0]])         # x/y asymmetry
		dfac = 4.0                                # Scale factor
		lp = np.sqrt( dw.pV[idx] / (dfac*patch.thick*dxdy) )  # Particle size
		S = getS( r, patch.dX, lp * np.diag( dw.dF[idx] ) )
		G = getG( r, patch.dX, lp * np.diag( dw.dF[idx] ) )
		return(S,G)
	
	@staticmethod
	def updateContrib( dw, patch, idx, cell ):
		S,G = updateSG( dw, patch, idx, cell )
		w = S[0]*S[1]
		grad = G * S[::-1]                        # Grad = Gx*Sy, Gy*Sx
		pCon = pContrib( idx, s, g )
		return pCon

	@staticmethod
	def updateContribList( dw, patch, matid ):
		