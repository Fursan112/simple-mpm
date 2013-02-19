import numpy as np
from element import Patch, Particle


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

	def uS(x,h,l):                    # Smooth spline interpolant
		r = np.abs(x)
		if( r < l ):      S = 1. - (r*r+l*l) / (2.*h*l)
		elif( r < h-l ):  S = 1. - r/h
		elif( r < h+l ):  S = (h+l-r) * (h+l-r) / (4.*h*l)
		else:             S = 0.
		return S

	def uG(x,h,l):                    # Derivative of S
		r = np.abs(x)
		sgnx = np.sign(x)
		if( r < l ):      G = -x/(h*l)
		elif( r < h-l ):  G = -sgnx/h
		elif( r < h+l ):  G = (h+l-r) / (-2.*sgnx*h*l)
		else:             G = 0.
		return G

	def getS( x, h, l ):              # Get S for a point vector
		S = np.zeros([length(x),1])
		for ii in range(length(x)):
			S[ii] = uS( x[ii], h[ii], l[ii] )
	
	def getG( x, h, l ):              # Get G for a point vector
		G = np.zeros([length(x),1])
		for ii in range(length(x)):
			G[ii] = uG( x[ii], h[ii], l[ii] )

	def getCell( part, patch ):                              # Gets lower left node of 4-cell block
		x_sc = (part.px - patch.X0)/patch.dX + patch.nGhost
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


