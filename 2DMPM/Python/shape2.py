import numpy as np
from element import Patch, Particle

class Shape:
	def __init__(self):
		self.dim = 2;
		self.S = np.zeros([self.dim,1])    # Value of Shape function
		self.G = np.zeros([self.dim,1])    # Value of Shape function derivative	

class GIMP(Shape):
	def __init__(self):
		self.nSupport = 9
		Shape.__init__(self)

	def uS(x,h,l):                    # Smooth spline interpolant
		r = np.abs(x)
		if( r < l ):      S = 1. - (r*r+l*l) / (2.*h*l)
		elif( r < h-l ):  S = 1. - r/h
		elif( r < h+l ):  S = (h+l-r) * (h+l-r) / (4.*h*l)
		else:             S = 0.
		return S

	def uG(x,h,l):
		r = np.abs(x)
		sgnx = np.sign(x)
		if( r < l ):      G = -x/(h*l)
		elif( r < h-l ):  G = -sgnx/h
		elif( r < h+l ):  G = (h+l-r) / (-2.*sgnx*h*l)
		else:             G = 0.
		return G

	def getS( x, h, l ):
		S = np.zeros([length(x),1])
		for ii in range(length(x)):
			S[ii] = uS( x[ii], h[ii], l[ii] )
	
	def getG( x, h, l ):
		G = np.zeros([length(x),1])
		for ii in range(length(x)):
			G[ii] = uG( x[ii], h[ii], l[ii] )



