import numpy

class Particle:
    def __init__(self,pos,mass,vol):
        dim = 2;
        
        Z1 = numpy.zeros([dim,1])
        Z2 = numpy.zeros([dim,dim])
    
        self.pCon = []      # Nodes and weights to which the particle contributes
        self.pF = Z2 * 1.   # Deformation Gradient
        self.pGv = Z2 * 1.  # Velocity Gradient
        self.pVS = Z2 * 1.  # Volume*Stress
        self.pX = pos       # Initial Position
        self.px = pos       # Position
        self.pv = Z1 * 1.   # Velocity
        self.pfe = Z1 * 1.  # External Force
        self.pw = Z1 * 1.   # Momentum
        self.pVI = Z1 * 1.  # Velocity increment
        self.pxI = Z1 * 1.  # Position increment
        self.pm = mass      # Mass
        self.pVol = vol     # Initial Volume
        self.pJ = 0.0       # Jacobian
        

class Node:
    def __init__(self,idx,pos,mass):
        dim = 2;
        Z1 = numpy.zeros([dim,1])  
        
        self.gIdx = idx     # Index
        self.gm = mass      # Mass        
        self.gx = pos       # Position
        self.gv = Z1 * 1.   # Velocity
        self.gw = Z1 * 1.   # Momentum        
        self.gfe = Z1 * 1.  # External Force
        self.gfi = Z1 * 1.  # Internal Force
        self.ga = Z1 * 1.   # Grid acceleration        