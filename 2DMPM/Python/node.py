import numpy

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