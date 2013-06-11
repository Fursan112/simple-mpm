import numpy as np

class Contact:
    def __init__( self, dwis ):
        self.dwis = dwis
        self.nodes = []

    def findIntersection( dw, patch ):
        pass
    

class StickyContact(Contact):
    def __init__( self, dwis ):
        Contact.__init__(self, dwis)
        
    def exchMomentumInterpolated( self, dw, patch, contacts ):
        pass
    
    def exchForceInterpolated( self, dw, patch, contacts ):
        pass

    def exchMomentumIntegrated( self, dw, patch, contacts ):
        pass