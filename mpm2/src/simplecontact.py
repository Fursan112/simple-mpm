import numpy as np

class SimpleContact:
    def __init__( self, dwis ):
        self.dwis = dwis
        self.nodes = []

    def findIntersection( self, dw ):
        self.findIntersectionSimple( dw )
    
    def findIntersectionSimple( self, dw ):
        # Assumes all materials share a common grid
        gm0 = dw.get('gm',self.dwis[0])
        gm1 = dw.get('gm',self.dwis[1])
        self.nodes = np.where( (gm0>0)*(gm1>0) == True )[0]
    

class FreeContact(SimpleContact):
    def __init__( self, dwis ):
        SimpleContact.__init__(self, dwis)
        
    def exchMomentumInterpolated( self, dw ):
        self.findIntersection( dw )
        if self.nodes.any():
            self.exchVals( 'gm', dw )
            self.exchVals( 'gw', dw )
         
    def exchForceInterpolated( self, dw ):
        if self.nodes.any():
            self.exchVals( 'gfi', dw )

    def exchMomentumIntegrated( self ):
        pass
    
    def exchVals( self, lbl, dw ):
        g0 = dw.get(lbl,self.dwis[0])
        g1 = dw.get(lbl,self.dwis[1])
       
        g0[self.nodes] += g1[self.nodes]
        g1[self.nodes] = g0[self.nodes]