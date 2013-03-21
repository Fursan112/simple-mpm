import numpy as np
import mpmutils as util
import materialmodel2d as mmodel


#===============================================================================
class Material:
    # Material - holds update functions - default is deformable
    # overridden by RigidMaterial for rigid materials
    def __init__(self, matid, props, model, exload):
        self.matid = matid
        self.props = props
        self.model = model
        self.exload = exload
        self.hasParts = False
        self.pIdx = []
        self.pCon = []        

        
    def getParticles( self, dw ):
        self.pIdx = dw.getMatIndex( self.matid )
        self.hasParts = ( len(self.pIdx) > 0 )
        if self.hasParts:
            self.pCon = dw.getData( 'pCon' )
            
    def setVelocity( self, dw, v ):
        pw = dw.getData( 'pw' )
        pm = dw.getData( 'pm' )
        for ii in self.pIdx:
            pw[ii] = v * pm[ii]
        
            
    def setExternalLoad( self, dw, fe ):
        pfe = dw.getData( 'pfe' )
        pm = dw.getData( 'pm' )        
        for ii in self.pIdx: 
            pfe[ii] = fe * pm[ii]


    def applyExternalLoads( self, dw, patch ):
        # Apply external loads to each material
        pp = dw.getData( 'pfe' )                         # External force
        gg = dw.getData( 'gfe')
        util.integrate( self.pCon, pp, gg, self.pIdx )

            
    def interpolateParticlesToGrid( self, dw, patch ):
        # Interpolate particle mass and momentum to the grid
        pp = dw.getData( 'pm' )                          # Mass
        gg = dw.getData( 'gm')
        util.integrate( self.pCon, pp, gg, self.pIdx )
        
        pp = dw.getData( 'pw' )                          # Momentum
        gg = dw.getData( 'gw')
        util.integrate( self.pCon, pp, gg, self.pIdx )        
     
    def computeStressTensor( self, dw, patch ):
        mm = mmodel.MaterialModel( self.model )
        pf  = dw.getData( 'pF' )                        # Deformation Gradient
        pvs = dw.getData( 'pVS' )                       # Volume * Stress
        pv  = dw.getData( 'pVol' )                      # Volume
        for ii in self.pIdx:
            S,Ja = mm.getStress( self.props, pf[ii] )   # Get stress and det(pf)
            pvs[ii] = S * (pv[ii] * Ja)                 # Stress * deformed volume     
            
    def computeInternalForce( self, dw, patch ):
        # Compute internal body forces - integrate divergence of stress to grid
        pp = dw.getData( 'pVS' )                          # Stress*Volume
        gg = dw.getData( 'gfi')
        util.divergence( self.pCon, pp, gg, self.pIdx )   
            
    def interpolateToParticlesAndUpdate( self, dw, patch ):
        pvI = dw.getData( 'pvI' )
        pxI = dw.getData( 'pxI' )
        pGv = dw.getData( 'pGv' )
        ga  = dw.getData( 'ga' )
        gv  = dw.getData( 'gv' )
        
        util.interpolate( self.pCon, pvI, ga, self.pIdx )
        util.interpolate( self.pCon, pxI, gv, self.pIdx )
        util.gradient( self.pCon, pGv, gv, self.pIdx )
        
        px = dw.getData( 'px' )
        pw = dw.getData( 'pw' )
        pm = dw.getData( 'pm'  )        
        pF = dw.getData( 'pF' )
        for ii in self.pIdx:
            pw[ii] += pvI[ii] * pm[ii] * patch.dt
            px[ii] += pxI[ii] * patch.dt
            pF[ii] += np.dot( pGv[ii], pF[ii] ) * patch.dt