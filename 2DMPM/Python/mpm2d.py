def timeAdvance( dw, patch, mats ):
    # Advance timestep
    applyExternalLoads( dw, patch, mats )
    interpolateParticlesToGrid( dw, patch, mats )
    computeInternalForce( dw, patch, mats )
    computeAndIntegrateAcceleration( dw, patch, mats )
    setGridBoundaryConditions( dw, patch, mats )
    computeStressTensor( dw, patch, mats )
    interpolateToParticlesAndUpdate( dw, patch, mats )
    
    dw.saveDataAndAdvance(patch.dt)

def applyExternalLoads( dw, patch, mats ):
    # Apply external loads to each material
    for mat in mats:
        mat.applyExternalLoads( dw, patch )
    
def interpolateParticlesToGrid( dw, patch, mats, sh ):
    # Interpolate particle mass and momentum to the grid
    for mat in mats:
            mat.interpolateParticlesToGrid( dw, patch )    
    
def computeInternalForce( dw, patch, mats ):
    # Compute internal body forces
    for mat in mats:
        mat.computeInternalForce( dw, patch )
    
def computeAndIntegrateAcceleration( dw, patch, tol ):
    # Integrate grid acceleration
    dt = patch.dt
    a_leap = 1                                        # Initializes leap-frog
    if( patch.it == 0 ): a_leap = 0.5
    gm = dw.getData( 'gm' )                           # Mass
    gw = dw.getData( 'gw' )                           # Momentum
    gfi = dw.getData( 'gfi' )                         # Internal Force
    gfe = dw.getData( 'gfe' )                         # External Force
    gv = dw.getData( 'gv' )                           # Velocity
    ga = dw.getData( 'ga' )
    
    for ii in range(len(gm)):
        gm[ii] += tol                                 # Make sure no divide by 0
        gv[ii] = gw[ii]/gm[ii]                        # v = momentum/mass
        ga[ii] = a_leap * (gfe[ii]+gfi[ii])/gm[ii]    # a = F/m
        gv[ii] += ga[ii] * dt                         # Integrate velocity
    

def setGridBoundaryConditions( dw, patch, bc ):
    bc.setBoundCond( dw, patch )


def computeStressTensor( dw, patch, mats ):
    for mat in mats:
        mat.computeStressTensor( dw, patch )

    
def interpolateToParticlesAndUpdate( dw, patch, mats ):
    for mat in mats:
        mat.computeStressTensor( dw, patch )