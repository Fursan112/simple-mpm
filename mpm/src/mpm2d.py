import time
import numpy as np

def timeAdvance( dw, patch, mats, sh ):
    # Advance timestep
    tm = np.zeros(9)
    tm[0] = time.time()    
    sh.updateContribList( dw, patch )            
    tm[1] = applyExternalLoads( dw, patch, mats )
    tm[2] = interpolateParticlesToGrid( dw, patch, mats )
    tm[3] = computeInternalForce( dw, patch, mats )
    tm[4] = computeAndIntegrateAcceleration( dw, patch, patch.tol )
    tm[5] = setGridBoundaryConditions( dw, patch, patch.tol )
    tm[6] = computeStressTensor( dw, patch, mats )
    tm[7] = interpolateToParticlesAndUpdate( dw, patch, mats, sh )
    tm[8] = updateMats( dw, patch, mats, sh )

    for ii in range(8,0,-1):
        tm[ii] = tm[ii] - tm[ii-1]

    return tm
    
def updateMats( dw, patch, mats, sh ):
    t0 = time.time()
    for mat in mats:
        mat.getParticles(dw)
    t2 = time.time()
    return time.time()
    
    
def applyExternalLoads( dw, patch, mats ):
    # Apply external loads to each material
    for mat in mats:
        mat.applyExternalLoads( dw, patch )
    return time.time()

    
def interpolateParticlesToGrid( dw, patch, mats ):
    # Interpolate particle mass and momentum to the grid
    for mat in mats:
        mat.interpolateParticlesToGrid( dw, patch )    
    return time.time()

    
def computeInternalForce( dw, patch, mats ):
    # Compute internal body forces
    for mat in mats:
        mat.computeInternalForce( dw, patch )
    return time.time()

    
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
        
    return time.time()
    

def setGridBoundaryConditions( dw, patch, tol ):
    for bc in patch.bcs:
        bc.setBoundCond( dw, patch, tol )
    return time.time()


def computeStressTensor( dw, patch, mats ):
    for mat in mats:
        mat.computeStressTensor( dw, patch )
    return time.time()

    
def interpolateToParticlesAndUpdate( dw, patch, mats, sh ):
    for mat in mats:
        mat.interpolateToParticlesAndUpdate( dw, patch )
        
    patch.stepTime()
    return time.time()