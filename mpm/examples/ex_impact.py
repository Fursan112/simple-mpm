import numpy as np
import time
from mpm_imports import *
import sys


#===============================================================================
def bcZero( x ):
    return( np.zeros(2) )


#===============================================================================
def init():
    #  Initialize Simulation
    # Result File Name
    fName = 'impact'
    fDir = 'test_data/impact'
    
    # Domain Constants
    x0 = np.array([0.0,0.0])                 #  Bottom left Corner
    x1 = np.array([2.,1.])                   #  Top right corner
    nN = np.array([64,32])                   #  Number of cells in each direction
    dx = (x1-x0)/nN
    nG = 2                                   #  Number of ghost nodes
    thick = 0.1                              #  Domain thickness
    ppe = 2                                  #  Particles per element 

    # Material Properties
    mProps1 = {'modulus':1.0e7, 'poisson':0.3, 'density':1.0e3, 
               'maxStress': sys.float_info.max}
    vw = np.sqrt( mProps1['modulus']/mProps1['density'] )   #  Wave speed
    modelName = 'planeStrainNeoHookean'
    mat1 = Material( matid=1, props=mProps1, model=modelName, exload=0)
    mat2 = Material( matid=2, props=mProps1, model=modelName, exload=0)
    mats = [mat1,mat2]
     
    # Time Constants
    t0 = 0.0
    CFL = 0.4
    dt = min(dx) * CFL / vw
    tf = 0.4
    t_out = tf/100.            
    
    # Create Data Warehouse
    dw = Dw( t=0.0, idx=0, sidx=0, tout = t_out, ddir=fDir )
    
    # Create Patch
    pch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, dw )
    
    # Create shape function
    sh = Shape()

    # Create boundary conditions
    bcVy = Bc( 'Y', 0.0, 'gv', bcZero )       # Zero velocity bc at y=0
    bcAy = Bc( 'Y', 0.0, 'ga', bcZero )       # Zero acceleration bc at y=0
    bcVx = Bc( 'X', 0.0, 'gv', bcZero )       # Zero velocity bc at y=0
    bcAx = Bc( 'X', 0.0, 'ga', bcZero )       # Zero acceleration bc at y=0
    pch.bcs = [bcVy,bcAy,bcVx,bcAx]           # BC list
 
    # Create Rectangles
    dx0 = dx[0]
    pt11 = np.array([0.2+2.0*dx0,0.3])
    pt12 = np.array([1.0+2.0*dx0,0.5])
    pt21 = np.array([1.0+4.0*dx0,0.0])
    pt22 = np.array([1.2+4.0*dx0,0.8])    
    r = np.array([0.0,0.2])
    matid1 = 1
    matid2 = 2
    geomutils.fillRectangle( 
        pt11, pt12, ppe, pch, dw, matid1, mProps1['density'])     # Impacter
    geomutils.fillRectangle( 
        pt21, pt22, ppe, pch, dw, matid2, mProps1['density'])     # Target
    
    v0 = np.array([15.0,0.0])
    v1 = np.array([0.0,0.0])
    gravity = np.array([0.0,-9.8])
    
    
    mpm.updateMats( dw, pch, mats, sh )
    mat1.setVelocity( dw, v0 )
    mat2.setVelocity( dw, v1 )
    mat1.setExternalLoad( dw, gravity )
    mat2.setExternalLoad( dw, gravity )
    
    dw.initArrays()
    
    print 'dt = ' + str(pch.dt)    
    return (dw, pch, mats, sh, fName )


#===============================================================================
def stepTime( dw, patch, mats, shape, saveName ):
    # Advance through time
    tbegin = time.time()
    try:
        while( (patch.t < patch.tf) and patch.allInPatch(dw.px) ):
            dw.resetNodes()     
            mpm.timeAdvance( dw, patch, mats, shape )
            dw.saveDataAndAdvance( patch.dt, saveName )
    except JacobError:
        print 'Negative Jacobian'
    
    tend = time.time()
    print (str(dw.idx) + ' iterations in: ' + readTime(tend-tbegin) 
           + ' t=' + str(patch.t) )

#===============================================================================            
def run():
    dw, pch, mats, sh, fname = init()
    stepTime( dw, pch, mats, sh, fname )
    return dw