import numpy as np
import time
from mpm_imports import *

# Define Zero Boundary Condition
def bcZero( x ):
    return( 0 )


def init():
    # Result File Name
    fName = 'test_data/mpm_test'
    
    # Domain Constants
    x0 = np.array([0.0,0.0])                 #  Bottom left Corner
    x1 = np.array([1.0,1.0])                 #  Top right corner
    nN = np.array([50,50])                   #  Number of nodes in each direction
    dx = (x1-x0)/nN
    nG = 2                                   #  Number of ghost nodes
    thick = 0.1                              #  Domain thickness
    ppe = 2                                  #  Particles per element 

    # Material Properties
    mProps1 = {'modulus':1.0e7, 'poisson':0.3, 'density':1.0e3 }
    vw = np.sqrt( mProps1['modulus']/mProps1['density'] )   #  Wave speed
    modelName = 'planeStrainNeoHookean'
    mat1 = Material( matid=1, props=mProps1, model=modelName, exload=0)
    mat2 = Material( matid=2, props=mProps1, model=modelName, exload=0)
    mats = [mat1,mat2]
    
    
    # Time Constants
    t0 = 0.0
    CFL = 0.2
    dt = min(dx) * CFL / vw
    tf = 2*dt
    
    # Create Data Warehouse
    dw = Dw( t=0.0, idx=0, ddir='test_data' )
    
    # Create Patch
    pch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, dw )
    
    # Create shape function
    sh = Shape()

    # Create boundary conditions
    bcV = Bc( 'Y', 0.0, 'gv', bcZero )            #  Zero velocity bc at y=0
    bcA = Bc( 'Y', 0.0, 'ga', bcZero )            #  Zero acceleration bc at y=0
    pch.bcs = [bcV,bcA]                           #  BC list
    
    # Create Rectangle
    pt1 = np.array([0.25,0])
    pt2 = np.array([0.75,0.5])
    matid1 = 1
    geomutils.fillRectangle( pt1,pt2,ppe,pch,dw,matid1,mProps1['density'] )
    
    return (dw, pch, mats, sh, fName )
    
    
def stepTime( dw, patch, mats, shape, saveName ):
    while( patch.t < patch.tf ):
        t0 = time.time()
        mpm.timeAdvance( dw, patch, mats, shape )
        t1 = time.time()
        dw.saveDataAndAdvance( patch.dt, saveName )
        t2 = time.time()
        print t1-t0, t2-t1
            
def run():
    dw, pch, mats, sh, fname = init()
    stepTime( dw, pch, mats, sh, fname )