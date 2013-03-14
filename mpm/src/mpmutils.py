import numpy as np
dim = 2
# Utils for moving data between particles and grid

def integrate( contrib, pp, gg, idx ):
    # Integrate particle values to grid (p->g)
    for ii in idx:
        for cc in contrib[ii]:
            gg[cc.idx] += pp[ii] * cc.w
    return gg

def interpolate( contrib, pp, gg, idx ):
    # Interpolate grid values to particles pp
    for ii in idx:
        pp[ii] = [0]
        for cc in contrib[ii]:
            pp[ii] += gg[cc.idx] * cc.w
    return pp

def gradient( contrib, pp, gg, idx ):
    # Interpolate a gradient
    for ii in idx:
        pp[ii] = [0]
        for cc in contrib[ii]:
            gR = np.reshape( gg[cc.idx], (dim,1) )
            cg = np.reshape( cc.grad, (1,dim) )
            pp[ii] += np.dot( gR, cg )
    return pp

def divergence( contrib, pp, gg, idx ):
    # Send divergence of particle field to the grid
    for ii in idx:
        for cc in contrib[ii]:
            cg = np.reshape( cc.grad, (dim,1) )            
            gg[cc.idx] -= np.reshape( np.dot( pp[ii], cg ), dim )
    return gg    