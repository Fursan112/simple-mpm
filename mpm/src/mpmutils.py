import numpy as np
# Utils for moving data between particles and grid

def integrate( contrib, pp, gg, idx ):
    # Integrate particle values to grid (p->g)
    for jj in range(len(idx)):
        ii = idx[jj]
        cc = contrib[ii]
        for kk in range(len(cc)):
            gg[cc[kk].idx] += pp[ii] * cc[kk].w
    return gg

def interpolate( contrib, pp, gg, idx ):
    # Interpolate grid values to particles pp
    for jj in range(len(idx)):
        ii = idx[jj]
        pp[ii] = [0]
        cc = contrib[ii]
        for kk in range(len(cc)):
            pp[ii] += gg[cc[kk].idx] * cc[kk].w
    return pp

def gradient( contrib, pp, gg, idx ):
    # Interpolate a gradient
    for jj in range(len(idx)):
        ii = idx[jj]
        pp[ii] = [0]
        cc = contrib[ii]
        dim = gg[cc[0].idx].size
        for kk in range(len(cc)):
            gR = np.reshape( gg[cc[kk].idx], (dim,1) )
            cg = np.reshape( cc[kk].grad, (1,dim) )
            pp[ii] += np.dot( gR, cg )
    return pp

def divergence( contrib, pp, gg, idx ):
    # Send divergence of particle field to the grid
    for jj in range(len(idx)):
        ii = idx[jj]
        cc = contrib[ii]
        dim = cc[0].grad.size
        for kk in range(len(cc)):
            cg = np.reshape( cc[kk].grad, (dim,1) )
            gg[cc[kk].idx] -= np.reshape( np.dot( pp[ii], cg ), dim )
    return gg    