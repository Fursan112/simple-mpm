import numpy as np

def fillRectangle( pt1, pt2, ppe, patch, dw, matid, density ):
    nn = pCeil( (pt2-pt1) / (patch.dX/ppe) )
    ps = (pt2-pt1)/nn
    vol = patch.thick * ps[0] * ps[1]
    for jj in range(int(nn[1])):
        for ii in range(int(nn[0])):
            ns = np.array([ii+0.5,jj+0.5])
            pt = pt1 + ps*ns
            if patch.inPatch( pt ):
                dw.addParticle( matid, pt, vol*density, vol )

def pCeil( x ):
    tol = 1.e-14
    return np.ceil(x-tol)