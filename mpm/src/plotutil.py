import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as copy
from matplotlib.widgets import Slider, Button, RadioButtons

def plot( data, domain=[0,1,0,1] ):
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0.25,bottom=0.25)
    t = sorted(data.keys())
    dwis = list(set( [jj for (ii,jj) in data[min(t)].dw.keys() ] ))
    
    pl = plt.scatter([0],[0],s=10,linewidths=(0,0,0))
    plt.axis(domain)        
    
    axTime = plt.axes([0.25, 0.1, 0.65, 0.03]) 
    sTime = Slider(axTime, 'Time', min(t), max(t), valinit=min(t))
    
    rax2 = plt.axes([0.05, 0.7, 0.15, 0.15])
    radio = RadioButtons(rax2, ('Vel', 'Stress'))
    
    def updateTime(val):
        t0 = sTime.val
        dt = abs(np.array(t)-t0)
        t_in = t[dt.argmin()]
        px = partVar(data,'Pos',t_in,dwis)
        pc = partVar(data,lbl,t_in,dwis)
        pl.set_offsets(px)  
        pl.set_array(pc)
        plt.draw()
            
    updateTime(t_in)
    sTime.on_changed(updateTime)
    
    def updateVar(lbl):
        px = partVar(data,'Pos',t_in,dwis)
        pc = partVar(data,lbl,t_in,dwis)
        pl.set_offsets(px)  
        pl.set_array(pc)
        plt.draw()    

    radio.on_clicked(updateVar)
        
    plt.show()
    
    
#======================================================================
def partVar( data, t, dwis, lbl='Vel' ):
    funcs = {'Pos':pPosition,'Vel':pVelocity,'Stress':pVonMises}
    func = funcs[lbl]
    out = func( data, t, dwis[0] )
    for ii in range(1,len(dwis)):
        out = np.append(out, func(data,t,dwis[ii]), axis=0 )
        
    return out

#======================================================================
def pPosition( data, t, dwi ):
    px = copy( data[t].get('px',dwi) )
    return px

#======================================================================
def pVelocity( data, t, dwi ):
    pv = copy( data[t].get('pxI',dwi) )
    return vnorm(pv)

#======================================================================
def pVonMises( data, t, dwi ):
    pVS,pVol = data[t].getMult( ['pVS','pVol'], dwi )
    pS = [pVS[ii]/pVol[ii] for ii in range(len(pVol))]
    ms = np.array([vonMises(pp) for pp in pS])
    return ms

#======================================================================    
def vnorm( x ):
    xx = np.copysign(np.sqrt((x*x).sum(axis=1)[:,np.newaxis]), x[:,0])
    return np.reshape( xx, xx.size )

#======================================================================    
def vonMises( S ):
    return np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] +
                    3.*S[1,0]*S[0,1] )