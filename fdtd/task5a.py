get_ipython().magic('matplotlib inline')
from matplotlib import pyplot as plt
import numpy as np
import math as m

fignum = 0
dpml = 200           # number of PML cells
dom = 400
nx = dom+(2*dpml)
wl = int(nx/10)
imp0 = 377.0                      #impedance

factor = 2
epshost = 1.5
epsmat = 2

nhost = np.sqrt(epshost)
nmat = np.sqrt(epsmat)

scale = nhost/nmat
thickness = wl/factor

#rescale thickness and nmat
nmat = (thickness*nhost)/int(thickness*scale)
thickness = int(thickness*scale)

epsmat = nmat**2
eps = np.zeros(nx)
eps[:] = epshost
eps[int(nx/2):int(nx/2)+thickness] = epsmat
        
c = 1/np.sqrt(epshost)

refeps = np.zeros(nx)
refeps[:] = epshost
refc = 1/np.sqrt(epshost)


srcori = int(nx/10)              #source origin
srcwid = wl*3.0*np.sqrt(max(epshost,epsmat))
srcdel = 10*srcwid              #source delay
nt = int(4*nx+srcdel)

ez = np.zeros(nx)
hy = np.zeros(nx)
x = np.arange(0,nx-1,1)

refez = np.zeros(nx)
refhy = np.zeros(nx)

R0 = 1e-6          # reflection factor
gra = 4             # order of polynomial grading


smax = -((gra+1)*m.log(R0))/(2*imp0*dpml)

es = np.zeros(nx)
hs = np.zeros(nx)

for i in range(dpml):
    #for the left side of the PML
    es[i+1] = smax*((dpml-i-0.5)/dpml)**gra
    hs[i] = smax*((dpml-i)/dpml)**gra  
    
    #for the right side of the PML
    es[nx-i-1] = smax*((dpml-i-0.5)/dpml)**gra 
    hs[nx-i-1] = smax*((dpml-i)/dpml)**gra

ea = np.exp(-es*imp0)-1
eb = np.exp(-es*imp0)

ha = np.exp(-hs*imp0)-1
hb = np.exp(-hs*imp0)

psihy = np.zeros(nx)
psiez = np.zeros(nx)

#reference

refes = np.zeros(nx)
refhs = np.zeros(nx)

for i in range(dpml):
    #for the left side of the PML
    refes[i+1] = smax*((dpml-i-0.5)/dpml)**gra
    refhs[i] = smax*((dpml-i)/dpml)**gra  
    
    #for the right side of the PML
    refes[nx-i-1] = smax*((dpml-i-0.5)/dpml)**gra 
    refhs[nx-i-1] = smax*((dpml-i)/dpml)**gra

refea = np.exp(-refes*imp0)-1
refeb = np.exp(-refes*imp0)

refha = np.exp(-refhs*imp0)-1
refhb = np.exp(-refhs*imp0)

refpsihy = np.zeros(nx)
refpsiez = np.zeros(nx)

amp = 0
for dt in range(0,nt):
    #modified field
    psihy[x] = hb[x]*psihy[x] + ha[x]*(ez[x+1] - ez[x])
    hy[x] = hy[x] + (ez[x+1] - ez[x])/imp0 + psihy[x]/imp0
    
    psiez[x+1] = eb[x+1]*psiez[x+1] + ea[x+1]*(hy[x+1]-hy[x])
    ez[x+1] = ez[x+1] + (hy[x+1]-hy[x])*imp0/eps[x+1] + psiez[x+1]*imp0/eps[x+1]

    if (dt > srcdel):
        amp = 1    
    else:
        amp = m.exp(-((dt-srcdel)*(dt-srcdel))/(srcwid*srcwid))
    ez[srcori] += amp/np.sqrt(epshost)*np.sin(2*np.pi*dt*c/wl)

    #reference field
    refpsihy[x] = refhb[x]*refpsihy[x] + refha[x]*(refez[x+1] - refez[x])
    refhy[x] = refhy[x] + (refez[x+1] - refez[x])/imp0 + refpsihy[x]/imp0
    
    refpsiez[x+1] = refeb[x+1]*refpsiez[x+1] + refea[x+1]*(refhy[x+1]-refhy[x])
    refez[x+1] = refez[x+1] + (refhy[x+1]-refhy[x])*imp0/refeps[x+1] + refpsiez[x+1]*imp0/refeps[x+1]

    if (dt > srcdel):
        amp = 1    
    else:
        amp = m.exp(-((dt-srcdel)*(dt-srcdel))/(srcwid*srcwid))
    refez[srcori] += amp/np.sqrt(epshost)*np.sin(2*np.pi*dt*refc/wl)
    
    plt.hold(True)
    if (dt % 1000 == 0 and dt > 1000 ):
        fignum = fignum + 1
        plt.figure(fignum)
#        plt.xlim(dpml,dom+dpml)
        plt.xlim(dpml,nx/2)
        plt.ylim(-0.0001,0.0001)
        plt.title("Field at t = " + str(dt))
        plt.ylabel("Field Amplitude")
        plt.xlabel("Position")
        plt.plot(ez-refez, label="E-Field")
#        plt.plot(hy*imp0, label="H-Field")
#        plt.legend()