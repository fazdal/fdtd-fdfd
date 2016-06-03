get_ipython().magic('matplotlib inline')
from matplotlib import pyplot as plt
import numpy as np
import math as m

fignum = 0

epsilon = 5                      #permittivity
imp0 = 377.0                      #impedance

nx = 400

srcori = int(nx/2)                #source is at center
srcwid = 30.0*np.sqrt(epsilon)    #source width
srcdel = 10*srcwid                #source delay

nt = int((nx+srcdel)*np.sqrt(epsilon))

ez = np.zeros(nx)
hy = np.zeros(nx)
x = np.arange(0,nx-1,1)

R0 = 1e-6          # reflection factor
gra = 4             # order of polynomial grading
dpml = 20           # number of PML cells

smax = -((gra+1)*m.log(R0))/(2*imp0*dpml)

es = np.zeros(nx)
hs = np.zeros(nx)

#polynomial gradng of the conductivity at the boundaries
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

for dt in range(0,nt):
    psihy[x] = hb[x]*psihy[x] + ha[x]*(ez[x+1] - ez[x])
    hy[x] = hy[x] + (ez[x+1] - ez[x])/imp0 + psihy[x]/imp0
    
    psiez[x+1] = eb[x+1]*psiez[x+1] + ea[x+1]*(hy[x+1]-hy[x])
    ez[x+1] = ez[x+1] + (hy[x+1]-hy[x])*imp0/epsilon + psiez[x+1]*imp0/epsilon

    ez[srcori] += m.exp(-((dt-srcdel)*(dt-srcdel))/(srcwid*srcwid))

    plt.hold(True)
    if (dt % 1000 == 0 and dt < 7000 ):
        fignum = fignum + 1
        plt.figure(fignum)
        plt.title("Field at t = " + str(dt))
        plt.ylabel("Field Amplitude")
        plt.xlabel("Position")
        plt.plot(ez, label="E-Field")
        plt.plot(hy*imp0, label="H-Field")
        plt.legend()