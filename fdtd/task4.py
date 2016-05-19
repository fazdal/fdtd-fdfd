#%matplotlib inline
from matplotlib import pyplot as plt

import numpy as np
import math as m

fignum = 0

nx = 300
imp0 = 337.0

eps = np.zeros(nx)
for i in range(nx-1):
    if (i >= nx/2):
        eps[i] = 1
    else:
        eps[i] = 9
        
        
n1 = m.sqrt(8)
n2 = m.sqrt(2)

srcori = int(nx/4)              #source origin
srcwid = 30*np.sqrt(eps[0])
srcdel = 10*srcwid              #source delay
nt = int(1.5*nx+srcdel)

ez = np.zeros(nx)
hy = np.zeros(nx)
x = np.arange(0,nx-1,1)

c = 1/np.sqrt(eps[0])
a = (c-1)/(c+1)
b = 2/(c + 1)

hwnp10, ewnp10 = 0,0 # W | ^{n+1} _{0}
hwnm11, ewnm11 = 0,0 # W | ^{n-1} _{1}
hwnp11, ewnp11 = 0,0 # W | ^{n+1} _{1}
hwnm10, ewnm10 = 0,0 # W | ^{n-1} _{0}
hwn0  , ewn0   = 0,0 # W | ^{n  } _{0}
hwn1  , ewn1   = 0,0 # W | ^{n  } _{1}

hwnp1im1, ewnp1im1 = 0,0 # W | ^{n+1} _{i-1}
hwnm1i  , ewnm1i   = 0,0 # W | ^{n-1} _{i  }
hwnp1i  , ewnp1i   = 0,0 # W | ^{n+1} _{i  }
hwnm1im1, ewnm1im1 = 0,0 # W | ^{n-1} _{i-1}
hwnim1  , ewnim1   = 0,0 # W | ^{n  } _{i-1}
hwni    , ewni     = 0,0 # W | ^{n  } _{i  }


emaxleft = np.zeros(nt)
emaxright = np.zeros(nt)
hmaxleft = np.zeros(nt)
hmaxright = np.zeros(nt)


for dt in range(0,nt):
    ######################
    #Magnetic field
    ######################
    hy[x] = hy[x] + (ez[x+1] - ez[x])/imp0

    #abc at left
    hwnp11 = hy[1]
    hwnp10 = -hwnm11 + a*(hwnp11 + hwnm10) + b*(hwn0 + hwn1)
    hy[0] = hwnp10
    hwnm11, hwnm10 = hwn1, hwn0
    hwn1, hwn0  = hwnp11, hwnp10
    
    #abc at right
    hwnp1im1 = hy[-2]
    hwnp1i = - hwnm1im1 + a*(hwnp1im1 + hwnm1i) + b*(hwnp1i + hwnim1)
    hy[-1] = hwnp1i
    hwnm1i, hwnm1im1 = hwni, hwnim1
    hwni, hwnim1  = hwnp1i, hwnp1im1
    
    ######################
    #Electric field
    ###################### 
    ez[x+1] = ez[x+1] + (hy[x+1]-hy[x])*imp0/eps[x]
    ez[srcori] += m.exp(-((dt-srcdel)*(dt-srcdel))/(srcwid*srcwid))

    #abc at left
    ewnp11 = ez[1]
    ewnp10 = -ewnm11 + a*(ewnp11 + ewnm10) + b*(ewn0 + ewn1)
    ez[0] = ewnp10
    ewnm11, ewnm10 = ewn1, ewn0
    ewn1, ewn0  = ewnp11, ewnp10
    
    #abc at right
    ewnp1im1 = ez[-2]
    ewnp1i = - ewnm1im1 + a*(ewnp1im1 + ewnm1i) + b*(ewnp1i + ewnim1)
    ez[-1] = ewnp1i
    ewnm1i, ewnm1im1 = ewni, ewnim1
    ewni, ewnim1  = ewnp1i, ewnp1im1

    plt.hold(True)
    if (dt % 100 == 0 or dt == 1140):
        fignum = fignum + 1
        plt.figure(fignum)
        plt.xlabel("Position")
        plt.ylabel("Amplitude")
        plt.title("Field at t = "+ str(dt))
        plt.plot(ez, label="E-field")
        plt.plot(hy*imp0, label="H-field")
        plt.legend()
    
#    emaxleft[dt] = np.amax(np.absolute(ez[0:((nx/2))]))
#    emaxright[dt] = np.amax(np.absolute(ez[((nx/2)):nx]))
#    hmaxleft[dt] = np.amax(imp0*np.absolute(hy[0:((nx/2))]))
#    hmaxright[dt] = np.amax(imp0*np.absolute(hy[((nx/2)):nx]))

    emaxleft[dt] = np.absolute(ez[nx/4])
    emaxright[dt] = np.absolute(ez[3*nx/4])
    hmaxleft[dt] = np.absolute(imp0*hy[nx/4])
    hmaxright[dt] = np.absolute(imp0*hy[3*nx/4])

eaveleft = np.amax(emaxleft)
eaveright = np.amax(emaxright)
haveleft = np.amax(hmaxleft)
haveright = np.amax(hmaxright)


fdtdratio = (haveright*eaveright)/(haveleft*eaveleft)