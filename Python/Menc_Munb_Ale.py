from astropy import constants as cons
from astropy import units as un
import numpy as np
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import sys
plt.rcParams.update({'font.size': 20})

folder = sys.argv[1]
simulation = sys.argv[2]
name = sys.argv[3]

x, y, z, vx, vy, vz, m, u = np.genfromtxt('/tank0/ballone/coll_set/{}/{}'.format(folder,simulation),
					usecols= (0,1,2,3,4,5,6,8), unpack=True)


# Calculate the position of the center of mass.

xcm = np.sum(x*m) / np.sum(m)
ycm = np.sum(y*m) / np.sum(m)
zcm = np.sum(z*m) / np.sum(m)

plt.scatter(x,y,s=1,color='k')
plt.scatter(xcm,ycm,marker='x',color='r')
plt.scatter(x[np.argmax(m)],y[np.argmax(m)],marker='o',color='b')
plt.xlim(-100.,100.)
plt.ylim(-100.,100.)
plt.show()

# Calculate the velocity of the center of mass.

vcmx = np.sum(vx*m) / np.sum(m)
vcmy = np.sum(vy*m) / np.sum(m)
vcmz = np.sum(vz*m) / np.sum(m)

# Ale's patch

xcm = x[np.argmax(m)]
ycm = y[np.argmax(m)]
zcm = z[np.argmax(m)]
vcmx = vx[np.argmax(m)]
vcmy = vy[np.argmax(m)]
vcmz = vz[np.argmax(m)]


# Calculate the radius and velocity wrt the cm.

r = np.sqrt((x-xcm)**2 + (y-ycm)**2 + (z-zcm)**2)
v = np.sqrt((vx-vcmx)**2 + (vy-vcmy)**2 + (vz-vcmz)**2)

# Index sorting the radius from lower to higher.

index = np.argsort(r)


# Reordering of all quantities based on Radius.

X = x[index]
Y = y[index]
Z = z[index]

VX = vx[index]
VY = vy[index]
VZ = vz[index]

R = r[index]
V = v[index]

U = u[index]
M = m[index]

# Generation of the enclosed mass profile.

Mencl = np.cumsum(M)


# Unit conversion for the gravitational constant.

G = (cons.G*(1/un.kg.to(un.M_sun))*(un.m.to(un.R_sun))).value

# Ale's patch

Gsec=((cons.G)/((un.R_sun.to(un.m)**3))*(un.M_sun.to(un.kg))*((1.8845e-2*86400)**2)).value

# Energy calculation for the unbound mass criteria. WITH THE RIGHT NORMALIZED G!

E = V**2 + U - Gsec*Mencl/R

unb = np.where(E > 0)[0]

Munb = M[unb]
percent = np.sum(Munb)*100/np.sum(M)
print np.sum(Munb)

fig = plt.figure(figsize=(10, 7))

plt.plot(R,Mencl,'.k')
plt.text(R[0] + 0.25, R[0] + 0.5,
        'Unbound Mass = {:0.3e} %'.format(percent),
        bbox={'facecolor': 'white', 'pad': 8}, fontsize=15)
plt.xlabel("$R$ [$R_{\odot}$]")
plt.ylabel("$M_{encl}$ [$M_{\odot}$]")
plt.semilogx()
plt.grid()

plt.tight_layout()
plt.show()
#plt.savefig('/tank0/ballone/coll_set/{}/'.format(folder)+name+'.png')
