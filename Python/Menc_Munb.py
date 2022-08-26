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

# Calculate the position and velocity of the center of mass.

indexm = np.where(m == np.amax(m))[0]
maxmass = len(indexm)

if maxmass == 1:
    xcm = x[np.argmax(m)]
    ycm = y[np.argmax(m)]
    zcm = z[np.argmax(m)]
    
    vcmx = vx[np.argmax(m)]
    vcmy = vy[np.argmax(m)]
    vcmz = vz[np.argmax(m)]
else:
    xcm = np.sum(x[indexm]*m[indexm]) / np.sum(m[indexm])
    ycm = np.sum(y[indexm]*m[indexm]) / np.sum(m[indexm])
    zcm = np.sum(z[indexm]*m[indexm]) / np.sum(m[indexm])
    
    vcmx = np.sum(vx[indexm]*m[indexm]) / np.sum(m[indexm])
    vcmy = np.sum(vy[indexm]*m[indexm]) / np.sum(m[indexm])
    vcmz = np.sum(vz[indexm]*m[indexm]) / np.sum(m[indexm])

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

G = ((cons.G)/((un.R_sun.to(un.m)**3))*(un.M_sun.to(un.kg))*((1.8845e-2*86400)**2)).value

# Energy calculation for the unbound mass criteria.

E = V**2 + U - G*Mencl/R
unb = np.where(E > 0)[0]

Munb = M[unb]
percent = np.sum(Munb)*100/np.sum(M)

Mencl_un = np.cumsum(Munb)
R_un = R[unb]

with open('/tank0/ballone/coll_set/{}/Munb_percent.txt'.format(folder), 'a') as myfile:
    myfile.write(str(round(percent,4))+'\t'+name+'\n')

# Mass enclosed profile.

Mass_array = np.array([Mencl, R,Mencl_un, R_un])

np.savetxt('/tank0/ballone/coll_set/{}/'.format(folder)+name+'.txt',Mass_array
           ,delimiter='\t')

fig = plt.figure(figsize=(10, 7))

plt.plot(R,Mencl,'.k', label='Total mass')
plt.plot(R_un,Mencl_un,'.b', label='Unbound mass')
plt.xlabel("$R$ [$R_{\odot}$]")
plt.ylabel("$M_{encl}$ [$M_{\odot}$]")
plt.semilogx()
plt.legend()
plt.grid()

plt.tight_layout()
plt.savefig('/tank0/ballone/coll_set/{}/'.format(folder)+name+'.png')
