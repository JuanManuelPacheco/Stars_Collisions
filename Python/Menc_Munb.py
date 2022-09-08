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

x, y, z, vx, vy, vz, m, d, u = np.genfromtxt('/tank0/ballone/coll_set/{}/{}'.format(folder,simulation),
                                          usecols= (0,1,2,3,4,5,6,7,8), unpack=True)

# Calculate the position and velocity of the center of mass.

xcm = x[np.argmax(d)]
ycm = y[np.argmax(d)]
zcm = z[np.argmax(d)]

vcmx = vx[np.argmax(d)]
vcmy = vy[np.argmax(d)]
vcmz = vz[np.argmax(d)]

# Calculate the radius and velocity wrt the cm.

r = np.sqrt((x-xcm)**2 + (y-ycm)**2 + (z-zcm)**2)
v = np.sqrt((vx-vcmx)**2 + (vy-vcmy)**2 + (vz-vcmz)**2)

# Index sorting the radius from lower to higher.

index = np.argsort(r)

# Reordering of all quantities based on Radius.

R = r[index]
V = v[index]

U = u[index]
M = m[index]

# Generation of the enclosed mass profile.

Mencl = np.cumsum(M)

# Unit conversion for the gravitational constant.

G = ((cons.G)/((un.R_sun.to(un.m)**3))*(un.M_sun.to(un.kg))*((1.8845e-2*86400)**2)).value

# Energy calculation for the unbound mass criteria.

E = V**2 - G*Mencl/R
unb = np.where(E > 0)[0]
bun = np.where(E <= 0)[0]

percent = np.sum(M[unb])*100/np.sum(M)

Mencl_un = np.cumsum(M[unb])
R_un = R[unb]
Mencl_bn = np.cumsum(M[bun])
R_bn = R[bun]

# Generation of the profiles.

np.savetxt('/tank0/ballone/coll_set/{}/'.format(folder)+name+'_R.txt', R)
np.savetxt('/tank0/ballone/coll_set/{}/'.format(folder)+name+'_Mencl.txt', Mencl)
np.savetxt('/tank0/ballone/coll_set/{}/'.format(folder)+name+'_R_un.txt', R_un)
np.savetxt('/tank0/ballone/coll_set/{}/'.format(folder)+name+'_Mencl_un.txt', Mencl_un)
np.savetxt('/tank0/ballone/coll_set/{}/'.format(folder)+name+'_R_bn.txt', R_bn)
np.savetxt('/tank0/ballone/coll_set/{}/'.format(folder)+name+'_Mencl_bn.txt', Mencl_bn)

with open('/tank0/ballone/coll_set/{}/Munb_percent.txt'.format(folder), 'a') as myfile:
    myfile.write(str(round(percent,4))+'\t'+name+'\n')

fig = plt.figure(figsize=(10, 7))

plt.plot(R,Mencl,'.k', label='Total mass')
plt.plot(R_un,Mencl_un,'.r', label='Unbound mass')
plt.plot(R_bn,Mencl_bn,'.b', label='Bound mass')
plt.xlabel("$R$ [$R_{\odot}$]")
plt.ylabel("$M_{encl}$ [$M_{\odot}$]")
plt.semilogx()
plt.legend()
plt.grid()

plt.tight_layout()
plt.savefig('/tank0/ballone/coll_set/{}/'.format(folder)+name+'.png')
