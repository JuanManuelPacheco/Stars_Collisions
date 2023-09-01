from astropy import constants as cons
from astropy import units as un
import numpy as np
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 50000
import matplotlib.pyplot as plt
import sys
plt.rcParams.update({'font.size': 20})

folder = sys.argv[1]
simulation = sys.argv[2]
name = sys.argv[3]

x, y, z, vx, vy, vz, m, d, u = np.genfromtxt('/home/pacheco/Collisions/{}/{}'.format(folder,simulation),
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

X = (x-xcm)[index]
Y = (y-ycm)[index]
Z = (z-zcm)[index]

U = u[index]
M = m[index]

# Generation of the enclosed mass profile.

Mencl = np.cumsum(M)

# Unit conversion for the gravitational constant.

G = ((cons.G)/((un.R_sun.to(un.m)**3))*(un.M_sun.to(un.kg))*((1.8845e-2*86400)**2)).value

# Optimized computation for the potential energy of each particle.

PE = np.zeros(len(R))

for i in range(len(R)):
    d = np.linalg.norm(np.array([X,Y,Z]).T - np.array([X,Y,Z])[:,i],axis=1)
    PE[i] = G*np.sum(np.delete(M,i)/np.delete(d,i))

# Energy calculation for the unbound mass criteria.

E = V**2 + U - PE

# Bound unbound criteria for the mass.

unb = np.where(E > 0)[0]
bun = np.where(E <= 0)[0] 

percent = np.sum(M[unb])*100/np.sum(M)

Mencl_un = np.cumsum(M[unb])
R_un = R[unb]
Mencl_bn = np.cumsum(M[bun])
R_bn = R[bun]

with open('/home/pacheco/Collisions/{}/Munb_percent_PE.txt'.format(folder), 'a') as myfile:
    myfile.write(str(round(percent,4))+'\t'+name+'\n')

fig = plt.figure(figsize=(10, 7))

plt.plot(R,G*Mencl/R,'.r', label = r'Aproximation')
plt.plot(R,PE,'.b', label = r'Theoretical')
plt.xlabel("$R$ [$R_{\odot}$]")
plt.ylabel("$Energy$")
plt.semilogx()
plt.vlines(1.5,np.max(PE),np.min(PE), colors='k',
           linestyles='dotted',linewidth = 3, label=r'$R = 1.5$ [$R_{\odot}$]')
plt.legend()
plt.tight_layout()
plt.savefig('/home/pacheco/Collisions/{}/'.format(folder)+name+'_PE.png')

# Energy calculation for the unbound mass criteria using aproximation

E2 = V**2 + U - G*Mencl/R

# Determination of the central radius where the potential energy aprox fails.

inner_R = np.where(R > 1.5)[0][0]

unb2 = np.where(E2[inner_R:] > 0)[0] + inner_R

bun_fin = np.where(E2[inner_R:] <= 0)[0] + inner_R
bun_ini = np.where(R < R[inner_R])[0]
bun2 = np.concatenate((bun_ini,bun_fin))

percent2 = np.sum(M[unb2])*100/np.sum(M)

with open('/home/pacheco/Collisions/{}/Munb_percent_PE.txt'.format(folder), 'a') as myfile:
    myfile.write(str(round(percent2,4))+'\t'+name+'_aprox'+'\n')
