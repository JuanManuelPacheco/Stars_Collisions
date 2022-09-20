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

new_x = (x-xcm)
new_y = (y-ycm)
new_z = (z-zcm)

new_vx = (vx-vcmx)
new_vy = (vy-vcmy)
new_vz = (vz-vcmz)

# Calculate the radius and velocity wrt the cm.

r = np.sqrt(new_x**2 + new_y**2 + new_z**2)
v = np.sqrt(new_vx**2 + new_vy**2 + new_vz**2)

# Index sorting the radius from lower to higher.

index = np.argsort(r)

# Reordering of all quantities based on Radius.

X = new_x[index]
Y = new_y[index]
Z = new_z[index]

VX = new_vx[index]
VY = new_vy[index]
VZ = new_vz[index]

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

# Determination of the central radius where the potential energy aprox fails.

inner_R = np.where(R > 1.5)[0][0]

# Bound unbound criteria for the mass.

unb = np.where(E[inner_R:] > 0)[0] + inner_R

bun_fin = np.where(E[inner_R:] <= 0)[0] + inner_R
bun_ini = np.where(R < R[inner_R])[0]
bun = np.concatenate((bun_ini,bun_fin))

# Calculation of the angular momentum for the bounded particles.

Pos_bn = np.array([X[bun],Y[bun],Z[bun]])
Vel_bn = np.array([VX[bun],VY[bun],VZ[bun]])
P_bn = M[bun]*Vel_bn

L_bn = np.cross(Pos_bn,P_bn, axis = 0)

Lt_bn = np.array([np.sum(L_bn[0]),np.sum(L_bn[1]),np.sum(L_bn[2])])

#Generation of the rotation matrices.

if Lt_bn[2] > 0:
    thx = np.arctan(Lt_bn[1]/Lt_bn[2])
    Rx = np.array([[1,0,0],[0,np.cos(thx),-np.sin(thx)],[0,np.sin(thx),np.cos(thx)]])
    thy = -1*np.arctan(Rx.dot(Lt_bn)[0]/Rx.dot(Lt_bn)[2])
    Ry = np.array([[np.cos(thy),0,np.sin(thy)],[0,1,0],[-np.sin(thy),0,np.cos(thy)]])
    
else:
    thx = np.arctan(Lt_bn[1]/Lt_bn[2])
    Rx = np.array([[1,0,0],[0,np.cos(thx),-np.sin(thx)],[0,np.sin(thx),np.cos(thx)]])
    thy = -1*np.arctan(Rx.dot(Lt_bn)[0]/Rx.dot(Lt_bn)[2]) + np.pi
    Ry = np.array([[np.cos(thy),0,np.sin(thy)],[0,1,0],[-np.sin(thy),0,np.cos(thy)]])

# Rotation to align z vector along the bounded angular momentum.

Pos_bn_rot = Ry.dot(Rx.dot(Pos_bn))
Vel_bn_rot = Ry.dot(Rx.dot(Vel_bn))
P_bn_rot = M[bun]*Vel_bn_rot

# Conversion to a cylindrical coordinate system.

R_bn_rot = np.sqrt(Pos_bn_rot[0]**2+Pos_bn_rot[1]**2)
Vrad = (Vel_bn_rot[0]*Pos_bn_rot[0]+Vel_bn_rot[1]*Pos_bn_rot[1])/R_bn_rot
Vtan = (Vel_bn_rot[1]*Pos_bn_rot[0]-Vel_bn_rot[0]*Pos_bn_rot[1])/R_bn_rot
Vz = Vel_bn_rot[2]

#Conversion to km/s units.

Vrad = Vrad*(un.R_sun.to(un.km))/(1.8845e-2*86400)
Vtan = Vtan*(un.R_sun.to(un.km))/(1.8845e-2*86400)
Vz = Vz*(un.R_sun.to(un.km))/(1.8845e-2*86400)

with open('/tank0/ballone/coll_set/{}/L_evo.txt'.format(folder), 'a') as myfile:
    myfile.write(str(np.mean(Vrad[1:]))+'\t'+ str(np.mean(Vtan[1:]))+'\t'+
                 str(np.mean(Vz[1:]))+'\t'+ str(np.std(Vrad[1:]))+'\t'+ 
                 str(np.std(Vtan[1:]))+'\t'+ str(np.std(Vz[1:]))+'\t'+name+'\n')

#Generation of the histogram for Vr and Vtan.

at = np.histogram(R_bn_rot, bins=100,  weights= Vtan*M[bun])
bt = np.histogram(R_bn_rot, bins=100,  weights= M[bun])
Velt = (at[0]/bt[0])
mean_bint = (at[1][1:] + at[1][:-1]) / 2

ar = np.histogram(R_bn_rot, bins=100,  weights= Vrad*M[bun])
br = np.histogram(R_bn_rot, bins=100,  weights= M[bun])
Velr = (ar[0]/br[0])
mean_binr = (ar[1][1:] + ar[1][:-1]) / 2

#Producing the final figure.

fig = plt.figure(figsize=(10, 7))

plt.plot(mean_binr,Velr,'-b', label = r'$V_{r}$')
plt.plot(mean_bint,Velt,'-r', label = r'$V_{\theta}$')

plt.xlabel("$R$ $[R_{\odot}]$")
plt.ylabel("$V$ $[km/s]$")
plt.grid()
plt.legend()

plt.tight_layout()
plt.savefig('/tank0/ballone/coll_set/{}/Vel_'.format(folder)+name+'.png')


