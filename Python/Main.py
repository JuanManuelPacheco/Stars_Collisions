from astropy import constants as cons
from astropy import units
import numpy as np
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 50000
import matplotlib.pyplot as plt
import sys
plt.rcParams.update({'font.size': 20})

import Menc_Munb as Men
import L_evo as Len

if __name__=="__main__":

    folder = sys.argv[1]  # Folder of each simulation.
    simulation = sys.argv[2] # Path of the ascii file.
    name = sys.argv[3]  # Snapshot name.

    x, y, z, vx, vy, vz, m, d, u = np.genfromtxt('/tank0/ballone/coll_set/{}/{}'.format(folder,simulation),
                                          usecols= (0,1,2,3,4,5,6,7,8), unpack=True)

    # Calculate the position and velocity wrt the center of mass.

    nw_x, nw_y, nw_z, nw_vx, nw_vy, nw_vz = Men.re_center(x,y,z,vx,vy,vz,d)

    # Reordering of all quantities based on Radius sorting from lower to higher.

    X, Y, Z, VX, VY, VZ, R, V, M, U = Men.re_order(nw_x,nw_y,nw_z,nw_vx,nw_vy,nw_vz,m,u)

    # Generation of the enclosed mass profile and total mass.

    Mt, M_enc = Men.mass_quantities(M)

    # Energy calculation for the unbound mass definition.

    e = Men.energy(V,U,R,M_enc)

    # Aplication of bound-unbound criteria.

    bn, un = Men.bound_unbound(R,e)

    # Calculating the porcentage of unbound mass in this snapshot.

    percent = np.sum(M[un])*100/Mt

    with open('/tank0/ballone/coll_set/{}/Munb_percent.txt'.format(folder), 'a') as myfile:
              myfile.write(str(round(percent,4))+'\t'+name+'\n')

    # Generation of the enclosed unbound and bound mass profile.

    Mencl_un = np.cumsum(M[un])
    R_un = R[un]
    
    Mencl_bn = np.cumsum(M[bn])
    R_bn = R[bn]

    # Due to the size of the array is needed the decomposition in order to save
    # all the mass profiles.

    np.savetxt('/tank0/ballone/coll_set/{}/Total_'.format(folder)+name+'.txt',
               np.array([R,M_enc]), delimiter='\t')

    np.savetxt('/tank0/ballone/coll_set/{}/Unb_'.format(folder)+name+'.txt',
               np.array([R_un,Mencl_un]), delimiter='\t')

    np.savetxt('/tank0/ballone/coll_set/{}/Bun_'.format(folder)+name+'.txt',
               np.array([R_bn,Mencl_bn]), delimiter='\t')

    fig = plt.figure(figsize=(10, 7))

    plt.plot(R,M_enc,'.k', label='Total mass')
    plt.plot(R_un,Mencl_un,'.r', label='Unbound mass')
    plt.plot(R_bn,Mencl_bn,'.g', label='Bound mass')
    plt.xlabel("$R$ [$R_{\odot}$]")
    plt.ylabel("$M_{encl}$ [$M_{\odot}$]")
    plt.semilogx()
    plt.legend()

    plt.tight_layout()
    plt.savefig('/tank0/ballone/coll_set/{}/'.format(folder)+name+'.png')

    # Calculation of the angular momentum for the bounded particles.

    pos_bn, vel_bn = np.array([X[bn],Y[bn],Z[bn]]), np.array([VX[bn],VY[bn],VZ[bn]])

    L_bn, Lt_bn = Len.angular_momemtum(pos_bn,vel_bn,M[bn])

    # Rotation to align z vector along the bounded angular momentum.

    pos_bn_rot, vel_bn_rot = Len.rotation(Lt_bn,pos_bn,vel_bn)

    # Conversion to a cylindrical coordinate system.

    R, Vrad, Vtan, Vz = Len.velocities(pos_bn_rot, vel_bn_rot)

    # Computation of the mean and standard deviation of all
    # components of the velocity.

    mean, std, w_mean, w_std = Len.mean_vel(M[bn],Vrad,Vtan,Vz)

    with open('/tank0/ballone/coll_set/{}/L_evo.txt'.format(folder), 'a') as myfile:
              myfile.write(str(mean[0])+'\t'+ str(mean[1])+'\t'+str(mean[2])+'\t'+
              str(std[0])+'\t'+ str(std[1])+'\t'+ str(std[2])+ str(w_mean[0])+'\t'+
              str(w_mean[1])+'\t'+str(w_mean[2])+'\t'+ str(w_std[0])+'\t'+
              str(w_std[1])+'\t'+ str(w_std[2])+'\t'+name+'\n')

    #Generation of the histogram for Vr and Vtan.

    mean_bin, Vr, Vt, Vz = Len.bining_vel(M[bn],R,Vrad,Vtan,Vz,150)

    np.savetxt('/tank0/ballone/coll_set/{}/Vel_'.format(folder)+name+'.txt',
               np.array([mean_bin,Vr,Vt,Vz]), delimiter='\t')

    #Producing the final figure.

    fig = plt.figure(figsize=(10, 7))
 
    plt.plot(mean_bin,Vr,'-b', label = r'$V_{r}$')
    plt.plot(mean_bin,Vt,'-r', label = r'$V_{\theta}$')
    plt.plot(mean_bin,Vz,'-g', label = r'$V_{z}$')

    plt.xlabel("$R$ $[R_{\odot}]$")
    plt.ylabel("$V$ $[km/s]$")
    plt.legend()

    plt.tight_layout()
    plt.savefig('/tank0/ballone/coll_set/{}/Vel_'.format(folder)+name+'.png')

