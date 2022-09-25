from astropy import units
import numpy as np
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 50000
import matplotlib.pyplot as plt
import sys
plt.rcParams.update({'font.size': 20})

import Menc_Munb as Men

def angular_momemtum(pos,vel,m):
    """
    Calculation of the angular momentum of each particle
    and the total angular momentum.

    Parameters
    ----------
    pos : numpy.ndarray
        Vector position of the particles.
    vel : numpy.ndarray
        Vector velocity of the particles.
    m : numpy.ndarray
        Mass of the particles.
  
    Returns
    -------
    L : numpy.ndarray
        Vector angular momentum of the particles.
    Lt : numpy.ndarray
        Total angular momentum.
    """
 
    p = m*vel
    L = np.cross(pos,p, axis = 0)
    Lt = np.array([np.sum(L[0]),np.sum(L[1]),np.sum(L[2])])
    
    return L, Lt

def rotation(Lt,pos,vel):
    """
    Rotation of position an velocity of particles to align
    the z-axis along the total angular momentum.

    Parameters
    ----------
    Lt : numpy.ndarray
        Total angular momentum.
    pos : numpy.ndarray
        Vector position of the particles.
    vel : numpy.ndarray
        Vector velocity of the particles.
  
    Returns
    -------
    pos_rot : numpy.ndarray
        Rotated vector position of the particles.
    vel_rot : numpy.ndarray
        Rotated vector velocity of the particles.
    """
    
    if Lt[2] > 0:
        thx = np.arctan(Lt[1]/Lt[2])
        Rx = np.array([[1,0,0],[0,np.cos(thx),-np.sin(thx)],[0,np.sin(thx),np.cos(thx)]])
        thy = -1*np.arctan(Rx.dot(Lt)[0]/Rx.dot(Lt)[2])
        Ry = np.array([[np.cos(thy),0,np.sin(thy)],[0,1,0],[-np.sin(thy),0,np.cos(thy)]])
    
    else:
        thx = np.arctan(Lt[1]/Lt[2])
        Rx = np.array([[1,0,0],[0,np.cos(thx),-np.sin(thx)],[0,np.sin(thx),np.cos(thx)]])
        thy = -1*np.arctan(Rx.dot(Lt)[0]/Rx.dot(Lt)[2]) + np.pi
        Ry = np.array([[np.cos(thy),0,np.sin(thy)],[0,1,0],[-np.sin(thy),0,np.cos(thy)]])
        
    pos_rot = Ry.dot(Rx.dot(pos))
    vel_rot = Ry.dot(Rx.dot(vel))
    
    return pos_rot, vel_rot

def velocities(pos,vel):
    """
    Conversion of cartesian coordinates to cylindrical
    coordinates.

    Parameters
    ----------
    pos : numpy.ndarray
        Vector position of the particles.
    vel : numpy.ndarray
        Vector velocity of the particles.
  
    Returns
    -------
    R : numpy.ndarray
        Radius of the particles.
    Vrad : numpy.ndarray
        Radial velocity of the particles.
    Vtan : numpy.ndarray
        Tangential velocity of the particles.
    Vz : numpy.ndarray
        Vertical velocity of the particles.
    """
    
    R = np.sqrt(pos[0]**2+pos[1]**2)
    Vrad = (vel[0]*pos[0]+vel[1]*pos[1])/R
    Vtan = (vel[1]*pos[0]-vel[0]*pos[1])/R
    Vz = vel[2]
        
    Vrad = Vrad*(units.R_sun.to(units.km))/(1.8845e-2*86400)
    Vtan = Vtan*(units.R_sun.to(units.km))/(1.8845e-2*86400)
    Vz = Vz*(units.R_sun.to(units.km))/(1.8845e-2*86400)
        
    return R, Vrad, Vtan, Vz

def bining_vel(M,R,Vrad,Vtan,Vz,points):
    """
    Bining of the components of velocity along the 
    radius.

    Parameters
    ----------
    M : numpy.ndarray
        Mass of the particles.
    R : numpy.ndarray
        Radius of the particles.
    Vrad : numpy.ndarray
        Radial velocity of the particles.
    Vtan : numpy.ndarray
        Tangential velocity of the particles.
    Vz : numpy.ndarray
        Vertical velocity of the particles.
    points : int
        Number of bins needed.
  
    Returns
    -------
    mean_bin : numpy.ndarray
        Mean value of the bins.
    Vr : numpy.ndarray
        Binned radial velocity of the particles.
    Vt : numpy.ndarray
        Binned tangential velocity of the particles.
    Vz : numpy.ndarray
        Binned vertical velocity of the particles.
    """    

    ar = np.histogram(R, bins=points, weights= Vrad*M)
    at = np.histogram(R, bins=points, weights= Vtan*M)
    az = np.histogram(R, bins=points, weights= Vz*M)
    
    b = np.histogram(R, bins=points, weights= M)
    
    Vr = (ar[0]/b[0])
    Vt = (at[0]/b[0])
    Vz = (az[0]/b[0])
    
    mean_bin = (b[1][1:] + b[1][:-1]) / 2
    
    return mean_bin, Vr, Vt, Vz

def mean_vel(M,Vrad,Vtan,Vz):
    """
    Computation of mean values and standard deviation of all
    velocity components, arithmetically and weigthed by the mass of 
    particles. 

    Parameters
    ----------
    M : numpy.ndarray
        Mass of the particles.
    R : numpy.ndarray
        Radius of the particles.
    Vrad : numpy.ndarray
        Radial velocity of the particles.
    Vtan : numpy.ndarray
        Tangential velocity of the particles.
    Vz : numpy.ndarray
        Vertical velocity of the particles.
  
    Returns
    -------
    mean : numpy.ndarray
        Mean value of the velocities.
    std : numpy.ndarray
        Standard deviation of the velocities.
    w_mean : numpy.ndarray
        Weigthed mean value of the velocities.
    w_std : numpy.ndarray
        Weigthed standard deviation of the velocities.
    """

    mn_Vr, mn_Vt, mn_Vz = np.nanmean(Vrad), np.nanmean(Vtan), np.nanmean(Vz)
    st_Vr, st_Vt, st_Vz = np.nanstd(Vrad), np.nanstd(Vtan), np.nanstd(Vz)
    
    mkVr = np.where(~np.isnan(Vrad) == True)[0]
    mkVt = np.where(~np.isnan(Vtan) == True)[0]
    mkVz = np.where(~np.isnan(Vz) == True)[0]
    
    vr, vt, vz = Vrad[mkVr], Vtan[mkVt], Vz[mkVz]
    Mr, Mt, Mz = M[mkVr], M[mkVt], M[mkVz]
    
    wm_Vr, wm_Vt, wm_Vz = np.average(vr, weights=Mr), np.average(vt, weights=Mt), np.average(vz, weights=Mz)
    
    ws_Vr = np.sqrt(np.average((vr-wm_Vr)**2, weights=Mr))
    ws_Vt = np.sqrt(np.average((vt-wm_Vt)**2, weights=Mt))
    ws_Vz = np.sqrt(np.average((vz-wm_Vz)**2, weights=Mz))
    
    mean = np.array([mn_Vr,mn_Vt,mn_Vz])
    std = np.array([st_Vr,st_Vt,st_Vz])
    
    w_mean = np.array([wm_Vr,wm_Vt,wm_Vz])
    w_std = np.array([ws_Vr,ws_Vt,ws_Vz])
    
    return mean, std, w_mean, w_std

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

    # Calculation of the angular momentum for the bounded particles.

    pos_bn, vel_bn = np.array([X[bn],Y[bn],Z[bn]]), np.array([VX[bn],VY[bn],VZ[bn]])

    L_bn, Lt_bn = angular_momemtum(pos_bn,vel_bn,M[bn])

    # Rotation to align z vector along the bounded angular momentum.

    pos_bn_rot, vel_bn_rot = rotation(Lt_bn,pos_bn,vel_bn)

    # Conversion to a cylindrical coordinate system.

    R, Vrad, Vtan, Vz = velocities(pos_bn_rot, vel_bn_rot)

    # Computation of the mean and standard deviation of all
    # components of the velocity.

    mean, std, w_mean, w_std = mean_vel(M[bn],Vrad,Vtan,Vz)

    with open('/tank0/ballone/coll_set/{}/L_evo.txt'.format(folder), 'a') as myfile:
              myfile.write(str(mean[0])+'\t'+ str(mean[1])+'\t'+str(mean[2])+'\t'+
              str(std[0])+'\t'+ str(std[1])+'\t'+ str(std[2])+ str(w_mean[0])+'\t'+
              str(w_mean[1])+'\t'+str(w_mean[2])+'\t'+ str(w_std[0])+'\t'+
              str(w_std[1])+'\t'+ str(w_std[2])+'\t'+name+'\n')

    #Generation of the histogram for Vr and Vtan.

    mean_bin, Vr, Vt, Vz = bining_vel(M[bn],R,Vrad,Vtan,Vz,150)

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
