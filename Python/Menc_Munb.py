from astropy import constants as cons
from astropy import units
import numpy as np
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 50000
import matplotlib.pyplot as plt
import sys
plt.rcParams.update({'font.size': 20})

def re_center(x,y,z,vx,vy,vz,d):
    """
    Redefinition of all coordinates and velocities wrt to
    the maximum density particle of each snapshot.

    Parameters
    ----------
    x, y, z : numpy.ndarray
        Position of the particles.
    vx, vy, vz : numpy.ndarray
        Velocity of the particles.
    d : numpy.ndarray
        Density of the particles.  
  
    Returns
    -------
    nw_x, nw_y, nw_z : numpy.ndarray
        Position of the particles recentered using
        the maximum density particle as cm.
    nw_vx, nw_vy, nw_vz : numpy.ndarray
        Velocity of the particles recentered using
        the maximum density particle as cm.
    """
    max_d = np.argmax(d)
    
    xcm = x[max_d]
    ycm = y[max_d]
    zcm = z[max_d]
    
    vcmx = vx[max_d]
    vcmy = vy[max_d]
    vcmz = vz[max_d]

    nw_x = (x-xcm)
    nw_y = (y-ycm)
    nw_z = (z-zcm)

    nw_vx = (vx-vcmx)
    nw_vy = (vy-vcmy)
    nw_vz = (vz-vcmz)
    
    return nw_x, nw_y, nw_z, nw_vx, nw_vy, nw_vz

def re_order(x,y,z,vx,vy,vz,m,u):
    """
    Reordering of all quantities extracted from the ascii
    files, sorting from lower to higher radius.

    Parameters
    ----------
    x, y, z : numpy.ndarray
        Position of the particles.
    vx, vy, vz : numpy.ndarray
        Velocity of the particles.
    m : numpy.ndarray
        Mass of the particles.  
    u : numpy.ndarray
        specific internal energy of the particles. 
  
    Returns
    -------
    ro_x, ro_y, ro_z : numpy.ndarray
        Position of the particles reordered with the
        radial criteria.
    ro_vx, ro_vy, ro_vz : numpy.ndarray
        Velocity of the particles reordered with the
        radial criteria.
    ro_r : numpy.ndarray
        Radius of the particles reordered with the
        radial criteria.
    ro_v : numpy.ndarray
        Velocity norm of the particles reordered with
        the radial criteria.
    ro_m : numpy.ndarray
        Mass of the particles reordered with the radial
        criteria.
    ro_u : numpy.ndarray
        Specific internal energy of the particles reordered
        with the radial criteria.
    """

    r = np.sqrt(x**2 + y**2 + z**2)
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    
    indx = np.argsort(r)
    
    ro_x = x[indx]
    ro_y = y[indx]
    ro_z = z[indx]

    ro_vx = vx[indx]
    ro_vy = vy[indx]
    ro_vz = vz[indx]

    ro_r = r[indx]
    ro_v = v[indx]

    ro_u = u[indx]
    ro_m = m[indx]
    
    return ro_x, ro_y, ro_z, ro_vx, ro_vy, ro_vz, ro_r, ro_v, ro_m, ro_u

def mass_quantities(m):
    """
    Computation of the total mass and enclosed mass profile.

    Parameters
    ----------
    m : numpy.ndarray
        Mass of the particles.  
  
    Returns
    -------
    mt : numpy.float64
        Total mass value.
    m_enc : numpy.ndarray
        Enclosed mass profile.
    """
    
    mt = np.sum(m)
    m_enc = np.cumsum(m)
    
    return mt, m_enc

def energy(v,u,r,m_enc):
    """
    Computation of the specific total energy of each particle.

    Parameters
    ----------
    v : numpy.ndarray
        Velocity norm of the particles reordered
        with the radial criteria.
    u : numpy.ndarray
        Specific internal energy of the particles reordered
        with the radial criteria.
    r : numpy.ndarray
        Radius of the particles reordered with the radial
        criteria.
    m_enc : numpy.ndarray
        Enclosed mass profile.
  
    Returns
    -------
    e : numpy.ndarray
        Specific total energy of the particles reordered
        with the radial criteria.
    """

    G = ((cons.G)/((units.R_sun.to(units.m)**3))*(units.M_sun.to(units.kg))*((1.8845e-2*86400)**2)).value
    
    e = v**2 + u - G*m_enc*(1/r)
    
    return e

def bound_unbound(r,e):
    """
    Definition of the bound and unbound particles based
    on the energy calculation, asuming that the inner
    particles (R < 1.5 [R_sun]) are automatically bound.

    Parameters
    ----------
    r : numpy.ndarray
        Radius of the particles reordered with the radial
        criteria.
    e : numpy.ndarray
        Specific total energy of the particles reordered
        with the radial criteria.
  
    Returns
    -------
    bn : numpy.ndarray
        Index of the bound particles.
    un : numpy.ndarray
        Index of the unbound particles.
    """
    
    inner_r = np.where(r > 1.5)[0][0]
    
    un = np.where(e[inner_r:] > 0)[0] + inner_r
    
    bn_fn = np.where(e[inner_r:] <= 0)[0] + inner_r
    bn_in = np.where(r < r[inner_r])[0]
    bn = np.concatenate((bn_in,bn_fn))
    
    return bn, un

if __name__=="__main__":

    folder = sys.argv[1]  # Folder of each simulation.
    simulation = sys.argv[2] # Path of the ascii file.
    name = sys.argv[3]  # Snapshot name.

    x, y, z, vx, vy, vz, m, d, u = np.genfromtxt('/tank0/ballone/coll_set/{}/{}'.format(folder,simulation),
                                          usecols= (0,1,2,3,4,5,6,7,8), unpack=True)

    # Calculate the position and velocity wrt the center of mass.

    nw_x, nw_y, nw_z, nw_vx, nw_vy, nw_vz = re_center(x,y,z,vx,vy,vz,d)

    # Reordering of all quantities based on Radius sorting from lower to higher.

    X, Y, Z, VX, VY, VZ, R, V, M, U = re_order(nw_x,nw_y,nw_z,nw_vx,nw_vy,nw_vz,m,u)

    # Generation of the enclosed mass profile and total mass.

    Mt, M_enc = mass_quantities(M)

    # Energy calculation for the unbound mass definition.

    e = energy(V,U,R,M_enc)

    # Aplication of bound-unbound criteria.

    bn, un = bound_unbound(R,e)

    # Calculating the porcentage of unbound mass in this snapshot.

    percent = np.sum(M[un])*100/Mt

    # Generation of the enclosed unbound and bound mass profile.

    Mencl_un = np.cumsum(M[un])
    R_un = R[un]
    
    Mencl_bn = np.cumsum(M[bn])
    R_bn = R[bn]

    with open('/tank0/ballone/coll_set/{}/Munb_percent.txt'.format(folder), 'a') as myfile:
              myfile.write(str(round(percent,4))+'\t'+name+'\n')

    fig = plt.figure(figsize=(10, 7))

    plt.plot(R,M_encl,'.k', label='Total mass')
    plt.plot(R_un,Mencl_un,'.r', label='Unbound mass')
    plt.plot(R_bn,Mencl_bn,'.g', label='Bound mass')
    plt.xlabel("$R$ [$R_{\odot}$]")
    plt.ylabel("$M_{encl}$ [$M_{\odot}$]")
    plt.semilogx()
    plt.legend()

    plt.tight_layout()
    plt.savefig('/tank0/ballone/coll_set/{}/'.format(folder)+name+'.png')
