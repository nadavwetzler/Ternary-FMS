
# Code author: nadav Wetzler
# Geological Survey of Israel, Jerusalem
# E-mail: nadavw@gsi.gov.il
# Revision 01/2025: The first version of the function.

import numpy as np
import matplotlib.pyplot as plb
from obspy.imaging.beachball import beachball
from obspy.imaging.beachball import beach
from matplotlib import path


def inpolygon(xq, yq, xv, yv):
    shape = xq.shape
    xq = xq.reshape(-1)
    yq = yq.reshape(-1)
    xv = xv.reshape(-1)
    yv = yv.reshape(-1)
    q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
    p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
    return p.contains_points(q).reshape(shape)


def convert2tri(strike, dip, rake):
    strike = np.array(strike) * np.pi / 180.
    dip    = np.array(dip) * np.pi / 180.
    rake   = np.array(rake) * np.pi / 180.

    N         = len(strike)
    n0        = np.zeros((N,3))
    u0        = np.zeros((N,3))
    P_azimuth = np.zeros(N)
    T_azimuth = np.zeros(N)
    P_theta   = np.zeros(N)
    T_theta   = np.zeros(N)
    dT        = np.zeros(N)
    dP        = np.zeros(N)
    dB        = np.zeros(N)

    # Normals and slip vectors
    n0[:,0] = -np.sin(dip) * np.sin(strike)
    n0[:,1] =  np.sin(dip) * np.cos(strike)
    n0[:,2] = -np.cos(dip)
    u0[:,0] =  np.cos(rake) * np.cos(strike) + np.cos(dip) * np.sin(rake) * np.sin(strike)
    u0[:,1] =  np.cos(rake) * np.sin(strike) - np.cos(dip) * np.sin(rake) * np.cos(strike)
    u0[:,2] = -np.sin(rake) * np.sin(dip)

    # PT-axes
    P_osa = (n0-u0) / np.rot90(np.tile(np.sqrt(np.sum((n0-u0)**2,axis=1)),(3,1)))
    T_osa = (n0+u0) / np.rot90(np.tile(np.sqrt(np.sum((n0+u0)**2,axis=1)),(3,1)))
    P_osa[P_osa[:,2]>0,:] = -P_osa[P_osa[:,2]>0,:]
    T_osa[T_osa[:,2]>0,:] = -T_osa[T_osa[:,2]>0,:]

    # Compute all angles
    for i in range(N):
        # Get azimuths and dip angles
        P_azimuth[i] = np.arctan2(P_osa[i,0],P_osa[i,1])
        P_theta[i]   = np.arccos(np.abs(P_osa[i,2]))

        T_azimuth[i] = np.arctan2(T_osa[i,0],T_osa[i,1])
        T_theta[i] = np.arccos(np.abs(T_osa[i,2]))

        # Get mechanism class
        dT[i] = np.pi/2 - T_theta[i]
        dP[i] = np.pi/2 - P_theta[i]
        dB[i] = np.arcsin(np.real(np.sqrt(1 - np.sin(dT[i])**2 - np.sin(dP[i])**2)))

    return dP*180/np.pi, dT*180/np.pi, dB*180/np.pi


def deg2hv(dP,dT,dB,dN):
    z = np.arctan(np.sin(dT)/np.sin(dP))-np.deg2rad(45)
    try:
        if np.isnan(z)==True:
            z = np.arctan(1) - np.deg2rad(45)
    except:
        pass

    h1 = (np.cos(dB)*np.sin(z)) / (np.sin(dN)*np.sin(dB)+np.cos(dN)*np.cos(dB)*np.cos(z))
    v1 = (np.cos(dN)*np.sin(dB)-np.sin(dN)*np.cos(dB)*np.cos(z))/(np.sin(dN)*np.sin(dB)+np.cos(dN)*np.cos(dB)*np.cos(z))
    return z, h1, v1


def TernaryFMS(dP0, dT0, dB0, ax=False):
    dP0 = np.deg2rad(dP0)
    dT0 = np.deg2rad(dT0)
    dB0 = np.deg2rad(dB0)
    dN = np.deg2rad(35.26)

    def mk_circle_r(x0, y0, r):
        theta = np.linspace(0, 2 * np.pi, 1000)  # Generate angles from 0 to 2π
        x = x0 + r * np.cos(theta)  # X-coordinates of the circle
        y = y0 + r * np.sin(theta)
        return x, y

    # Border
    zero0 = 1E-10
    [z, h1, v1] = deg2hv(np.deg2rad(90), zero0, zero0, dN)
    [z, h2, v2] = deg2hv(zero0, np.deg2rad(90), zero0, dN)
    [z, h3, v3] = deg2hv(zero0, zero0, np.deg2rad(90), dN)

    xc1, yc1 = mk_circle_r(h1, v1, 1.2)
    xc2, yc2 = mk_circle_r(h2, v2, 1.2)
    xc3, yc3 = mk_circle_r(h3, v3, 1.2)

    tri_px = np.array([h1, h2, h3, h1])
    tri_py = np.array([v1, v2, v3, v1])

    I1 = inpolygon(xc1, yc1, tri_px, tri_py)
    I2 = inpolygon(xc2, yc2, tri_px, tri_py)
    I3 = inpolygon(xc3, yc3, tri_px, tri_py)
    
    xc1i = np.append(xc1[I1], [h1, xc1[I1][0]])
    yc1i = np.append(yc1[I1], [v1, yc1[I1][0]])
    
    xc2i = np.append(xc2[I2], [h2, xc2[I2][0]])
    yc2i = np.append(yc2[I2], [v2, yc2[I2][0]])
    
    xc3i = np.append(xc3[I3], [h3, xc3[I3][0]])
    yc3i = np.append(yc3[I3], [v3, yc3[I3][0]])
    
    [z, h, v] = deg2hv(dP0, dT0, dB0, dN)

    Im1 = inpolygon(h, v, xc1i, yc1i)
    Im2 = inpolygon(h, v, xc2i, yc2i)
    Im3 = inpolygon(h, v, xc3i, yc3i)

    # Grid
    
    mClass = np.zeros(len(h))

    mClass[Im1] = 2
    mClass[Im2] = 3
    mClass[Im3] = 1

    Im4 = ~(Im1 + Im2 + Im3)

    if ax:
        dsz = 30
        ax.plot([h1, h2, h3, h1], [v1, v2, v3, v1], '-k')
        ax.plot(xc1[I1], yc1[I1], 'k')
        ax.plot(xc2[I2], yc2[I2], 'k')
        ax.plot(xc3[I3], yc3[I3], 'k')
        ax.scatter(h[Im1], v[Im1], dsz, 'g')  # Normal
        ax.scatter(h[Im2], v[Im2], dsz, 'b')  # Thrust
        ax.scatter(h[Im3], v[Im3], dsz, 'r')  # Strike-Slip
        ax.scatter(h[Im4], v[Im4], dsz, 'grey')

        dv = 0.2
        beach1 = beach([0, 45, -90], facecolor='g', edgecolor='k', xy=(h1, v1), width=0.20, linewidth=0.2, zorder=3)
        ax.add_collection(beach1)
        if sum(mClass == 2) > 0:
            ax.text(h1-dv, v1-dv, '%2.1f %s' % (sum(mClass == 2)/len(mClass)*100, '%'))
        else:
            ax.text(h1-dv, v1-dv, '0')

        beach2 = beach([0, 45, 90], facecolor='b', edgecolor='k', xy=(h2, v2), width=0.20, linewidth=0.2, zorder=3)
        ax.add_collection(beach2)
        if sum(mClass == 3) > 0:
            ax.text(h2, v2-dv, '%2.1f %s' % (sum(mClass==3)/len(mClass)*100, '%'))
        else:
            ax.text(h2, v2-dv, '0')

        beach3 = beach([45, 90, 180], facecolor='r', edgecolor='k', xy=(h3, v3), width=0.20, linewidth=0.2, zorder=3)
        ax.add_collection(beach3)
        if sum(mClass == 1) > 0:
            ax.text(h3, v3+dv, '%2.1f %s' % (sum(mClass==1)/len(mClass)*100, '%'))
        else:
            ax.text(h3, v3+dv, '0')

        ax.text(0, 0, '%2.1f' % (sum(mClass == 0)/len(mClass)*100))

        ax.set_aspect('equal', adjustable='box')

        ax.grid(False)
        # Hide axes ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
    return mClass