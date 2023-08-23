# function to compute the velocities given d, b and then theta, phi"
import numpy as np
from astropy.cosmology import default_cosmology
lcdm = default_cosmology.get()

def vel_turnaround(M1, M2, d, b, z, tmerge_ago):
    #Sarazin 2002 eqs 10 and 15
    #M_1, M_2   : in Msun
    #d0, b0     : in Mpc
    #tmerge_ago : in Gyr
    
    tobs = lcdm.age(z).to('Gyr').value
    tmerge = tobs - tmerge_ago
    d0 = 4.5 *((M1 + M2)/1e15)**(1/3) * (tmerge/10)**(2/3) #Mpc

    v = 2930 * np.sqrt((M1 + M2)/1e15) * np.sqrt(1/d) * np.sqrt((1 - d/d0)/(1-(b/d0)**2))
    return v #km/s

def pos_vel_3d(d, b, v, pos1, theta, phi):
    #d, b, pos1 : in Mpc
    #v          : in km/s
    #theta, phi : in degrees

    th = np.deg2rad(theta)
    ph = np.deg2rad(phi)
    n = np.array([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)])

    #the impact parameter is tilted 90deg wrt the merger
    n2 = np.array([np.sin(th+np.pi/2.)*np.cos(ph), np.sin(th+np.pi/2.)*np.sin(ph), np.cos(th+np.pi/2.)])
    
    d_merge = np.sqrt(d**2 - b**2)*np.sign(d)
    pos2 = pos1 + np.dot(n, d_merge)  + np.dot(n2, b)*np.sign(d)
    v3d = np.dot(n, v)
    return pos1*1000., pos2*1000., v3d

#now we can sample different merger orientations theta, phi that are all consistent with the turnaround criterion