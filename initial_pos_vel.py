# function to compute the velocities given d, b and then theta, phi"
import numpy as np
from astropy.cosmology import default_cosmology
import glob, os

lcdm = default_cosmology.get()

def vel_turnaround(M1, M2, d, b, z, tmerge_ago, delta_dex = 0):
    #Sarazin 2002 eqs 10 and 15
    #M_1, M_2   : in Msun
    #d0, b0     : in Mpc
    #tmerge_ago : in Gyr
    #delta_dex  : <= 0.2
    
    tobs = lcdm.age(z).to('Gyr').value
    tmerge = tobs - tmerge_ago
    d0 = 4.5 *((M1 + M2)/1e15)**(1/3) * (tmerge/10)**(2/3) #Mpc

    v = 2930 * np.sqrt((M1 + M2)/1e15) * np.sqrt(1/d) * np.sqrt((1 - d/d0)/(1-(b/d0)**2))
    return 10**(np.log10(v) + delta_dex) # km/s

def pos_vel_3d(d, b, v, pos1, theta, phi):
    #d, b, pos1 : in Mpc
    #v          : in km/s
    #theta, phi : in degrees

    th = np.deg2rad(theta)
    ph = np.deg2rad(phi)
    n = np.array([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)])

    #the impact parameter is tilted 90deg wrt the merger
    n2 = np.array([np.sin(th+np.pi)*np.cos(ph+np.pi/2), np.sin(th+np.pi)*np.sin(ph+np.pi/2), np.cos(th+np.pi)])
    
    d_merge = np.sqrt(d**2 - b**2)*np.sign(d)
    pos2 = pos1 + np.dot(n, d_merge)  + np.dot(n2, b)*np.sign(b)
    v3d = np.dot(n, v)
    return pos1*1000., pos2*1000., v3d

#now we can sample different merger orientations theta, phi that are all consistent with the turnaround criterion
#for NW1-S2-N
#first S2-N. S2 at center of box, N comes in from the south, NW1 from the southeast
v = vel_turnaround(2e14, 8e14, 3, 0.5, 0.308, 1, 0.2) #because N feels the mass of S2+NW1
pos1, pos2, vel2 = pos_vel_3d(-3, 0.5, v, np.array([7.,7.,7.]), 90,90)
print('S2-N: ', pos2.astype(int), vel2.astype(int))
#then NW1 wrt S2, similar start time but at an angle to line of sight
for d in [1, 2, 3]:
    print(d)
    v = vel_turnaround(2e14, 8e14, d, 0.1, 0.308, 1, 0.2) #because N feels the mass of S2+NW1
    for theta in [60, 30]:
        for phi in [15, 45, 75]:
            pos1, pos2, vel2 = pos_vel_3d(-d, 0.1, v, np.array([7.,7.,7.]), theta, phi)
            print(theta, phi, pos2.astype(int), vel2.astype(int))

def input_file(dir, pos2, vel2, pos3, vel3, pre1='s2', pre2='nw1', pre3='n'):
    geom = "#initial conditions \n \
    Merger_Coll_NumHalos            3       # number of halos \n \
    Merger_Coll_UseMetals           1                         \n \
                                                              \n    \
    #S2, initialised at center and rest frame of box \n\
    Merger_File_Prof1               %s_profile_gamer.h5     # profile table of cluster 1 \n\
    Merger_File_Par1                %s_gamerp_particles.h5  # particle file of cluster 1 \n\
    Merger_Coll_PosX1               7000.0                  # X-center of cluster 1 in kpc \n\
    Merger_Coll_PosY1               7000.0                  # Y-center of cluster 1 in kpc \n\
    Merger_Coll_PosZ1               7000.0                  # Z-center of cluster 1 in kpc \n\
    Merger_Coll_VelX1               0.0                     # X-velocity of cluster 1 in km/s \n\
    Merger_Coll_VelY1               0.0                     # Y-velocity of cluster 1 in km/s \n\
    Merger_Coll_VelZ1               0.0                     # Z-velocity of cluster 1 in km/s \n\
    \n\
    #NW1 Cluster, 4.7 kpc east of S2 at start, inclined 54º to line of sight \n\
    Merger_File_Prof2               %s_profile_gamer.h5    # profile table of cluster 2 \n\
    Merger_File_Par2                %s_gamerp_particles.h5 # particle file of cluster 2 \n\
    Merger_Coll_PosX2               %0.1f                  # X-center of cluster 2 in kpc \n\
    Merger_Coll_PosY2               %0.1f                  # Y-center of cluster 2 in kpc \n\
    Merger_Coll_PosZ2               %0.1f                  # Z-center of cluster 2 in kpc \n\
    Merger_Coll_VelX2               %0.1f                  # X-velocity of cluster 2 in km/s \n\
    Merger_Coll_VelY2               %0.1f                  # Y-velocity of cluster 2 in km/s \n\
    Merger_Coll_VelZ2               %0.1f                  # Z-velocity of cluster 2 in km/s \n\
    \n\
    #N Cluster, 2 Mpc south of S2 at start, moving 45º NW \n\
    Merger_File_Prof3               %s_profile_gamer.h5    # profile table of cluster 3 \n\
    Merger_File_Par3                %s_gamerp_particles.h5 # particle file of cluster 3 \n\
    Merger_Coll_PosX3               %0.1f                  # X-center of cluster 3 in kpc \n\
    Merger_Coll_PosY3               %0.1f                  # Y-center of cluster 3 in kpc \n\
    Merger_Coll_PosZ3               %0.1f                  # Z-center of cluster 3 in kpc \n\
    Merger_Coll_VelX3               %0.1f                  # X-velocity of cluster 3 in km/s \n\
    Merger_Coll_VelY3               %0.1f                  # Y-velocity of cluster 3 in km/s \n\
    Merger_Coll_VelZ3               %0.1f                  # Z-velocity of cluster 3 in km/s \
    " % (pre1, pre1, pre2, pre2, *pos2, *vel2, pre3, pre3, *pos3, *vel3)
    
    os.chdir(dir)
    with open('Input__TestProb', 'w') as f:
        f.writelines(geom)

#now NW1-S2-NW2
#first NW1-NW2. NW1 at center of box, NW2 comes in from the east, S2 from the northwest
v = vel_turnaround(2e14, 8e14, 3, 0.5, 0.308, 0, 0.2) #NW2 hasn't reached pericenter yet
pos1, pos2, vel2 = pos_vel_3d(-3, 0.5, v, np.array([7.,7.,7.]), 90,0)
print('NW1-NW2: ', pos2.astype(int), vel2.astype(int))
#then NW1 wrt S2, similar start time but at an angle to line of sight
pwd = os.getcwd()
for d in [1, 2, 3]:
    if not glob.glob('d2=%dMpc' % d):
        os.mkdir('d2=%dMpc' % d)
    v = vel_turnaround(2e14, 8e14, d, 0.1, 0.308, 1, 0.2) #because N feels the mass of S2+NW1
    for theta in [60, 30]:
        for phi in [15, 45, 75]:
            if not glob.glob('d2=%dMpc/theta=%d_phi=%d' % (d, theta, phi)):
                os.mkdir('d2=%dMpc/theta=%d_phi=%d' % (d, theta, phi))
            pos1, pos3, vel3 = pos_vel_3d(d, 0.1, -v, np.array([7.,7.,7.]), theta, phi)
            dirname = 'd2=%dMpc/theta=%d_phi=%d' % (d, theta, phi)
            input_file(dirname, pos2, vel2, pos3, vel3)
            os.chdir(pwd)