import numpy as np 
import matplotlib.pylab as plt 
from matplotlib import colors, cm
from scipy.optimize import curve_fit
from astropy.io import fits
import pickle, yt, glob
from astropy.constants import k_B
from skimage.feature import peak_local_max

kB = k_B.to('keV/K').value

def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def trimodal(x,mu1,sigma1,A1,mu2,sigma2,A2, mu3, sigma3, A3):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x,mu3,sigma3,A3)

def model_obs_velocities():
    basePath = '/Users/mila/Documents/Research/Postdoc/A2744/'
    furt = fits.getdata(basePath+'Furtak-JWST-lensing/A2744_cluster-members_v1.0 Lukas Furtak.fits')
    bcg2 = furt[1]
    fig, ax = plt.subplots()
    hist, bins = np.histogram(furt['z'], range=(0.28, 0.33), bins=50)
    p0 = (furt['z'][0], 0.002, 10, furt['z'][1],0.002,10, furt['z'][2], 0.002, 3)
    params,cov=curve_fit(trimodal, bins[:-1], hist, p0=p0)
    sigma = np.sqrt(cov)
    plt.plot(bins[:-1], trimodal(bins[:-1], *params), color='tab:purple', lw=3, label='Model')

#z_fit = np.array([0.31744283, 0.30067435, 0.30753658])
#sigma_fit = np.array([0.00296895, 0.00446729, 0.00094396])
#		   = np.array([ 890.68371263, 1340.18550331,  283.18850712])
#A_fit = np.array([10.06950447, 10.75791065,  3.08782618])

def catalog_regions(dir, npeaks=3):
    #Find 2 or 3 peak points and then find the kappa, sb, temp around them
    #this would be easiest to do in projection. you can find the 2-3 lensing peaks

    kappas = glob.glob(dir+'/kappa*fits'); kappas.sort()
    median = {}
    std = {}
    pos = {}
    keys = ['kappa', 'sb', 'temp']
    for key in keys:
        median[key] = {}
        std[key] = {}
        for ax in ['x', 'y', 'z']:
            median[key][ax] = {}
            std[key][ax] = {}
            pos[ax] = {}
            
    for k in kappas:
        snap = int(k.split('proj_')[1].split('_')[0])
        ax = k.split('_')[-1].split('.')[0]
        print(snap, ax)

        kappa = fits.getdata(k)
        sb = fits.getdata(k.replace('kappa', 'xray_photon_flux_0.5_7_keV'))
        temp = fits.getdata(k.replace('kappa', 'temperature')) * kB

        peaks = peak_local_max(kappa, num_peaks=npeaks)
        dx = 4000./512
        rad = 100/dx #kpc to pixels
        # masks = []
        X, Y = np.meshgrid(np.arange(kappa.shape[0]), np.arange(kappa.shape[1]))
        for (img, lab) in zip([kappa, sb, temp], keys):
            median[lab][ax][snap] = {}
            std[lab][ax][snap] = {}
            pos[ax][snap] = {}

        for peak in peaks:
            c = (peak[1], peak[0])
            dist = np.sqrt((X - c[0])**2 + (Y - c[1])**2)
            mask = (dist < rad)
            maskind = peaks.tolist().index(peak.tolist())
            pos[ax][snap][maskind] = c*dx

            for (img, lab) in zip([kappa, sb, temp], keys):
                median[lab][ax][snap][maskind] = np.nanmedian(img[mask])
                std[lab][ax][snap][maskind] = np.nanstd(img[mask])

            print('Peak %d done' % maskind)

    with open('median.pkl', 'wb') as f:
        pickle.dump(median, f)

    with open('std.pkl', 'wb') as f:
        pickle.dump(std, f)

    with open('pos.pkl', 'wb') as f:
            pickle.dump(pos, f)
            
    # with open('saved_dictionary.pkl', 'rb') as f:
    #     loaded_dict = pickle.load(f)

#ok this is good. Now what do I do with this? 
#I think I should rank order the peaks by the median kappa 
    #from obs, we have median and std. Honestly with the level of idealisation in this problem 
    #I think the median should just be within the observed median +/- std

def total_dens(field, data):
        #create a total density field
        return data['gas', 'density'] + data['particle_density_on_grid']
    
def find_bcg_velocities(snapshot, npeaks=2, min_distance=150, debug=False):
    ds = yt.load(snapshot)
    boxsize = (ds.parameters['BoxSize'][0], 'Mpc')
    #maximise FITS image resolution
    xmin = ds.find_min(("gamer","dx"))
    xmin = xmin[0].to('Mpc').value
    image_res = round(boxsize[0]/xmin)
    kappa = yt.FITSSlice(ds,'z', ('gamer','ParDens'), width=(boxsize), center='c', image_res = image_res)
    dx = kappa['ParDens'].header['CDELT1']/1000. #kpc to Mpc
    kappa = kappa['ParDens'].data
    peaks = peak_local_max(kappa, num_peaks = npeaks, min_distance=round(min_distance/(dx*1000.)))
    print(peaks)
    vpeaks = np.zeros((npeaks, 3))
    ppeaks = np.zeros((npeaks, 3))
    for i in range(len(peaks)):
        peak = peaks[i]
        c = [peak[1]*dx, peak[0]*dx, boxsize[0]/2.]
        box = ds.sphere(c, dx)
        v_box = [box[('nbody','ParVel%s' % axis)] for axis in ['X', 'Y', 'Z']]
        v_los = [np.mean(v_box[i].to('km/s').value) for i in range(3)]
        vpeaks[i] = v_los
        ppeaks[i] = c
    return vpeaks, ppeaks

def find_all_velocities(min_distance=150, save=True):
    snapshots = glob.glob('Data*')
    snapshots = [snap for snap in snapshots if 'png' not in snap]
    snapshots.sort()
    vpeaks = np.zeros((len(snapshots), 2, 3))
    ppeaks = np.zeros((len(snapshots), 2, 3))
    for i in range(len(snapshots)):
        print(snapshots[i])
        vpeaks[i], ppeaks[i] = find_bcg_velocities(snapshots[i], 2, min_distance=min_distance)
        if save:
            np.save('vpeaks.npy', vpeaks)
            np.save('ppeaks.npy', ppeaks)
        else:
            continue
    if not save:
        return vpeaks, ppeaks
    
def plot_vel_pos(dir = './', dt = 0.1, plot=True, ret=True):
    files = glob.glob(dir+'Data*')
    files.sort()
    times = np.array([int(file.split('_')[1]) for file in files]) * dt
    vpeaks = np.load(dir+'vpeaks.npy')
    ppeaks = np.load(dir+'ppeaks.npy')
    distance = np.sqrt(np.sum((ppeaks[:,0] - ppeaks[:,1])**2, axis=1))
    speed = np.sqrt(np.sum((vpeaks[:,0] - vpeaks[:,1])**2, axis=1))
    if plot:
        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        ax.plot(times, distance, color='tab:blue')
        ax2.plot(times, speed, color='tab:orange')
        ax.set_xlabel('Time (Gyr)')
        ax.set_ylabel('Distance (Mpc)', color='tab:blue')
        ax2.set_ylabel('Speed (km/s)', color='tab:orange')
        plt.tight_layout()
        plt.savefig(dir.split('/')[0]+'_vel_pos.png', dpi=192)
    if ret:
        return times, distance, speed

#how to match observations?
#s1-s2 is along the line of sight
#s2-n is plane of the sky
#s2-nw1 is inclined to the plane of the sky, angle TBD
#for s1-s2, we need (theta, phi) such that d_proj = d_3D * sin(theta) ~ 117 kpc
    #and v_los = v_3D * cos(theta) ~ 5640 km/s
    #current v_3D_max ~ 2600 km/s for v_i = 2000 km/s
    #so we're asking for v_i > 4300 km/s! That's so much!
    #t - t_peri ~ 0.1 Gyr. But with twice the velocity, it's 0.05 Gyr.

#s2-n has d_proj ~ 740 kpc, v_los ~ 370 km/s
    #at d_3D = 1Mpc, v_3D ~1000 km/s, i.e. theta = arccos(0.37) = 68 deg
    #so d_3D = 740/sin(68) ~ 800 kpc. Ok, pretty close. 
    #this is at t-t_peri ~ 0.2 Gyr. A result!!
    #although the cores need to be more disrupted, so b needs to be smaller. Try 150 kpc. 

#s2-nw1 has d_proj ~ 603 kpc, v_los ~ 2827 km/s
    #for such a high v_los, this merger has to be mostly along the line of sight
    #but for such a far separation, not entirely.
    #600 < d_proj < 1250 kpc for 0.3 < t-t_peri < 0.5 Gyr
    #so 0 < theta < 30 deg 
    #phi ~ 90 to make y-offset 0
    #2400 > v > 2300 km/s at these times
    #v_los = v_3D*0.87 --> v_3D ~ 2800 km/s
    #so v_i just has to be ~2300 km/s rather than 2000 km/s. 
    #this is at t-t_peri ~ 0.3-0.5 Gyr. But with 15% higher velocity, it's 0.3-0.4 Gyr.

#so we have a merger picture for s1, s2, n and nw1. let's set up these ICs. 
#Start s2-nw1 where they are now. 
    #theta = 30 deg, phi=?
#Start s2-n such that they get to their current position in 0.2 Gyr
    #theta = 68 deg. phi=? 
#Start s1-s2 such that they get to their current position in 0.4 Gyr 
    #theta ~ 0, phi ~ 0 
