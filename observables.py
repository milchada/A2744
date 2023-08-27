import numpy as np 
import matplotlib.pylab as plt 
from matplotlib import colors, cm
from scipy.optimize import curve_fit
from astropy.io import fits

basePath = '/Users/mila/Documents/Research/Postdoc/A2744/'
furt = fits.getdata(basePath+'Furtak-JWST-lensing/A2744_cluster-members_v1.0 Lukas Furtak.fits')
bcg2 = furt[1]

fig, ax = plt.subplots()
hist, bins = np.histogram(furt['z'], range=(0.28, 0.33), bins=50)

def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def trimodal(x,mu1,sigma1,A1,mu2,sigma2,A2, mu3, sigma3, A3):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x,mu3,sigma3,A3)

p0 = (furt['z'][0], 0.002, 10, furt['z'][1],0.002,10, furt['z'][2], 0.002, 3)
params,cov=curve_fit(trimodal, bins[:-1], hist, p0=p0)
sigma = np.sqrt(cov)
plt.plot(bins[:-1], trimodal(bins[:-1], *params), color='tab:purple', lw=3, label='Model')

#z_fit = np.array([0.31744283, 0.30067435, 0.30753658])
#sigma_fit = np.array([0.00296895, 0.00446729, 0.00094396])
#		   = np.array([ 890.68371263, 1340.18550331,  283.18850712])
#A_fit = np.array([10.06950447, 10.75791065,  3.08782618])

"Finding relative velocities between BCGs is so tricky!!"

import pickle 
from astropy.constants import k_B

kB = k_B.to('keV/K').value

def catalog_regions(dir, npeaks=3):
    #Find 2 or 3 peak points and then find the kappa, sb, temp around them
    
    kappas = glob.glob(dir+'/kappa*fits'); kappas.sort()
    median = {}
    std = {}
    keys = ['kappa', 'sb', 'temp']
    for key in keys:
        median[key] = {}
        std[key] = {}
        for ax in ['x', 'y', 'z']:
            median[key][ax] = {}
            std[key][ax] = {}
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

        for peak in peaks:
            c = (peak[1], peak[0])
            dist = np.sqrt((X - c[0])**2 + (Y - c[1])**2)
            mask = (dist < rad)
            maskind = peaks.tolist().index(peak.tolist())

            for (img, lab) in zip([kappa, sb, temp], keys):
                median[lab][ax][snap][maskind] = np.nanmedian(img[mask])
                std[lab][ax][snap][maskind] = np.nanstd(img[mask])

            print('Peak %d done' % maskind)

    with open('median.pkl', 'wb') as f:
        pickle.dump(median, f)

    with open('std.pkl', 'wb') as f:
        pickle.dump(median, f)
            
    # with open('saved_dictionary.pkl', 'rb') as f:
    #     loaded_dict = pickle.load(f)

#ok this is good. Now what do I do with this? 
#I think I should rank order the peaks by the median kappa 
    #from obs, we have median and std. Honestly with the level of idealisation in this problem 
    #I think the median should just be within the observed median +/- std
