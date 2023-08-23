import glob, os
from scipy.ndimage import rotate
from astropy.io import fits
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.constants import k_B
from astropy import units as u
import numpy as np
from skimage.feature import peak_local_max

z = 0.308
from astropy.constants import c, G
from astropy.cosmology import default_cosmology
lcdm = default_cosmology.get()
Dl = lcdm.angular_diameter_distance(z)
Sigma_crit = c**2/(4*np.pi*G*Dl)

#rotate(image, angle, reshape=False)
dirs = ['ncc-cc-r=1_8/','ncc-wcc-r=3/','wcc-cc-r=6_67/','ncc-ncc-r=5_2/', 'wcc-wcc-r=4/'] #'wcc-cc-r=4/',

kB = k_B.to('keV/K').value
props = ['kappa', 'photon']

def plot_rotation(dir, snap, xmin=100, xmax=400):
    sb = glob.glob(dir+'/xray_phot*%d*theta*phi*fits' % snap)
    temp = glob.glob(dir+'/temp*%d*theta*phi*fits' % snap)
    kappa = glob.glob(dir+'/kappa*%d*theta*phi*fits' % snap)
    sb.sort()
    temp.sort()
    kappa.sort()
    dx = fits.getheader(sb[0])['CDELT1']/1000. #Mpc per pix
    tix = np.arange(0, 14.1, 2)
    tpos = tix/dx
    tix -= 7
    for x in sb:
    	t = x.replace('xray_photon_flux_0.5_7_keV', 'temperature')
    	k = x.replace('xray_photon_flux_0.5_7_keV', 'kappa')
    	try:
	        fig, ax = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(12,4))
	        theta = x.split('theta')[1].split('_')[0]
	        phi = x.split('phi')[1].split('.')[0]

	        xray = fits.getdata(x)        
	        im1 = ax[0].imshow(xray, cmap=cm.magma, norm=colors.LogNorm(1,1e6))
	        plt.colorbar(im1, ax=ax[0], label='photons/cm**2/s', shrink=0.65, aspect=15)
	        plt.text(xmin+100, xmin+100, r'$\phi$=%s$^\circ$, $\theta$=%s$^\circ$' % (phi, theta))
	        kt = fits.getdata(t) * kB 
	        im2 = ax[1].imshow(kt, cmap=cm.afmhot, norm=colors.Normalize(1,22))
	        plt.colorbar(im2, ax=ax[1], label='keV', shrink=0.65, aspect=15)
	        
	        m = fits.getdata(k)
	        im3 = ax[2].imshow(m, cmap=cm.Greys, norm=colors.LogNorm(1e-4, .8))
	        plt.colorbar(im3, ax=ax[2], label = r'g/cm$^2$', shrink=0.65, aspect=15)
	        
	        ax[0].set_ylabel('Y (Mpc)')
	        for a in ax.flatten():
	        	a.set_xticks(tpos, tix)
		        a.set_yticks(tpos, tix)
		        a.set_xlabel('X (Mpc)')
		        a.set_xlim(xmin, xmax)
		        a.set_ylim(xmin, xmax)
	        
	        if len(phi) == 1:
	        	phi = '00'+phi
	        elif len(phi) == 2:
	        	phi = '0'+phi
	        if len(theta) == 1:
	        	theta = '0'+theta
	        plt.tight_layout()
	        fig.savefig('%s_sb_temp_kappa_snap%s_phi%s_theta%s.png' % (dir, snap, phi, theta), dpi=192)
	        plt.close()
	        print(phi, theta, 'done')
    	except FileError:
	    	continue

def plot_distance(dir,key=''):
	sb = glob.glob(dir+'/xray_phot*%s*fits' % key); sb.sort()
	dx = fits.getheader(sb[0])['CDELT1'] #kpc per pix
	dist_array = np.zeros((len(sb), 2))

	for i in range(len(sb)):
		x = sb[i]
		xray = fits.getdata(x)
		xray[np.isnan(xray)] = 0
		xpeak = peak_local_max(xray, min_distance=5, num_peaks=2)
		d = np.linalg.norm(xpeak[1] - xpeak[0])*dx
		snap = int(x.split('_')[-2].split('.')[0])
		dist_array[i] = snap, d
	return dist_array

def plot(dir, prop='kappa', key='', xmin=3000, xmax=9000, ymin=4000, ymax=10000):
	files = glob.glob(dir+'/*%s*%s*fits' % (prop,key))
	dx = fits.getheader(files[0])['CDELT1']
	files.sort()
	xtix = np.arange(xmin/1000., xmax/1000. + 0.1, 0.2)
	xmean = (xmin+xmax)/2000.
	xpos = xtix*1000./dx
	ytix = np.arange(ymin/1000., ymax/1000. + 0.1, 0.2)
	ypos = ytix*1000./dx

	area_pix = (dx*u.kpc.to('cm'))**2
	dL = lcdm.luminosity_distance(z).to('cm').value
	tcolor = 'w'
	for m in files:
		if 'kappa' in prop:
			img = fits.getdata(m)/Sigma_crit.to('g/cm**2').value
			cmap = cm.viridis
			bounds = 10**np.linspace(np.log10(0.086), np.log10(1.8), 10)
			norm = colors.BoundaryNorm(bounds, cmap.N)
		elif 'photon_flux' in prop:
			img = fits.getdata(m)#*area_pix/(4*np.pi*dL**2)
			#there is some unit conversion here that I am missing
			#the sim output is in photons/s/cm**2.. but this didn't know the redshift, did it?
			cmap = cm.afmhot
			norm = colors.PowerNorm(0.5, vmin=1.5e-8, vmax=4.8e-7)
		elif 'temp' in prop:
			img = fits.getdata(m)*kB
			cmap = cm.RdBu_r
			norm = colors.Normalize(5.5, 22)
			tcolor = 'k'
		plt.clf()
		plt.imshow(img, norm=norm, cmap=cmap)
		plt.xlim(xmin/dx, xmax/dx)
		plt.ylim(ymin/dx, ymax/dx)
		plt.colorbar()
		snap = m.split('_')[-2].split('.')[0]
		print(snap)
		plt.xticks(xpos, ['%0.1f' % t for t in (xtix - xmean)])
		plt.yticks(ypos, ['%0.1f' % t for t in (ytix - xmean)])
		plt.text(xpos[0]+10, ypos[0]+10, 't = %0.2f Gyr' % (float(snap)/100.), color=tcolor)
		plt.savefig(dir+'_%s_%s%s.png' % (prop, snap, key), dpi=192)

def plot_xray_mass(dir, nlevels=5, snapmin=0, xmin=0, xmax=512, key='', suffix='', pre=''):
    if pre == '':
        pre == dir
    sb = glob.glob(dir+'/xray_phot*%s*fits' % key)
    mass = glob.glob(dir+'/kappa*%s*fits' % key)
    sb.sort()
    mass.sort()
    dx = fits.getheader(sb[0])['CDELT1']/1000. #Mpc per pix
    tix = np.arange(0, 14.1, 2)
    tpos = tix/dx
    tix -= 7
    for (x, m) in zip (sb, mass):
    	snap = m.split('_')[-2].split('.')[0]
    	if int(snap) >= snapmin:
	        xray = fits.getdata(x)
	        xray[np.isnan(xray)] = 0
	        xpeak = peak_local_max(xray, min_distance=5, num_peaks=2)
	        lens = np.log10(fits.getdata(m))
	        plt.imshow(xray, cmap=cm.magma, norm=colors.LogNorm(1,1e6))
	        plt.colorbar(label='photons/cm**2/s')
	        plt.text(100, 50, 't = %0.1f Gyr' % (float(snap)/100), color='w')
	        for x in xpeak:
	        	plt.scatter(x[1], x[0], color='k', marker='x', linewidths=2)
	        plt.contour(lens, levels = nlevels, colors='w', linewidths = 1)
	        plt.xticks(tpos, tix)
	        plt.yticks(tpos, tix)
	        plt.xlabel('X (Mpc)')
	        plt.ylabel('Y (Mpc)')
	        plt.xlim(xmin, xmax)
	        plt.ylim(xmin, xmax)
	        plt.savefig('%s_xray_lensing_%s%s.png' % (dir, snap, suffix), dpi=192)
	        plt.close()
	        print(snap, 'done')