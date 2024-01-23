######################################################
## read simulation output and convert to FITS files ##
######################################################
import yt, glob, os, gc
import numpy as np

homedir = os.getcwd()
z = 0.308
from astropy.constants import c, G
from astropy.cosmology import default_cosmology
lcdm = default_cosmology.get()
Dl = lcdm.angular_diameter_distance(z)
Sigma_crit = c**2/(4*np.pi*G*Dl)
print(Sigma_crit.to('g/cm**2').value)

def total_dens(field, data):
        #create a total density field
        return data['gas', 'density'] + data['particle_density_on_grid']

def make_fits(file, mass = True, sb=True, temp=True, width=(10, 'Mpc'), suffix='', proj=True):
    thetas = np.arange(0, 91, 30)
    phis = np.arange(0, 181, 60)
    ds = yt.load(file, default_species_fields="ionized")
    snapmin = 0
    snapmax = 300
    filenum = int(1000*ds.current_time)
    print(filenum)
    if (filenum > snapmin) and (filenum < snapmax):
        pre=''
        if filenum < 10:
            pre = '000'
        elif filenum < 100:
            pre = '00'
        elif filenum < 1000:
            pre = '0'
        if proj:
            type = 'proj'
        else:
            type = 'slice'
        for theta in thetas:
            for phi in phis:
                print(theta, phi)
                n = [np.cos(np.deg2rad(theta))*np.sin(np.deg2rad(phi)), np.sin(np.deg2rad(theta))*np.sin(np.deg2rad(phi)), np.cos(np.deg2rad(phi))]
                suffix = '_theta%d_phi%d' % (theta, phi)
                if mass:
                    if not glob.glob("kappa_%s_%s%d%s.fits" % (type,pre, filenum, suffix)):
                        ds.add_field(("gamer", "total_density"), units="g/cm**3", function=total_dens, sampling_type="cell")
                        prj_fits = yt.FITSOffAxisProjection(ds, n, ('gamer','total_density'), weight_field=None, width=width, center=('min', 'Pote'))
                        prj_fits.writeto("kappa_%s_%s%d%s.fits" % (type,pre,filenum,suffix))
                        print("kappa done")
                        del(prj_fits)
                        gc.collect()

                if temp:
                    if not glob.glob("temperature_%s_%s%d%s.fits" % (type,pre, filenum, suffix)):
                        if proj:
                            prj_fits = yt.FITSOffAxisProjection(ds, n, ('gas','temperature'), weight_field='mazzotta_weighting', width=width, center=('min', 'Pote'))
                        else:
                            prj_fits = yt.FITSOffAxisSlice(ds, n, ('gas','temperature'),width=width, center=('min', 'Pote'))
                        prj_fits.writeto("temperature_%s_%s%d%s.fits" % (type,pre,filenum,suffix))
                        print("Temp done")
                        del(prj_fits)
                        gc.collect()

                if sb:
                    if not glob.glob("xray_photon_flux_0.5_7_keV_%s_%s%d%s.fits" % (type,pre,filenum,suffix)):
                        xray_fields = yt.add_xray_emissivity_field(ds, 0.5, 7, table_type='apec', metallicity=0.3)
                        prj1 = yt.FITSOffAxisProjection(ds, n, ('gas','xray_photon_emissivity_0.5_7_keV'), weight_field=None, width=width, center=('min', 'Pote'))
                        prj2 = yt.FITSOffAxisProjection(ds, n, ('gas','xray_emissivity_0.5_7_keV'),weight_field=None, width=width, center=('min', 'Pote'))
                        prj1.writeto("xray_photon_flux_0.5_7_keV_%s_%s%d%s.fits" % (type,pre,filenum,suffix), overwrite=True)
                        prj2.writeto("xray_sb_0.5_7_keV_%s_%s%d%s.fits" % (type,pre,filenum,suffix), overwrite=True)
                        print(" Xray SB done")
                        del(prj1, prj2)
                        gc.collect()

            print( 'Done for snap', filenum)

def parallelize_make_fits(snapmin=0, snapmax=1210):
    #run make_fits in parallel 
    from multiprocessing import Pool
    import itertools

    files = glob.glob('Data*')
    files.sort()

    if files:
        pool = Pool(processes=8)
        pool.map(make_fits, files)
        pool.close()
        pool.join()

if __name__ == '__main__':
    parallelize_make_fits()
  