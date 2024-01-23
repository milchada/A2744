import cluster_generator as cg
from yt.utilities.cosmology import Cosmology
import matplotlib.pyplot as plt
import numpy as np
from colossus.halo import concentration
from colossus.cosmology import cosmology
import unyt as u

# Redshift for cluster
z = 0.3072 

# Use Colossus cosmology class and create a yt cosmology class
ccosmo = cosmology.setCosmology('planck18')
cosmo = Cosmology(hubble_constant=ccosmo.h, omega_matter=ccosmo.Om0, omega_lambda=ccosmo.Ode0)

# masses for primary cluster
p_rmax = 5000.0

def make_clusters(run_dicts):
    fig, ax = plt.subplots(ncols = 2, sharex=True, figsize=(12, 8))
    colors = {0:'tab:blue', 1:'tab:orange', 2:'tab:red'}
    lss = {0.01:'dotted', 0.1:'dashdot', 0.5:'dashed', 1.0:'solid'}

    for run_dict in run_dicts:
        M200 = run_dict["M200"]
        rcore = run_dict["r_c"]
        alpha = run_dict["alpha"]

        # Use the concentration-mass relation to determine concentrations
        conc = concentration.modelDiemer19(M200/ccosmo.h, z)[0]

        # This helps us find r200 for the two clusters
        r200 = cg.find_overdensity_radius(M200, 200.0, z=z, cosmo=cosmo)

        # Compute parameters for truncated NFW profile
        r_s = r200/conc
        rho_s = cg.nfw_scale_density(conc, z=z)
        r_t = r200*2.0

        # This creates truncated NFW mass profile objects 
        # these are like functions, you can call them, e.g. Mt_1(r)
        Mt = cg.tnfw_mass_profile(rho_s, r_s, r_t)

        # This finds true radii and masses based on the profile
        r500, M500 = cg.find_radius_mass(Mt, z=z, delta=500.0, cosmo=cosmo)
        r2500, M2500 = cg.find_radius_mass(Mt, z=z, delta=2500.0, cosmo=cosmo)

        rhot = cg.tnfw_density_profile(rho_s, r_s, r_t)
        rhos = 0.9*cg.tnfw_density_profile(rho_s, 0.5*r_s, r_t).cutoff(20.0)
        f_g = 0.115
        rhog = cg.vikhlinin_density_profile(1.0, rcore*r2500, 1.1*r200,
                                            alpha, 0.67, 3)
        rhog = cg.rescale_profile_by_mass(rhog, f_g*M500, r500)
        hse = cg.HydrostaticEquilibrium.from_dens_and_tden(0.1, 15000.0, rhog, rhot,
                                                            stellar_density=rhos)
        
        c = colors[alpha]
        ls = lss[rcore]
        ax[0].loglog(hse["radius"], hse["density"], label=r'$r_c$ = %0.1f, $\alpha$ = %d' % (run_dict['r_c'], run_dict['alpha']), color=c, linestyle=ls)
        ax[0].set_title(f"Density")
        ax[0].set_ylabel(r'$g/cm^{-3}$')
        ax[0].set_ylim(1e4, 1e8)
        ax[1].plot(hse["radius"], hse["temperature"], color=c, linestyle=ls)
        ax[1].set_title('Temperature')
        ax[1].set_ylabel('keV')
        plt.xlim(10, 1000)
        plt.xlabel("r (kpc)")
        ax[1].set_ylim(0,13)
    ax[0].legend(ncols=3)
    plt.tight_layout()
    plt.savefig(f"rcore_alpha_%de14Msun.png" % (run_dict['M200']/1e14))


for M200 in [2e14, 4e14, 6e14]:
    run_dicts = []
    for alpha in [0, 1, 2]:
        for r_c in [0.01, 0.1, 0.5, 1.0]:
            run_dicts.append({'r_c': r_c, 'alpha': alpha, 'M200':M200})
            make_clusters(run_dicts)

