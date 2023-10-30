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

def make_clusters(run_dict):
    num = run_dict["num"]
    R = run_dict["R"]
    M200 = run_dict["M200"]
    M200_1 = M200*R/(1.0+R)
    M200_2 = M200/(1.0+R)
    
    # Use the concentration-mass relation to determine concentrations
    conc_1 = concentration.modelDiemer19(M200_1/ccosmo.h, z)[0]
    conc_2 = concentration.modelDiemer19(M200_2/ccosmo.h, z)[0]

    # This helps us find r200 for the two clusters
    r200_1 = cg.find_overdensity_radius(M200_1, 200.0, z=z, cosmo=cosmo)
    r200_2 = cg.find_overdensity_radius(M200_2, 200.0, z=z, cosmo=cosmo)

    # Compute parameters for truncated NFW profile
    r_s1 = r200_1/conc_1
    r_s2 = r200_2/conc_2
    rho_s1 = cg.nfw_scale_density(conc_1, z=z)
    rho_s2 = cg.nfw_scale_density(conc_2, z=z)
    r_t1 = r200_1*2.0
    r_t2 = r200_2*2.0

    # This creates truncated NFW mass profile objects 
    # these are like functions, you can call them, e.g. Mt_1(r)
    Mt_1 = cg.tnfw_mass_profile(rho_s1, r_s1, r_t1)
    Mt_2 = cg.tnfw_mass_profile(rho_s2, r_s2, r_t2)

    # This finds true radii and masses based on the profile
    r500_1, M500_1 = cg.find_radius_mass(Mt_1, z=z, delta=500.0, cosmo=cosmo)
    r2500_1, M2500_1 = cg.find_radius_mass(Mt_1, z=z, delta=2500.0, cosmo=cosmo)
    r500_2, M500_2 = cg.find_radius_mass(Mt_2, z=z, delta=500.0, cosmo=cosmo)
    r2500_2, M2500_2 = cg.find_radius_mass(Mt_2, z=z, delta=2500.0, cosmo=cosmo)
    rhot_1 = cg.tnfw_density_profile(rho_s1, r_s1, r_t1)
    rhot_2 = cg.tnfw_density_profile(rho_s2, r_s2, r_t2)
    rhos_1 = 0.9*cg.tnfw_density_profile(rho_s1, 0.5*r_s1, r_t1).cutoff(20.0)
    rhos_2 = 0.9*cg.tnfw_density_profile(rho_s2, 0.5*r_s2, r_t2).cutoff(20.0)
    f_g = 0.115
    rhog_1 = cg.vikhlinin_density_profile(1.0, run_dict["r_c1"]*r2500_1, 1.1*r200_1,
                                          run_dict["alpha1"], 0.67, 3)
    rhog_1 = cg.rescale_profile_by_mass(rhog_1, f_g*M500_1, r500_1)
    rhog_2 = cg.vikhlinin_density_profile(1.0, run_dict["r_c2"]*r2500_2, 1.1*r200_2, 
                                          run_dict["alpha2"], 0.67, 3)
    rhog_2 = cg.rescale_profile_by_mass(rhog_2, f_g*M500_2, r500_2)
    hse_1 = cg.HydrostaticEquilibrium.from_dens_and_tden(0.1, 15000.0, rhog_1, rhot_1, 
                                                        stellar_density=rhos_1)
    hse_2 = cg.HydrostaticEquilibrium.from_dens_and_tden(0.1, 15000.0, rhog_2, rhot_2, 
                                                        stellar_density=rhos_2)
    plt.figure(figsize=(10,10))
    plt.loglog(hse_1["radius"], hse_1["density"])
    plt.loglog(hse_2["radius"], hse_2["density"])
    plt.savefig(f"gas_density_{num}.png")
    plt.figure(figsize=(10,10))
    plt.loglog(hse_1["radius"], hse_1["dark_matter_density"])
    plt.loglog(hse_2["radius"], hse_2["dark_matter_density"])
    plt.loglog(hse_1["radius"], hse_1["stellar_density"])
    plt.loglog(hse_2["radius"], hse_2["stellar_density"])
    plt.savefig(f"particle_density_{num}.png")
    plt.figure(figsize=(10,10))
    plt.loglog(hse_1["radius"], hse_1["total_mass"], label="Total")
    plt.loglog(hse_1["radius"], hse_1["dark_matter_mass"], label="DM")
    plt.loglog(hse_1["radius"], hse_1["gas_mass"], label="Gas")
    plt.loglog(hse_1["radius"], hse_1["stellar_mass"], label="Stars")
    plt.loglog(hse_2["radius"], hse_2["total_mass"], label="Total")
    plt.loglog(hse_2["radius"], hse_2["dark_matter_mass"], label="DM")
    plt.loglog(hse_2["radius"], hse_2["gas_mass"], label="Gas")
    plt.xlim(1, 5000)
    plt.xlabel("r (kpc)")
    plt.ylabel("M (M$_\odot$)")
    plt.legend()
    plt.savefig(f"masses_{num}.png")
    masses1 = hse_1.mass_in_radius(p_rmax)
    masses2 = hse_2.mass_in_radius(p_rmax)
    f_star = (masses1["stellar"]+masses2["stellar"])/(masses1["dark_matter"]+masses2["dark_matter"])
    plt.figure(figsize=(10,10))
    plt.plot(hse_1["radius"], hse_1["gas_fraction"])
    plt.plot(hse_2["radius"], hse_2["gas_fraction"])
    plt.xscale("log")
    plt.xlim(1, 5000)
    plt.xlabel("r (kpc)")
    plt.ylabel("f$_{gas}$")
    plt.axhline(f_g)
    plt.legend()
    plt.savefig(f"gas_fraction_{num}.png")
    plt.figure(figsize=(10,10))
    plt.plot(hse_1["radius"], hse_1["temperature"])
    plt.plot(hse_2["radius"], hse_2["temperature"])
    plt.xscale('log')
    plt.xlim(1, 5000)
    plt.xlabel("r (kpc)")
    plt.ylabel("kT (keV)")
    plt.legend()
    plt.savefig(f"temperature_{num}.png")
    plt.figure(figsize=(10,10))
    plt.loglog(hse_1["radius"], hse_1["electron_number_density"])
    plt.loglog(hse_2["radius"], hse_2["electron_number_density"])
    plt.xlim(1,5000)
    #plt.ylim(1.0e-3, None)
    plt.xlabel("r (kpc)")
    plt.ylabel("n$_e$ (cm$^{-3}$)")
    plt.savefig(f"electron_density_{num}.png")
    plt.figure(figsize=(10,10))
    plt.loglog(hse_1["radius"], hse_1["entropy"])
    plt.loglog(hse_2["radius"], hse_2["entropy"])
    plt.xlim(1,5000)
    plt.xlabel("r (kpc)")
    plt.ylabel("S (keV cm$^{2}$)")
    plt.savefig(f"entropy_{num}.png")
    hse_1.write_model_to_h5(f"run{num}_profile_1.h5", overwrite=True)
    hse_2.write_model_to_h5(f"run{num}_profile_2.h5", overwrite=True)
    d = run_dict['d']
    b = run_dict["b"]
    vel = run_dict["v"]
    center = np.array([5000.0]*3)
    center1, center2 = cg.compute_centers_for_binary(center, d, b)
    velocity = np.array([(vel*u.km/u.s).to_value("kpc/Myr"), 0.0, 0.0])
    velocity1 = velocity/(1.0+R)
    velocity2 = -velocity*R/(1.0+R)
    n_dm = 5_000_000
    n_star = int(f_star*n_dm)
    num_particles = {"dm": n_dm, "star": n_star}
    print(num_particles)
   # print(r200_1, r200_2)
    ics = cg.ClusterICs(f"run{num}", 2, [f"run{num}_profile_1.h5", f"run{num}_profile_2.h5"], 
                        [center1, center2], [velocity1, velocity2], 
                        num_particles=num_particles, r_max=p_rmax)
    cg.setup_gamer_ics(ics, regenerate_particles=False)
run_dicts = [
        {"num": 1, "r_c1": 0.2, "r_c2": 0.1, "alpha1": 1.0, "alpha2": 2.0, "R" : 5.2, "M200" : 1.2e15, "b" : 150, "v" : 1000, "d":4000} #, #S1 - S2 merger
       # {"num": 2, "r_c1": 0.1, "r_c2": 0.1, "alpha1": 0.5, "alpha2": 2.0, "R": 1, "M200": 8e14, "b" : 100, "v" : 1512, "d": 1500} #N - S2 merger
 #       {"num": 3, "r_c1": 0.5, "r_c2": 0.1, "alpha1": 0.1, "alpha2": 1.0, "R" : 1, "M200" : 4e14, "b" : 150, "v" : 4000, "d": 3000} #S1, NW2; not currently merging
#        {"num": 4, "r_c1": 0.3, "r_c2": 0.2, "alpha1": 0.3, "alpha2": 1.0, "R" : 6.67, "M200" : 2.3e14, "b" : 150, "v" : 4500}]
#        {"num": 5, "r_c1": 0.3, "r_c2": 0.3, "alpha1": 0.3, "alpha2": 0.3, "R": 4, "M200" : 1e15, "b" : 2000, "v" : 4000}],
#        {"num": 6, "r_c1": 0.1, "r_c2": 0.1, "alpha1": 0.1, "alpha2": 1, "R" : 4, "M200" : 1e15, "b" : 1600, "v" : 4000}
]
for run_dict in run_dicts:
    make_clusters(run_dict)
