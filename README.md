# A2744
Abell 2744 is perhaps the most iconic example of a complex galaxy cluster merger. It seems that every time someone observers it anew, be it in the X-ray or gravitational lensing, they discover a new merging component, earning the system the nickname Pandora's Cluster. Merging clusters, due to their highly dynamic nature, are fantastic sites in which to explore questions of plasma physics. The appearance of shock fronts, cold fronts and radio relics all hold information about the strength and distribution of magnetic fields, thermal conduction, viscosity, and electron-ion equilibriation processes in the weakly magnetised intracluster medium (ICM). To measure these quantities, however, we first need a dynamical model of the system. What is moving where? 

Reconstructing the dynamics is complex enough for [binary mergers](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.1201C/abstract), with key free parameters including the cluster mass profiles, gas density and temperature profiles, relative velocities, angular momentum, and viewing direction (i.e. 10 parameters). Every additional merging component adds O(N) parameters, and exploring the parameter space blindly quickly becomes unfeasible. Here, we re-assess the available observations of the X-ray emitting ICM, galaxy velocities, and strong and weak gravitational lensing to simplify the problem.

## Observational constraints
First, we note that the X-ray emission suggests three main structures:
- A primary ICM, with strong lensing peaks around two southern BCGs S1 and S2.
- The S1-S2 system neighbours a sharp cold front and a bow shock ahead of it. This is a tell-tale signature of a recent merger between the S1 and S2 substructures; specifically, it suggests that S2 is less massive and has a strong cool core. 
- A northern extension that harbors a BCG (BCG-N), a cold front, and a shock front ahead of it. This means that whatever structure is centered on BCG-N recently had a pericenter passage with the cluster hosting either S2 or S1.
- A northwest interloper NW, which is close to - but offset from - two BCGs, NW-1 and NW-2.

Next, the galaxy velocities (compiled from Owers et al. 2011 and Furtak et al. 2022) suggest the following:
- S1 and S2 have a $\Delta z$ = 0.02, corresponding to a line-of-sight velocity offset of > 5000 km/s or, if the redshift is entirely cosmological, a separation of ~57 Mpc (with S2 behind S1). Given the X-ray cold front, however, we know that the separation cannot be that large. 
- BCG-S2 and BCG-N have no relative line-of-sight velocity, i.e. they can be treated as a plane-of-sky system. 
- NW-1 and NW-2 are also very close to the plane-of-sky.

The strong (Furtak et al 2022.) and weak (Medezinski et al 2016.) gravitational lensing maps show that:
- There are strong lensing peaks around all the aforementioned BCGs.
- The most massive component is centered on S1-S2.
- The NW-1 and NW-2 peaks are both offset to the east of the NW interloper, i.e. they are both closer to the primary ICM than NW. This has been particularly challenging to explain in the literature, although all the scenarios presented thus far have either considered only NW-2 or considered both BCGs to be associated with NW.

Last but not least, radio observations (Rajpurohit et al. 2022 and Paul et al. 2019) find three strong radio relics, i.e. extended features not associated with a point source. 
- R1 lies to ~ 1 Mpc the northeast of the cluster. Eckert et al. 2016 found an X-ray shock in this direction with XMM-Newton; as is common in the literature, the X-ray and radio Mach numbers differ. 
- R2 lies to the southeast of the cluster, ahead of the S1-S2 bow shock. It is not detected in the X-ray.
- R3 overlaps with NW. It is not concave with respect to the cluster center, as is characteristic of merger shocks, but convex. However, Böss et al. 2023 have reproduced exactly such a morphology when a shock front traveling outwards interacts with a heavy substructure falling in.
We posit that R2 and R3 are actually a shock pair, associated with an earlier merger of NW-1 with either S1 or S2. Neither S1 or S2 is preferred at this point, since they have an equal line-of-sight velocity offset with NW-1 and all three BCGs are aligned with the R2-R3 axis. The reverse shock R3 was then bent by the infall of the NW interloper, which in turn hosts NW-2.
- This means that, given the current trajectory of BCG-N, it would have previously interacted with the NW-1 - S1/S2 system. This interaction could have driven the R1 shock.

Thus, for the first time, we have a qualitative merger picture that explains all the BCG positions and velocities, gravitational lensing peaks, X-ray and radio features! 
- NW-1 is actually the BCG associated with the primary ICM.
- It underwent a merger with S2 some time ago, generating the R2 and R3 radio relics.
- Then, BCG-N passed both of them, forming the northern cold and shock front pair.
- Next, S1 and S2 underwent a merger close to, but not entirely along, the line of sight, producing the southern cold front and bow shock.
- Lastly, NW is undergoing a widely offset (high impact parameter/angular momentum) merger with the remaining components, leaving a trail of cool gas in its wake. At this phase of the merger, its BCG (NW-2) is leading the cold front. 

## Testing the picture with simulations
We start by simulating the following binary mergers:
- S1 and S2 merge close to the line of sight, producing a cold front and bow shock with the observed separations and a line-of-sight velocity offset of the BCGs ~ 5000 km/s. We pick the time of merger as the one that best reproduces the shock and cold front positions and surface brightness jumps as well as BCG velocity offset.
- NW-1 and S2 are oriented such that they produce a ~ 3000 km/s line-of-sight velocity offset at the present time, and produce shocks at the present locations of R2 and R3. We determine the time of merger based on the current distance and los velocity offset between NW-1 and S2, and the shock separation of R2 and R3. 
- The S2 and BCG-N merger is in the plane of the sky, and is timed to reproduce the BCG separation and approximate position of the northern cool core and shock front.
These binary simulations will allow us to constrain the 3D positions of all the BCGs with respect to one another.

Now that we know the relative timings of the above mergers, we set up the following triple mergers:
- NW-1 + (S2 + N) to produce R1, R2 and R3, in addition to the aforementioned cold and shock fronts.
- S2 + (NW-1 + NW-2) to produce the observed offset between NW-2 and the NW X-ray emission. 
Parantheses indicate plane-of-sky mergers. At this point, we expect the 3D orientations to differ slightly from the binary merger predictions due to 3-body effects.

Lastly, we will run the full system: {(NW-1 + [S2 + S1]) + N} + NW-2. 
Thus, we have six mergers to simulate - three binary, two triple, and one quadruple. Each simulation will have to be run a few different times to account for uncertainties in the inherent properties of the merging components. Different viewing directions do not need to be simulated; they are produced simply by projecting images along different lines of sight. 

## How to run the code
1 - Use [initial_conditions.py](initial_conditions.py) to generate the initial profiles of the sub-clusters. The free parameters here are cluster mass & dark matter concentration (these two are assumed to be initially along the mass-concentration relation of Diemer & Kravtsov 2015), gas cool core size (r_c) and strength (inner slope alpha). The gas extent and outer radius are kept fixed for now, although they could affect the strength of the radio relics, which are now all in the cluster outskirts. We do not, however, aim to reproduce their strengths at this point, because full modeling of the radio relics requires detailed modeling of the magnetic field and cosmic ray population. 

2 - For any binary, plane-of-sky merger, [initial_conditions.py](initial_conditions.py) can be used to output the (X,Y) positions of the cluster centers and velocities, used in the GAMER-2 setup file Input__TestProb. This assumes that the merger is in the X-Y plane and the Z-axis is the line of sight. If the merger axis is $\theta$ away from the z-axis (default line of sight), use the following:

$d_{XY} = d\cdot\sin(\theta)$

$d_Z = \sqrt{d^2 - d_{XY}^2}$

$v_{Z} = v\cdot\cos(\theta)$

$v_{XY} = \sqrt{v^2 - v_{Z}^2}$

3 - For triple mergers, we have to further decompose the plane-of-sky quantities into X and Y components, since these are no longer equivalent to different viewing directions. If a merger axis between two components is $\phi$ away from the x-axis:

$d_X = d_{XY}\cdot\cos{\phi}$

$d_Y = d_{XY}\cdot\sin{\phi}$

and similarly for ${v_X, v_Y}$.
