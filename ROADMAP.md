# ROADMAP

## Milestone 0: rewrite of the original implementation by Zhang, Rozenberg and Nájera

- [x] Port to MATLAB the [modern python implementation](legacy/PYTHON/real_ipt-text_v3.ipynb) by Óscar Nájera and crosscheck the basic output it provides.

- [x] Check the effect of the self-consistency condition by inspecting the spectral functions at one loop vs at convergence, with both the initial guesses for the bath. [it requires to implement a proper self-consistency control, possibly with linear mixing]

- [x] Analyze the U-driven MIT by extracting the quasiparticle weight from the self-energy, determining Uc<sub>2</sub> and Uc<sub>1</sub> points at some relevant temperatures and capturing first-order behaviour in between.

- [x] Complete the original tutorial by inspecting the T-driven MIT and defining the supercritical behaviour (bad metal and pseudogap phases).

#### <p align="right"> `> v0.1` </p>

## Milestone I: functional phase diagrams, pseudo-order-parameter zoo

- [x] Compute a MIvsFL marker based on the divergence of Im[Sigma(0)] (much sensible, very inexpensive), so to obtain sharper phase diagrams with respect to the Z-derived ones.  
 At the moment we implement the "strenght of correlations", S = norm[Sigma(0)-Sigma_HF], as defined in `PRL 114 185701 (2015)`, with actual neat results: the marker is almost zero accross the whole FL phase and starts increasing fast in the Mott insulator.

- [x] Compute the Luttinger Integral, as defined in `J. Phys.: Condens. Matter 28 (2016) 025601` and `PRB 102 081110(R) (2020)`. Since it appears to be quantized it could become the definitive _flag_ for phase diagrams; much better than Z or S for it is an _integer_.  
UPDATE: Luttinger Theorem currently works for very low temperatures only. It may well be an inherent limitation (in that case it would be a good marker for low temperature lines, but not phase diagrams). Also note that the quality of the first order step at the transition highly depends on the frequency resolution, so to have sharp boundaries you need to have at least `wres=2^14`, which makes the IPT solver far slower. [solved brilliantly with fast convolutions, see the `optimization` section below]

- [ ] Insert a "restarting" protocol for lines and full phase diagram spans. The gloc0=0 condition appers to be too unstable to obtain accurate UC1 values. Furthermore this would most probably speed up a lot the convergence (lowering the required number of iterations), so to enable a relevant slowdown of the single iteration (e.g. bigger frequency resolution).

- [ ] Implement an ergonomic 'full-roundtrip' protocol, so to enable suitable explorations of the coexistence region.

<!---
#### <p align="right"> `> v1.0` </p>
--->

## Milestone II: topological madness

- [ ] Try to reproduce the main result of `PRB 102 081110(R)`, namely the Lanczos tridiagonalization of the self-energy leading to the mapping of the quantum MIT to a generalized SSH SP-TQPT. The original result is achieved within DMFT/NRG and an insane bath dimension, so if we succeed this could be even a [ReScience](http://rescience.github.io) submission (given everything is Octave compatible).

<!---
#### <p align="right"> `> v2.0`</p>
--->

## Milestone III: imaginary-axis formalism and new entries to the zoo

- [ ] Write a Matsubara solver and overload the existing function to work with imaginary frequencies.

- [ ] Add relevant quantities that are accessible only via Mastubara summations [e.g. double occupancy, ref. to `Phys. Rev. B 93, 155162 (2016)`]

- [ ] Compute the Local Entanglement Entropy, as defined in `Rev. Mod. Phys. 80, 517 (2008)` (section V.F) and used in `Mod. Phys. Lett. B 2013 27:05` to characterize the MIT on the Bethe lattice.

<!---
#### <p align="right"> `> v3.0`</p>
--->

## Milestone IV: explorative fun... more lattices!

- [ ] Add a bunch of different particle-hole symmetric lattices, such as finite-coordination Bethe, honeycomb, 2d-square, 1d-chain, 3d-cubic, 3d-bcc... the main inspiration comes from the mighty [GFtools](https://github.com/DerWeh/gftools) by [DerWeh](https://github.com/DerWeh).

- [ ] Rewrite the code for the particle-hole broken case so to enable even more lattices, like kagome, fcc, triangular.

<!---
#### <p align="right"> `> v4.0`</p>
--->

# TODO

## Optimization

- [x] Make the SOPT runs faster, by optimizing the needed convolutions. [implemented a pow2-optimized FFTW-based custom algorithm that actually greatly improves the cpu-time for the `wres=2^15` calculation!]

- [ ] Start writing markdown [docs](./docs), with the aim to document the implementation details (usage instruction should go directly on the [README](./README.md) instead...)

