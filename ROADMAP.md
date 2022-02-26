# ROADMAP

## Milestone 0: rewrite of the original implementation by Zhang, Rozenberg and Nájera

- [x] Port to MATLAB the [modern python implementation](legacy/PYTHON/real_ipt-text_v3.ipynb) by Óscar Nájera and crosscheck the basic output it provides.

- [x] Check the effect of the self-consistency condition by inspecting the spectral functions at one loop vs at convergence, with both the initial guesses for the bath. [it requires to implement a proper self-consistency control, possibly with linear mixing]

- [x] Analyze the U-driven MIT by extracting the quasiparticle weight from the self-energy, determining Uc<sub>2</sub> and Uc<sub>1</sub> points at some relevant temperatures and capturing first-order behaviour in between.

- [x] Complete the original tutorial by inspecting the T-driven MIT and defining the supercritical behaviour (bad metal and pseudogap phases).

#### <p align="right"> [`> v0.1`](https://github.com/bellomia/MOTTlab/releases/tag/v0.1) </p>

## Milestone I: functional phase diagrams, pseudo-order-parameter zoo, performance

- [x] Compute a mottness-marker based on the divergence of the scattering rate (Im[Sigma(0)]... very sensible, basically inexpensive), so to obtain sharper phase diagrams with respect to the Z-derived ones. It should also give insight about the supercritical phases: bad metal and pseudogap. [**this yet to check**]   
  > At the moment we implement the "strenght of correlations", S = norm[Sigma(0)-Sigma(∞)], as defined in `PRL 114 185701 (2015)`, with actual neat results: the marker is almost zero accross the whole FL phase and starts increasing very fast in the Mott insulator.

- [x] Compute the Luttinger integral, as defined in `PRB 90 075150 (2014)`. Since it appears to be quantized _at very low temperatures_, it could become the definitive **flag** for _quantum_ phase diagrams; much better than Z or S for it is an **integer**. [→ easier phase-boundary recognition!]  
    > Luttinger Theorem currently works for very low temperatures only. It might well be an inherent limitation, restricting its domain to the quantum Mott transition. Also note that to have a sharp first order step-up at the transition a very highly frequency resolution is needed, so to make the IPT solver performance-critical! [solved brilliantly with fast convolutions, see the `solver-optimization` section below]

- [x] `SOLVER-OPTIMIZATION`: make the SOPT runs faster, by optimizing the needed convolutions. [implemented a pow2-optimized FFTW-based custom algorithm that actually greatly improves the cpu-time for the `wres=2^15` calculation: almost a x10 overall speedup!]

- [ ] `LOOP-OPTIMIZATION`: insert a "restarting" protocol for lines and full phase diagram spans. The gloc0=0 condition appers to be too unstable to obtain accurate UC1 values. Furthermore this would most probably speed up a lot the convergence, by lowering the required number of iterations.

- [x] `HPC-OPTIMIZATION`: configure an interface to cluster facilities and define the scheduling resources for optimal running [no distributed computing, just built-in handling of shared-memory parallelization]  

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

- [ ] Start writing markdown [docs](./docs), with the aim to document the implementation details (usage instruction should go directly on the [README](./README.md) instead...)

