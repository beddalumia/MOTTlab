## What `MOTTlab` is
An efficient implementation of (real axis) IPT-based dynamical mean-field theory for the half-filled Mott-Hubbard transition on the Bethe lattice, within pure MATLAB.

## What you can get out of it
Here we present a few examples of what you can obtain fiddling around the `MOTTlab`.

> Characterize the interaction-driven quantum Mott transition through different physical markers, such as the _timeless classic_ quasiparticle renormalization weight _Z<sub>F</sub>_ , or the _quantized_ Luttinger integral _I<sub>L</sub>_, getting some exotic hints of a topological interpretation.

![Mott-Transition-Markers-zeroT](./showcase/ZvsL_beta1e6_wres2e15.svg) 

> Truly visualize what happens at different relevant temperatures by means of pretty animations!

 Temperature  | DOS | SELF-ENERGY
:-------------------------:|:-------------------------:|:-------------------------:
`T -> 0: 2nd order MIT` |![Mott-Transition-AnimatedDOS-zeroT](./showcase/uDOS_zeroT.gif) | ![Mott-Transition-AnimatedSIGMA-zeroT](./showcase/uSigma_zeroT.gif) 
`Intermediate T: 1st order MIT`  |![Mott-Transition-AnimatedDOS-intermediateT](./showcase/uDOS_beta50.gif) | ![Mott-Transition-AnimatedSIGMA-intermediateT](./showcase/uSigma_beta50.gif)
`High T: supercritical MIT` |![Mott-Transition-AnimatedDOS-highT](./showcase/uDOS_beta1.gif) | ![Mott-Transition-AnimatedSIGMA-highT](./showcase/uSigma_beta1.gif)  

## Dependencies
I've tried hard to avoid exceeding with the infamous MathWorks add-ons. At the moment the only additional dependence is the [Signal Processing Toolbox](https://it.mathworks.com/products/signal.html), which provides an efficient implementation of the discrete Hilbert transforms, at the earth of the numerical evaluation of the Kramers-Kronig relationships. Nevertheless I provide also a FOSS options, through [GNU Octave](https://www.gnu.org/software/octave/) packages.  

First of all you need to install the main software. The easiest way that I know is through Homebrew (you won't need sudo rights, they just take care of everything):

```
 brew install octave
```

Then run an instance of Octave in your shell (just type `octave`) and install the [control package](https://octave.sourceforge.io/control/), by typing:

```
pkg install -forge control
```
Finally proceed to install the [signal package](https://octave.sourceforge.io/signal/index.html), which provides a FOSS implementation<sup>†</sup> of the `hilbert()` function.

<sub>† note that while the MATLAB `fft` family is always built on top of [FFTW](https://www.fftw.org), with Octave you could end up wrapping [FFTPACK](https://www.netlib.org/fftpack/) instead, depending on how your binaries have been built. </sub>

So, running

```
pkg install -forge signal
```

you should now be set... enjoy! 

## Licensing and legacy code
This code has been implemented taking inspiration from these two didactic sources:
1. http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf ([local copy](docs/haule_IPTtheory_rutgers.pdf))
2. https://www.cond-mat.de/events/correl19/manuscripts/rozenberg.pdf ([local copy](docs/rozenberg_review_julich.pdf))
and the hands-on material given therein, in particular a tutorial-intended jupyter notebook provided by Óscar Nájera (available [on his cloud](http://mycore.core-cloud.net/index.php/s/oAz0lIWuBM90Gqt), or [locally](legacy/PYTHON/real_ipt-text_v3.ipynb)) under the BSD 3-Clause License. Here we provide an extensive and efficient MATLAB rewrite, augmented by basic phase-diagram workflows, convergence and self-mixing control, various post-processing tools, finally some exploration of novel material with explicit references to the relevant research papers. Everything yet under the BSD 3-Clause License. You can read more about allowed use of this code in the [LICENSE](./LICENSE) file.

<kbd> © 2020 | Gabriele Bellomia 


