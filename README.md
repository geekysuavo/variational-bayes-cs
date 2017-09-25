
# Variational Bayesian Compressed Sensing

An extension of the variational relevance vector machine (VRVM) for
highly scalable sparse Bayesian learning in compressed sensing
recovery problems.

Technical details of the algorithm are _in preparation for submission to_:

> Worley, B., Malliavin, T. E., Nilges, M., _Scalable Sparse Bayesian
> Learning for Compressed Sensing Recovery_, IEEE Transactions on
> Signal Processing, 2017.

## Introduction

Compressed sensing deals with the recovery of signals from incomplete
measurements, using knowledge of the signal structure in some other
domain (e.g. frequency, wavelet, _etc._). While recovery problems
in compressed sensing have been analyzed from the perspective of
Bayesian inference, most of the resulting algorithms do not yield
fully Bayesian estimates, and fully Bayesian approaches scale
poorly.

The variational Bayesian compressed sensing (VBCS) algorithm extends
the VRVM to bring the computational requirements of fully Bayesian 
estimation (or rather its variational approximation) down to levels
that enable modeling of large problem instances.

### Source code

The main source code of **vbcs** is stored in [sources](sources), where
implementations of both VBCS and the VRVM are provided.

### Example problems

Numerical problems demonstrating the features of **vbcs** are stored in
subdirectories of [examples](examples).

### Figures

Scripts for generating the figures in the **vbcs** publication are stored
in subdirectories of [figures](figures).

## Licensing

The **vbcs** sources and example files released under the
[MIT license](https://opensource.org/licenses/MIT). See the
[LICENSE.md](LICENSE.md) file for the complete license terms.

*~ Brad.*

