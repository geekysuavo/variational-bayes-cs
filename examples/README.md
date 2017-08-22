
## variational-bayes-cs/examples

This directory contains a small selection of numerical example problems
that VBCS (and potentially the VRVM) may be used to solve. The list of
problems is:

 * [gauss](gauss): Unit impulse positive spikes, Gauss matrix measurements.
 * [nmr1d](nmr1d): 1D NMR trace, incomplete Fourier operator measurements.
 * [spikes-sm](spikes-sm): Small spikes instance, normal matrix measurements.
 * [spikes-lg](spikes-lg): Large spikes instance, normal matrix measurements.

### Running the scripts

Within each directory, there are two scripts: **instance.m** and **model.m**.
The instance script is used to generate an instance of the given problem.
Using octave, this is called as:

```bash
octave instance.m
```

The model script loads the current problem instance and either builds a
VBCS model or a VRVM model. Using octave, the models are built by calling:

```bash
octave model.m vbcs
octave model.m vrvm
```

The VRVM model script can take a very long time to complete, and typically
requires an order of magnitude more memory to work. Results can be cleaned
up by calling:

```bash
octave clean.m
```

