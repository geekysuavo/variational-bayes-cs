
## variational-bayes-cs/figures

This directory contains a set of scripts for summarizing some of the
numerical example problems in the [examples](../examples) directory.
The correspondence is:

 * [gauss](gauss) => [examples/gauss](../examples/gauss)
 * [nmr1d](nmr1d) => [examples/nmr1d](../examples/nmr1d)
 * [nmr2d](nmr2d) => [examples/nmr2d](../examples/nmr2d)
 * [spikes](spikes) => [examples/spikes-lg](../examples/spikes-lg)

### Generating the figures

Within each directory, there are two scripts: **prepare.m** and **render.gp**.
The prepare script is used to collect all required data from the corresponding
example directory into a gnuplot-readable text file. Using octave, this is
called as:

```bash
octave prepare.m
```

The render script reads the prepared text file and constructs two files:
an EPS file and a TEX file. The files are built by calling:

```bash
gnuplot render.gp
```

