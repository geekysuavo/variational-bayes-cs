
% set the input directory.
D = '../../examples/nmr2d/';
addpath(D);

% load the instance data.
load([D, 'instance.dat.gz']);

% load the models.
vbcs = load([D, 'model-vbcs.dat.gz']);
nesta = load([D, 'model-nesta.dat.gz']);

% build index vectors.
ii = vec(repmat([1 : n(1)]', n(2), 1));
jj = vec(repmat([1 : n(2)], n(1), 1));

% compute the ground truth spectrum.
x0 = mxifft2(y0, n)(1);
x0 = vec(fftshift(x0));

% compute the estimated spectra.
vbcs.xr = vec(fftshift(vbcs.x(1)));
nesta.xr = vec(fftshift(nesta.x(1)));

% collect the necessary data.
dat = [ii, jj, x0, vbcs.xr, nesta.xr];

% write the data to a text file.
save('-ascii', 'nmr2d.dat', 'dat');

