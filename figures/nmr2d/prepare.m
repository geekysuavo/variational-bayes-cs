
% set the input directory.
D = '../../examples/nmr2d/';
addpath(D);

% load the instance data.
load([D, 'instance.dat.gz']);

% load the models.
vbcs = load([D, 'model-vbcs.dat.gz']);
nesta = load([D, 'model-nesta.dat.gz']);

% set the contours.
clev0 = logspace(log10(0.01), log10(0.5), 20);
vbcs.clev = logspace(log10(0.02), log10(4), 20);
nesta.clev = logspace(log10(0.01), log10(2), 20);

% compute the ground truth spectrum.
x0 = fftshift(mxifft2(y0, n)(:,:,1));

% compute the estimated spectra.
vbcs.xr = fftshift(vbcs.x(:,:,1));
nesta.xr = fftshift(nesta.x(:,:,1));

% write contours from each spectrum.
cwrite(x0,       clev0,      'nmr2d-orig.dat');
cwrite(vbcs.xr,  vbcs.clev,  'nmr2d-vbcs.dat');
cwrite(nesta.xr, nesta.clev, 'nmr2d-nesta.dat');

