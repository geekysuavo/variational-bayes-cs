
% set the input directory.
D = '../../examples/nmr2d/';
addpath(D);

% load the instance data.
load([D, 'instance.dat.gz']);

% load the models.
vbcs = load([D, 'model-vbcs.dat.gz']);
nesta = load([D, 'model-nesta.dat.gz']);

% set the slice indices.
hset = [1453, 1584, 1657, 1736];
nset = [866, 867, 868];

% set the contours.
clev0 = logspace(log10(0.01), log10(0.5), 20);
vbcs.clev = logspace(log10(0.02), log10(4), 20);
nesta.clev = logspace(log10(0.01), log10(2), 20);

% set the vbcs variance contours.
vbcs.vlev = linspace(0.16, 0.26, 10);

% compute the ground truth spectrum.
x0 = fftshift(mxifft2(y0, n)(:,:,1));

% compute the estimated spectra.
vbcs.xr  = fftshift(vbcs.x(:,:,1));
nesta.xr = fftshift(nesta.x(:,:,1));

% compute the standard deviation.
vbcs.vr = sqrt(fftshift(vbcs.eta.v));

% compute the slices.
vbcs.sh  = [[1 : n(2)]', vbcs.xr(hset,:)', vbcs.vr(hset,:)'];
vbcs.sn  = [[1 : n(1)]', vbcs.xr(:,nset),  vbcs.vr(:,nset) ];
nesta.sh = [[1 : n(2)]', nesta.xr(hset,:)'];
nesta.sn = [[1 : n(1)]', nesta.xr(:,nset) ];

% write contours from each spectrum.
cwrite(x0,       clev0,      'nmr2d-orig.dat');
cwrite(vbcs.xr,  vbcs.clev,  'nmr2d-vbcs.dat');
cwrite(nesta.xr, nesta.clev, 'nmr2d-nesta.dat');

% write the variance contours.
cwrite(vbcs.vr,  vbcs.vlev,  'nmr2d-var.dat');

% write the slices.
dat = vbcs.sh;  save('-ascii', 'nmr2d-vbcs-h1.dat',   'dat');
dat = vbcs.sn;  save('-ascii', 'nmr2d-vbcs-n15.dat',  'dat');
dat = nesta.sh; save('-ascii', 'nmr2d-nesta-h1.dat',  'dat');
dat = nesta.sn; save('-ascii', 'nmr2d-nesta-n15.dat', 'dat');

