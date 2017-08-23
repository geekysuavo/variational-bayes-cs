
% set the input directory.
D = '../../examples/nmr2d/';
addpath(D);

% load the instance data.
load([D, 'instance.dat.gz']);

% load the models.
%vbcs = load([D, 'model-vbcs.dat.gz']);
%nesta = load([D, 'model-nesta.dat.gz']);

% build index vectors.
ii = vec(repmat([1 : n(1)]', n(2), 1));
jj = vec(repmat([1 : n(2)], n(1), 1));

% compute the ground truth spectrum.
x0 = sumsq(mxifft2(y0, n), 3);
x0 = vec(fftshift(x0));

% collect the necessary data.
dat = [ii, jj, x0];

% write the data to a text file.
save('-ascii', 'nmr2d.dat', 'dat');

