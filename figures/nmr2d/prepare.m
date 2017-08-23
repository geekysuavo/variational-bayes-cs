
% set the input directory.
D = '../../examples/nmr2d/';

% load the instance data.
load([D, 'instance.dat.gz']);

% load the models.
vbcs = load([D, 'model-vbcs.dat.gz']);
vrvm = load([D, 'model-vrvm.dat.gz']);

% collect the necessary data.
ax = linspace(-pi, pi, n)';
dat = [vrvm.x, vbcs.x];
dat = [ax, real(fftshift(dat, 1))];

% write the data to a text file.
save('-ascii', 'nmr2d.dat', 'dat');

