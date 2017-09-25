
% set the input directory.
D = '../../examples/gauss/';

% load the instance data.
load([D, 'instance.dat.gz']);

% load the models.
vbcs = load([D, 'model-vbcs.dat.gz']);
vrvm = load([D, 'model-vrvm.dat.gz']);

% collect the necessary data.
datx = [[1 : n]', x0, vrvm.x, vbcs.x];
daty = [[1 : m]', y, A * vrvm.x, A * vbcs.x];

% write the data to a text file.
save('-ascii', 'gauss-x.dat', 'datx');
save('-ascii', 'gauss-y.dat', 'daty');

