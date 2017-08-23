
% set the input directory.
D = '../../examples/spikes-lg/';

% load the instance data.
load([D, 'instance.dat.gz']);

% load the models.
vbcs = load([D, 'model-vbcs.dat.gz']);
vrvm = load([D, 'model-vrvm.dat.gz']);

% collect the necessary data.
dat = [x0, vrvm.x, vbcs.x];

% write the data to a text file.
save('-ascii', 'spikes.dat', 'dat');

