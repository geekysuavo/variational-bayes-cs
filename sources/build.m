
% get the command line arguments.
args = argv();

% set the model type to its default value.
mtype = 'vbcs';

% check if arguments were specified.
if (length(args) && ...
    (strcmp(args{end}, 'vbcs') || ...
     strcmp(args{end}, 'vrvm') || ...
     strcmp(args{end}, 'vrvm_lowrank') || ...
     strcmp(args{end}, 'ifsbl') || ...
     strcmp(args{end}, 'nesta')))
  % set the new model type.
  mtype = args{end};
end

% get the model function handle.
mfunc = str2func(mtype);

% read in the source files.
addpath('../../sources');

% load the instance data.
load('instance.dat.gz');

% execute a timed reconstruction.
tic;
[mu, obj, parms] = mfunc(y, A, At, nu0, lambda0, alpha0, beta0, iters);
runtime = toc;

% write the results.
save('-binary', '-z', ['model-', mtype, '.dat.gz']);

