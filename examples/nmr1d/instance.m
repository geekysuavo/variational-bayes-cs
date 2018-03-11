
% load the input data.
load('fid.dat.gz');
x0 = x0 ./ max(real(x0));

% set the problem sizes.
m = length(sched);
n = 2 * length(x0);

% set the prior parameters.
nu0 = 50;
lambda0 = 0.005;
alpha0 = 1e-3;
beta0 = 1e-3;

% set the iteration count.
iters = 500;

% define the final operators.
A = @(x) forward(x, sched, n);
At = @(y) adjoint(y, sched, n);

% construct the subsampled measurement.
y = x0(sched);

% save the current instance.
save('-binary', '-z', 'instance.dat.gz');

