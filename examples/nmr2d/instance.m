
% load the input data.
load('fid.dat.gz');

% set the problem sizes.
m = length(sched);
n = size(y0)(1:2);

% set the linearized sizes.
N = prod(n);

% build an index set for compacting to m-vectors.
S0 = sched(:,1) + (sched(:,2) - 1) * n(1);
S = bsxfun(@plus, S0, N * [0:3]);

% set the prior parameters.
mu0 = m;
sigma0 = 0.001;
lambda0 = m * sigma0^2;
alpha0 = 1e-3;
beta0 = 1e-3;

% set the iteration count.
iters = 500;

% define the final operators.
A = @(x) forward(x, S, n);
At = @(y) adjoint(y, S, n);

% construct the subsampled measurement.
y = y0(S);

% save the current instance.
save('-binary', '-z', 'instance.dat.gz');

