
% -*- texinfo -*-
% @deftypefn  {} {[@var{mu}, @dots{}] =} gssbl (@var{y}, @var{A}, @dots{})
% Recover a sparse signal @var{x} from incomplete data @var{y} and a
% measurement matrix (@var{A} or @var{A},@var{At}) using a Gibbs
% sampler for the sparse Bayesian learning model.
%
% See @ref{vrvm} for detailed usage information.
% @end deftypefn
%
function [mu, obj, parms] = ...
gssbl (y, A, At, nu0, lambda0, alpha0, beta0, iters, thin)
  % check for the minimum number of arguments.
  if (nargin < 2)
    error('at least two arguments required');
  end

  % check the required measurement vector.
  if (isempty(y) || !isvector(y) || !iscolumn(y))
    error('measurement must be a vector');
  end

  % check the required measurement operator.
  if (!isempty(A) && (nargin < 3 || isempty(At)))
    % matrix specification. check the data type.
    if (!ismatrix(A))
      error('invalid measurement: expected matrix');
    end

    % store copies of the measurement matrix and its gramian.
    B = A;
    G = A' * A;

    % create function handles for forward, inverse, and projection.
    A = @(x) B * x;
    At = @(y) B' * y;
    AtA = @(x) G * x;

  elseif (!isempty(A) && nargin >= 3 && !isempty(At))
    % function handle specification. check the data types.
    if (!is_function_handle(A) || !is_function_handle(At))
      error('invalid measurement: expected function handles');
    end

    % create a function handle for projection.
    AtA = @(x) At(A(x));
  end

  % check for a prior nu parameter.
  if (nargin < 4 || isempty(nu0))
    % none specified. use a default value.
    nu0 = 1e-3;
  end

  % check for a prior lambda parameter.
  if (nargin < 5 || isempty(lambda0))
    % none specified. use a default value.
    lambda0 = 1e-3;
  end

  % check for a prior alpha parameter.
  if (nargin < 6 || isempty(alpha0))
    % none specified. use a default value.
    alpha0 = 1e-3;
  end

  % check for a prior beta parameter.
  if (nargin < 7 || isempty(beta0))
    % none specified. use a default value.
    beta0 = 1e-3;
  end

  % check for an iteration count argument.
  if (nargin < 8 || isempty(iters))
    % none specified. use a default value.
    iters = 100;
  end

  % check for a thinning rate argument.
  if (nargin < 9 || isempty(thin))
    % none specified. use a default value.
    thin = 1;
  end

  % initialize the current sample.
  x = At(y);

  % get the problem sizes.
  m = length(y);
  n = length(x);

  % initialize the sample store.
  X = zeros(n, floor(iters / thin));

  % compute the constant updated parameters.
  alpha = alpha0 + (1/2);
  nu = nu0 + m/2;

  % initialize the objective vector.
  stdnrm = @(N) normrnd(0, 1, [N, 1]);
  obj = [];

  figure(1);
  drawnow();

  % iterate.
  for it = 1 : iters
    % update the tau rate parameter.
    lambda = lambda0 + 0.5 * norm(y - A(x))^2;

    % sample tau.
    tau = gamrnd(nu, 1 / lambda);

    % update the xi rate parameter.
    beta = beta0 + 0.5 .* abs(x).^2;

    % sample xi, with a lower bound.
    xi = bsxfun(@gamrnd, repmat(alpha, n, 1), 1 ./ beta);

    % sample independent normal vectors for sampling x.
    z1 = stdnrm(m);
    z2 = stdnrm(n);

    % build a sample from the distribution N(0, Sigma^-1).
    z3 = sqrt(tau) .* At(z1) + sqrt(xi) .* z2;

    % compute the right-hand vector.
    u = tau .* At(y) + z3;
    v = A(u ./ xi);

    % define the kernel matrix operator, and solve K*d = v.
    K = @(z) z ./ tau + A(At(z) ./ xi);
    d = K(eye(m)) \ v;

    % sample x.
    x = (u - At(d)) ./ xi;

    % store the sample.
    if (mod(it, thin) == 0)
      X(:,it/thin) = x;
    end
  end

  % compute the mean result.
  mu = mean(X, 2);

  % check if the complete set of variational parameters was requested.
  if (nargout >= 3)
    % store the parameters over x.
    parms.mu = mu;
    parms.sigma = var(X, [], 2);

    % store the samples.
    parms.X = X;
  end
end

