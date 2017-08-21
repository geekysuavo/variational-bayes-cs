
% -*- texinfo -*-
% @deftypefn  {} {@var{x} =} vrvm (@var{y}, @var{A})
% @deftypefnx {} {@var{x} =} vrvm (@var{y}, @var{A}, @var{At})
% @deftypefnx {} {@var{x} =} vrvm (@var{y}, @var{A}, @var{At}, @var{mu0})
% @deftypefnx {} {@var{x} =} vrvm (@var{y}, @var{A}, @var{At}, @var{mu0}, @
% @var{lambda0})
% @deftypefnx {} {@var{x} =} vrvm (@var{y}, @var{A}, @var{At}, @var{mu0}, @
% @var{lambda0}, @var{alpha0})
% @deftypefnx {} {@var{x} =} vrvm (@var{y}, @var{A}, @var{At}, @var{mu0}, @
% @var{lambda0}, @var{alpha0}, @var{beta0})
% @deftypefnx {} {@var{x} =} vrvm (@var{y}, @var{A}, @var{At}, @var{mu0}, @
% @var{lambda0}, @var{alpha0}, @var{beta0}, @var{iters})
% @deftypefnx {} {[@var{x}, @var{elbo}] =} vrvm (@dots{})
% @deftypefnx {} {[@var{x}, @var{elbo}, @var{eta}] =} vrvm (@dots{})
% Recover a sparse signal @var{x} from incomplete data @var{y} and a
% measurement matrix (@var{A} or @var{A},@var{At}) using the Variational
% Relevance Vector Machine method.
%
% See @ref{vbcs} for detailed usage information.
% @end deftypefn
%
function [x, elbo, eta] = ...
vrvm (y, A, At, mu0, lambda0, alpha0, beta0, iters)
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

  % check for a prior mu parameter.
  if (nargin < 4 || isempty(mu0))
    % none specified. use a default value.
    mu0 = 1e-3;
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

  % initialize the transformed data vector.
  h = At(y);

  % get the problem sizes.
  m = length(y);
  n = length(h);

  % if the projector was not defined, define it.
  if (!exist('G', 'var'))
    G = AtA(eye(n));
  end

  % initialize the weight means.
  x = zeros(n, 1);

  % initialize the precisions.
  xi = repmat(alpha0 ./ beta0, n, 1);

  % initialize the noise.
  tau = mu0 / lambda0;

  % compute the constant updated parameters.
  alpha = alpha0 + (1/2);
  mu = mu0 + m/2;

  % compute the constant portion of the lower bound.
  L0 = n * alpha0 * log(beta0) + mu0 * log(lambda0) ...
     + n * (lgamma(alpha) - lgamma(alpha0)) ...
     + (lgamma(mu) - lgamma(mu0)) ...
     + n * alpha + mu - (m/2) * log(2*pi) + n/2;

  % initialize the lower bound vector.
  elbo = repmat(L0, iters, 1);

  % iterate.
  for it = 1 : iters
    % update the variances.
    Sigma = inv(diag(xi) + tau .* G);

    % update the means.
    x = tau .* Sigma * h;

    % compute the trace term.
    trGxx = sum(vec(G .* (x * x' + Sigma)));
    x2 = conj(x) .* x + diag(Sigma);

    % update the precisions.
    beta = beta0 + 0.5 .* x2;
    xi = alpha ./ beta;

    % update the noise.
    lambda = lambda0 + 0.5 * y' * y - h' * x + 0.5 * trGxx;
    lambda = abs(lambda);
    tau = mu / lambda;

    % compute the new value of the lower bound.
    L = -(tau/2) * abs(y' * y - 2 * h' * x + trGxx) ...
      - sum(((alpha/2) * x2 + alpha * beta0) ./ beta) ...
      - mu * log(lambda) - mu * lambda0 / lambda - alpha * sum(log(beta)) ...
      + sum(log(diag(chol(Sigma))));

    % store the bound.
    elbo(it) += L;
  end

  % check if the complete set of variational parameters was requested.
  if (nargout >= 3)
    % store the parameters over x.
    eta.xhat = x;
    eta.v = diag(Sigma);

    % store the parameters over xi.
    eta.alpha = alpha;
    eta.beta = beta;

    % store the parameters over tau.
    eta.mu = mu;
    eta.lambda = lambda;
  end
end

