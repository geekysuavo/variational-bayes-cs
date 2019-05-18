
% -*- texinfo -*-
% @deftypefn  {} {@var{mu} =} vrvm_krylov (@var{y}, @var{A})
% @deftypefnx {} {@var{mu} =} vrvm_krylov (@var{y}, @var{A}, @var{At})
% @deftypefnx {} {@var{mu} =} vrvm_krylov (@var{y}, @var{A}, @var{At}, @
% @var{nu0})
% @deftypefnx {} {@var{mu} =} vrvm_krylov (@var{y}, @var{A}, @var{At}, @
% @var{nu0}, @var{lambda0})
% @deftypefnx {} {@var{mu} =} vrvm_krylov (@var{y}, @var{A}, @var{At}, @
% @var{nu0}, @var{lambda0}, @var{alpha0})
% @deftypefnx {} {@var{mu} =} vrvm_krylov (@var{y}, @var{A}, @var{At}, @
% @var{nu0}, @var{lambda0}, @var{alpha0}, @var{beta0})
% @deftypefnx {} {@var{mu} =} vrvm_krylov (@var{y}, @var{A}, @var{At}, @
% @var{nu0}, @var{lambda0}, @var{alpha0}, @var{beta0}, @var{iters})
% @deftypefnx {} {[@var{mu}, @var{obj}] =} vrvm_krylov (@dots{})
% @deftypefnx {} {[@var{mu}, @var{obj}, @var{parms}] =} vrvm_krylov (@dots{})
% Recover a sparse signal @var{x} from incomplete data @var{y} and a
% measurement matrix (@var{A} or @var{A},@var{At}) using an
% implementation of the VRVM that solves the marginalization
% subproblem of each iteration using lanczos iteration.
%
% See @ref{vrvm} for detailed usage information.
% @end deftypefn
%
function [mu, obj, parms] = ...
vrvm_krylov (y, A, At, nu0, lambda0, alpha0, beta0, iters)
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

  % initialize the transformed data vector.
  h = At(y);

  % get the problem sizes.
  m = length(y);
  n = length(h);

  % initialize the precisions.
  xi = repmat(alpha0 / beta0, n, 1);

  % initialize the noise.
  tau = nu0 / lambda0;

  % compute the constant updated parameters.
  alpha = alpha0 + (1/2);
  nu = nu0 + m/2;

  % compute the constant portion of the objective.
  phi0 = vrvm_const(m, n, nu0, lambda0, alpha0, beta0);

  % initialize the objective vector.
  obj = repmat(phi0, iters, 1);

  % initialize the marginals.
  mu = zeros(n, 1);
  sigma = zeros(n, 1);
  yhat = zeros(m, 1);

  % FIXME
  Tau = zeros(iters, 1);
  Xi = zeros(n, iters);
  Mu = zeros(n, iters);
  Sigma = zeros(n, iters);

  % iterate.
  for it = 1 : iters
    % solve the marginalization subproblem.
    [mu, sigma, f] = util_krylov(A, At, tau, xi, y, yhat);

    % FIXME
    Tau(it) = tau;
    Xi(:,it) = xi;
    Mu(:,it) = mu;
    Sigma(:,it) = sigma;

    % compute the log-determinant.
    lndetS = log(f) + m * log(tau) + sum(log(xi));

    % compute the trace term.
    trGS = (n - xi' * sigma) ./ tau;

    % update the precisions.
    beta = beta0 + 0.5 .* (abs(mu).^2 + sigma);
    xi = alpha ./ beta;

    % update the noise.
    yhat = A(mu);
    lambda = lambda0 + 0.5 * norm(y - yhat)^2 + 0.5 * trGS;
    tau = nu / lambda;

    % compute the objective function value.
    phi = nu * log(lambda) + alpha * sum(log(beta)) - 0.5 * lndetS;

    % store the objective.
    obj(it) += phi;
  end

  % check if the complete set of variational parameters was requested.
  if (nargout >= 3)
    % store the parameters over x.
    parms.mu = mu;
    parms.sigma = sigma;

    % store the parameters over tau.
    parms.nu = nu;
    parms.lambda = lambda;

    % store the parameters over xi.
    parms.alpha = alpha;
    parms.beta = beta;

    % FIXME
    parms.Tau = Tau;
    parms.Xi = Xi;
    parms.Mu = Mu;
    parms.Sigma = Sigma;
  end
end

