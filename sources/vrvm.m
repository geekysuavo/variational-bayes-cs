
% -*- texinfo -*-
% @deftypefn  {} {@var{mu} =} vrvm (@var{y}, @var{A})
% @deftypefnx {} {@var{mu} =} vrvm (@var{y}, @var{A}, @var{At})
% @deftypefnx {} {@var{mu} =} vrvm (@var{y}, @var{A}, @var{At}, @var{nu0})
% @deftypefnx {} {@var{mu} =} vrvm (@var{y}, @var{A}, @var{At}, @var{nu0}, @
% @var{lambda0})
% @deftypefnx {} {@var{mu} =} vrvm (@var{y}, @var{A}, @var{At}, @var{nu0}, @
% @var{lambda0}, @var{alpha0})
% @deftypefnx {} {@var{mu} =} vrvm (@var{y}, @var{A}, @var{At}, @var{nu0}, @
% @var{lambda0}, @var{alpha0}, @var{beta0})
% @deftypefnx {} {@var{mu} =} vrvm (@var{y}, @var{A}, @var{At}, @var{nu0}, @
% @var{lambda0}, @var{alpha0}, @var{beta0}, @var{iters})
% @deftypefnx {} {[@var{mu}, @var{obj}] =} vrvm (@dots{})
% @deftypefnx {} {[@var{mu}, @var{obj}, @var{parms}] =} vrvm (@dots{})
% Recover a sparse signal @var{x} from incomplete data @var{y} and a
% measurement matrix (@var{A} or @var{A},@var{At}) using a direct
% implementation of the Variational Relevance Vector Machine.
%
% For matrix-type measurements, only @var{A} is required. For custom
% measurements, the function handles @var{A} and @var{At} are both
% required.
%
% When function handles are supplied, they must be capable of handling
% both column vectors and matrices, where the operations are performed
% independently to each column of the matrix.
%
% The optional hyperparameters @var{nu0} and @var{lambda0}, which both
% default to 1e-3, control the prior distribution over the noise
% precision.
%
% The optional hyperparameters @var{alpha0} and @var{beta0}, which also
% both default to 1e-3, control the prior distribution over the signal
% coefficient precisions.
%
% The optional argument @var{iters} controls the iteration count, and
% defaults to 100 iterations.
%
% The function outputs @var{mu}, the mean estimate of the regression weights,
% and (optionally) the objective function values @var{obj} and the set of
% inferred parameters @var{parms}.
% @end deftypefn
%
function [mu, obj, parms] = ...
vrvm (y, A, At, nu0, lambda0, alpha0, beta0, iters)
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

  % if the projector was not defined, define it.
  if (!exist('G', 'var'))
    G = AtA(eye(n));
  end

  % initialize the weight means.
  mu = zeros(n, 1);

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

  % iterate.
  for it = 1 : iters
    % update the variances.
    U = chol(tau .* G + diag(xi));
    Sigma = chol2inv(U);

    % update the means.
    mu = tau .* Sigma * h;

    % compute the log-determinant.
    lndetS = -2 * sum(log(real(diag(U))));

    % compute the trace term.
    trGS = real(sum(vec(G .* Sigma)));

    % update the precisions.
    mu2 = conj(mu) .* mu + real(diag(Sigma));
    beta = beta0 + 0.5 .* mu2;
    xi = alpha ./ beta;

    % update the noise.
    lambda = lambda0 + 0.5 * norm(y - A(mu))^2 + 0.5 * trGS;
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
    parms.sigma = real(diag(Sigma));

    % store the parameters over tau.
    parms.nu = nu;
    parms.lambda = lambda;

    % store the parameters over xi.
    parms.alpha = alpha;
    parms.beta = beta;
  end
end

