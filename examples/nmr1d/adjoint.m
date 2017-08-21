
% define the adjoint operator.
function x = adjoint (y, schd, n)
  % create an infilled vector of y-values.
  f = zeros(n, columns(y));
  f(schd,:) = y;

  % compute the fourier transform.
  x = fft(f) ./ sqrt(n);
end

