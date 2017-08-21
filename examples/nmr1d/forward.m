
% define the forward operator.
function y = forward (x, schd, n)
  % inverse fourier transform and subsample.
  f = ifft(x) .* sqrt(n);
  y = f(schd,:);
end

