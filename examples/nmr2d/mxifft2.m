
% two-dimensional multicomplex inverse fft.
function y = mxifft2 (x, n)
  % compute the first-dimension inverse fourier transform.
  x12 = ifft(complex(x(:,:,1), x(:,:,2))) .* sqrt(n(1));
  x34 = ifft(complex(x(:,:,3), x(:,:,4))) .* sqrt(n(1));

  % repack.
  x(:,:,1) = real(x12);
  x(:,:,2) = imag(x12);
  x(:,:,3) = real(x34);
  x(:,:,4) = imag(x34);

  % compute the second dimension inverse fourier transform.
  x13 = ifft(complex(x(:,:,1), x(:,:,3)), [], 2) .* sqrt(n(2));
  x24 = ifft(complex(x(:,:,2), x(:,:,4)), [], 2) .* sqrt(n(2));

  % repack.
  x(:,:,1) = real(x13);
  x(:,:,3) = imag(x13);
  x(:,:,2) = real(x24);
  x(:,:,4) = imag(x24);

  % return.
  y = x;
end

