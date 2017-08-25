
% cwrite: write the contours of a matrix to a text file.
%
function cwrite (x, lev, fname)
  % compute the contours.
  C = contourc(x, lev);
  N = columns(C);
  ii = 1;

  % open the output file.
  fh = fopen(fname, 'w');

  % loop over each curve.
  while (ii <= N)
    % break each curve with a blank line.
    if (ii > 1)
      fprintf(fh, '\n');
    end

    % get the current curve size and height.
    z = C(1, ii);
    n = C(2, ii);

    % get the current curve points.
    x = C(1, ii+1 : ii+n);
    y = C(2, ii+1 : ii+n);

    % write the current curve.
    for jj = 1 : n
      fprintf(fh, '%e %e %e\n', x(jj), y(jj), z);
    end

    % increment the counter.
    ii += n + 1;
  end

  % close the output file.
  fclose(fh);
end

