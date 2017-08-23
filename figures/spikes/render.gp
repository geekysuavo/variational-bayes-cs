
set terminal epslatex color size 2.5in, 5in header \
 '\newcommand{\ft}[0]{\footnotesize}'

set output 'spikes.tex'

set lmargin 0.2
set rmargin 0.2

set multiplot layout 3, 1

unset key
set xlabel '$i$' off 0, 0.5
set xtics 5000 format '{\ft %g}'
set ytics 1 format '{\ft %g}'
set xrange [1 : 10000]
set yrange [-1.2 : 1.2]

set ylabel '$\m{x}^0$' off 3
p 'spikes.dat' u 1:2 w l lc rgb 'black'

set ylabel '$\xhat_{\textsc{vrvm}}$'
p 'spikes.dat' u 1:3 w l lc rgb 'black'

set ylabel '$\xhat_{\textsc{vbcs}}$'
p 'spikes.dat' u 1:4 w l lc rgb 'black'

unset multiplot

