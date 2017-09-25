
set terminal epslatex color size 2.5in, 2in header \
 '\newcommand{\ft}[0]{\footnotesize}'

set output 'gauss.tex'

set lmargin 0.2
set rmargin 0.2

unset key
set xtics format ''
set xrange [1 : 200]
set x2range [1 : 1000]
set xlabel '$i,j$' off 0, 0.5
set ylabel '$\m{y}^0$, $\m{y}_{\textsc{vbcs}}$' off 3

p 'gauss-y.dat' u 1:2 w p pt 6 ps 0.5 lc rgb 'black', \
  'gauss-y.dat' u 1:4 w l lc rgb 'black', \
  "<awk '$2>0.1' gauss-x.dat" u 1:(0) \
   w p ax x2y1 pt 1 lw 2.5 lc rgb 'black', \
  "<awk '$4>0.1' gauss-x.dat" u 1:(-0.007) \
   w p ax x2y1 pt 9 ps 0.8 lw 2 lc rgb 'black'

