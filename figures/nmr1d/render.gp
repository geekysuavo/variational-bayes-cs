
set terminal epslatex color size 2.5in, 5in header \
 '\newcommand{\ft}[0]{\footnotesize}'

set output 'nmr1d.tex'

set lmargin 0.2
set rmargin 0.2

set multiplot layout 3, 1

unset key
set xlabel '$\omega$' off 0, 0.5
set xtics ('{\ft $-\pi$}' -pi, '{\ft %g}' 0, '{\ft $\pi$}' pi)
set ytics ('{\ft %g}' 0)
unset ytics
set xrange [-pi : pi]

set yrange [-0.02 : 0.7]
set ylabel '$\m{\Phi}^* \m{y}^0$'
p 'nmr1d.dat' u 1:2 w l lc rgb 'black'

set yrange [-0.05 : 2.8]
set ylabel '$\xhat_{\textsc{nesta}}$'
p 'nmr1d.dat' u 1:3 w l lc rgb 'black'

set yrange [-0.1 : 4.4]
set ylabel '$\xhat_{\textsc{vbcs}}$'
p 'nmr1d.dat' u 1:5 w l lc rgb 'black'

unset multiplot

