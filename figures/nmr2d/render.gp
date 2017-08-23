
set terminal epslatex color size 6in, 2in header \
 '\newcommand{\ft}[0]{\footnotesize}'

lx = 0.05
ly = 0.9

set output 'nmr2d.tex'

#set lmargin 0.2
set rmargin 0.2

set multiplot layout 1, 3

unset key
set xlabel '$\omega_1$' off 0, 0.5
set ylabel '$\omega_2$'
unset xtics
unset ytics

set label 1 '$\m{\Phi}^* \m{y}^0$' at graph lx, ly
p 'nmr2d.dat' u 1:2 w l lc rgb 'black'

set label 1 '$\xhat_{\textsc{nesta}}$' at graph lx, ly
p 'nmr2d.dat' u 1:3 w l lc rgb 'black'

set label 1 '$\xhat_{\textsc{vbcs}}$' at graph lx, ly
p 'nmr2d.dat' u 1:4 w l lc rgb 'black'

unset multiplot

