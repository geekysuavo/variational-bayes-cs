
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
unset colorbox
set palette gray

set xrange [-1475 : -415]
set yrange [200 : 710]

set label 1 '$\m{\Phi}^* \m{y}^0$' at graph lx, ly
p 'nmr2d-orig.dat' u (-$1):2:3 w l lc pal

set xrange [595 : 1630]
set yrange [-1855 : -1344]

set label 1 '$\xhat_{\textsc{nesta}}$' at graph lx, ly
p 'nmr2d-nesta.dat' u 1:(-$2):3 w l lc pal

set label 1 '$\xhat_{\textsc{vbcs}}$' at graph lx, ly
p 'nmr2d-vbcs.dat' u 1:(-$2):3 w l lc pal

unset multiplot

