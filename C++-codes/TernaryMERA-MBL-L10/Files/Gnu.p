set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "t"   font ",35" textcolor rgb "black"
set ylabel "S(t)"  font ",35"  textcolor rgb "black"

set ytics nomirror (0.5,0.4,0.3,0.2,0.1,0.0);
set xtics nomirror (35,36,37,38,39,40);
set tics scale 4

set key box t r  
set key font ",30"
set key spacing 1


set output "QRLinear.eps"
set multiplot

p [34:39] [-0.05:0.5] "Entexact.text" u ($1):($2)  t "exact" with line  linetype 2 lw 4 lc rgb "#204a87","Ent1.text" u ($1):($2)  t "{/Symbol s}=10^{-1}" with line  dashtype 4 lw 4  lc rgb "#f57900","Ent5.text" u ($1):($2)  t "{/Symbol s}=10^{-3}" with line  dashtype 6 lw 4  lc rgb "#f5006a"


set origin .17, .4800
set size 0.49,0.350
clear
unset key
unset xlabel
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set ylabel ""  font ",20"  textcolor rgb "black"
set xtics ( 10, 11, 12)
set ytics nomirror (0.1,0.2,0.3,0.4,0.5);
set tics scale 2

#set key box t c
#set key font ",15"
#set key spacing 1.0
#set key width 0.5


p [10:12] [0.04:0.3] "Entexact1.text" u ($1):($2)  t "" with line  linetype 2 lw 4 lc rgb "#204a87","Ent12.text" u ($1):($2)  t "" with line  dashtype 4 lw 4  lc rgb "#f57900","Ent64.text" u ($1):($2)  t "" with line  dashtype 6 lw 4  lc rgb "#f5006a"


unset multiplot






