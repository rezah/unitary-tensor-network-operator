set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "n"   font ",40" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"

#set xtics ("1/5" .2, "1/8" .125,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics (-0.4886,-0.4950,-0.49580, -0.4968,-0.4980);


set key font ",35"
#set key spacing 6

set key box b r  
set output "Spectrum.eps"
p [127-10:127+10] "EnergyIter8.txt" u 1:3 t "ED"   with points  pointtype 7 lw 2 ps 6.0 lc rgb "#729fcf","EnergyIter81.txt" u 1:2 t "vTN: lay=1"   with points  pointtype 9 lw 2 ps 6.0 lc rgb "#a40000","EnergyIter82.txt" u 1:2 t "vTN: lay=2"   with points  pointtype 11 lw 2 ps 6.0 lc rgb "#75507b","EnergyIter83.txt" u 1:2 t "vTN: lay=3"   with points  pointtype 13 lw 2 ps 6.0 lc rgb "#f57900"#,"EnergyIter84.txt" u 1:2 t "vTN: lay=4"   with points  pointtype 15 lw 2 ps 6.0 lc rgb "#c17d11"
