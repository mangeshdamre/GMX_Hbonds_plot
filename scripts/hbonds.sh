toilet -f mono12 -F metal HBONDS
echo "
|===========================================|
| Run this script in MD_Traj_Protein Folder |
|===========================================|";

echo "===================";
echo "Name of EM gro file"; 
echo "===================";
ls *.gro
echo "===================";
read GRO;

echo "===================";
echo "Name of XTC file"; 
echo "===================";
ls *.xtc
echo "===================";
read XTC;

echo "===================";
echo "Name of EM TPR file"; 
echo "===================";
ls *.tpr
echo "===================";
read TPR;

mkdir Hbonds
cd Hbonds		
gmx make_ndx -f ../$GRO -o hbonds.ndx << EOL
a 2708-2792
name 10 A

a 5852-5936
name 11 B

a 8996-9080
name 12 C

a 12140-12224
name 13 D

a 15284-15368
name 14 E

a 18428-18512
name 15 F

q
EOL

for chain in A B C D E F
do
mkdir $chain
if [ ! -f $chain/hbnum$chain.xvg ]; then
	gmx hbond -f ../$XTC \
	-s ../../$TPR \
	-n hbonds.ndx \
	-num $chain/hbnum$chain.xvg \
	-g $chain/hbond$chain.log \
	-ac $chain/hbac$chain.xvg \
	-dist $chain/hbdist$chain.xvg \
	-hx $chain/hbhelix$chain.xvg \
	-hbn $chain/hbond$chain.ndx \
	-hbm $chain/hbmap$chain.xpm \
	-a 20 -r 0.3 <<EOL
$chain
$chain
EOL
fi
done

gnuplot <<EOL
reset
set terminal pngcairo  background "#ffffff" enhanced font "Times-New-Roman-Bold,10" fontscale 1.0 size 500, 500 
set key on b r inside horizontal
set output 'hbdist.png'
set title "Hydrogen Bond Distribution"

set xlabel "Donor - Acceptor Distance (\305)" rotate parallel
set ylabel "" rotate parallel

set style line 1  lt 1 dt 1 pi 0 ps 1 # lt = color; pt = point type; dt = dash type; ps = size  
set style line 2  lt 2 dt 1 pi 0 ps 1
set style line 3  lt 4 dt 1 pi 0 ps 1
set style line 4  lt 5 dt 1 pi 0 ps 1
set style line 5  lt 6 dt 1 pi 0 ps 1
set style line 6  lt 7 dt 1 pi 0 ps 1
set style line 7  lt 8 dt 1 pi 0 ps 1
set style line 8  lt 10 dt 1 pi 0 ps 1

set grid xtics ytics mxtics mytics

p [2:4][] "A/hbdistA.xvg" u (column(1)*10):2 w l ls 1 t "A", \
  "B/hbdistB.xvg" u (column(1)*10):2 w l ls 2 lw 0.9 t "B", \
  "C/hbdistC.xvg" u (column(1)*10):2 w l ls 3 lw 0.9 t "C", \
  "D/hbdistD.xvg" u (column(1)*10):2 w l ls 4 lw 0.9 t "D", \
  "E/hbdistE.xvg" u (column(1)*10):2 w l ls 5 lw 0.9 t "E", \
  "F/hbdistF.xvg" u (column(1)*10):2 w l ls 6 lw 0.9 t "F"
EOL

gnuplot <<EOL
reset
set terminal pngcairo  background "#ffffff" enhanced font "Times-New-Roman-Bold,10" fontscale 1.0 size 500, 500 
set key on b r inside horizontal
set output 'hbnum.png'
set multiplot layout 3,2
set title "Hydrogen Bonds"

#set xlabel "Time (ns)" rotate parallel
set ylabel "Number" rotate parallel

set style line 1  lt 1 dt 1 pi 0 ps 1 # lt = color; pt = point type; dt = dash type; ps = size  
set style line 2  lt 2 dt 1 pi 0 ps 1
set style line 3  lt 4 dt 1 pi 0 ps 1
set style line 4  lt 5 dt 1 pi 0 ps 1
set style line 5  lt 6 dt 1 pi 0 ps 1
set style line 6  lt 7 dt 1 pi 0 ps 1
set style line 7  lt 8 dt 1 pi 0 ps 1
set style line 8  lt 10 dt 1 pi 0 ps 1

set yrange [0:12]

set grid xtics ytics mxtics mytics
unset xtics
set label 10 "A" c at graph 0.90,0.89
set ylabel " "
p "A/hbnumA.xvg" u (column(1)/1000):2 w l ls 1 lw 0.9 t ""

set title " "
set label 10 "D" c at graph 0.90,0.89
p   "D/hbnumD.xvg" u (column(1)/1000):2 w l ls 2 lw 0.9 t ""

unset title
set ylabel "# Bonds"
set label 10 "B" c at graph 0.90,0.89
p   "B/hbnumB.xvg" u (column(1)/1000):2 w l ls 3 lw 0.9 t ""

set ylabel " "
set label 10 "E" c at graph 0.90,0.89
p   "E/hbnumE.xvg" u (column(1)/1000):2 w l ls 4 lw 0.9 t ""

set ylabel " "
set xtics 10
set label 10 "C" c at graph 0.905,0.89
p   "C/hbnumC.xvg" u (column(1)/1000):2 w l ls 5 lw 0.9 t ""

set ylabel " "
set label 10 "F" c at graph 0.90,0.89
p   "F/hbnumF.xvg" u (column(1)/1000):2 w l ls 6 lw 0.9 t ""
unset multiplot
EOL

gnuplot <<EOL
reset
set terminal pngcairo  background "#ffffff" enhanced font "Times-New-Roman-Bold,10" fontscale 1.0 size 500, 500 
set key on b r inside horizontal
set output 'hbnum_within0.3nm.png'
set title "Hydrogen Bonds"

set xlabel "Time (ns)" rotate parallel
set ylabel "Number" rotate parallel

set style line 1  lt 1 dt 1 pi 0 ps 1 # lt = color; pt = point type; dt = dash type; ps = size  
set style line 2  lt 2 dt 1 pi 0 ps 1
set style line 3  lt 4 dt 1 pi 0 ps 1
set style line 4  lt 5 dt 1 pi 0 ps 1
set style line 5  lt 6 dt 1 pi 0 ps 1
set style line 6  lt 7 dt 1 pi 0 ps 1
set style line 7  lt 8 dt 1 pi 0 ps 1
set style line 8  lt 10 dt 1 pi 0 ps 1

set grid xtics ytics mxtics mytics

p "A/hbnumA.xvg" u (column(1)/1000):3 w l ls 1 lw 0.9 t "A", \
  "B/hbnumB.xvg" u (column(1)/1000):3 w l ls 2 lw 0.9 t "B", \
  "C/hbnumC.xvg" u (column(1)/1000):3 w l ls 3 lw 0.9 t "C", \
  "D/hbnumD.xvg" u (column(1)/1000):3 w l ls 4 lw 0.9 t "D", \
  "E/hbnumE.xvg" u (column(1)/1000):3 w l ls 5 lw 0.9 t "E", \
  "F/hbnumF.xvg" u (column(1)/1000):3 w l ls 6 lw 0.9 t "F"
EOL

gnuplot <<EOL
reset
set terminal pngcairo  background "#ffffff" enhanced font "Times-New-Roman-Bold,10" fontscale 1.0 size 500, 500 
set key on b r inside horizontal
set output 'hbnum_smoothBezier.png'
set multiplot layout 3,2
set title "Hydrogen Bonds"

#set xlabel "Time (ns)" rotate parallel
set ylabel "Number" rotate parallel

set style line 1  lt 1 dt 1 pi 0 ps 1 # lt = color; pt = point type; dt = dash type; ps = size  
set style line 2  lt 2 dt 1 pi 0 ps 1
set style line 3  lt 4 dt 1 pi 0 ps 1
set style line 4  lt 5 dt 1 pi 0 ps 1
set style line 5  lt 6 dt 1 pi 0 ps 1
set style line 6  lt 7 dt 1 pi 0 ps 1
set style line 7  lt 8 dt 1 pi 0 ps 1
set style line 8  lt 10 dt 1 pi 0 ps 1

set yrange [0:8]

set grid xtics ytics mxtics mytics
unset xtics
set label 10 "A" c at graph 0.90,0.89
set ylabel " "
p "A/hbnumA.xvg" u (column(1)/1000):2 smooth bezier ls 1 lw 2 t ""

set title " "
set label 10 "D" c at graph 0.90,0.89
p "D/hbnumD.xvg" u (column(1)/1000):2 smooth bezier ls 2 lw 2 t ""

unset title
set ylabel "# Bonds"
set label 10 "B" c at graph 0.90,0.89
p "B/hbnumB.xvg" u (column(1)/1000):2 smooth bezier ls 3 lw 2 t ""

set ylabel " "
set label 10 "E" c at graph 0.90,0.89
p "E/hbnumE.xvg" u (column(1)/1000):2 smooth bezier ls 4 lw 2 t ""

set ylabel " "
set xtics 10
set label 10 "C" c at graph 0.905,0.89
p "C/hbnumC.xvg" u (column(1)/1000):2 smooth bezier ls 5 lw 2 t ""

set ylabel " "
set label 10 "F" c at graph 0.90,0.89
p "F/hbnumF.xvg" u (column(1)/1000):2 smooth bezier ls 6 lw 2 t ""

unset multiplot
EOL

gnuplot <<EOL
reset
set terminal pngcairo  background "#ffffff" enhanced font "Times-New-Roman-Bold,10" fontscale 1.0 size 500, 500 
set key on b r inside horizontal
set output 'hbnum_smoothBezier_within0.3nm.png'
set multiplot layout 3,2
set title "Hydrogen Bonds"

#set xlabel "Time (ns)" rotate parallel
set ylabel "Number" rotate parallel

set style line 1  lt 1 dt 1 pi 0 ps 1 # lt = color; pt = point type; dt = dash type; ps = size  
set style line 2  lt 2 dt 1 pi 0 ps 1
set style line 3  lt 4 dt 1 pi 0 ps 1
set style line 4  lt 5 dt 1 pi 0 ps 1
set style line 5  lt 6 dt 1 pi 0 ps 1
set style line 6  lt 7 dt 1 pi 0 ps 1
set style line 7  lt 8 dt 1 pi 0 ps 1
set style line 8  lt 10 dt 1 pi 0 ps 1

set yrange [30:50]
set grid xtics ytics mxtics mytics
unset xtics
set label 10 "A" c at graph 0.90,0.89
set ylabel " "
p "A/hbnumA.xvg" u (column(1)/1000):3 smooth bezier ls 1 lw 2 t ""

set title " "
set label 10 "D" c at graph 0.90,0.89
p "D/hbnumD.xvg" u (column(1)/1000):3 smooth bezier ls 2 lw 2 t ""

unset title
set ylabel "# Bonds"
set label 10 "B" c at graph 0.90,0.89
p "B/hbnumB.xvg" u (column(1)/1000):3 smooth bezier ls 3 lw 2 t ""

set ylabel " "
set label 10 "E" c at graph 0.90,0.89
p "E/hbnumE.xvg" u (column(1)/1000):3 smooth bezier ls 4 lw 2 t ""

set ylabel " "
set xtics 10
set label 10 "C" c at graph 0.905,0.89
p "C/hbnumC.xvg" u (column(1)/1000):3 smooth bezier ls 5 lw 2 t ""

set ylabel " "
set label 10 "F" c at graph 0.90,0.89
p "F/hbnumF.xvg" u (column(1)/1000):3 smooth bezier ls 6 lw 2 t ""

unset multiplot
EOL
