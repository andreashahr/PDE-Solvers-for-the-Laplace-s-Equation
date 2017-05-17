set output 'NEW_mv.png'
set terminal png size 800,500 enhanced font "Helvetica,16"
stats "filedata.out" matrix

set title "multigrid: gridSize=12, numIter=1252"
unset key
set xrange [0:*]
set yrange [0:*]
set zrange [0:STATS_max]
unset ytics
unset xtics
unset ztics
set cbtics (STATS_min,STATS_max)
set xyplane at 0
set border 4095
set format z "%1.4g"
set format cb "%1.4g"
set hidden3d front
set pm3d at b
set palette grey

splot "filedata.out" matrix lt -1
