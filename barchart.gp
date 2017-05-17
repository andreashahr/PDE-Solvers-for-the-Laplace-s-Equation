set terminal png size 800,500 enhanced font "Helvetica,12"
set output 'barchart.png'

set style line 2 lc rgb 'black' lt 1 lw 1
set style data histogram
set style histogram cluster gap 1.4
set style fill pattern border -1
set boxwidth 0.9
set xtics format ""
set yrange [0:50]
set grid ytics

plot "bardata" using 2:xtic(1) title "sequential" ls 2, \
            '' using 3 title "1 thread" ls 2, \
            '' using 4 title "2 threads" ls 2, \
            '' using 5 title "3 threads" ls 2, \
	    '' using 6 title "4 threads" ls 2
