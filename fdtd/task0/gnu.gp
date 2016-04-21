set terminal gif animate delay 5 loop 0
set output 'ehplot.gif'
set xrange[0:299]
set yrange[-1:1]
do for [i=0:318] {
    set title sprintf('t = %d', i)
    plot "eplot" every :::i::i using 2:3 title "E-field" with lines lw 2,     "hplot" every :::i::i using 2:3 title "H-field" with lines
}
set output
set output 'eplotfinal.gif'
set autoscale y
i = 318
set title 'at t = 319'
    plot "eplot" every :::i::i using 2:3 title "E-field" with lines
set output