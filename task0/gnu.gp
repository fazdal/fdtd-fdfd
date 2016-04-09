set terminal gif animate delay 5 loop 0
set output 'eplot.gif'
set xrange[0:199]
set yrange[-1:1]
do for [i=0:499] {
    set title sprintf('t = %d', i)
    plot "eplot" every :::i::i using 2:3 title "E-field" with lines,     "hplot" every :::i::i using 2:3 title "H-field" with lines
}
set output
