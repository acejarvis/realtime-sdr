# example.gnuplut : configuration for plotting (change as needed)

reset                                   # reset
set size ratio 0.2                      # set relative size of plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 1,1 scale 1.0,1.0  # set two plots for this figure

set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 2

# time domain
#set ylabel 'Sample value'               # set y-axis label
#set xlabel 'Sample #'                   # set x-axis label
#set yrange [-0.007:0.007]                       # set y plot range
#set xrange [0:608]                      # set x plot range
#plot '../data/rrc_output_before.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle

#plot '../data/rrc_output_i.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle


# freq domain (Fourier)
#set ylabel 'Spectrum (Mag)'              # set y-axis label
#set xlabel 'Frequency bin'               # set x-axis label
#set yrange [-0.1:0.1]                    # set y plot range
#set xrange [0:255]                       # set x plot range
#plot '../data/demod_freq.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle

# rrc_out

 set yrange [-1:1]                       # set y plot range
 set xrange [-1:1]                       # set x plot range
# add your own .dat file for PSD as part of the take-home
plot '../data/rrc_output.dat' with points pt 7 ps 1

#plot '../data/sample.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle

#set xrange [300:400]                       # set x plot range
#plot '../data/sample.dat' with lines linestyle 1, \
 '../data/sample1.dat' with lines linestyle 2

