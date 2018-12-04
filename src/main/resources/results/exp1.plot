set xlabel 'k'
set ylabel 'log(num clusters)'

set key top right
set key box

set logscale y

plot 'rnn_d31_ari.txt'  using 1:2 title 'd31' with lines, \
     'rnn_r15_ari.txt'  using 1:2 title 'r15' with lines, \
     'rnn_aggr_ari.txt' using 1:2 title 'aggregation' with lines, \
     'rnn_path_ari.txt' using 1:2 title 'pathbased' with lines
