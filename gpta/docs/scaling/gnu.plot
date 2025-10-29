set term pngcairo enha lw 3 size 1440,1050 font "Arial,26"

set xlabel 'Total number of CPUs used'
set ylabel 'Time [s]'
set cblabel 'log_2 (N_{omp})'
unset key

set g
set log y

#==> res_3 <==
#  #Proc     MPI  OpenMP Reading    Comm  Action   Neigh   Total
#      1       1       1  0.124   0.000   0.450   4.295   4.967
t0 = 4.295
rc=3
set o 'fig_'.rc.'.png'
set title "R_{cut} = ".rc."\305"
p t0/x lc 0 , t0/(0.5*x+0.5) lc 0 dt "-" , t0/(0.25*x+0.75) lc 0 dt "." \
, 'res_'.rc u 1:($2==1 ? $7 : 1/0) lc "red", '' u 1:7:(log($3)/log(2)) w p pt 7 ps 4 palette ,'' u 1:4:(log($3)/log(2)) w p pt 2 ps 3 palette 

#  #Proc     MPI  OpenMP Reading    Comm  Action   Neigh   Total
#      1       1       1  0.110   0.000   0.854   6.002   7.062
t0 = 6.002
rc=4
set o 'fig_'.rc.'.png'
set title "R_{cut} = ".rc."\305"
p t0/x lc 0 , t0/(0.5*x+0.5) lc 0 dt "-" , t0/(0.25*x+0.75) lc 0 dt "." \
, 'res_'.rc u 1:($2==1 ? $7 : 1/0) lc "red", '' u 1:7:(log($3)/log(2)) w p pt 7 ps 4 palette ,'' u 1:4:(log($3)/log(2)) w p pt 2 ps 3 palette 

#==> res_6 <==
#  #Proc     MPI  OpenMP Reading    Comm  Action   Neigh   Total
#      1       1       1  0.221   0.000   4.115  19.311  23.752
t0 = 19.311
rc=6
set o 'fig_'.rc.'.png'
set title "R_{cut} = ".rc."\305"
p t0/x lc 0 , t0/(0.5*x+0.5) lc 0 dt "-" , t0/(0.25*x+0.75) lc 0 dt "." \
, 'res_'.rc u 1:($2==1 ? $7 : 1/0) lc "red", '' u 1:7:(log($3)/log(2)) w p pt 7 ps 4 palette ,'' u 1:4:(log($3)/log(2)) w p pt 2 ps 3 palette 

#==> res_10 <==
#  #Proc     MPI  OpenMP Reading    Comm  Action   Neigh   Total
#      1       1       1  0.474   0.000  15.493  60.947  77.027
t0 = 60.947
rc=10
set o 'fig_'.rc.'.png'
set title "R_{cut} = ".rc."\305"
p t0/x lc 0 , t0/(0.5*x+0.5) lc 0 dt "-" , t0/(0.25*x+0.75) lc 0 dt "." \
, 'res_'.rc u 1:($2==1 ? $7 : 1/0) lc "red", '' u 1:7:(log($3)/log(2)) w p pt 7 ps 4 palette ,'' u 1:4:(log($3)/log(2)) w p pt 2 ps 3 palette 

#==> res_20 <==
#  #Proc     MPI  OpenMP Reading    Comm  Action   Neigh   Total
#      1       1       1  0.316   0.000  98.146 288.104 386.695
t0 = 288.104
rc=20
set o 'fig_'.rc.'.png'
set title "R_{cut} = ".rc."\305"
p t0/x lc 0 , t0/(0.5*x+0.5) lc 0 dt "-" , t0/(0.25*x+0.75) lc 0 dt "." \
, 'res_'.rc u 1:($2==1 ? $7 : 1/0) lc "red", '' u 1:7:(log($3)/log(2)) w p pt 7 ps 4 palette ,'' u 1:4:(log($3)/log(2)) w p pt 2 ps 3 palette 
