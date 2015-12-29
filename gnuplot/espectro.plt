#!/bin/gnuplot
set grid
set xlabel "numero"
set ylabel "autovalor"

#set yrange[*:0.0001]
#set log y

# tensao total
#set title "Matriz simetrica" 
#plot 'exemplo1_mesh0_dt1_sym_0_autovalor.txt' using 1:2 with line title 'dt1',\
#     'exemplo1_mesh0_dt2_sym_0_autovalor.txt' using 1:2 with line title 'dt2',\
#     'exemplo1_mesh0_dt3_sym_0_autovalor.txt' using 1:2 with line title 'dt3'

#pause -1

#set title "Matriz simetrica x nap simetrica" 
#plot 'exemplo1_mesh0_dt1_un_0_autovalor.txt'  using 1:2 with line title 'unSym dt1',\
#     'exemplo1_mesh0_dt1_sym_0_autovalor.txt' using 1:2 with line title 'sym dt1'

#pause -1


#set log y
#set title "Matriz nao-simetrica" 
#plot 'exemplo1_mesh0_dt1_un_0_autovalor.txt' using 1:2 with line title 'dt1',\
#     'exemplo1_mesh0_dt2_un_0_autovalor.txt' using 1:2 with line title 'dt2',\
#     'exemplo1_mesh0_dt3_un_0_autovalor.txt' using 1:2 with line title 'dt3'

#pause -1

set format y '%1.0e'
set yrange[1.e-10:*]
set log y
set ylabel "valores singulares"
set title "Matriz nao-simetrica" 
plot 'exemplo1_mesh0_dt1_un_0_svd.txt' using 1:2 with points title 'dt1',\
     'exemplo1_mesh0_dt2_un_0_svd.txt' using 1:2 with points title 'dt2',\
     'exemplo1_mesh0_dt3_un_0_svd.txt' using 1:2 with points title 'dt3',\
     'exemplo1_mesh0_dt4_un_0_svd.txt' using 1:2 with points title 'dt4',\
     'exemplo1_mesh0_dt5_un_0_svd.txt' using 1:2 with points title 'dt5'

pause -1

set ylabel "valores singulares"
set title "Matriz simetrica" 
plot 'exemplo1_mesh0_dt1_sym_0_svd.txt' using 1:2 with points title 'dt1',\
     'exemplo1_mesh0_dt2_sym_0_svd.txt' using 1:2 with points title 'dt2',\
     'exemplo1_mesh0_dt3_sym_0_svd.txt' using 1:2 with points title 'dt3',\
     'exemplo1_mesh0_dt4_sym_0_svd.txt' using 1:2 with points title 'dt4',\
     'exemplo1_mesh0_dt5_sym_0_svd.txt' using 1:2 with points title 'dt5'

pause -1


set ylabel "valores singulares"
set title "Matriz dt1" 
plot 'exemplo1_mesh0_dt1_un_0_svd.txt'  using 1:2 with points title 'unsym',\
     'exemplo1_mesh0_dt1_sym_0_svd.txt' using 1:2 with points title 'sym  '  
pause -1

set ylabel "valores singulares"
set title "Matriz dt2" 
plot 'exemplo1_mesh0_dt2_un_0_svd.txt'  using 1:2 with points title 'unsym',\
     'exemplo1_mesh0_dt2_sym_0_svd.txt' using 1:2 with points title 'sym  '  
pause -1

set ylabel "valores singulares"
set title "Matriz dt3" 
plot 'exemplo1_mesh0_dt3_un_0_svd.txt'  using 1:2 with points title 'unsym',\
     'exemplo1_mesh0_dt3_sym_0_svd.txt' using 1:2 with points title 'sym  '  
pause -1

set ylabel "valores singulares"
set title "Matriz dt4" 
plot 'exemplo1_mesh0_dt4_un_0_svd.txt'  using 1:2 with points title 'unsym',\
     'exemplo1_mesh0_dt4_sym_0_svd.txt' using 1:2 with points title 'sym  '  
pause -1

set ylabel "valores singulares"
set title "Matriz dt5" 
plot 'exemplo1_mesh0_dt5_un_0_svd.txt'  using 1:2 with points title 'unsym',\
     'exemplo1_mesh0_dt5_sym_0_svd.txt' using 1:2 with points title 'sym  '  
pause -1
