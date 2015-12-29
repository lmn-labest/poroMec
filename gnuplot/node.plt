#!/bin/gnuplot
set grid
set xlabel "Dias"
set ylabel "MPa"

# tensao total
#plot 'exemplo1_mesh0_t1_stress_node_44.txt'   using ($2/86400):($5/1000000) with line title 'sigmaZZ '

#pause -1

#plot 'exemplo1_mesh0_t1_stress_node_44.txt'   using ($2/86400):($3/1000000) with line title 'sigmaXX '

#pause -1

#tensao efetiva de terzaghi
#plot 'exemplo1_mesh0_t1_stressE_node_44.txt'   using ($2/86400):($5/1000000) with line title 'sigmaZZ '

#pause -1

#plot 'exemplo1_mesh0_t1_stressE_node_44.txt'   using ($2/86400):($3/1000000) with line title 'sigmaXX '

#pause -1

#tensao efetiva de biot
#plot 'exemplo1_mesh0_t1_stressB_node_44.txt'   using ($2/86400):($5/1000000) with line title 'sigmaZZ '

#pause -1

#plot 'exemplo1_mesh0_t1_stressB_node_44.txt'   using ($2/86400):($3/1000000) with line title 'sigmaXX '

#pause -1

#exit 0

set ylabel "Desloc(m)"
plot 'solo2_bcg_up_node_9041.txt'   using ($2/86400):5 with line title 'z bcg',\
     'solo2_up_node_9041.txt'       using ($2/86400):5 with points title 'z  cg',\
     'solo1_up_node_1271.txt'       using ($2/86400):5 with line title 'z  cg'

pause -1

set ylabel "P(MPa)"
plot 'solo2_bcg_up_node_220.txt'  using ($2/86400):6 with line title 'p bcg',\
     'solo2_up_node_220.txt'      using ($2/86400):6 with points title 'p  cg',\
     'solo1_up_node_61.txt'       using ($2/86400):6 with line title 'p  cg'
pause -1

exit 0
#set yrange[-3.0:3.0]
#plot 'exemplo1_mesh0_t1_node_44.txt'   using ($2/86400):($3/10000000) with line title 'pres 1'
#plot 'exemplo1_mesh0_t1_node_44.txt'   using ($2/86400):($6/1000000) with line title 'pres 1'
#plot 'exemplo1_mesh0_t1_node_44.txt'   using ($2/86400):($6/10000000) with line title 'pres 1',\
#     'exemplo1_mesh0_t2_node_44.txt' using ($2/86400):($6/10000000) with line title 'pres 2',\
#     'exemplo1_mesh0_t3_node_44.txt' using ($2/86400):($6/10000000) with line title 'pres 3'


      
pause -1
