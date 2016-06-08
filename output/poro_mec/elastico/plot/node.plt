#!/bin/gnuplot
set grid
set xlabel "Dias"

set title 'Deslocamento no topo'
set ylabel "z (metros)"
plot 'exemplo1_mesh0_t2_up_node_1.txt'         using ($2/86400):5 with line title 'z dominio 1',\
     'solo1_up_node_1271.txt'                  using ($2/86400):5 with line title 'z dominio 2',\
     'solo_cilindro_hexa_1_up_node_67.txt'     using ($2/86400):5 with line title 'z dominio 3'
pause -1

set title 'Delta de pressão na base'
set ylabel "deltaP (MPa)"
plot 'exemplo1_mesh0_t2_up_node_44.txt'     using ($2/86400):($6/1000000) with line title 'dominio 1',\
     'solo1_up_node_61.txt'               using ($2/86400):($6/1000000) with line title 'dominio 2',\
     'solo_cilindro_hexa_1_up_node_633.txt' using ($2/86400):($6/1000000) with line title 'dominio 3'


pause -1

# tensao total
set title 'Tensão Total na base'
set ylabel "Tensão Total sigmaZZ (MPa)"
plot 'exemplo1_mesh0_t2_stress_node_44.txt'     using ($2/86400):($5/1000000) with line title 'dominio 1',\
     'solo1_stress_node_61.txt'                 using ($2/86400):($5/1000000) with line title 'dominio 2',\
     'solo_cilindro_hexa_1_stress_node_633.txt' using ($2/86400):($5/1000000) with line title 'dominio 3'

pause -1

set ylabel "Tensão Total sigmaXX (MPa)"
plot 'exemplo1_mesh0_t2_stress_node_44.txt'       using ($2/86400):($3/1000000) with line title 'dominio 1',\
     'solo1_stress_node_61.txt'                   using ($2/86400):($3/1000000) with line title 'dominio 2',\
     'solo_cilindro_hexa_1_stress_node_633.txt'   using ($2/86400):($3/1000000) with line title 'dominio 3'

pause -1

#tensao efetiva de Terzaghi
set title 'Tensão de Tergaghi na base'
set ylabel "Tensão efetiva de Terzaghi sigmaZZ (MPa)"
plot 'exemplo1_mesh0_t2_stressE_node_44.txt'       using ($2/86400):($5/1000000) with line title 'dominio 1',\
     'solo1_stressE_node_61.txt'                   using ($2/86400):($5/1000000) with line title 'dominio 2',\
     'solo_cilindro_hexa_1_stressE_node_633.txt'   using ($2/86400):($5/1000000) with line title 'dominio 3'

pause -1


set ylabel "Tensão efetiva de Terzaghi sigmaXX (MPa)"
plot 'exemplo1_mesh0_t2_stressE_node_44.txt'        using ($2/86400):($3/1000000) with line title 'dominio 1',\
     'solo1_stressE_node_61.txt'                    using ($2/86400):($3/1000000) with line title 'dominio 2',\
     'solo_cilindro_hexa_1_stressE_node_633.txt'    using ($2/86400):($3/1000000) with line title 'dominio 3'

pause -1

#tensao efetiva de biot
set title 'Tensão de Biot na base'
set ylabel "Tensão efetiva de Biot sigmaZZ(MPa)"
plot 'exemplo1_mesh0_t2_stressB_node_44.txt'       using ($2/86400):($5/1000000) with line title 'dominio 1',\
     'solo1_stressB_node_61.txt'                   using ($2/86400):($5/1000000) with line title 'dominio 2',\
     'solo_cilindro_hexa_1_stressB_node_633.txt'   using ($2/86400):($5/1000000) with line title 'dominio 3'

pause -1

set ylabel "Tensão efetiva de Biot sigmaXX(MPa)"
plot 'exemplo1_mesh0_t2_stressB_node_44.txt'       using ($2/86400):($3/1000000) with line title 'dominio 1',\
     'solo1_stressB_node_61.txt'                   using ($2/86400):($3/1000000) with line title 'dominio 2',\
     'solo_cilindro_hexa_1_stressB_node_633.txt'   using ($2/86400):($3/1000000) with line title 'dominio 3'

pause -1



exit 0
