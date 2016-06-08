#!/bin/gnuplot
set term png
set output 'desloc.png'
set grid
set xlabel "Dias"

set title 'Deslocamento do topo'
set ylabel "z (metros)"
plot 'exemplo1_mesh0_t2_up_node_1.txt'        using ($2/86400):5 with line title 'h1',\
     'solo1_up_node_1271.txt'                 using ($2/86400):5 with line title 'h2',\
     'solo_cilindro_hexa_1_up_node_67.txt'    using ($2/86400):5 with line title 'h3',\
	 'exemplo1_mesh_t2_tetra_up_node_1.txt'   using ($2/86400):5 with line title 't1',\
	 'solo1_tetra_up_node_1271.txt'           using ($2/86400):5 with line title 't2',\
	 'solo_cilindro_tetra_1_up_node_610.txt'    using ($2/86400):5 with line title 't3'
pause -1


exit 0
