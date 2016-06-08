#!/bin/gnuplot
set term png
set output 'pres.png'
set grid
set xlabel "Dias"

set title 'Pressão na base'
set ylabel "p(MPa)"
plot 'exemplo1_mesh0_t2_up_node_44.txt'        using ($2/86400):6 with line title 'h1',\
     'solo1_up_node_61.txt'                    using ($2/86400):6 with line title 'h2',\
     'solo_cilindro_hexa_1_up_node_633.txt'    using ($2/86400):6 with line title 'h3',\
	 'exemplo1_mesh_t2_tetra_up_node_44.txt'     using ($2/86400):6 with line title 't1',\
	 'solo1_tetra_up_node_1271.txt'             using ($2/86400):6 with line title 't2',\
	 'solo_cilindro_tetra_1_up_node_747.txt'    using ($2/86400):6 with line title 't3'

exit 0
