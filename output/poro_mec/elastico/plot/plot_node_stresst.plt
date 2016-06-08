#!/bin/gnuplot
set term png
set grid
set xlabel "Dias"

#set output 'sigmazz_total.png'
set title 'Tensão Total na base'
set ylabel "Tensão Total SigmaZZ (MPa)"
plot 'exemplo1_mesh0_t2_stress_node_44.txt'        using ($2/86400):5 with line title 'h1',\
     'solo1_stress_node_61.txt'                    using ($2/86400):5 with line title 'h2',\
     'solo_cilindro_hexa_1_stress_node_633.txt'    using ($2/86400):5 with line title 'h3',\
	 'exemplo1_mesh_t2_tetra_stress_node_44.txt'   using ($2/86400):5 with line title 't1',\
	 'solo1_tetra_stress_node_61.txt'              using ($2/86400):5 with line title 't2',\
	 'solo_cilindro_tetra_1_stress_node_747.txt'   using ($2/86400):5 with line title 't3'
	 
pause -1

	 
#set output 'sigmaxx_total.png'
set title 'Tensão Total na base'
set ylabel "Tensão Total SigmaXX (MPa)"
plot 'exemplo1_mesh0_t2_stress_node_44.txt'        using ($2/86400):4 with line title 'h1',\
     'solo1_stress_node_61.txt'                    using ($2/86400):4 with line title 'h2',\
     'solo_cilindro_hexa_1_stress_node_633.txt'    using ($2/86400):4 with line title 'h3',\
	 'exemplo1_mesh_t2_tetra_stress_node_44.txt'   using ($2/86400):4 with line title 't1',\
	 'solo1_tetra_stress_node_61.txt'              using ($2/86400):4 with line title 't2',\
	 'solo_cilindro_tetra_1_stress_node_747.txt'   using ($2/86400):4 with line title 't3'
	 
exit 0
