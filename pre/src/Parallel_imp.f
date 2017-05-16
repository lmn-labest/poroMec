c *********************************************************************
c * DIV_WORK: divide os trabalhos no paralelo                         *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * ierr  - variavel do mpi                                           *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * ierr  - 1 (sem MPI)                                               *  
c *********************************************************************
      subroutine div_work(MPIW,ierr,nprcs,nnos,div,part,loop,rank)
      implicit none
      integer MPIW,ierr
      integer nprcs,nnos,part,div,loop
      integer rank
c =====================================================================
c
c === with MPI
#ifdef MPI
      include 'mpif.h' 
      integer status(MPI_STATUS_SIZE)
c ===      
      call MPI_div_work(MPIW,status,ierr,nprcs,nnos,div,part,loop,rank)
c =====================================================================
c
c === without MPI
#else
c ===
      div  = nprcs
      loop = nprcs
      part = 0
      ierr = 1
c =====================================================================      
#endif      
      return
      end
c *********************************************************************
c
c *********************************************************************
c * FINALIZE: finaliza as operacoes em paralelo                       *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * ierr  - variavel do mpi                                           *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * ierr  - variavel do mpi (saida de erro) (com MPI)                 *  
c * ierr  - 1 (sem MPI)                                               *  
c *********************************************************************
      subroutine finalize()
      implicit none
c === with MPI     
#ifdef MPI 
      include 'mpif.h'
#endif
c =====================================================================
c
c ===
       integer ierr
       ierr = 1
c =====================================================================
c
c === with MPI
#ifdef MPI
      call MPI_Finalize(ierr)
#endif
      stop   
c =====================================================================
      return
      end
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 14/05/2017                                    * 
c * COMM_MEF :  comnicacao do paralelo                                 *  
c * ------------------------------------------------------------------ *  
c * Parametros de entrada :                                            *  
c * ------------------------------------------------------------------ *  
c * MPIW    -  variavel dp mpi                                         *  
c * nnode   -  numero de nos                                           *  
c * numel   -  numero de elementos                                     *  
c * nen     -  numero de nos por elementos                             *  
c * numat   -  numero de materiais                                     *  
c * ndf     -  graus de liberdade do problema mecanico                 *  
c * ndft    -  graus de liberdade do problema termico                  *  
c * ndm     -  numero de dimensoes da malha                            *  
c * npi     -                                                          *  
c * lines   -  linhas de macros comandos depois do end mesh            *  
c * nlines  -  numero linhas de macros comandos depois do end mesh     *  
c * i_ix    -  ponteiro conectividade da malha com material            *  
c * i_ie    -                                                          *  
c * i_num   -                                                          *  
c * i_e     -                                                          *  
c * i_x     -  ponteiro coordenadas da malha com material              *  
c * i_id    -                                                          *  
c * i_nload -                                                          *  
c * i_eload -                                                          *  
c * i_f     - ponteiro para forcas do problema mecanico                *  
c * i_u     - ponteiro para deslocamentos iniciais mecanico            *  
c * i_v     -                                                          *  
c * i_a     -                                                          *  
c * i_idt   -                                                          *  
c * i_nloadt-                                                          *  
c * i_eloadt-                                                          *  
c * i_ft    - ponteiro para forcas do problema termico                 *  
c * i_ut0   - ponteiro para temperaturas inicias                       *  
c * i_vt    -                                                          *  
c * i_w     -                                                          *  
c * rank    - id do processo                                           *  
c * ------------------------------------------------------------------ *  
c * Parametros de saida :                                              *  
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * Obs:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine comm_mesh_mef(MPIW,nnode,numel,nen,numat,ndf,ndft
     1                        ,ndm,npi,lines,nlines,i_ix,i_ie,i_e,i_x
     2                        ,i_id,i_nload,i_eload,i_eloadp
     3                        ,i_f,i_u,i_v,i_a
     4                        ,i_idt,i_nloadt,i_eloadt
     5                        ,i_ft,i_ut0,i_vt,i_w
     6                        ,i_xl,i_lel,i_idl,i_nloadl
     7                        ,i_eloadl,i_eloadpl
     8                        ,i_fl,i_idtl,i_nloadtl,i_eloadtl,i_ftl
     9                        ,i_ut0l
     1                        ,i_nodedist,i_elmdist
     2                        ,rank,npes,rload,ncont)
c ===
      implicit none
c ... Mpi      
      integer MPIW
      integer rank,npes
      integer*8 i_nodedist,i_elmdist
c ... malha      
      integer nnode,numel,nen,numat,ndm,npi
      integer*8 i_ix,i_ie,i_e,i_x
      integer*8 i_xl,i_lel
c ... problema 1
      integer ndf  
      integer*8 i_id,i_nload,i_eload,i_eloadp,i_f,i_u,i_v,i_a
      integer*8 i_idl,i_nloadl,i_eloadl,i_eloadpl,i_fl
c ... problema 2
      integer ndft
      integer*8 i_idt,i_nloadt,i_eloadt,i_ft,i_ut0,i_vt,i_w
      integer*8 i_idtl,i_nloadtl,i_eloadtl,i_ftl,i_ut0l
c ... arquivo do mefpar      
      character lines(*)
      integer nlines,ncont
      logical rload(ncont)
c =====================================================================
c
c ===
#ifdef MPI
      if(npes .gt. 1) then
        call MPI_comunicate_mesh_mef(MPIW,nnode,numel,nen,numat,ndf,ndft
     .                              ,ndm,npi,lines,nlines,i_ix,i_ie
     .                              ,i_e,i_x,i_id,i_nload,i_eload 
     .                              ,i_f,i_u,i_v,i_a,i_idt,i_nloadt
     .                              ,i_eloadt,i_ft,i_ut0,i_vt,i_w
     .                              ,i_xl,i_lel,i_idl,i_nloadl,i_eloadl
     .                              ,i_fl
     .                              ,i_idtl,i_nloadtl,i_eloadtl,i_ftl
     .                              ,i_ut0l
     .                              ,i_nodedist,i_elmdist
     .                              ,rank,npes,rload,ncont)
       else
         i_xl       = i_x
         i_lel      = i_ix
         i_idl      = i_id 
         i_nloadl   = i_nload
         i_eloadl   = i_eload
         i_fl       = i_f    
         i_idtl     = i_idt  
         i_nloadtl  = i_nloadt 
         i_eloadtl  = i_eloadt 
         i_ftl      = i_ft   
         i_ut0l     = i_ut0  
         i_nodedist = 1
         i_elmdist  = 1
       endif
#else
      i_xl       = i_x
      i_lel      = i_ix
      i_idl      = i_id 
      i_nloadl   = i_nload
      i_eloadl   = i_eload
      i_eloadpl  = i_eloadp
      i_fl       = i_f    
      i_idtl     = i_idt  
      i_nloadtl  = i_nloadt 
      i_eloadtl  = i_eloadt 
      i_ftl      = i_ft   
      i_ut0l     = i_ut0  
      i_nodedist = 1
      i_elmdist  = 1
#endif     
c =====================================================================
c
c ===      
      return
      end
c *********************************************************************
c * COMM_METIS: comunicacao do paralelo                               *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * MPIW    -  variavel dp mpi                                        *  
c * nnode   -  numero de nos                                          *  
c * numel   -  numero de elementos                                    *  
c * nen     -  numero de nos por elementos                            *  
c * ntria3  -  numero de elementos tria3                              *  
c * nquad4  -  numero de elementos quad4                              *  
c * nquad4  -  numero de elementos tetra4                             *  
c * nheaxa8 -  numero de elementos hexa8                              *  
c * i_ix0   -  ponteiro para as conectividades sem material           *  
c * rank    -  id do processo                                         *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine comm_metis(MPIW,nnode,numel,nen,ntria3,nquad4
     .                     ,ntetra4,nhexa8,i_ix0,rank)
c === Variaveis
      implicit none
      integer MPIW,nnode,numel,nen
      integer ntria3,nquad4,ntetra4,nhexa8
      integer*8 i_ix0
      integer rank
c === with MPI      
#ifdef MPI  
       call Mpi_comunicate_metis(MPIW,nnode,numel,nen,ntria3,nquad4
     .                          ,ntetra4,nhexa8,i_ix0,rank)
#endif
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c * COMM_OPT : comunicacao do paralelo                                *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * MPIW    -  variavel dp mpi                                        *  
c * ierr    -  variavel dp mpi (saida de erro)                        *  
c * filein  -  nome do arquivo de entrada                             *  
c * fileout -  nome do arquivo de saida                               *  
c * ovlp    -  metodo de particionamento                              *  
c * vtkp    -  escrita do vtk                                         *  
c * mefp    -  escrita do mef par                                     *  
c * nprcs   -  numero de divisoes da malha                            *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine comm_opt(MPIW,filein,fileout,ovlp,vtkp,mefp,nprcs,ierr)
c === Variaveis
      implicit none
      character filein(*),fileout(*)
      logical ovlp,vtkp,mefp
      integer MPIW,rank,nnos,ierr,nprcs
c === with MPI      
#ifdef MPI      
      call MPI_comunicate_opt(MPIW,filein,fileout,ovlp,vtkp,mefp,nprcs
     .                       ,ierr)
#endif
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * START    : Inicia variaveis utilizados no programa                *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * MPIW  -  variavel dp mpi                                          *  
c * rank  -  variavel do mpi (id do processo)                         *  
c * nnos  -  variavel do mpi (numero de processos)                    *  
c * ierr  -  variavel do mpi (saida de erro)                          *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine start(MPIW,rank,nnos,ierr)
c === Variaveis
      implicit none
      integer MPIW,rank,nnos,ierr
c === with MPI      
#ifdef MPI      
      call Mpi_start(MPIW,rank,nnos,ierr)
c =====================================================================
c
c === without MPI
#else
      MPIW = 0
      rank = 0
      nnos = 1
      ierr = 0
#endif
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
