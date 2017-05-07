#ifdef MPI
c *********************************************************************
c * MPI_INIT :                                                        *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine Mpi_start(MPIW,rank,nnos,ierr)
c === Variaveis
      implicit none
      integer MPIW,rank,nnos,ierr
      include "mpif.h"
c =====================================================================
c
c === inicializando MPI
      call MPI_Init(ierr)
c ... obtendo meu rank
      MPIW = MPI_COMM_WORLD
      call MPI_Comm_rank(MPIW,rank,ierr)
c ... numero de prcs Mips      
      call MPI_Comm_size(MPIW,nnos,ierr)
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * MPI_COMUNICATE_MESH_MEF : Distribui as informacoes da malha entre *  
c * os processos                                                      *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * MPIW       - comunicador do mpi                                   *  
c * nnode      - numero de nos                                        *  
c * nen        - numero do nos por elemento                           *  
c * numat      - numero de materiais                                  *  
c * ndf        - grau de liberdade mecanico                           *  
c * ndft       - grau de liberdade termico                            *  
c * ndm        - numero de dimensoes                                  *  
c * npi        -                                                      *  
c * lines      - linhas lidas apos a macro end mesh                   *  
c * nlines     - numero de linha lidas apos a macro end mesh          *  
c * i_ix       - ponteiro para as conectividades                      *  
c * i_ie       - ponteiro para o tipo de elemento                     *  
c * i_e        - ponteiro para as propriendades                       *  
c * i_x        - ponteiro para as coordenadas                         *  
c * i_id       - ponteiro para as restricoes mecanicas                *  
c * i_nload    - ponteiro para cargas no nos                          *  
c * i_eload    - ponteiro para cargas nos elementos                   *  
c * i_f        - ponteiro para cargas                                 *  
c * i_u        - ponteiro para o deslocamento incial                  *  
c * i_v        - ponteiro para derivada 1 no tempo                    *  
c * i_a        - ponteiro para derivada 2 no tempo                    *  
c * i_idt      - ponteiro para restricoes termicas                    *  
c * i_nloadt   - ponteiro para cargas termicas                        *  
c * i_eloadt   - ponteiro para cargas termicas nos elementos          *  
c * i_ft       - ponteiro para fontes termicas nodais                 *  
c * i_ut0      - ponteiro para temperaturas iniciais                  *  
c * i_vt       - ponteiro para derivada primeira                      *  
c * i_w        - ponteiro para o campo de velicidades(conveccao)      *  
c * i_xl       -                                                      *  
c * i_lel      -                                                      *  
c * i_idl      -                                                      *  
c * i_nloadl   -                                                      *  
c * i_eloadl   -                                                      *  
c * i_fl       -                                                      *  
c * i_idtl     -                                                      *  
c * i_nloadtl  -                                                      *  
c * i_eloadtl  -                                                      *  
c * i_ftl      -                                                      *  
c * i_utl0     -                                                      *  
c * i_nodedist -                                                      *  
c * my_rank    - numero do meu processo                               *  
c * npes       - numero de processos mpi                              *  
c * rload      - restricoes lidas na rdat                             *  
c * ncont      - numero de restricoes lidas na rdat                   *  
c * ----------------------------------------------------------------- *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * i_xl       - ponteiro para coordenadas distribuidas               *  
c * i_lel      - ponteiro para conectividades distribuidas            *  
c * i_idl      - ponteiro para restricoes  distribuidas               *  
c * i_nloadl   - ponteiro para cargas nodais distribuidas             *  
c * i_eloadl   - ponteiro para cargas elmt   distribuidas             *  
c * i_fl       - ponteiro para cargas  distribuidas                   *  
c * i_idtl     - ponteiro para restricoes  distribuidas               *  
c * i_nloadtl  - ponteiro para cargas nodais distribuidas             *  
c * i_eloadtl  - ponteiro para cargas elmt   distribuidas             *  
c * i_ftl      - ponteiro para cargas distribuidas                    *  
c * i_utl0     - ponteiro para temperaturas iniciais distribuidas     *  
c * i_nodedist - ponteiro para distribuidas das propriedades nodais   *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine Mpi_comunicate_mesh_mef(MPIW,nnode,numel,nen,numat,ndf
     .                                  ,ndft,ndm,npi,lines,nlines,i_ix
     .                                  ,i_ie,i_e,i_x,i_id,i_nload
     .                                  ,i_eload,i_f,i_u,i_v,i_a
     .                                  ,i_idt,i_nloadt,i_eloadt
     .                                  ,i_ft,i_ut0,i_vt,i_w,i_xl,i_lel
     .                                  ,i_idl,i_nloadl,i_eloadl
     .                                  ,i_fl
     .                                  ,i_idtl,i_nloadtl,i_eloadtl
     .                                  ,i_ftl,i_ut0l
     .                                  ,i_nodedist,i_elmdist
     .                                  ,my_rank,npes,rload,ncont)
c === Variaveis  
      use Malloc
      implicit none
      include "mpif.h"
      include "termprop.fi"
      include "elementos.fi"
      include "string.fi"
      include "parallel.fi"
      include "load.fi"
c ... MPI      
      integer MPIW
      integer my_rank,npes,nno,nel,i
      real*8 df
c ... malha
      integer nnode,numel,nen
      integer ndf,ndft,ndm,numat,npi
c ... ponteiros       
      integer*8 i_ix,i_id,i_ie,i_e,i_x
      integer*8 i_nodedist,i_elmdist
      integer*8 i_xl,i_lel
c ... mecanico
      integer*8 i_nload,i_eload
      integer*8 i_f,i_u,i_v,i_a
      integer*8 i_idl,i_nloadl,i_eloadl,i_fl
c ... termico     
      integer*8 i_idt,i_nloadt,i_eloadt,i_ft,i_ut0,i_vt,i_w
      integer*8 i_idtl,i_nloadtl,i_eloadtl,i_ftl,i_ut0l
c ...      
      integer nlines
      character*200 lines(*)
      integer ncont
      logical rload(*)
c =====================================================================      
c
c === difundido dados inicias master -> escravos
      call MPI_Bcast(nlines ,1           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(lines  ,nlines*200  ,MPI_CHARACTER,0,MPIW,ierr)
      call MPI_Bcast(macros ,100*15      ,MPI_CHARACTER,0,MPIW,ierr)
      call MPI_Bcast(rload  ,ncont       ,MPI_LOGICAL  ,0,MPIW,ierr)
      call MPI_Bcast(nnode  ,1           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(numel  ,1           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(nen    ,1           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(numat  ,1           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(ndf    ,1           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(ndft   ,1           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(ndm    ,1           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(npi    ,1           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(nbar2  ,2           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(ntria3 ,2           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(ntetra4,2           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(nquad4 ,2           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(nhexa8 ,2           ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(load   ,2*numload   ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(eprop  ,prop*maxmate,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(fload,npar*nterms*numload ,MPI_DOUBLE_PRECISION
     .              ,0,MPIW,ierr)
      call MPI_Bcast(nprop,prop*ntterms*maxmate,MPI_DOUBLE_PRECISION
     .              ,0,MPIW,ierr)
      call MPI_Bcast(tad,maxpol*maxmate,MPI_DOUBLE_PRECISION
     .              ,0,MPIW,ierr)
      call MPI_Bcast(tempad,maxpol*maxmate,MPI_DOUBLE_PRECISION
     .              ,0,MPIW,ierr)
c =====================================================================
c
c === difundindo a malha master->escravos
      if(my_rank .ne. 0) then
        i_ie    = alloc_4('ie      ',      1,numat)
        i_e     = alloc_8('e       ',maxprop,numat)
      endif
c ...
      call MPI_Bcast(ia(i_ie)   ,numat        ,MPI_INTEGER         ,0
     .                ,MPIW,ierr)
c ...
      call MPI_Bcast(ia(i_e)    ,maxprop*numat ,MPI_DOUBLE_PRECISION,0
     .                ,MPIW,ierr)
c =====================================================================
c
c ===
c ... divisao das propriedades do elementos pelos processos
      i_elmdist = alloc_4('elmdist ',  1,npes+1)
      call mpi_split_prop(numel,nel,ia(i_elmdist),MPIW,npes,my_rank)
c     print*, (ia(i_elmdist+i),i=0,npes)
c ... elementos  
      i_lel     = alloc_4('lel     ',nen+1,nel)
      call mpi_dist_prop_int(ia(i_ix),ia(i_lel),ia(i_elmdist)
     .                      ,MPIW,nen+1,npes,my_rank)
c =====================================================================
c
c ===
c ... divisao das propriedades nodais pelos processos
      i_nodedist = alloc_4('nodedist',  1,npes+1)
      call mpi_split_prop(nnode,nno,ia(i_nodedist),MPIW,npes,my_rank)
c     print*,nnode,nno
c     print*, (ia(i_nodedist+i),i=0,npes)
c ... coordenadas
      i_xl       = alloc_8('lx      ',ndm,nno)
      call mpi_dist_prop_node_double(ia(i_x),ia(i_xl),ia(i_nodedist)
     .                              ,MPIW,ndm,npes,my_rank)
c .....................................................................
c
c ... restricoes mecanico
      if(rload(1))then
        i_id   =  locate('id      ')
        if(my_rank .gt. 0) i_id = 1
        i_idl  = alloc_4('idl     ',ndf,nno)
        call mpi_dist_prop_int(ia(i_id),ia(i_idl),ia(i_nodedist)
     .                             ,MPIW,ndf,npes,my_rank)
      endif 
c .....................................................................
c
c ... nodalloads mecanico
      if(rload(2))then
        i_nload =  locate('nload   ')
        if(my_rank .gt. 0) i_nload = 1
        i_nloadl= alloc_4('nloadl  ',ndf,nno)
        call mpi_dist_prop_int(ia(i_nload),ia(i_nloadl)
     .                             ,ia(i_nodedist),MPIW,ndf,npes
     .                             ,my_rank)
      endif
c .....................................................................
c
c ... eload mecanico
      if(rload(3))then
        i_eload =  locate('eload   ')
        if(my_rank .gt. 0) i_eload = 1
        i_eloadl= alloc_4('eloadl  ',  7,nel)
        call mpi_dist_prop_int(ia(i_eload),ia(i_eloadl)
     .                             ,ia(i_elmdist),MPIW,7,npes
     .                             ,my_rank)
      endif
c .....................................................................
c
c ... nodalforces mecanico
      if(rload(4))then
        i_f     =  locate('f       ')
        if(my_rank .gt. 0) i_f = 1
        i_fl    = alloc_8('fl      ',ndf,nno)
        call mpi_dist_prop_node_double(ia(i_f),ia(i_fl),ia(i_nodedist)
     .                                ,MPIW,ndf,npes,my_rank)
      endif
c .....................................................................
c
c ... restricoes termico
      if(rload(8))then
        i_idt   =  locate('idt     ')
        if(my_rank .gt. 0) i_idt = 1
        i_idtl  = alloc_4('idtl    ',ndft,nno)
        call mpi_dist_prop_int(ia(i_idt),ia(i_idtl),ia(i_nodedist)
     .                             ,MPIW,ndft,npes,my_rank)
      endif
c .....................................................................
c
c ...
      if(rload(9))then
        i_nloadt  =  locate('nloadt  ')
        if(my_rank .gt. 0) i_nloadt = 1
        i_nloadtl = alloc_4('nloadtl ',ndft,nno)
        call mpi_dist_prop_int(ia(i_nloadt),ia(i_nloadtl)
     .                        ,ia(i_nodedist),MPIW,ndft,npes
     .                        ,my_rank)
      endif
c .....................................................................
c
c ... carga termicas nos elementos
      if(rload(10))then
        i_eloadt= locate('eloadt  ')
        if(my_rank .gt. 0) i_eloadt = 1
        i_eloadtl = alloc_4('eloadtl ',7,nel)
        call mpi_dist_prop_int(ia(i_eloadt),ia(i_eloadtl)
     .                        ,ia(i_elmdist),MPIW,7,npes
     .                        ,my_rank)
      endif
c .....................................................................
c
c ... fontest termicas             
      if(rload(11))then
        i_ft    = locate('ft      ')
        if(my_rank .gt. 0) i_ft = 1
        i_ftl   = alloc_8('ftl     ',ndft,nno)
        call mpi_dist_prop_node_double(ia(i_ft),ia(i_ftl)
     .                                ,ia(i_nodedist),MPIW,ndft,npes
     .                                ,my_rank)
      endif
c .....................................................................
c
c ... condicoes iniciais termico
      if(rload(12))then
        i_ut0   = locate('ut0     ')
        if(my_rank .gt. 0) i_ut0 = 1
        i_ut0l  = alloc_8('ut0l    ',ndft,nno)
        call mpi_dist_prop_node_double(ia(i_ut0),ia(i_ut0l)
     .                                ,ia(i_nodedist)
     .                                ,MPIW,ndft,npes
     .                                ,my_rank)
      endif 
c .....................................................................
c
c ... liberando memoria
      if(my_rank .eq. 0) then
        i_ix    = dealloc('ix      ')
        i_x     = dealloc('x       ')
        if(rload(1))  i_id    = dealloc('id      ')
        if(rload(2))  i_nload = dealloc('nload   ')
        if(rload(3))  i_eload = dealloc('eload   ')
        if(rload(4))  i_f     = dealloc('f       ')
        if(rload(8))  i_idt   = dealloc('idt     ')
        if(rload(9))  i_nloadt= dealloc('nloadt  ')
        if(rload(10)) i_eloadt= dealloc('eloadt  ')
        if(rload(11)) i_ft    = dealloc('ft      ')
        if(rload(12)) i_ut0   = dealloc('ut0     ')
      endif
c ... recalculando vetores      
      i_nodedist = locate('nodedist')
      i_elmdist  = locate('elmdist ')
      i_xl       = locate('lx      ')
      i_lel      = locate('lel     ')
      i_ie       = locate('ie      ')
      i_e        = locate('e       ')
      if(rload(1))  i_idl      = locate('idl     ')
      if(rload(2))  i_nloadl   = locate('nloadl  ')
      if(rload(3))  i_eloadl   = locate('eloadl  ')                
      if(rload(4))  i_fl       = locate('fl      ')
      if(rload(8))  i_idtl     = locate('idtl    ')
      if(rload(9))  i_nloadtl  = locate('nloadtl ')
      if(rload(10)) i_eloadtl  = locate('eloadtl ')
      if(rload(11)) i_ftl      = locate('ftl     ')
      if(rload(12)) i_ut0l     = locate('ut0l    ')
c .....................................................................
c
c === alocar memoria para os esvravos
      if(my_rank .ne. 0) then
c       if( rload(3)) i_eload  = alloc_4('eload   ',     7,numel)
c       if(rload(10)) i_eloadt = alloc_4('eloadt  ',     7,numel)
        if(rload(13)) i_vt     = alloc_8('vt      ',  ndft,nnode)
      endif
c =====================================================================
c
c === difundindo a malha master->escravos
c ...
c      if(rload(3)) then
c        i_eload = locate('eload   ')
c        call MPI_Bcast(ia(i_eload),7*numel,MPI_INTEGER,0
c     .                           ,MPIW,ierr)
c       endif
c ...
c     if(rload(10))then
c       i_eloadt= locate('eloadt  ')
c       call MPI_Bcast(ia(i_eloadt),7*numel,MPI_INTEGER,0
c    .                            ,MPIW,ierr)
c     endif
c ... condicoes iniciais(derivada primeira)
      if(rload(13))then
        i_vt    = locate('vt      ')
        call MPI_Bcast(ia(i_vt),ndft*nnode
     .                ,MPI_DOUBLE_PRECISION,0,MPIW,ierr)
      endif
c =====================================================================      
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * MPI_COMUNICATE_OPT :                                              *  
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
      subroutine Mpi_comunicate_opt(MPIW,filein,fileout,ovlp,vtkp,mefp
     .                             ,nprcs,ierr)
c ===
      implicit none
      include "mpif.h"
c ... MPI      
      integer MPIW
      integer ierr
c ... opcoes do progema
      character filein(*),fileout(*)
      integer nprcs          
      logical vtkp,mefp,ovlp
c =====================================================================
c
c ===
      call MPI_Bcast(filein ,80    ,MPI_CHARACTER,0,MPIW,ierr)
      call MPI_Bcast(fileout,80    ,MPI_CHARACTER,0,MPIW,ierr)
      call MPI_Bcast(ovlp   ,1     ,MPI_LOGICAL  ,0,MPIW,ierr)
      call MPI_Bcast(vtkp   ,1     ,MPI_LOGICAL  ,0,MPIW,ierr)
      call MPI_Bcast(mefp   ,1     ,MPI_LOGICAL  ,0,MPIW,ierr)
      call MPI_Bcast(nprcs  ,1     ,MPI_INTEGER  ,0,MPIW,ierr)
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * MPI_TEMP : gera a soma dos tempos em todos os processos           *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * MPIW      - comunicados mpi                                       *  
c * tcomunica - tempo da comunicao por processo                       *  
c * tmetis    - tempo do metis por processo                           *  
c * tpart     - tempo do particionamento por processo                 *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * tc - soma do tempo da comunicacao                                 *  
c * tm - soma do tempo do metis                                       *  
c * tp - soma do tempo do particionamento                             *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine Mpi_temp(MPIW,tcomunica,tmetis,tpart,tc
     .                   ,tm,tp,ierr)
c ===
      implicit none
      include "mpif.h"
c ... MPI      
      integer MPIW
      integer ierr
c ...                    
      real*8  tcomunica,tmetis,tpart,tc,tm,tp
c =====================================================================
c
c ===
      call MPI_Reduce(tc       ,tcomunica,1,MPI_DOUBLE_PRECISION
     .               ,MPI_SUM,0,MPIW,ierr)
      call MPI_Reduce(tm       ,tmetis   ,1,MPI_DOUBLE_PRECISION
     .               ,MPI_SUM,0,MPIW,ierr)
      call MPI_Reduce(tp       ,tpart    ,1,MPI_DOUBLE_PRECISION
     .               ,MPI_SUM,0,MPIW,ierr)
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * MPI_COMUNICATE_METIS :                                            *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine Mpi_comunicate_metis(MPIW,nnode,numel,nen,ntria3,nquad4
     .                                ,ntetra4,nhexa8,my_rank)
c ===
      use Malloc
      implicit none
      include "mpif.h"
c ... MPI      
      integer MPIW
      integer my_rank,ierr
c ... malha
      integer nnode,numel,nen
      integer ntria3(*),nquad4(*),ntetra4(*),nhexa8(*)
c =====================================================================       
c
c ====
      call MPI_Bcast(nnode  ,1     ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(numel  ,1     ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(nen    ,1     ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(ntria3 ,2     ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(nquad4 ,1     ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(ntetra4,1     ,MPI_INTEGER  ,0,MPIW,ierr)
      call MPI_Bcast(nhexa8 ,1     ,MPI_INTEGER  ,0,MPIW,ierr)
c =====================================================================       
c
c ===
      return
      end
c *********************************************************************      
c
c *********************************************************************
c * MPI_DIV_WORK : divide o particionamentos pelo porcessos           *
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * MPIW  -                                                           *  
c * status-                                                           *  
c * ierr  -                                                           *  
c * ndiv  - numero de divisoes                                        *  
c * npes  - numero de processo MPI                                    *  
c * div   -                                                           *  
c * part  -                                                           *  
c * mype  - numero do processo MPI                                    *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * div   - numero de divisoes atribuido ao processo                  *  
c * part  - numero da particao inicial do processo                    *  
c * loop  - numero de loop no mapa                                    *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine MPI_div_work(MPIW,status,ierr,ndiv,npes,div,part,loop
     .                       ,mype)
c ===
      use Malloc
      implicit none
      include 'mpif.h'
c ... MPI      
      integer MPIW,ierr,tag,source
      integer status(*)
c ... Particionamento      
      integer ndiv,npes,part,div,loop
      integer mype
      integer*8 i_dpn,i_pp
      integer dv,md
      integer i
c ...
      data tag/10/
c .....................................................................
c =====================================================================
c
c ===
      div  = 0
      part = 0
      if(mype .eq. 0)then
c ...      
        i_dpn = alloc_4('dpn     ',1,npes)
        i_pp  = alloc_4('pp      ',1,npes)
c .....................................................................
c
c ...
        dv = ndiv/npes
        md = mod(ndiv,npes)
c .....................................................................
c
c ... div de particoes por no
        do i = 1,npes
          ia(i_dpn+i-1) = dv
        enddo
        do while(md .gt. 0)
          do i = 1,npes
            ia(i_dpn+i-1) = ia(i_dpn+i-1) + 1
                       md = md - 1
            if(md .eq. 0) goto 100          
          enddo
        enddo
c .....................................................................
 100    continue
c ... div dos ranks por no
        ia(i_pp) = 0
        do i = 2,npes
          ia(i_pp+i-1) = ia(i_dpn+i-2) + ia(i_pp+i-2)
                    md = md - 1
        enddo
c       do i = 1, nnos
c         print*,ia(i_dpn+i-1),ia(i_pp+i-1)
c       enddo
        div   = ia(i_dpn)
        part  = ia(i_pp)
c ... envindo para os outros nos        
        do i = 1,npes-1
          call MPI_Send(ia(i_dpn+i),1,MPI_INTEGER,i,tag,MPIW,ierr)
          call MPI_Send(ia(i_pp+i),1,MPI_INTEGER,i,tag,MPIW,ierr)
        enddo
c .....................................................................        
        i_dpn = dealloc('dpn     ')
        i_pp  = dealloc('pp      ')
c ... recebendo        
      else
        source = 0
        call MPI_Recv(div,1,MPI_INTEGER,source,tag,MPIW,status
     .               ,ierr)
        call MPI_Recv(part,1,MPI_INTEGER,source,tag,MPIW,status
     .               ,ierr)
      endif
c ... 
      call Mpi_Allreduce(div,loop,1,MPI_INTEGER,MPI_MAX,MPIW,ierr)
c .....................................................................
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * MPI_SPLIT_VECTOR_BCAST : divide um vetor para mpi_bcast           *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * x     - vetor a ser comunicado                                    *  
c * n     - dimensao total do vetor                                   *  
c * nparts- numero de partes q o vetor ira ser enviado                *  
c * cod   - 1 - integer 2 - double                                    *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine mpi_split_vector_bcast(ix,dx,n,nparts,MPIW,cod)
      implicit none
      include 'mpif.h'
      integer n
      integer ix(*)
      real*8  dx(*)
      integer MPIW,i,ierr,dv,md,nparts
      integer cod
c ===  
      dv  = n/nparts
      md  = mod(n,nparts)
      do i = 0,nparts-1
        if( cod .eq. 1) then
          call MPI_Bcast(ix(dv*i+1)   ,dv,MPI_INTEGER,0,MPIW,ierr)
        elseif( cod .eq. 2) then 
          call MPI_Bcast(dx(dv*i+1)   ,dv,MPI_DOUBLE_PRECISION,0,MPIW
     .                  ,ierr)
        endif  
      enddo
c ... resto do vetor      
      if( md .gt. 0) then
        if( cod .eq. 1) then
          call MPI_Bcast(ix(dv*nparts+1)   ,md,MPI_INTEGER,0,MPIW,ierr)
        elseif( cod .eq. 2) then 
          call MPI_Bcast(dx(dv*nparts+1)   ,md,MPI_DOUBLE_PRECISION,0
     .                  ,MPIW,ierr)
        endif                
      endif
c =====================================================================
c
c ===  
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * MPI_SPLIT_PROP: distribuicao das propriedade entre  os processos  *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * n    - numero total das propriedades a serem distribuidos         *  
c * ln   -                                                            *  
c * dist -                                                            *  
c * MPIW -                                                            *  
c * mype - numero do meu processo(MPI)                                *  
c * npes - numero de processos (MPI)                                  *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * dist - distribuicao da propriedades ( todos o proc)               *  
c * ln   - numero local de propriedades                               *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine mpi_split_prop(n,ln,dist,MPIW,npes,mype)
      implicit none
      include 'mpif.h'
      include 'time.fi'
      integer npes,mype
      integer n,ln
      integer dist(*)
      integer MPIW,i,j,k,ierr
c ===  
      if(mype .eq. 0) then
        dist(1) = 0
        j = n
        do i = 1, npes
                  k = j/(npes-i+1)
          dist(i+1) = dist(i) + k
                  j = j - k
        enddo
      endif
c .....................................................................
c
c ... comunicando a distribuicao
c
      tcomunica = getime() - tcomunica
      call MPI_Bcast(dist,npes+1,MPI_INTEGER,0,MPIW,ierr)
      tcomunica = getime() - tcomunica
      ln = dist(mype+2)-dist(mype+1)
c =====================================================================
c
c ===  
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * MPI_SPLIT_PROP_NODE_DOUBLE: distribuicao das propriedade nodais   *
c * entre os processos (double)                                       *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * x        - coordenadas                                            *  
c * xl       -                                                        *  
c * nodedist - distribuicao dos propriedade nodais pelos processos    *  
c * MPIW     -                                                        *  
c * mype     - numero do meu processo(MPI)                            *  
c * npes     - numero de processos (MPI)                              *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * xl       - distribuicao das coordenadas                           *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine mpi_dist_prop_node_double(x,xl,nodedist,MPIW,ndm,npes
     .                                    ,mype)
      implicit none
      include 'mpif.h'
      include 'time.fi'
      integer status(MPI_STATUS_SIZE)
      integer npes,mype
      integer nnode,ndm
      integer nodedist(*)
      real*8 x(ndm,*),xl(ndm,*)
      integer MPIW,nno,ierr,your_nno,nnoi,i,j,pe
c ===
c     call MPI_Barrier(MPIW,ierr)
      tcomunica = getime() - tcomunica
      nno = nodedist(mype+2)-nodedist(mype+1)
      if(mype .eq. 0) then
        do i = 1, nno
          do j = 1, ndm
            xl(j,i) = x(j,i)
          enddo
        enddo
        do pe = 2, npes
          your_nno = nodedist(pe+1)-nodedist(pe)
          nnoi = nodedist(pe) + 1
          call MPI_Send(x(1,nnoi),your_nno*ndm,MPI_DOUBLE_PRECISION
     .                 ,pe-1,0,MPIW,ierr) 
        enddo
      else
        call MPI_Recv(xl,nno*ndm,MPI_DOUBLE_PRECISION,0,0,MPIW
     .               ,status,ierr) 
      endif
      tcomunica = getime() - tcomunica
c ===  
c     do i = 1, nno 
c       print*,mype,i,(xl(j,i),j=1,ndm)  
c     enddo  
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * MPI_SPLIT_PROP_INT: distribuicao das propriedade entre            *
c * os processos (INT)                                                *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * x    - propriedades globais                                       *  
c * xl   -                                                            *  
c * dist - distribuicao dos propriedade pelos processos               *  
c * MPIW -                                                            *  
c * mype - numero do meu processo(MPI)                                *  
c * npes - numero de processos (MPI)                                  *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * xl       - distribuicao das coordenadas                           *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine mpi_dist_prop_int(x,xl,dist,MPIW,l,npes,mype)
      implicit none
      include 'mpif.h'
      include 'time.fi'
      integer status(MPI_STATUS_SIZE)
      integer npes,mype
      integer l
      integer dist(*)
      integer x(l,*),xl(l,*)
      integer MPIW,nxl,ierr,your_nxl,nxli,i,j,pe
c ===
c     call MPI_Barrier(MPIW,ierr)
      tcomunica = getime() - tcomunica
      nxl = dist(mype+2)-dist(mype+1)
      if(mype .eq. 0) then
        do i = 1, nxl
          do j = 1, l
            xl(j,i) = x(j,i)
          enddo
        enddo
        do pe = 2, npes
          your_nxl = dist(pe+1)-dist(pe)
          nxli = dist(pe) + 1
          call MPI_Send(x(1,nxli),your_nxl*l,MPI_INTEGER
     .                 ,pe-1,0,MPIW,ierr) 
        enddo
      else
        call MPI_Recv(xl,nxl*l,MPI_INTEGER,0,0,MPIW,status,ierr) 
      endif
      tcomunica = getime() - tcomunica
c === 
c     do i = 1, nno
c        print*,mype,i,(xl(j,i),j=1,ndm)  
c     enddo  
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * GLOBAL_DOUBLE:gera o vetor global de real*8                       *
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * g        -                                                        *  
c * l        - vetor distribuido                                      *  
c * dist     - distribuicao dos propriedade nodais pelos processos    *  
c * MPIW     - comunicador mpi                                        *  
c * mype     - numero do meu processo(MPI)                            *  
c * npes     - numero de processos (MPI)                              *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * g        - vetor global                                           *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine global_double(g,l,dist,nl,nc,npes,mype,MPIW)
      implicit none
      include 'mpif.h'
      include 'time.fi'
      integer npes,mype
      integer nl,nc
      integer dist(*)
      real*8 g(*),l(*)
c     integer MPIW,nno,ierr,your_nno,nnoi,i,j,pe,nnoa,nnod
      integer MPIW,nn,ierr,i
c ... !!alocacao dinamica!!     
      integer dspl(npes),rcont(npes)
c ===
      tcomunica = getime() - tcomunica
      nn = dist(mype+2) - dist(mype+1)
      do i = 1, npes
        rcont(i) = (dist(i+1)-dist(i))*nc
      enddo
      dspl(1) = 0
      do i = 1, npes - 1 
        dspl(i+1) = dspl(i) + rcont(i)
      enddo
      call Mpi_Allgatherv(l,nn*nc,MPI_DOUBLE_PRECISION,g,rcont,dspl
     .                   ,MPI_DOUBLE_PRECISION,MPIW,ierr)  
      tcomunica = getime() - tcomunica
c ===
      return
      end
c =======================================================================
c *********************************************************************
c
c *********************************************************************
c * GLOBAL_INT:gera o vetor global de integer                         *
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * g        -                                                        *  
c * l        - vetor distribuido                                      *  
c * dist     - distribuicao dos propriedade nodais pelos processos    *  
c * MPIW     - comunicador mpi                                        *  
c * mype     - numero do meu processo(MPI)                            *  
c * npes     - numero de processos (MPI)                              *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * g        - vetor global                                           *  
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine global_int(g,l,dist,nl,nc,npes,mype,MPIW)
      implicit none
      include 'mpif.h'
      include 'time.fi'
      integer npes,mype
      integer nl,nc
      integer dist(*)
      integer g(*),l(*)
      integer MPIW,nn,ierr,i,j
c ... !!alocacao dinamica!!     
      integer dspl(npes),rcont(npes)
c .....................................................................      
c ===
      tcomunica = getime() - tcomunica
      nn = dist(mype+2) - dist(mype+1)
      do i = 1, npes
        rcont(i) = (dist(i+1)-dist(i))*nc
      enddo
      dspl(1) = 0
      do i = 1, npes - 1 
        dspl(i+1) = dspl(i) + rcont(i)
      enddo
      call Mpi_Allgatherv(l,nn*nc,MPI_INTEGER,g,rcont,dspl,MPI_INTEGER
     .                   ,MPIW,ierr)  
      tcomunica = getime() - tcomunica
c .....................................................................
c     do i = 1,nl  
c        if(mype .eq. 0)  print*,mype,"el",i,(g((i-1)*nc+j),j=1,nc)  
c       if(mype .eq. 1)  print*,mype,"el",i,(g(j,i),j=1,nc)  
c     enddo  
c =======================================================================
c
c ===
c     call finalize()
      return
      end
c =======================================================================
c *********************************************************************
#else
      subroutine global_double(g,l,dist,nl,nc,npes,mype,MPIW)
      return
      end
      subroutine global_int(g,l,dist,nl,nc,npes,mype,MPIW)
      return
      end
#endif
