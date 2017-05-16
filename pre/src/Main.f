c **********************************************************************
      program PREPAR
c *********************************************************************
c *                                                                   *
c *   Pre-processador para programa de elementos finitos              *
c *      paralelizado para memoria distriburanka com MPI              *
c *                                                                   *
c *     baseado na rotina rdat.f do programa Mef                      *
c * ----------------------------------------------------------------- *
c * Variaveis da malha :                                              *
c * nnode   - numero de nos                                           *
c * numel   - numero de elementos                                     *
c * numat   - numero de materias                                      *
c * maxno   - numero de nos nos elementos                             *
c * ndf     - graus de liberdade mecanico                             *
c * ndft    - graus de liberdade termico                              *
c * i_ix    - ponteiro das conectivrankades                           *
c * i_ie    -                                                         *
c * i_e     -                                                         *
c * i_x     - ponterio das coordenadas                                *
c * i_f     -                                                         *
c * i_nload -                                                         *
c * i_eload -                                                         *
c * i_num   -                                                         *
c * i_u     -                                                         *
c * i_v     -                                                         *
c * i_a     -                                                         *
c * i_nloadt-                                                         *
c * i_eloadt-                                                         *
c * i_ut    -                                                         *
c * i_vt    -                                                         *
c * i_w     -                                                         *
c * Variaveis do Mpi :                                                *
c * ierr    - variavel de erro do Mpi                                 *
c * rank    - rank do processo                                        *
c * nprcs   - numero de precesso                                      *
c *********************************************************************
      use Malloc
      implicit none
c .....................................................................
c
c ... Inclusao de headers:
c
      include 'elementos.fi'
      include 'load.fi'
      include 'string.fi'
      include 'parallel.fi'
      include 'time.fi'
      include 'termprop.fi'
c
c ... Variaveis de arquivos
c
      character*80 filein,fileout
      integer nin,nout,naux,nelemtload
      data nin /1/  nout /2/ naux /3/ nelemtload /4/
c
c ... Variaveis descritivas do problema:
c
      integer nnodev,nnode,numel,numat,maxno,maxnov,ndf,ndft,ndm
      integer npi
      integer*8 i_ix,i_nen,i_ie,i_e,i_x
c ... mecanico - poromec
      integer*8 i_f,i_u,i_v,i_a,i_tx0,i_id,i_nload,i_eload,i_eloadp
c ... termico
      integer*8 i_ft,i_ut,i_ut0,i_du,i_vt,i_w,i_idt,i_nloadt,i_eloadt
      integer nlines,pnlines
      character*200 lines(2000),plines(100)
      integer ncont
      parameter (ncont = 100) 
      logical rload(ncont)
      character*80 file_prop(maxmate)
c
c ... Variaveis de particionamento:
c
      integer*8 i_ep,i_np,i_ix0,i_nnof,i_elGL
      integer*8 i_xl,i_lel,i_rank
      integer*8 i_idl,i_nloadl,i_eloadl,i_eloadpl,i_fl
      integer*8 i_idtl,i_nloadtl,i_eloadtl,i_ftl,i_ut0l
      integer*8 i_nodedist,i_elmdist
      integer my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
      integer my_nnode,my_numel_nov,my_numel_ov,my_numel
      integer part
      logical partmesh
c
c ... opcoes no saida para o prompt
c
      logical vrdat
c
c ... opcoes de saida 
c
      logical vtkf,vtkp,mefp,bvtk,contp,vtkgeores
      logical legacy
c .....................................................................
c
c ... Mpi
c
c ... variaveis      
      integer MPIW
      integer rank,nnos,i,j,k,div,loop
c .....................................................................
c
c ... auxilires            
c
      integer dum
c
c ... inicializacao dos ponteiros
      i_ix       = 1 
      i_ie       = 1
      i_e        = 1
      i_x        = 1
c ... mecanico - poromec
      i_id       = 1
      i_nload    = 1
      i_eload    = 1
      i_eloadp   = 1
      i_u        = 1
      i_f        = 1
      i_v        = 1
      i_a        = 1
      i_tx0      = 1
c ... termico
      i_idt      = 1
      i_nloadt   = 1
      i_eloadt   = 1
      i_ft       = 1
      i_ut       = 1
      i_ut0      = 1
      i_du       = 1
      i_vt       = 1
      i_w        = 1 
c
      i_ep       = 1
      i_np       = 1
      i_ix0      = 1
      i_nen      = 1
      i_xl       = 1
      i_lel      = 1
      i_idl      = 1
      i_nloadl   = 1
      i_eloadl   = 1
      i_eloadpl  = 1
      i_fl       = 1
      i_idtl     = 1
      i_nloadtl  = 1
      i_eloadtl  = 1
      i_ftl      = 1
      i_ut0l     = 1
      i_nodedist = 1
      i_elmdist  = 1
c.....................................................................
c
c ...
c     modo verbose para rdata                                  
      vrdat      = .true.
c     escrita do particionamento malha global no formato vtk
      vtkf       = .false.
c     escrita do particionamento no formato vtk
      vtkp       = .false.
c     escrita no formato do mefpar
      mefp       = .false.
c     escrita da geometria com as restricoes malha global vtk
      vtkgeores  = .false.
c     escrita do vtk binario  
      bvtk       = .false.
c     escrita do vtk no formato legacy       
      legacy     = .true.
c.....................................................................
c
c ... gera arquivos apeans com o map
      contp      = .false.
c.....................................................................
c 
c ... tipo geometrico dos elementos
      geom_el_type(1)  = type_barra2  !barra!
      geom_el_type(2)  = type_tria3   !tria3!
      geom_el_type(3)  = type_tria3   !tria3!
      geom_el_type(4)  = type_quad4   !quad4!
      geom_el_type(5)  = type_quad4   !quad4!
      geom_el_type(6)  = type_tetra4  !tetra4!
      geom_el_type(7)  = type_hexa8   !hexa8!
      geom_el_type(10) = type_prism6  !prisma-interface!
      geom_el_type(11) = type_hexa8   !hexa8-interface!
      geom_el_type(12) = type_hexa8   !hexa8-grampo!
      geom_el_type(13) = type_prism6  !prisma-interface!
      geom_el_type(15) = type_tetra10 !tetra10!
      geom_el_type(16) = type_tetra10 !tetra10-vprop!
      geom_el_type(17) = type_hexa20  !hexa20!
      geom_el_type(18) = type_hexa20  !hexa20-vprop!
      geom_el_type(35) = type_tetra10 !tetra10-plastic!
      geom_el_type(36) = type_tetra10 !tetra10-plastic-vprop!
      geom_el_type(37) = type_hexa20  !hexa20-plastic!
      geom_el_type(38) = type_hexa20  !hexa20-plastic-vprop!
c......................................................................
c
c ... variavel de controle de tempo
      trmetis    = 0.d0
      tcomunica  = 0.d0
      tmetis     = 0.d0
      tpart      = 0.d0
      trmesh     = 0.d0
c.....................................................................
c =====================================================================
c
c === iniciando
      call start(MPIW,rank,nnos,ierr)
c =====================================================================
c
c === lendo dados necessario para o metis 
      if(rank.eq.0)then
        call getinfoterm(maxmem   ,filein,fileout,nprcs,ovlp
     1                  ,vtkf     ,vtkp  ,mefp   ,bvtk ,vtkgeores
     2                  ,nin)
      endif
c =====================================================================
c
c === iniciando malloc
      call init_malloc()
c =====================================================================
c 
c === 
      if(rank.eq.0)then
        print*,"Reading file for metis ..."
        trmetis = getime() 
        call read_metis(i_ix0,i_nen,nnodev,nnode,numel
     .                 ,maxnov,maxno,nin,.true.)
        trmetis = getime() - trmetis
        print*,"End of the reading."
      endif
c =====================================================================
c
c === comunicando so no MPI
c     call comm_opt(MPIW,filein,fileout,ovlp,vtkp,mefp,nprcs,ierr)
c =====================================================================
c
c === comunicando dados para o metis no MPI
c     call comm_metis(MPIW,nnode,numel,maxno,ntria3,nquad4,ntetra4
c    .               ,nhexa8,rank)
c =====================================================================
c
c
c === gerando o grafo da malha
      tmetis = getime()
      call call_metis(i_np,i_ep,i_ix0,i_nen,nnode,numel,maxno,rank
     .               ,nnos,nprcs,MPIW)
      tmetis = getime() - tmetis
c =====================================================================
c
c ===  
c ... leitura da malha
      if((rank .eq. 0) )then
        open(nin, file= filein, status= 'old')
        print*, 'Reading mesh ...'
        trmesh  = getime()
        call read_mef(nnodev,nnode,numel,numat,maxnov,maxno
     1                ,ndf   ,ndft ,ndm  ,npi  
     2                ,i_ix  ,i_ie ,i_e  ,i_x
     3                ,i_id  ,i_nload    ,i_eload ,i_eloadp
     4                ,i_f   ,i_u        ,i_tx0   ,i_v    ,i_a
     5                ,i_idt ,i_nloadt   ,i_eloadt
     6                ,i_ut  ,i_ut0      ,i_du    ,i_vt   ,i_w
     7                ,lines ,nlines     ,plines  ,pnlines
     8                ,vrdat      ,rload   ,file_prop,ncont,nin)
        close(nin)
        trmesh  = getime()- trmesh
        print*,"End of the reading."
      endif
c      call printvetorint(ia(i_eloadt),numel,7)
c ####################################################################  
c ######################PARA MALHA ESTRUTURADAS#######################  
c     i_np  = alloc_4('np      ',  1,nnode)  
c     i_ep  = alloc_4('ep      ',  1,numel)  
c     call struct_mesh(ia(i_ix),ia(i_x),ia(i_np),ia(i_ep),nnode,numel
c    .                ,maxnov,maxno,ndm,nprcs)
c     call struct_cubo(ia(i_ix),ia(i_x),ia(i_np),ia(i_ep),nnode,numel
c    .                ,maxnov,maxno,ndm,nprcs)
c ####################################################################  
c .....................................................................
c =====================================================================
c
c ... escrita da malha 
      if((rank .eq. 0) )then
        print*,"Writing mesh ..."
c ... particionamento
        if(vtkf) then      
          print*,"Partition ..."
          call writemesh(ia(i_np),ia(i_ep),ia(i_ix),ia(i_x),nnode
     .                  ,numel,maxno,ndm,fileout,bvtk,legacy,nout)
        endif
c .....................................................................          
c
c ... carregamentos
        if(vtkgeores) then 
          print*,"loading ..."
          call write_mesh_geo(ia(i_ix) ,ia(i_x),ia(i_ie)
     1                 ,ia(i_id) ,ia(i_nload)  ,ia(i_eload),ia(i_eloadp)
     2                 ,ia(i_f)  ,ia(i_u)      ,ia(i_tx0)
     3                 ,ia(i_idt),ia(i_nloadt) ,ia(i_eloadt)
     4                 ,ia(i_ft) ,ia(i_ut0)
     5                 ,nnode    ,numel ,ndf    ,ndft
     6                 ,maxno    ,ndm   ,fileout,bvtk
     7                 ,macros   ,legacy,nout   ,nelemtload)
c .....................................................................
          print*,"End of the writing."
        endif
      endif
c .....................................................................
c ===
c
c === comunicar a malha para o mef dat
      tcomunica = getime()
      call comm_mesh_mef(MPIW,nnode,numel,maxno,numat,ndf,ndft
     1                  ,ndm,npi,lines,nlines,i_ix,i_ie,i_e,i_x
     2                  ,i_id,i_nload,i_eload,i_eloadp
     3                  ,i_f,i_u,i_v,i_a
     4                  ,i_idt,i_nloadt,i_eloadt
     5                  ,i_ft,i_ut0,i_vt,i_w
     6                  ,i_xl,i_lel,i_idl,i_nloadl
     7                  ,i_eloadl,i_eloadpl
     8                  ,i_fl,i_idtl,i_nloadtl,i_eloadtl,i_ftl
     9                  ,i_ut0l
     1                  ,i_nodedist,i_elmdist
     2                  ,rank,nnos,rload,ncont)
       tcomunica = getime() - tcomunica
c =====================================================================
c
c ===
c... Div das part pelos processos com/sem Mpi
      call div_work(MPIW,ierr,nprcs,nnos,div,part,loop,rank)
c     print*,"Main",div,part,loop
c =====================================================================
c
c === particinamento
      j          = part  
      k          = div  
      partmesh   = .false.
c     do i = part, part+div-1
      do i = 1 , loop
        if(k.gt.0)print*,"mype",rank,"partition",j, "generating map ..."
        tpart  = getime() - tpart
        call partdmesh(i_np,i_ep,i_lel
     .                ,i_elGL,i_elLG,i_noGL,i_noLG
     .                ,i_fmap,i_dspl,i_rcvs,i_nnof,i_ranks,i_sizes  
     .                ,i_elmdist,numel,nnode,maxno
     .                ,my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
     .                ,my_nnode,my_numel,my_numel_nov,my_numel_ov
     .                ,nprcs,j,rank,nnos,MPIW,ovlp,partmesh,k)
        partmesh = .true.
        tpart  = getime() - tpart
        if(k.gt.0)print*,"mype",rank,"partition",j, "map generated"
c       call use_mem("MB")
c =====================================================================
c 
c ... ecreve a malha particonada
c
c ... apenas o mapa paralelo
        if(contp .and. k .gt. 0)then
          call controlepre(ia(i_noLG),ia(i_elLG),ia(i_ranks),ia(i_sizes)
     .                    ,ia(i_rcvs)
     .                    ,ia(i_dspl),ia(i_fmap),ia(i_nnof)
     .                    ,my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
     .                    ,my_nnode,my_numel_nov,my_numel_ov
     .                    ,nnode,numel,ovlp,j,nprcs,fileout,nout)
        endif
c ... vtk
        if(vtkp) then
          call writemeshpart(ia(i_np),ia(i_ep),ia(i_lel),ia(i_xl)
     .                      ,ia(i_ie)                        
     .                      ,ia(i_elLG),ia(i_noLG),ia(i_noGL)
     .                      ,ia(i_nodedist),ia(i_elmdist),nnode,numel
     .                      ,my_nnode,my_numel,my_numel_nov,my_numel_ov
     .                      ,maxno,ndm,fileout,bvtk
     .                      ,legacy,nout,j,k,rank,nnos,MPIW)
        endif 
c ... dat 
        if(mefp)then
          call writemydatpart(ia(i_lel)  ,ia(i_xl)     ,ia(i_elLG)
     1         ,ia(i_noLG) ,ia(i_noGL)   ,ia(i_ranks)
     2         ,ia(i_sizes),ia(i_fmap)   ,ia(i_nnof)
     3         ,ia(i_rcvs) ,ia(i_dspl)   ,ia(i_ie),ia(i_e)
     4         ,ia(i_idl)  ,ia(i_nloadl) ,ia(i_eloadl),ia(i_eloadpl)
     5         ,ia(i_fl)   ,ia(i_u)      ,ia(i_tx0)
     6         ,ia(i_idtl) ,ia(i_nloadtl),ia(i_eloadtl)
     7         ,ia(i_ftl)  ,ia(i_ut0l)
     8         ,ia(i_nodedist),ia(i_elmdist)
     9         ,nnodev,nnode,numel,maxnov,maxno
     1         ,ndm,numat,ndf,ndft,npi ,file_prop
     2         ,my_nno1,my_nno2,my_nno3,my_nno4
     3         ,my_nno1a,my_nnode
     4         ,my_numel,my_numel_nov,my_numel_ov
     5         ,lines,nlines,plines,pnlines
     6         ,macros,fileout,nprcs,ovlp
     7         ,j,k,rank,nnos,MPIW,nout)
        endif
c ...       
        k   =   k - 1
        j   =   1 + j
c .....................................................................        
      enddo
c     print*,'exit',rank
c =====================================================================
c
c ===
 200  continue
      call mytemp(MPIW,fileout,nnos,rank,nprcs,naux)
c =====================================================================
c
c ===
      call finalize()
      end
c *********************************************************************
c 
c *********************************************************************
c
c *********************************************************************
c * USE_MEM : estimando a memoria usada no programa.                  *
c *                                                                   *
c * Parametro de entrada :                                            *
c * ----------------------                                            *
c *                                                                   *
c * Parametros de saranka:                                            *
c * --------------------                                              *
c * use_mem - retorna o numero de bytes usasos do vetor ia            *
c *                                                                   *
c *********************************************************************
      subroutine use_mem(memo)
      use Malloc
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
      character memo*2
      character*8 arname(maxnpts)
      integer*8 i,nalp,ip(maxnpts)
      common /malloc_info/ arname,ip,nalp
      i=ip(nalp+1)
      if(memo .eq. ' B') then  
        print*,'Memoria usada B',i*4.0
      else if(memo .eq. 'KB') then  
        print*,'Memoria usada KB',i*4.0/1024
      else if(memo .eq. 'MB') then  
        print*,'Memoria usada MB',i*4.0/1024**2
      else if(memo .eq. 'GB') then  
        print*,'Memoria usada GB',i*4.0/1024**3
      endif  
      return
      end
c *********************************************************************
 
c *********************************************************************
c * TES_MEM : estimando a memoria usada no programa.                  *
c *                                                                   *
c * Parametro de entrada :                                            *
c * ----------------------                                            *
c *                                                                   *
c * Parametros de saranka:                                            *
c * --------------------                                              *
c * use_mem - retorna o numero de bytes usasos do vetor ia            *
c *********************************************************************
      subroutine tes_mem(i_i,iname)
      implicit none
#ifdef MPI     
      include 'mpif.h'
      integer ierr
#endif      
      integer*8 i_i
      character*9 iname 
      if(i_i.lt.0)then
        write(*,'(a,a,i16)')"Erro na alocacao do vetor: ",iname,i_i
#ifdef MPI     
        call MPI_Finalize(ierr)   
#endif      
        stop
      endif
      return
      end
c *********************************************************************
