c*****************************Svn***************************************
c*$Date: 2015-06-09 16:23:50 -0300 (Tue, 09 Jun 2015) $                 
c*$Rev: 971 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************  
      program Mef
c **********************************************************************
c *                                                                    *
c *   Metodo dos elementos finitos para problemas termo-mecanicos      *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'string.fi'
      include 'transiente.fi'
      include 'parallel.fi'
      include 'elementos.fi'
      include 'time.fi'
      include 'openmp.fi'
c ----------------------------------------------------------------------
c
c ... Variaveis da estrutura interna de macro-comandos:
c
      character*8 mc,macro(40),lmacro(50)
      character*30 string
      character*80 str
      integer  iloop,imacro,nmacro,nmc,loop,j
c ----------------------------------------------------------------------
c
c ... Variaveis para controle de arquivos:
c
      character*80 prename,fname,name,filein
      character*80 pnodename
      integer nin,nplot,ngid,nout,logsolv,fconf
      integer totfiles,openflag
      integer*8 i_no,i_nfile
      integer num_pnode
c ... arquivo de impressao nos nos ( temp,flux,disp,stress,...)  
      integer nfiles
      parameter ( nfiles = 4)
      logical new_file(nfiles),flag_pnd
c      logical cont1
c ----------------------------------------------------------------------
c
c ... Variaveis de controle de solucao:
c
      integer maxit,maxnlit,tmaxnlit,ngram,stge,solver,istep
      integer ilib,ntn,code
      real*8  tol,solvtol,resid,resid0
      logical reordf,unsym
c ----------------------------------------------------------------------
c
c ... Variaveis descritivas do problema:
c
      integer nnodev,nnode,numel,numat,nen,nenv,ndf,ndm,nst
      integer neq,nequ,neqp,nad,nadpu
      logical dualCsr
c ----------------------------------------------------------------------
c
c ... Variaveis da interface de linha de comando
      integer nargs
      character arg*80
c ----------------------------------------------------------------------
c
c ... Variaveis locais:
c
      integer i,k
      real*8  dot_par
c ----------------------------------------------------------------------
c
c ... Variaveis de medicao de tempo:
c
      real*8 timei      
c ----------------------------------------------------------------------
c
c ... Variaveis Colormesh:
c
      integer numcolors 
c ... Ponteiros:
      integer*8 i_colorg,i_elcolor    
c ----------------------------------------------------------------------
c
c ... Ponteiros:
c
c ... malha
      integer*8 i_ix,i_id,i_ie,i_nload,i_eload,i_e,i_x,i_xq,i_inum
c ... arranjos locais ao elemento
      integer*8 i_xl,i_ul,i_vl,i_zl,i_wl,i_pl,i_sl,i_ld
      integer*8 i_f,i_f0
      integer*8 i_u,i_u0
c ... sistema de equacoes
      integer*8 i_ia,i_ja,i_ad,i_au,i_al,i_m,i_b,i_b0
c ... arranjos globais (MPI - escrita)
      integer*8 i_g,i_g1,i_g2
c ----------------------------------------------------------------------
c
c ... Variaveis de controle do MPI:
c
      integer status(MPI_STATUS_SIZE)      
c ----------------------------------------------------------------------
c
c ... Macro-comandos disponiveis:
c
      data nmc /40/
      data macro/'loop    ','        ','mesh    ','solv    ','dt      ',
     .'pgeo    ','pgeoquad','        ','        ','        ','        ',
     .'        ','        ','        ','        ','pdisp   ','pstress ',
     .'        ','        ','        ','        ','        ','        ',
     .'        ','        ','maxnlit ','tmaxnlit','nltol   ','        ',
     .'        ','        ','setpnode','        ','        ','pndisp  ',
     .'pnstress','        ','maxit   ','solvtol ','stop    '/
c ----------------------------------------------------------------------
c
c ... Arquivos de entrada e saida:
c
      data nin /1/, nplot /3/, ngid /4/, logsolv /10/, nout /15/
      data fconf /5/
      data flag_pnd /.false./ 
c     arquivo de impressao de nos associados aos seguintes inteiros
c     nfile = 50,51,52,...,num_pnode
c ----------------------------------------------------------------------      
c
c ... Inicializacao MPI:
c
      call MPI_INIT( ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprcs, ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
c ----------------------------------------------------------------------
c
c ... Inicializacao de variaveis da estrutura interna de macro-comandos:
c
      iloop   = 0
      imacro  = 0
      nmacro  = 0
c ----------------------------------------------------------------------
c
c ... Inicializacao de variaveis de controle de solucao:
c
c.......................................................................
      istep   =  0
      t       =  0.d0
      dt      =  1.d0
      alfa    =  1.d0
      beta    =  1.d0
c ... reordf  =  true -> reordenacao Cuthill-McKee
      reordf  = .false.
c ... maxit   =  numero max. de iteracoes do solver iterativo
c ... solvtol =  tolerancia para o solver iterativo
c ... maxnlit =  numero max. de iteracoes nao-lineares
c ... tol     =  tolerancia do algoritmo nao-linear
c ... ngram   =  base de Krylov (utilizada somente no gmres)
      maxit   =  50000
      solvtol =  1.d-10
      maxnlit =  20
      tmaxnlit=  20
      tol     =  1.d-04
      ngram   =  15
c ... unsym   =  true -> matriz de coeficientes nao-simetrica      
c ... solver  =  1 (pcg), 2 (gmres), 3 (gauss / LDU), 4 (bicgstab)
c ... stge    =  1 (csr), 2 (edges), 3 (ebe), 4 (skyline)
      unsym   = .false.
      solver  =  1
      stge    =  1
      dualCsr = .true.
      resid0  =  0.d0
c ... ilib    =  1 define a biblioteca padrão ( default = termico )
      ilib    =  1
c     OpenMP 
      openmp  = .false.
c .....................................................................
c
c ... Inicializacao da estrutura de dados de gerenciamento de memoria:
c
      call init_malloc(maxmem)
c ----------------------------------------------------------------------
c ...
      call init_openmp(num_threads,openmp,my_id)
c ......................................................................
c
c ... Abertura de arquivos:    
      nargs = iargc()
   10 continue
      if (my_id .eq. 0) then
c ... intervace de linha de comando        
        if(nargs .gt. 0) then
          call getarg(1,arg)
          filein = arg
        else
          print*, 'Arquivo de dados:'
          read(*,'(a)') filein
        endif 
      endif      
      if (nprcs .eq. 1) then
         open(nin, file= filein, status= 'old', err=15, action= 'read')
         goto 20
   15    continue
         print*, 'Arquivo nao existente !'
         nargs = 0
         goto 10   
      else
c ...    Passa nome do arquivo de entrada do processo 0 para os demais:
         call MPI_BCAST(filein,20,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c ...    Nome do arquivo de dados do processo my_id:
         fname = name(filein,my_id,13)
c ...    Abre o arquivo de dados de cada processo:
         totfiles = 0
         openflag = 0
         open(nin, file= fname, status= 'old', err=11, action= 'read')
         openflag = 1
   11    continue
c ...    Testa se todos abriram sem erro:
         call MPI_ALLREDUCE(openflag,totfiles,1,MPI_INTEGER,MPI_SUM,
     .                      MPI_COMM_WORLD,ierr)
         if (totfiles .ne. nprcs) then
            if (my_id .eq. 0) print*, '*** Erro na abertura de arquivos'
            call stop_mef()
         endif
      endif
   20 continue
c ......................................................................   
      if (my_id .eq. 0) then   
c ... intervace de linha de comando        
        if(nargs .eq. 2) then
          call getarg(2,arg)
          prename = arg
        else
          print*, 'Arquivo de saida: '
          read(*,'(a)') prename
        endif 
         fname = name(prename,nprcs,15)
         open(logsolv,file=fname)
         write(logsolv,'(a)') 'Solver control flop.'
c ...    Arquivo de resultados para o GID:
c         fname = name(prename,istep,2)
c         open(ngid,file=fname)
c         write(ngid,'(a)') 'GID post results file 1.0'
      endif
c ......................................................................      
      call MPI_barrier(MPI_COMM_WORLD,ierr)
c-----------------------------------------------------------------------
c
c ... Leitura dos macro-comandos:
c
c ......................................................................
   50 continue
      if (iloop .eq. 0) then
         call readmacro(nin,.true.)
         write(mc,'(8a)') (word(i),i=1,8)
      else
         if (imacro .eq. 0 .or. imacro .eq. nmacro) then
            imacro = 1
         else
            imacro = imacro + 1
         endif
         mc  = lmacro(imacro)
         iloop = iloop - 1
      endif
      do 60 j = 1, nmc
         if (mc .eq. macro(j)) goto 70
   60 continue
      goto 50
   70 continue
      goto (100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,
     .     1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,
     .     2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,
     ,     3800,3900,5000) j
c ----------------------------------------------------------------------
c
c ... Execucao dos macro-comandos:
c
c ----------------------------------------------------------------------
c
c ... Macro-comando LOOP:
c ......................................................................
  100 continue
      call readmacro(nin,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 120,end = 120) loop
      nmacro = 0
      imacro = 0
      iloop  = 0
  110 continue
      call readmacro(nin,.true.)
      write(mc,'(8a)') (word(i),i=1,8)        
      if (mc .eq. 'next ') goto 50
      nmacro = nmacro + 1
      iloop = loop*nmacro
      lmacro(nmacro) = mc
      goto 110
  120 continue
      print*,'Erro na leitura da macro (LOOP) !'
      goto 5000            
c ----------------------------------------------------------------------
c
c ... Macro-comando NCONC: 
c
c ......................................................................
  200 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando MESH:
c
c ......................................................................
  300 continue
      if(my_id.eq.0)print*, 'Macro MESH'  
c
c.... Leitura de dados:
c
c      call rdat(nnode,numel,numat,nen,ndf,ndm,nst,i_ix,i_id,i_ie,
c     .          i_nload,i_eload,i_inum,i_e,i_x,i_f,i_u,i_v,i_a,nin)
c
      call rdat(nnode ,nnodev ,numel  ,numat  
     .         ,nen   ,nenv ,ndf    ,ndm    ,nst
     .         ,i_ix  ,i_ie   ,i_inum ,i_e    ,i_x
     .         ,i_id  ,i_nload,i_eload,i_f    ,i_f0
     .         ,i_u   ,i_u0   ,nin ) 
c
c     -----------------------------------------------------------------
c     | ix | id | ie | nload | eload | inum | e | x | f | f0 | u | u0 |
c     -----------------------------------------------------------------
c
c ......................................................................      
c
c.... Controle de tempos:
c
      soltime    = 0.d0
      elmtime    = 0.d0
      vectime    = 0.d0
      ovhtime    = 0.d0
      dottime    = 0.d0
      gvtime     = 0.d0
      allgtime   = 0.d0
      allrtime   = 0.d0
      sendtime   = 0.d0
      matvectime = 0.d0
      colortime  = 0.d0
      pmatrixtime= 0.d0
      writetime  = 0.d0
      totaltime  = MPI_Wtime()
c ......................................................................      
c
c.... Otimizacao da largura de banda:
c
      timei = MPI_Wtime()
      if (nprcs .gt. 1) then
       call reord(ia(i_ix),ia(i_inum),nno1-nno1a,nnode,numel,nen,reordf)
      else
       call reord(ia(i_ix),ia(i_inum),nnode,nnode,numel,nen,reordf)
      endif
      reordtime = MPI_Wtime()-timei
c ......................................................................
c
c.... Numeracao nodal das equacoes:
c
      timei = MPI_Wtime()
      if(ndf .gt. 0) then
c ... primero numera os deslocamento e depois as pressoes
        if(dualCsr) then
          call numeqpmec1(ia(i_id),ia(i_inum),ia(i_id),
     .                    nnode,nnodev,ndf,neq,nequ,neqp)
        else
c ... primero numera os deslocamento e pressoes ao mesmo tempo
          call numeqpmec2(ia(i_id),ia(i_inum),ia(i_id),
     .                    nnode,nnodev,ndf,neq)
          nequ = 0
          neqp = 0
        endif
      endif
      numeqtime = MPI_Wtime()-timei
c ......................................................................
c
c ... Front Mpi
      timei = MPI_Wtime()
      call init_front(i_noLG,i_noGL,nno1,nno2,nno3,nnofi
     .                ,nno_pload
     .                ,nnoG,nelG,nnode,numel,ovlp,novlp,nprcs,nviz1
     .                ,nviz2,i_rreqs,i_sreqs,i_rcvs0i,i_dspl0i
     .                ,i_fmap0i,mpi)
c
c.... Mapa de equacoes de fronteira:
c
      if (ndf  .gt. 0) call frontb(ndf,ia(i_id),neq,neq1
     .                ,neq2,neq3,neq4,neq1a,neqf1,neqf2,neq32
     .                ,neq_dot
     .                ,i_fmap,i_rcvs,i_dspl,i_xf
     .                ,'fmap     ','rcvs    ','dspl    ','xf      ',0)
      frontime = MPI_Wtime()-timei
c -------------------------------------------------------------------
c | noLG | noGL | elLG | fmap | rcvs | dspl | fmapt | rcvst | dsplt |
c -------------------------------------------------------------------
c
c                     ----------------------------                    
c                     | fmap0i | rcvs0i | dspl0i |
c                     ----------------------------                    
c ......................................................................
c
c.... Arranjos locais de elemento:
c
      i_xl = alloc_8('xl      ',ndm,nenv)
      i_ul = alloc_8('ul      ',1  ,nst)
      i_pl = alloc_8('pl      ',1  ,nst)
      i_sl = alloc_8('sl      ',nst,nst)
      i_ld = alloc_4('ld      ',  1,nst)
c
c     -----------------------------------------
c     | xl | ul | pl | sl | ld |
c     -----------------------------------------
c
c ......................................................................
c
c ... Memoria para a estrutura de dados do sistema de equacoes:
c
      timei = MPI_Wtime()
      if (ndf .gt. 0) then
         call datastruct(ia(i_ix),ia(i_id),ia(i_inum),nnode,nnodev,
     .   numel,nen,ndf,nst,neq,nequ,neqp,stge,unsym,nad,nadpu,
     .   i_ia,i_ja,i_au,i_al,i_ad,
     .   'ia      ','ja      ','au      ','al      ','ad      ',
     .   ovlp,dualCsr)
      endif
      dstime = MPI_Wtime()-timei
c
c ... colorir a malha (openmp)
c
c     colortime = MPI_Wtime()
c     call coloredmesh(ia(i_ix),nnode,numel,nen,numcolors,i_colorg,
c    .                 i_elcolor)           
c     colortime = MPI_Wtime()-colortime
c ......................................................................
c
c ......................................................................
c
c ... Memoria para o vetor de forcas e solucao:
c
c ......................................................................
      if (ndf .gt. 0) then
         i_b0  = alloc_8('b0      ',    1,neq+neq3+neq4)
         i_b   = alloc_8('b       ',    1,neq+neq3+neq4)
         call azero(ia(i_b0),neq+neq3+neq4)      
         call azero(ia(i_b ),neq+neq3+neq4)
c ...    Memoria para o precondicionador diagonal:
         i_m   = alloc_8('m       ',    1,neq+1)
      endif
c ......................................................................
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando SOLV:
c
c ......................................................................
  400 continue
      if(my_id.eq.0)print*, 'Macro SOLV '
         ilib = 1
         i    = 0
c ... forcas de volume(f)                         
         call pformpmec(ia(i_ix),ia(i_eload),ia(i_ie),ia(i_e),
     .               ia(i_x),ia(i_id),ia(i_ia),ia(i_ja),
     .               ia(i_au),ia(i_al),ia(i_ad),ia(i_b0),ia(i_u),
     .               ia(i_xl),ia(i_ul),ia(i_pl),ia(i_sl),ia(i_ld),
     .               numel,nen,nenv,ndf,
     .               ndm,nst,neq,nequ,nad,nadpu,.false.,.true.,unsym,
     .               stge,3,ilib,i,dualCsr)
c .....................................................................
c
c ... calculo da matriz de coeficientes ( Ku = f )
        call pformpmec(ia(i_ix),ia(i_eload),ia(i_ie),ia(i_e),
     .               ia(i_x),ia(i_id),ia(i_ia),ia(i_ja),
     .               ia(i_au),ia(i_al),ia(i_ad),ia(i_b0),ia(i_u),
     .               ia(i_xl),ia(i_ul),ia(i_pl),ia(i_sl),ia(i_ld),
     .               numel,nen,nenv,ndf,
     .               ndm,nst,neq,nequ,nad,nadpu,.true.,.true.,unsym,
     .               stge,2,ilib,i,dualCsr)
c .....................................................................
        call testeMatvec(ia(i_ia),ia(i_ja),ia(i_ad),ia(i_al)
     .                  ,ia(i_m) ,ia(i_b0)
     .                  ,neq,nequ,nad,nadpu)
      goto 50 
c ----------------------------------------------------------------------
c
c ... Macro-comando: DT
c
c ......................................................................
  500 continue
      if(my_id.eq.0)print*, 'Macro DT   '
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err = 510,end = 510) dt
      goto 50
  510 continue
      print*,'Erro na leitura da macro (DT) !'
      goto 5000      
c ----------------------------------------------------------------------
c
c ... Macro-comando: PGEO
c
c ......................................................................
  600 continue
      call global_ix(nenv+1,numel_nov,i_ix,i_g,'ixg     ')
      call global_v(ndm,nno_pload,i_x,i_g2,'xg      ')
      if (my_id .eq. 0) then
         print*, 'Macro PGEO'
c ...    Geometria: 
         writetime = writetime + MPI_Wtime()-timei
         call writeMeshGeo(ia(i_g)   ,ia(i_g2),nnodev ,numel
     .                       ,nenv ,ndm     ,prename,.false.
     .                       ,.true. ,nplot)
         writetime = writetime + MPI_Wtime()-timei
c ......................................................................         
      endif
      if (nprcs .gt. 1) then
         i_g2 = dealloc('xg      ')
         i_g  = dealloc('ixg     ')
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierr)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PGEOQ
c
c ......................................................................
  700 continue
      if(my_id.eq.0)  then  
        print*, 'Macro PGEOQ'
        i_xq = alloc_8('xq      ',ndm,nnode)
        call mkCoorQuad(ia(i_x) ,ia(i_xq)
     .                 ,ia(i_ix)
     .                 ,numel   ,nen  
     .                 ,nnode   ,nnodev  ,ndm)   
        nhexa20(1:4) = nhexa8(1:4)   
        nhexa8(1:4)  = 0 
        call writeMeshGeo(ia(i_ix)   ,ia(i_xq),nnode   ,numel
     .                   ,nen    ,ndm     ,prename,.false.
     .                   ,.true. ,nplot)
        nhexa8(1:4) = nhexa20(1:4)   
        i_xq = dealloc('xq      ')      
      endif
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PFLUX
c
c ......................................................................
  800 continue
c      soma = 0.d0  
      if(my_id.eq.0)print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
  900 continue
      if (my_id .eq. 0)   print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 1000 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PDEF
c
c ......................................................................
 1100 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SCALE
c
c ......................................................................
 1200 continue
      if(my_id.eq.0) print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
 1300 continue
      if(my_id.eq.0) print*, 'Macro '
      goto 50     
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
 1400 continue
      if(my_id.eq.0)print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
 1500 continue
      if(my_id.eq.0)print*,'Macro SOLVM'
      goto 50 
c ----------------------------------------------------------------------
c
c ... Macro-comando: PDISP
c
c ......................................................................
 1600 continue
      if(my_id.eq.0)print*,'Macro '
      goto 50     
c ----------------------------------------------------------------------
c
c ... Macro-comando: PSTRESS
c
c ......................................................................
 1700 continue
      if(my_id.eq.0)print*, 'Macro Macro'
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 1800 continue
      if (my_id .eq. 0 )   print*, 'Macro  '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 1900 continue
      if(my_id .eq. 0)   print*, 'Macro '
      goto 50     
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2000 continue
      print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2100 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2200 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2300 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2400 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2500 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2600 continue
      if(my_id.eq.0)print*, 'Macro MAXNLIT    '
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2610,end =2610) maxnlit
      write(*,'(a,i10)')' Set max mechanic nonlinear it for ',maxnlit
      goto 50
 2610 continue
      print*,'Erro na leitura da macro (MAXNLIT) !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2700 continue
      if(my_id.eq.0)print*, 'Macro TMAXNLIT    '
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2710,end =2710) tmaxnlit
      write(*,'(a,i10)')' Set max Thermical-Chemical nonlinear it for '
     .                  ,tmaxnlit
      goto 50
 2710 continue
      print*,'Erro na leitura da macro (TMAXNLIT) !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2800 continue
      if(my_id.eq.0)print*, 'Macro NLTOL'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2810,end =2810) tol
      write(*,'(a,d10.2)')' Set noliner tol for ',tol
      goto 50
 2810 continue
      print*,'Erro na leitura da macro (NLTOL) !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2900 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 3000 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 3100 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SETPNODE impressao de grandezas por no no tempo
c
c ......................................................................
 3200 continue
      if(my_id.eq.0) print*, 'Macro SETPNODE'
      if(my_id.eq.0) then
        call readmacro(nin,.false.)
        write(str,'(80a)') (word(i),i=1,80)
        read(str,*,err=2310,end = 2310) pnodename
        goto 2320
c ... problema no arquivo auxiliar        
 2310   continue
        print*,'Erro na leitura da macro (SETPNODE)'
        flag_pnd = .false.
        goto 2330
c ... leitura normal 
 2320   continue     
        call readpnode(pnodename,i_no,i_nfile,num_pnode,flag_pnd,nout)
        new_file(1:nfiles) = .true.
 2330   continue
      endif
      call MPI_BCAST(flag_pnd,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c ... erro na letura do nome do arquivo auxiliar      
      if( flag_pnd .eqv. .false.) call stop_mef()
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:                                                    
c ......................................................................
 3300 continue
      if(my_id.eq.0) print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PNTEMP impressao do fluxo de Calor por no no tempo
c     (SETPNODE)                                                   
c ......................................................................
 3400 continue
      if(my_id.eq.0) print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PNDISP impressao do deslocamento por no no tempo
c     (SETPNODE)                                                   
c ......................................................................
 3500 continue
      if(my_id.eq.0) print*, 'Macro PNDISP    '
      if(flag_pnd.eqv..false.) then
        if(my_id.eq.0)print*,'Nemhum no de impressao para PNDISP!'
        call stop_mef()
      endif
c ... codigo para o arquivo _disp.txt      
      code = 29
c .....................................................................
      call global_v(ndf,nno_pload,i_u,i_g1,'dispG   ')
      string = 'Deselocamentos'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          call printnode(ia(i_g1),ia(i_no+j-1),ndf,istep,dt
     .                 ,string,prename,ia(i_nfile+j-1),code,new_file(3))
        enddo
        new_file(3) = .false.
      endif
      if( nprcs .gt. 1) then
        i_g1      = dealloc('dispG   ')
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierr)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:                                                    
c ......................................................................
 3600 continue
      if(my_id.eq.0) print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 3700 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: Set maxit solver
c
c ......................................................................
 3800 continue
      if(my_id.eq.0)print*, 'Macro MAXIT '
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =3910,end =3910) maxit
      write(*,'(a,i10)')' Set max it solver for ',maxit
      goto 50
 3810 continue
      print*,'Erro na leitura da macro (MAXIT) !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando: Set Solver Tol
c
c ......................................................................
 3900 continue
      if(my_id.eq.0)print*, 'Macro SOLVTOL'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =3910,end =3910) solvtol
      write(*,'(a,d10.2)')' Set solver tol for ',solvtol
      goto 50
 3910 continue
      print*,'Erro na leitura da macro (SOLVTOL) !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando STOP:
c
c ......................................................................
 5000 continue
      if(my_id.eq.0)print*, 'Macro STOP'
c ... fecha o arquivo de entrada de dados
      close(nin)
c ... fecha o arquivo gid , view3d  e log do solv
c     close(ngid)
c     close(nplot)
      close(logsolv)
c .....................................................................
c
c ...
      totaltime = MPI_Wtime() - totaltime
c
c ... arquivo de tempo      
      call write_log_file(nnode ,numel      ,numel_nov,numel_ov
     .                   ,ndf   ,neq        ,neq1     ,neq2    
     .                   ,neq32 ,neq4       ,neq1a    ,neqf1
     .                   ,neqf2 ,nad        ,nad1
     .                   ,openmp,num_threads,numcolors,prename
     .                   ,my_id ,nprcs      ,nout)
c .....................................................................
c
c ... desalocando a memoria do vetor de trabalho  
      call common_finalize()
c .....................................................................
c
c ...
      call stop_mef()
c .....................................................................
      end
c **********************************************************************      
      subroutine testeMatvec(ia ,ja,ad,al
     .                      ,x  ,b
     .                      ,neq,nequ,nad,nadpu)
      implicit none
      integer ia(*),ja(*),dum
      real*8 ad(*),al(*),x(*),b(*)
      integer neq,nequ,nad,nadpu,i
      integer*8 idum
      call matvec_csrcb(neq,nequ,ia,ja,ia(neq+2),ja(nad+1)
     .                 ,ad ,al  ,al(nad+1),b,x 
     .                 ,dum,dum,idum,idum,idum,idum)
c ...
      do i = 1, neq
        print*,i,x(i),b(i)
      enddo
c .....................................................................
      return
      end


