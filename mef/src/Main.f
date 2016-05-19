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
      include 'gravity.fi'
      include 'elementos.fi'
      include 'time.fi'
      include 'openmp.fi'
c ......................................................................
c
c ... Variaveis da estrutura interna de macro-comandos:
c
      character*8 mc,macro(40),lmacro(50)
      character*30 string
      character*80 str
      integer  iloop,imacro,nmacro,nmc,loop,j
      logical flag_macro_mesh
c ......................................................................
c
c ... Variaveis para controle de arquivos:
c
      character*80 prename,fname,name,filein
      character*80 pnodename
      integer nin,nplot,nout,nout_face
      integer logsolv,fconf,logsolvd
      integer totfiles,openflag
      integer*8 i_no,i_nfile
      integer num_pnode
c ... arquivo de impressao nos nos ( pu,stress,stressE,stressB,flux,...)  
      integer nfiles,ifiles
      parameter ( nfiles = 5)
      logical new_file(nfiles),flag_pnd
c      logical cont1
c ......................................................................
c
c ... Variaveis de controle de solucao:
c
      integer maxit,maxnlit,tmaxnlit,ngram,stge,solver,istep
      integer ilib,ntn,code
      real*8  tol,solvtol,resid,resid0
      logical reordf,unsym
c ... pcg duplo
      integer cmaxit
      real*8  ctol,alfap,alfau
c.......................................................................
c
c ... Variaveis descritivas do problema:
c
      integer nnodev,nnode,numel,numat,nen,nenv,ndf,ndm,nst
      logical fporomec,fmec
c ......................................................................
c
c ... Variaveis do sistema de equacoes:
      integer neq,nequ,neqp,nad,naduu,nadpp,nadpu
      integer n_blocks_pu
      logical block_pu
      character*8 sia,sja,sau,sal,sad
c .....................................................................
c .......................................................................
c
c ... precondicionador
      integer precond  
c .......................................................................
c
c ... Variaveis da interface de linha de comando
      integer nargs
      character arg*80
c ......................................................................
c
c ... Variaveis locais:
c
      integer i,k
      real*8  dot_par
c ......................................................................
c 
c ...
      logical bvtk
c ......................................................................
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
c ......................................................................
c
c ...
      real*8 smachn
c ......................................................................
c
c ... Ponteiros:
c
c ... malha
      integer*8 i_ix,i_id,i_ie,i_nload,i_eload,i_e,i_x,i_xq,i_inum
      integer*8 i_ic
c ... arranjos locais ao elemento
      integer*8 i_xl,i_ul,i_pl,i_sl,i_ld,i_dpl,i_txl,i_txnl
c ... forcas e graus de liberdade 
      integer*8 i_f
      integer*8 i_u,i_u0,i_tx0,i_dp
      integer*8 i_tx,i_txb,i_txe,i_flux
c ... sistema de equacoes
      integer*8 i_ia,i_ja,i_ad,i_au,i_al,i_b,i_b0,i_x0,i_bst0
c ... precondicionador
      integer*8 i_m
c ... arranjos globais (MPI - escrita)
      integer*8 i_g,i_g1,i_g2
c ... coo
      integer*8 i_lin,i_col,i_acoo
      integer*8 i_linuu,i_coluu,i_acoouu
      integer*8 i_linpp,i_colpp,i_acoopp
c ......................................................................
c
c ... tensoes iniciais
      logical fstress0,fcstress0
c ......................................................................
c
c ... Variaveis de controle do MPI:
c
      integer status(MPI_STATUS_SIZE)      
c ......................................................................
c
c ... Macro-comandos disponiveis:
c
      data nmc /40/
      data macro/'loop    ','hextotet','mesh    ','solv    ','dt      ',
     .'pgeo    ','pgeoquad','block_pu','gravity ','pardiso ','gmres   ',
     .'deltatc ','pcoo    ','bcgs    ','pcg     ','pres    ','spcgm   ',
     .'solvm   ','pmecres ','bcgsl2  ','minres  ','        ','        ',
     .'        ','        ','maxnlit ','        ','nltol   ','        ',
     .'        ','        ','setpnode','        ','        ','pnup    ',
     .'pnsf    ','config  ','maxit   ','solvtol ','stop    '/
c ......................................................................
c
c ... Arquivos de entrada e saida:
c
      data nin /1/, nplot /3/, logsolv /10/, nout /15/,
     .     logsolvd /16/, nout_face /17/ 
      data fconf /5/
      data flag_pnd /.false./ 
c     arquivo de impressao de nos associados aos seguintes inteiros
c     nfile = 50,51,52,...,num_pnode
c ......................................................................      
c
c ... Inicializacao MPI:
c
      call MPI_INIT( ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprcs, ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
c ......................................................................
c
c ... Inicializacao de variaveis da estrutura interna de macro-comandos:
c
      iloop   = 0
      imacro  = 0
      nmacro  = 0
c ......................................................................
c
c ... Inicializacao de variaveis de controle de solucao:
c
c ......................................................................
      istep   =  0
      t       =  0.d0
      dt      =  1.d0
      alfa    =  1.d0
      beta    =  1.d0
c ... tipo do problema
c ... fporomec  = problema poromecanico                    
c ... fmec      = problema mecanico              
      fporomec = .false.
      fmec     = .false.
c ... tensoes iniciais
c ... fstress0  = tensoes iniciais como condicao inicial
c ... fcstress0 = tensoes iniciais utilizadas       
      fstress0  = .false.
      fcstress0 = .false.
c ... reordf  =  true -> reordenacao Cuthill-McKee
      reordf  = .true.
c ... maxit   =  numero max. de iteracoes do solver iterativo
c ... solvtol =  tolerancia para o solver iterativo
c ... maxnlit =  numero max. de iteracoes nao-lineares
c ... tol     =  tolerancia do algoritmo nao-linear
c ... ngram   =  base de Krylov (utilizada somente no gmres)
c ... precond =  1 - NONE , 2 - diag, 3 - iLDLt(0), 4 - iC(0)
c                5 - diagm
      maxit   =  50000
      solvtol =  1.d-11
      maxnlit =  2 
      tol     =  1.d-04
      ngram   =  50
      precond =  2
c ... cmaxit  =  numero max. de iteracoes do ciclo externo do pcg duplo
c ... ctol    =  tolerancia do ciclo externo do pcg duplo
      cmaxit  =  200
      ctol    =  1.d-6
c ... unsym   =  true -> matriz de coeficientes nao-simetrica      
c ... solver  =  1 (pcg)       , 2 (gmres)       , 3 (gauss / LDU)
c                4 (bicgstab)  , 5 (block_pcg_it), 6 (bicgstabl2) 
c                7 (minres)                      , 9 (pardiso)
c ... stge    =  1 (csr), 2 (edges), 3 (ebe), 4 (skyline), 6 (csr3)
      unsym   = .false.
      solver  =  1
      stge    =  1
c     block_pu= .true.
      n_blocks_pu = 0 
      block_pu    = .false.
      resid0      =  0.d0
c ... ilib    =  1 define a biblioteca padrão ( default = poromec )
      ilib    =  1
c ... campo gravitacional (Padrao)
      gravity(1)  =   0.0d0
      gravity(2)  =   0.0d0
      gravity(3)  = -9.81d0
      gravity_mod = dsqrt(gravity(1)*gravity(1)
     .                   +gravity(2)*gravity(2)
     .                   +gravity(3)*gravity(3))
c ...
      flag_macro_mesh = .false.
c ... 
      bvtk = .false.
c ... OpenMP
      omp_elmt = .false.
      omp_solv = .false.
      nth_elmt = 1
      nth_solv = 1
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
      endif
c ......................................................................      
      call MPI_barrier(MPI_COMM_WORLD,ierr)
c ......................................................................
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
c ......................................................................
c
c ... Execucao dos macro-comandos:
c
c ......................................................................
c
c ... Macro-comando LOOP:
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
c ......................................................................
c
c ... Macro-comando HEXTOTET:
c
c ......................................................................
  200 continue
      if(my_id.eq.0)print*, 'Macro HEXTOTET'  
      call hexa_to_tetra(ia(i_ix),numel,nen,prename,nplot)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando MESH:
c
c ......................................................................
  300 continue
      if(my_id.eq.0)print*, 'Macro MESH'  
c ... Inicializacao da estrutura de dados de gerenciamento de memoria:
c
      call init_malloc(maxmem)
c ......................................................................
c
c ...
      call init_openmp(omp_elmt,omp_solv,nth_elmt,nth_solv,my_id)
c ......................................................................
c
c
      flag_macro_mesh = .true.
c
c.... Leitura de dados:
c
c      call rdat(nnode,numel,numat,nen,ndf,ndm,nst,i_ix,i_id,i_ie,
c     .          i_nload,i_eload,i_inum,i_e,i_x,i_f,i_u,i_v,i_a,nin)
c
      call rdat_pm(nnode  ,nnodev  ,numel  ,numat  
     .         ,nen       ,nenv   
     .         ,ndf       ,ndm     ,nst    ,i_ix  
     .         ,i_ie      ,i_inum  ,i_e    ,i_x
     .         ,i_id      ,i_nload ,i_eload,i_f  
     .         ,i_u       ,i_u0    ,i_tx0  ,i_dp
     .         ,fstress0  ,fporomec,fmec
     .         ,nin ) 
c
c    -----------------------------------------------------------------
c    | ix | id | ie | nload | eload | inum | e | x | f | u | u0 | tx0 |
c    -----------------------------------------------------------------
c
c    -----------------------------------------------------------------
c    | dp |                                                            
c    -----------------------------------------------------------------
c ......................................................................      
c
c ... calculo das forcas internas devidos as tensoes inicias tensoes
c     habilitado
      if(fstress0) fcstress0 = .true.
c ......................................................................   
c
c ... desabilita o csrc blocado em problemas mecanicos
      if(fmec) then
         block_pu    = .false. 
         n_blocks_pu = 0 
      endif
c ......................................................................   
c
c.... Controle de tempos:
c
      soltime      = 0.d0
      elmtime      = 0.d0
      tformtime    = 0.d0
      vectime      = 0.d0
      ovhtime      = 0.d0
      dottime      = 0.d0
      gvtime       = 0.d0
      allgtime     = 0.d0
      allrtime     = 0.d0
      sendtime     = 0.d0
      matvectime   = 0.d0
      colortime    = 0.d0
      pmatrixtime  = 0.d0
      writetime    = 0.d0
      precondtime  = 0.d0
      ifatsolvtime = 0.d0
      totaltime    = MPI_Wtime()
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
c ... mecanico
      if(fmec) then
        call numeq(ia(i_id),ia(i_inum),ia(i_id),nnode,ndf,neq)
          nequ = 0
          neqp = 0
c ......................................................................
c
c ... poro mecanico
      else if(fporomec) then
c ... primero numera os deslocamento e depois as pressoes
        if(block_pu) then
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
c ......................................................................
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
c ... mecanico
      if(fmec) then
        i_xl  = alloc_8('xl      ',ndm,nenv)
        i_ul  = alloc_8('ul      ',1  ,nst)
        i_pl  = alloc_8('pl      ',1  ,nst)
c ...      
        i_txnl= alloc_8('txnl    ',  6,nen)
c ...      
        i_sl  = alloc_8('sl      ',nst,nst)
        i_ld  = alloc_4('ld      ',  1,nst)
c .....................................................................
c
c     -----------------------------------------
c     | xl | ul | pl | sl | ld |
c     -----------------------------------------
c
c ... poro mecanico
      else if(fporomec) then
        i_xl  = alloc_8('xl      ',ndm,nenv)
        i_ul  = alloc_8('ul      ',1  ,nst)
        i_dpl = alloc_8('dpl     ',1  ,nenv)
        i_pl  = alloc_8('pl      ',1  ,nst)
c ... 6 - tensoes totais, 6 - tensoes de biot , 3 - fluxo de darcy
        i_txl = alloc_8('txl     ', 15,nenv)
c ...      
        i_txnl= alloc_8('txnl    ',  6,nen)
c ...      
        i_sl  = alloc_8('sl      ',nst,nst)
        i_ld  = alloc_4('ld      ',  1,nst)
      endif 
c .......................................................................
c
c     ---------------------------------------------
c     | xl | ul | dpl | pl | txl | txnl | sl | ld |
c     ---------------------------------------------
c
c ... Memoria para a estrutura de dados do sistema de equacoes:
c
      timei = MPI_Wtime()
c ... poromecanico
      if (fporomec) then
         sia = 'ia' 
         sja = 'ja' 
         sau = 'au' 
         sal = 'al' 
         sad = 'ad' 
         call datastruct_pm(ia(i_ix),ia(i_id),ia(i_inum),nnode,nnodev
     .                     ,numel,nen,ndf,nst,neq,nequ,neqp,stge,unsym
     .                     ,nad  ,naduu,nadpp,nadpu
     .                     ,i_ia ,i_ja,i_au  ,i_al,i_ad 
     .                     ,sia  ,sja ,sau   ,sal ,sad    
     .                     ,ovlp ,n_blocks_pu,block_pu  )
c .....................................................................
c
c ... mec
      else if(fmec) then
         sia = 'ia' 
         sja = 'ja' 
         sau = 'au' 
         sal = 'al' 
         sad = 'ad'  
         call datastruct(ia(i_ix),ia(i_id),ia(i_inum),nnode
     .                  ,numel   ,nen     ,ndf       ,nst
     .                  ,neq     ,stge    ,unsym     ,nad  ,nad1
     .                  ,i_ia    ,i_ja    ,i_au      ,i_al ,i_ad 
     .                  ,sia     ,sja     ,sau       ,sal  ,sad         
     .                  ,ovlp     )
      endif
      dstime = MPI_Wtime()-timei
c
c ... colorir a malha (openmp)
c
      colortime = MPI_Wtime()
      call coloredmesh(ia(i_ix),nnode,nnodev,numel,nenv,nen,numcolors
     .               ,i_colorg,i_elcolor)     
      colortime = MPI_Wtime()-colortime
c ......................................................................
c
c ......................................................................
c
c ... Memoria para o vetor de forcas e solucao:
c
c ......................................................................
      if (ndf .gt. 0) then
         i_x0  = alloc_8('x0      ',    1,neq+neq3+neq4)
         i_bst0= alloc_8('bstress0',    1,neq+neq3+neq4)
         i_b0  = alloc_8('b0      ',    1,neq+neq3+neq4)
         i_b   = alloc_8('b       ',    1,neq+neq3+neq4)
         call azero(ia(i_x0)  ,neq+neq3+neq4)
         call azero(ia(i_bst0),neq+neq3+neq4) 
         call azero(ia(i_b0)  ,neq+neq3+neq4)      
         call azero(ia(i_b )  ,neq+neq3+neq4)
c ...
         i_m   = 1
c ...  Memoria para o precondicionador diagonal:
         if(precond .eq. 2 .or. precond .eq. 5) then 
           i_m   = alloc_8('m       ',    1,neq)
           call azero(ia(i_m),neq)
c ......................................................................
c
c ...  Memoria para o precondicionador iLDLt e iLLT (cholesky)
         else if( precond .eq. 3 .or.  precond .eq. 4) then
           i_m   = alloc_8('m       ',    1,neq+nad)
           call azero(ia(i_m),neq+nad)
c ..................................................................... 
         endif       
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
c ...
      ilib   = 1
      i      = 1
      istep  = istep + 1
      resid0 = 0.d0
      t      = t + dt
c .....................................................................
c
c ...
      if(my_id.eq.0)write(*,'(a,i8,a,f15.5,a,f15.5,a,f15.5)')
     .                                  ,' STEP '      ,istep
     .                                  ,' Time(s)'    ,t
     .                                  ,' Time(horas)',t/3600.d0
     .                                  ,' Time(dias)' ,t/86400.d0
c .....................................................................
c
c ... Cargas nodais e valores prescritos no tempo t+dt:
      timei = MPI_Wtime()
      call pload_pm(ia(i_id),ia(i_f),ia(i_u),ia(i_b0),ia(i_nload)
     .            ,nnode,nnodev,ndf)
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ... forcas internas devidos as tensoes inicias tensoes
      if(fstress0 .and. fcstress0) then
          timei = MPI_Wtime()
          call pform_pm(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e)  
     .                 ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
     .                 ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_bst0) 
     .                 ,ia(i_u0)    ,ia(i_dp)     ,ia(i_tx0)
     .                 ,ia(i_xl)    ,ia(i_ul)     ,ia(i_dpl),ia(i_pl)
     .                 ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl) 
     .                 ,numel       ,nen          ,nenv     ,ndf 
     .                 ,ndm         ,nst          
     .                 ,neq         ,nequ         ,neqp 
     .                 ,nad         ,naduu        ,nadpp    ,nadpu 
     .                 ,.false.     ,.true.       ,unsym 
     .                 ,stge        ,5            ,ilib     ,i
     .                 ,ia(i_colorg),ia(i_elcolor),numcolors,.true.
     .                 ,block_pu    ,n_blocks_pu)
          elmtime = elmtime + MPI_Wtime()-timei
          fcstress0= .false.
      endif 
c .....................................................................      
c
c ... forcas de volume e superficie do tempo t+dt e graus de liberade 
c     do passo t:  
      timei = MPI_Wtime()
      call pform_pm(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e) 
     .             ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja) 
     .             ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_b0) 
     .             ,ia(i_u0)    ,ia(i_dp)     ,ia(i_tx0)
     .             ,ia(i_xl)    ,ia(i_ul)     ,ia(i_dpl),ia(i_pl)
     .             ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl) 
     .             ,numel       ,nen          ,nenv     ,ndf 
     .             ,ndm         ,nst          
     .             ,neq         ,nequ         ,neqp 
     .             ,nad         ,naduu        ,nadpp    ,nadpu 
     .             ,.false.     ,.true.       ,unsym 
     .             ,stge,4      ,ilib         ,i
     .             ,ia(i_colorg),ia(i_elcolor),numcolors,.false.
     .             ,block_pu    ,n_blocks_pu)
      elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ... tensao inicial
      if(fstress0) then
        call vsum(ia(i_b0),ia(i_bst0),neq,ia(i_b0))
      endif  
c .....................................................................
c
c ... Preditor: du(n+1,0) = u(n+1) - u(n)  
      timei = MPI_Wtime()
      call predict_pm(nnode,nnodev,ndf,ia(i_u),ia(i_u0))
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ---------------------------------------------------------------------
c loop nao linear:
c ---------------------------------------------------------------------
  410 continue
c ... Loop multi-corretor:      
      if(my_id.eq.0) print*,'nonlinear iteration ',i
c ...
      timei = MPI_Wtime()
      call aequalb(ia(i_b),ia(i_b0),neq)
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ... Residuo: b = F - K.du(n+1,i)
      timei = MPI_Wtime()
      call pform_pm(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e)
     .             ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
     .             ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_b)
     .             ,ia(i_u)     ,ia(i_dp)     ,ia(i_tx0)
     .             ,ia(i_xl)    ,ia(i_ul)     ,ia(i_dpl),ia(i_pl)
     .             ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl)
     .             ,numel       ,nen          ,nenv     ,ndf
     .             ,ndm         ,nst          
     .             ,neq         ,nequ         ,neqp
     .             ,nad         ,naduu        ,nadpp    ,nadpu
     .             ,.true.      ,.true.       ,unsym
     .             ,stge        ,2            ,ilib     ,i
     .             ,ia(i_colorg),ia(i_elcolor),numcolors,.false.
     .             ,block_pu,n_blocks_pu)
      elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ......................................................................
      resid = dsqrt(dot_par(ia(i_b),ia(i_b),neq_dot))
      if(i .eq. 1) resid0 = max(resid0,resid)
      print*,'resid/resid0',resid/resid0,'resid',resid
      if ((resid/resid0) .lt. tol) goto 420     
c ......................................................................            
c
c ... solver (Kdu(n+1,i+1) = b; du(t+dt) )
      timei = MPI_Wtime()
      call get_res(ia(i_u),ia(i_x0),ia(i_id),nnode,nnodev,ndf)
      call solv_pm(neq  ,nequ    ,neqp  
     .         ,nad     ,naduu   ,nadpp      
     .         ,ia(i_ia),ia(i_ja),ia(i_ad)   ,ia(i_al)
     .         ,ia(i_m) ,ia(i_b) ,ia(i_x0)   ,solvtol,maxit
     .         ,ngram   ,block_pu,n_blocks_pu,solver,istep
     .         ,cmaxit  ,ctol    ,alfap      ,alfau ,precond
     .         ,fmec    ,fporomec  
     .         ,neqf1   ,neqf2   ,neq3       ,neq4  ,neq_dot
     .         ,i_fmap  ,i_xf    ,i_rcvs     ,i_dspl)
      soltime = soltime + MPI_Wtime()-timei
c .....................................................................
c
c ... atualizacao :      du(n+1,i+1) = du(n+1,i)      + dv(n+1,i+1)
      timei = MPI_Wtime()
      call update_pm(nnode,nnodev,ndf,ia(i_id),ia(i_u),ia(i_x0))
      vectime = vectime + MPI_Wtime()-timei
c     call update2(nnode,ndf,ia(i_id),ia(i_u),ia(i_u),ia(i_u)
c    .            ,ia(i_x0)  ,0.d0,0.0d0,0.d0)
c .....................................................................
c
c ...
      if (i .ge. maxnlit) goto 420
      i = i + 1
      goto 410 
c .....................................................................
c
c ---------------------------------------------------------------------
c fim do loop nao linear:
c ---------------------------------------------------------------------
c
c ...
  420 continue 
c ... calculo da solucao no tempo t + dt e atualizacao dos valores do 
c     passo de tempo anterior: 
c     u(n+1)  = u(n) + du(n+1) 
c     u(n)    = u(n+1) 
c     dp(n+1) = p(n) - p(0)
      timei = MPI_Wtime()
c     call update_res(nnode  ,nnodev  ,ndf
c    .               ,ia(i_u),ia(i_u0),ia(i_dp),ia(i_pres0))
      call update_res_v2(nnode  ,nnodev  ,ndf
     .                  ,ia(i_u),ia(i_u0),ia(i_dp))
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
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
c ...
      print*, 'Macro PGEO'
      ntn   = 6
c ... Geometria:
      writetime = writetime + MPI_Wtime()-timei 
      call write_mesh_geo_pm(ia(i_ix)   ,ia(i_x)    ,ia(i_ie)
     .                      ,ia(i_id)   ,ia(i_f)    ,ia(i_u) 
     .                      ,ia(i_tx0)  ,ia(i_nload),ia(i_eload)
     .                      ,nnodev     ,numel      ,ndf     ,ntn
     .                      ,nen        ,ndm        ,prename
     .                      ,bvtk       ,macros     ,.true.
     .                      ,nplot      ,nout_face)
      writetime = writetime + MPI_Wtime()-timei
c ......................................................................
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PGEOQ
c
c ......................................................................
  700 continue
      if(my_id.eq.0)  then  
        print*, 'Macro PGEOQ'
c ... tetraedros de 10 nos     
        if( nen .eq. 10 ) then
          i_xq = alloc_8('xq      ',ndm,nnode)
          call mkCoorQuad(ia(i_x) ,ia(i_xq)
     .                   ,ia(i_ix)
     .                   ,numel   ,nen  
     .                   ,nnode   ,nnodev  ,ndm)   
          ntetra10(1:4) = ntetra4(1:4)   
          ntetra4(1:4)  = 0 
          call write_mesh_geo(ia(i_ix)   ,ia(i_xq),nnode   ,numel
     .                     ,nen    ,ndm     ,prename,.false.
     .                     ,.true. ,nplot)
          ntetra4(1:4)  = ntetra10(1:4)
          ntetra10(1:4) = 0
          i_xq = dealloc('xq      ') 
c ... hexaedros de 20 nos     
        else if( nen .eq. 20 ) then
          i_xq = alloc_8('xq      ',ndm,nnode)
          call mkCoorQuad(ia(i_x) ,ia(i_xq)
     .                   ,ia(i_ix)
     .                   ,numel   ,nen  
     .                   ,nnode   ,nnodev  ,ndm)   
          nhexa20(1:4) = nhexa8(1:4)   
          nhexa8(1:4)  = 0 
          call write_mesh_geo(ia(i_ix)   ,ia(i_xq),nnode   ,numel
     .                     ,nen    ,ndm     ,prename,.false.
     .                     ,.true. ,nplot)
          nhexa8(1:4)  = nhexa20(1:4)
          nhexa20(1:4) = 0
          i_xq = dealloc('xq      ')      
        endif
      endif
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: BLOCK_PU
c
c ......................................................................
  800 continue
      if(my_id.eq.0)print*, 'Macro Block_pu'
      if(flag_macro_mesh) then
        print*,'Macro so pode ser utilizada antes da macro mesh'
        goto 5000
      endif
      block_pu = .true.
c ... n_blocks_pu
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =801,end =801) n_blocks_pu   
      goto 50
  801 continue
      print*,'Erro na leitura da macro (BLOCK_PU) n_blocks_pu   !'
      goto 5000
c ......................................................................
c
c ... Macro-comando: GRAVITY
c
c ......................................................................
  900 continue
      if (my_id .eq. 0)   print*, 'Macro Gravity'
c ... gx
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =910,end =910) gravity(1)    
c ... gy
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =920,end =920) gravity(2)    
c ... gz
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =930,end =930) gravity(3)    
      gravity_mod = dsqrt(gravity(1)*gravity(1)
     .                   +gravity(2)*gravity(2)
     .                   +gravity(3)*gravity(3))
      goto 50
  910 continue
      print*,'Erro na leitura da macro (GRAVITY) gx !'
      goto 5000
  920 continue
      print*,'Erro na leitura da macro (GRAVITY) gy !'
      goto 5000
  930 continue
      print*,'Erro na leitura da macro (GRAVITY) gz !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando:       
c
c ......................................................................
 1000 continue
      if(my_id.eq.0)print*, 'Macro PARDISO'
      solver = 9 
      stge   = 6
      goto 50
c ......................................................................
c
c ... Macro-comando: GMRES 
c
c ......................................................................
 1100 continue
      if(my_id.eq.0)print*, 'Macro GMRES'
      solver = 2       
c ... numero de base de krylov
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1101,end =1101) ngram    
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,i10)')'Krylov subspace:',ngram
      endif
c ......................................................................
c
c ... numero maximo de iteracoes
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1102,end =1102) maxit    
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,i10)')'Set max it solver for:',maxit
      endif
c ......................................................................
c
c ... tolerancia 
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1103,end =1103) solvtol   
      if( solvtol .eq. 0.d0) solvtol = smachn()
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,e10.3)')'Set solver tol for:', solvtol  
      endif
c ......................................................................
      goto 50
 1101 continue
      print*,'Erro na leitura da macro (GMRES) ngram !'
      goto 5000
 1102 continue
      print*,'Erro na leitura da macro (GMRES) maxit !'
      goto 5000
 1103 continue
      print*,'Erro na leitura da macro (GMRES) solvtol !'
      goto 5000
c ......................................................................
c
c ... Macro-comando: DELTATC
c
c ......................................................................
 1200 continue
      if(my_id.eq.0) print*, 'Macro DELTATC'
      if(fporomec) then
          call deltat_critico(ia(i_ix),ia(i_eload),ia(i_ie),ia(i_e)
     .                       ,ia(i_x), ia(i_xl)
     .                       ,numel,nen,nenv,ndf
     .                       ,ndm,nst,1,ilib)
      endif
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
 1300 continue
      if(my_id.eq.0) print*, 'Macro PCOO'
c ... csr divido em bloco Kuu, kpp e kpu
      if(block_pu) then
c ...
        i_lin    = alloc_4('lincoo  ',1,neq+2*(nad+nadpu)) 
        i_col    = alloc_4('colcoo  ',1,neq+2*(nad+nadpu)) 
        i_acoo   = alloc_8('acoo    ',1,neq+2*(nad+nadpu))
c ... Kuu
        i_linuu  = alloc_4('lincoouu',1,nequ+2*naduu) 
        i_coluu  = alloc_4('colcoouu',1,nequ+2*naduu) 
        i_acoouu = alloc_8('acoouu  ',1,nequ+2*naduu)
c ... Kpp
        i_linpp  = alloc_4('lincoopp',1,neqp+2*nadpp) 
        i_colpp  = alloc_4('colcoopp',1,neqp+2*nadpp) 
        i_acoopp = alloc_8('acoopp  ',1,neqp+2*nadpp)
c ......................................................................
c
c ...
        call csr_block_to_coo(ia(i_lin)  ,ia(i_col)  ,ia(i_acoo)
     .                       ,ia(i_linuu),ia(i_coluu),ia(i_acoouu)
     .                       ,ia(i_linpp),ia(i_colpp),ia(i_acoopp)
     .                       ,ia(i_ia)   ,ia(i_ja)
     .                       ,ia(i_al)   ,ia(i_ad)
     .                       ,neq        ,nequ   ,nad, nadpu            
     .                       ,.false.)
c ......................................................................
c
c ...
        fname   = name(prename,0,50)
        call write_coo(ia(i_lin),ia(i_col),ia(i_acoo)
     .               ,ia(i_b  ),neq      ,neq+2*(nad+nadpu)  
     .               ,fname    ,nout   
     .               ,.true.   ,.true.   )
c ......................................................................
c
c ...
        fname   = name(prename,0,51)
        call write_coo(ia(i_linuu),ia(i_coluu),ia(i_acoouu)
     .               ,ia(i_b  )   ,nequ       ,nequ+2*(naduu)  
     .               ,fname       ,nout   
     .               ,.false.     ,.true.     )
c ......................................................................
c
c ...
        fname   = name(prename,0,52)
        call write_coo(ia(i_linpp),ia(i_colpp),ia(i_acoopp)
     .               ,ia(i_b  )   ,neqp       ,neqp+2*(nadpp)  
     .               ,fname       ,nout   
     .               ,.false.     ,.true.     )
c ......................................................................
c
c ...
        i_acoopp = dealloc('acoopp  ')
        i_colpp  = dealloc('colcoopp') 
        i_linpp  = dealloc('lincoopp')
c 
        i_acoouu = dealloc('acoouu  ')
        i_coluu  = dealloc('colcoouu') 
        i_linuu  = dealloc('lincoouu') 
c
        i_acoo   = dealloc('acoo    ')
        i_col    = dealloc('colcoo  ') 
        i_lin    = dealloc('lincoo  ') 
c ......................................................................
c
c ...
      else
c ...
        fname   = name(prename,0,50)
c ......................................................................
c
c ...
        i_lin  = alloc_4('lincoo  ',1,neq+2*nad) 
        i_col  = alloc_4('colcoo  ',1,neq+2*nad) 
        i_acoo = alloc_8('acoo    ',1,neq+2*nad)
c ......................................................................
c
c ...
        call csr_to_coo_pm(ia(i_lin) ,ia(i_col) ,ia(i_acoo)
     .                    ,ia(i_ia)  ,ia(i_ja)
     .                    ,ia(i_al)  ,ia(i_ad)  
     .                    ,neq       ,nad       ,.false.)
c ......................................................................
c
c ...
        call write_coo(ia(i_lin),ia(i_col),ia(i_acoo)
     .                   ,ia(i_b  ),neq      ,neq+2*nad  
     .                   ,fname    ,nout   
     .                   ,.true.   ,.true.   )
c ......................................................................
c
c ...
        i_acoo = dealloc('acoo    ')
        i_col  = dealloc('colcoo  ') 
        i_lin  = dealloc('lincoo  ') 
c ......................................................................
      endif
c ......................................................................
c
c ...      
      goto 50     
c ----------------------------------------------------------------------
c
c ... Macro-comando: BICGSTAB
c
c ......................................................................
 1400 continue
      if(my_id.eq.0)print*, 'Macro BICGSTAB'
      solver = 4 
c ... numero maximo de iteracoes
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1402,end =1402) maxit    
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,i10)')'Set max it solver for:',maxit
      endif
c ......................................................................
c
c ... tolerancia 
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1403,end =1403) solvtol   
      if( solvtol .eq. 0.d0) solvtol = smachn()
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,e10.3)')'Set solver tol for:', solvtol  
      endif
c ......................................................................
c
c ... precondicionador
      call readmacro(nin,.false.)
      write(string,'(6a)') (word(i),i=1,6)
      call set_precond(word,precond,nin,my_id)  
c ......................................................................
      goto 50
 1402 continue
      print*,'Erro na leitura da macro (BICGSTAB) maxit !'
      goto 5000
 1403 continue
      print*,'Erro na leitura da macro (BICGSTAB) solvtol !'
      goto 5000
c ----------------------------------------------------------------------
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PCG
c
c ......................................................................
 1500 continue
      if(my_id.eq.0)print*,'Macro PCG'
      solver = 1       
c ... numero maximo de iteracoes
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1502,end =1502) maxit    
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,i10)')'Set max it solver for:',maxit
      endif
c ......................................................................
c
c ... tolerancia 
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1503,end =1503) solvtol  
      if( solvtol .eq. 0.d0) solvtol = smachn() 
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,e10.3)')'Set solver tol for:', solvtol  
      endif     
c ......................................................................
c
c ... precondicionador
      call readmacro(nin,.false.)
      write(string,'(6a)') (word(i),i=1,6)
      call set_precond(word,precond,nin,my_id)  
c ......................................................................
      goto 50
 1502 continue
      print*,'Erro na leitura da macro (PCG) maxit !'
      goto 5000
 1503 continue
      print*,'Erro na leitura da macro (PCG) solvtol !'
      goto 5000
c ----------------------------------------------------------------------
      goto 50 
c ----------------------------------------------------------------------
c
c ... Macro-comando: PRES     
c
c ......................................................................
 1600 continue
      if(my_id.eq.0)print*,'Macro PRES'
c ... calculo da tensoes, tensoes efetivas e fluxo de darcy nos vertices.
      ntn   = 6
c .....................................................................
      i_tx  = alloc_8('tx      ',  ntn,nnodev)
      i_txe = alloc_8('txe     ',  ntn,nnodev)
      i_txb = alloc_8('txb     ',  ntn,nnodev)
      i_flux= alloc_8('flux    ',  ndm,nnodev)
      i_ic  = alloc_4('ic      ',    1,nnodev)
c .....................................................................
c
c ...
      timei = MPI_Wtime()
      call tform_pm(ia(i_ix) ,ia(i_x)  ,ia(i_e)  ,ia(i_ie)
     .             ,ia(i_ic) ,ia(i_xl) ,ia(i_ul) ,ia(i_dpl)
     .             ,ia(i_txl),ia(i_u)  ,ia(i_dp) ,ia(i_tx0)
     .             ,ia(i_tx) ,ia(i_txb),ia(i_txe),ia(i_flux)
     .             ,nnodev   ,numel   ,nen       ,nenv
     .             ,ndm      ,ndf     ,nst       ,ntn
     .             ,3        ,ilib)
      tformtime = tformtime + MPI_Wtime()-timei
c ......................................................................
c
c ...
      fname = name(prename,istep,2)
      open(nplot,file=fname)
      call write_mesh_res_pm(ia(i_ix),ia(i_x)  ,ia(i_u)  ,ia(i_dp)
     .                      ,ia(i_tx),ia(i_txb),ia(i_txe),ia(i_flux)
     .                      ,nnodev  ,numel    ,istep    ,t
     .                      ,nen     ,ndm      ,ndf      ,ntn
     .                      ,fname   ,.false.,.true.     ,nplot)
      close(nplot)  
c ......................................................................
c
c ...
      i_ic  = dealloc('ic      ')
      i_flux= dealloc('flux    ')
      i_txe = dealloc('txb     ')
      i_txe = dealloc('txe     ')
      i_tx  = dealloc('tx      ')
c ......................................................................
      goto 50     
c ......................................................................
c
c ... Macro-comando: SPCGM        
c
c ......................................................................
 1700 continue
      if(my_id.eq.0)print*,'Macro SPCGM'
c ...
      fname = name(prename,nprcs,53)
      open(logsolvd,file=fname)
c ......................................................................
      solver = 5       
c ... fator de sobre-relaxação p
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1703,end =1701) alfap    
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,e10.3)')'Set alfaP :',alfap
      endif
c ......................................................................
c
c ... fator de sobre-relaxação p
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1702,end =1702) alfau    
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,e10.3)')'Set alfaU :',alfau
      endif
c ......................................................................
c
c ... numero maximo de iteracoes
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1703,end =1703) maxit    
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,i10)')'Set max it solver for:',maxit
      endif
c ......................................................................
c
c ... tolerancia 
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1704,end =1704) solvtol   
      if( solvtol .eq. 0.d0) solvtol = smachn()
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,e10.3)')'Set solver tol for:', solvtol  
      endif
c ......................................................................
c
c
c ... tolerancia no ciclo
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1705,end =1705) cmaxit   
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,i10)')'Set max cycles it solver for:',cmaxit
      endif
c ......................................................................
c
c ... tolerancia no ciclo
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1702,end =1702) ctol   
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,e10.3)')'Set cycle tol for:', ctol 
      endif
c ......................................................................
      goto 50
 1701 continue
      print*,'Erro na leitura da macro (SPCG) alfaP   !'
      goto 5000
 1702 continue
      print*,'Erro na leitura da macro (SPCG) alfaU!'
      goto 5000
 1703 continue
      print*,'Erro na leitura da macro (SPCG) maxit !'
      goto 5000
 1704 continue
      print*,'Erro na leitura da macro (SPCG) solvtol !'
      goto 5000
 1705 continue
      print*,'Erro na leitura da macro (SPCG) cmaxit !'
      goto 5000
 1706 continue
      print*,'Erro na leitura da macro (SPCG) ctol !'
      goto 5000
c ----------------------------------------------------------------------
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SOLVM
c
c ......................................................................
 1800 continue
      if (my_id .eq. 0 ) print*, 'Macro  SOLVM'
c ...
      ilib   = 1
      i      = 1
      istep  = istep + 1
      resid0 = 0.d0
c .....................................................................
c
c ... Cargas nodais e valores prescritos no tempo t+dt:
      timei = MPI_Wtime()
      call pload_mec(ia(i_id),ia(i_f),ia(i_u),ia(i_b0),ia(i_nload)
     .              ,nnode   ,ndf)
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ... forcas internas devidos as tensoes inicias tensoes
      if(fstress0 .and. fcstress0) then
          timei = MPI_Wtime()
          call pform_mec(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e) 
     .                 ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
     .                 ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_bst0) 
     .                 ,ia(i_u0)    ,ia(i_tx0)
     .                 ,ia(i_xl)    ,ia(i_ul)     ,ia(i_pl)
     .                 ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl) 
     .                 ,numel       ,nen          ,nenv     ,ndf 
     .                 ,ndm         ,nst          ,neq      ,nad        
     .                 ,.false.     ,.true.       ,unsym 
     .                 ,stge        ,5            ,ilib     ,i
     .                 ,ia(i_colorg),ia(i_elcolor),numcolors,.true.)
          elmtime = elmtime + MPI_Wtime()-timei
          fcstress0= .false.
      endif 
c .....................................................................      
c
c ... forcas de volume e superficie do tempo t+dt :  
      timei = MPI_Wtime()
      call pform_mec(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e) 
     .              ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
     .              ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_b0) 
     .              ,ia(i_u0)    ,ia(i_tx0)
     .              ,ia(i_xl)    ,ia(i_ul)     ,ia(i_pl)
     .              ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl) 
     .              ,numel       ,nen          ,nenv     ,ndf 
     .              ,ndm         ,nst          
     .              ,neq         ,nad           
     .              ,.false.     ,.true.       ,unsym 
     .              ,stge,4      ,ilib         ,i
     .              ,ia(i_colorg),ia(i_elcolor),numcolors,.false.)
      elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ... tensao inicial
      if(fstress0) then
        call vsum(ia(i_b0),ia(i_bst0),neq,ia(i_b0))
      endif  
c .....................................................................
c
c ---------------------------------------------------------------------
c loop nao linear:
c ---------------------------------------------------------------------
 1810 continue
c ... Loop multi-corretor:      
      if(my_id.eq.0) print*,'nonlinear iteration ',i
c ...
      timei = MPI_Wtime()
      call aequalb(ia(i_b),ia(i_b0),neq)
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ... Residuo: b = F - K.u(n+1,i)
      timei = MPI_Wtime()
      call pform_mec(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e)
     .              ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
     .              ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_b)
     .              ,ia(i_u)     ,ia(i_tx0)
     .              ,ia(i_xl)    ,ia(i_ul)     ,ia(i_pl)
     .              ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl)
     .              ,numel       ,nen          ,nenv     ,ndf
     .              ,ndm         ,nst          
     .              ,neq         ,nad        
     .              ,.true.      ,.true.       ,unsym
     .              ,stge        ,2            ,ilib     ,i
     .              ,ia(i_colorg),ia(i_elcolor),numcolors,.false.)
      elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ......................................................................
      resid = dsqrt(dot_par(ia(i_b),ia(i_b),neq_dot))
      if(i .eq. 1) resid0 = max(resid0,resid)
      print*,'resid/resid0',resid/resid0,'resid',resid
      if ((resid/resid0) .lt. tol) goto 1820     
c ......................................................................            
c
c ... solver (Ku(n+1,i+1) = b; u(t+dt) )
      timei = MPI_Wtime()
      call solv_pm(neq     ,nequ    ,neqp  
     .            ,nad     ,naduu   ,nadpp      
     .            ,ia(i_ia),ia(i_ja),ia(i_ad)    ,ia(i_al)
     .            ,ia(i_m) ,ia(i_b) ,ia(i_x0)    ,solvtol,maxit
     .            ,ngram   ,block_pu,n_blocks_pu ,solver ,istep
     .            ,cmaxit  ,ctol    ,alfap       ,alfau  ,precond 
     .            ,fmec    ,fporomec  
     .            ,neqf1   ,neqf2   ,neq3        ,neq4   ,neq_dot
     .            ,i_fmap  ,i_xf    ,i_rcvs      ,i_dspl)
      soltime = soltime + MPI_Wtime()-timei
c .....................................................................
c
c ... atualizacao :      u(n+1,i+1) = x(n+1,i)
      timei = MPI_Wtime()
      call update_mec(nnode,ndf,ia(i_id),ia(i_u),ia(i_x0))
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ...
      if (i .ge. maxnlit) goto 1820
      i = i + 1
      goto 1810 
c .....................................................................
c
c ---------------------------------------------------------------------
c fim do loop nao linear:
c ---------------------------------------------------------------------
c
c ...
 1820 continue 
c .....................................................................
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PMECRES
c
c ......................................................................
 1900 continue
      if(my_id.eq.0)print*,'Macro PMECRES'
c ... calculo da tensoes, tensoes efetivas e fluxo de darcy nos vertices.
      ntn   = 6
c .....................................................................
      i_tx  = alloc_8('tx      ',  ntn,nnodev)
      i_ic  = alloc_4('ic      ',    1,nnodev)
c .....................................................................
c
c ...
      timei = MPI_Wtime()
      call tform_mec(ia(i_ix)   ,ia(i_x)  ,ia(i_e)   ,ia(i_ie)
     .               ,ia(i_ic)  ,ia(i_xl) ,ia(i_ul) 
     .               ,ia(i_txnl),ia(i_u)  ,ia(i_tx0),ia(i_tx) 
     .               ,nnodev    ,numel   ,nen       ,nenv
     .               ,ndm       ,ndf     ,nst       ,ntn
     .               ,3         ,ilib)
      tformtime = tformtime + MPI_Wtime()-timei
c ......................................................................
c
c ...
      fname = name(prename,istep,2)
      open(nplot,file=fname)
      call write_mesh_res_mec(ia(i_ix),ia(i_x)  ,ia(i_u),ia(i_tx)
     .                       ,nnodev  ,numel
     .                       ,nen     ,ndm      ,ndf   ,ntn
     .                       ,fname   ,.false.,.true.  ,nplot)
      close(nplot)  
c ......................................................................
c
c ...
      i_ic  = dealloc('ic      ')
      i_tx  = dealloc('tx      ')
c ......................................................................
      goto 50     
c ----------------------------------------------------------------------
c
c ... Macro-comando: BICGSTABL2
c
c ......................................................................
 2000 continue
      if(my_id.eq.0)print*, 'Macro BICGSL2'
      solver = 6 
c ... numero maximo de iteracoes
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2002,end =2002) maxit    
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,i10)')'Set max it solver for:',maxit
      endif
c ......................................................................
c
c ... tolerancia 
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2003,end =2003) solvtol   
      if( solvtol .eq. 0.d0) solvtol = smachn()
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,e10.3)')'Set solver tol for:', solvtol  
      endif
c ......................................................................
c
c ... precondicionador
      call readmacro(nin,.false.)
      write(string,'(6a)') (word(i),i=1,6)
      call set_precond(word,precond,nin,my_id)  
c ......................................................................
      goto 50
 2002 continue
      print*,'Erro na leitura da macro (BICGSTABL2) maxit !'
      goto 5000
 2003 continue
      print*,'Erro na leitura da macro (BICGSTABL2) solvtol !'
      goto 5000
c ----------------------------------------------------------------------
      goto 50
c ----------------------------------------------------------------------
c
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2100 continue
      if(my_id.eq.0)print*, 'Macro MINRES'
      solver = 7 
c ... numero maximo de iteracoes
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2102,end =2102) maxit    
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,i10)')'Set max it solver for:',maxit
      endif
c ......................................................................
c
c ... tolerancia 
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2103,end =2103) solvtol   
      if( solvtol .eq. 0.d0) solvtol = smachn()
      if(my_id.eq.0) then
        write(*,'(1x,a25,1x,e10.3)')'Set solver tol for:', solvtol  
      endif
c ......................................................................
c
c ... precondicionador
      call readmacro(nin,.false.)
      write(string,'(6a)') (word(i),i=1,6)
      call set_precond(word,precond,nin,my_id)  
c ......................................................................
      goto 50
 2102 continue
      print*,'Erro na leitura da macro (MINRES) maxit !'
      goto 5000
 2103 continue
      print*,'Erro na leitura da macro (MINRES) solvtol !'
      goto 5000
c ----------------------------------------------------------------------
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
      write(*,'(a,i10)')' Set max nonlinear it for ', maxnlit
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
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: NLSTATIC
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
      if(my_id.eq.0) print*, 'Macro PNUP'      
      if(flag_pnd.eqv..false.) then
        if(my_id.eq.0)print*,'Nemhum no de impressao para PNUP!'  
        call stop_mef()
      endif
c ... codigo para o arquivo up_node.txt      
      code   = 30
      ifiles = 1
c .....................................................................
      call global_v(ndf,nno_pload,i_u,i_g1,'dispG   ')
c     call global_v(  1,nno_pload,i_dp,i_g1,'dispG   ')
      string = 'DeslocAndPress'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          call printnode(ia(i_g1),ia(i_no+j-1),ndf            ,istep,dt
     .                  ,string  ,prename     ,ia(i_nfile+j-1)
     .                  ,code    ,new_file(ifiles))
        enddo
        new_file(ifiles) = .false.
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
      if(my_id.eq.0) print*, 'Macro PNSF    '
      if(flag_pnd.eqv..false.) then
        if(my_id.eq.0)print*,'Nemhum no de impressao para PNUP!'  
        call stop_mef()
      endif
c ... calculo da tensoes, tensoes efetivas e fluxo de darcy nos vertices.
      ntn   = 6
c .....................................................................
      i_tx  = alloc_8('tx      ',  ntn,nnodev)
      i_txe = alloc_8('txe     ',  ntn,nnodev)
      i_txb = alloc_8('txb     ',  ntn,nnodev)
      i_flux= alloc_8('flux    ',  ndm,nnodev)
      i_ic  = alloc_4('ic      ',    1,nnodev)
c .....................................................................
c
c ...
      timei = MPI_Wtime()
      call tform_pm(ia(i_ix) ,ia(i_x)  ,ia(i_e)  ,ia(i_ie)
     .             ,ia(i_ic) ,ia(i_xl) ,ia(i_ul) ,ia(i_dpl)
     .             ,ia(i_txl),ia(i_u)  ,ia(i_dp) ,ia(i_tx0)
     .             ,ia(i_tx) ,ia(i_txb),ia(i_txe),ia(i_flux)
     .             ,nnodev   ,numel    ,nen      ,nenv
     .             ,ndm      ,ndf      ,nst      ,ntn
     .             ,3        ,ilib)
      tformtime = tformtime + MPI_Wtime()-timei
c ......................................................................
c
c ... codigo para o arquivo stress_node.txt      
      code   = 31
      ifiles = 2
      string = 'stressTotal'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          call printnode(ia(i_tx),ia(i_no+j-1),ntn            ,istep,dt
     .                 ,string   ,prename     ,ia(i_nfile+j-1)
     .                 ,code     ,new_file(ifiles))
        enddo
        new_file(ifiles) = .false.
      endif
c ......................................................................
c
c ... codigo para o arquivo stressE_node.txt      
      code   = 32
      ifiles = 3
      string = 'stressE'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          call printnode(ia(i_txe),ia(i_no+j-1),ntn            ,istep,dt
     .                 ,string   ,prename     ,ia(i_nfile+j-1)
     .                 ,code     ,new_file(ifiles))
        enddo
        new_file(ifiles) = .false.
      endif
c ......................................................................
c
c ... codigo para o arquivo stressB_node.txt      
      code   = 33
      ifiles = 4
      string = 'stressBiot'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          call printnode(ia(i_txb),ia(i_no+j-1),ntn            ,istep,dt
     .                 ,string   ,prename     ,ia(i_nfile+j-1)
     .                 ,code     ,new_file(ifiles))
        enddo
        new_file(ifiles) = .false.
      endif
c ......................................................................
c
c ... codigo para o arquivo flux_node.txt      
      code   = 34
      ifiles = 5
      string = 'stressFlux'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          call printnode(ia(i_flux),ia(i_no+j-1),ndm          ,istep,dt
     .                 ,string   ,prename     ,ia(i_nfile+j-1)
     .                 ,code     ,new_file(ifiles))
        enddo
        new_file(ifiles) = .false.
      endif
c ......................................................................
c
c ...
      i_ic  = dealloc('ic      ')
      i_flux= dealloc('flux    ')
      i_txe = dealloc('txb     ')
      i_txe = dealloc('txe     ')
      i_tx  = dealloc('tx      ')
c ......................................................................
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: CONFIG
c
c ......................................................................
 3700 continue
      print*, 'Macro CONFIG    '
      if(flag_macro_mesh) then
        print*,'Macro so pode ser utilizada antes da macro mesh'
        goto 5000
      endif
      call read_config(maxmem
     .                ,omp_elmt,omp_solv
     .                ,nth_elmt,nth_solv
     .                ,reordf  ,bvtk 
     .                ,nin)
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
      close(logsolv)
      if(solver .eq. 5 ) close(logsolvd)
c .....................................................................
c
c ...
      totaltime = MPI_Wtime() - totaltime
c
c ... arquivo de tempo      
      call write_log_file(nnode    ,numel   ,numel_nov,numel_ov,ndf 
     .                   ,neq      ,nequ    ,neqp     ,neq1    ,neq2
     .                   ,neq32    ,neq4    ,neq1a    ,neqf1   ,neqf2 
     .                   ,nad      ,naduu   ,nadpp    ,nadpu   ,nad1
     .                   ,omp_elmt ,nth_elmt,omp_solv ,nth_solv
     .                   ,fporomec ,fmec    ,numcolors,prename
     .                   ,my_id    ,nprcs   ,nout)
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
c     call matvec_csrcb(neq,nequ,ia,ja,ia(neq+2),ja(nad+1)
c    .                 ,ad ,al  ,al(nad+1),b,x 
c    .                 ,dum,dum,idum,idum,idum,idum)
c ...
      do i = 1, neq
        print*,i,x(i),b(i)
      enddo
c .....................................................................
      return
      end
c      
      subroutine printv (v,n)
      implicit none 
      real*8 v(*)
      integer i,n
      do i = 1, n
      print*,i,v(i)
      enddo
      return
      end

