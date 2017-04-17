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
      include 'precond.fi'
      include 'transiente.fi'
      include 'parallel.fi'
      include 'gravity.fi'
      include 'elementos.fi'
      include 'time.fi'
      include 'openmp.fi'
      include 'termprop.fi'
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
      character*80 pnodename,ppiname
      integer nin,nincl,naux,nplot,nout,nout_face,nout_nonlinear
      integer logsolv,fconf,logsolvd,log_hist_solv
      integer totfiles,openflag
      integer*8 i_no,i_nfile,i_pi,i_nfile_pi
      integer print_nnode
      integer num_pnode,num_pel
c ... arquivo de impressao nos nos ( pu,stress,stressE,stressB,flux,...)  
      integer nfiles,ifiles,num_print,n_opt_print_flag
      parameter ( nfiles = 10,n_opt_print_flag = 20)
      logical new_file(nfiles),new_file_pi(nfiles),flag_pnd,flag_ppi
      logical print_flag(n_opt_print_flag),fhist_log,legacy_vtk
c ......................................................................
c
c ... solver
      logical fprint
c ......................................................................
c
c ... propriedades variaveis
      logical vprop(4)
c ......................................................................
c
c ... Variaveis de controle de solucao:
c
      integer maxit,maxnlit,tmaxnlit,ngram,stge,solver,istep
      integer ilib,ntn,code,istop,stop_crit
      real*8  tol,solvtol,resid,resid0
      logical reordf,unsym,lhs,newton_raphson
c ... pcg duplo
      integer cmaxit
      real*8  ctol,alfap,alfau
c.......................................................................
c
c ... Variaveis descritivas do problema:
c
      integer nnodev,nnode,numel,numat,nen,nenv,ndf,ndm,nst
      logical fporomec,fmec,fplastic
c ......................................................................
c
c ... Variaveis do sistema de equacoes:
      integer neq,nequ,neqp,naduu,nadpp,nadpu
      integer*8 nad
      integer n_blocks_pu
      logical block_pu,block_pu_sym,fneqs
      character*8 sia,sja,sau,sal,sad
c .....................................................................
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
      integer*8 i_ic,i_fnno
c ... arranjos locais ao elemento
      integer*8 i_xl,i_ul,i_pl,i_sl,i_ld,i_dpl,i_txl,i_txnl,i_plasticl
      integer*8 i_tx1pl,i_tx2pl,i_depsl,i_porosity,i_dporo,i_p0l
      integer*8 i_vpropell
c ... forcas e graus de liberdade 
      integer*8 i_f
      integer*8 i_u,i_u0,i_tx0,i_tx1p,i_tx2p,i_plastic,i_depsp,i_dp
      integer*8 i_tx,i_txb,i_txe,i_flux,i_pc,i_elplastic
c ... sistema de equacoes
      integer*8 i_ia,i_ja,i_ad,i_au,i_al,i_b,i_b0,i_x0,i_bst0
c ... precondicionador
      integer*8 i_m
c ... variacao das propriedades por elemento
      integer*8 i_vpropel
c ... arranjos globais (MPI - escrita)
      integer*8 i_g,i_g1,i_g2,i_g3,i_g4,i_g5,i_g6,i_g7,i_g8,i_g9,i_g10
      integer*8 i_g11
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
      data macro/'loop    ','hextotet','mesh    '
     1          ,'solv    ','dt      ','pgeo    '
     2          ,'presolv ','block_pu','gravity '
     3          ,'conseq  ','solver  ','deltatc '
     4          ,'pcoo    ','        ','        '
     5          ,'pres    ','        ','solvm   '
     6          ,'pmecres ','        ','        '
     7          ,'        ','        ','        '
     8          ,'        ','nl      ','        '
     9          ,'        ','        ','return  '
     1          ,'setpnode','setprint','setpi   '
     2          ,'ppi     ','pnup    ','pnsf    '
     3          ,'config  ','neqs    ','run     '
     4          ,'stop    '/
c ......................................................................
c
c ... Arquivos de entrada e saida:
c
      data nin /1/, nplot /3/, nout_nonlinear /4/ , nincl / 7/
     .    , logsolv /10/, nout /15/,logsolvd /16/, nout_face /17/
     .    ,log_hist_solv /18/ 
      data fconf /5/
      data flag_pnd /.false./,flag_ppi /.false./
c     arquivo de impressao de nos associados aos seguintes inteiros
c     nfile    = 51 ,51 , 52,...,num_pnode,...,100
c     nfile_pi = 101,102,103,...,num_pel  ,...,150
c ......................................................................      
c
c ... Inicializacao MPI:
c
      mpi = .false.
      call MPI_INIT( ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprcs, ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
      if(nprcs .gt. 1) mpi = .true.
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
c ... escrita de variavies no vtk
c ... malha quadratica (1)
c ... deslocamento     (2)
c ... pressao          (3)
c ... delta pressa     (4)
c ... stress Total     (5)
c ... stress Biot      (6)
c ... stress Terzaghi  (7)
c ... fluxo de darcy   (8)
c ... delta prosidade  (9)
c ... pconsolidation  (10)
c ... plastic elemt   (11)
      print_flag(1) = .false. 
      print_flag(2) = .true.
      print_flag(3) = .true.  
      print_flag(4) = .false.  
      print_flag(5) = .false.  
      print_flag(6) = .false.  
      print_flag(7) = .false.
      print_flag(8) = .false. 
      print_flag(9) = .false.
      print_flag(10)= .false. 
      print_flag(11)= .false. 
c ... tipo do problema
c ... fporomec  = problema poromecanico                    
c ... fmec      = problema mecanico              
c ... fplastic  = plasticidade
c ... vprop     = propriedades variaveis
c             1 - prop por pontos de integracao (true|false)                          
c             2 - konzey-Caraman                (true|false)                          
c             3 - mecanico                      (true|false)                          
c             4 - massa especifica              (true|false)                          
      fporomec   = .false.
      fmec       = .false.
      fplastic   = .false.
      vprop(1)   = .false.
      vprop(2)   = .false.
      vprop(3)   = .false.
      vprop(4)   = .false.
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
c ... resid_type = tipo do criterio de para do algoritmo nao-linear
c                  1 - |R|/|R0| < tol
c                  2 - |Ru|/|Ru0| < tol
c                      |Rp|/|Rp0| < tol
c                  3 - |du|/|du0| < tol
c                      |dp|/|dp0| < tol
c ... ngram   =  base de Krylov (utilizada somente no gmres)
c ... precond =  1 - NONE , 2 - diag, 3 - iLDLt(0), 4 - iC(0)
c                5 - diagm, 6 - bdiag, 7 -diagS
      maxit          =  50000
      solvtol        =  1.d-11
      maxnlit        =  500
      tol            =  1.d-04
      stop_crit      =  2
      ngram          =  150
      precond        =  2
      fhist_log      = .false.
      newton_raphson = .false.
c ... cmaxit  =  numero max. de iteracoes do ciclo externo do pcg duplo
c ... ctol    =  tolerancia do ciclo externo do pcg duplo
      cmaxit  =  200
      ctol    =  1.d-6
c ... unsym   =  true -> matriz de coeficientes nao-simetrica      
c ... solver  =  1 (pcg)       , 2 (gmres)       , 3 (gauss / LDU)
c                4 (bicgstab)  , 5 (block_pcg_it), 6 (bicgstabl2) 
c                7 (minres)    , 8 (pcr)         , 9 (symmlq)
c               10 (pardiso)   ,11 (sqmr)
c ... stge    =  1 (csr), 2 (edges), 3 (ebe), 4 (skyline), 6 (csr3)
      unsym   = .false.
      solver  =  11
      stge    =  1
c ... Matriz blocada 
c     block_pu = true e n_blocks_up = 1 monta o bloco 
c     kuu e kpp separados e nao monta o bloco kup           
c                                                                     
c     block_pu = true e n_blocks_up = 2 monta o bloco 
c     kuu e kpp juntos e o bloco kup separado              
c                                                                    
c     block_pu = true e n_blocks_up = 3 monta o bloco 
c     kuu, kpp e kup separados                                   
c 
c     block_pu_sym = true  monta o bloco kuu, kpp e kup juntos na 
c     forma simetrica
c     ( primeiras equacoes u e depois as esquacoes de pressao)  
c 
c     block_pu_sym = false e block_pu = false                               
c     monta o bloco kuu, kpp e kup juntos sem considera a estrutura      
c     blocada, i.e., graus de liberda juntos                      
      n_blocks_pu  = 0 
      block_pu     = .false.
      block_pu_sym = .false.
c
      resid0  =  0.d0
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
c ... sai da paraview
c     bvtk = true vtk binario
c     leagacy_vtk = true (.vtk) | false (.vtu) 
      bvtk       = .false.
      legacy_vtk = .false.
c ... OpenMP
      omp_elmt = .false.
      omp_solv = .false.
      nth_elmt = 1
      nth_solv = 1
c ... calculo do numero de equacoes
      fneqs = .false.
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
          print*, 'Input file:'
          read(*,'(a)') filein
        endif 
      endif      
      if (nprcs .eq. 1) then
         open(nin, file= filein, status= 'old', err=15, action= 'read')
         goto 20
   15    continue
         print*,'File ',trim(filein),' not found !'
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
            if (my_id .eq. 0) print*,'File ',trim(fname),' not found !'
            call stop_mef()
         endif
      endif
   20 continue
c ......................................................................   
      if (my_id .eq. 0) then   
c ... interface de linha de comando        
        if(nargs .eq. 2) then
          call getarg(2,arg)
          prename = arg
        else
          print*, 'Output file: '
          read(*,'(a)') prename
        endif
c ... log do solver 
        fname = name(prename,nprcs,15)
        open(logsolv,file=fname)
        write(logsolv,'(a)') 'Solver control flop.'
c ... log do nao linear 
        fname = name(prename,nprcs,16)
        open(nout_nonlinear,file=fname)
        write(nout_nonlinear,'(a)') 'Non-linear control flop.'
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
      goto (100 , 200, 300 !'loop    ','hextotet','mesh    '
     1     ,400 , 500, 600 !'solv    ','dt      ','pgeo    '
     2     ,700 , 800, 900 !'        ','block_pu','gravity '
     3     ,1000,1100,1200 !'conseq  ','solver  ','deltatc '
     4     ,1300,1400,1500 !'pcoo    ','        ','        '
     5     ,1600,1700,1800 !'pres    ','        ','solvm   '
     6     ,1900,2000,2100 !'pmecres ','        ','        '
     7     ,2200,2300,2400 !'        ','        ','        '
     8     ,2500,2600,2700 !'        ','nl      ','        '
     9     ,2800,2900,3000 !'        ','        ','return  '
     1     ,3100,3200,3300 !'setpnode','setprint','setpg   '
     2     ,3400,3500,3600 !'ppi     ','pnup    ','pnsf    '
     3     ,3700,3800,3900 !'config  ','neqs    ','run     '
     4     ,5000) j        !'stop    ' 
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
      print*,'Error reading macro (LOOP) !'
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
      call rdat_pm(nnode  ,nnodev  ,numel      ,numat  
     1         ,nen       ,nenv    ,ntn        ,ndf
     1         ,ndm       ,nst     ,i_ix       ,i_ie   
     2         ,i_inum    ,i_e    ,i_x
     3         ,i_id      ,i_nload ,i_eload    ,i_f  
     4         ,i_u       ,i_u0    ,i_tx0      ,i_dp
     5         ,i_tx1p    ,i_tx2p  ,i_depsp    ,i_plastic
     6         ,i_porosity,i_fnno  ,i_elplastic,i_vpropel
     7         ,fstress0  ,fporomec,fmec       ,print_flag(1)
     8         ,fplastic  ,vprop   ,nin )
c    -----------------------------------------------------------------
c    | ix | id | ie | nload | eload | inum | e | x | f | u | u0 | tx0 |
c    -----------------------------------------------------------------
c
c    -----------------------------------------------------------------
c    | dp | tx1p | tx2p | epsp | plastic | porosity | fnno | elpastic | 
c    -----------------------------------------------------------------
c
c    -----------------------------------------------------------------
c    | vpropel |                                                
c    -----------------------------------------------------------------
c ......................................................................      
c
c ... calculo das forcas internas devidos as tensoes inicias tensoes
c     habilitado
      if(fstress0) fcstress0 = .true.
c ......................................................................   
c
c ... estrutura de dados para o pardiso 
      if(solver .eq. 10) then
        stge         = 6
        block_pu_sym = .false.
        block_pu     = .false.
      endif
c .....................................................................
c
c ... estrutura de dados para o bloclo iterativo pcg (poro_mec)
      if(solver .eq. 5) then
        fname = name(prename,nprcs,53)
        open(logsolvd,file=fname)
        block_pu_sym = .false.
        block_pu     = .true.
        n_blocks_pu  = 3
      endif
c .....................................................................
c 
c ... desabilita o csrc blocado em problemas mecanicos
      if(fmec) then
        block_pu        = .false. 
        block_pu_sym    = .false. 
        n_blocks_pu     = 0 
      endif
c ......................................................................   
c
c ... desabilita a impressao da pressao de consolidacao no problemas 
c     elasticos
      if(.not.fplastic) then
        print_flag(10) = .false.
        print_flag(11) = .false.
      endif  
c .....................................................................
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
      prebdiagtime = 0.d0
      totaltime    = MPI_Wtime()
c ......................................................................      
c
c.... Otimizacao da largura de banda:
c
      timei = MPI_Wtime()
      if (mpi) then
       call reord(ia(i_ix),ia(i_inum),nno1-nno1a,nnode,numel,nen,reordf)
      else
       call reord(ia(i_ix),ia(i_inum),nnode,nnode,numel,nen,reordf)
      endif
      reordtime = MPI_Wtime()-timei
c ......................................................................
c
c ... calcula o numero de equacoes
      if(fneqs) then 
        call numeq_pmec_e(ia(i_id),ia(i_fnno),nnode,ndf,my_id) 
        goto 5000 
      endif
c ......................................................................
c
c ... Numeracao nodal das equacoes:
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
        if(block_pu .or. block_pu_sym) then
          call numeqpmec1(ia(i_id),ia(i_inum),ia(i_id),ia(i_fnno),
     .                    nnode,nnodev,ndf,neq,nequ,neqp)
        else
c ... numera os deslocamento e pressoes ao mesmo tempo
          call numeqpmec2(ia(i_id),ia(i_inum),ia(i_id),ia(i_fnno),
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
     1                ,nno_pload
     2                ,nnovG,nnoG,nelG,nnodev,nnode
     3                ,numel,ovlp,novlp,nprcs,nviz1
     4                ,nviz2,i_rreqs,i_sreqs,i_rcvs0i,i_dspl0i
     5                ,i_fmap0i,mpi)
c
c.... Mapa de equacoes de fronteira:
c
      if (ndf  .gt. 0) call frontb(ndf,ia(i_id),neq,neq1
     1                ,neq2,neq3,neq4,neq1a,neqf1,neqf2,neq32
     2                ,neq_dot
     3                ,i_fmap,i_rcvs,i_dspl,i_xf
     4                ,'fmap     ','rcvs    ','dspl    ','xf      ',0)
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
        i_p0l = alloc_8('p0l     ',1  ,nenv)
        i_pl  = alloc_8('pl      ',1  ,nst)
c ... 6 - tensoes totais, 6 - tensoes de biot , 3 - fluxo de darcy
        i_txl      = alloc_8('txl     ', 15,nenv)
c ...      
        i_txnl     = alloc_8('txnl    ',  6,nen)
c ...
        i_tx1pl     = 1 
        i_tx2pl     = 1
        i_depsl     = 1
        i_plasticl  = 1        
        if(fplastic) then 
c ... tensoes nos pontos de integracao
          i_tx1pl = alloc_8('tx1l     ',ntn,npi)
          i_tx2pl = alloc_8('tx2l     ',ntn,npi)
c     delta deformacao e pressoes entre iteracoes nao lineares
          i_depsl = alloc_8('depsl    ',ntn+1,npi)
c ... dilatacao volumetrica plastica
c    , dilatacao volumetrica plastica do passo de tempo anterior
c     e paramentro de endurecimento
          i_plasticl = alloc_8('plasticl', 3,npi)
        endif
c ...
        i_vpropell  = 1 
        if(vprop(1)) then 
c ... propriedades variaveis nos pontos de integracao
          i_vpropell = alloc_8('vpropell ',nvprop,npi)
        endif
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
c ... -emoria para a estrutura de dados do sistema de equacoes:
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
     1                     ,numel   ,nen     ,ndf       ,nst
     2                     ,neq     ,nequ    ,neqp      ,stge ,unsym 
     3                     ,nad     ,naduu   ,nadpp   ,nadpu,nadr
     4                     ,i_ia ,i_ja,i_au  ,i_al,i_ad 
     5                     ,sia  ,sja ,sau   ,sal ,sad    
     6                     ,ovlp ,n_blocks_pu,block_pu ,block_pu_sym) 
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
     1                  ,numel   ,nen     ,ndf       ,nst
     2                  ,neq     ,stge    ,unsym     ,nad  ,nadr
     3                  ,i_ia    ,i_ja    ,i_au      ,i_al ,i_ad 
     4                  ,sia     ,sja     ,sau       ,sal  ,sad         
     5                  ,ovlp     )
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
c   vetor que usam a subrotina comunicate necessitam ter a dimensao        
c   neq+neq3+neq4                                                       
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
         if(precond .eq. 2 .or. precond .eq. 5 
     .     .or.  precond .eq. 7  ) then 
           i_m   = alloc_8('m       ',    1,neq)
           call azero(ia(i_m),neq)
c .....................................................................
c
c ...  Memoria para o precondicionador iLDLt e iLLT (cholesky)
         else if( precond .eq. 3 .or.  precond .eq. 4) then
           i_m   = alloc_8('m       ',    1,neq+nad)
           call azero(ia(i_m),neq+nad)
c ..................................................................... 
c
c ...  Memoria para o precondicionador block diagonal
         else if( precond .eq. 6) then
           i_m   = alloc_8('m       ',iparam(3),neq)
           call azero(ia(i_m),iparam(3)*neq)
c ..................................................................... 
         endif       
      endif
c ......................................................................
c
c ...
      print_nnode = nnovG   
      if(print_flag(1)) print_nnode = nnoG  
c ......................................................................
c
c ... porosidade nodal inicial
      i_ic        = alloc_4('ic      ',    1,nnode)
      call initial_porosity(ia(i_ix)      ,ia(i_e)  ,ia(i_ie),ia(i_ic) 
     1                     ,ia(i_porosity),ia(i_fnno)
     2                     ,nnode         ,numel     ,nen ,nenv
     3                     ,ndm           ,ndf       ,i_xf,novlp)     
      i_ic       = dealloc('ic      ') 
c ..................................................................... 
c
c ... 
      if( vprop(1) ) then                    
        call initial_prop(ia(i_ix),ia(i_e),ia(i_ie),ia(i_vpropel) 
     1                   ,numel   ,nen    ,numat   ,npi
     2                   ,vprop)
      endif
c ..................................................................... 
c
c ...
      if(fstress0 .and. fcstress0) then
        timei = MPI_Wtime()
        call initial_stress(ia(i_ix)     ,ia(i_ie)  ,ia(i_e)  
     1       ,ia(i_x)     ,ia(i_id)      ,ia(i_bst0) 
     2       ,ia(i_u0)    ,ia(i_tx1p)    ,ia(i_tx0) ,ia(i_dp)
     3       ,ia(i_xl)    ,ia(i_ul)      ,ia(i_dpl) ,ia(i_pl)
     4       ,ia(i_ld)    ,ia(i_txnl)    ,ia(i_tx1pl) 
     5       ,numel       ,nen           ,nenv     ,ndf 
     6       ,ndm         ,nst           ,npi      ,ntn
     7       ,neq         ,stge          ,ilib      
     8       ,block_pu    ,fplastic)
        elmtime = elmtime + MPI_Wtime()-timei
c ... inicializa as tensoes tx1p -> tx2p
        if(fplastic) then 
          call aequalb(ia(i_tx2p),ia(i_tx1p),ntn*npi*numel) 
        endif
        fcstress0= .false.
      endif 
c ..................................................................... 
c
c ... inicializa o paramentro de encruamento
      if(fplastic) call initial_pc0(ia(i_ix),ia(i_ie) ,ia(i_e)
     1                   ,ia(i_x) ,ia(i_tx0),ia(i_u)   ,ia(i_plastic)
     3                   ,ia(i_xl),ia(i_ul) ,ia(i_txnl),ia(i_plasticl)  
     5                   ,numel    ,nen     ,nenv     ,ndf 
     6                   ,ndm      ,nst     ,npi      ,ntn       
     7                   ,ilib    )
c .....................................................................
c
c ...
      go to 50
c ......................................................................
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
      t      = t + dt
c .....................................................................
c
c ...
      if(my_id.eq.0) then
        write(*,'(a,i8,a,f15.1,a,f15.5,a,f15.5)')
     .                                  ,' STEP '      ,istep
     .                                  ,' Time(s)'    ,t
     .                                  ,' Time(hours)',t/3600.d0
     .                                  ,' Time(days)' ,t/86400.d0
        write(nout_nonlinear,'(a,i9,e13.6)')'step',istep,t
        if(fhist_log)write(log_hist_solv,'(a,i9,e13.6)')'step',istep,t
      endif
c .....................................................................
c
c ... Cargas nodais e valores prescritos no tempo t+dt:
      timei = MPI_Wtime()
      call pload_pm(ia(i_id),ia(i_f),ia(i_u),ia(i_b0),ia(i_nload)
     .             ,ia(i_fnno),nnode,ndf)
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................   
c
c ... forcas de volume e superficie do tempo t+dt e graus de liberade 
c     do passo t:  
      timei = MPI_Wtime()
      call pform_pm(ia(i_ix)  ,ia(i_eload),ia(i_ie)  ,ia(i_e) 
     1           ,ia(i_x)     ,ia(i_id)   ,ia(i_ia)  ,ia(i_ja) 
     2           ,ia(i_au)    ,ia(i_al)   ,ia(i_ad)  ,ia(i_b0) 
     3           ,ia(i_u0)    ,ia(i_u)    ,ia(i_tx1p),ia(i_tx2p)
     4           ,ia(i_depsp) ,ia(i_dp)   ,ia(i_plastic),ia(i_elplastic)
     5           ,ia(i_vpropel)   
     6           ,ia(i_xl)    ,ia(i_ul)   ,ia(i_p0l)  ,ia(i_dpl)
     7           ,ia(i_pl)    ,ia(i_sl)   ,ia(i_ld)   ,ia(i_txnl)
     8           ,ia(i_tx1pl) ,ia(i_tx2pl),ia(i_depsl),ia(i_plasticl)
     9           ,ia(i_vpropell)   
     1           ,numel       ,nen        ,nenv       ,ndf 
     2           ,ndm         ,nst        ,npi        ,ntn
     3           ,neq         ,nequ       ,neqp 
     4           ,nad         ,naduu      ,nadpp      ,nadpu,nadr 
     5           ,.false.     ,.true.     ,unsym 
     6           ,stge        ,4          ,ilib     ,i
     7           ,ia(i_colorg),ia(i_elcolor),numcolors
     8           ,block_pu    ,n_blocks_pu  ,fplastic,vprop)
      elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ... tensao inicial elastico
      if(fstress0 .and. ( .not. fplastic) ) then
        call vsum(ia(i_b0),ia(i_bst0),neq,ia(i_b0))
      endif  
c .....................................................................
c
c ... Preditor: du(n+1,0) = u(n+1) - u(n)  
      timei = MPI_Wtime()
      call delta_predict_pm(nnode,ndf,ia(i_u),ia(i_u0),ia(i_fnno))
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ---------------------------------------------------------------------
c loop nao linear:
c ---------------------------------------------------------------------
      lhs = .true.  
  410 continue
c ...
      timei = MPI_Wtime()
      call aequalb(ia(i_b),ia(i_b0),neq)
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c      
c ... Residuo:
c ... plastic
c              Fu = Fu - int(BeT*sigma*dv)
c              bp = Fp - K.dp(n+1,i)
c ... elastic
c              Fu = Fu - K.du(n+1,i)
c              bp = Fp - K.dp(n+1,i)
      timei = MPI_Wtime()
      call pform_pm(ia(i_ix)    ,ia(i_eload) ,ia(i_ie)  ,ia(i_e)
     1         ,ia(i_x)     ,ia(i_id)    ,ia(i_ia)  ,ia(i_ja)
     2         ,ia(i_au)    ,ia(i_al)    ,ia(i_ad)  ,ia(i_b)
     3         ,ia(i_u0)    ,ia(i_u)     ,ia(i_tx1p),ia(i_tx2p)
     4         ,ia(i_depsp) ,ia(i_dp)    ,ia(i_plastic),ia(i_elplastic) 
     5         ,ia(i_vpropel)   
     6         ,ia(i_xl)    ,ia(i_ul)    ,ia(i_p0l)  ,ia(i_dpl)
     7         ,ia(i_pl)    ,ia(i_sl)    ,ia(i_ld)   ,ia(i_txnl)
     8         ,ia(i_tx1pl) ,ia(i_tx2pl) ,ia(i_depsl),ia(i_plasticl)
     9         ,ia(i_vpropell)   
     1         ,numel       ,nen         ,nenv     ,ndf
     2         ,ndm         ,nst         ,npi      ,ntn
     3         ,neq         ,nequ        ,neqp
     4         ,nad         ,naduu       ,nadpp    ,nadpu,nadr
     5         ,lhs         ,.true.      ,unsym
     6         ,stge        ,2           ,ilib     ,i
     7         ,ia(i_colorg),ia(i_elcolor),numcolors
     8         ,block_pu    ,n_blocks_pu  ,fplastic,vprop)
      elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ... Comunicacao do residuo para o caso non-overlapping:
      if (novlp) call communicate(ia(i_b),neqf1,neqf2,i_fmap,i_xf,
     .                            i_rcvs,i_dspl)
c ......................................................................
c
c ...
      if( stop_crit .eq. 1 .or. stop_crit .eq. 2) then
        call cal_residuo_pm(ia(i_id) ,ia(i_b)   ,ia(i_x0)
     1                  ,nno_pload,ndf       ,neq_dot,i     
     3                  ,istop    ,stop_crit ,tol
     4                  ,my_id    ,mpi       ,nout_nonlinear)
        if (istop .eq. 2) goto 420 
      endif
      if(fhist_log)write(log_hist_solv,'(a,i7)')'nlit',i       
c ......................................................................            
c
c ... solver (Kdu(n+1,i+1) = b; du(t+dt) )
      timei = MPI_Wtime()
      call solv_pm(neq  ,nequ    ,neqp  
     1         ,nad     ,naduu   ,nadpp      
     2         ,ia(i_ia),ia(i_ja),ia(i_ad)   ,ia(i_al)
     3         ,ia(i_m) ,ia(i_b) ,ia(i_x0)   ,solvtol,maxit
     4         ,ngram   ,block_pu,n_blocks_pu,solver,istep
     5         ,cmaxit  ,ctol    ,alfap      ,alfau ,precond
     6         ,fmec    ,fporomec,fhist_log  ,fprint 
     7         ,neqf1   ,neqf2   ,neq3       ,neq4  ,neq_dot
     8         ,i_fmap  ,i_xf    ,i_rcvs     ,i_dspl)
      soltime = soltime + MPI_Wtime()-timei
c .....................................................................
c
c ... atualizacao :      du(n+1,i+1) = du(n+1,i)      + dv(n+1,i+1)
      timei = MPI_Wtime()
      call delta_update_pm(nnode   ,ndf
     1                    ,ia(i_id),ia(i_u)
     2                    ,ia(i_x0),ia(i_fnno))
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ... atualizacao da propriedades nos pontos de integracao sem atualizar
c     as porosidades
      timei = MPI_Wtime()
      if(vprop(1) .and. newton_raphson) then
        call update_prop(ia(i_ix),ia(i_x),ia(i_e),ia(i_ie),ia(i_vpropel)
     1                  ,ia(i_u) ,ia(i_plastic)
     2                  ,ia(i_xl),ia(i_ul),ia(i_plasticl),ia(i_vpropell)
     3                  ,numel,nen ,nenv
     4                  ,ndm  ,ndf ,nst      ,npi 
     5                  ,10   ,ilib,,fplastic,vprop,.false.)
      endif
      upproptime = upproptime + MPI_Wtime()-timei
c ....................................................................
c
c ...
      if( stop_crit .eq. 3 .or. stop_crit .eq. 4) then
        call cal_residuo_pm(ia(i_id) ,ia(i_b)   ,ia(i_x0)
     1                  ,nno_pload,ndf       ,neq_dot,i     
     3                  ,istop    ,stop_crit ,tol
     4                  ,my_id    ,mpi       ,nout_nonlinear)
        if (istop .eq. 2) goto 420   
      endif
c .....................................................................
c
c ...
      if (i .ge. maxnlit)then
        if(my_id.eq.0) then
          print*,'Newton-Raphson: no convergence reached after '
     .          ,i,' iteretions !'      
        endif
        call stop_mef()
c       goto 420
      endif
      i = i + 1
c .....................................................................
c
c ... matriz K global calculada apenas na primeira iteracao
      if (.not. newton_raphson) lhs = .false.
c .....................................................................
      goto 410 
c .....................................................................
c
c .....................................................................
c fim do loop nao linear:
c .....................................................................
c
c ...
  420 continue
c ... atualizacao da propriedades nos pontos de integracao e as
c     porosidades
      timei = MPI_Wtime()
      if(vprop(1)) then
        call update_prop(ia(i_ix),ia(i_x),ia(i_e),ia(i_ie),ia(i_vpropel)
     1                  ,ia(i_u) ,ia(i_plastic)
     2                  ,ia(i_xl),ia(i_ul),ia(i_plasticl),ia(i_vpropell)
     3                  ,numel,nen ,nenv
     4                  ,ndm  ,ndf ,nst     ,npi 
     5                  ,10   ,ilib,fplastic,vprop,.true.)
      endif
      upproptime = upproptime + MPI_Wtime()-timei
c ....................................................................
c
c ... calculo da porosidade nodal
      i_ic    = alloc_4('ic      ',    1,nnode)
      i_dporo = alloc_8('dporo   ',    1,nnode)
      call porosity_form(ia(i_ix),ia(i_x) ,ia(i_e) ,ia(i_ie)
     1     ,ia(i_u) ,ia(i_porosity),ia(i_plastic),ia(i_dporo)
     2     ,ia(i_vpropel),ia(i_ic),ia(i_fnno)
     3     ,ia(i_xl),ia(i_ul),ia(i_pl),ia(i_plasticl),ia(i_vpropell)
     4     ,nnode   ,numel   ,nen    ,nenv
     5     ,ndm     ,ndf     ,nst    ,npi 
     6     ,7       ,ilib    ,i_xf   ,novlp
     7     ,fplastic,vprop)
      i_dporo = dealloc('dporo   ')
      i_ic    = dealloc('ic      ') 
c .....................................................................
c
c ... calculo da solucao no tempo t + dt e atualizacao dos valores do 
c     passo de tempo anterior: 
c     u(n+1)  = u(n) + du(n+1) 
c     u(n)    = u(n+1) 
c     dp(n+1) = p(n) - p(0)
      timei = MPI_Wtime()
c     call update_res(nnode  ,nnodev  ,ndf
c    .               ,ia(i_u),ia(i_u0),ia(i_dp),ia(i_pres0))
      call update_res_v2(nnode   ,ndf
     1                  ,ia(i_u) ,ia(i_u0)
     2                  ,ia(i_dp),ia(i_fnno))
c .....................................................................
c
c ... 
      if(fplastic) then
c ... atualizacoes as tensoes do passo de tempo anterior tx2p -> tx1p
        call aequalb(ia(i_tx1p),ia(i_tx2p),ntn*npi*numel)   
c ... atualizacoes as deformacoes volumetricas plastica 
c     do passo de tempo anterior 1 -> 2
        call update_plastic(ia(i_ix),ia(i_e),ia(i_plastic),ia(i_vpropel)
     .                     ,nen,npi,numel,vprop)
      endif
c .....................................................................
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
      print*,'Error reading macro (DT) !'
      goto 5000      
c ----------------------------------------------------------------------
c
c ... Macro-comando: PGEO - PRINT GEOMETRY
c
c ......................................................................
  600 continue
c ... numero do tensor de tensoes
c ... | sxx syy szz sxy  0 0 0|
      if( ndm .eq. 2) then
        ntn = 4
c ... | sxx syy szz  sxy syz sxz |
      else if(ndm .eq. 3) then
        ntn = 6
      endif
c .....................................................................
c
c ... Geometria:
      if(mpi) then       
        writetime = writetime + MPI_Wtime()-timei 
        call global_ix(nen+1,numel_nov,i_ix,i_g,'ixg     ')
        call global_v(ndm   ,nno_pload,i_x ,i_g1,'xg      ')
c .....................................................................
c
c ...
        if( my_id .eq. 0 ) then
          print*, 'Macro PGEO'
          call write_mesh_geo(ia(i_g)    ,ia(i_g1),print_nnode,nelG
     1                       ,nen        ,ndm     ,prename    ,bvtk
     2                       ,legacy_vtk ,nplot)
        endif        
c .....................................................................
c
c ...        
        i_g1 = dealloc('xg      ')
        i_g  = dealloc('ixg     ')
        writetime = writetime + MPI_Wtime()-timei
c ......................................................................
c
c ...
      else
        print*, 'Macro PGEO'
        writetime = writetime + MPI_Wtime()-timei 
        call write_mesh_geo_pm(ia(i_ix)   ,ia(i_x)    ,ia(i_ie)
     1                      ,ia(i_id)     ,ia(i_f)    ,ia(i_u) 
     2                      ,ia(i_tx0)    ,ia(i_nload),ia(i_eload)
     3                      ,print_nnode  ,numel      ,ndf     ,ntn
     4                      ,nen          ,ndm        ,prename
     5                      ,bvtk         ,macros     ,legacy_vtk
     6                      ,print_flag(1),nplot      ,nout_face)
        writetime = writetime + MPI_Wtime()-timei
      endif
c .....................................................................
      goto 50
c .....................................................................
c
c ... Macro-comando:           
c
c ......................................................................
  700 continue
c ...                                             
      if(my_id.eq.0) print*, 'Macro '
      goto 50
c .....................................................................
c
c ... Macro-comando: BLOCK_PU
c
c ......................................................................
  800 continue
      if(my_id.eq.0) print*, 'Macro Block_pu'
      if(flag_macro_mesh) then
        print*,'This macro can only be used before macro mesh'
        goto 5000
      endif
      block_pu = .true.
c ... n_blocks_pu
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =801,end =801) n_blocks_pu   
      goto 50
  801 continue
      print*,'Error reading macro (BLOCK_PU) n_blocks_pu   !'
      goto 5000
c ......................................................................
c
c ... Macro-comando: GRAVITY - aceleracao da gravidade
c
c ......................................................................
  900 continue
      if (my_id .eq. 0) print*, 'Macro Gravity'
      call read_gravity(gravity,gravity_mod,nin)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:       
c
c ......................................................................
 1000 continue
      if(my_id.eq.0)print*, 'Macro CONSEQ'
      if(flag_macro_mesh) then
        print*,'This macro can only be used before macro mesh'
        goto 5000
      endif
      call read_constitutive_equation(fplastic,vprop,my_id,nin)
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: SOLVER configuracao do solver
c
c ......................................................................
 1100 continue
      if(my_id.eq.0)print*, 'Macro SOLVER'
      if(flag_macro_mesh) then
        print*,'This macro can only be used before macro mesh'
        goto 5000
      endif
      call read_solver_config_pm(solver   ,solvtol
     1                          ,maxit    ,precond
     2                          ,ngram    ,fhist_log
     3                          ,fprint   
     4                          ,prename  ,log_hist_solv
     5                          ,alfap    ,alfau
     6                          ,ctol     ,cmaxit
     7                          ,nprcs    ,my_id,nin)
      goto 50
c ......................................................................
c
c ... Macro-comando: DELTATC - deltat critico
c
c ......................................................................
 1200 continue
      if(my_id.eq.0) print*, 'Macro DELTATC'
      if(fporomec) then
        call deltat_critico(ia(i_ix)   ,ia(i_eload),ia(i_ie),ia(i_e)
     1                     ,ia(i_x)    ,ia(i_xl)
     2                     ,numel      ,nen         ,nenv,ndf
     3                     ,ndm  ,nst  ,1           ,ilib
     4                     ,my_id,mpi)
      endif
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PCOO - impressao da matriz no formato COO
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
        i_lin  = alloc_4('lincoo  ',1,neq+nad) 
        i_col  = alloc_4('colcoo  ',1,neq+nad) 
        i_acoo = alloc_8('acoo    ',1,neq+nad)
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
     .                   ,ia(i_b  ),neq      ,neq+nad  
     .                   ,fname    ,nout   
     .                   ,.false.  ,.false.  )
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
c ... Macro-comando:            
c
c ......................................................................
 1400 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
 1500 continue
      goto 50 
c ----------------------------------------------------------------------
c
c ... Macro-comando: PRES - impressao dos resultados    
c
c ......................................................................
 1600 continue
      if(my_id.eq.0) print*,'Macro PRES'
c ... print_flag (true| false)
c     2  - desloc
c     3  - pressao 
c     4  - delta pressa 
c ... 5  - stress Total
c ... 6  - stress Biot
c ... 7  - stress Terzaghi
c ... 8  - fluxo de darcy 
c ... 9  - delta porosidade
c ...10  - pressao de consolidacao
c
c ... calculo da tensoes, tensoes efetivas e fluxo de darcy nos vertices.
c
c ...
      i_tx   = 1
      i_txe  = 1
      i_txb  = 1
      i_flux = 1
      i_ic   = 1
      i_g    = 1
      i_g1   = 1
      i_g2   = 1
      i_g3   = 1
      i_g4   = 1
      i_g5   = 1
      i_g6   = 1
      i_g7   = 1
      i_g8   = 1
      i_g9   = 1
      i_g10  = 1
c ......................................................................
c
c ...
      i_ic        = alloc_4('ic      ',    1,nnode)
      i_tx        = alloc_8('tx      ',  ntn,nnode)
      i_txe       = alloc_8('txe     ',  ntn,nnode)
      i_txb       = alloc_8('txb     ',  ntn,nnode)
      i_flux      = alloc_8('flux    ',  ndm,nnode)
      i_pc        = alloc_8('pc      ',    1,nnode)
c .....................................................................
c
c ... 
      if(mpi) then
c ...
        if(print_flag(5) .or. print_flag(6) .or. print_flag(7) 
     .     .or. print_flag(8)) then
          timei = MPI_Wtime()
          call tform_pm(ia(i_ix)   ,ia(i_x)  ,ia(i_e)  ,ia(i_ie)
     1        ,ia(i_ic)   ,ia(i_xl) ,ia(i_ul),ia(i_dpl),ia(i_tx1pl)
     2        ,ia(i_vpropell) 
     3        ,ia(i_txl)  ,ia(i_u)  ,ia(i_dp),ia(i_tx1p),ia(i_vpropel) 
     4        ,ia(i_tx)   ,ia(i_txb),ia(i_flux),ia(i_fnno) 
     5        ,nnode      ,numel   ,nen       ,nenv
     6        ,ndm        ,ndf     ,nst       ,ntn  ,npi 
     7        ,3          ,ilib    ,i_xf      ,novlp,fplastic,vprop)
          tformtime = tformtime + MPI_Wtime()-timei
        endif
c ......................................................................
c
c ...
        if(print_flag(10))then
          timei = MPI_Wtime()
          call consolidation_pressure(ia(i_ix),ia(i_ie),ia(i_ic)
     1                    ,ia(i_plasticl)    ,ia(i_pl) 
     2                    ,ia(i_pc)          ,ia(i_plastic),ia(i_fnno)  
     3                    ,nnode      ,numel  ,nen      ,nenv
     4                    ,ndm        ,ndf    ,nst      ,npi
     5                    ,8          ,ilib   ,i_xf     ,novlp)
          tformtime = tformtime + MPI_Wtime()-timei
        endif
c ......................................................................
c
c ... comunicao ( nno_pload = no1 + no2 )
        call global_v(ndf   ,nno_pload,i_u        ,i_g ,'upG     ')
        call global_ix(nen+1,numel_nov,i_ix       ,i_g1,'ixG     ')
        call global_v(ndm   ,nno_pload,i_x        ,i_g2,'xG      ')
        call global_v(1     ,nno_pload,i_dp       ,i_g3,'dpG     ')
        call global_v(ntn   ,nno_pload,i_tx0      ,i_g4,'tx0G    ')
c ... tensao e fluxo 
        if(print_flag(5) .or. print_flag(6) .or. print_flag(7) ) then
          if(novlp) then
            call global_v(ntn,nnode,i_tx ,i_g5 ,'txG     ')
            call global_v(ntn,nnode,i_txb,i_g6 ,'txbG    ')
            i_g7        = alloc_8('txeG    ',ntn,nnovG)
            call global_v(ndm,nnode,i_flux,i_g8,'fluxG   ')
          else          
            call global_v(ntn,nno_pload,i_tx ,i_g5,'txG     ')
            call global_v(ntn,nno_pload,i_txb,i_g6,'txbG    ')
            i_g7        = alloc_8('txeG    ',ntn,nnovG)
            call global_v(ndm,nno_pload,i_flux,i_g8,'fluxG   ')
          endif
c ... add tensao incial
          if(.not. fplastic) then
            call vsum(ia(i_g5),ia(i_g4),nnovG*ntn,ia(i_g5))
            call vsum(ia(i_g6),ia(i_g4),nnovG*ntn,ia(i_g6))
          endif
          call effective_stress(ia(i_g7),ia(i_g5),ia(i_g) 
     .                         ,nnovG   ,ntn     ,ndf)
c ......................................................................
        endif
c ......................................................................
c
c ... delta porosidade 
        if(print_flag(9)) then
          if(novlp) then
            call global_v(1,nnode,i_porosity,i_g9 ,'poroG   ')
          else          
            call global_v(1,nno_pload,i_porosity,i_g9 ,'poroG   ')
          endif
        endif
c ......................................................................
c
c ... pressao de consolidacao
        if(print_flag(10)) then
          if(novlp) then
            call global_v(1,nnode,i_pc,i_g10 ,'pcG     ')
          else          
            call global_v(1,nno_pload,i_pc,i_g10 ,'pcG     ')
          endif
        endif
c ......................................................................
c
c ... indica se o elemento plastificou ou nao ( 0 ou 1 )
        if(print_flag(11)) then
          call global_v_elm(1 ,numel_nov,i_elplastic,i_g11,'eplG    ',1)
        endif
c ......................................................................
c
c ...
        if(my_id .eq. 0) then
          call write_mesh_res_pm(ia(i_g1),ia(i_g2) ,ia(i_g)   ,ia(i_g3)
     1               ,ia(i_g10)        ,ia(i_g11)
     2               ,ia(i_g9)         ,ia(i_g5)   ,ia(i_g6)  ,ia(i_g7)
     3               ,ia(i_g8)         ,print_nnode,nelG      ,istep   
     4               ,t                ,nen        ,ndm       ,ndf   
     5               ,ntn              ,fname      ,prename   
     6               ,bvtk             ,legacy_vtk ,print_flag,nplot)
        endif
c ......................................................................
c
c ...
          if(print_flag(11)) i_g11 = dealloc('eplG    ')
c
          if(print_flag(10)) i_g10 = dealloc('pcG     ')
c   
          if(print_flag(9)) i_g9  = dealloc('poroG   ')
c
          if(print_flag(5) .or. print_flag(6) .or. print_flag(7) ) then
            i_g8  = dealloc('fluxG   ')
            i_g7  = dealloc('txeG    ')
            i_g6  = dealloc('txbG    ')
            i_g5  = dealloc('txG     ')
          endif
c
          i_g4  = dealloc('tx0G    ')
          i_g3  = dealloc('dpG     ')
          i_g2  = dealloc('xG      ')
          i_g1  = dealloc('ixG     ')
          i_g   = dealloc('upG     ')
c .......................................................................  
c
c ...
      else
c ...
        if( print_flag(5) .or. print_flag(6) .or. print_flag(7) 
     .     .or. print_flag(8)) then
          timei = MPI_Wtime()
          call tform_pm(ia(i_ix)   ,ia(i_x)    ,ia(i_e)  ,ia(i_ie)
     1       ,ia(i_ic)   ,ia(i_xl)  ,ia(i_ul)  ,ia(i_dpl),ia(i_tx1pl)
     2       ,ia(i_vpropell) 
     3       ,ia(i_txl)  ,ia(i_u)   ,ia(i_dp)  ,ia(i_tx1p),ia(i_vpropel)
     4       ,ia(i_tx)   ,ia(i_txb) ,ia(i_flux),ia(i_fnno) 
     5       ,nnode      ,numel     ,nen       ,nenv
     6       ,ndm        ,ndf       ,nst       ,ntn  ,npi 
     7       ,3          ,ilib      ,i_xf      ,novlp,fplastic,vprop)
          tformtime = tformtime + MPI_Wtime()-timei
c ......................................................................
c
c ... add tensao inicial
          if(.not. fplastic) then
            call vsum(ia(i_tx) ,ia(i_tx0),nnodev*ntn,ia(i_tx))
            call vsum(ia(i_txb),ia(i_tx0),nnodev*ntn,ia(i_txb))
          endif
          call effective_stress(ia(i_txe),ia(i_tx),ia(i_u) 
     .                         ,nnodev  ,ntn     ,ndf)
c ......................................................................
        endif
c ......................................................................
c
c ...
        if(print_flag(10))then
          timei = MPI_Wtime()
          call consolidation_pressure(ia(i_ix),ia(i_ie),ia(i_ic)
     1                    ,ia(i_plasticl)    ,ia(i_pl) 
     2                    ,ia(i_pc)          ,ia(i_plastic),ia(i_fnno)  
     3                    ,nnode      ,numel  ,nen      ,nenv
     4                    ,ndm        ,ndf    ,nst      ,npi
     5                    ,8          ,ilib   ,i_xf     ,novlp)
          tformtime = tformtime + MPI_Wtime()-timei
        endif
c ......................................................................
c
c ...
        call write_mesh_res_pm(ia(i_ix),ia(i_x)    ,ia(i_u)  ,ia(i_dp)
     1              ,ia(i_pc)       ,ia(i_elplastic)
     2              ,ia(i_porosity) ,ia(i_tx)   ,ia(i_txb),ia(i_txe)
     3              ,ia(i_flux)     ,print_nnode,numel    ,istep   
     4              ,t              ,nen        ,ndm      ,ndf   
     5              ,ntn            ,fname      ,prename   
     6              ,bvtk           ,legacy_vtk ,print_flag,nplot)
c ......................................................................
      endif
c ......................................................................
c
c ...
      i_pc       = dealloc('pc      ')
      i_flux     = dealloc('flux    ')
      i_txb      = dealloc('txb     ')
      i_txe      = dealloc('txe     ')
      i_tx       = dealloc('tx      ')
      i_ic       = dealloc('ic      ')
c ......................................................................
      goto 50     
c ......................................................................
c
c ... Macro-comando:       
c
c ......................................................................
 1700 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SOLVM
c
c ......................................................................
 1800 continue
      if (my_id .eq. 0 ) print*, 'Macro  SOLVM'
c ...
c     ilib   = 1
c     i      = 1
c     istep  = istep + 1
c     resid0 = 0.d0
c .....................................................................
c
c ... Cargas nodais e valores prescritos no tempo t+dt:
c     timei = MPI_Wtime()
c     call pload_mec(ia(i_id),ia(i_f),ia(i_u),ia(i_b0),ia(i_nload)
c    .              ,nnode   ,ndf)
c     vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ... forcas internas devidos as tensoes inicias tensoes
c     if(fstress0 .and. fcstress0) then
c         timei = MPI_Wtime()
c         call pform_mec(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e) 
c    1                 ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
c    2                 ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_bst0) 
c    3                 ,ia(i_u0)    ,ia(i_tx0)
c    4                 ,ia(i_xl)    ,ia(i_ul)     ,ia(i_pl)
c    5                 ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl) 
c    6                 ,numel       ,nen          ,nenv     ,ndf 
c    7                 ,ndm         ,nst          ,neq      ,nad,nadr   
c    8                 ,.false.     ,.true.       ,unsym 
c    9                 ,stge        ,5            ,ilib     ,i
c    1                 ,ia(i_colorg),ia(i_elcolor),numcolors,.true.)
c         elmtime = elmtime + MPI_Wtime()-timei
c         fcstress0= .false.
c     endif 
c .....................................................................      
c
c ... forcas de volume e superficie do tempo t+dt :  
c     timei = MPI_Wtime()
c     call pform_mec(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e) 
c    1              ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
c    2              ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_b0) 
c    3              ,ia(i_u0)    ,ia(i_tx0)
c    4              ,ia(i_xl)    ,ia(i_ul)     ,ia(i_pl)
c    5              ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl) 
c    6              ,numel       ,nen          ,nenv     ,ndf 
c    7              ,ndm         ,nst          ,neq      ,nad,nadr  
c    8              ,.false.     ,.true.       ,unsym 
c    9              ,stge,4      ,ilib         ,i
c    1              ,ia(i_colorg),ia(i_elcolor),numcolors,.false.)
c     elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ... tensao inicial
c     if(fstress0) then
c       call vsum(ia(i_b0),ia(i_bst0),neq,ia(i_b0))
c     endif  
c .....................................................................
c
c ---------------------------------------------------------------------
c loop nao linear:
c ---------------------------------------------------------------------
c1810 continue
c ... Loop multi-corretor:      
c     if(my_id.eq.0) print*,'nonlinear iteration ',i
c ...
c     timei = MPI_Wtime()
c     call aequalb(ia(i_b),ia(i_b0),neq)
c     vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ... Residuo: b = F - K.u(n+1,i)
c     timei = MPI_Wtime()
c     call pform_mec(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e)
c    1              ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
c    2              ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_b)
c    3              ,ia(i_u)     ,ia(i_tx0)
c    4              ,ia(i_xl)    ,ia(i_ul)     ,ia(i_pl)
c    5              ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl)
c    6              ,numel       ,nen          ,nenv     ,ndf
c    7              ,ndm         ,nst          ,neq      ,nad ,nadr    
c    8              ,.true.      ,.true.       ,unsym
c    9              ,stge        ,2            ,ilib     ,i
c    1              ,ia(i_colorg),ia(i_elcolor),numcolors,.false.)
c     elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ... Comunicacao do residuo para o caso non-overlapping:
c     if (novlp) call communicate(ia(i_b),neqf1,neqf2,i_fmap,i_xf,
c    .                            i_rcvs,i_dspl)
c ......................................................................
c
c ......................................................................
c     resid = dsqrt(dot_par(ia(i_b),ia(i_b),neq_dot))
c     if(i .eq. 1) resid0 = max(resid0,resid)
c     if(my_id .eq. 0 ) print*,'resid/resid0',resid/resid0,'resid',resid
c     if ((resid/resid0) .lt. tol) goto 1820     
c ......................................................................            
c
c ... solver (Ku(n+1,i+1) = b; u(t+dt) )
c     timei = MPI_Wtime()
c     call solv_pm(neq     ,nequ    ,neqp  
c    .            ,nad     ,naduu   ,nadpp      
c    .            ,ia(i_ia),ia(i_ja),ia(i_ad)    ,ia(i_al)
c    .            ,ia(i_m) ,ia(i_b) ,ia(i_x0)    ,solvtol,maxit
c    .            ,ngram   ,block_pu,n_blocks_pu ,solver ,istep
c    .            ,cmaxit  ,ctol    ,alfap       ,alfau  ,precond 
c    .            ,fmec    ,fporomec,fhist_log
c    .            ,neqf1   ,neqf2   ,neq3        ,neq4   ,neq_dot
c    .            ,i_fmap  ,i_xf    ,i_rcvs      ,i_dspl)
c     soltime = soltime + MPI_Wtime()-timei
c .....................................................................
c
c ... atualizacao :      u(n+1,i+1) = x(n+1,i)
c     timei = MPI_Wtime()
c     call update_mec(nnode,ndf,ia(i_id),ia(i_u),ia(i_x0))
c     vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ...
c     if (i .ge. maxnlit) goto 1820
c     i = i + 1
c     goto 1810 
c .....................................................................
c
c ---------------------------------------------------------------------
c fim do loop nao linear:
c ---------------------------------------------------------------------
c
c ...
c1820 continue 
c .....................................................................
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PMECRES
c
c ......................................................................
 1900 continue
      if(my_id.eq.0) print*,'Macro PMECRES'
c ... print_flag (true| false)
c     2  - desloc
c     10 - stress
c
c ... numero do tensor de tensoes
c ... | sxx syy szz sxy  0 0 0|
c     if( ndm .eq. 2) then
c       ntn = 4
c ... | sxx syy szz  sxy syz sxz |
c     else if(ndm .eq. 3) then
c       ntn = 6
c     endif
c ......................................................................
c
c ... 
c     i_tx = 1
c     i_g3 = 1
c     if(mpi) then
c ... comunicao
c       call global_v(ndf   ,nno_pload,i_u   ,i_g ,'dispG   ')
c       call global_ix(nen+1,numel_nov,i_ix  ,i_g1,'ixG     ')
c       call global_v(ndm   ,nno_pload,i_x   ,i_g2,'xG      ')
c       if(print_flag(10)) then
c         call global_v(ntn   ,nno_pload,i_tx0 ,i_g3,'tx0G    ')
c       endif
c ......................................................................
c
c ...
c       if(my_id.eq.0) then
c ...
c         if(print_flag(10)) then
c           i_tx  = alloc_8('tx      ',  ntn,print_nnode)
c           i_ic  = alloc_4('ic      ',    1,print_nnode)
c           call azero(ia(i_tx)    ,print_nnode*ntn)
c           call mzero(ia(i_ic)    ,print_nnode)
c .....................................................................
c
c ...
c           timei = MPI_Wtime()
c           call tform_mec(ia(i_g1) ,ia(i_g2),ia(i_e)   ,ia(i_ie)
c    .                   ,ia(i_ic)  ,ia(i_xl),ia(i_ul) 
c    .                   ,ia(i_txnl),ia(i_g) ,ia(i_g3),ia(i_tx) 
c    .                   ,nnovG     ,nelG    ,nenv      ,nen
c    .                   ,ndm       ,ndf     ,nst  ,ntn
c    .                   ,3         ,ilib)
c           tformtime = tformtime + MPI_Wtime()-timei
c ......................................................................
c         endif
c ......................................................................
c
c ...
c         call write_mesh_res_mec(ia(i_g1) ,ia(i_g2) ,ia(i_g),ia(i_tx)
c    .                         ,print_nnode,nelG
c    .                         ,nen        ,ndm ,ndf        ,ntn
c    .                         ,prename    ,istep
c    .                         ,bvtk       ,legacy_vtk,print_flag,nplot)
c ......................................................................
c
c ...
c         if(print_flag(10)) then
c           i_ic  = dealloc('ic      ')
c           i_tx  = dealloc('tx      ')
c         endif
c ......................................................................
c       endif
c ......................................................................
c
c ...
c       if(print_flag(10)) i_g3  = dealloc('tx0G    ')
c       i_g2  = dealloc('xG      ')
c       i_g1  = dealloc('ixG     ')
c       i_g   = dealloc('dispG   ')
c ......................................................................
c
c ...
c     else
c ...
c       if(print_flag(10)) then
c         i_tx  = alloc_8('tx      ',  ntn,print_nnode)
c         i_ic  = alloc_4('ic      ',    1,print_nnode)
c         call azero(ia(i_tx)    ,print_nnode*ntn)
c         call azero(ia(i_ic)    ,print_nnode)
c .....................................................................
c
c ...
c         timei = MPI_Wtime()
c         call tform_mec(ia(i_ix)  ,ia(i_x)  ,ia(i_e)  ,ia(i_ie)
c    1                  ,ia(i_ic)  ,ia(i_xl) ,ia(i_ul) 
c    2                  ,ia(i_txnl),ia(i_u)  ,ia(i_tx0),ia(i_tx) 
c    3                  ,nnodev    ,numel    ,nenv     ,nen
c    5                  ,ndm       ,ndf      ,nst      ,ntn
c    6                  ,3         ,ilib)
c         tformtime = tformtime + MPI_Wtime()-timei
c ......................................................................
c       endif
c ......................................................................
c
c ...
c       call write_mesh_res_mec(ia(i_ix) ,ia(i_x)    ,ia(i_u)  ,ia(i_tx)
c    1                       ,print_nnode,numel
c    2                       ,nen        ,ndm        ,ndf      ,ntn
c    3                       ,prename    ,istep
c    4                       ,bvtk       ,legacy_vtk,print_flag,nplot)
c ......................................................................
c
c ...
c       if(print_flag(10)) then
c         i_ic  = dealloc('ic      ')
c         i_tx  = dealloc('tx      ')
c       endif
c ......................................................................
c     endif
c ......................................................................
      goto 50     
c ......................................................................
c
c ... Macro-comando:
c
c ......................................................................
 2000 continue
      goto 50
c ......................................................................
c
c ... Macro-comando:
c
c ......................................................................
 2100 continue
      goto 50
c ......................................................................
c
c ... Macro-comando: 
c
c ......................................................................
 2200 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
 2300 continue
      goto 50 
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2400 continue
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
      if(my_id.eq.0)print*, 'Macro NL    '
c
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2610,end =2610) maxnlit
      if(my_id.eq.0)write(*,'(a,i10)')' Setting max nonlinear it for '
     .                      ,maxnlit
c
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2650,end =2650) tol
      if(my_id.eq.0) write(*,'(a,d10.2)')' Setting nonlinear tol for '
     .                                    ,tol
c
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2655,end =2655) stop_crit 
      if(my_id.eq.0) write(*,'(a,i2)')' Setting stoping criterion for '
     .              ,stop_crit
c
      goto 50
 2610 continue
      print*,'Error reading maxit in macro (NL) !'
      goto 5000
 2650 continue
      print*,'Error reading tol in macro (NL) !'
      goto 5000
 2655 continue
      print*,'Error reading stop_crit macro (NL) !'
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
c ... Macro-comando: 
c
c ......................................................................
 2800 continue
      goto 50
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
      if(my_id.eq.0) print*, 'Macro RETURN'
      nin   = naux 
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SETPNODE impressao de grandezas por no no tempo
c
c ......................................................................
 3100 continue
      if(my_id.eq.0) then
        print*, 'Macro SETPNODE'
        call readmacro(nin,.false.)
        write(pnodename,'(80a)') (word(i),i=1,80)
        goto 3120
c ... problema no arquivo auxiliar        
 3110   continue
        print*,'Error reading macro (SETPNODE)'
        flag_pnd = .false.
        goto 3130
c ... leitura normal 
 3120   continue     
        call readpnode(pnodename,i_no,i_nfile,num_pnode,flag_pnd,nout)
        new_file(1:nfiles) = .true.
 3130   continue
      endif
      call MPI_BCAST(flag_pnd,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c ... erro na letura do nome do arquivo auxiliar      
      if( flag_pnd .eqv. .false.) call stop_mef()
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SETPRINT impressao de grandezas no arquivo vtk  
c
c ......................................................................
 3200 continue
      if(my_id.eq.0) print*, 'Macro SETPRINT'
      call set_print_vtk_pm(print_flag,my_id,nin)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: impressao das grandesas no pontos de integracao   
c
c ......................................................................
 3300 continue
      if(my_id.eq.0) then
        print*, 'Macro SETPI'
        call readmacro(nin,.false.)
        write(str,'(80a)') (word(i),i=1,80)
        read(str,*,err=3310,end = 3310) ppiname
        goto 3320
c ... problema no arquivo auxiliar        
 3310   continue
        print*,'Error reading macro (SETPI)'
        flag_ppi = .false.
        goto 3330
c ... leitura normal 
 3320   continue     
        call readppi(ppiname,i_pi,i_nfile_pi,num_pel,flag_ppi,nout)
        new_file_pi(1:nfiles) = .true.
 3330   continue
      endif
      call MPI_BCAST(flag_ppi,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c ... erro na letura do nome do arquivo auxiliar      
      if( flag_ppi .eqv. .false.) call stop_mef()
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PPG impressao das variaveis por pontos de gauss
c     no tempo (SETPG)                                                   
c ......................................................................
 3400 continue
      if(my_id.eq.0) print*, 'Macro PGI'  
      if(.not. fplastic) goto 50    
      if(flag_pnd.eqv..false.) then
        if(my_id.eq.0)print*,'No cell to print for PGP!'  
        call stop_mef()
      endif
c ... codigo para o arquivo stress_pi.txt      
      code   = 35
      ifiles = 1
c .....................................................................
      call global_v_elm(ntn*npi,numel_nov,i_tx1p,i_g1,'stressG ',2)
      string = 'Stress'
      if( my_id .eq. 0) then
        do j = 1, num_pel  
          call printpi(ia(i_g1),ia(i_pi+j-1),ntn,npi     ,istep,t
     1                ,string  ,prename       ,ia(i_nfile_pi+j-1)
     2                ,code    ,new_file_pi(ifiles))
        enddo
        new_file_pi(ifiles) = .false.
      endif
      if(mpi) then
        i_g1 = dealloc('stressG ')
      endif
c .....................................................................
c
c ... codigo para o arquivo pc_pi.txt      
      code   = 36
      ifiles = 2
c .....................................................................
      call global_v_elm(3*npi,numel_nov,i_plastic,i_g1,'plasticG',2)
      string = 'Plastic'
      if( my_id .eq. 0) then
        do j = 1, num_pel  
          call printpi(ia(i_g1),ia(i_pi+j-1),  3,npi     ,istep,t
     1                ,string  ,prename       ,ia(i_nfile_pi+j-1)
     2                ,code    ,new_file_pi(ifiles))
        enddo
        new_file_pi(ifiles) = .false.
      endif
      if(mpi) then
        i_g1 = dealloc('plasticG')
      endif
c ......................................................................
c
c ...
      call MPI_barrier(MPI_COMM_WORLD,ierr)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PNDISP impressao do deslocamento por no no tempo
c     (SETPNODE)                                                   
c ......................................................................
 3500 continue
      if(my_id.eq.0) print*, 'Macro PNUP'      
      if(flag_pnd.eqv..false.) then
        if(my_id.eq.0)print*,'No node to print for PNUP!'  
        call stop_mef()
      endif
c ... codigo para o arquivo up_node.txt      
      code   = 30
      ifiles = 1
c .....................................................................
      call global_v(ndf,nno_pload,i_u,i_g1,'dispG   ')
      string = 'DeslocAndPress'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          call printnode(ia(i_g1),ia(i_no+j-1),ndf            ,istep,t
     1                  ,string  ,prename     ,ia(i_nfile+j-1)
     2                  ,code    ,new_file(ifiles))
        enddo
        new_file(ifiles) = .false.
      endif
      if(mpi) then
        i_g1      = dealloc('dispG   ')
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierr)
      goto 50
c ......................................................................
c
c ... Macro-comando:                                                    
c ......................................................................
 3600 continue
      if(my_id.eq.0) print*, 'Macro PNSF'
      if(flag_pnd.eqv..false.) then
        if(my_id.eq.0)print*,'No node to print for PNSF!'   
        call stop_mef()
      endif
c .....................................................................
c
c ... calculo da tensoes, tensoes efetivas e fluxo de darcy nos vertices.
      i_tx  = alloc_8('tx      ',  ntn,nnode)
      i_txe = alloc_8('txe     ',  ntn,nnode)
      i_txb = alloc_8('txb     ',  ntn,nnode)
      i_flux= alloc_8('flux    ',  ndm,nnode)
      i_ic  = alloc_4('ic      ',    1,nnode)
c .....................................................................
c
c ...
      timei = MPI_Wtime()
      call tform_pm(ia(i_ix)   ,ia(i_x)  ,ia(i_e)   ,ia(i_ie)
     1        ,ia(i_ic)   ,ia(i_xl) ,ia(i_ul)  ,ia(i_dpl),ia(i_tx1pl)
     2        ,ia(i_vpropell) 
     3        ,ia(i_txl)  ,ia(i_u)  ,ia(i_dp),ia(i_tx1p),ia(i_vpropel)
     4        ,ia(i_tx)   ,ia(i_txb),ia(i_flux),ia(i_fnno) 
     5        ,nnode      ,numel   ,nen        ,nenv
     6        ,ndm        ,ndf     ,nst        ,ntn   ,npi
     7        ,3          ,ilib    ,i_xf       ,novlp,fplastic,vprop)
      tformtime = tformtime + MPI_Wtime()-timei
c ......................................................................
c
c ... comunicao
      if(mpi) then
        call global_v(ndf   ,nno_pload,i_u   ,i_g ,'upG     ')
        call global_ix(nen+1,numel_nov,i_ix  ,i_g1,'ixG     ')
        call global_v(ndm   ,nno_pload,i_x   ,i_g2,'xG      ')
        call global_v(1     ,nno_pload,i_dp  ,i_g3,'dpG     ')
        call global_v(ntn   ,nno_pload,i_tx0 ,i_g4,'tx0G    ')
        if(novlp) then
          call global_v(ntn,nnode,i_tx ,i_g5 ,'txG     ')
          call global_v(ntn,nnode,i_txb,i_g6 ,'txbG    ')
          i_g7        = alloc_8('txeG    ',ntn,nnovG)
          call global_v(ndm,nnode,i_flux,i_g8,'fluxG   ')
        else          
          call global_v(ntn,nno_pload,i_tx ,i_g5,'txG     ')
          call global_v(ntn,nno_pload,i_txb,i_g6,'txbG    ')
          i_g7        = alloc_8('txeG    ',ntn,nnovG)
          call global_v(ndm,nno_pload,i_flux,i_g8,'fluxG   ')
        endif
c ... add tensao incial
        if(.not. fplastic) then
          call vsum(ia(i_g5),ia(i_g4),nnovG*ntn,ia(i_g5))
          call vsum(ia(i_g6),ia(i_g4),nnovG*ntn,ia(i_g6))
        endif
        call effective_stress(ia(i_g7),ia(i_g5),ia(i_g) 
     .                       ,nnovG   ,ntn     ,ndf)
      else
        if(.not. fplastic) then
          call vsum(ia(i_tx) ,ia(i_tx0),nnovG*ntn,ia(i_tx))
          call vsum(ia(i_txb),ia(i_tx0),nnovG*ntn,ia(i_txb))
        endif
        call effective_stress(ia(i_txe),ia(i_tx),ia(i_u) 
     .                       ,nnovG    ,ntn     ,ndf)      
      endif
c .....................................................................
c
c
c ... codigo para o arquivo stress_node.txt      
      code   = 31
      ifiles = 2
      string = 'totalStress'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          if(mpi) then
           call printnode(ia(i_g5),ia(i_no+j-1),ntn          ,istep,t
     1                   ,string   ,prename     ,ia(i_nfile+j-1)
     2                   ,code     ,new_file(ifiles))
          else
            call printnode(ia(i_tx),ia(i_no+j-1),ntn          ,istep,t
     1                   ,string   ,prename     ,ia(i_nfile+j-1)
     2                   ,code     ,new_file(ifiles))
          endif  
        enddo
        new_file(ifiles) = .false.
      endif
c ......................................................................
c
c ... codigo para o arquivo stressE_node.txt      
      code   = 32
      ifiles = 3
      string = 'eSstress'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
         if(mpi) then
           call printnode(ia(i_g7),ia(i_no+j-1),ntn          ,istep,t
     1                 ,string   ,prename     ,ia(i_nfile+j-1)
     2                 ,code     ,new_file(ifiles))
         else
           call printnode(ia(i_txe),ia(i_no+j-1),ntn          ,istep,t
     1                 ,string   ,prename     ,ia(i_nfile+j-1)
     2                 ,code     ,new_file(ifiles))
          endif
        enddo
        new_file(ifiles) = .false.
      endif
c ......................................................................
c
c ... codigo para o arquivo stressB_node.txt      
      code   = 33
      ifiles = 4
      string = 'biotStress'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          if(mpi) then
            call printnode(ia(i_g6),ia(i_no+j-1),ntn         ,istep,t
     1                   ,string   ,prename     ,ia(i_nfile+j-1)
     2                   ,code     ,new_file(ifiles))
          else
            call printnode(ia(i_txb),ia(i_no+j-1),ntn         ,istep,t
     1                   ,string   ,prename     ,ia(i_nfile+j-1)
     2                   ,code     ,new_file(ifiles))
          endif
        enddo
        new_file(ifiles) = .false.
      endif
c ......................................................................
c
c ... codigo para o arquivo flux_node.txt      
      code   = 34
      ifiles = 5
      string = 'darcyFlux'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          if(mpi) then
            call printnode(ia(i_g8),ia(i_no+j-1),ndm         ,istep,t
     1                     ,string   ,prename     ,ia(i_nfile+j-1)
     2                     ,code     ,new_file(ifiles))
          else
            call printnode(ia(i_flux),ia(i_no+j-1),ndm        ,istep,t
     1                 ,string   ,prename     ,ia(i_nfile+j-1)
     2                 ,code     ,new_file(ifiles))
          endif
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
c
c ...
      if(mpi) then
        i_g8      = dealloc('fluxG   ')
        i_g7      = dealloc('txeG    ')
        i_g6      = dealloc('txbG    ')
        i_g5      = dealloc('txG     ')
        i_g4      = dealloc('tx0G    ')
        i_g3      = dealloc('dpG     ')
        i_g2      = dealloc('xG      ')
        i_g1      = dealloc('ixG     ')
        i_g       = dealloc('upG     ')
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierr)
      goto 50
c ......................................................................
c
c ... Macro-comando: CONFIG
c
c ......................................................................
 3700 continue
      if(my_id .eq. 0 )print*, 'Macro CONFIG    '
      if(flag_macro_mesh) then
        print*,'This macro can only be used before macro mesh'
        goto 5000
      endif
      call read_config(maxmem
     1                ,omp_elmt ,omp_solv
     2                ,nth_elmt ,nth_solv
     3                ,reordf   ,newton_raphson
     4                ,bvtk     ,legacy_vtk 
     5                ,mpi      ,nprcs
     6                ,my_id    ,nin)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: NEQS
c
c ......................................................................
 3800 continue
      if(my_id .eq. 0 )print*, 'Macro NEQS'
      fneqs = .true.
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: run
c
c ......................................................................
 3900 continue
      if(my_id .eq. 0 )print*, 'Macro RUN'
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      naux = nin
      open(nincl, file= fname,status= 'old',err=3951,action='read')
      nin = nincl  
      goto 50
 3951 continue
      print*,'File ',trim(fname),' not found !'
      call stop_mef()
c ----------------------------------------------------------------------
c
c ... Macro-comando STOP:
c
c ......................................................................
 5000 continue
      if(my_id.eq.0)print*, 'Macro STOP'
c ... fecha o arquivo de entrada de dados
      close(nin)
c ... fecha log do solv e do loop nao linear
      close(logsolv)
      close(nout_nonlinear)
      if(fhist_log) close(log_hist_solv)
c ... fecha arquivo extra dp log do solv block_pcg_it
      if(solver .eq. 5 ) close(logsolvd)
c .....................................................................
c
c ...
      totaltime = MPI_Wtime() - totaltime
c
c ... arquivo de tempo      
      call write_log_file(nnode    ,numel   ,numel_nov,numel_ov,ndf 
     1                   ,neq      ,nequ    ,neqp     ,neq1    ,neq2
     2                   ,neq32    ,neq4    ,neq1a    ,neqf1   ,neqf2 
     3                   ,nad      ,naduu   ,nadpp    ,nadpu   ,nadr
     4                   ,omp_elmt ,nth_elmt,omp_solv ,nth_solv
     5                   ,fporomec ,fmec    ,numcolors,prename
     6                   ,my_id    ,nprcs   ,nout)
c .....................................................................
c
c ... media do tempo mpi 
      if(mpi) then    
        call mpi_log_mean_time(nnovG,nnoG,nelG
     1                        ,omp_elmt ,nth_elmt
     2                        ,omp_solv ,nth_solv
     3                        ,fporomec ,fmec 
     4                        ,numcolors,prename
     5                        ,my_id    ,nprcs   ,nout)
      endif
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

