c **********************************************************************
c *                                                                    *
c *   SOLV_PM.F                                            06/12/2015  *
c *                                                                    *
c *   Metodos iterativos de solucao:                                   *
c *                                                                    *
c *   pcg                                                              *
c *   gmres                                                            *
c *   pbcgstab                                                         *
c *   bi_pcg_it                                                        *
c *                                                                    *
c **********************************************************************
      subroutine solv_pm(neq    ,nequ    ,neqp
     .                  ,nad    ,naduu   ,nadpp
     .                  ,ip     ,ja      ,jat        ,iat    ,kat     
     .                  ,ad     ,al 
     .                  ,m      ,b       ,x          ,tol    ,maxit
     .                  ,ngram  ,block_pu,n_blocks_up,solver ,istep
     .                  ,cmaxit ,ctol    ,alfap      ,alfau  ,precond 
     .                  ,neqf1i ,neqf2i  ,neq3i      ,neq4i  ,neq_doti
     .                  ,i_fmapi,i_xfi   ,i_rcvsi    ,i_dspli)
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'openmp.fi'
      include 'time.fi'
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      integer*8 i_z,i_r,i_s,i_c,i_h,i_g,i_y,i_a
      integer*8 i_ilu,i_jat,i_iat,i_kat
c ...
      integer*8 i_afull
c ......................................................................
      integer neq3i,neq4i,neq_doti
      integer ip(*),ja(*),neq,nequ,neqp,nad,naduu,nadpp
      integer maxit,solver,ngram,istep,n_blocks_up
      real*8  ad(*),al(*),m(*),x(*),b(*),tol,energy
c ... pcg duplo
      integer cmaxit
      real*8  ctol,alfap,alfau
c ......................................................................
      logical block_pu
c ... precondicionador
      logical diag
      integer precond
      integer jat(*),iat(*),kat(*)
c ...................................................................... 
      integer neqovlp
      external dot,dot_par
      external matvec_csrc_pm,matvec_csrc_sym_pm        
c     OpenMP'ed subroutines
      external dot_par_omp,dot_par_omp_loopwise
      external matvec_csrc_sym_pm_omp,matvec_csrc_pm_omp
c ......................................................................
c ... numero total de equacoes na particao overlapping:
c    (neqovlp = neq, no sequencial e no non-overlapping)
      neqovlp = neq+neq3i+neq4i
      if (omp_solv) then
c ... matriz aramazenada em csrc blocado (Kuu,Kpp,Kpu)
         if(block_pu) then
           pmatrixtime = Mpi_Wtime() - pmatrixtime 
           i_threads_y = alloc_8('buffer_y',nth_solv,neq)
           call partition_matrix(ip,ja,neq,nequ,nad,ovlp,block_pu)
           pmatrixtime = Mpi_Wtime() - pmatrixtime
c .......................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
         else   
           pmatrixtime = Mpi_Wtime() - pmatrixtime 
           i_threads_y = alloc_8('buffer_y',nth_solv,neq)
           call partition_matrix(ip,ja,neq,nequ,nad,ovlp,block_pu)
           pmatrixtime = Mpi_Wtime() - pmatrixtime
         endif
      endif
c ......................................................................
c
c ... PCG
      if (solver .eq. 1) then
c ...
         i_z = alloc_8('zsolver ',1,neq)
         i_r = alloc_8('rsolver ',1,neq)
c ......................................................................
c
c ...    precondicionador diagonal:
         if(precond .eq. 2) then 
           call pre_diag(m,ad,neq)
c .....................................................................
         else if(precond .eq. 3) then
c ...
           precondtime = Mpi_Wtime() - precondtime 
           call ildlt(neq,ip,ja,al,ad,m,ia(i_z),0.0d0,.false.)
           precondtime = Mpi_Wtime() - precondtime 
c ...
         endif           
c ...    Comunicacao da diagonal para o caso non-overlapping:
         if (novlp) call communicate(m,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                               i_rcvsi,i_dspli)
c .....................................................................
c

c
c ... matriz aramazena em csrc blocado (Kuu,Kpp,Kpu)
         if(block_pu) then
           print*,"PCG não disponivel para a matriz blocada"
           stop   
c .....................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
         else
c ... omp
           if(omp_solv) then
             call pcg_omp(neq    ,nequ   ,nad,ip   ,ja
     .                   ,ad     ,al     ,al ,m    ,b ,x
     .                   ,ia(i_z),ia(i_r),tol,maxit
c ... matvec comum:
     .                   ,matvec_csrc_sym_pm_omp,dot_par_omp 
     .                   ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                   ,i_xfi ,i_rcvsi,i_dspli,ia(i_threads_y)
     .                   ,.true.)
c .....................................................................
c
c ... sequencial
           else
             call call_cg(neq      ,nequ  ,nad   ,ip      ,ja
     .                   ,ad       ,al    ,m    
     .                   ,jat      ,iat   ,kat     
     .                   ,b        ,x     ,ia(i_z)  ,ia(i_r) 
     .                   ,tol      ,maxit ,precond
     .                   ,my_id    ,neqf1i,neqf2i,neq_doti,i_fmapi
     .                   ,i_xfi ,i_rcvsi,i_dspli)  
c .....................................................................
c
           endif      
c .....................................................................
         endif
c .....................................................................
c
c ...
         i_r = dealloc('rsolver ')
         i_z = dealloc('zsolver ')
c ......................................................................
c
c ... Gmres com precondicionador diagonal:
      elseif(solver .eq. 2) then
         diag = .true.
c ...    precondicionador diagonal:
         if(diag) then 
           call pre_diag(m,ad,neq)
c ...    
         else                           
           m(1:neq) = 1.d0  
         endif           
c ...    Comunicacao da diagonal para o caso non-overlapping:
         if (novlp) call communicate(m,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                               i_rcvsi,i_dspli)
c .....................................................................
c
c ...
         i_g = alloc_8('gsolver ',neq    ,ngram+1)         
         i_h = alloc_8('hsolver ',ngram+1,ngram)
         i_y = alloc_8('ysolver ',1,ngram)
         i_c = alloc_8('csolver ',1,ngram)
         i_s = alloc_8('ssolver ',1,ngram)
         i_r = alloc_8('rsolver ',1,ngram+1)
c .....................................................................
c
c ... matriz aramazena em csrc blocado (Kuu,Kpp,Kpu)
         if(block_pu) then
c ... omp
           if(omp_solv) then
             call gmres_omp(neq,nequ,nad,ip,ja 
     .                 ,ad ,al  ,al ,m ,b ,x,ngram,ia(i_g) 
     .                 ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r) 
     .                 ,tol    ,maxit   
c ... matvec comum:
     .                 ,matvec_csrc_pm_omp,dot_par_omp 
     .                 ,neqovlp 
     .                 ,my_id  ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     .                 ,i_xfi  ,i_rcvsi,i_dspli  
     .                 ,ia(i_threads_y),.true.)
c .....................................................................
c
c ... sequencial
           else 
             call gmres(neq,nequ,nad,ip,ja 
     .                 ,ad ,al  ,al ,m ,b ,x,ngram,ia(i_g) 
     .                 ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r) 
     .                 ,tol    ,maxit   
c ... matvec comum:
     .                 ,matvec_csrc_pm  ,dot_par 
     .                 ,neqovlp 
     .                 ,my_id   ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     .                 ,i_xfi   ,i_rcvsi,i_dspli  
     .                 ,.true.)
           endif
c .....................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
         else
c ... omp
           if(omp_solv) then
             call gmres_omp(neq,nequ,nad,ip,ja 
     .                 ,ad ,al  ,al ,m ,b ,x,ngram,ia(i_g) 
     .                 ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r) 
     .                 ,tol    ,maxit   
c ... matvec comum:
     .                 ,matvec_csrc_sym_pm_omp,dot_par_omp 
     .                 ,neqovlp 
     .                 ,my_id   ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                 ,i_xfi   ,i_rcvsi,i_dspli 
     .                 ,ia(i_threads_y),.true.)
c .....................................................................
c
c ... sequencial
           else 
             call gmres(neq,nequ,nad,ip,ja 
     .                 ,ad ,al  ,al ,m ,b ,x,ngram,ia(i_g) 
     .                 ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r) 
     .                 ,tol    ,maxit   
c ... matvec comum:
     .                 ,matvec_csrc_sym_pm,dot_par 
     .                 ,neqovlp
     .                 ,my_id  ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     .                 ,i_xfi  ,i_rcvsi,i_dspli 
     .                 ,.true.)
           endif 
c .....................................................................
         endif
c .....................................................................
c
c ......................................................................
         i_r = dealloc('rsolver ')
         i_s = dealloc('ssolver ')
         i_c = dealloc('csolver ')
         i_y = dealloc('ysolver ')
         i_h = dealloc('hsolver ')
         i_g = dealloc('gsolver ')
c ......................................................................
c
c ...                                         
      elseif(solver .eq. 3) then
        print*,'Solver GUASS nao disponivel para o poromecanico !!'
c ......................................................................
c
c ... BICGSTAB com precondicionador diagonal:
      else if (solver .eq. 4) then
c ...
        if(n_blocks_up .eq. 1 ) then
          print*,"BICGSTAB não disponivel para a matriz",
     .           " blocada [Kuu Kpp]."
          stop   
        else if( n_blocks_up .eq. 3 ) then
          print*,"BICGSTAB não disponivel para a matriz",
     .           " blocada [Kuu] [Kpp] [Kpu]."
          stop   
        endif
c ...
         diag = .true.
c ...    precondicionador diagonal:
         if(diag) then 
           call pre_diag(m,ad,neq)
c ...    
         else                           
           m(1:neq) = 1.d0  
         endif           
c ...    Comunicacao da diagonal para o caso non-overlapping:
         if (novlp) call communicate(m,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                               i_rcvsi,i_dspli)
c .....................................................................
c
c ... alocacao dos arronjos auxiliares (5neq)    
         i_c = alloc_8('tsolver ',1,neq)
         i_h = alloc_8('hsolver ',1,neq)
         i_r = alloc_8('rsolver ',1,neq)
         i_s = alloc_8('psolver ',1,neq)
         i_z = alloc_8('zsolver ',1,neq)
c .....................................................................
c
c ... matriz aramazena em csrc blocado (Kuu,Kpp,Kpu)
         if(block_pu) then
c ... omp
           if(omp_solv) then
             call pbicgstab_omp(neq   ,nequ,nad    ,ip   ,ja      
     .                     ,ad     ,al     ,al     ,m      ,b      ,x  
     .                     ,ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z)
     .                     ,tol    ,maxit 
c ... matvec comum:
     .                     ,matvec_csrc_pm_omp,dot_par_omp 
     .                     ,my_id  ,neqf1i  ,neqf2i  ,neq_doti   
     .                     ,i_fmapi,i_xfi  ,i_rcvsi,i_dspli
     .                     ,ia(i_threads_y),.true.)
c .....................................................................
c
c ... sequencial
           else 
             call pbicgstab(neq   ,nequ   ,nad  ,ip     ,ja      
     .                    ,ad     ,al     ,al   ,m      ,b      ,x  
     .                    ,ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z)
     .                    ,tol    ,maxit 
c ... matvec comum:
     .                    ,matvec_csrc_pm ,dot_par 
     .                    ,my_id    ,neqf1i,neqf2i ,neq_doti 
     .                    ,i_fmapi ,i_xfi  ,i_rcvsi,i_dspli
     .                    ,.true.)
           endif 
c .....................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
         else
c ... omp
           if(omp_solv) then
             call pbicgstab_omp(neq   ,nequ,nad    ,ip   ,ja      
     .                     ,ad     ,al     ,al     ,m      ,b      ,x  
     .                     ,ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z)
     .                     ,tol    ,maxit 
c ... matvec comum:
     .                     ,matvec_csrc_sym_pm_omp,dot_par_omp 
     .                     ,my_id  ,neqf1i  ,neqf2i  ,neq_doti   
     .                     ,i_fmapi,i_xfi  ,i_rcvsi,i_dspli
     .                     ,ia(i_threads_y),.true.)
c .....................................................................
c
c ... sequencial
           else 
             call pbicgstab(neq   ,nequ   ,nad    ,ip   ,ja      
     .                     ,ad    ,al     ,al     ,m    ,b      ,x  
     .                     ,ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z)
     .                     ,tol    ,maxit 
c ... matvec comum:
     .                     ,matvec_csrc_sym_pm,dot_par 
     .                     ,my_id        ,neqf1i  ,neqf2i  ,neq_doti   
     .                     ,i_fmapi      ,i_xfi  ,i_rcvsi,i_dspli
     .                     ,.true.)
           endif      
c .....................................................................
         endif
c .....................................................................
c
c ... 
         i_z = dealloc('zsolver ')     
         i_s = dealloc('psolver ')
         i_r = dealloc('rsolver ')
         i_h = dealloc('hsolver ')
         i_c = dealloc('tsolver ')
c ......................................................................         
c
c ...                                                                     
      else if (solver .eq. 5) then
         diag = .true.
c ...    precondicionador diagonal:
         if(diag) then 
           call pre_diag(m        ,ad        ,nequ)
           call pre_diag(m(nequ+1),ad(nequ+1),neqp)
c ...    
         else                           
           m(1:neq) = 1.d0  
         endif           
c ......................................................................
c
c ... 
         if(block_pu) then
            if(n_blocks_up .eq. 1 ) then
              print*,"PCG_BLOCK indisponivel para a matriz",
     .               " blocada [Kuu Kpp]."
              stop   
            else if( n_blocks_up .eq. 2 ) then
              print*,"PCG_BLOCK indisponivel para a matriz",
     .               " blocada [Kuu Kpp] [Kpu]."
              stop   
            endif
c ......................................................................
c
c ...
           i_z  = alloc_8('zsolver ',1,nequ)
           i_r  = alloc_8('rsolver ',1,nequ)
           i_s  = alloc_8('ssolver ',1,nequ)
           i_c  = alloc_8('csolver ',1,neqp)
           i_h  = alloc_8('hsolver ',1,nequ)
           i_g  = alloc_8('gsolver ',1,neqp)
           i_y  = alloc_8('ysolver ',1,nequ)
           i_a  = alloc_8('asolver ',1,neqp)
c ......................................................................
c
c ...
            call pcg_block_it(neq    ,nequ       ,neqp  
     .                    ,nad       ,naduu      ,nadpp
     .                    ,ip        ,ja     
     .                    ,ip(nequ+2),ja(naduu+1)
     .                    ,ip(neq+3) ,ja(naduu+nadpp+1)         
     .                    ,ad        ,ad(nequ+1)
     .                    ,al        ,al(naduu+1),al(naduu+nadpp+1)  
     .                    ,m         ,m(nequ+1)  ,b               ,x
     .                    ,ia(i_z),ia(i_r),ia(i_s)  ,ia(i_c)
     .                    ,ia(i_h),ia(i_g),ia(i_y)  ,ia(i_a)
     .                    ,tol    ,ctol   ,maxit    ,cmaxit
     .                    ,alfap  ,alfau 
     .                    ,.true.,istep
     .                    ,my_id  ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                    ,i_xfi  ,i_rcvsi,i_dspli) 
c ......................................................................
c
c ... 
         else 
           print*,"PCG_BLOCK não disponivel para a matriz nao blocada"
           stop   
         endif 
c ......................................................................
c
c ...
         i_a  = dealloc('asolver ')        
         i_y  = dealloc('ysolver ')        
         i_g  = dealloc('gsolver ')        
         i_h  = dealloc('hsolver ')        
         i_c  = dealloc('csolver ')        
         i_s  = dealloc('ssolver ')         
         i_r  = dealloc('rsolver ')        
         i_z  = dealloc('zsolver ') 
c ......................................................................
      endif
c ......................................................................         
c
c ...                                                                     
      if (omp_solv) then
        pmatrixtime = Mpi_Wtime() - pmatrixtime 
        i_threads_y = dealloc('buffer_y')
        pmatrixtime = Mpi_Wtime() - pmatrixtime
      endif
c ......................................................................         
c
c ...   
      print*,precondtime,ifatsolvtime,time                            
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 11/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_CG : chama a versao do gradiente conjudado desejado           *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * al(*)    - parte triangular inferior de A                          *
c * b(neq)   - vetor de forcas                                         *
c * m(*)     - precondicionador                                        *
c * jat(*)   - arranjo auxiliar para precondicionador iLDLt            *
c * iat(*)   - arranjo auxiliar para precondicionador iLDLt            *
c * kat(*)   - arranjo auxiliar para precondicionador iLDLt            *
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * precond  - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diag                                                *
c *            3 - iLDLt                                               *
c *            4 -                                                     *
c *            5 -                                                     *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Arranjos jat,iat e kat são utilizados na retrosubstituizao do      *
c * solver iLDLt                                                       *
c **********************************************************************  
      subroutine call_cg(neq      ,nequ  ,nad   ,ia      ,ja
     .                  ,ad       ,al    ,m    
     .                  ,jat      ,iat   ,kat     
     .                  ,b        ,x     ,z     ,r     
     ,                  ,tol      ,maxit ,precond
     .                  ,my_id    ,neqf1i,neqf2i,neq_doti,i_fmapi
     .                  ,i_xfi    ,i_rcvsi,i_dspli)
      implicit none
      include 'time.fi'
c ... mpi
      integer my_id
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c .....................................................................
      integer neq,nequ,nad,neq_doti 
      integer ia(*),ja(*)
      real*8  ad(*),al(*),m(*),x(*),b(*),z(*),r(*)
c ...
      real*8  tol
      integer maxit  
c ... precondicionador
      logical diag
      integer precond
      integer jat(*),iat(*),kat(*)
c ...
      external dot_par
      external matvec_csrc_sym_pm        
c ......................................................................

c ... cg
      if(precond .eq. 1) then
c ...  
        call cg(neq    ,nequ   ,nad,ia   ,ja
     .         ,ad     ,al     ,al ,b    ,x
     .         ,z      ,r      ,tol,maxit
c ... matvec comum:
     .         ,matvec_csrc_sym_pm,dot_par 
     .         ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .         ,i_xfi ,i_rcvsi,i_dspli
     .         ,.true.,.true.)
c .....................................................................
c
c ... pcg - cg com precondicionador diagonal
      else if(precond .eq. 2) then
c ...  
        call pcg(neq    ,nequ   ,nad,ia   ,ja
     .          ,ad     ,al     ,al ,m    ,b ,x
     .          ,z      ,r      ,tol,maxit
c ... matvec comum:
     .          ,matvec_csrc_sym_pm,dot_par 
     .          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .          ,i_xfi ,i_rcvsi,i_dspli
     .          ,.true.,.true.)
c .....................................................................
c
c ...
      elseif(precond .eq. 3 ) then
        call icg(neq      ,nequ  ,nad   ,ia      ,ja
     .          ,ad       ,al    ,al    
     .          ,m    
     .          ,jat      ,iat   ,kat     
     .          ,b        ,x    
     .          ,z        ,r     ,tol      ,maxit
c ... matvec comum:
     .          ,matvec_csrc_sym_pm,dot_par 
     .          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .          ,i_xfi ,i_rcvsi,i_dspli
     .          ,.true.)
      endif  
c .....................................................................
      return
      end    
c ************************************************************************
c
c
c *************************************************************************
      subroutine get_res(u,x,id,nnode,nnodev,ndf)
      implicit none
      integer nnode,nnodev,ndf
      integer id(ndf,*),i,j,k
      real*8 u(ndf,*),x(*)
c ... loop nos
      do i = 1, nnode
        do j = 1, ndf - 1
          k    = id(j,i)
          if( k .gt. 0 ) then
            x(k) = u(j,i)
          endif
        enddo
      enddo
c .....................................................................
c
c ... loop nos vertices
      do i = 1, nnodev
        k    = id(ndf,i)
        if( k .gt. 0 ) then
          x(k) = u(ndf,i)
        endif
      enddo
c .....................................................................
c
c ...
      return
      end
c *************************************************************************
c
c **********************************************************************
c *                                                                    *
c *   SOLV.F                                             31/08/2005    *
c *                                                                    *
c *   Metodos iterativos de solucao:                                   *
c *                                                                    *
c *   pcg                                                              *
c *   gmres                                                            *
c *   prediag                                                          *
c *   pbcgstab                                                         *
c *                                                                    *
c **********************************************************************
c     subroutine solv(neq,nequ,nad,ip,ja,ad,au,al,m,b,x,tol,maxit,ngram,
c    .               unsym,solver,neqf1i,neqf2i,neq3i,neq4i,neq_doti,
c    .               i_fmapi,i_xfi,i_rcvsi,i_dspli)
c     use Malloc
c     implicit none
c     include 'mpif.h'
c     include 'parallel.fi'
c     include 'openmp.fi'
c     include 'time.fi'
c     integer neqf1i,neqf2i
c ... ponteiros      
c     integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c     integer*8 i_z,i_r,i_g,i_h,i_y,i_c,i_s
c ......................................................................
c     integer neq3i,neq4i,neq_doti
c     integer ip(*),ja(*),neq,nequ,nad
c     integer maxit,ngram,solver
c     real*8  ad(*),au(*),al(*),m(*),x(*),b(*),tol,energy
c     logical unsym
c     integer neqovlp
c     external matvec_csrsym1
c     external dot,dot_par
c     external matvec_csrc,matvec_csrcr,matvec_csrcsym,matvec_csrcrsym
c     external matvec_csrc1,matvec_csrcr1
c     external matvec_csrcsym1,matvec_csrcrsym1
c     OpenMP'ed subroutines
c     external dot_par_omp,dot_par_omp_loopwise
c     external matvec_csrc_omp,matvec_csrcsym_omp,
c    .         matvec_csrcr_omp,matvec_csrcrsym_omp
c ......................................................................
c ... numero total de equacoes na particao overlapping:
c    (neqovlp = neq, no sequencial e no non-overlapping)
c     neqovlp = neq+neq3i+neq4i
c     if (openmp) then
c        pmatrixtime = Mpi_Wtime() - pmatrixtime 
c        i_threads_y = alloc_8('buffer_y',nth_solv,neq)
c        call partition_matrix(ip,ja,neq,ovlp)
c        pmatrixtime = Mpi_Wtime() - pmatrixtime
c     endif
c ......................................................................
c
c ... Gradientes conjugados com precondicionador diagonal:
c     if (solver .eq. 1) then
c        if (unsym) then
c           print*,'Solver 1 nao disponivel para matriz nao-simetrica !'
c           call stop_mef()
c        endif
c        i_z = alloc_8('z       ',1,neq)
c         i_z = alloc_8('z       ',1,neq+1) + 2
c        i_r = alloc_8('r       ',1,neq)
c ...    precondicionador diagonal:
c        call aequalb(m,ad,neq)
c ...    Comunicacao da diagonal para o caso non-overlapping:
c        if (novlp) call communicate(m,neqf1i,neqf2i,i_fmapi,i_xfi,
c    .                   i_rcvsi,i_dspli)
c ......................................................................
c        if (ovlp) then
c ......... Overlapping:
c           if (openmp) then
c              call pcg_omp(neq,ip,ja,ad,au,al,m,b,x,ia(i_z),
c    .                  ia(i_r),tol,maxit,matvec_csrcrsym_omp,
c    .                  dot_par_omp,my_id,neqf1i,neqf2i,neq_doti,
c    .                  i_fmapi,i_xfi,i_rcvsi,i_dspli,ia(i_threads_y))
c           else
c              call pcg(neq,ip,ja,ad,au,al,m,b,x,ia(i_z),ia(i_r),
c    .                  tol,maxit,
c ... matvec comum:
c    .                  matvec_csrcrsym,dot_par,
c ... matvec desenrolado:
c    .                  matvec_csrcrsym1,dot_par,
c    .                  my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
c    .                  i_rcvsi,i_dspli)
c ......................................................................
c           endif
c        else
c ......... Sequencial e non-overlapping:
c           if (openmp) then
c              call pcg_omp(neq,ip,ja,ad,au,al,m,b,x,ia(i_z),ia(i_r),
c    .                  tol,maxit,matvec_csrcsym_omp,dot_par_omp,my_id,
c    .                  neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,i_rcvsi,
c    .                  i_dspli,ia(i_threads_y))
c           else
c              call pcg(neq,ip,ja,ad,au,al,m,b,x,ia(i_z),ia(i_r),
c    .                  tol,maxit,
c ... matvec comum:
c    .                  matvec_csrcsym,dot_par,
c ... matvec desenrolado:
c    .                  matvec_csrcsym1,dot_par,
c    .                  my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
c    .                  i_rcvsi,i_dspli)
c           endif
c        endif
c ......................................................................
c        i_r = dealloc('r       ')
c        i_z = dealloc('z       ')
c ......................................................................
c
c ... Gmres com precondicionador diagonal:
c     elseif(solver .eq. 2) then
c        i_g = alloc_8('g       ',neqovlp,ngram+1)         
c        i_h = alloc_8('h       ',ngram+1,ngram)
c        i_y = alloc_8('y       ',1,ngram)
c        i_c = alloc_8('c       ',1,ngram)
c        i_s = alloc_8('s       ',1,ngram)
c        i_r = alloc_8('r       ',1,ngram+1)
c ...... precondicionador diagonal:
c        call aequalb(m,ad,neq)
c ...... Comunicacao da diagonal para o caso non-overlapping:
c        if (novlp) call communicate(m,neqf1i,neqf2i,i_fmapi,i_xfi,
c    .                   i_rcvsi,i_dspli)
c        if(unsym) then
c ......................................................................
c           if(ovlp) then
c ............ Matriz nao-simetrica, overlapping:
c              if (openmp) then
c                 call gmres_omp(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
c    .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c    .                       tol,maxit,matvec_csrcr_omp,dot_par_omp,
c    .                       neqovlp,my_id,neqf1i,neqf2i,neq_doti,
c    .                       i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .                       ia(i_threads_y))
c              else
c                 call gmres(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
c    .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c ... matvec comum:
c     .                       maxit,matvec_csrcr,dot_par,neqovlp)
c ... matvec desenrolado:
c    .                       tol,maxit,matvec_csrcr1,dot_par,neqovlp,
c    .                       my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
c    .                       i_rcvsi,i_dspli)
c              endif
c              call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,
c    .                          i_rcvsi,i_dspli)
c ......................................................................
c           else
c ............ Matriz nao-simetrica, sequencial e non-overlapping:
c              if (openmp) then
c                 call gmres_omp(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
c    .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c    .                       tol,maxit,matvec_csrc_omp,dot_par_omp,
c    .                       neqovlp,my_id,neqf1i,neqf2i,neq_doti,
c    .                       i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .                       ia(i_threads_y))
c              else
c                 call gmres(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
c    .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c ... matvec comum:
c     .                       maxit,matvec_csrc,dot_par,neqovlp)
c ... matvec desenrolado:
c    .                       tol,maxit,matvec_csrc1,dot_par,neqovlp,
c    .                       my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
c    .                       i_rcvsi,i_dspli)
c             endif
c           endif
c ......................................................................
c        else
c ......................................................................
c           if (ovlp) then
c ............ Matriz simetrica, overlapping:
c              if (openmp) then
c                 call gmres_omp(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
c    .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c    .                       tol,maxit,matvec_csrcrsym_omp,dot_par_omp,
c    .                       neqovlp,my_id,neqf1i,neqf2i,neq_doti,
c    .                       i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .                       ia(i_threads_y))
c              else
c                 call gmres(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
c    .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c ... matvec comum:
c     .                       maxit,matvec_csrcrsym,dot_par,neqovlp)
c ... matvec desenrolado:
c    .                       tol,maxit,matvec_csrcrsym1,dot_par,neqovlp,
c    .                       my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
c    .                       i_rcvsi,i_dspli)
c              endif
c              call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
c    .                          i_dspli)
c ......................................................................
c           else
c ............ Matriz simetrica, sequencial e non-overlapping:
c              if (openmp) then
c                 call gmres_omp(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
c    .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c    .                       tol,maxit,matvec_csrcsym_omp,dot_par_omp,
c    .                       neqovlp,my_id,neqf1i,neqf2i,neq_doti,
c    .                       i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .                       ia(i_threads_y))
c             else
c                call gmres(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
c    .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c ... matvec comum:
c     .                       maxit,matvec_csrcsym,dot_par,neqovlp)
c ... matvec desenrolado:
c    .                       tol,maxit,matvec_csrcsym1,dot_par,neqovlp,
c    .                       my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
c    .                       i_rcvsi,i_dspli)
c              endif
c           endif
c ......................................................................
c        endif
c ......................................................................
c        i_r = dealloc('r       ')
c        i_s = dealloc('s       ')
c        i_c = dealloc('c       ')
c        i_c = dealloc('y       ')
c        i_h = dealloc('h       ')
c        i_g = dealloc('g       ')
c ......................................................................
c
c ... Gauss:
c     elseif(solver .eq. 3) then
c        time0 = MPI_Wtime()
c        call dtri(ad,au,al,ja,neq,unsym)
c        call dsolv(ad,au,al,b,ja,neq,energy,.true.)
c        time = MPI_Wtime()
c        print*,'CPU time (s) ',time-time0
c        x(1:neq) = b(1:neq)
c ......................................................................
c
c ... BICGSTAB com precondicionador diagonal:
c     elseif(solver .eq. 4) then
c ...    precondicionador diagonal:
c        call aequalb(m,ad,neq)      
c ...    Comunicacao da diagonal para o caso non-overlapping:
c        if (novlp) call communicate(m,neqf1i,neqf2i,i_fmapi,i_xfi,
c    .                               i_rcvsi,i_dspli)
c .....................................................................
c
c ...     
c        i_c = alloc_8('tsolver ',1,neq)
c        i_h = alloc_8('hsolver ',1,neq)
c        i_r = alloc_8('rsolver ',1,neq)
c        i_s = alloc_8('psolver ',1,neq)
c        i_z = alloc_8('zsolver ',1,neq)
c .....................................................................
c
c ............ Matriz nao-simetrica 
c        if(unsym) then
c ............ Matriz nao-simetrica, overlapping:
c          if(ovlp) then
c ... Openmp            
c            if(openmp) then
c              call pbicgstab_omp_loopwise(neq,ip,ja,ad,au,al,m,b,x,
c              call pbicgstab_omp(neq,ip,ja,ad,au,al,m,b,x,           
c    .                            ia(i_c),ia(i_h),   
c    .                            ia(i_r),ia(i_s),ia(i_z),tol,maxit,
c ... dot_par_omp_loopwise:     
c    .                            matvec_csrcr_omp,dot_par_omp_loopwise,
c ... dot_par_omp      
c    .                            matvec_csrcr_omp,dot_par_omp,     
c    .                            my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
c    .                            i_xfi,i_rcvsi,i_dspli,
c    .                            ia(i_threads_y))
c ... sem openmp     
c            else
c              call pbicgstab(neq,ip,ja,ad,au,al,m,b,x,
c              call bicgstab(neq,ip,ja,ad,au,al,m,b,x,
c    .                      ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z),
c    .                      tol,maxit,
c ... matvec comum:
c    .                      matvec_csrcr,dot_par,
c ... matvec desenrolado:
c    .                      matvec_csrcr1,dot_par,
c    .                      my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
c    .                      i_xfi,i_rcvsi,i_dspli)
c            endif
c            call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,
c    .                        i_rcvsi,i_dspli)
c .....................................................................
c
c ............ Matriz nao-simetrica, sequencial e non-overlapping:  
c          else
c ... Openmp            
c            if(openmp) then
c               call pbicgstab_omp_loopwise(neq,ip,ja,ad,au,al,m,b,x,
c               call pbicgstab_omp(neq,ip,ja,ad,au,al,m,b,x,           
c    .                            ia(i_c),ia(i_h),   
c    .                            ia(i_r),ia(i_s),ia(i_z),tol,maxit,
c ... dot_par_omp_loopwise:     
c    .                            matvec_csrc_omp,dot_par_omp_loopwise,
c ... dot_par_omp      
c    .                            matvec_csrc_omp,dot_par_omp,     
c    .                            my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
c    .                            i_xfi,i_rcvsi,i_dspli,
c    .                            ia(i_threads_y))
c .....................................................................
c     
c ... sem openmp     
c            else
c              call pbicgstab(neq,nequ,nad,ip,ja,ad,au,al,m,b,x,
c              call bicgstab(neq,ip,ja,ad,au,al,m,b,x,
c    .                      ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z),
c    .                      tol,maxit,
c ... matvec comum:
c     .                      matvec_csrc,dot_par,
c ... matvec desenrolado:
c    .                      matvec_csrc1,dot_par,
c    .                      my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
c    .                      i_xfi,i_rcvsi,i_dspli)
c           endif
c .....................................................................
c
c ....................................................................     
c          endif
c ....................................................................             
c
c ............ Matriz simetrica 
c        else       
c ............ Matriz simetrica, overlapping:
c          if(ovlp) then
c ... Openmp       
c            if(openmp) then
c              call pbicgstab_omp_loopwise(neq,ip,ja,ad,au,al,m,b,x,
c              call pbicgstab_omp(neq,ip,ja,ad,au,al,m,b,x,           
c    .                         ia(i_c),ia(i_h),   
c    .                         ia(i_r),ia(i_s),ia(i_z),tol,maxit,
c ... dot_par_omp_loopwise:     
c    .                         matvec_csrcrsym_omp,dot_par_omp_loopwise,
c ... dot_par_omp      
c    .                         matvec_csrcrsym_omp,dot_par_omp,     
c    .                         my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
c    .                         i_xfi,i_rcvsi,i_dspli,
c    .                         ia(i_threads_y))
c .....................................................................
c     
c ... sem openmp     
c            else
c              call pbicgstab(neq,ip,ja,ad,au,al,m,b,x,
c    .                      ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z),
c    .                      tol,maxit,
c ... matvec comum:
c    .                      matvec_csrcrsym,dot_par,
c ... matvec desenrolado:
c    .                      matvec_csrcrsym1,dot_par,
c    .                      my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
c    .                      i_xfi,i_rcvsi,i_dspli)
c            endif
c .....................................................................
c
c ...
c            call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
c    .                        i_dspli)
c .....................................................................
c
c ............ Matriz simetrica, sequencial e non-overlapping:
c          else
c ... Openmp       
c            if(openmp) then
c              call pbicgstab_omp_loopwise(neq,ip,ja,ad,au,al,m,b,x, 
c              call pbicgstab_omp(neq,ip,ja,ad,au,al,m,b,x,           
c    .                          ia(i_c),ia(i_h),   
c    .                          ia(i_r),ia(i_s),ia(i_z),tol,maxit,
c ... dot_par_omp_loopwise:     
c    .                          matvec_csrcsym_omp,dot_par_omp_loopwise,
c ... dot_par_omp      
c    .                          matvec_csrcsym_omp,dot_par_omp,     
c    .                          my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
c    .                          i_xfi,i_rcvsi,i_dspli,
c    .                          ia(i_threads_y))
c .....................................................................
c     
c ... sem openmp     
c            else
c              call pbicgstab(neq   ,nequ   ,nad    ,ip     ,ja     ,
c    .                       ad     ,au     ,al     ,m      ,b      ,x, 
c    .                       ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z),
c    .                       tol    ,maxit,
c ... matvec comum:
c    .                       matvec_csrcb  ,dot_par,
c    .                     matvec_csrcsym,dot_par,
c ... matvec desenrolado:
c    .                     matvec_csrcsym1,dot_par,
c    .                       my_id   ,neqf1i  ,neqf2i,
c    .                       neq_doti,i_fmapi,i_xfi  ,i_rcvsi,i_dspli)
c            endif
c .....................................................................
c          endif  
c .....................................................................       
c        endif
c .....................................................................
c
c ...             
c        i_z = dealloc('zsolver ')     
c        i_s = dealloc('psolver ')
c        i_r = dealloc('rsolver ')
c        i_h = dealloc('hsolver ')
c        i_c = dealloc('tsolver ')
c ......................................................................         
c     endif
c ......................................................................
c     if (openmp) then
c        pmatrixtime = Mpi_Wtime() - pmatrixtime 
c        i_threads_y = dealloc('buffer_y')
c        pmatrixtime = Mpi_Wtime() - pmatrixtime
c     endif
c     return
c     end
c **********************************************************************
