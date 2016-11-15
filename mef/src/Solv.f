c **********************************************************************
c *                                                                    *
c *   SOLV_PM.F                                            06/12/2015  *
c *                                                        28/06/2016  *
c *   Metodos iterativos de solucao:                                   *
c *                                                                    *
c *   cg                                                               *
c *   pcg                                                              *
c *   bcg                                                              *
c *   iccg                                                             *              
c *   minres                                                           *
c *   pminres                                                          *
c *   cr                                                               * 
c *   pcr                                                              * 
c *   sqrm                                                             *
c *   rsqrm                                                            *
c *   lsqrm                                                            *
c *   bicgstab                                                         *
c *   pbicgstab                                                        * 
c *   icbicgstab                                                       *                
c *   bicgstabl2                                                       *
c *   pbicgstabl2                                                      * 
c *   gmres(m)                                                         *
c *   gmres2(m)                                                        *
c *   block_it_pcg                                                     *
c **********************************************************************
      subroutine solv_pm(neq    ,nequ    ,neqp
     .                  ,nad    ,naduu   ,nadpp
     .                  ,ip     ,ja      ,ad         ,al 
     .                  ,m      ,b       ,x          ,tol    ,maxit
     .                  ,ngram  ,block_pu,n_blocks_up,solver ,istep
     .                  ,cmaxit ,ctol    ,alfap      ,alfau  ,precond
     .                  ,fmec   ,fporomec,fhist_solv 
     .                  ,neqf1i ,neqf2i  ,neq3i      ,neq4i  ,neq_doti
     .                  ,i_fmapi,i_xfi   ,i_rcvsi    ,i_dspli)
      use Malloc
      implicit none
      include 'precond.fi'
      include 'mpif.h'
      include 'parallel.fi'
      include 'openmp.fi'
      include 'time.fi'
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      integer*8 i_a,i_b,i_c,i_d,i_g,i_h,i_y,i_z,i_r,i_s
c ......................................................................
      integer neq3i,neq4i,neq_doti
      integer ip(*),ja(*),neq,nequ,neqp,nad,naduu,nadpp
      integer maxit,solver,ngram,istep,n_blocks_up
      real*8  ad(*),al(*),m(*),x(*),b(*),tol,energy
c ... pcg duplo
      integer cmaxit
      real*8  ctol,alfap,alfau
c ......................................................................
      logical block_pu,fmec,fporomec,fhist_solv
c ... precondicionador
      logical diag
      integer precond
      real*8  max_block_a(max_block*max_block)
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
        i_s = alloc_8('psolver ',1,neqovlp)
c ......................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_z)
     2                  ,precond,neq   ,nequ
     3                  ,neqf1i ,neqf2i,i_fmapi
     4                  ,i_xfi ,i_rcvsi,i_dspli
     5                  ,novlp,my_id)
c .....................................................................
c
c ... matriz aramazena em csrc blocado (Kuu,Kpp,Kpu)
        if(block_pu) then
           print*,"PCG nao disponivel para a matriz", 
     .            " blocada nao simetrica !!"
          stop   
c .....................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
        else
c ... omp
          if(omp_solv) then
            call call_cg_omp(neq   ,nequ    ,nad   ,ip      ,ja
     1                   ,ad       ,al     ,m      ,b       ,x   
     2                   ,ia(i_z)  ,ia(i_r),ia(i_s) 
     3                   ,tol      ,maxit  ,precond,iparam ,fhist_solv 
     4                   ,my_id    ,neqf1i ,neqf2i,neq_doti,i_fmapi
     5                   ,i_xfi ,i_rcvsi   ,i_dspli,ia(i_threads_y)
     6                   ,nprcs   ,ovlp    ,mpi)  
c .....................................................................
c
c ... sequencial (cg, pcg e iccg)
          else
            call call_cg(neq      ,nequ   ,nad    ,ip      ,ja
     1                   ,ad      ,al     ,m      ,b       ,x   
     2                   ,ia(i_z) ,ia(i_r),ia(i_s) 
     3                   ,tol     ,maxit  ,precond,iparam ,fhist_solv 
     4                   ,my_id   ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     5                   ,i_xfi   ,i_rcvsi,i_dspli
     6                   ,nprcs   ,ovlp   ,mpi)  
c .....................................................................
c
          endif      
c .....................................................................
        endif
c .....................................................................
c
c ...
        i_s = dealloc('psolver ')
        i_r = dealloc('rsolver ')
        i_z = dealloc('zsolver ')
c ......................................................................
c
c ... GMRES :
      elseif(solver .eq. 2) then
c ...
        if(block_pu) then
          if(n_blocks_up .eq. 1 ) then
            print*,"GMRES nao disponivel para a matriz",
     .             " blocada [Kuu Kpp]."
            stop   
          else if( n_blocks_up .eq. 3 ) then
            print*,"GMRES nao disponivel para a matriz",
     .           " blocada [Kuu] [Kpp] [Kpu]."
            stop   
          endif
        endif
c ......................................................................
c
c ...
         diag = .true.
c ...    precondicionador diagonal:
         if(diag) then 
           precondtime = Mpi_Wtime() - precondtime  
           call pre_diag(m,ad,neq,.false.)
           precondtime = Mpi_Wtime() - precondtime  
c ...    
         else                           
           m(1:neq) = 1.d0  
         endif           
c ......................................................................
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
     1                 ,ad ,al  ,al ,m ,b ,x,ngram,ia(i_g) 
     2                 ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r) 
     3                 ,tol    ,maxit   
c ... matvec comum:
     1                 ,matvec_csrc_pm_omp,dot_par_omp 
     2                 ,neqovlp 
     3                 ,my_id  ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     4                 ,i_xfi  ,i_rcvsi,i_dspli  
     5                 ,ia(i_threads_y),.true.)
c .....................................................................
c
c ... sequencial
           else 
             call gmres2(neq,nequ,nad,ip,ja 
     1                 ,ad ,al  ,al ,m ,b ,x,ngram,ia(i_g) 
     2                 ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r) 
     3                 ,tol    ,maxit   
c ... matvec comum:
     1                 ,matvec_csrc_pm  ,dot_par 
     2                 ,neqovlp 
     3                 ,my_id   ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     4                 ,i_xfi   ,i_rcvsi,i_dspli  
     5                 ,.true.)
           endif
c .....................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
         else
c ... omp
           if(omp_solv) then
             call gmres_omp(neq,nequ,nad,ip,ja 
     1                 ,ad ,al  ,al ,m ,b ,x,ngram,ia(i_g) 
     2                 ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r) 
     3                 ,tol    ,maxit   
c ... matvec comum:
     1                 ,matvec_csrc_sym_pm_omp,dot_par_omp 
     2                 ,neqovlp 
     3                 ,my_id   ,neqf1i ,neqf2i,neq_doti,i_fmapi
     4                 ,i_xfi   ,i_rcvsi,i_dspli 
     5                 ,ia(i_threads_y),.true.)
c .....................................................................
c
c ... sequencial
           else 
             call gmres2(neq,nequ,nad,ip,ja 
     1                 ,ad ,al  ,al ,m ,b ,x,ngram,ia(i_g) 
     2                 ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r) 
     3                 ,tol    ,maxit   
c ... matvec comum:
     1                 ,matvec_csrc_sym_pm,dot_par 
     2                 ,neqovlp
     3                 ,my_id  ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     4                 ,i_xfi  ,i_rcvsi,i_dspli 
     5                 ,.true.)
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
c ... BICGSTAB :
      else if (solver .eq. 4) then
c ...
        if(block_pu) then
          if(n_blocks_up .eq. 1 ) then
            print*,"BICGSTAB nao disponivel para a matriz",
     .             " blocada [Kuu Kpp]."
            stop   
          else if( n_blocks_up .eq. 3 ) then
            print*,"BICGSTAB nao disponivel para a matriz",
     .           " blocada [Kuu] [Kpp] [Kpu]."
            stop   
          endif
        endif
c ......................................................................
c
c ... alocacao dos arronjos auxiliares (6neq)    
         i_c = alloc_8('tsolver ',1,neq)
         i_h = alloc_8('hsolver ',1,neq)
         i_r = alloc_8('rsolver ',1,neq)
         i_s = alloc_8('psolver ',1,neq)
         i_z = alloc_8('zsolver ',1,neq)
         i_y = alloc_8('ysolver ',1,neq)
c .....................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_z)
     2                  ,precond,neq   ,nequ
     3                  ,neqf1i ,neqf2i,i_fmapi
     4                  ,i_xfi ,i_rcvsi,i_dspli
     5                  ,novlp,my_id)
c .....................................................................
c
c ... matriz aramazena em csrc blocado (Kuu,Kpp,Kpu)
         if(block_pu) then
c ... omp
           if(omp_solv) then
             call pbicgstab_omp(neq   ,nequ,nad    ,ip   ,ja      
     1                     ,ad     ,al     ,al     ,m      ,b      ,x  
     2                     ,ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z)
     3                     ,tol    ,maxit 
c ... matvec comum:
     1                     ,matvec_csrc_pm_omp,dot_par_omp 
     2                     ,my_id  ,neqf1i  ,neqf2i  ,neq_doti   
     3                     ,i_fmapi,i_xfi  ,i_rcvsi,i_dspli
     4                     ,ia(i_threads_y),.true.)
c .....................................................................
c
c ... sequencial
           else 
             call pbicgstab(neq   ,nequ   ,nad  ,ip     ,ja      
     1                ,ad     ,al     ,al   ,m      ,b      ,x  
     2                ,ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z),ia(i_y)
     3                ,tol    ,maxit 
c ... matvec comum:
     1                ,matvec_csrc_pm ,dot_par 
     2                ,my_id    ,neqf1i,neqf2i ,neq_doti 
     3                ,i_fmapi ,i_xfi  ,i_rcvsi,i_dspli
     4                ,.true.  ,.true.,.true.)
           endif 
c .....................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
         else
c ... omp
           if(omp_solv) then
             call pbicgstab_omp(neq   ,nequ,nad    ,ip   ,ja      
     1                     ,ad     ,al     ,al     ,m      ,b      ,x  
     2                     ,ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z)
     3                     ,tol    ,maxit 
c ... matvec comum:
     1                     ,matvec_csrc_sym_pm_omp,dot_par_omp 
     2                     ,my_id  ,neqf1i  ,neqf2i  ,neq_doti   
     3                     ,i_fmapi,i_xfi  ,i_rcvsi,i_dspli
     4                     ,ia(i_threads_y),.true.)
c .....................................................................
c
c ... sequencial (bigstab, pbigstab e icbigstab)
           else
             call call_bicgstab(neq      ,nequ   ,nad    ,ip,ja
     1                    ,ad       ,al  ,m      ,b      ,x ,ia(i_y)
     2                    ,ia(i_c)  ,ia(i_h),ia(i_r),ia(i_s),ia(i_z)   
     3                    ,tol      ,maxit  ,precond
     4                    ,my_id    ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     5                    ,i_xfi    ,i_rcvsi,i_dspli) 
           endif      
c .....................................................................
         endif
c .....................................................................
c
c ... 
         i_y = dealloc('ysolver ')   
         i_z = dealloc('zsolver ')     
         i_s = dealloc('psolver ')
         i_r = dealloc('rsolver ')
         i_h = dealloc('hsolver ')
         i_c = dealloc('tsolver ')
c ......................................................................         
c
c ... PCG_BLOCK_IT                                                                   
      else if (solver .eq. 5) then
c ...
         if(omp_solv) then 
           print*,"PCG_BOLCK_IT: openmp nao disponivel !!"
         endif 
c .....................................................................
c
c ...
         diag = .true.
c ...    precondicionador diagonal:
         if(diag) then 
           precondtime = Mpi_Wtime() - precondtime    
           call pre_diag(m        ,ad        ,nequ,.false.)
           call pre_diag(m(nequ+1),ad(nequ+1),neqp,.false.)
           precondtime = Mpi_Wtime() - precondtime   
c ...    
         else                           
           m(1:neq) = 1.d0  
         endif           
c ......................................................................
c
c ... BLOCK_IT_PCG :
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
           i_d  = alloc_8('dsolver ',1,nequ)
           i_c  = alloc_8('csolver ',1,neqp)
           i_h  = alloc_8('hsolver ',1,nequ)
           i_g  = alloc_8('gsolver ',1,neqp)
           i_y  = alloc_8('ysolver ',1,nequ)
           i_a  = alloc_8('asolver ',1,neqp)
c ......................................................................
c
c ...
           call pcg_block_it(neq     ,nequ       ,neqp  
     1                    ,nad       ,naduu     ,nadpp
     2                    ,ip        ,ja     
     3                    ,ip(nequ+2),ja(naduu+1)
     4                    ,ip(neq+3) ,ja(naduu+nadpp+1)         
     5                    ,ad        ,ad(nequ+1)
     6                    ,al        ,al(naduu+1),al(naduu+nadpp+1)  
     7                    ,m         ,m(nequ+1)  ,b               ,x
     8                    ,ia(i_z)   ,ia(i_r),ia(i_s) 
     9                    ,ia(i_d)   ,ia(i_c),ia(i_h),ia(i_g)
     1                    ,ia(i_y)   ,ia(i_a)
     2                    ,tol       ,ctol   ,maxit    ,cmaxit
     3                    ,alfap     ,alfau 
     4                    ,.true.    ,istep
     5                    ,my_id     ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                    ,i_xfi     ,i_rcvsi,i_dspli) 
c ......................................................................
c
c ... 
         else 
           print*,"PCG_BLOCK indisponivel para a matriz nao blocada"
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
         i_d  = dealloc('dsolver ')       
         i_s  = dealloc('ssolver ')         
         i_r  = dealloc('rsolver ')        
         i_z  = dealloc('zsolver ') 
c ......................................................................
c
c ... BICGSTAB(2):
      else if(solver .eq. 6 ) then
c ...
         if(omp_solv) then 
           print*,"BICGSTAB(2): openmp nao disponivel !!"
         endif 
c .....................................................................
c
c ...
        if(block_pu) then
          if(n_blocks_up .eq. 1 ) then
            print*,"BICGSTAB(2) nao disponivel para a matriz",
     .             " blocada [Kuu Kpp]."
            stop   
          else if( n_blocks_up .eq. 3 ) then
            print*,"BICGSTAB(2) nao disponivel para a matriz",
     .           " blocada [Kuu] [Kpp] [Kpu]."
            stop   
          endif
        endif
c ......................................................................
c
c ... alocacao dos arronjos auxiliares (8neq)    
        i_c = alloc_8('tsolver ',1,neq)
        i_h = alloc_8('hsolver ',1,neq)
        i_r = alloc_8('rsolver ',1,neq)
        i_s = alloc_8('psolver ',1,neq)
        i_z = alloc_8('zsolver ',1,neq)
        i_y = alloc_8('ysolver ',1,neq)
        i_a = alloc_8('asolver ',1,neq)
        i_d = alloc_8('dsolver ',1,neq)
        i_g = alloc_8('gsolver ',1,neqovlp)
        i_b = alloc_8('bsolver ',1,neqovlp)
c ....................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_z)
     2                  ,precond,neq   ,nequ
     3                  ,neqf1i ,neqf2i,i_fmapi
     4                  ,i_xfi ,i_rcvsi,i_dspli
     5                  ,novlp,my_id)
c .....................................................................
c
c ... matriz aramazena em csrc blocado (Kuu,Kpp,Kpu)
        if(block_pu) then
           call call_bicgstabl2_v2(neq   ,nequ   ,nad,ip ,ja 
     1          ,ad      ,al     ,m      ,b      ,x   
     2          ,ia(i_c) ,ia(i_h),ia(i_r),ia(i_s),ia(i_z)
     3          ,ia(i_y) ,ia(i_a),ia(i_d),ia(i_b),ia(i_g)
     4          ,tol     ,maxit  ,precond
     5          ,my_id ,neqf1i   ,neqf2i,neq_doti,i_fmapi
     6          ,i_xfi ,i_rcvsi  ,i_dspli)
c .....................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
        else
           call call_bicgstabl2(neq       ,nequ   ,nad  ,ip     ,ja
     1                  ,ad       ,al     ,m      ,b    ,x    
     2                  ,ia(i_c) ,ia(i_h) ,ia(i_r),ia(i_s),ia(i_z)
     3                  ,ia(i_y) ,ia(i_a) ,ia(i_d),ia(i_b),ia(i_g)
     4                  ,tol      ,maxit  ,precond
     5                  ,my_id    ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                  ,i_xfi    ,i_rcvsi,i_dspli
     7                  ,nprcs    ,ovlp   ,mpi)  
        endif
c ......................................................................
c
c ... 
         i_b = dealloc('bsolver ')
         i_g = dealloc('gsolver ')
         i_d = dealloc('dsolver ')
         i_a = dealloc('asolver ')
         i_y = dealloc('ysolver ')   
         i_z = dealloc('zsolver ')     
         i_s = dealloc('psolver ')
         i_r = dealloc('rsolver ')
         i_h = dealloc('hsolver ')
         i_c = dealloc('tsolver ')
c ...................................................................... 
c
c ... MINRES:
      else if(solver .eq. 7 ) then
c ...
         if(omp_solv) then 
           print*,"MINRES: openmp nao disponivel !!"
         endif 
c .....................................................................
c
c ... matriz aramazena em csrc blocado nao simentrico (Kuu,Kpp,Kpu)
        if(block_pu) then
          print*,"MINRES nao disponivel para a matriz", 
     .           " blocada nao simetrica !!"
          stop 
        endif    
c .....................................................................
        
c ...
        i_c = alloc_8('tsolver ',1,neq)
        i_h = alloc_8('hsolver ',1,neq)
        i_r = alloc_8('rsolver ',1,neq)
        i_s = alloc_8('psolver ',1,neq)
        i_z = alloc_8('zsolver ',1,neq)
        i_y = alloc_8('ysolver ',1,neq)
        i_a = alloc_8('asolver ',1,neq)
        i_g = alloc_8('gsolver ',1,neq)
c .....................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_z)
     2                  ,precond,neq   ,nequ
     3                  ,neqf1i ,neqf2i,i_fmapi
     4                  ,i_xfi ,i_rcvsi,i_dspli
     5                  ,novlp,my_id)
c .....................................................................
c
c ...
        call call_minres(neq      ,nequ  ,nad   ,ip      ,ja
     1               ,ad       ,al    ,m     ,b       ,x    
     2               ,ia(i_c)  ,ia(i_h),ia(i_r),ia(i_s),ia(i_z),ia(i_y) 
     3               ,ia(i_a)  ,ia(i_g)  
     4               ,tol      ,maxit ,precond,iparam,fhist_solv 
     5               ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     6               ,i_xfi    ,i_rcvsi,i_dspli)
c .....................................................................
c
c ...
        i_g = dealloc('gsolver ')
        i_a = dealloc('asolver ')
        i_y = dealloc('ysolver ')   
        i_z = dealloc('zsolver ')     
        i_s = dealloc('psolver ')
        i_r = dealloc('rsolver ')
        i_h = dealloc('hsolver ')
        i_c = dealloc('tsolver ')
c .....................................................................
c
c ... CR - Conjugate Residual
      else if(solver .eq. 8 ) then
c ...
         if(omp_solv) then 
           print*,"CR: openmp nao disponivel !!"
         endif 
c .....................................................................
c
c ... matriz aramazena em csrc blocado nao simentrico (Kuu,Kpp,Kpu)
        if(block_pu) then
          print*,"CR nao disponivel para a matriz", 
     .           " blocada nao simetrica !!"
          stop 
        endif    
c .....................................................................
c
c ...
        i_c = alloc_8('tsolver ',1,neq)
        i_h = alloc_8('hsolver ',1,neq)
        i_r = alloc_8('rsolver ',1,neq)
        i_s = alloc_8('psolver ',1,neq)
        i_z = alloc_8('zsolver ',1,neq)
        i_y = alloc_8('ysolver ',1,neq)
c .....................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_z)
     2                  ,precond,neq   ,nequ
     3                  ,neqf1i ,neqf2i,i_fmapi
     4                  ,i_xfi ,i_rcvsi,i_dspli
     5                  ,novlp,my_id)
c .....................................................................
c
c ...   
        call call_cr(neq      ,nequ  ,nad,ip      ,ja
     1              ,ad       ,al    ,m     ,b       ,x    
     2              ,ia(i_c),ia(i_h),ia(i_r) ,ia(i_s),ia(i_z)
     3              ,ia(i_y)  
     4              ,tol      ,maxit ,precond,iparam  ,fhist_solv 
     5              ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     6              ,i_xfi    ,i_rcvsi,i_dspli)
c .....................................................................
c
c ...
        i_y = dealloc('ysolver ')   
        i_z = dealloc('zsolver ')     
        i_s = dealloc('psolver ')
        i_r = dealloc('rsolver ')
        i_h = dealloc('hsolver ')
        i_c = dealloc('tsolver ')
c .....................................................................
c
c ... SYMMLQ:
      else if(solver .eq. 9 ) then
c ...
         if(omp_solv) then 
           print*,"SYMMLQ: openmp nao disponivel !!"
         endif 
c .....................................................................
c
c ... matriz aramazena em csrc blocado nao simentrico (Kuu,Kpp,Kpu)
        if(block_pu) then
          print*,"SYMMLQ nao disponivel para a matriz", 
     .           " blocada nao simetrica !!"
          call stop_mef() 
        endif    
c .....................................................................
c
c ...
        i_c = alloc_8('tsolver ',1,neq)
        i_h = alloc_8('hsolver ',1,neq)
        i_r = alloc_8('rsolver ',1,neq)
        i_s = alloc_8('psolver ',1,neq)
        i_z = alloc_8('zsolver ',1,neq)
        i_y = alloc_8('ysolver ',1,neq)
c .....................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_z)
     2                  ,precond,neq   ,nequ
     3                  ,neqf1i ,neqf2i,i_fmapi
     4                  ,i_xfi ,i_rcvsi,i_dspli
     5                  ,novlp,my_id)
c .....................................................................
c
c ...
        call call_symmlq(neq      ,nequ  ,nad   ,ip      ,ja
     1                  ,ad       ,al    ,m     ,b       ,x    
     2                  ,ia(i_c)  ,ia(i_h)      ,ia(i_r),ia(i_s)
     3                  ,ia(i_z)  ,ia(i_y) 
     4                  ,tol      ,maxit ,precond,iparam ,fhist_solv
     5                  ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     6                  ,i_xfi    ,i_rcvsi,i_dspli)
c .....................................................................
c
c ...
        i_y = dealloc('ysolver ')   
        i_z = dealloc('zsolver ')     
        i_s = dealloc('psolver ')
        i_r = dealloc('rsolver ')
        i_h = dealloc('hsolver ')
        i_c = dealloc('tsolver ')
c .....................................................................
c
c ... mkl_pardiso
      else if(solver .eq. 10 ) then
        i_z  = alloc_8('zsolver ',1,neq)
        if(fmec) then
          call call_mkl_pardiso(neq,ip,ja,ad,b,x,ia(i_z),2)
        else if(fporomec) then
          call call_mkl_pardiso(neq,ip,ja,ad,b,x,ia(i_z),-2)
        endif
        i_z  = dealloc('zsolver ') 
c ......................................................................
c
c ... SQRM - QRM simetrico
      else if(solver .eq. 11) then
c ... matriz aramazena em csrc blocado nao simentrico (Kuu,Kpp,Kpu)
        if(block_pu) then
          print*,"SQRM nao disponivel para a matriz", 
     .           " blocada nao simetrica !!"
          stop 
        endif    
c .....................................................................
c
c ...
        i_z = alloc_8('zsolver ',1,neq)
        i_h = alloc_8('hsolver ',1,neq)
        i_r = alloc_8('rsolver ',1,neqovlp)
        i_s = alloc_8('psolver ',1,neq)
c .....................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_z)
     2                  ,precond,neq   ,nequ
     3                  ,neqf1i ,neqf2i,i_fmapi
     4                  ,i_xfi ,i_rcvsi,i_dspli
     5                  ,novlp ,my_id)
c .....................................................................
c
c ...
        if(omp_solv)then
          call call_sqrm_omp(neq  ,nequ  ,nad   ,ip      ,ja
     1                    ,ad     ,al    ,m     ,b       ,x    
     2                    ,ia(i_z),ia(i_h),ia(i_r) ,ia(i_s)     
     3                    ,tol    ,maxit ,precond,iparam ,fhist_solv
     4                    ,my_id  ,neqf1i,neqf2i ,neq_doti,i_fmapi
     5                    ,i_xfi  ,i_rcvsi,i_dspli,ia(i_threads_y)
     6                    ,nprcs  ,ovlp   ,mpi) 
c .....................................................................
c
c ...
        else
          call call_sqrm(neq    ,nequ  ,nad     ,ip      ,ja
     1                ,ad       ,al    ,m       ,b       ,x    
     2                ,ia(i_z),ia(i_h),ia(i_r)  ,ia(i_s) 
     3                ,tol      ,maxit ,precond ,iparam  ,fhist_solv
     4                ,my_id    ,neqf1i,neqf2i  ,neq_doti,i_fmapi
     5                ,i_xfi    ,i_rcvsi,i_dspli
     6                ,nprcs    ,ovlp   ,mpi) 
       endif
c .....................................................................
c
c ...
        i_s = dealloc('psolver ')
        i_r = dealloc('rsolver ')
        i_h = dealloc('hsolver ')
        i_z = dealloc('zsolver ')
c .....................................................................
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
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 11/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_CG : chama a versao do gradiente conjudado desejada           *    
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
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * pc       - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * my_id    - MPI                                                     *
c * neqf1i   - MPI                                                     *
c * neqf2i   - MPI                                                     *
c * neq_doti - MPI                                                     *
c * i_fmap   - MPI                                                     *
c * i_xfi    - MPI                                                     *
c * i_rvcs   - MPI                                                     *
c * i_dspli  - MPI                                                     *
c * mpi      - true|false                                              *
c * ovlp     - overllaping                                             *
c * nprcs    - numero de processos mpi                                 *  
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_cg(neq      ,nequ  ,nad   ,ia      ,ja
     1                  ,ad       ,al    ,m     ,b       ,x    
     2                  ,z        ,r     ,s    
     3                  ,tol      ,maxit ,precond,iparam ,fhist_log
     4                  ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     5                  ,i_xfi    ,i_rcvsi,i_dspli
     6                  ,nprcs    ,ovlp   ,mpi    )
      implicit none
      include 'time.fi'
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,nequ,nad,neq_doti 
      integer ia(*),ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 z(*),r(*),s(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log
c ... precondicionador
      integer precond,iparam(*)
      real*8 m(*)
c ...
      external dot_par
      external matvec_csrc_sym_pm,matvec_csrcr_sym_pm  
      external ildlt_solv,illt_solv  
c ......................................................................

c ... cg
      if(precond .eq. 1) then
c ...  
        call cg(neq    ,nequ   ,nad,ia   ,ja
     1         ,ad     ,al     ,al ,b    ,x
     2         ,z      ,r      ,s  ,tol  ,maxit
c ... matvec comum:
     3         ,matvec_csrc_sym_pm,dot_par 
     4         ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     5         ,i_xfi ,i_rcvsi,i_dspli
     6         ,.true.,.true. ,fhist_log ,.true.)
c .....................................................................
c
c ... pcg - cg com precondicionador diagonal
      else if(precond .eq. 2 .or. precond .eq. 5 ) then
c ... overllaping
         if(ovlp) then  
           call pcg(neq    ,nequ   ,nad,ia   ,ja
     1             ,ad     ,al     ,al ,m    ,b
     2             ,x      ,z      ,r  ,s
     3             ,tol    ,maxit
c ... matvec comum:
     4             ,matvec_csrcr_sym_pm,dot_par 
c
     5             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6             ,i_xfi ,i_rcvsi,i_dspli
     7             ,.true.,.true. ,fhist_log,.true.
     8             ,nprcs ,mpi)
c .....................................................................
c
c ...non-overllaping
         else 
           call pcg(neq    ,nequ   ,nad,ia   ,ja
     1             ,ad     ,al     ,al ,m    ,b
     2             ,x      ,z      ,r  ,s
     3             ,tol    ,maxit
c ... matvec comum:
     4             ,matvec_csrc_sym_pm,dot_par 
c
     5             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6             ,i_xfi ,i_rcvsi,i_dspli
     7             ,.true.,.true. ,fhist_log,.true.
     8             ,nprcs ,mpi)
         endif
c .....................................................................
c
c ... iccg - cg com precondicionador LDLT(0) imcompleto
      elseif(precond .eq. 3 ) then
        call iccg(neq      ,nequ  ,nad   ,ia      ,ja
     1           ,ad       ,al    ,al    ,m       ,b    
     2           ,x        ,z     ,r     ,s   
     3           ,tol      ,maxit
c ... matvec comum:
     4           ,matvec_csrc_sym_pm,dot_par,ildlt_solv 
     5           ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6           ,i_xfi ,i_rcvsi,i_dspli
     7           ,.true.,.true. ,fhist_log,.true.)
c .....................................................................
c
c ... iccg - cg com precondicionador LLT(0) imcompleto
      elseif(precond .eq. 4 ) then
        call iccg(neq      ,nequ  ,nad   ,ia      ,ja
     1           ,ad       ,al    ,al    ,m       ,b    
     2           ,x        ,z     ,r     ,s   
     3           ,tol      ,maxit
c ... matvec comum:
     4           ,matvec_csrc_sym_pm,dot_par,illt_solv
     5           ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6           ,i_xfi ,i_rcvsi,i_dspli
     7           ,.true.,.true. ,fhist_log,.true.)
c .....................................................................
c
c ... bpcg - cg com bloco diagonal 
      elseif(precond .eq. 6 ) then
        call bpcg(neq      ,nequ  ,nad   ,ia      ,ja
     1           ,ad       ,al    ,al    ,m       ,b    
     2           ,x        ,z     ,r     ,s   
     3           ,tol      ,maxit 
c ... matvec comum:
     4           ,matvec_csrc_sym_pm,dot_par
     5           ,my_id ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     6           ,i_xfi ,i_rcvsi,i_dspli
     7           ,.true.,.true. ,fhist_log ,.true.)
      endif  
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 23/09/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_CG_OMP : chama a versao do gradiente conjudado desejada       *    
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
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * pc       - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * my_id    - MPI                                                     *
c * neqf1i   - MPI                                                     *
c * neqf2i   - MPI                                                     *
c * neq_doti - MPI                                                     *
c * i_fmap   - MPI                                                     *
c * i_xfi    - MPI                                                     *
c * i_rvcs   - MPI                                                     *
c * i_dspli  - MPI                                                     *
c * thread_y - buffer de equacoes para o vetor y (openmp)              *  
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_cg_omp(neq      ,nequ  ,nad   ,ia      ,ja
     1                  ,ad       ,al    ,m     ,b       ,x    
     2                  ,z        ,r     ,s    
     3                  ,tol      ,maxit ,precond,iparam ,fhist_log
     4                  ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     5                  ,i_xfi    ,i_rcvsi,i_dspli,thread_y
     6                  ,nprcs    ,ovlp   ,mpi    )
      implicit none
      include 'time.fi'
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,nequ,nad,neq_doti 
      integer ia(*),ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 z(*),r(*),s(*)
c ... buffer do CSR omp
      real*8 thread_y(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log
c ... precondicionador
      integer precond,iparam(*)
      real*8 m(*)
c ...
      external dot_par_omp
      external matvec_csrc_sym_pm_omp,matvec_csrcr_sym_pm_omp  
      external ildlt_solv,illt_solv  
c ......................................................................

c ... cg
      if(precond .eq. 1) then
        print*,"CG: openmp nao disponivel !!"  
c
c ... pcg - cg com precondicionador diagonal
      else if(precond .eq. 2 .or. precond .eq. 5 ) then
c ... overllaping
        if(ovlp) then  
          call pcg_omp(neq    ,nequ   ,nad,ia   ,ja
     1                ,ad     ,al     ,al ,m    ,b 
     2                ,x      ,z      ,r  ,s
     3                ,tol    ,maxit
c ... matvec comum:
     4                ,matvec_csrcr_sym_pm_omp,dot_par_omp
c             
     5                ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7                ,.true.,.true. ,fhist_log,.true.
     8                ,nprcs ,mpi)
c .....................................................................
c
c ...non-overllaping
        else 
          call pcg_omp(neq    ,nequ   ,nad,ia   ,ja
     1                ,ad     ,al     ,al ,m    ,b 
     2                ,x      ,z      ,r  ,s
     3                ,tol    ,maxit
c ... matvec comum:
     4                ,matvec_csrc_sym_pm_omp,dot_par_omp
c             
     5                ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7                ,.true.,.true. ,fhist_log,.true.
     8                ,nprcs ,mpi)
        endif
c .....................................................................
c
c ... iccg - cg com precondicionador LDLT(0) imcompleto
      elseif(precond .eq. 3 ) then
        print*,"ICCG: openmp nao disponivel !!" 
c .....................................................................
c
c ... iccg - cg com precondicionador LLT(0) imcompleto
      elseif(precond .eq. 4 ) then
        print*,"ICCG: openmp nao disponivel !!" 
c .....................................................................
c
c ... bpcg - cg com bloco diagonal 
      elseif(precond .eq. 6 ) then
        print*,"BCCG: openmp nao disponivel !!" 
      endif  
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 27/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_SYMMLQ : chama a versao do SYMMLQ                             *    
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
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * y(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * pc       - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_symmlq(neq      ,nequ  ,nad   ,ia      ,ja
     1                  ,ad       ,al    ,m     ,b       ,x    
     2                  ,c        ,h     ,r     ,s       ,z  ,y 
     3                  ,tol      ,maxit ,precond,iparam ,fhist_log
     4                  ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     5                  ,i_xfi    ,i_rcvsi,i_dspli)
      implicit none
c ... mpi
      integer my_id
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c .....................................................................
      integer neq,nequ,nad,neq_doti 
      integer ia(*),ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 c(*),h(*),r(*),s(*),z(*),y(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log
c ... precondicionador
      integer precond,iparam(*)
      real*8 m(*)
c ...
      external dot_par
      external matvec_csrc_sym_pm  
      external ildlt_solv,illt_solv  
c ......................................................................
c
c ...   
        if(precond .eq. 1 ) then
          call symmlq(neq  ,nequ   ,nad  ,ia     ,ja
     1          ,ad     ,al     ,al      ,b      ,x
     2          ,c      ,h      ,r       ,s      
     3          ,tol    ,maxit
c ... matvec comum:
     4          ,matvec_csrc_sym_pm,dot_par 
     5          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6          ,i_xfi ,i_rcvsi,i_dspli
     7          ,.true.,.true.,.true.)
c .....................................................................
c
c ... precondicionador diagonal:
        else if(precond .eq. 2 .or. precond .eq. 5 
     1         .or. precond .eq. 7) then
          call psymmlq(neq  ,nequ   ,nad ,ia   ,ja
     1          ,ad     ,al  ,al     ,b  ,m    ,x
     2          ,c      ,h   ,r      ,s  ,z    ,y      
     3          ,tol    ,maxit
c ... matvec comum:
     4          ,matvec_csrc_sym_pm,dot_par 
     5          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6          ,i_xfi ,i_rcvsi,i_dspli
     7          ,.true.,.true.,.true.)
c .....................................................................
        endif
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 27/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_CR : chama a versao do CR                                     *    
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
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * y(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * pc       - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_cr(neq      ,nequ  ,nad   ,ia      ,ja
     1                  ,ad       ,al    ,m     ,b       ,x    
     2                  ,c        ,h     ,r     ,s       ,z  ,y 
     3                  ,tol      ,maxit ,precond,iparam ,fhist_log
     4                  ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     5                  ,i_xfi    ,i_rcvsi,i_dspli)
      implicit none
c ... mpi
      integer my_id
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c .....................................................................
      integer neq,nequ,nad,neq_doti 
      integer ia(*),ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 c(*),h(*),r(*),s(*),z(*),y(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log
c ... precondicionador
      integer precond,iparam(*)
      real*8 m(*)
c ...
      external dot_par
      external matvec_csrc_sym_pm  
      external ildlt_solv,illt_solv  
c ......................................................................
c
      if(precond .eq. 1 ) then
        call cr(neq  ,nequ   ,nad        ,ia     ,ja
     1          ,ad     ,al     ,al      ,b      ,x
     2          ,c      ,h      ,r       ,s      ,z
     3          ,tol    ,maxit
c ... matvec comum:
     4          ,matvec_csrc_sym_pm,dot_par 
     5          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6          ,i_xfi ,i_rcvsi,i_dspli
     7          ,.true.,.true.,.true.)
c .....................................................................
c
c ...
       else if(precond .eq. 2 
     1         .or. precond .eq. 5 .or.  precond .eq. 7) then
         call pcr(neq  ,nequ   ,nad    ,ia    ,ja
     1          ,ad     ,al    ,al     ,b     ,m  ,x
     2          ,c      ,h     ,r      ,s     ,z  ,y  
     3          ,tol    ,maxit
c ... matvec comum:
     4          ,matvec_csrc_sym_pm,dot_par 
     5          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6          ,i_xfi ,i_rcvsi,i_dspli
     7          ,.true.,.true.,.true.)
c .....................................................................
        endif
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/06/2016                                    *
c * Data de modificaco : 28/10/2016                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_SQRM:chama a versao do metodo QMR simetrico                   *    
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
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * pc       - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * my_id    - MPI                                                     *
c * neqf1i   - MPI                                                     *
c * neqf2i   - MPI                                                     *
c * neq_doti - MPI                                                     *
c * i_fmap   - MPI                                                     *
c * i_xfi    - MPI                                                     *
c * i_rvcs   - MPI                                                     *
c * i_dspli  - MPI                                                     *
c * mpi      - true|false                                              *
c * ovlp     - overllaping                                             *
c * nprcs    - numero de processos mpi                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_sqrm(neq      ,nequ  ,nad      ,ia      ,ja
     1                    ,ad       ,al    ,m        ,b       ,x    
     2                    ,c        ,h     ,r        ,s    
     3                    ,tol      ,maxit ,pc       ,iparam ,fhist_log
     4                    ,my_id    ,neqf1i,neqf2i   ,neq_doti,i_fmapi
     5                    ,i_xfi    ,i_rcvsi,i_dspli
     6                    ,nprcs    ,ovlp   ,mpi    )
      implicit none
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,nequ,nad,neq_doti 
      integer ia(*),ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 c(*),h(*),r(*),s(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log
c ... precondicionador
      integer pc,iparam(*)
      real*8 m(*)
c ...
      external dot_par
      external matvec_csrc_sym_pm,matvec_csrcr_sym_pm  
      external ildlt_solv,illt_solv  
c ......................................................................
c
      if(pc .eq. 1 ) then
        call sqrm(neq    ,nequ   ,nad    ,ia    ,ja
     1           ,ad     ,al     ,al     ,m     ,b   ,x 
     2           ,c      ,h      ,r      ,s 
     3           ,tol   ,maxit
     4           ,matvec_csrc_sym_pm,dot_par
     5           ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6           ,i_xfi ,i_rcvsi,i_dspli
     7           ,.true.,.true. ,fhist_log,.true.)
c .....................................................................
c
c ...
       else if(pc .eq. 2 .or. pc .eq. 5 .or.  pc .eq. 7) then
c ... overllaping
         if(ovlp) then
           call rpsqrm(neq    ,nequ   ,nad    ,ia    ,ja
     1                ,ad     ,al     ,al      ,m    ,b   ,x 
     2                ,c      ,h      ,r       ,s 
     3                ,tol   ,maxit
     4                ,matvec_csrcr_sym_pm,dot_par
     5                ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                ,i_xfi ,i_rcvsi,i_dspli
     7                ,.true.,.true. ,fhist_log,.true.
     8                ,nprcs ,mpi)
c .....................................................................
c
c ...non-overllaping
         else 
           call rpsqrm(neq    ,nequ   ,nad    ,ia    ,ja
     1                ,ad     ,al     ,al      ,m     ,b   ,x 
     2                ,c      ,h      ,r       ,s 
     3                ,tol   ,maxit
     4                ,matvec_csrc_sym_pm,dot_par
     5                ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                ,i_xfi ,i_rcvsi,i_dspli
     7                ,.true.,.true. ,fhist_log,.true.
     8                ,nprcs ,mpi)
          endif
c .....................................................................
        endif
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 22/09/2016                                    *
c * Data de modificaco : 01/11/2016                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_SQRM_OMP:chama a versao do metodo QMR simetrico               *    
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
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * pc       - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * thread_y - buffer de equacoes para o vetor y (openmp)              * 
c * mpi      - true|false                                              *
c * ovlp     - overllaping                                             *
c * nprcs    - numero de processos mpi                                 *     
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_sqrm_omp(neq      ,nequ  ,nad   ,ia      ,ja
     1                    ,ad       ,al    ,m     ,b       ,x    
     2                    ,c        ,h     ,r     ,s    
     3                    ,tol      ,maxit ,pc    ,iparam ,fhist_log
     4                    ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     5                    ,i_xfi    ,i_rcvsi,i_dspli,thread_y
     6                    ,nprcs    ,ovlp   ,mpi    )
      implicit none
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,nequ,nad,neq_doti 
      integer ia(*),ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 c(*),h(*),r(*),s(*)
c ... buffer do CSR omp
      real*8 thread_y(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log
c ... precondicionador
      integer pc,iparam(*)
      real*8 m(*)
c ...
      external dot_par_omp
      external matvec_csrc_sym_pm_omp,matvec_csrcr_sym_pm_omp
      external ildlt_solv,illt_solv  
c ......................................................................
c
      if(pc .eq. 1 ) then
        print*,"SQRM: openmp nao disponivel !!"  
c .....................................................................
c
c ...
      else if(pc .eq. 2 .or. pc .eq. 5 .or.  pc .eq. 7) then
c ... overllaping
        if(ovlp) then
          call rpsqrm_omp(neq    ,nequ   ,nad    ,ia    ,ja
     1             ,ad     ,al     ,al      ,m     ,b   ,x 
     2             ,c      ,h      ,r       ,s 
     3             ,tol   ,maxit
     4             ,matvec_csrcr_sym_pm_omp,dot_par_omp
     5             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6             ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7             ,.true.,.true. ,fhist_log,.true.
     8             ,nprcs ,mpi)
c .....................................................................
c
c ...non-overllaping
        else 
          call rpsqrm_omp(neq    ,nequ   ,nad    ,ia    ,ja
     1             ,ad     ,al     ,al      ,m     ,b   ,x 
     2             ,c      ,h      ,r       ,s 
     3             ,tol   ,maxit
     4             ,matvec_csrc_sym_pm_omp,dot_par_omp
     5             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6             ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7             ,.true.,.true. ,fhist_log,.true.
     8             ,nprcs ,mpi)
        endif
c .....................................................................
      endif
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 27/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_MINRES: chama a versao do MINRES                              *    
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
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * y(neq)   - arranjo local de trabalho                               *
c * a(neq)   - arranjo local de trabalho                               *
c * g(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * precond  - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diagonal                                            *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_minres(neq      ,nequ  ,nad   ,ia      ,ja
     1                  ,ad       ,al    ,m     ,b       ,x    
     2                  ,c        ,h     ,r     ,s       ,z  ,y 
     3                  ,a        ,g  
     4                  ,tol      ,maxit ,precond,iparam ,fhist_log
     5                  ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     6                  ,i_xfi    ,i_rcvsi,i_dspli)
      implicit none
c ... mpi
      integer my_id
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c .....................................................................
      integer neq,nequ,nad,neq_doti 
      integer ia(*),ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 c(*),h(*),r(*),s(*),z(*),y(*),a(*),g(*)
c ...
      real*8  tol
      integer maxit,nrestart 
      logical fhist_log
c ... precondicionador
      integer precond,iparam(*)
      real*8 m(*)
c ...
      external dot_par
      external matvec_csrc_sym_pm  
      external ildlt_solv,illt_solv  
c ......................................................................
c
c ...   
      nrestart = 10
c ......................................................................
c
c ...
      if(precond .eq. 1 ) then
        call minres(neq    ,nequ,nad  ,ia  ,ja
     1          ,ad     ,al     ,al     ,b   ,x
     2          ,c      ,h      ,r      ,s   ,z  ,y
     3          ,tol    ,maxit
c ... matvec comum:
     4          ,matvec_csrc_sym_pm,dot_par 
     5          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6          ,i_xfi ,i_rcvsi,i_dspli
     7          ,.true.,.true.,.true.)
c .....................................................................
c
c ...
      else if(precond .eq. 2 .or. precond .eq. 5
     1          .or. precond .eq. 7) then
        call pminres(neq      ,nequ   ,nad   ,ia     ,ja
     1          ,ad     ,al     ,al     ,b      ,x    ,m
     2          ,c      ,h      ,r      ,s      ,z
     3          ,y      ,a      ,g
     4          ,tol    ,maxit  ,nrestart
c ... matvec comum:
     5          ,matvec_csrc_sym_pm,dot_par 
     6          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     7          ,i_xfi ,i_rcvsi,i_dspli
     8          ,.true.,.true.,.true.)
c .....................................................................
      endif
c .....................................................................
c
c ...
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 14/04/2016                                    *
c * Data de modificaco : 23/05/2016                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_BICGSTABCG : chama a versao do Bbicgstab desejada             *    
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
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * y(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * precond  - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diag                                                *
c *            3 - iLDLt                                               *
c *            4 -                                                     *
c *            5 - modulo da diagonal                                  *
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
c **********************************************************************  
      subroutine call_bicgstab(neq      ,nequ  ,nad   ,ia      ,ja
     1                        ,ad       ,al    ,m     ,b       ,x    
     2                        ,c        ,h     ,r     ,s       ,z  ,y   
     3                        ,tol      ,maxit ,precond
     4                        ,my_id    ,neqf1i,neqf2i,neq_doti,i_fmapi
     5                        ,i_xfi    ,i_rcvsi,i_dspli)
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
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 z(*),r(*),s(*),c(*),h(*),y(*)
c ...
      real*8  tol
      integer maxit  
c ... precondicionador
      logical diag
      integer precond
      real*8  m(*)
c ...
      external dot_par
      external matvec_csrc_sym_pm        
c ......................................................................

c ... bicgstab
      if(precond .eq. 1) then
c ...  
        call bicgstab(neq   ,nequ   ,nad    ,ia   ,ja      
     1                ,ad    ,al     ,al    ,b    ,x  
     2                ,c     ,h      ,r     ,s    ,z
     3                ,tol   ,maxit 
c ... matvec comum:
     4                ,matvec_csrc_sym_pm,dot_par 
     5                ,my_id        ,neqf1i  ,neqf2i  ,neq_doti   
     6                ,i_fmapi      ,i_xfi  ,i_rcvsi,i_dspli
     7                ,.true.       ,.true. ,.true.)
c .....................................................................
c
c ... pbicgstab - bicgstab com precondicionador diagonal
      else if(precond .eq. 2 .or. precond .eq. 5) then
c ...  
        call pbicgstab(neq   ,nequ   ,nad    ,ia   ,ja      
     1                ,ad    ,al     ,al     ,m    ,b      ,x  
     2                ,c     ,h      ,r      ,s    ,z      ,y
     3                ,tol   ,maxit 
c ... matvec comum:
     4                ,matvec_csrc_sym_pm,dot_par 
     5                ,my_id        ,neqf1i  ,neqf2i  ,neq_doti   
     6                ,i_fmapi      ,i_xfi  ,i_rcvsi,i_dspli
     7                ,.true.       ,.true. ,.true.)
c .....................................................................
c
c ...  icbicgstab - bicgstab com precondicionador LDLT(0) imcompleto
      elseif(precond .eq. 3 ) then
        call icbicgstab(neq      ,nequ  ,nad   ,ia      ,ja
     1          ,ad       ,al    ,al    ,m     ,b       ,x    
     2          ,c        ,h     ,r     ,s     ,z       ,y
     3          ,tol      ,maxit
c ... matvec comum:
     4          ,matvec_csrc_sym_pm,dot_par 
     5          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6          ,i_xfi ,i_rcvsi,i_dspli
     7          ,.true.,.true. ,.true.)
c .....................................................................
c
c ...
      elseif(precond .eq. 4 ) then
        print*,'iLLt nao implementado para bicgstab!!!'
        stop
      endif  
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 14/04/2016                                    *
c * Data de modificaco : 11/11/2016                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_BICGSTABCGL2 : chama a versao do bicgstabl2 desejada          *    
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
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * y(neq)   - arranjo local de trabalho                               *
c * a(neq)   - arranjo local de trabalho                               *
c * d(neq)   - arranjo local de trabalho                               *
c * e(neq)   - arranjo local de trabalho                               *
c * g(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * precond  - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diag                                                *
c *            3 - iLDLt                                               *
c *            4 -                                                     *
c *            5 - modulo da diagonal                                  *
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
c * mpi      - true|false                                              *
c * ovlp     - overllaping                                             *
c * nprcs    - numero de processos mpi                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_bicgstabl2(neq    ,nequ   ,nad  ,ia        ,ja
     1                        ,ad       ,al     ,m      ,b       ,x    
     2                        ,c        ,h      ,r      ,s       ,z     
     3                        ,y        ,a      ,d      ,e       ,g 
     4                        ,tol      ,maxit  ,precond
     5                        ,my_id    ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                        ,i_xfi    ,i_rcvsi,i_dspli
     7                        ,nprcs    ,ovlp   ,mpi) 
      implicit none
      include 'time.fi'
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
c .....................................................................
      integer neq,nequ,nad,neq_doti 
      integer ia(*),ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 z(*),r(*),s(*),c(*),h(*),y(*),a(*),d(*),g(*),e(*)
c ...
      real*8  tol
      integer maxit  
c ... precondicionador
      logical diag
      integer precond
      real*8  m(*)
c ...
      external dot_par
      external matvec_csrc_sym_pm,matvec_csrcr_sym_pm        
c ......................................................................

c ... bicgstabl2
      if(precond .eq. 1) then
c ...  
        call bicgstabl2(neq     ,nequ   ,nad,ia ,ja 
     1          ,ad      ,al     ,al     ,b  ,x   
     2          ,c       ,h      ,r      ,s  ,z
     3          ,y       ,a 
     4          ,tol     ,maxit  
c ... matvec comum:
     5          ,matvec_csrc_sym_pm,dot_par 
     6          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     7          ,i_xfi ,i_rcvsi,i_dspli
     8          ,.true.,.true. ,.true.)
c .....................................................................
c
c ... pbicgstabl2 - bicgstabl2 com precondicionador diagonal
      else if(precond .eq. 2 .or. precond .eq. 5 ) then
c ... overllaping
         if(ovlp) then
           call pbicgstabl2(neq     ,nequ   ,nad,ia ,ja 
     1          ,ad      ,al     ,al ,m  ,b  ,x   
     2          ,c       ,h      ,r      ,s      ,z      
     3          ,y       ,a      ,d      ,e      ,g      
     4          ,tol     ,maxit  
c ... matvec comum:
     5          ,matvec_csrcr_sym_pm,dot_par 
     6          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     7          ,i_xfi ,i_rcvsi,i_dspli
     8          ,.true.,.true. ,.true.
     9          ,nprcs ,mpi)
c .....................................................................
c
c ...non-overllaping
         else 
           call pbicgstabl2(neq     ,nequ   ,nad,ia ,ja 
     1          ,ad      ,al     ,al ,m  ,b  ,x   
     2          ,c       ,h      ,r      ,s      ,z      
     3          ,y       ,a      ,d      ,e      ,g      
     4          ,tol     ,maxit  
c ... matvec comum:
     5          ,matvec_csrc_sym_pm,dot_par 
     6          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     7          ,i_xfi ,i_rcvsi,i_dspli
     8          ,.true.,.true. ,.true.
     9          ,nprcs ,mpi)
         endif
c .....................................................................
c
c ...
      elseif(precond .eq. 3 ) then
        print*,'iLDLt nao implementado para bicgstab(2)!!!'
        stop
c ...
      elseif(precond .eq. 4 ) then
        print*,'iLLt nao implementado para bicgstab(2)!!!'
        stop
      endif  
c .....................................................................
      return
      end    
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 15/11/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_BICGSTABCGL2_V2 : chama a versao do bicgstabl2 com bloco pu   * 
c * desejada                                                           *    
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
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * y(neq)   - arranjo local de trabalho                               *
c * a(neq)   - arranjo local de trabalho                               *
c * d(neq)   - arranjo local de trabalho                               *
c * e(neq)   - arranjo local de trabalho                               *
c * g(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * precond  - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diag                                                *
c *            3 - iLDLt                                               *
c *            4 -                                                     *
c *            5 - modulo da diagonal                                  *
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
c * Arranjos jat,iat e kat sao utilizados na retrosubstituicao do      *
c * solver iLDLt                                                       *
c **********************************************************************  
      subroutine call_bicgstabl2_v2(neq    ,nequ  ,nad  ,ia        ,ja
     1                        ,ad       ,al    ,m      ,b       ,x    
     2                        ,c        ,h     ,r      ,s       ,z     
     3                        ,y        ,a     ,d      ,e       ,g 
     4                        ,tol      ,maxit ,precond
     5                        ,my_id    ,neqf1i,neqf2i,neq_doti,i_fmapi
     6                        ,i_xfi    ,i_rcvsi,i_dspli)
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
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 z(*),r(*),s(*),c(*),h(*),y(*),a(*),d(*),g(*),e(*)
c ...
      real*8  tol
      integer maxit  
c ... precondicionador
      logical diag
      integer precond
      real*8  m(*)
c ...
      external dot_par
      external matvec_csrc_pm        
c ......................................................................

c ... bicgstabl2
      if(precond .eq. 1) then
c ...  
        call bicgstabl2(neq     ,nequ   ,nad,ia ,ja 
     1          ,ad      ,al     ,al     ,b  ,x   
     2          ,c       ,h      ,r      ,s  ,z
     3          ,y       ,a 
     4          ,tol     ,maxit  
c ... matvec comum:
     6          ,matvec_csrc_pm,dot_par 
     7          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     8          ,i_xfi ,i_rcvsi,i_dspli
     9          ,.true.,.true. ,.true.)
c .....................................................................
c
c ... pbicgstabl2 - bicgstabl2 com precondicionador diagonal
      else if(precond .eq. 2 .or. precond .eq. 5 ) then
         call pbicgstabl2(neq     ,nequ   ,nad,ia ,ja 
     1          ,ad      ,al     ,al ,m  ,b  ,x   
     2          ,c       ,h      ,r      ,s      ,z      
     3          ,y       ,a      ,d      ,e      ,g      
     4          ,tol     ,maxit  
c ... matvec comum:
     5          ,matvec_csrc_pm,dot_par 
     6          ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     7          ,i_xfi ,i_rcvsi,i_dspli
     8          ,.true.,.true. ,.true.)
c .....................................................................
c
c ...
      elseif(precond .eq. 3 ) then
        print*,'iLDLt nao implementado para bicgstab(2)!!!'
        stop
c ...
      elseif(precond .eq. 4 ) then
        print*,'iLLt nao implementado para bicgstab(2)!!!'
        stop
      endif  
c .....................................................................
      return
      end    
c *********************************************************************
c *********************************************************************
      subroutine get_res(u,x,id,fnno,nnode,ndf)
      implicit none
      integer nnode,ndf
      integer id(ndf,*),fnno(*),i,j,k
      real*8 u(ndf,*),x(*)
c ... loop nos
      do i = 1, nnode
        do j = 1, ndf - 1
          k    = id(j,i)
          if( k .gt. 0 ) then
            x(k) = u(j,i)
          endif
        enddo
        if(fnno(i) .eq. 1) then
          k    = id(ndf,i)
          if( k .gt. 0 ) then
            x(k) = u(ndf,i)
          endif
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
c * Data de criacao    : 18/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * SET_PRECOND : escolhe o precondicionandor                          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * macro   - precondicionador escolhido                               *
c * solver  - nao definido                                             *
c * nin     - aqruivo de entrada                                       *
c * my_id   - id do processo do mpi                                    *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * solver  - solver escolhido                                         *
c *         1 - CG                                                     *
c *         2 - GMRES                                                  *
c *         3 -                                                        *
c *         4 - BICGSTAB                                               *
c *         5 -                                                        *
c *         6 - BICGSTABL2                                             *
c *         7 - MINRES                                                 *
c *         8 - PCR                                                    *
c *         9 - SYMMLQ                                                 *
c *        10 - pardiso                                                *
c *        11 - SQRM                                                   *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c ********************************************************************** 
      subroutine set_solver(macro,solver,nin,my_id)
      implicit none
      include 'string.fi'
      include 'precond.fi'
      character macro(maxstrl)
      character*8 macros(12),string
      integer solver,nin,my_id
      integer i,nmc 
      data macros/'cg      ','sqrm    ','symmlq  '
     .           ,'cr      ','minres  ','bicgsl2 '
     .           ,'bicgs   ','gmres   ','pardiso '
     .           ,'block_it','        ','        '/
      data nmc /12/
c ...
      write(string,'(8a)') (word(i),i=1,8)
c ... CG
      if( string .eq. macros(1)) then
        solver = 1
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(1))
        endif
c .....................................................................
c
c ... SQRM
      else if( string .eq. macros(2)) then
        solver = 11
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(2))
        endif
c .....................................................................
c
c ... symmlq
      else if( string .eq. macros(3)) then
        solver =  9
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(3))
        endif
c .....................................................................
c
c ... cr     
      else if( string .eq. macros(4)) then
        solver =  8
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(4))
        endif
c .....................................................................
c
c ... minres
      else if( string .eq. macros(5)) then
        solver =  7
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(5))
        endif
c .....................................................................
c
c ... BICGSTABL2
      else if( string .eq. macros(6)) then
        solver =  6
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(6))
        endif
c .....................................................................
c
c ... BICGSTAB
      else if( string .eq. macros(7)) then
        solver =  4
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(7))
        endif
c .....................................................................
c
c ... GMRES   
      else if( string .eq. macros(8)) then
        solver =  2
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :', adjustr(macros(8))
        endif
c .....................................................................
c
c ... PARDISO 
      else if( string .eq. macros(9)) then
        solver =  10
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :', adjustr(macros(9))
        endif
c .....................................................................
c
c ... BLOK_IT_PCG
      else if( string .eq. macros(10)) then
        solver =  5
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :', adjustr(macros(9))
        endif
c .....................................................................
c
c ...                         
      else
        print*,'Erro na leitura da macro set solver !'
        print*,'Solver disponiveis:'
        do i = 1, nmc
          if(my_id.eq.0) print*,'Solver: ',macros(i)
        enddo
        call stop_mef()
      endif 
c .....................................................................
c
c ...
      return
      end
c ********************************************************************** 
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
