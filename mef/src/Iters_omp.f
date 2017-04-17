c *********************************************************************
c * Metodos iterativos para solucao de sistemas lineares (OpenMP)     *
c * ----------------------------------------------------------------- *
c * simetricos:                                                       *
c * ----------------------------------------------------------------- *
c *                                                                   *
c * PCG_OMP - gradiente conjugados com precondicionador diagonal      *
c *                                                                   *
c * RSQRM_OMP - QRM simetrico com precondicionador diagonal a direta  *
c *                                                                   *    
c * ----------------------------------------------------------------- *
c * nao-simetricos:                                                   *
c * ----------------------------------------------------------------- *
c *                                                                   *                                                       *
c * pbicgstab - gradiente bi-conjugados estabilizados  com            * 
c * precondicionador diagonal                                         *
c *                                                                   *
c * gmres2_omp(m) - GMRES com precondicionador diagonal               *
c *                                                                   *
c * gmres_omp(m) - GMRES com precondicionador diagonal                *
c *                                                                   *
c * ----------------------------------------------------------------- *
c ********************************************************************* 
      subroutine pcg_omp(neq   ,nequ,nad,ia ,ja
     1                  ,ad    ,au  ,al ,m  ,b
     2                  ,x     ,z   ,r  ,p
     3                  ,tol,maxit
     4                  ,matvec,dot
     5                  ,my_id ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     6                  ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7                  ,fprint,flog   ,fhist  ,fnew
     8                  ,nprcs ,mpi)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 15/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * PCG_OMP : Solucao de sistemas de equacoes pelo metodo dos          *
c * gradientes conjugados com precondicionador diagonal para matrizes  *
c * simetricas                                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    - MPI                                                     *
c * neqf1i   - MPI                                                     *
c * neqf2i   - MPI                                                     *
c * neq_doti - MPI                                                     *
c * i_fmap   - MPI                                                     *
c * i_xfi    - MPI                                                     *
c * i_rvcs   - MPI                                                     *
c * i_dspli  - MPI                                                     *
c * thread_y - buffer de equacoes para o vetor y (openmp)              *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fhist    - log dos resuduos por iteracao                           *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
c * mpi      - true|false                                              *
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
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
c ... mpi
      logical mpi        
      integer neqf1i,neqf2i,neq_doti,nprcs,ierr
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................
      integer*8 ia(*),nad   
      integer neq,nequ,maxit,i,j,jj
      integer ja(*),my_id
      real*8  ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8  r(*),z(*),p(*)
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp,norm_b
      real*8  norm_r,norm_m_r
      real*8  time0,time
      real*8  thread_y(*)
      logical flog,fprint,fnew,fhist
c ...
      real*8 flop_cg
      real*8  mflops,vmean
c .....................................................................
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime() 
c......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i,ad(i)
          call stop_mef()
        endif 
   5  continue
c ......................................................................
c
c ......................................................................
c$omp parallel default(none) 
c$omp.private(i,j,jj,conv,alpha,beta,xkx,norm,d,di,tmp)
c$omp.private(norm_b,norm_r,norm_m_r)
c$omp.shared(neq,nequ,nad,ia,ja,al,ad,au,b,x,m,p,z,r,tol,maxit,thread_y)
c$omp.shared(neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,neq_doti)
c$omp.shared(flog,fprint,fnew,fhist,my_id,time,time0,mpi,mflops)
c$omp.shared(vmean,nprcs,ierr)
c$omp.num_threads(nth_solv)                                          
c ......................................................................
c
c ... Chute inicial:
c
      if(fnew) then  
c$omp do
        do 10 i = 1, neq
          x(i) = 0.d0
   10   continue
c$omp end do
      endif  
c ......................................................................
c
c ... conv = tol * |(M-1)b|m = tol *(b,M-1b)
c$omp do
      do 15 i = 1, neq
         z(i) = b(i) * m(i)
   15 continue
c$omp end do
      d      = dot(b,z,neq_doti)
      norm_b = dsqrt(dabs(d))  
      conv   = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0 
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp do
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
c ... z0 = (M-1)r0
         z(i) = r(i) * m(i)
c ... p0 = r0
         p(i) = z(i)
  100 continue
c$omp end do
      d    = dot(r,z,neq_doti)
c ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )
      d    = dot(r,z,neq_doti)
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,p,z
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,thread_y)
c .....................................................................
c
c ... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
c$omp do
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
c ... z  = (M-1)r0
            z(i) = r(i) * m(i)
  210    continue
c$omp end do
c .....................................................................
c
c ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
         di   = dot(r,z,neq_doti) 
c ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
         beta = di / d
c .....................................................................
c
c ...
c$omp do
         do 220 i = 1, neq
c ... p(j+1) = (M-1)r(j+1) + beta*p(j) = z + beta*p(j)
            p(i) = z(i) + beta * p(i)
  220    continue
c$omp end do
c ......................................................................
c
c ...
         if(fhist) write(18,1500),j,dsqrt(dabs(d))/norm_b,sqrt(dabs(d)) 
c .....................................................................
c
c ...
         d =  di
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
c
c ...
c$omp master
         if( jj .eq.1000) then
           jj = 0
           if(my_id .eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c$omp end master
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
c$omp single
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
c$omp end single
  300 continue
c
c ... Energy norm:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,thread_y)
      xkx = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =M(-1)(b - Ax) (calculo do residuo explicito)
c$omp do
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
        z(i) = r(i)*m(i)
  310 continue
c$omp end do
      norm_m_r = dot(r,z,neq_doti)
      norm_m_r = dsqrt(dabs(norm_m_r))
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. conv ) then
c$omp single
         if(my_id .eq.0)then
           write(*,1400) norm_m_r,conv
         endif 
c$omp end single
      endif
c ......................................................................
c$omp single
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
c 
c ...    
      if(mpi) then
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call mpi_mean(vmean,time,nprcs) 
        time   = vmean        
      endif    
      mflops = (flop_cg(neq,nad,j,2,mpi)*1.d-06)/time  
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
       if(mpi) then
          write(*,1110)tol,conv,j,xkx,norm,norm_r,norm_m_r,mflops,time
        else
          write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,norm_m_r
     .                ,mflops,time
        endif
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'PCG_OMP: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
c$omp end single
c$omp end parallel
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG_OMP:',/,5x,'Diagonal coefficient ' 
     . '- equation ',i9)
 1100 format(' (PCG_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1110 format(' (PCG_OMP_MPI) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,d20.10
     .        ' iterations !',/)
 1300 format (' PCG_OMP:',5x,'It',i7,5x,2d20.10)
 1400 format (' PCG_OMP:',1x,'Explicit residual > tol * ||b||| :'
     .       ,1x,d20.10,1x,d20.10)
 1500 format ( 5x,i7,5x,2es20.10)
      end
c *********************************************************************  
c
c *********************************************************************  
      subroutine rpsqmr_omp(neq   ,nequ   ,nad   ,ia  ,ja
     1                 ,ad    ,au     ,al    ,m   ,b  ,x  
     2                 ,t     ,r     ,q   ,d   
     3                 ,tol   ,maxit
     4                 ,matvec,dot
     5                 ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                 ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7                 ,fprint,flog   ,fhist  ,fnew
     8                 ,nprcs ,mpi    ,fsqmr)
c **********************************************************************
c * Data de criacao    : 22/09/2016                                    *
c * Data de modificaco : 02/02/2017                                    * 
c * ------------------------------------------------------------------ *   
c * RPSQRM : Solucao de sistemas de equacoes pelo metodo QMR simetrico *
c * diagonal a direita                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * r(neq)   - arranjo local de trabalho                               *
c * q(neq)   - arranjo local de trabalho                               *
c * t(neq)   - arranjo local de trabalho                               *
c * d(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    - MPI                                                     *
c * neqf1i   - MPI                                                     *
c * neqf2i   - MPI                                                     *
c * neq_doti - MPI                                                     *
c * i_fmap   - MPI                                                     *
c * i_xfi    - MPI                                                     *
c * i_rvcs   - MPI                                                     *
c * i_dspli  - MPI                                                     *
c * thread_y - buffer de equacoes para o vetor y (openmp)              *                                                 *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fhist    - log dos resuduos por iteracao                           *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
c * mpi      - true|false                                              *
c * nprcs    - numero de processos mpi                                 * 
c * fsqmr    - nao definido                                            *  
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * fsqmr    - true para falha do sqmr                                 *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * Fonte: A New Krylov-subspace method for symmetric indefinite       * 
c * linear systems                                                     *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
c ... mpi
      logical mpi        
      integer neqf1i,neqf2i,neq_doti,nprcs,ierr
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,maxit,i,j,jj
      integer*8 ia(*),nad
      integer ja(*),my_id
      real*8 ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8 r(*),t(*),q(*),d(*)
      real*8 dot,tol,conv,xkx,norm,alpha,beta,tmp1,tmp2,tau,ro,vn,v0
      real*8 sigma,cn,norm_b 
      real*8 norm_r
      real*8 time0,time
      real*8 thread_y(*)
      logical flog,fprint,fnew,fhist,fsqmr
c ...
      real*8 flop_sqrm
      real*8  mflops,vmean
c .....................................................................
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
          call stop_mef()
        endif 
   5  continue
c ......................................................................
c
c ......................................................................
c$omp parallel default(none) 
c$omp.private(i,j,jj,conv,alpha,beta,xkx,norm_b,norm_r,norm)  
c$omp.private(tau,ro,v0,sigma,vn,cn,tmp1,tmp2)
c$omp.shared(neq,nequ,nad,ia,ja,al,ad,au,b,x,m,t,r,q,d)
c$omp.shared(tol,maxit,thread_y)
c$omp.shared(neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,neq_doti,flog)
c$omp.shared(my_id,time,time0,fnew,fhist,fprint,mpi,mflops)
c$omp.shared(vmean,nprcs,ierr,fsqmr)
c$omp.num_threads(nth_solv)                                          
c ......................................................................
c
c ... Chute inicial:
c
      if(fnew) then 
c$omp do 
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
c$omp end do
      endif 
c .......................................................................
c
c ... conv = tol * |b| 
      norm_b = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(norm_b))
c .......................................................................
c  
c ... Ax0                                                            
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,t,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,thread_y)
c .......................................................................
c
c ...
c$omp do
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - t(i)
c ... q = t
         q(i) = r(i) * m(i)
c ... d = 0.0
         d(i) = 0.d0
  100 continue
c$omp end do
c ... ( r,r ) 
      tau = dsqrt(dot(r,r,neq_doti))
c ... ( r,q ) 
      ro  = dot(r,q,neq_doti)
c ......................................................................
c
c ...
      v0 = 0.0d0
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... t = Aq(j-1)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,q,t,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,thread_y)
c .....................................................................
c
c ... sigma = ( q(j-1),t)
         sigma = dot(q,t,neq_doti)
         if( sigma .eq. 0.0) then
c$omp single
           print*,"RSQMR_OMP fail (sigma)!"
c$omp end single
           call stop_mef()  
         endif  
c .....................................................................
c
c ... alpha(j-1) = ro(j-1)/sigma(j-1)
         alpha = ro/sigma
c .....................................................................
c
c ...
c$omp do
         do 210 i = 1, neq
c ... r(j) = r(j-1) - alpha(j-1)*t
            r(i) = r(i) - alpha * t(i)
  210    continue
c$omp end do
c .....................................................................
c
c ... v(j) = ||r||/tau(j-1)
         vn   = dsqrt(dot(r,r,neq_doti))/tau
c ... c(j) = 1/sqrt(1+v(j)*v(j) 
         cn   = 1.0d0/dsqrt(1.d0+vn*vn)
c ... tau(j) = tau(j-1)*v(j)*c(j)
         tau = tau*vn*cn
c .....................................................................
c
c ... tau(j) = (c(j)*c(j)*v(j-1)*v(j-1)) d(j-1)
c            + c(j)*c(j)*alpha(j-1) * q(j-1)
         tmp1 = cn*cn*v0*v0
         tmp2 = cn*cn*alpha
c$omp do 
         do 215 i = 1, neq
c ... d(j) = (cj*cj*vj*vj) d(j-1) + cj*cj*alpha(j-1) * q(j-1)
           d(i) = tmp1*d(i) + tmp2*q(i) 
c ... x(j) = x(j-1) + d(j)
           x(i) = x(i) + d(i) 
c ... u(j) = M(-1)r
           t(i) = r(i)*m(i)
  215    continue 
c$omp end do
c .....................................................................
c
c ...
         v0 = vn
c ......................................................................
c
c ... 
         if( ro .eq. 0.0) then
c$omp single
           print*,"RSQMR_OMP fail (ro)!"
c$omp end single
           call stop_mef()  
         endif  
c .....................................................................
c
c ... (r,u) 
         tmp1 = dot(r,t,neq_doti) 
c ... beta = (r,u)/ro(j-1)
         beta = tmp1/ro
c ...
         ro = tmp1
c .....................................................................
c
c ...
c$omp do
         do 220 i = 1, neq
c ... q(j+1) = (M-1)r(j) + beta*q(j-1) = t + beta*q(j-1)
            q(i) = t(i) + beta * q(i)
  220    continue
c$omp end do
c .....................................................................
c
c ...
         norm_r = dsqrt(dot(r,r,neq_doti))
c$omp master
         if(fhist) write(18,1500),j,norm_r/norm_b,norm_r 
c$omp end master
c .....................................................................
c
c ...
         if (norm_r .lt. conv) goto 300
c ......................................................................
c$omp master
         if( jj .eq. 1000) then
           jj = 0
           if(my_id .eq.0 .and. fprint) write(*,1300),j,norm_r,conv 
         endif  
         jj = jj + 1
c$omp end master
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
c$omp master
      if(my_id .eq. 0 .and. fprint ) then 
        write(*,1200) maxit,norm_r
        write(*,*)' *** Switching to Gmres(l) !'
        if(flog) write(10,1200) maxit,norm_r
      endif
c$omp end master
      fsqmr = .true.
      goto 400
  300 continue
c
c ... Energy norm:  x*Kx
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .           ,x,t,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,thread_y)
      xkx = dot(x,t,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =(b - Ax) (calculo do residuo explicito)
c$omp do
      do 310 i = 1, neq
        r(i) = b(i) - t(i)
  310 continue
c$omp end do
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_r .gt. conv ) then
c$omp single
         if(my_id .eq.0  .and. fprint )then
           write(*,1400) norm_r,conv
         endif 
c$omp end single
      endif
c ......................................................................
c$omp single
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
c 
c ...    
      if(mpi) then
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call mpi_mean(vmean,time,nprcs) 
        time   = vmean        
      endif    
      mflops = (flop_sqrm(neq,nad,j,2,mpi)*1.d-06)/time  
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        if(mpi) then
          write(*,1110)tol,conv,j,xkx,norm,norm_r,mflops,time
        else
          write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r
     .                ,mflops,time
        endif
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,3(a,d20.10),2(a,f20.2))')
     .       'RPSQMR: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' Mflops ',mflops,' time ',time
        endif
      endif
c$omp end single
c ......................................................................
c
c ...
  400 continue
c$omp end parallel
c ...
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA LPSMRQ:',/,5x,'Diagonal coefficient ' 
     . '- equacao ',i9)
 1100 format(' (RPSQMR_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||          = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
1110  format(' (RPSQMR_MPI_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||          = ',d20.6/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,d20.10,
     .        ' iterations !',/)
 1300 format (' RPSQMR_OMP:',5x,'It',i7,5x,2d20.10)
 1400 format (' RPSQMR_OMP:',1x,'Explicit residual > tol * ||b|| :'
     .       ,1x,d20.10,1x,d20.10)
 1500 format ( 'RPSQMR_OMP: ',5x,i7,5x,2es20.10)
      end
c **********************************************************************
c
c *********************************************************************  
      subroutine gmres2_omp(neq   ,nequ ,nad,ia,ja
     1                 ,ad    ,au   ,al ,m ,b ,x,k
     2                 ,g     ,h    ,y  ,c ,s ,e
     3                 ,tol   ,maxit
     4                 ,matvec,dot
     5                 ,neqovlp
     6                 ,my_id ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     7                 ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7                 ,fprint,flog   ,fhist  ,fnew   
     8                 ,nprcs ,mpi)
c **********************************************************************
c * Data de criacao    : 29/12/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * GMRES: Solucao iterativa de sistemas simetricos e nao-simetricos   *
c *        pelo metodo GMRES com precondicionador diagonal.            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * k        - base de Krylov                                          *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *]
c * thread_y - buffer de equacoes para o vetor y (openmp)              *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fhist    - log dos resuduos por iteracao                           *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
c * mpi      - true|false                                              *
c * nprcs    - numero de processos mpi                                 *  
c * ------------------------------------------------------------------ * 
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(*),ad(*),al(*),au(*) - inalterados                               *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Arranjos locais de trabalho:                                       *
c *                                                                    *
c * g(neqovlp+1,k+1)                                                   *
c * h(k+1,k)                                                           *
c * y(k)                                                               *
c * c(k)                                                               *
c * s(k)                                                               *
c * e(k+1)                                                             *
c *                                                                    *
c * versao com solver triagular superio (Ry = e) versao coluna         *
c * versao com matriz vetor geral ( x = Vy ) versao coluna             *
c * versao com refined modified gram-schmidt                           *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
c ... mpi
      logical mpi        
      integer neqf1i,neqf2i,neq_doti,nprcs,ierr
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ
      integer*8 ia(*),nad
      integer k,maxit,ja(*),neqovlp,nit,i,j,jj,l,ni,ic,nadr
      real*8  ad(*),au(*),al(*),m(*),b(*),x(*)
      real*8  g(neqovlp,k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  xkx,econv,norm,dot,r,aux1,aux2,beta,inorm,norm_r,norm_m_r
      real*8 tau,kesp
      real*8  time0,time
      real*8 dum1
      real*8  thread_y(*)
      logical flog,fprint,fnew,fhist
c ...
      real*8  flop_gmres
      real*8  mflops,vmean
c .....................................................................
      external matvec,dot
      integer my_id
      parameter (kesp = 0.25d0)
c ......................................................................
      time0 = MPI_Wtime()
c ......................................................................
c
c.... Chute inicial:
c
c$omp parallel private(aux1) num_threads(nth_solv)
c$omp do
      do 10 i = 1, neq
         if(fnew) x(i) = 0.d0
c ...    pre-condicionador diagonal:                  
         g(i,1) = b(i)*m(i)
   10 continue
c$omp enddo
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
      aux1  = dot(g(1,1),g(1,1),neq_doti)
c$omp single
      econv = tol*dsqrt(aux1)
c$omp end single
c$omp end parallel
c ----------------------------------------------------------------------      
c
c ... Ciclos GMRES:
c
      nit = 0
      jj  = 0
      do 1000 l = 1, maxit
c
c ...... Residuo g(1) = b - A x:
c
c$omp parallel private(aux1) num_threads(nth_solv) 
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,x,g(1,1)
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi 
     .              ,i_rcvsi,i_dspli,thread_y)
c
c ...... Residuo com precondicionador diagonal:
c
c$omp do
         do 200 i = 1, neq
            g(i,1) = (b(i) - g(i,1))*m(i)
  200    continue
c$omp enddo
c
c ...... Norma do residuo:
c
         aux1 = dot(g(1,1),g(1,1),neq_doti)
c$omp single
         e(1) = dsqrt(aux1) 
c$omp end single
c
c ...... Normalizacao de g1:
c
         inorm = 1.d0/e(1) 
c$omp do
         do 210 i = 1, neq
            g(i,1) = g(i,1)*inorm
  210    continue
c$omp enddo
c$omp end parallel
c
c ...... Iteracoes GMRES:
c
         ni = 0
         do 400 i = 1, k
            nit = nit + 1
            ni  = ni  + 1
c
c ......... Produto g(i+1) = A.g(i):
c
c$omp parallel num_threads(nth_solv) 
            call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al
     .                 ,al(nad+1)
     .                 ,g(1,i),g(1,i+1)
     .                 ,neqf1i,neqf2i 
     .                 ,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c
c ......... Precondicionador diagonal:
c
c$omp do
            do 300 j = 1, neq
                g(j,i+1) = g(j,i+1)*m(j)
  300       continue
c$omp enddo
c
c ......... Ortogonalizacao (Gram-Schmidt modificado com refinamento):
c
c ...
            aux1 =  dot(g(1,i+1),g(1,i+1),neq_doti)
c$omp single
            tau = dsqrt(aux1) 
c$omp end single
c$omp end parallel
c .....................................................................
            do 320 j = 1, i
c$omp parallel private(aux1) num_threads(nth_solv) 
               aux1 = dot(g(1,i+1),g(1,j),neq_doti)
c$omp single
               beta = aux1
c$omp end single
c$omp do
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
c$omp end do
c$omp end parallel
               h(j,i) = beta
  320       continue
c ......................................................................
c
c ....
c$omp parallel private(aux1) num_threads(nth_solv)
            aux1 = dot(g(1,i+1),g(1,i+1),neq_doti)
c$omp single
            norm = dsqrt(aux1)
c$omp end single
c$omp end parallel
            if( norm .le. kesp*tau) then
              do 321 j = 1, i
c$omp parallel private(aux1) num_threads(nth_solv) 
                  aux1 = dot(g(1,i+1),g(1,j),neq_doti)
c$omp single
                  beta = aux1
c$omp end single
c$omp do
                  do 311 ic = 1, neq
                    g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  311             continue
c$omp end do
c$omp end parallel
                  h(j,i) = h(j,i) + beta
  321          continue
            endif
c ......................................................................
c
c ......... Norma de g(i+1):
c
c$omp parallel private(aux1) num_threads(nth_solv) 
            aux1 = dot(g(1,i+1),g(1,i+1),neq_doti)
c$omp single
            norm = dsqrt(aux1)
c$omp end single
            h(i+1,i) = norm
c
c ......... Normalizacao de g(i+1):
c
            inorm = 1.d0/norm
c$omp do
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)*inorm
  330       continue
c$omp end do
c$omp end parallel
c ...........................................................
c
c ... Givens QR Methods
            do 340 j = 1, i-1
               aux1 =  c(j) * h(j,i) + s(j) * h(j+1,i)
               aux2 = -s(j) * h(j,i) + c(j) * h(j+1,i)
               h(j,i)   = aux1
               h(j+1,i) = aux2
  340       continue
c  
            call sym_ortho2(h(i,i),h(i+1,i),c(i),s(i),r)
c
            h(i,i)   = r
            h(i+1,i) = 0.d0
            e(i+1) = -s(i) * e(i)
            e(i)   =  c(i) * e(i)
            if (dabs(e(i+1)) .le. econv) goto 500
  400    continue
  500    continue
c
c ...... Resolve o sistema r y = e :
c
c ...
         y(1:ni) = e(1:ni)
         do 520 j = ni,2,-1
            y(j) = y(j)/h(j,j)     
            r    = y(j)
            do 510 i = 1 , j - 1
               y(i) = y(i) - h(i,j)*r
  510       continue
  520    continue
         y(1) = y(1)/h(1,1) 
c .....................................................................
c
c ...... Atualizacao de x:
c
         do 600 j = 1, ni
           r = y(j) 
c$omp parallel do num_threads(nth_solv) 
           do 610 i = 1, neq
             x(i) = x(i) + r * g(i,j)
 610       continue
c$omp end parallel do
 600     continue
c ......................................................................
c
c ...
         jj = jj + 1
         if( jj .eq. 10) then
           jj = 0
           if(my_id .eq.0 .and. fprint) then
             write(*,2300),l,nit,dabs(e(ni+1)),econv
           endif
         endif
c ......................................................................
c
c
c ...
         if(fhist) then
           if(my_id .eq.0) write(18,2500)l,dabs(e(ni+1))/norm
     .                          ,dabs(e(ni+1))
         endif  
c .....................................................................
c
c ...... Verifica a convergencia:
c
c         nii(l)=ni
         if (dabs(e(ni+1)) .le. econv) goto 1100
c ......................................................................
 1000 continue
c ......................................................................
 1100 continue
c
c ... Norma da solucao: x*Kx
c
c$omp parallel  num_threads(nth_solv) 
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,
     .           al(nad+1),x     ,g(1,1)  ,neqf1i,neqf2i,
     .           i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
      aux1 = dot(x,g,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      aux2 = dot(x,x,neq_doti)
c$omp single
      xkx  = aux1
      norm = dsqrt(aux2)
c$omp end single  
c ......................................................................
c
c ... r =M(-1)(b - Ax) (calculo do residuo explicito)
c$omp do
      do 1200 i = 1, neq
        g(i,2) = b(i) - g(i,1)
        g(i,3) = g(i,2)*m(i)
 1200 continue
c$omp end do
      aux1 = dot(g(1,2),g(1,2),neq_doti)
      aux2 = dot(g(1,3),g(1,3),neq_doti)
c$omp single
      norm_r   = dsqrt(aux1)
      norm_m_r = dsqrt(aux2)  
      if(  norm_m_r .gt. econv ) then
         if(my_id .eq.0 )then
           write(*,2400)  norm_m_r,econv
         endif 
      endif
c$omp end single
c$omp end parallel  
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if (dabs(e(ni+1)) .gt. econv) then
         if(my_id .eq. 0) then
           write(*,2100) maxit,k,nit
           if(flog) write(10,2100) maxit,k,nit
         endif 
         call stop_mef()
      endif
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        if(mpi) then
          write(*,2010) tol,econv,l,nit,dabs(e(ni+1))
     .                 ,xkx,norm,norm_r,norm_m_r,time
        else
          write(*,2000) tol,econv,neq,l,nit,dabs(e(ni+1))
     .                 ,xkx,norm,norm_r,norm_m_r,time
        endif
      endif
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,i9,a,d20.10,a,f20.2)')
     .         'GMRES2: ',' it ',nit, ' x * Kx ',xkx,' ||x|| ',norm,
     .         ' nKylov ',k,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ......................................................................
 2000 format(' (GMRES2_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.10/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2010 format(' (GMRES2_MPI_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.10/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for '
     .      ,i9,' cycles !',5x,i7,' nKylov',5x,' It ',i7/)
 2300 format (' GMRES2:',5x,'cycles',i7,5x,'It',i7,5x,2d20.10)
 2400 format (' GMRES2:',1x,'Explicit residual > tol * ||b||| :'
     .       ,1x,d20.10,1x,d20.10)
 2500 format ( 5x,i7,5x,2es20.10)
      end
c **********************************************************************
c
c **********************************************************************
      subroutine gmres_omp(neq    ,nequ ,nad,ia ,ja
     .                    ,ad     ,au   ,al ,m  ,b
     .                    ,x      ,k    ,g  ,h  ,y
     .                    ,c      ,s    ,e  ,tol,maxit
     .                    ,matvec ,dot
     .                    ,neqovlp
     .                    ,my_id  ,neqf1i,neqf2i,neq_doti,i_fmapi
     .                    ,i_xfi  ,i_rcvsi,i_dspli,thread_y,flog)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 12/12/2015                                    * 
c * ------------------------------------------------------------------ *   
c * GMRES_OMP:Solucao iterativa de sistemas simetricos e nao-simetricos*
c * pelo metodo GMRES com precondicionador diagonal.                   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * k        - base de Krylov                                          *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * thread_y - buffer de equacoes para o vetor y (openmp)              *
c * flog     - log do arquivo de saida                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq),ad(*),al(*),au(*) - modificados                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Arranjos locais de trabalho:                                       *
c *                                                                    *
c * g(neq+1,k+1)                                                       *
c * h(k+1,k)                                                           *
c * y(k)                                                               *
c * c(k)                                                               *
c * s(k)                                                               *
c * e(k+1)                                                             *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,k,maxit,ia(*),ja(*)
      integer neqovlp,nit,i,j,jj,l,ni,ic,nad,nad1
      real*8  ad(neq),au(*),al(*),m(*),b(*),x(*)
      real*8  g(neqovlp,1:k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  energy,econv,norm,dot,ddot,r,aux1,aux2,beta
      real*8  time0,time
      real*8  thread_y(*)
      logical flog
      external matvec,dot
      integer my_id
c ......................................................................
      time0 = MPI_Wtime()
c ......................................................................
c
c.... Chute inicial:
c
c$omp parallel do num_threads(nth_solv)
      do 10 i = 1, neq
        x(i) = 0.d0
c ...    pre-condicionador diagonal:                           
        g(i,1) = b(i) * m(i)
   10 continue
c$omp end parallel do
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
c$omp parallel private(aux1) num_threads(nth_solv) 
      aux1  = dot(g(1,1),g(1,1),neq_doti)
c$omp single
      econv = tol*dsqrt(aux1)
c$omp end single
c$omp end parallel
c ----------------------------------------------------------------------
c
c ... Ciclos GMRES:
c
      nit = 0
      jj  = 0
      do 1000 l = 1, maxit
c
c ...... Residuo g(1) = b - A x:
c
c$omp parallel num_threads(nth_solv) 
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,x,g(1,1)
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,thread_y)
c$omp end parallel
c
c ...... Residuo com precondicionador diagonal:
c
c$omp parallel do num_threads(nth_solv) 
         do 200 i = 1, neq
            g(i,1) = (b(i) - g(i,1))*m(i)
  200    continue
c$omp end parallel do
c
c ...... Norma do residuo:
c
c$omp parallel private(aux1) num_threads(nth_solv) 
         aux1 = dot(g(1,1),g(1,1),neq_doti)
c$omp single
         e(1) = dsqrt(aux1)
c$omp end single
c$omp end parallel
c
c ...... Normalizacao de g1:
c
c$omp parallel do num_threads(nth_solv) 
         do 210 i = 1, neq
            g(i,1) = g(i,1)/e(1)
  210    continue
c$omp end parallel do
c
c ...... Iteracoes GMRES:
c
         ni = 0
         do 400 i = 1, k
            nit = nit + 1
            ni  = ni  + 1
c
c ......... Produto g(i+1) = A.g(i):
c
c$omp parallel num_threads(nth_solv) 
            call matvec(neq,nequ
     .                 ,ia     ,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .                 ,g(1,i) ,g(1,i+1)
     .                 ,neqf1i ,neqf2i
     .                 ,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
c
c ......... Precondicionador diagonal:
c
c$omp parallel do num_threads(nth_solv) 
            do 300 j = 1, neq
               g(j,i+1) = g(j,i+1)*m(j)
  300       continue
c$omp end parallel do
c
c ......... Ortogonalizacao (Gram-Schmidt modificado):
c
            do 320 j = 1, i
c$omp parallel private(aux1) num_threads(nth_solv) 
               aux1 = dot(g(1,i+1),g(1,j),neq_doti)
c$omp single
               beta = aux1
c$omp end single
c$omp do 
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
c$omp end do
c$omp end parallel
               h(j,i) = beta
  320       continue
c
c ......... Norma de g(i+1):
c
c$omp parallel private(aux1) num_threads(nth_solv) 
            aux1 = dot(g(1,i+1),g(1,i+1),neq_doti)
c$omp single
            norm = dsqrt(aux1)
c$omp end single
c$omp end parallel
c
            h(i+1,i) = norm
c
c ......... Normalizacao de g(i+1):
c
c$omp parallel do num_threads(nth_solv) 
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)/norm
  330       continue
c$omp end parallel do
c
            do 340 j = 1, i-1
               aux1 =  c(j) * h(j,i) + s(j) * h(j+1,i)
               aux2 = -s(j) * h(j,i) + c(j) * h(j+1,i)
               h(j,i)   = aux1
               h(j+1,i) = aux2
  340       continue
            r = dsqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
            c(i) = h(i,i)/r
            s(i) = h(i+1,i)/r
            h(i,i)   = r
            h(i+1,i) = 0.d0
            e(i+1) = -s(i) * e(i)
            e(i)   =  c(i) * e(i)
            if (dabs(e(i+1)) .le. econv) goto 500
  400    continue
  500    continue
c
c ...... Resolve o sistema h y = e :
c
         y(ni) = e(ni) / h(ni,ni)
         do 520 i = ni-1, 1, -1
            y(i) = 0.d0
            do 510 j = i+1, ni
               y(i) = y(i) - h(i,j)*y(j)
  510       continue
            y(i) = ( y(i) + e(i)) / h(i,i)
  520    continue
c
c ...... Atualizacao de x:
c
c$omp parallel do num_threads(nth_solv) 
         do 610 i = 1, neq
            do 600 j = 1, ni
               x(i) = x(i) + y(j) * g(i,j)
  600       continue
  610    continue
c$omp end parallel do
c ......................................................................
c
c ...
         jj = jj + 1
         if( jj .eq. 10) then
           jj = 0
           write(*,2300),l,nit,dabs(e(ni+1)),econv
         endif
c ......................................................................
c
c ...... Verifica a convergencia:
c
         if (dabs(e(ni+1)) .le. econv) goto 1100
 1000 continue
c ----------------------------------------------------------------------
 1100 continue
c
c ... Norma de energia da solucao:
c
c$omp parallel  num_threads(nth_solv) 
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .           ,x     ,g(1,1)  
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,thread_y)
      aux1 = dot(x,g(1,1),neq_doti)
c$omp single
      energy = aux1
c$omp end single
c$omp end parallel
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if (dabs(e(ni+1)) .gt. econv) then
         if(my_id .eq. 0) then
            write(*,2100) maxit,k,nit
           if(flog) write(10,2100) maxit,k,nit
         endif 
         call stop_mef()
      endif
c .....................................................................
      if(my_id.eq.0)write(*,2000) tol,neq,l,nit,dabs(e(ni+1)),energy
     .                           ,time
c ......................................................................
      if(flog) then
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10,a,f20.2)')
     .         "GMRES_OMP: "," it ",nit, " energy norm ",energy
     .        ," tol ",tol," time ",time
      endif
c .....................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.10/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for '
     .      ,i9,' cycles !',5x,i7,' nKylov',5x,' It ',i7/)
 2300 format (' GMRES_OMP:',5x,'cycles',i7,5x,'It',i7,5x,2d20.10)
      end
c **********************************************************************
      subroutine pbicgstab_omp(neq   ,nequ   ,nad    ,ia,ja
     .                        ,ad    ,au     ,al     ,m ,b
     .                        ,x     ,t      ,v      ,r ,p ,z
     .                        ,tol   ,maxit
     .                        ,matvec,dot
     .                        ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                        ,i_xfi ,i_rcvsi,i_dspli,thread_y,flog)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 12/12/2015                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTAB_OMP : Solucao de sistemas de equacoes pelo metodo dos    * 
c * gradientes biconjugados com precondicionador diagonal para         *
c * matrizes nao-simetricas.                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                   *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * t(neq)   - arranjo local de trabalho                               *
c * v(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
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
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
c ... Mpi      
      integer ierr
c .....................................................................      
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli      
c .....................................................................      
      integer neq,nequ,maxit,nad,i,j,jj,k
      integer ia(*),ja(*),my_id
      real*8 ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),t(*),v(*),z(*)
      real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
      real*8  time0,time,dottmp1,dottmp2
      real*8 thread_y(*) 
      logical flog
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
c$omp parallel default(none) 
c$omp.private(i,j,jj,d,conv,beta,alpha,rr0,w,energy) 
c$omp.shared(neq,nequ,nad,ia,ja,al,ad,au)
c$omp.shared(m,x,r,p,t,v,z,tol,maxit,thread_y)
c$omp.shared(neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,neq_doti,flog)
c$omp.shared(my_id,time,time0)
c$omp.num_threads(nth_solv)                                          
c
c ... Chute inicial:
c
c$omp do
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end do
c ----------------------------------------------------------------------
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp do
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         p(i) = r(i)
         b(i) = p(i)
         z(i) = p(i)*m(i) 
  100 continue
c$omp end do
      d = dot(r,z,neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      jj = 1
      do 230 j = 1, maxit
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,z,v
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,thread_y)
         rr0   = dot(b,r,neq_doti)
         alpha = rr0 / dot(v,r,neq_doti)
c$omp do          
         do 210 i = 1, neq
            x(i) = x(i) + alpha * z(i)
            b(i) = b(i) - alpha * v(i)
            z(i) = b(i) * m(i)
  210    continue
c$omp end do
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,z,t
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,thread_y)
         w = dot(t,b,neq_doti)/dot(t,t,neq_doti)
c$omp do         
         do 220 i = 1, neq
             x(i) = x(i) + w*z(i)
             b(i) = b(i) - w*t(i)
  220    continue
c$omp end do
         d = dot(b,z,neq_doti)
         if (dsqrt(dabs(d)) .lt. conv) goto 300
         beta = (dot(r,b,neq_doti) / rr0)*(alpha/w)
c$omp do         
         do 225 i = 1, neq
              p(i) = b(i) + beta*(p(i)-w*v(i))
              z(i) = p(i)*m(i)
  225    continue
c$omp end do  
c ......................................................................
c$omp master
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c$omp end master
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
c$omp single      
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
c$omp end single      
  300 continue
c
c ... Energy norm:
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,thread_y)
      energy   = dot(x,z,neq_doti)
c ......................................................................
c$omp single
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100) tol,neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10,a,f20.2)')
     .               "PBICGSTAB_OMP: ","it",j, " energy norm ",energy,
     .               " tol ",tol," time ",time
c ......................................................................
c$omp end single
c$omp end parallel
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PBICGSTAB:',/,5x,'Coeficiente da diagonal
     . nulo ou negativo - equacao ',i7)
 1100 format(' (PBICGSTAB_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BICGSTAB:',5x,'It',i7,5x,2d20.10)
      end
c **********************************************************************
      subroutine pcg_omp_loopwise(neq,ia,ja,ad,au,al,m,b,x,z,r,tol,
     .              maxit,matvec,dot,my_id,neqf1i,neqf2i,neq_doti,
     .              i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c **********************************************************************
c *                                                                    *
c *   Subroutine PCG                                                   *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   conjugados com precondicionador diagonal para matrizes           *
c *   simetricas.                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   energy - nao definido                                            *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - inalterados                                  *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................
      integer neq,maxit,i,j,nad
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),z(*),b(*)
      real*8  dot,tol,conv,energy,d,alpha,beta
      real*8  time0,time
      real*8  dottmp
      real*8 thread_y(*)
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
      nad = ia(neq+1)-1
      if(my_id.eq.0)print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
c$omp parallel do
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end parallel do
c ----------------------------------------------------------------------
c$omp parallel
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
c$omp parallel do
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         z(i) = r(i) / m(i)
         b(i) = z(i)
  100 continue
c$omp end parallel do
      d    = dot(r(1),z(1),neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do 230 j = 1, maxit
c$omp parallel
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               b,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
c$omp end parallel
         dottmp = dot(b,z,neq_doti)
         alpha = d / dottmp
c$omp parallel do
         do 210 i = 1, neq
            x(i) = x(i) + alpha * b(i)
            r(i) = r(i) - alpha * z(i)
            z(i) = r(i) / m(i)
  210    continue
c$omp end parallel do
         dottmp = dot(r,z,neq_doti)
         beta = dottmp / d
c$omp parallel do
         do 220 i = 1, neq
            b(i) = z(i) + beta * b(i)
  220    continue
c$omp end parallel do
         d = beta * d
         if (dsqrt(dabs(d)) .lt. conv) goto 300
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
c$omp parallel
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
      energy = dot(x(1),z(1),neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100) tol,neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10)')
     .               "PCG_OMP_LOOPWISE: ",
     .               "it",j, " energy norm ",energy," tol ",tol
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Coeficiente da diagonal nulo
     .ou negativo - equacao ',i7)
 1100 format(' (PCG_OMP_LOOPWISE) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
      end
      subroutine gmres_omp_loopwise(neq,ia,ja,ad,au,al,m,b,x,k,g,h,y,c,
     .             s,e,tol,maxit,matvec,dot,neqovlp,my_id,
     .             neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .             thread_y) 
c **********************************************************************
c *                                                                    *
c *   GMRES: Solucao iterativa de sistemas simetricos e nao-simetricos *
c *          pelo metodo GMRES com precondicionador diagonal.          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Arranjos locais de trabalho:                                     *
c *                                                                    *
c *      g(neq+1,k+1)                                                  *
c *      h(k+1,k)                                                      *
c *      y(k)                                                          *
c *      c(k)                                                          *
c *      s(k)                                                          *
c *      e(k+1)                                                        *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq),ad(*),al(*),au(*) - modificados                           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      include 'openmp.fi'
      integer neq,k,maxit,ia(*),ja(*),neqovlp,nit,i,j,l,ni,ic,nad,nad1
      real*8  ad(neq),au(*),al(*),m(*),b(*),x(*)
      real*8  g(neqovlp,1:k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  energy,econv,norm,dot,r,aux1,aux2,beta, gnorm
      real*8  time0,time
      real*8 thread_y(*)
      external matvec,dot
      integer nii(maxit),my_id
      integer ns
c ......................................................................
      time0 = MPI_Wtime()
c ......................................................................
c     call omp_set_num_threads(num_threads)
      nad = ia(neq+1)-1
      if(my_id.eq.0)print*,nad
c ----------------------------------------------------------------------
c
c.... Chute inicial:
c
c$omp parallel do
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end parallel do
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
      norm  = dsqrt(dot(b(1),b(1),neq_doti))
      econv = tol*norm
c ----------------------------------------------------------------------
c
c ... Ciclos GMRES:
c
      nit = 0
      do 1000 l = 1, maxit
c
c ...... Residuo g(1) = b - A x:
c
c$omp parallel
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               x,g(1,1),neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .               i_dspli,thread_y)
c$omp end parallel
c
c ...... Residuo com precondicionador diagonal:
c
c$omp parallel do
         do 200 i = 1, neq
            g(i,1) = (b(i) - g(i,1))/m(i)
  200    continue
c$omp end parallel do
c
c ...... Norma do residuo:
c
         e(1) = dsqrt(dot(g(1,1),g(1,1),neq_doti))
c
c ...... Normalizacao de g1:
c
c$omp parallel do
         do 210 i = 1, neq
            g(i,1) = g(i,1)/e(1)
  210    continue
c$omp end parallel do
c
c ...... Iteracoes GMRES:
c
         ni = 0
         do 400 i = 1, k
            nit =nit + 1
            ni  = ni  + 1
c
c ......... Produto g(i+1) = A.g(i):
c
c$omp parallel
            call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,
     .                  al(nad+1),g(1,i),g(1,i+1),neqf1i,neqf2i,i_fmapi,
     .                  i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
c
c ......... Precondicionador diagonal:
c
c$omp parallel do
            do 300 j = 1, neq
               g(j,i+1) = g(j,i+1)/m(j)
  300       continue
c$omp end parallel do
c
c ......... Ortogonalizacao (Gram-Schmidt modificado):
c
            do 320 j = 1, i
               beta = dot(g(1,i+1),g(1,j),neq_doti)
c$omp parallel do
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
c$omp end parallel do
               h(j,i) = beta
  320       continue
c
c ......... Norma de g(i+1):
c
            norm = dsqrt(dot(g(1,i+1),g(1,i+1),neq_doti))
c
            h(i+1,i) = norm
c
c ......... Normalizacao de g(i+1):
c
c$omp parallel do
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)/norm
  330       continue
c$omp end parallel do
  400    continue
c
c ..... Givens QR factorization of the Hessenberg matrix
c
         r = dsqrt(h(1,1)*h(1,1) + h(2,1)*h(2,1))
         c(1)   = h(1,1)/r
         s(1)   = h(2,1)/r
         h(1,1) = r
         h(2,1) = 0.d0
         e(2) = -s(1) * e(1)
         e(1) =  c(1) * e(1)
         ns = 1
         if (dabs(e(2)) .le. econv) goto 500
         do 450 i = 2, k
c$omp parallel do
            do 340 j = i, k
               aux1 =  c(i-1) * h(i-1,j) + s(i-1) * h(i,j)
               aux2 = -s(i-1) * h(i-1,j) + c(i-1) * h(i,j)
               h(i-1,j)   = aux1
               h(i,j)     = aux2
  340       continue
c$omp end parallel do
c ......... New Givens matrix elements
            r = dsqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
            c(i) = h(i,i)/r
            s(i) = h(i+1,i)/r
            h(i,i)   = r
            h(i+1,i) = 0.d0
            e(i+1) = -s(i) * e(i)
            e(i)   =  c(i) * e(i)
            ns = i
            if (dabs(e(i+1)) .le. econv) goto 500
  450    continue
  500    continue
c
c ...... Resolve o sistema h y = e :
c
         y(ns) = e(ns) / h(ns,ns)
         do 520 i = ns-1, 1, -1
            y(i) = 0.d0
            do 510 j = i+1, ns
               y(i) = y(i) - h(i,j)*y(j)
  510       continue
            y(i) = (y(i) + e(i)) / h(i,i)
  520    continue
c
c ...... Atualizacao de x:
c
c TODO: junior dense matvec with columns arrangement
c$omp parallel do
         do 610 i = 1, neq
            do 600 j = 1, ns
               x(i) = x(i) + y(j) * g(i,j)
  600       continue
  610    continue
c$omp end parallel do
c
c ...... Verifica a convergencia:
c
         nii(l)=ni
         if (dabs(e(ns+1)) .le. econv) goto 1100
 1000 continue
c ----------------------------------------------------------------------
 1100 continue
c
c ... Norma de energia da solucao:
c
      energy = dot(x,b,neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,2000) tol,neq,l,nit,dabs(e(ns+1)),energy
     .                           ,time
      if (dabs(e(ns+1)) .gt. econv) then
         write(*,2100) maxit
c         stop
      endif
c ......................................................................
c     Controle de flops
      if(my_id.eq.0)write(10,'(999(i4,1x))') l,nit,(nii(j),j=1,l)
c ......................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES_OMP_LOOPWISE) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.6/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for ',i9,' cycles !',
     . /)
      end
c **********************************************************************
      subroutine pbicgstab_omp_loopwise(neq,ia,ja,ad,au,al,m,b,x,t,v,r,
     .            p,z,tol,maxit,matvec,dot,my_id,neqf1i,neqf2i,
     .            neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c **********************************************************************
c *                                                                    *
c *   Subroutine PBICGSTAB                                             *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   biconjugados com precondicionador diagonal para matrizes         *
c *   simetricas.                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   t(neq) - arranjo local de trabalho                               *
c *   v(neq) - arranjo local de trabalho
c *   r(neq) - arranjo local de trabalho                               *
c *   p(neq) - arranjo local de trabalho
c *   z(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   energy - nao definido                                            *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - inalterados                                  *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
c ... Mpi      
      integer ierr
c .....................................................................      
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli      
c .....................................................................      
      integer neq,maxit,nad,i,j,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),t(*),v(*),z(*)
      real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
      real*8  time0,time
      real*8  thread_y(*)
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
      nad = ia(neq+1)-1
      if(my_id.eq.0)print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
c$omp parallel do private(i)
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end parallel do
c ----------------------------------------------------------------------
c$omp parallel
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
c
c$omp parallel do private(i)
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         p(i) = r(i)
         b(i) = p(i)
         z(i) = p(i)/m(i) 
  100 continue
c$omp end parallel do
      d    = dot(r(1),z(1),neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do 230 j = 1, maxit
c$omp parallel      
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,v,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
c$omp end parallel     
         rr0 = dot(b,r,neq_doti)
         alpha = rr0/dot(v,r,neq_doti)
c$omp parallel do private(i)         
         do 210 i = 1, neq
            x(i) = x(i) + alpha * z(i)
            b(i) = b(i) - alpha * v(i)
            z(i) = b(i) / m(i)
  210    continue
c$omp end parallel do   
c$omp parallel  
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,t, neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
c$omp end parallel     
         w = dot(t,b,neq_doti) / dot(t,t,neq_doti)
c$omp parallel do private(i)        
         do 220 i = 1, neq
            x(i) = x(i) + w*z(i)
            b(i) = b(i) - w*t(i)
  220    continue
c$omp end parallel do
         d = dot(b,z,neq_doti)
         if (dsqrt(dabs(d)) .lt. conv) goto 300   
         beta = (dot(r,b,neq_doti) / rr0)*(alpha/w)
c$omp parallel do private(i)         
         do 225 i = 1, neq
             p(i) = b(i) + beta*(p(i)-w*v(i))
             z(i) = p(i)/m(i)
  225    continue
c$omp end parallel do  

  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
c$omp parallel
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
      energy   = dot(x(1),z(1),neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100) tol,neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10)')
     .               "PBICGSTAB_OMP: ","it",j, " energy norm ",energy,
     .               " tol ",tol
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PBICGSTAB:',/,5x,'Coeficiente da diagonal
     . nulo ou negativo - equacao ',i7)
 1100 format(' (PBICGSTAB_OMP_LOOPWISE) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
      end
c *********************************************************************

