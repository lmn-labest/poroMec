c*****************************Svn***************************************      
c*$Date: 2011-12-02 22:10:15 -0200 (Fri, 02 Dec 2011) $                 
c*$Rev: 959 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************      
      subroutine pcg_omp(neq,ia,ja,ad,au,al,m,b,x,z,r,tol,maxit,
     .              matvec,dot,my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .              i_xfi,i_rcvsi,i_dspli,thread_y)
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
c ... Mpi      
      integer ierr
c .....................................................................      
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c ......................................................................      
      integer neq,maxit,i,j,nad
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),z(*),b(*)
      real*8  dot,tol,conv,energy,d,alpha,beta
      real*8  time0,time
      real*8  dottmp
      real*8  thread_y(*)
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
      nad = ia(neq+1)-1
      if(my_id.eq.0)print *, 'nad :',nad
c$omp parallel private(i,j,d,conv,alpha,beta,dottmp,energy)
!$      thread_id = omp_get_thread_num()
c ......................................................................
c
c ... Chute inicial:
c
c$omp do
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end do
c ----------------------------------------------------------------------
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp do
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         z(i) = r(i) / m(i)
         b(i) = z(i)
  100 continue
c$omp end do
      d    = dot(r(1),z(1),neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do 230 j = 1, maxit
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               b,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
         dottmp = dot(b,z,neq_doti)
         alpha = d / dottmp
c$omp do
         do 210 i = 1, neq
            x(i) = x(i) + alpha * b(i)
            r(i) = r(i) - alpha * z(i)
            z(i) = r(i) / m(i)
  210    continue
c$omp end do
         dottmp = dot(r,z,neq_doti)
         beta = dottmp / d
c$omp do
         do 220 i = 1, neq
            b(i) = z(i) + beta * b(i)
  220    continue
c$omp end do
         d = beta * d
         if (dsqrt(dabs(d)) .lt. conv) goto 300
  230 continue
c ----------------------------------------------------------------------
c$omp master
      write(*,1200) maxit
      call stop_mef()
c$omp end master
  300 continue
c
c ... Energy norm:
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
      energy = dot(x(1),z(1),neq_doti)
c ......................................................................
c$omp single
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100) tol,neq,j,energy,time
c ......................................................................
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10)')"PCG_OMP: ",
     .               "it",j, " energy norm ",energy," tol ",tol
c$omp end single
c$omp end parallel
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Coeficiente da diagonal nulo
     .ou negativo - equacao ',i7)
 1100 format(' (PCG_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
      end
      subroutine gmres_omp(neq,ia,ja,ad,au,al,m,b,x,k,g,h,y,c,s,e,
     .               tol,maxit,matvec,dot,neqovlp,my_id,neqf1i,
     .               neqf2i,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
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
      include 'openmp.fi'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,k,maxit,ia(*),ja(*),neqovlp,nit,i,j,l,ni,ic,nad,nad1
      real*8  ad(neq),au(*),al(*),m(*),b(*),x(*)
      real*8  g(neqovlp,1:k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  energy,econv,norm,dot,ddot,r,aux1,aux2,beta
      real*8  time0,time
      real*8  thread_y(*)
      external matvec,dot
      integer my_id
c      integer nii(maxit)
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
c ...    pre-condicionador diagonal:                           
         g(i,1) = b(i) / m(i)
   10 continue
c$omp end parallel do
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
c$omp parallel private(aux1)
      aux1  = dot(g(1,1),g(1,1),neq_doti)
c$omp single
      econv = tol*dsqrt(aux1)
c$omp end single nowait
c$omp end parallel
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
     .               x,g(1,1),neqf1i,neqf2i,i_fmapi,i_xfi,
     .               i_rcvsi,i_dspli,thread_y)
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
c$omp parallel private(aux1)
         aux1 = dot(g(1,1),g(1,1),neq_doti)
c$omp single
         e(1) = dsqrt(aux1)
c$omp end single
c$omp end parallel
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
            nit = nit + 1
            ni  = ni  + 1
c
c ......... Produto g(i+1) = A.g(i):
c
c$omp parallel
            call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,
     .                  al(nad+1),g(1,i),g(1,i+1),neqf1i,neqf2i,
     .                  i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
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
c$omp parallel private(aux1)            
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
c$omp parallel private(aux1)
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
c$omp parallel do
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
            aux1 = 0.d0
c$omp parallel do reduction(-:aux1)            
            do 510 j = i+1, ni
               aux1 = aux1 - h(i,j)*y(j)
  510       continue
c$omp end parallel do
            y(i) = (aux1 + e(i)) / h(i,i)
  520    continue
c
c ...... Atualizacao de x:
c
         do 610 i = 1, neq
            aux1 = x(i)
c$omp parallel do reduction(+:aux1)
            do 600 j = 1, ni
               aux1 = aux1 + y(j) * g(i,j)
  600       continue
c$omp end parallel do
            x(i) = aux1
  610    continue
c
c ...... Verifica a convergencia:
c
c         nii(l)=ni
         if (dabs(e(ni+1)) .le. econv) goto 1100
 1000 continue
c ----------------------------------------------------------------------
 1100 continue
c
c ... Norma de energia da solucao:
c
c$omp parallel private(aux1)
      aux1 = dot(x(1),b(1),neq_doti)
c$omp single
      energy = aux1
c$omp end single
c$omp end parallel
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,2000) tol,neq,l,nit,dabs(e(ni+1)),energy
     .                           ,time
      if (dabs(e(ni+1)) .gt. econv) then
         write(*,2100) maxit
c         stop
      endif
c ......................................................................
c     Controle de flops
c      if(my_id.eq.0)write(10,'(999(i4,1x))') l,nit,(nii(j),j=1,l)
c ......................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES_OMP) solver:'/
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
      subroutine pbicgstab_omp(neq,ia,ja,ad,au,al,m,b,x,t,v,r,
     .            p,z,tol,maxit,matvec,dot,my_id,neqf1i,neqf2i,
     .            neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c **********************************************************************
c *                                                                    *
c *   Subroutine PBICGSTAB_OMP                                         *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   biconjugados com precondicionador diagonal para matrizes         *
c *   naosimetricas.                                                   *
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
      include 'omp_lib.h'
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
      real*8  time0,time,dottmp1,dottmp2
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
c$omp parallel private(i,j,d,conv,beta,alpha,rr0,w,energy) 
c$omp do
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end do
c ----------------------------------------------------------------------
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp do
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         p(i) = r(i)
         b(i) = p(i)
         z(i) = p(i)/m(i) 
  100 continue
c$omp end do
      d = dot(r(1),z(1),neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do 230 j = 1, maxit
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,v,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
         rr0   = dot(b,r,neq_doti)
         alpha = rr0 / dot(v,r,neq_doti)
c$omp do          
         do 210 i = 1, neq
            x(i) = x(i) + alpha * z(i)
            b(i) = b(i) - alpha * v(i)
            z(i) = b(i) / m(i)
  210    continue
c$omp end do
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,t, neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
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
              z(i) = p(i)/m(i)
  225    continue
c$omp end do  
  230 continue
c ----------------------------------------------------------------------
c$omp single      
      write(*,1200) maxit
      call stop_mef()
c$omp end single      
  300 continue
c
c ... Energy norm:
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
      energy   = dot(x,z,neq_doti)
      
c ......................................................................
c$omp single
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
      end
c **********************************************************************
c
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
      energy = dot(x(1),b(1),neq_doti)
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
