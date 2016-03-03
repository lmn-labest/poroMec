      subroutine pcg(neq   ,nequ  ,nad   ,ia      ,ja
     .              ,ad    ,au    ,al    ,m        ,b      
     .              ,x     ,z     ,r     ,tol      ,maxit
     .              ,matvec,dot
     .              ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .              ,i_xfi ,i_rcvsi,i_dspli
     .              ,flog)
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
c *   nequ   - numero de equacoes no bloco Kuu                         *
c *   nad    - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
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
c *   flog   - log do arquivo de saida                                 *
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
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,nad,maxit,i,j,jj
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),z(*),b(*)
      real*8  dot,ddot,tol,conv,energy,d,alpha,beta
      real*8  time0,time
      real*8 dum1
      logical flog
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
c
c ... Chute inicial:
c
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c ----------------------------------------------------------------------
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1), 
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         z(i) = r(i) * m(i)
         b(i) = z(i)
  100 continue
      d    = dot(r,z,neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      jj = 1
      do 230 j = 1, maxit
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .               b,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               dum1)
         alpha = d / dot(b,z,neq_doti)
         do 210 i = 1, neq
            x(i) = x(i) + alpha * b(i)
            r(i) = r(i) - alpha * z(i)
            z(i) = r(i) * m(i)
  210    continue
         beta = dot(r,z,neq_doti) / d
         do 220 i = 1, neq
            b(i) = z(i) + beta * b(i)
  220    continue
         d = beta * d
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .           x,z,
     .           neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      energy = dot(x,z,neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100)tol,neq,nad,j,energy,time
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10,a,f20.2)')
     .             "PCG: "," it ",j, " energy norm ",energy," tol ",tol,
     .             " time ",time
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Coeficiente da diagonal nulo
     .ou negativo - equacao ',i7)
 1100 format(' (PCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PCG:',5x,'It',i7,5x,2d20.10)
      end
      subroutine gmres(neq,nequ,nad,ia,ja,ad,au,al,m,b,x,k,g,h,y,c,s,e,
     .              tol,maxit,matvec,dot,neqovlp,my_id,neqf1i,neqf2i,
     .              neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,flog)
c **********************************************************************
c *                                                                    *
c *   GMRES: Solucao iterativa de sistemas simetricos e nao-simetricos *
c *          pelo metodo GMRES com precondicionador diagonal.          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   nequ   - numero de equacoes no bloco Kuu                         *
c *   nad    - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
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
c *   flog   - log do arquivo de saida                                 *
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
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,nad
      integer k,maxit,ia(*),ja(*),neqovlp,nit,i,j,jj,l,ni,ic,nad1
      real*8  ad(*),au(*),al(*),m(*),b(*),x(*)
      real*8  g(neqovlp,1:k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  energy,econv,norm,dot,ddot,r,aux1,aux2,beta
      real*8  time0,time
      real*8 dum1
      logical flog
      external matvec,dot
      integer my_id
c      integer nii(maxit)
c ......................................................................
      time0 = MPI_Wtime()
c ......................................................................
c     nad = ia(neq+1)-1
c     if(my_id.eq.0)print*,nad
c ----------------------------------------------------------------------
c
c.... Chute inicial:
c
      do 10 i = 1, neq
         x(i) = 0.d0
c ...    pre-condicionador diagonal:                  
         g(i,1) = b(i)*m(i)
   10 continue
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
      norm  = dsqrt(dot(g(1,1),g(1,1),neq_doti))
      econv = tol*norm
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
c        call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
c    .               x,g(1,1),neqf1i,neqf2i,i_fmapi,i_xfi,
c    .               i_rcvsi,i_dspli,dum1)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .               x,g(1,1),neqf1i,neqf2i,i_fmapi,i_xfi,
     .               i_rcvsi,i_dspli,dum1)
c
c ...... Residuo com precondicionador diagonal:
c
         do 200 i = 1, neq
            g(i,1) = (b(i) - g(i,1))*m(i)
  200    continue
c
c ...... Norma do residuo:
c
         e(1) = dsqrt(dot(g(1,1),g(1,1),neq_doti))
c
c ...... Normalizacao de g1:
c
         do 210 i = 1, neq
            g(i,1) = g(i,1)/e(1)
  210    continue
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
c           call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,
c    .                  al(nad+1),g(1,i),g(1,i+1),neqf1i,neqf2i,
c    .                  i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
            call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,
     .                  al(nad+1),g(1,i),g(1,i+1),neqf1i,neqf2i,
     .                  i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c
c ......... Precondicionador diagonal:
c
            do 300 j = 1, neq
                g(j,i+1) = g(j,i+1)*m(j)
  300       continue
c
c ......... Ortogonalizacao (Gram-Schmidt modificado):
c
            do 320 j = 1, i
               beta = dot(g(1,i+1),g(1,j),neq_doti)
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
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
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)/norm
  330       continue
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
            y(i) = (y(i) + e(i)) / h(i,i)
  520    continue
c
c ...... Atualizacao de x:
c
         do 610 i = 1, neq
            do 600 j = 1, ni
               x(i) = x(i) + y(j) * g(i,j)
  600       continue
  610    continue
c
c ...... Verifica a convergencia:
c
c         nii(l)=ni
         if (dabs(e(ni+1)) .le. econv) goto 1100
c ......................................................................
         jj = jj + 1
         if( jj .eq. 200) then
           jj = 0
           write(*,2300),l,nit,dabs(e(ni+1)),econv
         endif
c ......................................................................
 1000 continue
c ......................................................................
 1100 continue
c
c ... Norma de energia da solucao:
c
      energy = dot(x,b,neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if (dabs(e(ni+1)) .gt. econv) then
         if(my_id .eq. 0) then
           write(*,2100) maxit
           if(flog) write(10,2100) maxit
         endif 
         call stop_mef()
      endif
c ......................................................................
      if(my_id.eq.0)write(*,2000) tol,neq,l,nit,dabs(e(ni+1)),energy
     .                           ,time
c ......................................................................
c     Controle de flops
c      if(my_id.eq.0)write(10,'(999(i4,1x))') l,nit,(nii(j),j=1,l)
      if(flog) then
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10,a,f20.2)')
     .         "GMRES: "," it ",nit, " energy norm ",energy," tol ",tol,
     .         " time ",time
      endif
c ......................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.6/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for ',i9,' cycles !',
     . /)
 2300 format (' GMRES:',5x,'cycles',i7,5x,'It',i7,5x,2d20.10)
      end      
      subroutine bicgstab(neq,ia,ja,ad,au,al,m,b,x,y,z,p,r,s,tol,maxit,
     .                    matvec,dot,my_id,neqf1i,neqf2i,neq_doti,
     .                    i_fmapi,i_xfi,i_rcvsi,i_dspli)
c **********************************************************************
c *                                                                    *
c *   Subroutine BICGSTAB                                              *
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
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - modificado                                   *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................
      integer neq,maxit,nad,i,j,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),y(*),z(*),s(*)
      real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
      real*8  time0,time
      real*8 dum1
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
      do 10 i = 1, neq
         x(i) = 0.d0
c ...... pre-condicionador diagonal:         
         b(i)  = b(i)/m(i)
         ad(i) = ad(i)/m(i)
         do 5 k = ia(i), ia(i+1)-1
            j = ja(k)
            al(k) = al(k) / m(i)
            au(k) = au(k) / m(j)
   5     continue      
  10  continue
c ----------------------------------------------------------------------
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,r,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)    
      do 100 i = 1, neq
         r(i) = b(i) - r(i)
         p(i) = r(i)
         b(i) = r(i)
  100 continue
      d    = dot(r(1),r(1),neq)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do 230 j = 1, maxit
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               dum1)
         rr0 = dot(r,b,neq)
         alpha = rr0/dot(z,b,neq)
          do 210 i = 1, neq
            s(i) = r(i) - alpha * z(i)
  210    continue
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               s,y, neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               dum1)
         w = dot(y,s,neq) / dot(y,y,neq)
         do 220 i = 1, neq
            x(i) = x(i) + alpha*p(i) + w * s(i)
            r(i) = s(i) - w*y(i)
  220    continue
         beta = (dot(r,b,neq) / rr0)*(alpha/w)
         do 225 i = 1, neq
             p(i) = r(i) + beta*(p(i)-w*z(i))
  225    continue
         d = dot(r,r,neq)  
         if (dsqrt(dabs(d)) .lt. conv) goto 300
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
c      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
c     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
      b(1:neq) = b(1:neq)*m(1:neq)
      energy   = dot(x(1),b(1),neq)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100) neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0)write(10,'(a,a,i9,a,d20.10,a,d20.10,f20.2)')
     .              "BICGSTAB: ", "it",j, " energy norm ",energy,
     .              " tol ",tol,"time",time
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA BICGSTAB:',/,5x,'Coeficiente da diagonal
     . nulo ou negativo - equacao ',i7)
 1100 format(' (BICGSTAB) solver:'/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.6/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i4,
     .        ' iterations !',/)
      end
c **********************************************************************
      subroutine pbicgstab(neq     ,nequ  ,nad,ia ,ja,
     .                     ad      ,au    ,al ,m  ,b ,  x,  
     .                     t       ,v     ,r  ,p ,  z,
     .                     tol     ,maxit ,
     .                     matvec  ,dot   ,
     .                     my_id   ,neqf1i,neqf2i,
     .                     neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,flog)
c **********************************************************************
c *                                                                    *
c *   Subroutine PBICGSTAB                                             *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   biconjugados com precondicionador diagonal para matrizes         *
c *   nao-simetricas.                                                  *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   nequ   - numero de equacoes no bloco Kuu                         *
c *   nad    - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   t(neq) - arranjo local de trabalho                               *
c *   v(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   p(neq) - arranjo local de trabalho                               *
c *   z(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   flog   - log do arquivo de saida                                 *
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
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................   
      integer neq,nequ,nad
      integer maxit,i,j,jj,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),t(*),v(*),z(*)
      real*8  dot,tol,conv,energy,d,alpha,beta,rr0,w,vi
      real*8  time0,time
      real*8  dum1 
      logical flog
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
c     if(my_id.eq.0) print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c 
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c ----------------------------------------------------------------------
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         p(i) = r(i)
         b(i) = p(i)
         z(i) = p(i)*m(i) 
  100 continue
      d    = dot(r,z,neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      jj = 1
      do 230 j = 1, maxit
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .               z,v,
     .               neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
         rr0   = dot(b,r,neq_doti)
         alpha = rr0/dot(v,r,neq_doti)
         do 210 i = 1, neq
            x(i) = x(i) + alpha * z(i)
            b(i) = b(i) - alpha * v(i)
            z(i) = b(i) * m(i)
  210    continue
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .               z,t,
     .               neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
         w = dot(t,b,neq_doti) / dot(t,t,neq_doti)
         do 220 i = 1, neq
            x(i) = x(i) + w*z(i)
            b(i) = b(i) - w*t(i)
  220    continue
         d = dot(b,z,neq_doti)
c        if ( j .eq. 300) goto 300 
         if (dsqrt(dabs(d)) .lt. conv) goto 300
         beta = (dot(r,b,neq_doti) / rr0)*(alpha/w)
         do 225 i = 1, neq
             p(i) = b(i) + beta*(p(i)-w*v(i))
             z(i) = p(i)*m(i)
  225    continue
c ......................................................................
         if( jj .eq. 200) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      energy   = dot(x,z,neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100)tol,neq,nad,j,energy,time
c ......................................................................
c     Controle de flops
      if(flog) then
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10,a,f20.2)')
     .               "PBICGSTAB: ","it",j, " energy norm ",energy,
     .               " tol ",tol," time ",time
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA BICGSTAB:',/,5x,'Coeficiente da diagonal
     . nulo ou negativo - equacao ',i7)
 1100 format(' (PBICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BICGSTAB:',5x,'It',i7,5x,2d20.10)
      end
c *********************************************************************
c
c *********************************************************************      
      subroutine pcg_block_it(neq   ,nequ  ,neqp  ,nad  ,naduu,nadpp
     .                       ,iau   ,jau   ,iap   ,jap  ,iapu ,japu    
     .                       ,adu   ,adp   ,alu   ,alp  ,alpu
     .                       ,mu    ,mp    ,b     ,x     
     .                       ,z     ,r     
     .                       ,bu    ,bp    ,bu0   ,bp0
     .                       ,u     ,p
     .                       ,tol   ,ctol  ,maxit ,cmaxit,alfap ,alfau 
     .                       ,fnew  ,istep
     .                       ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                       ,i_xfi ,i_rcvsi,i_dspli)
c **********************************************************************
c *                                                                    *
c *   Subroutine PCG_BLOCK_IT                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   nequ   - numero de equacoes no bloco Kuu                         *
c *   neqp   - numero de equacoes no bloco Kuu                         *
c *   nad    - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c *   naduu  - numero de termos nao nulos no bloco Kuu                 *  
c *   nadpp  - numero de termos nao nulos no bloco Kpp                 *
c *   ia(*)  - ponteiro do formato CSR (bloco Kuu e Kpp)               *
c *   ja(*)  - ponteiro das colunas no formato CSR (bloco Kuu e Kpp)   *
c *   iapu(*)- ponteiro do formato CSR (bloco Kpu )                    *
c *   japu(*)- ponteiro das colunas no formato CSR (bloco Kpu)         *
c *   ad(neq)- diagonal da matriz A  (bloco Kuu e Kpp)                 *
c *   al(*)  - parte triangular inferior de A (bloco Kuu e Kpp)        *
c *   alpu(*)- parte triangular inferior de A (bloco Kpu)              *
c *   m(*)   - precondicionador diagonal (bloco Kuu e Kpp)             *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - valores do passo anterior                               *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   bu(nequ) - arranjo local de trabalho                             *
c *   bp(neqp) - arranjo local de trabalho                             *
c *   bu0(nequ)- arranjo local de trabalho                             *
c *   bp0(neqp)- arranjo local de trabalho                             *
c *   u(nequ) - arranjo local de trabalho                              *
c *   p(neq0) - arranjo local de trabalho                              *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   ctol   - tolerancia de convergencia do ciclo externo             *
c *   cmaxit - numero maximo de iteracoes do ciclo externo             *
c *   fnew   - .true. chute inicial nulo                               * 
c *            .false. passo anterior                                  *
c *   istep  - numero do passo de tempo                                *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   energy - nao definido                                            *
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
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,neqp,nad,naduu,nadpp,maxit,i,j,jj,istep
      integer cmaxit
      integer iau(*),jau(*),iap(*),jap(*),iapu(*),japu(*)
      integer my_id,idum
      real*8 adu(*),adp(*),alu(*),alp(*),alpu(*)
      real*8 mu(*),mp(*),x(*),r(*),z(*),b(*)
      real*8 bp(*),bu(*),bp0(*),bu0(*)
      real*8 u(*),p(*)
      real*8 dot,tol,ctol,energy,d,alpha,beta
      real*8 resid_u,resid_p,p_conv,u_conv,alfap,alfau
      real*8 time0,time,time_csr
      real*8 dum1
      logical l_u_conv,l_p_conv,fnew
      external matvec_csrc_sym_pm
      external dot_par,dot
c ======================================================================
      time0    = MPI_Wtime()
      time_csr = 0.d0
c ......................................................................
c
c ...
      time0 = MPI_Wtime()
c ... 
c     alfap = 0.1d0
c     alfau = 0.2d0
c.......................................................................
c 
c ...
      if(fnew) then
        do j = 1, nequ 
          u(j) = 0.d0 
        enddo
        do j = 1, neqp
          p(j) = 0.d0        
        enddo
c.......................................................................
c 
c ...
      else
        do j = 1, nequ 
          u(j) = x(j) 
        enddo
        do j = 1, neqp
          p(j) = x(nequ+j) 
        enddo
      endif
c.......................................................................
c 
c ... 
      l_u_conv = .false.
      l_p_conv = .false.
c.......................................................................
c 
c ... bu 
      call aequalb(bu0,b,nequ)
c ... bp
      call aequalb(bp0,b(nequ+1),neqp)
c.......................................................................
c
c ... 
c     print*,ctol,cmaxit,maxit,tol
      do i = 1, cmaxit
c
c ...  r = Fp-kpu*U
        time_csr = MPI_Wtime() - time_csr
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,u,r,.false.)
        call aminusb(bp0,r,bp,neqp) 
        time_csr = MPI_Wtime() - time_csr
c ......................................................................
c
c ... P = inv(Kpp)*(Fp - kpu*U)
        call pcg(neqp      ,neqp       ,nadpp
     .          ,iap       ,jap
     .          ,adp       ,alp        ,alp            
     .          ,mp        ,bp         ,x               
     .          ,z         ,r          ,tol,maxit
     .          ,matvec_csrc_sym_pm,dot_par 
     .          ,my_id ,neqf1i ,neqf2i,neqp    ,i_fmapi
     .          ,i_xfi ,i_rcvsi,i_dspli,.false.)
c ......................................................................
c
c ... x - > p
       do j = 1, neqp
         p(j) = (1.d0-alfap)*p(j) + alfap*x(j)
       enddo
c ......................................................................
c
c ...  r = Fu-kup*P
        time_csr = MPI_Wtime() - time_csr
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,p,r,.true.)
        call aminusb(bu0,r,bu,nequ)
        time_csr = MPI_Wtime() - time_csr
c ......................................................................
c
c ... U = inv(Kuu)*(Fu - kup*P)
        call pcg(nequ   ,nequ   ,naduu
     .          ,iau    ,jau
     .          ,adu    ,alu    ,alu
     .          ,mu     ,bu     ,x
     .          ,z      ,r      
     .          ,tol    ,maxit
     .          ,matvec_csrc_sym_pm,dot_par 
     .          ,my_id ,neqf1i ,neqf2i,nequ    ,i_fmapi
     .          ,i_xfi ,i_rcvsi,i_dspli,.true.)
c ......................................................................
c
c ... x - > u
       do j = 1, nequ 
         u(j) = (1.d0-alfau)*u(j) + alfau*x(j) 
       enddo
c ......................................................................
c            
c ... Kpp*P
        time_csr = MPI_Wtime() - time_csr
        call matvec_csrc_sym_pm(neqp      ,neqp
     .                        ,iap       ,jap ,idum     ,idum        
     .                        ,adp       ,alp ,alp       
     .                        ,p         ,z
     .                        ,neqf1i,neqf2i
     .                        ,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ... Kpu*U
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,u,r,.false.)
c ...
        call aminusb(bp0,r,bp,neqp)
        call aminusb(bp ,z,bp,neqp)
        time_csr = MPI_Wtime() - time_csr
        resid_p = dsqrt(dot(bp,bp,neqp))
c .....................................................................
c            
c ... Kuu*U
        time_csr = MPI_Wtime() - time_csr
        call matvec_csrc_sym_pm(nequ      ,nequ
     .                        ,iau       ,jau ,idum     ,idum         
     .                        ,adu       ,alu ,alu              
     .                        ,u         ,z
     .                        ,neqf1i,neqf2i
     .                        ,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ... Kpu*P
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,p,r,.true.)
c ...
        call aminusb(bu0,r,bu,nequ)
        call aminusb(bu,z,bu,nequ)
        time_csr = MPI_Wtime() - time_csr
        resid_u = dsqrt(dot(bu,bu,nequ))
c ... 
        if( i .eq. 1) then
          u_conv = dsqrt(dot(b,b,neq))*ctol
          p_conv = dsqrt(dot(b,b,neq))*ctol
c ......................................................................
c 
c ... 
        else 
          if( resid_u .lt. u_conv ) l_u_conv = .true.
          if( resid_p .lt. p_conv ) l_p_conv = .true.
        endif
c ......................................................................
        print*,'it',i,'t',resid_p ,p_conv,resid_u,u_conv
        if( l_u_conv .and. l_p_conv ) goto 100
        if( i .eq. 500) then
          print*, 'MAXIT'
          stop
        endif
      enddo
  100 continue   
c ......................................................................
      time = MPI_Wtime()
c ......................................................................
      write(16,*) istep,i,time-time0,time_csr
c ...
      do j = 1, nequ 
        x(j) = u(j) 
      enddo
      do j = 1, neqp
        x(nequ+j) = p(j) 
      enddo
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      real*8 function smachn()
c **********************************************************************
c *                                                                    *
c *   SMACHN: calcula a precisao da maquina para real*8                *
c *                                                                    *
c **********************************************************************
      implicit none
      smachn = 1.0d0
100   smachn = smachn*0.5d0
      if(smachn+1.d0 .ne. 1.d0) go to 100
c     smachn = 2.0d0*smachn
      return
      end
c ***********************************************************************
