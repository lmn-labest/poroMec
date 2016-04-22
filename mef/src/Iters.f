c *********************************************************************
c * Metodos iterativos para solucao de sistemas lineares              *
c * ----------------------------------------------------------------- *
c * simetricos:                                                       *
c * ----------------------------------------------------------------- *
c * CG   - gradiente conjugados                                       *
c *                                                                   *
c * PCG  - gradiente conjugados com precondicionador diagonal         *
c *                                                                   *
c * ICCG - gradiente conjugados com precondicionador de fatoracoes    *              
c * incompletas                                                       *
c *                                                                   *
c * ----------------------------------------------------------------- *
c * nao-simetricos:                                                   *
c * ----------------------------------------------------------------- *
c * bicgstab - gradiente bi-conjugados estabilizados                  *
c *                                                                   *
c * pbicgstab - gradiente bi-conjugados estabilizados  com            * 
c * precondicionador diagonal                                         *
c *                                                                   *
c * icbicgstab - gradiente bi-conjugados estabilizados fatoracoes     *              
c * incompletas                                                       *   
c *                                                                   *
c * gmres(m) - GMRES com precondicionador diagonal                    *
c *                                                                   *
c * pcg_split - resolucao iterativa do problema poro mecanico com     *
c * matriz blocada | Kuu Kup | onde kup = -kpu                        * 
c *                | kpu Kpp |                                        *
c * ----------------------------------------------------------------- *
c *********************************************************************  
      subroutine cg(neq   ,nequ  ,nad   ,ia      ,ja
     .             ,ad    ,au    ,al    ,b       ,x
     .             ,z     ,r     ,p     ,tol     ,maxit
     .             ,matvec,dot
     .             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .             ,i_xfi ,i_rcvsi,i_dspli
     .             ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 20/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * Subroutine CG : Solucao de sistemas de equacoes pelo metodo dos    *    
c * gradientes conjugados                                              *
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
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
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
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ad(*),al(*),au(*) e b - inalterados                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
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
      real*8  ad(*),au(*),al(*),x(*),b(*)
      real*8  r(*),z(*),p(*)
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp
      real*8  time0,time
      real*8 dum1
      logical flog,fprint,fnew
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
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif  
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)  
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
c ... p0 = r0
         p(i) = r(i)
  100 continue
c ... ( r(0),r(0) )
      d    = dot(r,r,neq_doti)
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),r(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
  210    continue
c .....................................................................
c
c ...    
         di   = dot(r,r,neq_doti) 
c ... beta = ( r(j+1),r(j+1) ) / ( r(j),r(j) )
         beta = di / d
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... p(j+1) = r(j+1) + beta*p(j)
            p(i) = r(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         d =  di
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
c
c ...
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
       xkx = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. conv ) then
        write(*,1400) tmp
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "CG: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA CG:',/,1x,'Coeficiente da diagonal nulo '
     .,i9)
 1100 format(' (CG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' CG:',5x,'It',i7,5x,2d20.10)
 1400 format (' CG:',1x,'Residuo exato > conv ',1x,d20.10)
      end
c *********************************************************************  
c
c *********************************************************************  
      subroutine pcg(neq   ,nequ   ,nad   ,ia       ,ja
     .              ,ad    ,au     ,al    ,m        ,b      
     .              ,x     ,z      ,r     ,p     
     .              ,tol   ,maxit
     .              ,matvec,dot
     .              ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .              ,i_xfi ,i_rcvsi,i_dspli
     .              ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 20/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * Subroutine PCG : Solucao de sistemas de equacoes pelo metodo dos   *
c * gradientes conjugados com precondicionador diagonal para matrizes  *
c * simetricas.                                                        *
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
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
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
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
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
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,nad,maxit,i,j,jj
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8  r(*),z(*),p(*)
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp
      real*8  time0,time
      real*8 dum1
      logical flog,fprint,fnew
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
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif  
c .......................................................................
c
c ... conv = tol * |(M-1)b|
      do 15 i = 1, neq
         z(i) = b(i) * m(i)
   15 continue
      d    = dot(z,z,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0                                                            
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
c ... z0 = (M-1)r0
         z(i) = r(i) * m(i)
c ... p0 = r0
         p(i) = z(i)
  100 continue
c ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )
      d    = dot(r,z,neq_doti)
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
c ... z  = (M-1)r0
            z(i) = r(i) * m(i)
  210    continue
c .....................................................................
c
c ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
         di   = dot(r,z,neq_doti) 
c ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
         beta = di / d
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... p(j+1) = (M-1)r(j+1) + beta*p(j) = z + beta*p(j)
            p(i) = z(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         d =  di
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
c ... Energy norm:  x*Kx
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .           ,x,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'PCG: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Coeficiente da diagonal ',i9)
 1100 format(' (PCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PCG:',5x,'It',i7,5x,2d20.10)
 1400 format (' PCG:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c *********************************************************************  
c
c *********************************************************************  
      subroutine gmres(neq,nequ,nad,ia,ja,ad,au,al,m,b,x,k,g,h,y,c,s,e,
     .              tol,maxit,matvec,dot,neqovlp,my_id,neqf1i,neqf2i,
     .              neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,flog)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 22/04/2016                                    * 
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
c * i_dspli  -                                                         *
c * flog     - log do arquivo de saida                                 *
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
c * g(neq+1,k+1)                                                       *
c * h(k+1,k)                                                           *
c * y(k)                                                               *
c * c(k)                                                               *
c * s(k)                                                               *
c * e(k+1)                                                             *
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
      real*8  xkx,econv,norm,dot,r,aux1,aux2,beta
      real*8  time0,time
      real*8 dum1
      logical flog
      external matvec,dot
      integer my_id
c ......................................................................
      time0 = MPI_Wtime()
c ......................................................................
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
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,x,g(1,1)
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi 
     .              ,i_rcvsi,i_dspli,dum1)
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
            call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al
     .                 ,al(nad+1)
     .                 ,g(1,i),g(1,i+1)
     .                 ,neqf1i,neqf2i 
     .                 ,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
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
c         nii(l)=ni
         if (dabs(e(ni+1)) .le. econv) goto 1100
c ......................................................................
 1000 continue
c ......................................................................
 1100 continue
c
c ... Norma da solucao: x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,
     .           al(nad+1),x     ,g(1,1)  ,neqf1i,neqf2i,
     .           i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx = dot(x,g,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 1200 i = 1, neq
        g(i,2) = b(i) - g(i,1)
 1200 continue
      aux1 = dot(g(1,2),g(1,2),neq_doti)
      aux1 = dsqrt(aux1)
      if( aux1 .gt. 3.16d0*econv ) then
         if(my_id .eq.0 )then
           write(*,2400) aux1,econv
         endif 
      endif
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
      if(my_id.eq.0)write(*,2000) tol,neq,l,nit,dabs(e(ni+1)),xkx,norm
     .                           ,time
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,i9,a,d20.10,a,f20.2)')
     .         'GMRES: ',' it ',nit, ' x * Kx ',xkx,' ||x|| ',norm,
     .         ' nKylov ',k,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.10/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for '
     .      ,i9,' cycles !',5x,i7,' nKylov',5x,' It ',i7/)
 2300 format (' GMRES:',5x,'cycles',i7,5x,'It',i7,5x,2d20.10)
 2400 format (' GMRES:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c **********************************************************************
c
c **********************************************************************
      subroutine bicgstab(neq      ,nequ  ,nad,ia ,ja 
     .                    ,ad      ,au    ,al,b  ,x   
     .                    ,t       ,v     ,r ,p  ,r0
     .                    ,tol     ,maxit  
     .                    ,matvec  ,dot    
     .                    ,my_id   ,neqf1i,neqf2i 
     .                    ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                    ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 15/04/2016                                    *
c * Data de modificaco : 20/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * BICGSTAB  : Solucao de sistemas de equacoes pelo metodo dos        * 
c * gradientes biconjugados para matrizes nao-simetricas.              *                                         *
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
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) e b - inalterados                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *                                                                      *
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
      real*8  ad(*),au(*),al(*),x(*),b(*)
      real*8  r(*),p(*),t(*),v(*),r0(*)
      real*8  dot,tol,conv,xkx,norm,d,alpha,beta,rr0,w,tmp
      real*8  time0,time
      real*8  dum1 
      logical flog,fprint,fnew
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
c     if(my_id.eq.0) print *, 'nad :',nad
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
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif  
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ...
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,p,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - p(i)
c ... r = r0
         r(i)  = r0(i)
c ... p = r0
         p(i)  = r0(i)
  100 continue
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit
c ... v = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,v
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... alpha = ( r(j),r0 ) / ( Ap(j), r0 ))
         rr0   = dot(r,r0,neq_doti)
         alpha = rr0/dot(v,r0,neq_doti)
c .......................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... s(j)   = r(j) - alpha*Ap(j)
            r(i) = r(i) - alpha * v(i)
  210    continue
c ........................................................................
c
c ... (s,s)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... t = Ar(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,r,t 
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ........................................................................
c
c ... w = ( Ar(j),r(j) ) / ( Ar(j), Ar(j) ))
         w = dot(t,r,neq_doti) / dot(t,t,neq_doti)
c ........................................................................
c
c ... 
         do 220 i = 1, neq
c ... x(j+1) = x(j) + w*r(j)
            x(i) = x(i) + w*r(i)
c ... r(j+1) = s(j) - w*As(j)
            r(i) = r(i) - w*t(i)
  220    continue
c .......................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c .......................................................................
c
c ... beta = ( r(j+1),r0 ) / ( r(j), r0 )) * (alpha/w) 
         beta = (dot(r,b,neq_doti) / rr0)*(alpha/w)
c .......................................................................
c
c ...
         do 225 i = 1, neq
             p(i) = r(i) + beta*(p(i)-w*v(i))
  225    continue
c .......................................................................
c
c ...
         if( jj .eq. 500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
c
c ...
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c ......................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,p 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,p,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - p(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. conv ) then
        write(*,1400) tmp
      endif
c ......................................................................
c
c ...
      time = MPI_Wtime()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "BICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA BICCSTAB:',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (BICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BICGSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' BICGSTAB:',1x,'Residuo exato > conv ',1x,d20.10)
      end
c *********************************************************************
c
c **********************************************************************
      subroutine pbicgstab(neq     ,nequ  ,nad,ia ,ja 
     .                    ,ad      ,au    ,al ,m  ,b ,x   
     .                    ,t       ,v     ,r  ,p  ,z ,r0 
     .                    ,tol     ,maxit  
     .                    ,matvec  ,dot    
     .                    ,my_id   ,neqf1i,neqf2i 
     .                    ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                    ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 20/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTAB : Solucao de sistemas de equacoes pelo metodo dos        * 
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
c * r0(neq)  - arranjo local de trabalho                               *
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
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ad(*),al(*),au(*) e b(*)  inalterados                              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * A(M-1)y=b precondicionador direita                                 *
c * ------------------------------------------------------------------ *                                                                      *
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
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*)
      real*8  r(*),p(*),t(*),v(*),z(*),r0(*)
      real*8  dot,tol,conv,d,alpha,beta,rr0,w,xkx,norm,tmp
      real*8  time0,time
      real*8  dum1 
      logical flog,fprint,fnew
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
c ... Chute inicial:
c 
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - z(i)
c ... r = r0
         p(i)  = r0(i)
c ... p = r0
         r(i)  = r0(i)
c ... z = M(-1)p
         z(i)  = p(i)*m(i) 
  100 continue
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit
c ... v = Az(j) = AM(-1)p(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,z,v 
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... alpha = ( r(j),r0 ) / ( AM(-1)p(j), r0 ))
         rr0   = dot(r,r0,neq_doti)
         alpha = rr0/dot(v,r0,neq_doti)
c .......................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*M(-1)p
            x(i) = x(i) + alpha * z(i)
c ... s(j)   = r(j) - alpha*AM(-1)p(j)
            r(i) = r(i) - alpha * v(i)
c ... z = M(-1)s
            z(i) = r(i) * m(i)
  210    continue
c ........................................................................
c
c ... (s,s)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... t = Az = AM(-1)s(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,z,t 
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ........................................................................
c
c ... w = ( AM(-1)s(j) ,s(j) ) / ( AM(-1)s(j), AM(-1)s(j) )
         w = dot(t,r,neq_doti) / dot(t,t,neq_doti)
c ........................................................................
c
c ... 
         do 220 i = 1, neq
c ... x(j+1) = x(j) + w*M(-1)s
            x(i) = x(i) + w*z(i)
c ... r(j+1) = s(j) - w*AM(-1)s
            r(i) = r(i) - w*t(i)
  220    continue
c ........................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... beta = ( r(j+1),r0 ) / ( r(j), r0 )) * (alpha/w) 
         beta = (dot(r,r0,neq_doti) / rr0)*(alpha/w)
c .......................................................................
c
c ...
         do 225 i = 1, neq
c ... p(j+1) = r(i) + beta*(p(j)-w*v(i))
             p(i) = r(i) + beta*(p(i)-w*v(i))
c ... z = M(-1)p
             z(i) = p(i)*m(i)
  225    continue
c .......................................................................
c
c ...
         if( jj .eq. 500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
c
c ...
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c ......................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,z,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
        write(*,1400) tmp,conv
      endif
c ......................................................................
c
c ...
      time = MPI_Wtime()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "BICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA PBICCSTAB:',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (PBICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PBICGSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' PBICGSTAB:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c *********************************************************************
c
c **********************************************************************
      subroutine icbicgstab(neq     ,nequ  ,nad,ia ,ja 
     .                     ,ad      ,au    ,al ,m  ,b     ,x   
     .                     ,t       ,v     ,r  ,p  ,z     ,r0
     .                     ,tol     ,maxit  
     .                     ,matvec  ,dot    
     .                     ,my_id   ,neqf1i,neqf2i 
     .                     ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                     ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 15/04/2016                                    *
c * Data de modificaco : 22/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTAB : Solucao de sistemas de equacoes pelo metodo dos        * 
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
c * r0(neq)  - arranjo local de trabalho                               *
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
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * A(M-1)y=b precondicionador direita                                 *
c * ------------------------------------------------------------------ *                                                                      *
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
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*)
      real*8  r(*),p(*),t(*),v(*),z(*),r0(*)
      real*8  dot,tol,conv,xkx,norm,d,alpha,beta,rr0,w,tmp
      real*8  time0,time
      real*8  dum1 
      logical flog,fprint,fnew
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
c ... Chute inicial:
c 
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - z(i)
c ... r = r0
         p(i) = r0(i)
c ... p = r0
         r(i) = r0(i)
  100 continue
c .......................................................................
c
c ... Mz=p  
      call ildlt_solv(neq,ia,ja,m,m(neq+1),p,z)
c .......................................................................      
c
c ...
      jj = 1
      do 230 j = 1, maxit
c ... v = Az(j) = AM(-1)p(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .               z,v,
     .               neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... alpha = ( r(j),r0 ) / ( AM(-1)p(j), r0 ))
         rr0   = dot(r,r0,neq_doti)
         alpha = rr0/dot(v,r0,neq_doti)
c .......................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*M(-1)p
            x(i) = x(i) + alpha * z(i)
c ... s(j)   = r(j) - alpha*AM(-1)p(j)
            r(i) = r(i) - alpha * v(i)
  210    continue
c ........................................................................
c
c
c ... (s,s)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... Mz=r  
         call ildlt_solv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ...  t = Az = AM(-1)s(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .               z,t,
     .               neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ........................................................................
c
c ... w = ( AM(-1)s(j) ,s(j) ) / ( AM(-1)s(j), AM(-1)s(j) )
         w = dot(t,r,neq_doti) / dot(t,t,neq_doti)
c ........................................................................
c
c ... 
         do 220 i = 1, neq
c ... x(j+1) = x(j) + w*M(-1)s
            x(i) = x(i) + w*z(i)
c ... r(j+1) = s(j) - w*AM(-1)s
            r(i) = r(i) - w*t(i)
  220    continue
c ........................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... beta = ( r(j+1),r0 ) / ( r(j), r0 )) * (alpha/w) 
         beta = (dot(r,r0,neq_doti) / rr0)*(alpha/w)
c .......................................................................
c
c ...
         do 225 i = 1, neq
c ... p(j+1) = r(i) + beta*(p(j)-w*v(i))
             p(i) = r(i) + beta*(p(i)-w*v(i))
  225    continue
c .......................................................................
c
c ... Mz=p  
         call ildlt_solv(neq,ia,ja,m,m(neq+1),p,z)
c .......................................................................
c
c ...
         if( jj .eq. 500) then
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
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
       xkx = dot(x,z,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. conv ) then
        write(*,1400) tmp
      endif
c ......................................................................
c
c ...
      time = MPI_Wtime()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "ICBICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif

      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA ICBICCSTAB:',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (ICBICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' ICBICCSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' ICBICCSTAB:',1x,'Residuo exato > conv ',1x,d20.10)
      end
c *********************************************************************
c
c *********************************************************************      
      subroutine pcg_block_it(neq   ,nequ  ,neqp  
     .                       ,nad   ,naduu ,nadpp
     .                       ,iau   ,jau   
     .                       ,iap   ,jap  
     .                       ,iapu  ,japu    
     .                       ,adu   ,adp  
     .                       ,alu   ,alp  ,alpu
     .                       ,mu    ,mp    ,b     ,x     
     .                       ,z     ,r     ,s
     .                       ,bu    ,bp    ,bu0   ,bp0
     .                       ,u     ,p
     .                       ,tol   ,ctol  ,maxit ,cmaxit,alfap ,alfau 
     .                       ,fnew  ,istep
     .                       ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                       ,i_xfi ,i_rcvsi,i_dspli)
c **********************************************************************
c * Data de criacao    : 12/12/2015                                    *
c * Data de modificaco : 22/04/2016                                    * 
c * ------------------------------------------------------------------ *                                                                   *
c * Subroutine PCG_BLOCK_IT                                            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                    *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * neqp     - numero de equacoes no bloco Kpp                         *
c * nad      - numero de termos nao nulos no bloco Kuu, kpp e kpu      *
c * naduu    - numero de termos nao nulos no bloco Kuu                 *  
c * nadpp    - numero de termos nao nulos no bloco Kpp                 *
c * iau(*)   - ponteiro do formato CSR (bloco Kuu)                     *
c * jau(*)   - ponteiro das colunas no formato CSR (bloco Kuu)         *
c * iap(*)   - ponteiro do formato CSR (bloco Kp)                      *
c * jap(*)   - ponteiro das colunas no formato CSR (bloco Kpp)         *
c * iapu(*)  - ponteiro do formato CSR (bloco Kpu )                    *
c * japu(*)  - ponteiro das colunas no formato CSR (bloco Kpu)         *
c * adu(*)   - diagonal da matriz A  (bloco Kuu)                       *
c * alu(*)   - parte triangular inferior de A (bloco Kuu)              *
c * adp(*)   - diagonal da matriz A  (bloco Kpp)                       *
c * alp(*)   - parte triangular inferior de A (bloco Kpp)              *
c * alpu(*)  - parte triangular inferior de A (bloco Kpu)              *
c * mu(*)    - precondicionador diagonal (bloco Kuu)                   *
c * mp(*)    - precondicionador diagonal (bloco Kuu)                   *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - valores do passo anterior                               *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * bu(nequ) - arranjo local de trabalho                               *
c * bp(neqp) - arranjo local de trabalho                               *
c * bu0(nequ)- arranjo local de trabalho                               *
c * bp0(neqp)- arranjo local de trabalho                               *
c * u(nequ)  - arranjo local de trabalho                               *
c * p(neqp)  - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * ctol     - tolerancia de convergencia do ciclo externo             *
c * cmaxit   - numero maximo de iteracoes do ciclo externo             *
c * fnew     - .true. chute inicial nulo                               * 
c *          .false. passo anterior                                    *
c * istep    - numero do passo de tempo                                *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * energy   - nao definido                                            *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *  
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
      real*8 mu(*),mp(*),x(*),r(*),z(*),b(*),s(*)
      real*8 bp(*),bu(*),bp0(*),bu0(*)
      real*8 u(*),p(*),xkx,norm
      real*8 dot,tol,ctol,d,alpha,beta
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
      jj = 0
      do i = 1, cmaxit
c
c ...  rp= Fp-kpu*U
        time_csr = MPI_Wtime() - time_csr
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,u,r,.false.)
        call aminusb(bp0,r,bp,neqp) 
        time_csr = MPI_Wtime() - time_csr
c ......................................................................
c
c ... P = inv(Kpp)*(Fp - kpu*U)= inv(Kpp)*rp
        call pcg(neqp      ,neqp       ,nadpp,iap       ,jap
     .          ,adp       ,alp        ,alp  ,mp        ,bp  
     .          ,x(nequ+1) ,z          ,r    ,s     
     .          ,tol       ,maxit
     .          ,matvec_csrc_sym_pm,dot_par 
     .          ,my_id     ,neqf1i     ,neqf2i,neqp    ,i_fmapi
     .          ,i_xfi     ,i_rcvsi    ,i_dspli
     .          ,.false.   ,.false.    ,.false.)
c ......................................................................
c
c ... x - > p
       do j = 1, neqp
         p(j) = (1.d0-alfap)*p(j) + alfap*x(nequ+j)
       enddo
c ......................................................................
c
c ...  ru= Fu-kup*P
        time_csr = MPI_Wtime() - time_csr
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,p,r,.true.)
        call aminusb(bu0,r,bu,nequ)
        time_csr = MPI_Wtime() - time_csr
c ......................................................................
c
c ... U = inv(Kuu)*(Fu - kup*P)=inv(Kuu)*ru
        call pcg(nequ   ,nequ   ,naduu,iau    ,jau
     .          ,adu    ,alu    ,alu  ,mu     ,bu   
     .          ,x      ,z      ,r    ,s      
     .          ,tol    ,maxit
     .          ,matvec_csrc_sym_pm,dot_par 
     .          ,my_id  ,neqf1i ,neqf2i,nequ    ,i_fmapi
     .          ,i_xfi  ,i_rcvsi,i_dspli
     .          ,.false.,.true.  ,.false.)
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
          u_conv = dsqrt(dot(bu,bu,nequ))*ctol
          p_conv = dsqrt(dot(bp,bp,neqp))*ctol
c ......................................................................
c 
c ... 
        else 
          if( resid_u .lt. u_conv ) l_u_conv = .true.
          if( resid_p .lt. p_conv ) l_p_conv = .true.
        endif
c ......................................................................
        jj = jj + 1
        if( jj . eq. 10 ) then
          write(*,200),i,resid_p ,p_conv,resid_u,u_conv
          jj = 0
        endif
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
c ... calculo de x*(kx)
c
c ... x*kx = u*bu + p*bp
      xkx = dot(u,bu0,nequ) + dot(p,bp0,neqp)
c ......................................................................
c
c ... norm-2 = || x || = || u || + || p ||
      norm = dsqrt(dot(u,u,nequ)) + dsqrt(dot(p,p,neqp))
c ......................................................................
c
c ...
      write(*,300)ctol,neq,i,xkx,norm,time-time0
c ......................................................................
c
c ...
      write(16,*) istep,i,time-time0,time_csr
c ......................................................................
c
c ...
      do j = 1, nequ 
        x(j) = u(j) 
      enddo
      do j = 1, neqp
        x(nequ+j) = p(j) 
      enddo
c ......................................................................
      return
 200  format (1x,'It',1x,i4,1x,'rp',1x,2d10.2,1x,'ru',1x,2d10.2)
 300  format(' (BLOCKPCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
      end
c **********************************************************************
c
c *********************************************************************
      subroutine iccg(neq   ,nequ  ,nad   ,ia      ,ja
     .               ,ad    ,au    ,al    ,m       ,b       
     .               ,x     ,z     ,r     ,p   
     .               ,tol   ,maxit
     .               ,matvec,dot
     .               ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .               ,i_xfi ,i_rcvsi,i_dspli
     .               ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 20/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * IcCG : Solucao de sistemas de equacoes pelo metodo dos gradientes   *
c * conjugados com precondicionador incompleto                         *
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
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
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
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
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
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,nad,maxit,i,j,jj
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),x(*),b(*),r(*),z(*),p(*)
      real*8  m(*)
      real*8  dot,ddot,tol,conv,xkx,norm,d,di,alpha,beta,tmp
      real*8  time0,time
      real*8 dum1
      logical flog,fnew,fprint
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
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif  
c ......................................................................
c
c ... conv = tol * |(M-1)b|
      call ildlt_solv(neq,ia,ja,m,m(neq+1),b,z)
      d    = dot(z,z,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)  
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................      
c
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
  100 continue
c .......................................................................      
c
c ... Mz=r  
      call ildlt_solv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................      
c
c ...
      do 105 i = 1, neq
c ... p0 = r0
         p(i) = z(i)
  105 continue
c .......................................................................      
c
c ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )
      d    = dot(r,z,neq_doti)
c .......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
  210    continue
c ......................................................................
c
c ... Mz=r  
         call ildlt_solv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................

c ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
         di   = dot(r,z,neq_doti) 
c ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
         beta = di / d
c .....................................................................
c
c ...         
         do 220 i = 1, neq
            p(i) = z(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         d = di           
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
c ... Energy norm: x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .          ,x,z 
     .          ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx    = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. conv ) then
        write(*,1400) tmp
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "ICCG: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Coeficiente da diagonal nulo
     .ou negativo - equacao ',i7)
 1100 format(' (ICG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format ( 'ICCG: ',5x,'It',i7,5x,2d20.10)
 1400 format (' ICCG:',1x,'Residuo exato > conv ',1x,d20.10)
      end
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
c     subroutine bicgstab(neq,ia,ja,ad,au,al,m,b,x,y,z,p,r,s,tol,maxit,
c    .                    matvec,dot,my_id,neqf1i,neqf2i,neq_doti,
c    .                    i_fmapi,i_xfi,i_rcvsi,i_dspli)
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
c     implicit none
c     include 'mpif.h'
c     integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
c     integer*8 i_fmapi,i_xfi
c     integer*8 i_rcvsi,i_dspli
c .....................................................................
c     integer neq,maxit,nad,i,j,k
c     integer ia(*),ja(*),my_id
c     real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),y(*),z(*),s(*)
c     real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
c     real*8  time0,time
c     real*8 dum1
c     external matvec,dot
c ======================================================================
c     time0 = MPI_Wtime() 
c ......................................................................
c     nad = ia(neq+1)-1
c     if(my_id.eq.0)print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
c     do 10 i = 1, neq
c        x(i) = 0.d0
c ...... pre-condicionador diagonal:         
c        b(i)  = b(i)/m(i)
c        ad(i) = ad(i)/m(i)
c        do 5 k = ia(i), ia(i+1)-1
c           j = ja(k)
c           al(k) = al(k) / m(i)
c           au(k) = au(k) / m(j)
c  5     continue      
c 10  continue
c ----------------------------------------------------------------------
c     call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,r,
c    .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)    
c     do 100 i = 1, neq
c        r(i) = b(i) - r(i)
c        p(i) = r(i)
c        b(i) = r(i) 
c 100 continue
c     d    = dot(r(1),r(1),neq)
c     conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
c     do 230 j = 1, maxit
c        call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
c    .               p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .               dum1)
c        rr0 = dot(r,b,neq)
c        alpha = rr0/dot(z,b,neq)
c         do 210 i = 1, neq
c           s(i) = r(i) - alpha * z(i)
c 210    continue
c        call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
c    .               s,y, neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .               dum1)
c        w = dot(y,s,neq) / dot(y,y,neq)
c        do 220 i = 1, neq
c           x(i) = x(i) + alpha*p(i) + w * s(i)
c           r(i) = s(i) - w*y(i)
c 220    continue
c        beta = (dot(r,b,neq) / rr0)*(alpha/w)
c        do 225 i = 1, neq
c            p(i) = r(i) + beta*(p(i)-w*z(i))
c 225    continue
c        d = dot(r,r,neq)  
c        if (dsqrt(dabs(d)) .lt. conv) goto 300
c 230 continue
c ----------------------------------------------------------------------
c     write(*,1200) maxit
c     call stop_mef()
c 300 continue
c
c ... Energy norm:
c
c      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
c     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c     b(1:neq) = b(1:neq)*m(1:neq)  
c     energy   = dot(x(1),b(1),neq)
c ......................................................................
c     time = MPI_Wtime()
c     time = time-time0
c ----------------------------------------------------------------------
c     if(my_id.eq.0)write(*,1100) neq,j,energy,time
c ......................................................................
c     Controle de flops
c     if(my_id.eq.0)write(10,'(a,a,i9,a,d20.10,a,d20.10,f20.2)')
c    .              "BICGSTAB: ", "it",j, " energy norm ",energy,
c    .              " tol ",tol,"time",time
c ......................................................................
c     return
c ======================================================================
c1000 format (//,5x,'SUBROTINA BICGSTAB:',/,5x,'Coeficiente da diagonal
c    . nulo ou negativo - equacao ',i7)
c1100 format(' (BICGSTAB) solver:'/
c    . 5x,'Number of equations  = ',i20/
c    . 5x,'Number of iterations = ',i20/
c    . 5x,'Energy norm          = ',d20.6/
c    . 5x,'CPU time (s)         = ',f20.2/)
c1200 format (' *** WARNING: No convergence reached after ',i4,
c    .        ' iterations !',/)
c     end