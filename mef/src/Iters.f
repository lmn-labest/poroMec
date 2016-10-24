c *********************************************************************
c * Metodos iterativos para solucao de sistemas lineares              *
c * ----------------------------------------------------------------- *
c * simetricos:                                                       *
c * ----------------------------------------------------------------- *
c * CG   - gradiente conjugados                                       *
c *                                                                   *
c * PCG  - gradiente conjugados com precondicionador diagonal         *
c *                                                                   *
c * BCG  - gradiente conjugados com precondicionador bloco diagonal   *
c *                                                                   *
c * ICCG - gradiente conjugados com precondicionador de fatoracoes    *   
c * incompletas  (LLT e LDLt)                                         *
c *                                                                   *
c * MINRES - matriz indefinida                                        *
c *                                                                   *
c * PMINRES - MINRES com precondicionador diagonal M=D(1/2)D(1/2)     *
c *                                                                   *
c * CR -Residuos conjugados                                           *
c *                                                                   *
c * PCR -Residuos conjugados com precondicionador diagonal            *
c *                                                                   *
c * SQRM - QRM simetrico                                              *
c *                                                                   *
c * RSQRM - QRM simetrico com precondicionador diagonal a direita     *
c *                                                                   *
c * LSQRM - QRM simetrico com precondicionador diagonal a esquerda    *
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
c * incompletas (LLT e LDLt)                                          *  
c *                                                                   *
c * bicgstabl2 - gradiente bi-conjugados estabilizados (l=2)          *
c *                                                                   *
c * pbicgstabl2- gradiente bi-conjugados estabilizados (l=2) com      *
c * precondicionador diagonal                                         *
c *                                                                   *
c * gmres(m) - GMRES com precondicionador diagonal                    *
c *                                                                   *
c * gmres2(m) - GMRES com precondicionador diagonal                   *
c *   (Ry = e) versao coluna                                          *
c *   (x = Vy) versao coluna                                          *
c *   versao com refined modified gram-schmidt                        * 
c *                                                                   *
c * block_it_pcg - resolucao iterativa do problema poro mecanico com  *
c * matriz blocada | Kuu Kup | onde kup = -kpu                        * 
c *                | kpu Kpp |                                        *
c * ----------------------------------------------------------------- *
c *********************************************************************  
      subroutine cg(neq   ,nequ  ,nad     ,ia      ,ja
     .             ,ad    ,au    ,al      ,b       ,x
     .             ,z     ,r     ,p       ,tol     ,maxit
     .             ,matvec,dot
     .             ,my_id ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     .             ,i_xfi ,i_rcvsi,i_dspli
     .             ,fprint,flog   ,fhist  ,fnew)
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
c * fhist    - log dos resuduos por iteracao                           *
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
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp,norm_b
      real*8  time0,time
      real*8 dum1
      logical flog,fprint,fnew,fhist
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
      d      = dot(b,b,neq_doti)
      norm_b = dsqrt(d)
      conv   = tol*dsqrt(d)
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
         if(fhist) write(18,1500),j,dsqrt(d)/norm_b 
c .....................................................................
c
c ...
         d =  di
         if (dsqrt(d) .lt. conv) goto 300
c ......................................................................
c
c ...
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(d),conv 
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
      if( tmp .gt. 3.16d0*conv ) then
        write(*,1400) tmp,conv
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,tmp,time
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
 1000 format (//,5x,'SUBROTINA CG:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (CG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' CG:',5x,'It',i7,5x,2d20.10)
 1400 format (' CG:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
 1500 format (' CG: ',5x,i7,5x,2es20.10)
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
     .              ,fprint,flog   ,fhist  ,fnew)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 18/06/2016                                    * 
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
c * fhist    - log dos resuduos por iteracao                           *
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
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp,norm_b
      real*8  norm_r,norm_m_r
      real*8  time0,time
      real*8 dum1
      logical flog,fprint,fnew,fhist
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
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
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif 
c .......................................................................
c
c ... conv = tol * |(M-1)b|m = tol *(b,M-1b)
      do 15 i = 1, neq
         z(i) = b(i) * m(i)
   15 continue
      d      = dot(b,z,neq_doti)
      norm_b = dsqrt(dabs(d))  
      conv   = tol*dsqrt(dabs(d))
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
         if(fhist) write(18,1500),j,dsqrt(dabs(d))/norm_b 
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
c ... r =M(-1)(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
        z(i) = r(i)*m(i)
  310 continue
      norm_m_r = dot(r,z,neq_doti)
      norm_m_r = dsqrt(dabs(norm_m_r))
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_m_r,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,norm_m_r,time
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
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Diagonal coefficient ' 
     . '- equation ',i9)
 1100 format(' (PCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,d20.10
     .        ' iterations !',/)
 1300 format (' PCG:',5x,'It',i7,5x,2d20.10)
 1400 format (' PCG:',1x,'Explicit residual > tol * ||b||| :'
     .       ,1x,d20.10,1x,d20.10)
 1500 format ( 5x,i7,5x,2es20.10)
      end
c *********************************************************************  
c
c *********************************************************************
      subroutine iccg(neq   ,nequ   ,nad     ,ia      ,ja
     .               ,ad    ,au     ,al      ,m       ,b       
     .               ,x     ,z      ,r       ,p   
     .               ,tol   ,maxit
     .               ,matvec,dot    ,triasolv
     .               ,my_id ,neqf1i ,neqf2i  ,neq_doti,i_fmapi
     .               ,i_xfi ,i_rcvsi,i_dspli
     .               ,fprint,flog   ,fhist   ,fnew)
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
c * triasolv - nome da funcao externa para o o solver triangular       *
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
c * fhist    - log dos resuduos por iteracao                           *
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
      real*8  norm_r,norm_m_r,norm_b
      real*8  time0,time
      real*8 dum1
      logical flog,fnew,fprint,fhist
      external matvec,dot,triasolv
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
c ... conv = tol * |(M-1)b|m = tol *(b,M-1b)
      call triasolv(neq,ia,ja,m,m(neq+1),b,z)
      d      = dot(b,z,neq_doti)
      norm_b = dsqrt(dabs(d))
      conv   = tol*dsqrt(dabs(d))
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
      call triasolv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................      
c
c ...
      do 105 i = 1, neq
c ... p0 = r0
         p(i) = z(i)
  105 continue
c .......................................................................      
c
c ... ( z(0),z(0) )m = ( r(0), z0 ) (r(0), (M-1)r0 )
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
         call triasolv(neq,ia,ja,m,m(neq+1),r,z)
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
         if(fhist) write(18,1500),j,dsqrt(dabs(d))/norm_b 
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
c ... r =(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
c ... Mz=r  
      call triasolv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................
      norm_m_r = dot(r,z,neq_doti)
      norm_m_r = dsqrt(norm_m_r)
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_m_r,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,norm_m_r,time
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
 1000 format (//,5x,'SUBROTINA ICCG:',/,5x,'Coeficiente da diagonal nulo
     . - equacao ',i9)
 1100 format(' (ICCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' ICCG: ',5x,'It',i7,5x,2d20.10)
 1400 format (' ICCG:',1x,'Residuo exato > 3.16*conv ',1x,d20.10
     .       ,1x,d20.10)
 1500 format ( 'ICCG ',5x,i7,5x,2es20.10)
      end
c **********************************************************************
c
c *********************************************************************
      subroutine bpcg(neq   ,nequ   ,nad     ,ia      ,ja
     .               ,ad    ,au     ,al      ,m       ,b       
     .               ,x     ,z      ,r       ,p   
     .               ,tol   ,maxit  
     .               ,matvec,dot    
     .               ,my_id ,neqf1i ,neqf2i  ,neq_doti,i_fmapi
     .               ,i_xfi ,i_rcvsi,i_dspli
     .               ,fprint,flog   ,fhist   ,fnew)
c **********************************************************************
c * Data de criacao    : 17/06/2016                                    *
c * Data de modificaco : 18/06/2016                                    * 
c * ------------------------------------------------------------------ *   
c * BPCG : Solucao de sistemas de equacoes pelo metodo dos gradientes  *
c * conjugados com precondicionador bloco diagonal                     *
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
c * m(*)     - precondicionador bloco diagonal                         *
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
c * fhist    - log dos resuduos por iteracao                           *
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
      include 'precond.fi'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,nad,maxit,i,j,jj
      integer ia(*),ja(*),my_id,bsize
      real*8  ad(*),au(*),al(*),x(*),b(*),r(*),z(*),p(*)
      real*8  m(*),max_block_a(max_block*max_block)
      real*8  dot,ddot,tol,conv,xkx,norm,d,di,alpha,beta,tmp
      real*8  norm_r,norm_m_r,norm_b
      real*8  time0,time
      real*8 dum1
      logical flog,fnew,fprint,fhist
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
c ... conv = tol * |(M-1)b|m = tol *(b,M-1b)
      call op_block_precond(m,b,z,max_block_a,iparam)
      d      = dot(b,z,neq_doti)
      norm_b = dsqrt(dabs(d))
      conv   = tol*dsqrt(dabs(d))
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
c ... z=M(-1)r  
      call op_block_precond(m,r,z,max_block_a,iparam)   
c .......................................................................      
c
c ...
      do 105 i = 1, neq
c ... p0 = r0
         p(i) = z(i)
  105 continue
c .......................................................................      
c
c ... ( z(0),z(0) )m = ( r(0), z0 ) (r(0), (M-1)r0 )
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
c ... z=M(-1)r 
         call op_block_precond(m,r,z,max_block_a,iparam)    
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
         if(fhist) write(18,1500),j,dsqrt(dabs(d))/norm_b 
c .....................................................................
c
c ...
         d = di
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
         if( jj .eq.500)then
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
c ... r =M(-1)(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
c ... z=M(-1)r
      call op_block_precond(m,r,z,max_block_a,iparam) 
c
      norm_m_r = dot(r,z,neq_doti)
      norm_m_r = dsqrt(dabs(norm_m_r))
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_m_r,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      bsize = iparam(3)*iparam(1)+ iparam(2) 
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,iparam(4),bsize
     .              ,xkx,norm,norm_r,norm_m_r,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,i3,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'BPCG: ',' it ',j,' block ',iparam(4), ' x * Kx '
     .      ,xkx,' ||x|| ',norm,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA BPCG:',/,5x,'Coeficiente da diagonal nulo
     . - equacao ',i9)
 1100 format(' (BPCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Size of blocks       = ',i20/
     . 5x,'Size of M            = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BPCG: ',5x,'It',i7,5x,2d20.10)
 1400 format (' BPCG:',1x,'Residuo exato > 3.16*conv ',1x,d20.10
     .       ,1x,d20.10)
 1500 format ( 'BPCG: ',5x,i7,5x,2es20.10)
      end
c **********************************************************************
c
c **********************************************************************  
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
c *********************************************************************  
      subroutine gmres2(neq,nequ,nad,ia,ja,ad,au,al,m,b,x,k,g,h,y,c,s,e,
     .              tol,maxit,matvec,dot,neqovlp,my_id,neqf1i,neqf2i,
     .              neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,flog)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 30/05/2016                                    * 
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
c *                                                                    *
c * versao com solver triagular superio (Ry = e) versao coluna         *
c * versao com matriz vetor geral ( x = Vy ) versao coluna             *
c * versao com refined modified gram-schmidt                           *
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
      real*8  xkx,econv,norm,dot,r,aux1,aux2,beta,inorm
      real*8 tau,kesp
      real*8  time0,time
      real*8 dum1
      logical flog
      external matvec,dot
      integer my_id
      parameter (kesp = 0.25d0)
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
         inorm = 1.d0/e(1) 
         do 210 i = 1, neq
            g(i,1) = g(i,1)*inorm
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
c ......... Ortogonalizacao (Gram-Schmidt modificado com refinamento):
c
c ...
            tau =  dsqrt(dot(g(1,i+1),g(1,i+1),neq_doti))
c .....................................................................
            do 320 j = 1, i
               beta = dot(g(1,i+1),g(1,j),neq_doti)
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
               h(j,i) = beta
  320       continue
c ......................................................................
c
c ....
            norm = dsqrt(dot(g(1,i+1),g(1,i+1),neq_doti))
            if( norm .le. kesp*tau) then
              do 321 j = 1, i
                  beta = dot(g(1,i+1),g(1,j),neq_doti)
                  do 311 ic = 1, neq
                    g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  311             continue
                  h(j,i) = h(j,i) + beta
  321          continue
            endif
c ......................................................................
c
c ......... Norma de g(i+1):
c
            norm = dsqrt(dot(g(1,i+1),g(1,i+1),neq_doti))
c
            h(i+1,i) = norm
c
c ......... Normalizacao de g(i+1):
c
            inorm = 1.d0/norm
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)*inorm
  330       continue
c ...........................................................
c
c ... Givenxs QR Methods
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
           do 610 i = 1, neq
             x(i) = x(i) + r * g(i,j)
 610       continue
 600     continue
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
     .         'GMRES2: ',' it ',nit, ' x * Kx ',xkx,' ||x|| ',norm,
     .         ' nKylov ',k,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES2) solver:'/
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
 2300 format (' GMRES2:',5x,'cycles',i7,5x,'It',i7,5x,2d20.10)
 2400 format (' GMRES2:',1x,'Residuo exato > 3.16d0*conv '
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
      real*8  dot,tol,conv,xkx,norm,d,alpha,beta,rr0,w,tmp,norm_r
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
      norm_r = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
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
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,time
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
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BICGSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' BICCSTAB:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c *********************************************************************
c
c **********************************************************************
      subroutine bicgstabl2(neq     ,nequ  ,nad,ia ,ja 
     .                     ,ad      ,au    ,al ,b  ,x   
     .                     ,t       ,v     ,r  ,u  ,r0
     .                     ,w       ,s 
     .                     ,tol     ,maxit  
     .                     ,matvec  ,dot    
     .                     ,my_id   ,neqf1i,neqf2i 
     .                     ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                     ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 01/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * BICGSTABL2: Solucao de sistemas de equacoes pelo metodo dos        * 
c * gradientes biconjugados(2) para matrizes nao-simetricas.           *                                         *
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
c * u(neq)   - arranjo local de trabalho                               *
c * r0(neq)  - arranjo local de trabalho                               *
c * w(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
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
c * Versão do livro Iterative krylov Method for large linear Systems   *                                                                      *
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
      real*8  r(*),u(*),t(*),v(*),r0(*),w(*),s(*)
      real*8  dot,tol,conv,xkx,norm,d
      real*8  alpha,beta,rr0,rr1,omega1,omega2,mi,ni,gamma,tau,tmp
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
c ... Ax
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,t,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - t(i)
c ... r = r0
         r(i)  = r0(i)
c ... u = 0.0d0
         u(i)  = 0.d0 
  100 continue
c .......................................................................
c
c ...
      rr0     = 1.0d0
      alpha   = 0.d0
      omega2  = 1.d0
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit, 1
c ... rro = -w2*rr0
         rr0 = -omega2*rr0 
c ... even BiCG step:
c ... rr1 = (r,r0)
         rr1  = dot(r,r0,neq_doti)
c ... beta = alpha * rr1/ rr0
         beta = alpha*rr1/rr0
c ... rr0 = rr1
         rr0  = rr1
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j)   = r(j) - beta*u(i)
           u(i) = r(i) - beta * u(i)
         enddo        
c .......................................................................
c
c ... v = Au(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,u,v
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...  gamma = (v,r0)
         gamma = dot(v,r0,neq_doti)
c ... alpha = rr0 / gamma
         alpha = rr0/gamma
c .......................................................................
c
c ...
         do i = 1, neq
c ... r(j)   = r(j) - alpha*v(i)
           r(i) = r(i) - alpha * v(i)
         enddo        
c .......................................................................
c
c ... s = Ar(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,r,s
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
         do i = 1, neq
c ... x(j) = x(j) + alpha*u
            x(i) = x(i) + alpha * u(i)
         enddo   
c ........................................................................
c
c ... odd BiCG step:
c ... rr1 = (s,r0)
         rr1  = dot(s,r0,neq_doti)
c ... beta = alpha*rr1/rr0
         beta = alpha*rr1/rr0
c ... rr0 = rr1
         rr0  = rr1
c .......................................................................
c
c ...
         do i = 1, neq
c ... v(j) = s(j) - beta*v
            v(i) = s(i) - beta * v(i)
         enddo   
c ........................................................................
c
c ... w = Av(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,v,w
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... gamma = (w,r0) 
         gamma = dot(w,r0,neq_doti)
c ... alpha = rr0/gamma
         alpha = rr0/gamma
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j)   = r(j) - beta*u(i)
           u(i) = r(i) - beta * u(i)
c ... r(j)   = r(j) - alpha*v(i)
           r(i) = r(i) - alpha * v(i)
c ... s(j)   = s(j) - alpha*w(i)
           s(i) = s(i) - alpha * w(i)
         enddo        
c .......................................................................
c
c ... t = As(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,s,t
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... GCR(2)-part
c ... w1 = (r,s)
         omega1 = dot(r,s,neq_doti)
c ... mi = (s,s)
         mi = dot(s,s,neq_doti)
c ... ni = (s,t)
         ni = dot(s,t,neq_doti)
c ... tau = (t,t)
         tau = dot(t,t,neq_doti)
c ... w2
         omega2 = dot(r,t,neq_doti)
c ... tau = tau - ni*ni/mi
         tau = tau - ni*ni/mi
c ... w2 = (w2 - ni*w1/mi)/tau
         omega2 = (omega2 - ni*omega1/mi)/tau
c ... w1 = (w1 - ni*w2)/mi
         omega1 =(omega1 - ni*omega2)/mi
c .......................................................................
c
c ...
         do i = 1, neq
c ... x(j+2) = x(i) + w1 * r(i) + w2 * s(i) + alpha * u(i)
           x(i) = x(i) + omega1 * r(i) + omega2 * s(i) + alpha * u(i)
c ... r(j+2) = r(j) - w1 * s(i) - w2 * t(i)
           r(i) = r(i) - omega1 * s(i) - omega2 * t(i)
         enddo        
c .......................................................................
c
c ...
         d  = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j+1) = u(j) - w1 * v(i) - w2 * w(i)
           u(i) = u(i) - omega1 * v(i) - omega2 * w(i)
         enddo    
c .......................................................................
c
c ...
         if( jj .eq.500) then
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
     .           ,x,t 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,t,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - t(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
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
     .       "BICGSTAB(2): "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA BICCSTAB(2):',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (BICGSTAB(2)) solver:'/
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
 1300 format (' BICGSTAB(2):',5x,'It',i7,5x,2d20.10)
 1400 format (' BICCSTAB(2):',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c **********************************************************************
c
c **********************************************************************
      subroutine pbicgstabl2(neq     ,nequ  ,nad,ia ,ja 
     .                     ,ad      ,au    ,al  ,m  ,b  ,x   
     .                     ,t       ,v     ,r  ,u   ,r0
     .                     ,w       ,s     ,p  ,h   ,z
     .                     ,tol     ,maxit  
     .                     ,matvec  ,dot    
     .                     ,my_id   ,neqf1i,neqf2i 
     .                     ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                     ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 01/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTABl2  : Solucao de sistemas de equacoes pelo metodo dos     * 
c * gradientes biconjugados para matrizes nao-simetricas com           *
c * prendicionador diagonal                                            *
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
c * u(neq)   - arranjo local de trabalho                               *
c * r0(neq)  - arranjo local de trabalho                               *
c * w (neq)  - arranjo local de trabalho                               *
c * s (neq)  - arranjo local de trabalho                               *
c * p (neq)  - arranjo local de trabalho                               *
c * h (neq)  - arranjo local de trabalho                               *
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
c * Versão do livro Iterative krylov Method for large linear Systems   * 
c * A(M-1)y=b precondicionador a direita                               *                                                                      *
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
      real*8  r(*),u(*),t(*),v(*),r0(*),w(*),s(*),p(*),z(*),h(*)
      real*8  dot,tol,conv,xkx,norm,d,morm_r
      real*8  alpha,beta,rr0,rr1,omega1,omega2,mi,ni,gamma,tau,tmp
      real*8  time0,time
      real*8  dum1 
      real*8  breaktol,btol
      parameter (btol = 1.d-32)
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
c ... Ax
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,t,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - t(i)
c ... r = r0
         r(i)  = r0(i)
c ... u = 0.0d0
         u(i)  = 0.d0 
  100 continue
c .......................................................................
c
c ...
      rr0     = 1.0d0
      alpha   = 0.d0
      omega2  = 1.d0
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit, 1
c ... rro = -w2*rr0
         rr0 = -omega2*rr0 
         if( dabs(rr0) .lt. breaktol) then
           write(*,1510)dabs(rr0)
           call stop_mef() 
         endif  
c ... even BiCG step:
c ... rr1 = (r,r0)
         rr1  = dot(r,r0,neq_doti)
c ... beta = alpha * rr1/ rr0
         beta = alpha*rr1/rr0
c ... rr0 = rr1
         rr0  = rr1
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j) = r(j) - beta*u(i)
           u(i) = r(i) - beta * u(i)
c ... p    = M(-1)u
           z(i)= u(i)*m(i) 
         enddo        
c .......................................................................
c
c ... v = Ap(j) = AM(-1)u
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,z,v
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...  gamma = (v,r0)
         gamma = dot(v,r0,neq_doti)
c ... alpha = rr0 / gamma
         alpha = rr0/gamma
c .......................................................................
c
c ...
         do i = 1, neq
c ... r(j)   = r(j) - alpha*v(i)
           r(i) = r(i) - alpha * v(i)
c ... p    = M(-1)r
           p(i) = r(i)*m(i) 
         enddo        
c .......................................................................
c
c ... s = Ap(j) = AM(-1)r
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,s
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
         do i = 1, neq
c ... x(j) = x(j) + alpha*(M-1)u
            x(i) = x(i) + alpha * z(i)
         enddo   
c ........................................................................
c
c ... odd BiCG step:
c ... rr1 = (s,r0)
         rr1  = dot(s,r0,neq_doti)
c ... beta = alpha*rr1/rr0
         beta = alpha*rr1/rr0
c ... rr0 = rr1
         rr0  = rr1
c .......................................................................
c
c ...
         do i = 1, neq
c ... v(j) = s(j) - beta*v
           v(i) = s(i) - beta * v(i)
c ... p    = M(-1)v
           p(i) = v(i)*m(i) 
         enddo   
c ........................................................................
c
c ... w = Av(j) = AM(-1)v
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,w
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... gamma = (w,r0) 
         gamma = dot(w,r0,neq_doti)
c ... alpha = rr0/gamma
         alpha = rr0/gamma
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j)   = r(j) - beta*u(i)
           u(i) = r(i) - beta * u(i)
c ... r(j)   = r(j) - alpha*v(i)
           r(i) = r(i) - alpha * v(i)
c ... s(j)   = s(j) - alpha*w(i)
           s(i) = s(i) - alpha * w(i)
c ... p    = M(-1)s
           p(i) = s(i)*m(i)
c ... z    = M(-1)u
           z(i) = u(i)*m(i) 
c ... h    = M(-1)r
           h(i) = r(i)*m(i)  
         enddo        
c .......................................................................
c
c ... t = Ap(j) = AM(-1)s
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,t
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... GCR(2)-part
c ... w1 = (r,s)
         omega1 = dot(r,s,neq_doti)
c ... mi = (s,s)
         mi = dot(s,s,neq_doti)
c ... ni = (s,t)
         ni = dot(s,t,neq_doti)
c ... tau = (t,t)
         tau = dot(t,t,neq_doti)
c ... w2
         omega2 = dot(r,t,neq_doti)
c ... tau = tau - ni*ni/mi
         tau = tau - ni*ni/mi
c ... w2 = (w2 - ni*w1/mi)/tau
         omega2 = (omega2 - ni*omega1/mi)/tau
c ... w1 = (w1 - ni*w2)/mi
         omega1 =(omega1 - ni*omega2)/mi
c .......................................................................
c
c ...
         do i = 1, neq
c ... x(j+2) = x(i) + w1 * r(i) + w2 * s(i) + alpha * u(i)
           x(i) = x(i) + omega1 * h(i) + omega2 * p(i) + alpha * z(i)
c ... r(j+2) = r(j) - w1 * s(i) - w2 * t(i)
           r(i) = r(i) - omega1 * s(i) - omega2 * t(i)
         enddo        
c .......................................................................
c
c ...
         d  = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j+1) = u(j) - w1 * v(i) - w2 * w(i)
           u(i) = u(i) - omega1 * v(i) - omega2 * w(i)
         enddo    
c .......................................................................
c
c ...
         if( jj .eq.500) then
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
     .           ,x,t 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,t,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - t(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      morm_r = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
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
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,morm_r,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "PBICGSTAB(2): "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA BICCSTAB(2):',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (PBICGSTAB(2)) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PBICGSTAB(2):',5x,'It',i7,5x,2d20.10)
 1400 format (' PBICCSTAB(2):',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
 1510 format (' PBICGSTAB:',1x,'Breakdown:',1x,'(r0)',1x,d20.10)
      end
c **********************************************************************
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
c * A(M-1)y=b precondicionador a direita                               *
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
      real*8  dot,tol,conv,d,alpha,beta,rr0,w,xkx,norm,tmp,norm_r
      real*8  time0,time
      real*8  breaktol,btol
      parameter (btol = 1.d-32)
      real*8  dum1 
      logical flog,fprint,fnew,f
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
      conv = tol*dsqrt(d)
      breaktol = btol*dsqrt(d)
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
c ... p = r0
         p(i)  = r0(i)
c ... r = r0
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
         if( dsqrt(dabs(rr0)) .lt. breaktol) then
           write(*,1510)dabs(rr0)
           call stop_mef() 
         endif 
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
         if( dabs(w) .lt. breaktol) then
           print*,breaktol
           write(*,1515),dabs(w)
           call stop_mef()
         endif 
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
      norm_r = dsqrt(tmp)
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
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "PBICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
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
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PBICGSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' PBICGSTAB:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
 1510 format (' PBICGSTAB:',1x,'Breakdown:',1x,'(r,r0)',1x,d20.10)
 1515 format (' PBICGSTAB:',1x,'Breakdown:',1x,'w',1x,d20.10)
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
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
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
 1400 format (' ICBICCSTAB:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
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
c * Data de modificaco : 18/07/2016                                    * 
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
     .          ,.false.   ,.false.    ,.false.,.false.)
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
     .          ,.false.,.true.  ,.false.,.false.)
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
c
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
      write(16,*) istep,i,time-time0,time_csr,xkx,norm
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
      subroutine minres(neq   ,nequ  ,nad   ,ia      ,ja
     .             ,ad    ,au    ,al    ,b     ,x
     .             ,z     ,v0    ,v     ,w     ,w0    ,w00    
     .             ,tol   ,maxit
     .             ,matvec,dot
     .             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .             ,i_xfi ,i_rcvsi,i_dspli
     .             ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 16/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * MINRES  : Solucao de sistemas de equacoes pelo metodo MINRES       *
c * (matriz simetrica geral)                                           *
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
c * v0(neq)  - arranjo local de trabalho                               *
c * w(neq)   - arranjo local de trabalho                               *
c * w0(neq)  - arranjo local de trabalho                               *
c * w00(neq) - arranjo local de trabalho                               *
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
c * fonte: Iterative Krylov Methods for Large Linear Systems           *
c * Henk A. van de Vorst                                               *
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
      real*8  v(*),v0(*),w0(*),w00(*),w(*),z(*)
      real*8  dot,tol,conv,xkx,norm,tmp
      real*8  neta,beta_old,beta,c_old,c,s_old,s,normr
      real*8  ro1,ro2,ro3,alpha,delta
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
      tmp  = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(tmp))
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
c ... v1 = b - Ax0
         v(i) = b(i) - z(i)
c ... v0 = 0.0
         v0(i) = 0.d0
c ... w0 = 0.0
         w0(i) = 0.d0
c ... w1 = 0.0
         w00(i) = 0.d0
  100 continue
c ... beta  = ||v1||
      beta  = dsqrt(dot(v,v,neq_doti))
      normr = beta
c ......................................................................
c
c ...
      neta  = beta
      c_old = 1.d0 
      c     = 1.d0
      s_old = 0.d0
      s     = 0.d0
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... The Lancos recurrence:
        tmp = 1.d0/beta
        do 210 i = 1, neq
c ... v(j) = v(j)/beta
          v(i) = v(i)*tmp
  210   continue
c .....................................................................
c
c ... z = Av(j)
        call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .             ,v,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .             ,dum1)
c .....................................................................
c
c ... alpha = ( v(j), z  ) = ( v(j), Av(j))
        alpha = dot(v,z,neq_doti)
c .....................................................................
c
c ...
        do 215 i = 1, neq
c ... v(j+1) = Av(j) - alpha(j)*v(j) - beta(j)*v(j-1)
          z(i) = z(i) - alpha * v(i) - beta * v0(i)
  215   continue
c .....................................................................
c
c ... beta(j+1) = ||v(j+1)||
       beta_old  = beta
       beta      = dsqrt(dot(z,z,neq_doti))
c ......................................................................
c
c ... QR part:
c ... delta = c(j)*alfa(j) - c(j-1)*s(j)*beta(j)
        delta = c*alpha - c_old*s*beta_old
c ... ro2 = s(j)*alpha(j) + c(j-1)*c(j)*beta(j)
        ro2 = s*alpha + c_old*c*beta_old
c ... ro3 = s(j-1)*beta(j)
        ro3 = s_old*beta_old 
c .......................................................................
c
c ...
        c_old = c
        s_old = s
c .......................................................................
c
c ... New Givens rotation for subdiag element:
      call sym_ortho(delta,beta,c,s,ro1)
c ... ro1 = raiz( delta^2 + beta(j+1)^2)
c       ro1 = dsqrt(delta*delta + beta*beta)
c ... c(j+1) = delta*ro1
c       c = delta/ro1
c ... s(j+1) = beta(j+1)/ro1
c       s = beta/ro1
c .........................................................................
c
c ... Update of solution (W = VR^-1) 
        tmp = 1.d0/ro1
        do 220 i = 1, neq
c ... w(j) = (v(j) - ro3*w(j-2) - ro2*w(j-1))/ro1
          w(i) = (v(i) - ro3*w00(i) - ro2*w0(i))*tmp
c ... x(j) = x(j-1) + c(j+1)*neta*w(i)
          x(i)  = x(i) + c*neta*w(i)
c ...
          v0(i) = v(i)
          v(i)  = z(i)
c ...
          w00(i) = w0(i)
          w0(i)  = w(i)   
  220   continue
c .......................................................................
c
c .......................................................................
c
c ... ||r(j) || = |c(j+1)| || r(j-1) ||
        normr = dabs(s)*normr 
        if ( normr .lt. conv) goto 300  
c .......................................................................
c
c ... neta = - s(j+1) * neta
        neta = - s*neta
c .......................................................................
c
c ...
        if( jj .eq. 500) then
          jj = 0
          write(*,1300),j,normr,conv 
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
        v(i) = b(i) - z(i)
  310 continue
      tmp  = dot(v,v,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. 3.16*conv ) then
        write(*,1400) tmp,conv
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
     .       "MINRES: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA MINRES:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (MINRES) solver:'/
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
 1300 format (' MINRES:',5x,'It',i7,5x,2d20.10)
 1400 format (' MINRES:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c *********************************************************************
c
c *********************************************************************  
      subroutine pminres(neq   ,nequ  ,nad   ,ia    ,ja
     .             ,ad    ,au    ,al  ,b     ,x     ,m
     .             ,v0    ,v     ,w   ,w0    ,w00
     .             ,z     ,z0    ,p
     .             ,tol   ,maxit ,nrestart
     .             ,matvec,dot
     .             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .             ,i_xfi ,i_rcvsi,i_dspli
     .             ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 16/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * PMINRES  : Solucao de sistemas de equacoes pelo metodo MINRES      *     
c * precondicionado diagonal M=D(1/2)D(1/2)                            *
c * (matriz simetrica geral)                                           *
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
c * v0(neq)  - arranjo local de trabalho                               *
c * w(neq)   - arranjo local de trabalho                               *
c * w0(neq)  - arranjo local de trabalho                               *
c * w00(neq) - arranjo local de trabalho                               *
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
c * fonte:                                                             *
c * 1. Iterative Krylov Methods for Large Linear Systems               *
c * Henk A. van de Vorst                                               *
c * 2. Precondicionador Iterative methods for singular linear          *  
c * equations and lest-squares problems                                * 
c *( Sou-Cheng (Terrya) Choi - 2006                                    *
c * D(-1/2)AD(-1/2)x = D(-1/2)b                                        *
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
      integer ia(*),ja(*),my_id,it,nrestart
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*)
      real*8  v(*),v0(*),w0(*),w00(*),w(*),z(*),z0(*),p(*)
      real*8  dot,tol,conv,xkx,norm,tmp1,tmp2,tmp3,tmp4
      real*8  neta,beta_old,beta,c_old,c,s_old,s,normr
      real*8  ro1,ro2,ro3,alpha,delta,norm_m_r,norm_r
      real*8  time0,time
      real*8 dum1
      logical flog,fprint,fnew
      external matvec,dot
c ======================================================================
      time0    = MPI_Wtime()
      it       = 0
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
c ... conv = tol * |M(-1/2)b|
      do 15 i = 1, neq
        z(i)  = b(i)*dsqrt(dabs(m(i)))
   15 continue
      tmp1 = dot(z,z,neq_doti)
      conv = tol*dsqrt(dabs(tmp1))
c .......................................................................
c  
   50 continue
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)  
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... v1 = b - Ax0
         v(i) = b(i) - z(i)
c ... v0 = 0.0
         v0(i) = 0.d0
c ... w0 = 0.0
         w0(i) = 0.d0
c ... w1 = 0.0
         w00(i) = 0.d0
c ... z  =M(-1)v
         z0(i) = v(i)*m(i) 
  100 continue
c ... beta  = raiz( vt*z )
      beta  = dsqrt(dabs(dot(v,z0,neq_doti)))
      normr = beta
c ......................................................................
c
c ...
      neta     = beta
      beta_old = beta
      c_old    = 1.d0 
      c        = 1.d0
      s_old    = 0.d0
      s        = 0.d0
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
        it = it + 1
c ... p = Az(j)
        call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .             ,z0,p,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .             ,dum1)
c .....................................................................
c
c ... alpha = ( z(j), p  ) / beta^2        
        tmp1  = 1.d0/(beta*beta)
        alpha = dabs(dot(z0,p,neq_doti))*tmp1
c .....................................................................
c
c ...
        tmp1 =  1.0d0/beta
        tmp2 = -alpha/beta
        tmp3 = -beta/beta_old
        do 215 i = 1, neq
c ... v(j+1) = Az(j) - alpha(j)*v(j) - beta(j)*z(j-1)
          p(i) = tmp1*p(i) + tmp2* v(i) + tmp3 * v0(i)
c ... z  =M(-1)v
          z(i) = p(i)*m(i)
  215   continue
c .....................................................................
c
c ... beta(j+1) = raiz( (zt(j+1),v(j+1)) )
       beta_old  = beta
       tmp1      = dot(z,p,neq_doti)
       beta      = dsqrt(dabs(tmp1))
c ......................................................................
c
c ... QR part:
c ... delta = c(j)*alfa(j) - c(j-1)*s(j)*beta(j)
        delta = c*alpha - c_old*s*beta_old
c ... ro2 = s(j)*alpha(j) + c(j-1)*c(j)*beta(j)
        ro2 = s*alpha + c_old*c*beta_old
c ... ro3 = s(j-1)*beta(j)
        ro3 = s_old*beta_old 
c .......................................................................
c
c ...
        c_old = c
        s_old = s
c .......................................................................
c
c ... New Givens rotation for subdiag element:
        call sym_ortho(delta,beta,c,s,ro1)
c ... ro1 = raiz( delta^2 + beta(j+1)^2)
c       ro1 = dsqrt(delta*delta + beta*beta)
c ... c(j+1) = delta*ro1
c       c = delta/ro1
c ... s(j+1) = beta(j+1)/ro1
c       s = beta/ro1
c .........................................................................
c
c ... Update of solution (W = VR^-1) 
        tmp1 = 1.d0/(beta_old*ro1)
        tmp2 = -ro2/ro1        
        tmp3 = -ro3/ro1
        tmp4 = c*neta
        do 220 i = 1, neq
c ... w(j) = (v(j)/beta(j) - ro2*w(j-1) - ro3*w(j-2))/ro1
c         w(i) = (tmp2*z0(i) - ro2*w0(i) - ro3*w00(i))*tmp1
          w(i) = tmp1*z0(i) + tmp2*w0(i) + tmp3*w00(i)
c ... x(j) = x(j-1) + c(j+1)*neta*w(i)
          x(i)  = x(i) + tmp4*w(i)
c ...
          v0(i) = v(i)
          v(i)  = p(i)
c ...
          z0(i) = z(i)
c ...
          w00(i) = w0(i)
          w0(i)  = w(i)   
  220   continue
c .......................................................................
c
c .......................................................................
c
c ... ||r(j) || = |c(j+1)| || r(j-1) ||
        normr = dabs(s)*normr 
        if ( normr .lt. conv) goto 300  
c .......................................................................
c
c ... neta = - s(j+1) * neta
        neta = - s*neta
c .......................................................................
c
c ...
        if( jj .eq. 500) then
          jj = 0
          write(*,1300)nrestart,it,normr,conv 
        endif  
        jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
      if( nrestart .ne. 0 ) then
        nrestart = nrestart - 1
        goto 50
      endif
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
        v(i) = b(i) - z(i)
        z(i) = v(i)*dsqrt(dabs(m(i)))
  310 continue
      tmp1     = dot(v,v,neq_doti)
      norm_r   = dsqrt(tmp1)
      tmp1     = dot(z,z,neq_doti)
      norm_m_r = dsqrt(tmp1)
      if( norm_m_r .gt. 3.16*conv ) then
        write(*,1400) norm_m_r,conv
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,it,xkx,norm,norm_r,time
      endif
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "PMINRES: "," it ",it, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA MINRES:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (PMINRES) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || M(-1/2)b || = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PMINRES:','nrestart',i7,5x,'It',i7,5x,2d20.10)
 1400 format (' PMINRES:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c **********************************************************************
c 
c *********************************************************************  
      subroutine cr(neq   ,nequ   ,nad   ,ia       ,ja
     .              ,ad    ,au    ,al    ,b        ,x
     .              ,z     ,r     ,p     ,ap       ,ar
     .              ,tol   ,maxit
     .              ,matvec,dot
     .              ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .              ,i_xfi ,i_rcvsi,i_dspli
     .              ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 22/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CR : Solucao de sistemas de equacoes pelo metodo dos               *
c * residuos conjugados para matrizes simetricas.                      *                                                        *
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
c * ar(neq)  - arranjo local de trabalho                               *
c * ap(neq)  - arranjo local de trabalho                               *
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
      real*8  ad(*),au(*),al(*),b(*),x(*)
      real*8  r(*),z(*),p(*),ar(*),ap(*)
      real*8  dot,tol,conv,xkx,norm,d,rar,alpha,beta,tmp
      real*8  norm_r
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
     .           ,x,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
c ... p0 = r0
         p(i) = r(i)
  100 continue
c ......................................................................
c
c ... Ar(0)                                
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,r,ar,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,dum1)
c ......................................................................
c
c ... Ap(0) = Ar(0)
      do 110 i = 1, neq
         ap(i) = ar(i)
  110 continue
c ......................................................................
c
c ... loop         
      jj = 1
      do 230 j = 1, maxit
c ... alpha = (r(j),Ar(j))/(Ap(j),Ap(j))
         alpha = dot(r,ar,neq_doti) / dot(ap,ap,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alfa(j)*p(j)
           x(i) = x(i) + alpha * p(i)
  210    continue
c .....................................................................
c
c ... rar = (r(j),Ar(j))
         rar = dot(r,ar,neq_doti) 
c .....................................................................
c
c ...
         do 215 i = 1, neq
c ...r(j+1) = r(j) - alfa(j)*Ap(j)
           r(i) = r(i) - alpha * ap(i)
  215    continue
c .....................................................................
c
c ... Ar(j+1)                                
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,r,ar,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,dum1)
c ......................................................................
c
c ... beta(j) = (r(j+1),Ar(j+1))/rAr   
          beta = dot(r,ar,neq_doti) / rAr
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... p(j+1) = r(j+1) + beta(j)p(i)
            p(i) = r(i) + beta * p(i)
c ... Ap(j+1) = Ar(j+1) + betaAp(j)
            ap(i) = ar(i) + beta*ap(i)
  220    continue
c .....................................................................
c
c ...
         d = dot(r,r,neq_doti)
         if (dsqrt(d) .lt. conv) goto 300
c ......................................................................
         if( jj .eq.1000) then
           jj = 0
           write(*,1300),j,dsqrt(d),conv 
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
c ... r =(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_r .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_r,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'CR: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA CR:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (CR) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||          = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' CR:',5x,'It',i7,5x,2d20.10)
 1400 format (' CR:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c ********************************************************************* 
c
c *********************************************************************  
      subroutine pcr(neq   ,nequ  ,nad   ,ia   ,ja
     .              ,ad    ,au    ,al    ,b    ,m  ,x
     .              ,z     ,r     ,p     ,ap   ,az
     .              ,t  
     .              ,tol   ,maxit
     .              ,matvec,dot
     .              ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .              ,i_xfi ,i_rcvsi,i_dspli
     .              ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 22/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * PCR : Solucao de sistemas de equacoes pelo metodo dos              *
c * residuos conjugados com precondicionador diagonal para matrizes    *
c * simetricas.                                                        *                                                        *
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
c * az(neq)  - arranjo local de trabalho                               *
c * ap(neq)  - arranjo local de trabalho                               *
c * t(neq)   - arranjo local de trabalho                               *
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
c * M(-1)Ax=M(-1)b onde M tem que simetrica positiva definida para     *
c * gerar um produto interno proprio conjugado com M, ou seja:         *
c * (x,y)M = (x,My)                                                    *
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
      real*8  ad(*),au(*),al(*),b(*),x(*),m(*)
      real*8  r(*),z(*),p(*),az(*),ap(*),t(*)
      real*8  dot,tol,conv,xkx,norm,d,zaz,alpha,beta,tmp
      real*8  norm_r,norm_m_r
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
c ... conv = tol * |M(-1)b|M = tol * (M(-1)b,b) 
      do 20 i = 1, neq
        z(i) = b(i)*m(i)
   20 continue
      d    = dot(b,z,neq_doti)
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
c ... z = M(-1)r0
        z(i) = r(i)*m(i)  
c ... p0 = r0
        p(i) = z(i)
  100 continue
c ......................................................................
c
c ... Az(0)                                
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,z,az,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,dum1)
c ......................................................................
c
c ... Ap(0) = Az(0)
      do 110 i = 1, neq
         ap(i) = az(i)
  110 continue
c ......................................................................
c
c ... loop         
      jj = 1
      do 230 j = 1, maxit
c ... t = M(-1)Ap(j)
        do 205 i = 1, neq
          t(i) = ap(i)*m(i)
  205   continue
c .....................................................................
c
c ... alpha = (z(j),Az(j))/(Ap(j),M(-1)Ap(j))
         alpha = dot(z,az,neq_doti) / dot(ap,t,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alfa(j)*p(j)
           x(i) = x(i) + alpha * p(i)
  210    continue
c .....................................................................
c
c ... zaz = (z(j),Az(j))
         zaz = dot(z,az,neq_doti) 
c .....................................................................
c
c ...
         do 215 i = 1, neq
c ... r(j+1) = r(j) - alfa(j)*Ap(j)
           r(i) = r(i) - alpha * ap(i)
c ... z(j+1) = (-1)r(j+1)
           z(i) = r(i)*m(i)  
  215    continue
c .....................................................................
c
c ... Az(j+1)                                
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,z,az,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,dum1)
c ......................................................................
c
c ... beta(j) = (z(j+1),Az(j+1))/zAz   
          beta = dot(z,az,neq_doti) / zaz
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... p(j+1) = z(j+1) + beta(j)p(i)
            p(i) = z(i) + beta * p(i)
c ... Ap(j+1) = Az(j+1) + betaAp(j)
            ap(i) = az(i) + beta*ap(i)
  220    continue
c .....................................................................
c
c ...
         d = dot(r,z,neq_doti)
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
         if( jj .eq.1000) then
           jj = 0
           write(*,1300),j,sqrt(dabs(d)),conv 
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
c ... r =M(-1)(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
        z(i) = r(i)*m(i)
  310 continue
      norm_m_r = dot(r,z,neq_doti)
      norm_m_r = dsqrt(dabs(norm_m_r))
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_m_r,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,norm_m_r,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'PCR: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCR:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (PCR) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PCR:',5x,'It',i7,5x,2d20.10)
 1400 format (' PCR:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c ********************************************************************* 
c
c *********************************************************************  
      subroutine symmlq(neq   ,nequ  ,nad   ,ia      ,ja
     .             ,ad    ,au    ,al    ,b  ,x
     .             ,v0    ,v     ,y     ,wb   
     .             ,tol   ,maxit
     .             ,matvec,dot
     .             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .             ,i_xfi ,i_rcvsi,i_dspli
     .             ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 30/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * SYMMLQ  : Solucao de sistemas de equacoes pelo metodo SYMMLQ       *
c * (matriz simetrica geral)                                           *
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
c * v0(neq)  - arranjo local de trabalho                               *
c * v(neq)   - arranjo local de trabalho                               *
c * y(neq)   - arranjo local de trabalho                               *
c * wb(neq)  - arranjo local de trabalho                               *
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
c * fonte: http://web.stanford.edu/group/SOL/software/symmlq/          *
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
      real*8  v(*),v0(*),wb(*),y(*)
      real*8  dot,tol,conv,xkx,tmp,alpha,gamma,delta,epsln,norm,norm_r
      real*8  beta1,beta_old,beta,c,s,normr,zi,t1,t2
      real*8  rhs1,rhs2,snprod,tnorm,ynorm2,Anorm,gbar,dbar  
      real*8  epsa,eps,epsr,diag,lqnorm,cgnorm,ynorm,qrnorm
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
      tmp  = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(tmp))
c .......................................................................
c  
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)  
     .           ,x,y 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... v1 = b - Ax0
         v0(i) = b(i) - y(i)
  100 continue
c ... beta  = ||v1||
      beta1  = dsqrt(dot(v0,v0,neq_doti))
      beta_old = beta1  
      beta     = beta_old
c ......................................................................
c
c ... primeiro vetor de lanczos
c     v1 = v1 / ||v1||
      tmp = 1.0d0/beta_old
      do 110 i = 1, neq
        v0(i) = v0(i)*tmp
  110 continue
c ... z = Av
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)  
     .           ,v0,y 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      alpha = dot(y,v0,neq_doti)
c ... vetor de lanczos
      do i = 1, neq
        y(i) = y(i) - alpha*v0(i) 
      enddo
      beta    = dsqrt(dot(y,y,neq_doti))
c .....................................................................
c
c ...
      do 115 i = 1, neq
        wb(i) = v0(i)
  115 continue
c ......................................................................
c
c ...
	rhs1     = beta1
      rhs2     = 0.d0
	snprod   = 1.0d0
      tnorm    = alpha*alpha + beta*beta
      ynorm2   = 0.d0
      Anorm    = dsqrt(tnorm)
      gbar     = alpha
      dbar     = beta
c .....................................................................
      jj = 1
      do 230 j = 1, maxit
c ... The Lancos recurrence:
        tmp = 1.d0/beta
        do 210 i = 1, neq
c ... v(j) = v(j)/beta
          v(i) = y(i)*tmp
  210   continue
c .....................................................................
c
c ... y = Av(j)
        call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .             ,v,y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .             ,dum1)
c .....................................................................
c
c ... alpha = ( v(j), y  ) = ( v(j), Av(j))
        alpha = dot(v,y,neq_doti)
c .....................................................................
c
c ...
        do 215 i = 1, neq
c ... v(j+1) = Av(j) - alpha(j)*v(j) - beta(j)*v(j-1)
          y(i) = y(i) - alpha * v(i) - beta * v0(i)
  215   continue
c .....................................................................
c
c ... beta(j+1) = ||v(j+1)||
       beta_old  = beta
       beta      = dsqrt(dot(y,y,neq_doti))
c ......................................................................
c
c ...
        tnorm  = tnorm  +  alpha* alpha  +  beta_old*beta_old  
     .         +  beta*beta
c .....................................................................
c
c ... proximo plano de rotacao para Q. 
c       gamma = dsqrt(gbar*gbar+beta_old*beta_old)
c       c     =     gbar/gamma
c       s     = beta_old/gamma
        call sym_ortho(gbar,beta_old,c,s,gamma)
c ...
        delta  = c * dbar  +  s * alpha
        gbar   = s * dbar  -  c * alpha
        epsln  =   s * beta
        dbar   = - c * beta
c ... News Givens rotation 
        zi = rhs1 /gamma
        t1 = zi*c
        t2 = zi*s
c .....................................................................
c
c ...  
        do 220 i = 1, neq
c ... 
          x(i)   = x(i) + t1*wb(i)  + t2*v(i)
c ... 
          wb(i)  = s*wb(i) - c*v(i)
c ...
          v0(i) = v(i)
c ...
  220   continue
c .......................................................................
c
c .......................................................................
c ...
        snprod = snprod * s
        ynorm2 = zi*zi   +  ynorm2
        rhs1   = rhs2  -  delta * zi
        rhs2   =       -  epsln * zi
c .....................................................................
c
c ...
        Anorm  = sqrt( tnorm  )
        ynorm  = sqrt( ynorm2 )
        epsa   = Anorm * eps
c       epsr   = Anorm * ynorm * tol
        diag   = gbar
        epsa   = Anorm * eps
        if (diag .eq. 0.d0) diag = epsa
        lqnorm = dsqrt( rhs1*rhs1 + rhs2*rhs2)
        qrnorm = snprod * beta1
        cgnorm = qrnorm*beta/dabs(diag)
        if( cgnorm .le. conv) goto 300
c .....................................................................
c
c ...
        if( jj .eq. 2000) then
          jj = 0
          write(*,1300),j,cgnorm,conv
        endif  
        jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c ...
      if (cgnorm .le. lqnorm ) then
        t1 = rhs1/diag
 	  do i = 1, neq
          x(i)   = x(i) + t1*wb(i)
        enddo
 	endif
c .....................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,y 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
       xkx = dot(x,y,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        v(i) = b(i) - y(i)
  310 continue
      tmp  = dot(v,v,neq_doti)
      norm_r = dsqrt(tmp)
      if( norm_r .gt. 3.16*conv ) then
        write(*,1400) norm_r,conv
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,time
      endif
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "SYMMLQ: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA SYMMLQ:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (SYMMLQ) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' SYMMLQ:',5x,'It',i7,5x,2d20.10)
 1400 format (' SYMMLQ:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c *********************************************************************
c
c *********************************************************************  
      subroutine psymmlq(neq   ,nequ  ,nad ,ia ,ja
     .             ,ad    ,au    ,al  ,b   ,m  ,x   
     .             ,v     ,r1    ,r2  ,y   ,wb ,z   
     .             ,tol   ,maxit
     .             ,matvec,dot
     .             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .             ,i_xfi ,i_rcvsi,i_dspli
     .             ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 30/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * SYMMLQ  : Solucao de sistemas de equacoes pelo metodo SYMMLQ       *
c * (matriz simetrica geral)                                           *
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
c * v(neq)   - arranjo local de trabalho                               *
c * r1(neq)  - arranjo local de trabalho                               *
c * r2(neq)  - arranjo local de trabalho                               *
c * y(neq)   - arranjo local de trabalho                               *
c * wb(neq)  - arranjo local de trabalho                               *
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
c * ad(*),al(*),au(*) e b - inalterados                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * fonte: http://web.stanford.edu/group/SOL/software/symmlq/          *
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
      real*8  ad(*),au(*),al(*),x(*),m(*),b(*)
      real*8  v(*),r1(*),r2(*),wb(*),y(*),z(*)
      real*8  dot,tol,conv,xkx,tmp,alpha,gamma,delta,epsln,norm_r_m
      real*8  beta1,beta_old,beta,c,s,normr,zi,t1,t2,norm,norm_r
      real*8  rhs1,rhs2,snprod,tnorm,ynorm2,Anorm,gbar,dbar  
      real*8  epsa,eps,epsr,diag,lqnorm,cgnorm,ynorm,qrnorm
      real*8  time0,time
      real*8 dum1
      real*8 smachn
      logical flog,fprint,fnew
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
      eps = smachn()
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
c ... conv = tol * (b,M(-1)b)
      do i = 1, neq
c ... y = M(-1)b
        y(i) = m(i)*b(i)
      enddo
      tmp  = dot(b,y,neq_doti)
      conv = tol*dsqrt(dabs(tmp))
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
c ... v1 = b - Ax0
        r1(i) = b(i) - z(i)
c ... y = M(-1)v1
        y(i)  = r1(i)*m(i)
  100 continue
c ... ||v1||M = (v1,M(-1)v1)
      beta1    = dsqrt(abs(dot(r1,y,neq_doti)))
      beta_old = beta1  
      beta     = beta_old
c ......................................................................
c
c ... primeiro vetor de lanczos
c     v1 = M(-1)v1 / ||v1||M
      tmp = 1.0d0/beta_old
      do 110 i = 1, neq
        v(i) = y(i)*tmp
  110 continue
c ... z = Av
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)  
     .           ,v,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      alpha = dot(v,z,neq_doti)
c ... vetor de lanczos
      do i = 1, neq
        y(i) = z(i) - (alpha/beta1)*r1(i) 
      enddo
      beta    = dsqrt(dot(y,y,neq_doti))
c .....................................................................
c
c ... p = M(-1)y
      do 115 i = 1, neq
        r2(i) = y(i)
        y(i)  = r2(i)*m(i)
        wb(i) = v(i)
  115 continue
c ... ||v2||M
      beta    = dsqrt(dabs(dot(r2,y,neq_doti)))
c ......................................................................
c
c ...
	rhs1     = beta1
      rhs2     = 0.d0
	snprod   = 1.0d0
      tnorm    = alpha*alpha + beta*beta
      ynorm2   = 0.d0
      Anorm    = dsqrt(tnorm)
      gbar     = alpha
      dbar     = beta
c .....................................................................
      jj = 1
      do 230 j = 1, maxit
c ... The Lancos recurrence:
        tmp = 1.d0/beta
        do 210 i = 1, neq
c ... v(j) = v(j)/beta
          v(i) = y(i)*tmp
  210   continue
c .....................................................................
c
c ... y = Av(j)
        call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .             ,v,y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .             ,dum1)
c .....................................................................
c
c ... y = Av(j) - (beta(j-1)/beta(j-2))*v(j-1) 
        t1 = - beta/beta_old
        do 215 i = 1, neq
          y(i) = y(i) + t1*r1(i) 
  215   continue
c .....................................................................
c
c ... alpha = ( v(j), y  ) = ( v(j), Av(j))
        alpha = dot(v,y,neq_doti)
c .....................................................................
c
c ... v(j+1) = y - (alfa(j)/beta(j-1))*v(j) 
        t1 = - alpha/beta
        do i = 1, neq
          y(i) = y(i) + t1*r2(i)  
c 
          r1(i) = r2(i)
c
          r2(i) = y(i)
c ... p = M(-1)y
          y(i) = r2(i)*m(i)
        enddo
c .....................................................................
c
c ... beta(j+1) = ||v(j+1)||
        beta_old  = beta
        beta      = dsqrt(dabs(dot(r2,y,neq_doti)))
c ......................................................................
c
c ...
        tnorm  = tnorm  +  alpha* alpha  +  beta_old*beta_old  
     .         +  beta*beta
c .....................................................................
c
c ... proximo plano de rotacao para Q. 
c       gamma = dsqrt(gbar*gbar+beta_old*beta_old)
c       c     =     gbar/gamma
c       s     = beta_old/gamma
        call sym_ortho(gbar,beta_old,c,s,gamma)
c ...
        delta  = c * dbar  +  s * alpha
        gbar   = s * dbar  -  c * alpha
        epsln  =   s * beta
        dbar   = - c * beta
c ... News Givens rotation 
        zi = rhs1 /gamma
        t1 = zi*c
        t2 = zi*s
c .....................................................................
c
c ...  
        do 220 i = 1, neq
c ... 
          x(i)   = x(i) + t1*wb(i)  + t2*v(i)
c ... 
          wb(i)  = s*wb(i) - c*v(i)
  220   continue
c .......................................................................
c
c .......................................................................
c ...
        snprod = snprod * s
        ynorm2 = zi*zi   +  ynorm2
        rhs1   = rhs2  -  delta * zi
        rhs2   =       -  epsln * zi
c .....................................................................
c
c ...
        Anorm  = dsqrt( tnorm  )
        ynorm  = dsqrt( ynorm2 )
        epsa   = Anorm * eps
c       epsr   = Anorm * ynorm * tol
        diag   = gbar
        epsa   = Anorm * eps
        if (diag .eq. 0.d0) diag = epsa
        lqnorm = dsqrt( rhs1*rhs1 + rhs2*rhs2)
        qrnorm = snprod * beta1
        cgnorm = qrnorm*beta/dabs(diag)
        if( cgnorm .le. conv) goto 300
c .....................................................................
c
c ...
        if( jj .eq. 2000) then
          jj = 0
          write(*,1300),j,cgnorm,conv
        endif  
        jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c ...
      if (cgnorm .le. lqnorm ) then
        t1 = rhs1/diag
 	  do i = 1, neq
          x(i)   = x(i) + t1*wb(i)
        enddo
 	endif
c .....................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,y 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
       xkx = dot(x,y,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        v(i) = b(i) - y(i)
c ... y = M(-1)v1
        y(i)  = v(i)*m(i)
  310 continue
      tmp      = dot(v,v,neq_doti)
      norm_r   = dsqrt(tmp)
      tmp      = dot(v,y,neq_doti)
      norm_r_m = dsqrt(abs(tmp))
      if( norm_r .gt. 3.16*conv ) then
        write(*,1400) norm_r,conv
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,norm_r_m,time
      endif
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "PSYMMLQ: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PSYMMLQ:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (PSYMMLQ) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||M       = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||M        = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PSYMMLQ:',5x,'It',i7,5x,2d20.10)
 1400 format (' PSYMMLQ:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c *********************************************************************
c
c *********************************************************************  
      subroutine sqrm(neq   ,nequ   ,nad   ,ia  ,ja
     .                 ,ad    ,au     ,al    ,m   ,b  ,x  
     .                 ,t     ,r     ,q   ,d   
     .                 ,tol   ,maxit
     .                 ,matvec,dot
     .                 ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                 ,i_xfi ,i_rcvsi,i_dspli
     .                 ,fprint,flog   ,fhist  ,fnew)
c **********************************************************************
c * Data de criacao    : 28/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
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
c * fhist    - log dos resuduos por iteracao                           *
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
c * Fonte: A New Krylov-subspace method for symmetric indefinite       * 
c * linear systems                                                     *
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
      real*8 ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8 r(*),t(*),q(*),d(*)
      real*8 dot,tol,conv,xkx,norm,alpha,beta,tmp1,tmp2,tau,ro,vn,v0
      real*8 sigma,cn,norm_b 
      real*8 norm_r,norm_m_r
      real*8 time0,time
      real*8 dum1
      logical flog,fprint,fnew,fhist
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
      norm_b = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(norm_b))
c .......................................................................
c  
c ... Ax0                                                            
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,t,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - t(i)
c ... q = t
         q(i) = r(i) 
c ... d = 0.0
         d(i) = 0.d0
  100 continue
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
     .              ,dum1)
c .....................................................................
c
c ... sigma = ( q(j-1),t)
         sigma = dot(q,t,neq_doti)
         if( sigma .eq. 0.0) then
           print*,"RSQRM fail (sigma)!"
           call stop_mef()  
         endif  
c .....................................................................
c
c ... alpha(j-1) = ro(j-1)/sigma(j-1)
         alpha = ro/sigma
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... r(j) = r(j-1) - alpha(j-1)*t
            r(i) = r(i) - alpha * t(i)
  210    continue
c .....................................................................
c
c ... v(j) = ||r||/tau(j-1)
         vn   = dsqrt(dot(r,r,neq_doti))/tau
c ... c(j) = 1/sqrt(1+v(j)*v(j) 
         cn   = 1.0d0/dsqrt(1+vn*vn)
c ... tau(j) = tau(j-1)*v(j)*c(j)
         tau = tau*vn*cn
c .....................................................................
c
c ... tau(j) = (c(j)*c(j)*v(j-1)*v(j-1)) d(j-1)
c            + c(j)*c(j)*alpha(j-1) * q(j-1)
         tmp1 = cn*cn*v0*v0
         tmp2 = cn*cn*alpha 
         do 215 i = 1, neq
c ... d(j) = (cj*cj*vj*vj) d(j-1) + cj*cj*alpha(j-1) * q(j-1)
           d(i) = tmp1*d(i) + tmp2*q(i) 
c ... x(j) = x(j-1) + d(j)
           x(i) = x(i) + d(i) 
  215    continue 
c .....................................................................
c
c ...
         v0 = vn
c ......................................................................
c
c ... 
         if( ro .eq. 0.0) then
           print*,"RSQRM fail (ro)!"
           call stop_mef()  
         endif  
c .....................................................................
c
c ... (r,u) 
         tmp1 = dot(r,r,neq_doti) 
c ... beta = (r,u)/ro(j-1)
         beta = tmp1/ro
c ...
         ro = tmp1
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... q(j+1) = r(j) + beta*q(j-1)
            q(i) = r(i) + beta * q(i)
  220    continue
c .....................................................................
c
c ...
         norm_r = dsqrt(dot(r,r,neq_doti))
         if(fhist) write(18,1500),j,norm_r/norm_b 
c .....................................................................
c
c ...
         if (norm_r .lt. conv) goto 300
c ......................................................................
         if( jj .eq.1000) then
           jj = 0
           write(*,1300),j,norm_r,conv 
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
     .           ,x,t,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx = dot(x,t,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - t(i)
  310 continue
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_r,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'SMRQ: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA SMRQ:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (SMRQ) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||          = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' SMRQ:',5x,'It',i7,5x,2d20.10)
 1400 format (' SMRQ:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
 1500 format ( 'SMRQ: ',5x,i7,5x,2es20.10)
      end
c **********************************************************************
c
c *********************************************************************  
      subroutine lpsqrm(neq   ,nequ   ,nad   ,ia  ,ja
     .                 ,ad    ,au     ,al    ,m   ,b  ,x  
     .                 ,t     ,r     ,q   ,d   
     .                 ,tol   ,maxit
     .                 ,matvec,dot
     .                 ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                 ,i_xfi ,i_rcvsi,i_dspli
     .                 ,fprint,flog   ,fhist  ,fnew)
c **********************************************************************
c * Data de criacao    : 28/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * LPSQRM : Solucao de sistemas de equacoes pelo metodo QMR simetrico *
c * diagonal a esquerda                                                *
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
c * fhist    - log dos resuduos por iteracao                           *
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
c * Fonte: A New Krylov-subspace method for symmetric indefinite       * 
c * linear systems                                                     *
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
      real*8 ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8 r(*),t(*),q(*),d(*)
      real*8 dot,tol,conv,xkx,norm,alpha,beta,tmp1,tmp2,tau,ro,vn,v0
      real*8 sigma,cn,norm_r,norm_m_r,norm_b
      real*8 time0,time
      real*8 dum1
      logical flog,fprint,fnew,fhist
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
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
         t(i) = b(i) * m(i)
   15 continue
      norm_b = dot(t,t,neq_doti)
      conv   = tol*dsqrt(norm_b)
c .......................................................................
c  
c ... Ax0                                                            
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,t,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - t(i)
c ... t = (M-1)r0
         t(i) = r(i) * m(i)
c ... q = t
         q(i) = t(i)
c ... d = 0.0
         d(i) = 0.d0
  100 continue
c ... ( t,t ) 
      tau = dsqrt(dot(t,t,neq_doti))
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
     .              ,dum1)
c .....................................................................
c
c ... sigma = ( q(j-1),t)
         sigma = dot(q,t,neq_doti)
         if( sigma .eq. 0.0) then
           print*,"lSQRM fail (sigma)!"
           stop  
         endif  
c .....................................................................
c
c ... alpha(j-1) = ro(j-1)/sigma(j-1)
         alpha = ro/sigma
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... r(j) = r(j-1) - alpha(j-1)*t
            r(i) = r(i) - alpha * t(i)
c ... t  = (M-1)r(j)
            t(i) = r(i) * m(i)
  210    continue
c .....................................................................
c
c ... v(j) = ||t||/tau(j-1)
         vn   = dsqrt(dot(t,t,neq_doti))/tau
c ... c(j) = 1/sqrt(1+v(j)*v(j) 
         cn   = 1.0d0/dsqrt(1+vn*vn)
c ... tau(j) = tau(j-1)*v(j)*c(j)
         tau = tau*vn*cn
c .....................................................................
c
c ... tau(j) = (c(j)*c(j)*v(j-1)*v(j-1)) d(j-1)
c            + c(j)*c(j)*alpha(j-1) * q(j-1)
         tmp1 = cn*cn*v0*v0
         tmp2 = cn*cn*alpha 
         do 215 i = 1, neq
c ... d(j) = (cj*cj*vj*vj) d(j-1) + cj*cj*alpha(j-1) * q(j-1)
           d(i) = tmp1*d(i) + tmp2*q(i) 
c ... x(j) = x(j-1) + d(j)
           x(i) = x(i) + d(i) 
  215    continue 
c .....................................................................
c
c ...
         v0 = vn
c ......................................................................
c
c ... 
         if( ro .eq. 0.0) then
           print*,"lSQRM fail (ro)!"
           stop  
         endif  
c .....................................................................
c
c ... (r,t) 
         tmp1 = dot(r,t,neq_doti) 
c ... beta = (r,t)/ro(j-1)
         beta = tmp1/ro
c ...
         ro = tmp1
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... q(j+1) = (M-1)r(j) + beta*q(j-1) = t + beta*q(j-1)
            q(i) = t(i) + beta * q(i)
  220    continue
c .....................................................................
c
c ...
         norm_r = dsqrt(dot(r,r,neq_doti))
         if(fhist) write(18,1500),j,norm_r/norm_b
c .....................................................................
c
c ...
         if (norm_r .lt. conv) goto 300
c ......................................................................
         if( jj .eq.1000) then
           jj = 0
           write(*,1300),j,norm_r,conv 
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
     .           ,x,t,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx = dot(x,t,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =M(-1)(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - t(i)
        q(i) = r(i)*m(i)
  310 continue
      norm_m_r = dot(q,q,neq_doti)
      norm_m_r = dsqrt(norm_m_r)
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_m_r,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,norm_m_r,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'LPSMRQ: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA LPSMRQ:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9,d20.10)
 1100 format(' (LPSMRQ) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||M(-1)b||     = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| M(-1)(b - Ax) ||  = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' LPSMRQ:',5x,'It',i7,5x,2d20.10)
 1400 format (' LPSMRQ:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
 1500 format ( 5x,i7,5x,2es20.10)
      end
c ********************************************************************* 
c
c *********************************************************************  
      subroutine rpsqrm(neq   ,nequ   ,nad   ,ia  ,ja
     .                 ,ad    ,au     ,al    ,m   ,b  ,x  
     .                 ,t     ,r     ,q   ,d   
     .                 ,tol   ,maxit
     .                 ,matvec,dot
     .                 ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                 ,i_xfi ,i_rcvsi,i_dspli
     .                 ,fprint,flog   ,fhist  ,fnew)
c **********************************************************************
c * Data de criacao    : 28/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
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
c * fhist    - log dos resuduos por iteracao                           *
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
c * Fonte: A New Krylov-subspace method for symmetric indefinite       * 
c * linear systems                                                     *
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
      real*8 ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8 r(*),t(*),q(*),d(*)
      real*8 dot,tol,conv,xkx,norm,alpha,beta,tmp1,tmp2,tau,ro,vn,v0
      real*8 sigma,cn,norm_b 
      real*8 norm_r
      real*8 time0,time
      real*8 dum1
      logical flog,fprint,fnew,fhist
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
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
      norm_b = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(norm_b))
c .......................................................................
c  
c ... Ax0                                                            
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,t,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - t(i)
c ... q = t
         q(i) = r(i) * m(i)
c ... d = 0.0
         d(i) = 0.d0
  100 continue
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
     .              ,dum1)
c .....................................................................
c
c ... sigma = ( q(j-1),t)
         sigma = dot(q,t,neq_doti)
         if( sigma .eq. 0.0) then
           print*,"RSQRM fail (sigma)!"
           call stop_mef()  
         endif  
c .....................................................................
c
c ... alpha(j-1) = ro(j-1)/sigma(j-1)
         alpha = ro/sigma
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... r(j) = r(j-1) - alpha(j-1)*t
            r(i) = r(i) - alpha * t(i)
  210    continue
c .....................................................................
c
c ... v(j) = ||r||/tau(j-1)
         vn   = dsqrt(dot(r,r,neq_doti))/tau
c ... c(j) = 1/sqrt(1+v(j)*v(j) 
         cn   = 1.0d0/dsqrt(1+vn*vn)
c ... tau(j) = tau(j-1)*v(j)*c(j)
         tau = tau*vn*cn
c .....................................................................
c
c ... tau(j) = (c(j)*c(j)*v(j-1)*v(j-1)) d(j-1)
c            + c(j)*c(j)*alpha(j-1) * q(j-1)
         tmp1 = cn*cn*v0*v0
         tmp2 = cn*cn*alpha 
         do 215 i = 1, neq
c ... d(j) = (cj*cj*vj*vj) d(j-1) + cj*cj*alpha(j-1) * q(j-1)
           d(i) = tmp1*d(i) + tmp2*q(i) 
c ... x(j) = x(j-1) + d(j)
           x(i) = x(i) + d(i) 
c ... u(j) = M(-1)r
           t(i) = r(i)*m(i)
  215    continue 
c .....................................................................
c
c ...
         v0 = vn
c ......................................................................
c
c ... 
         if( ro .eq. 0.0) then
           print*,"RSQRM fail (ro)!"
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
         do 220 i = 1, neq
c ... q(j+1) = (M-1)r(j) + beta*q(j-1) = t + beta*q(j-1)
            q(i) = t(i) + beta * q(i)
  220    continue
c .....................................................................
c
c ...
         norm_r = dsqrt(dot(r,r,neq_doti))
         if(fhist) write(18,1500),j,norm_r/norm_b 
c .....................................................................
c
c ...
         if (norm_r .lt. conv) goto 300
c ......................................................................
         if( jj .eq.1000) then
           jj = 0
           write(*,1300),j,norm_r,conv 
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
     .           ,x,t,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx = dot(x,t,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - t(i)
  310 continue
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_r .gt. conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_r,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'RPSMRQ: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA LPSMRQ:',/,5x,'Diagonal coefficient ' 
     . '- equation ',i9,d20.10)
 1100 format(' (RPSMRQ) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||          = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' RPSMRQ:',5x,'It',i7,5x,2d20.10)
 1400 format (' RPSMRQ:',1x,'Explicit residual > tol * ||b|| :'
     .       ,1x,d20.10,1x,d20.10)
 1500 format ( 5x,i7,5x,2es20.10)
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
c **********************************************************************
      subroutine sym_ortho(a,b,c,s,r)
c **********************************************************************
c * Data de criacao    : 18/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * SYM_ORTHO: Givens rotation                                        * 
c * (versa com melhor comportamento numerico)                          *      
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * a   - paramentro                                                   *
c * b   - paramentro s da diagonal principal                           *
c * c   - nao definido                                                 *
c * s   - nao definido                                                 *
c * r   - nao definido                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * c   - cos(teta)                                                    *
c * s   - seno(teta)                                                   *
c * r   - raiz(a^2 + b^2)                                              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * | c  s | | a |    | raiz(a^2 + b^2) |   | r |                      *
c   |      | |   |  = |                 | = |   |                      *
c * | s -c | | b |    |        0        |   | 0 |                      *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
      real*8 a,b,c,s,r,t,ma,mb,sa,sb
      real*8 sign1
c ...      
      ma = dabs(a)
      mb = dabs(b)
      sa = sign1(a)
      sb = sign1(b)
c .....................................................................
c
c ...        
      if(b .eq. 0.d0) then
        s = 0.d0
        r = ma
        if ( a .eq. 0.d0) then
          c = 1.0d0
        else
          c = sa
        endif
      else if( a .eq. 0.d0) then
        c = 0.d0
        s = sb
        r = mb
      else if( mb .gt. a) then
        t = a/b
        s = sb/dsqrt(1.d0 + t*t)
        c = s*t
        r = b/s
      else if( ma .gt. mb ) then
        t = b/a
        c = sa/dsqrt(1.d0+t*t)
        s = c*t
        r = a/c
      endif
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 18/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * SIGN1 : retorna o sinal de a                                       * 
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * a   - paramentro                                                   *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * 0.0d0, 1.d0 ou -1.d0
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      real*8 function sign1(a)
      implicit none
      real*8 a
      sign1 = 0.0d0
      if( a .gt. 0.d0) then
        sign1 = 1.d0
      else if( a .lt. 0.d0) then
        sign1 = -1.d0
      endif
      return
      end
c ***********************************************************************
c
c ***********************************************************************
c * Data de criacao    : 30/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * SYM_ORTHO: Givens rotation                                         * 
c * (versa com melhor comportamento numerico)                          *      
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * a   - paramentro                                                   *
c * b   - paramentro s da diagonal principal                           *
c * c   - nao definido                                                 *
c * s   - nao definido                                                 *
c * r   - nao definido                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * c   - cos(teta)                                                    *
c * s   - seno(teta)                                                   *
c * r   - raiz(a^2 + b^2)                                              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * | c  s | | a |    | raiz(a^2 + b^2) |   | r |                      *
c   |      | |   |  = |                 | = |   |                      *
c * |-s  c | | b |    |        0        |   | 0 |                      *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine sym_ortho2(a,b,c,s,r)
      implicit none
      real*8 a,b,c,s,r,t,ma,mb
      real*8 sign1
c ...      
      ma = dabs(a)
      mb = dabs(b)
c .....................................................................
c
c ...     
      c = 1.d0
      s = 0.d0   
      r = a
      if(b .ne. 0.d0) then
        if( mb .gt. ma) then
          t = a/b
          s = 1.d0/dsqrt(1.d0 + t*t)
          c = s*t
          r = b/s
        else 
          t = b/a
          c = 1.d0/dsqrt(1.d0+t*t)
          s = c*t
          r = a/c
        endif
      endif
      return
      end
c **********************************************************************
 
