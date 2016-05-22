      subroutine cyclic_jacobi(a,x,w,ip,n,tol,maxciclos,cflag)
c **********************************************************************
c * Data de criacao    : 14/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * CYCLIC_JACOBI : Calculo de autovalores e autovetores de matrizes   * 
c * simetricas pelo metodo de JACOBI.                                  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * a(n,n) - matriz simetrica                                          *
c * x(n,n) - nao definido                                              *
c * w(n)   - nao definido                                              *
c * ip(n)  - arranjo auxiliar (nao definido)                           *
c * n      - dimensao da matriz                                        *
c * tol    - tolerancia de convergencia                                *
c * maxciclos - numero maximo de ciclos                                *
c * cflag     -.true. ordena os autovalores e vetores em ordem         *
c *             crescente                                              *
c *            .false. ordena os autovalores e vetores em ordem        *
c *             decrescente                                            *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * a - matriz espectral: a(i,i) = i-esimo autovalor                   *
c * x - matriz de autovetores: x(i,j),(i=1,n) = autovetor              *
c *                            associado ao j-esimo autovalor          *
c * w - w(i) = (a(i,i) = i-esimo autovalor)                            *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Versao: Cyclic-by-row Jacobi                                       *
c * fonte: Matrix Computation Golub e Van Loan pag. 480 Ed.4           * 
c **********************************************************************
      implicit none
      integer k
      integer ip(n),n,p,q,maxciclos
      real*8 a(n,n),x(n,n),w(n),tol,conv,off
      real*8 c,s,d
      real*8 forbenius_norm_sym,off_forbenius_norm_sym
      logical cflag
c ......................................................................
c
c ... Inicializa a matriz x (matriz identidade):
c
      call initx(x,n)
c
c ...
      off = off_forbenius_norm_sym(a,n)
      if( off .eq. 0.d0) goto 2000
c .....................................................................
c
c ... Ciclos de Jacobi:
c
      conv = tol*forbenius_norm_sym(a,n)
c ...
      do 1000 k = 1, maxciclos
c ...
        do  110 p = 1, n
c ...
          do  100 q = p+1, n
c
c ... Calculo da rotacao para o coeficiente a(p,q): 
c ... 2-by-2 Symmetric Schur Decomposition
c
            call sym_schur2(a(q,q),a(p,p),a(p,q),c,s)
c .....................................................................
c                   
c ... Produto Xt A X:
c
            call xax(a,c,s,n,p,q)
c
c ... Produto (Xi)(Xi+1):
c                   
            call xx(x,c,s,n,p,q)
  100     continue
  110   continue
        off =  off_forbenius_norm_sym(a,n)
        if(off .lt. conv) goto 2000
 1000 continue
      print*,' CYCLIC_JACOBI: Sem convergecia apos '
     .      ,maxciclos,' ciclos !'
      stop
 2000 continue
c .....................................................................
c
c ...    Ordena os autovalores e autovetores:
c
c ... ordem crescente
      if(cflag) then
        call ordena(a,x,w,ip,n)
c .....................................................................
c
c ... decrescente
      else
        call dec_ordena(a,x,w,ip,n)
      endif
c .....................................................................
      return
      end
c ********************************************************************** 
c
c ********************************************************************** 
      subroutine jacobi(a,x,w,ip,n,tol,maxciclos)
c **********************************************************************
c *                                                                    *
c *   JACOBI                                                           *
c *   ------                                                           *
c *                                                                    *
c *   Calculo de autovalores e autovetores de matrizes simetricas      *
c *   pelo metodo de JACOBI.                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   a(n,n) - matriz simetrica                                        *
c *   x(n,n) - nao definido                                            *
c *   w(n)   - nao definido                                            *
c *   ip(n)  - arranjo auxiliar (nao definido)                         *
c *   n      - dimensao da matriz                                      *
c *   tol    - tolerancia de convergencia                              *
c *   maxciclos - nmumero maximo de ciclos                             *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *   a - matriz espectral: a(i,i) = i-esimo autovalor                 *
c *   x - matriz de autovetores: x(i,j),(i=1,n) = autovetor            *
c *                              associado ao j-esimo autovalor        *
c *   w - w(i) = (a(i,i) = i-esimo autovalor                           *
c *                                                                    *
c **********************************************************************
      implicit none
      integer k
      integer ip(n),n,p,q,maxciclos
      real*8 a(n,n),x(n,n),w(n),tol,diag,off,conv
      real*8 t,c,s,d,r,r1,r2
c ......................................................................
c
c ... Inicializa a matriz x (matriz identidade):
c
      call initx(x,n)
c
c ... Ciclos de Jacobi:
c
      call diagnorm(a,n,diag)
      do 1000 k = 1, maxciclos
        do  110 p = 1, n
          do  100 q = p+1, n
c
c ... Calculo da rotacao para o coeficiente a(p,q):
c
            c = 1.d0
            s = 0.d0
            if(a(p,q) .ne. 0.d0) then
              r = (a(q,q)-a(p,p))/(2.d0*a(p,q))
              d = dsqrt(1.d0+r*r)
              r1 = -r + d
              r2 = -r - d
              t = r1
              if(dabs(r2).lt.dabs(r1)) then
                t = r2
              endif
              c = (1.d0/dsqrt(1.d0+t*t))
              s = (t*c)
            endif
c                   t
c ... Produto X A X:
c
            call xax(a,c,s,n,p,q)
c
c ... Produto (Xi)(Xi+1):
c                   
            call xx(x,c,s,n,p,q)
  100     continue
  110   continue
        call offnorm(a,n,off)
        conv = off/diag
        if(conv .lt. tol) goto 2000
 1000 continue
      print*,'JACOBI: Sem convergecia apos ',maxciclos,' ciclos !'
      stop
 2000 continue
c
c ...    Ordena os autovalores e autovetores:
c
      call ordena(a,x,w,ip,n)
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine xax(a,c,s,n,p,q)
c **********************************************************************
c *                                                                    *
c *      XAX - Produto Xt(AX)                                          *
c *   ------                                                           *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j
      integer n,p,q
      real*8  a(n,n),c,s,a1,a2
c ......................................................................
c
c ... Produto AX:
c
      do 100 i = 1, n
        a1 = a(i,p)
        a2 = a(i,q)
        a(i,p) = a1*c - a2*s
        a(i,q) = a1*s + a2*c
  100 continue
c              
c ... Produto Xt(AX):
c
      do 200 j = 1, n
        a1 = a(p,j)
        a2 = a(q,j)
        a(p,j) = c*a1 - s*a2
        a(q,j) = s*a1 + c*a2
  200 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine xx(x,c,s,n,p,q)
c **********************************************************************
c *                                                                    *
c *      XX  Produto XX                                                *
c *   --                                                               *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      integer n,p,q
      real*8  x(n,n),c,s,a1,a2
c ......................................................................
      do 100 i = 1, n
        a1 = x(i,p)
        a2 = x(i,q)
        x(i,p) =  a1*c - a2*s
        x(i,q) =  a1*s + a2*c
  100 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine initx(x,n)
c **********************************************************************
c *                                                                    *
c *      INITX                                                         *
c *   -----                                                            *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j
      integer n
      real*8  x(n,n)
c ......................................................................
      do 110 i = 1, n
         x(i,i) = 1.d0
         do 100 j = i+1, n
            x(i,j) = 0.d0
            x(j,i) = 0.d0
  100    continue
  110 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine ordena(a,x,w,ip,n)
c **********************************************************************
c *                                                                    *
c *   ORDENA: ordem crescente                                          *
c *   ------                                                           *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j,k
      integer ip(n),n,iaux
      real*8  a(n,n),x(n,n),w(n),aux
      logical itroca
c ......................................................................
      do 100 i = 1, n
         ip(i) = i
         w(i)  = a(i,i)
  100 continue
  200 continue
      itroca = .false.
      do 210 i = 2, n
         j = i-1
         if(w(j).gt.w(i)) then
            aux   = w(i)
            w(i) = w(j)
            w(j) = aux
            iaux  = ip(i)
            ip(i) = ip(j)
            ip(j) = iaux
            do 110 k = 1, n
               aux    = x(k,i)
               x(k,i) = x(k,j)
               x(k,j) = aux
  110       continue
            itroca = .true.
         endif
  210 continue
      if(itroca) goto 200
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine dec_ordena(a,x,w,ip,n)
c **********************************************************************
c *                                                                    *
c *   DEC_ORDENA : ordem decrescente                                   *
c *   ----------                                                       *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j,k
      integer ip(n),n,iaux
      real*8  a(n,n),x(n,n),w(n),aux
      logical itroca
c ......................................................................
      do 100 i = 1, n
         ip(i) = i
         w(i)  = a(i,i)
  100 continue
  200 continue
      itroca = .false.
      do 210 i = 2, n
         j = i-1
         if(w(j).lt.w(i)) then
            aux   = w(i)
            w(i) = w(j)
            w(j) = aux
            iaux  = ip(i)
            ip(i) = ip(j)
            ip(j) = iaux
            do 110 k = 1, n
               aux    = x(k,i)
               x(k,i) = x(k,j)
               x(k,j) = aux
  110       continue
            itroca = .true.
         endif
  210 continue
      if(itroca) goto 200
      return
      end
c **********************************************************************
c
c **********************************************************************
       subroutine diagnorm(a,n,diag)
c **********************************************************************
c *                                                                    *
c *      DIAGNORM                                                      *
c *   --------                                                         *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      integer n
      real*8  a(n,n),diag
c ......................................................................
      diag = 0.d0
      do 100 i = 1, n
         diag = diag + a(i,i)*a(i,i)
  100 continue
      return
      end
c **********************************************************************
      real*8 function forbenius_norm_sym(a,n)
c **********************************************************************
c * Data de criacao    : 14/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * FORBENIUS_NORM_SYM : norma de Forbenius                            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * a - matriz simetrica                                               *
c * n - numero de linhas                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      integer i,j
      integer n
      real*8  a(n,n),norm
c ......................................................................
      norm = 0.d0
      do 100 j = 1, n
        do 110 i = j+1,n
          norm = norm + a(i,j)*a(i,j)
  110   continue
  100 continue
      forbenius_norm_sym = dsqrt(2.d0*norm)
      return
      end
c **********************************************************************
      subroutine offnorm(a,n,off)
c **********************************************************************
c *                                                                    *
c *      OFFNORM                                                       *
c *   --------                                                         *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j
      integer n
      real*8  a(n,n),off
c ......................................................................
      off = 0.d0
      do 110 i = 1, n
        do 100 j = i+1,n
          off = off + a(i,j)*a(i,j)
  100   continue
  110 continue
      off = off*2.d0
      return
      end
c **********************************************************************
c
c **********************************************************************
      real*8 function off_forbenius_norm_sym(a,n)
c **********************************************************************
c * Data de criacao    : 14/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * OFF_FORBENIUS_NORM_SYM : norma de Forbenius dos elementos fora     *    
c * da diagonal                                                        *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * a - matriz simetrica                                               *
c * n - numero de linhas                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      integer i,j
      integer n
      real*8  a(n,n),off
c ......................................................................
      off = 0.d0
      do 110 i = 1, n
        do 100 j = i+1,n
          off = off + a(i,j)*a(i,j)
  100   continue
  110 continue
      off = off*2.d0
      off_forbenius_norm_sym = dsqrt(off)
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine sym_schur2(aqq,app,apq,c,s)
c **********************************************************************
c * Data de criacao    : 19/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * SYM_SCHUR2 :                                                       *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * aqq - matriz simetrica                                             *
c * app - matriz simetrica                                             *
c * apq - matriz simetrica                                             *
c * c   - nao definido                                                 *
c * s   - nao defindido                                                *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * c   - nao definido                                                 *
c * s   - nao definido                                                 *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      real*8 aqq,app,apq,c,s,r,d,t
      c = 1.d0
      s = 0.d0
      if(apq .ne. 0.d0) then
        r = 0.5d0*(aqq-app)/apq
        d = dsqrt(1.d0+r*r)
        if( r .ge. 0.d0) then
          t = 1.d0/(r+d)
        else
          t = 1.d0/(r-d)
        endif
        c = (1.d0/dsqrt(1.d0+t*t))
        s = t*c
      endif
      return
      end
c **********************************************************************
      subroutine xaxold(a,c,s,n,p,q)
c **********************************************************************
c *                                                                    *
c *      XAX                                                           *
c *   ---                                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j
      integer n,p,q
      real*8  a(n,n),c,s,a1,a2
c ......................................................................
c
c ... Produto AX:
c
      do 100 i = 1, n
        a1 = a(i,p)*c - a(i,q)*s
        a2 = a(i,p)*s + a(i,q)*c
        a(i,p) = a1
        a(i,q) = a2
  100 continue
c                   t
c ... Produto X(AX):
c
      do 200 j = 1, n
        a1 = c*a(p,j) - s*a(q,j)
        a2 = s*a(p,j) + c*a(q,j)
        a(p,j) = a1
        a(q,j) = a2
  200 continue
      return
      end
      subroutine xxold(x,c,s,n,p,q)
c **********************************************************************
c *                                                                    *
c *      XX                                                            *
c *   --                                                               *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      integer n,p,q
      real*8  x(n,n),c,s,a1,a2
c ......................................................................
      do 100 i = 1, n
        a1 =  x(i,p)*c - x(i,q)*s
        a2 =  x(i,p)*s + x(i,q)*c
        x(i,p) = a1
        x(i,q) = a2
  100 continue
      return
      end

