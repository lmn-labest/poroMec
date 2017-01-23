c **********************************************************************
c *                                                                    *
c *   MATVEC.F                                           31/08/2005    *
c *                                                                    *
c *   Este arquivo contem subrotinas para produto matriz-vetor e       *
c *   produto escalar de vetores                                       *
c *                                                                    *
c *   matvec_csr                                                       *
c *   matvec_csr_sym_v3                                                *
c *   matvec_csrsym                                                    *
c *   matvec_csrsym1*         (loops aninhados)                        *
c *   matvec_csrc                                                      *
c *   matvec_csrc1                                                     *
c *   matvec_csrc2                                                     *
c *   matvec_csrcsym                                                   *
c *   matvec_csrcsym1                                                  *
c *   matvec_csrcsym2                                                  *
c *   matvec_csrcr                                                     *
c *   matvec_csrcr1                                                    *
c *   matvec_csrcrsym                                                  *
c *   matvec_csrcrsym1                                                 *
c * --------------------------Poromecanico --------------------------- *
c *   matvec_csrc_pm                                                   *
c *   matvec_csrcr_sym_pm                                              *
c *   matvec_csrc_sym_pm                                               *
c *   matvec_csr_sym_pm                                                *
c * ------------------------------------------------------------------ *
c *   saxpb                                                            *
c *   dot                                                              *
c *   dot_par                                                          *
c *   ddot                                                             *
c *   aequalb                                                          *
c *   matvec_sym                                                       *
c *                                                                    *
c *   1 = loop interno desenrolado                                     *
c *   2 = loops desenrolados                                           *
c *                                                                    *
c **********************************************************************
      subroutine matvec_csr(neq,ia,ja,ad,a,al,x,y,neqf1i,neqf2i,
     .                      i_fmapi,i_xfi,i_rcvsi,i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSR: produto matriz-vetor y = Ax  (A nao-simetrica),      *
c *   ----------                      coef. de A no formato CSR        *
c *                                    e grafo nao-simetrico.          *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro       *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor a                               *
c *   ad(neq) - diagonal da matriz A                                   *
c *   a(nad)  - coef. de A, sem a diagonal                             *
c *   al(*)   - nao utilizado                                          *
c *   x(neq+1)- vetor a ser multiplicado                               *
c *   y(neq+1)- nao definido                                           *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y=Ax              *
c *                                                                    *
c **********************************************************************      
      implicit none
      integer neq,ia(*),ja(*),i,k
      integer neqf1i,neqf2i
c ... ponteiros        
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      real*8  ad(*),a(*),al(*),x(*),y(*),t
      real*8 dum4
c ......................................................................
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         t = ad(i)*x(i)
c
c ...    Produto da linha i pelo vetor x:
c
         do 100 k = ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
c ......................................................................
      return
      end                
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csr_sym_v3(neq,ia  ,ja,a  ,x   ,y,flag)
c **********************************************************************
c * Data de criacao    : 01/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ * 
c * MATVEC_CSR_SYM_V3 : produto matriz-vetor CSR padrao simetrico      *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq       - numero de equacoes                                     *
c * ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro        *
c *                   coeficiente nao-nulo da linha   i                *
c * ja(*)     - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor au                                *
c * a(*)   - coeficientes                                              *
c * x(neq) - vetor a ser multiplicado                                  *
c * y(neq) - nao definido                                              *
c * flag   - .true.  triangular superior                               *
c *        - .false. triangular inferior                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *                                                                    *
c * ------------------------------------------------------------------ *
c * y(neq) - vetor contendo o resultado do produto y = Ax              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * CSR padrao com tem vetores (ia,ja,a)                               *
c * ja em ordem crescente                                              *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),i,k,kk,jak
      real*8  a(*),x(*),y(*),xi,t,s
      logical flag
c ......................................................................
c
c ...
      if(flag) then
        y(1:neq) = 0.0d0
        do i = 1, neq
          kk   = ia(i)
          xi   = x(i)
          y(i) = y(i) + a(kk)*xi
c ...
          do k = kk+1, ia(i+1)-1
            jak   = ja(k)
            s     = a(k)
c ... parte superior
            y(i)   =  y(i)  + s*x(jak)
c ... parte inferior
            y(jak) = y(jak) + s*xi
          enddo
c .....................................................................
        enddo
c .....................................................................
c
c ...
      else
       do i = 1, neq
          xi   = x(i)
c ...
          do k = ia(i), ia(i+1)-2
            jak   = ja(k)
            s     = a(k)
c ... parte superior
            t     =  t      + s*x(jak)
c ... parte inferior
            y(jak) = y(jak) + s*xi
          enddo
c .....................................................................
          y(i) = t + a(k)*xi
        enddo
c .....................................................................
      endif
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csrcsym(neq,ia,ja,dum0,dum1,ad,al,dum2,dum3,x,y,
     .                          neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                          i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCSYM: produto matriz-vetor y = Ax  (A simetrica),      *
c *                   coef. de A no formato CSRC.                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),dum0,dum1,dum2,dum3,i,k,jak
      real*8  ad(*),al(*),x(*),y(*),s,t,xi
      real*8 dum4
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
c
c ... Comunicacao do vetor y no sistema non-overlapping:
c
      if (novlp) call communicate(y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli)
c ......................................................................
      return
      end
      subroutine matvec_csrc(neq,ia,ja,dum0,dum1,ad,al,au,dum2,x,y,
     .                       neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                       i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRC: produto matriz-vetor y = Ax  (A nao-simetrica),     *
c *                coef. de A no formato CSRC e grafo simetrico.       *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor al                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A no formatp CSR           *
c *   au(*)  - parte triangular superior de A no formatp CSC           *
c *            al(k) = Aij,  au(k) = Aji, i < j                        *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),dum0,dum1,dum2,i,k,jak
      real*8  ad(*),al(*),au(*),x(*),y(*),t,xi
      real*8 dum4   
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + al(k)*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + au(k)*xi
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
c
c ... Comunicacao do vetor y no sistema non-overlapping:
c
      if (novlp) call communicate(y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli)
c ......................................................................
      return
      end
      subroutine matvec_csrcrsym(neq,ia,ja,ia1,ja1,ad,al,dum0,ar,x,y,
     .                           neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                           i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCRSYM: produto matriz-vetor y = Ax  (A simetrica),     *
c *                    coef. de A no formato CSRC e grafo simetrico.   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor al                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A no formatp CSR           *
c *   au(*)  - parte triangular superior de A no formatp CSC           *
c *            al(k) = Aij,  au(k) = Aji, i < j                        *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),ia1(*),ja1(*),i,k,jak,dum0
      real*8  ad(*),al(*),ar(*),x(*),y(*),t,xi,s
      real*8 dum4   
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
c
c ... Comunicacao do vetor x no sistema overlapping:
c
      call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
      do 200 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            s = al(k)
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
  100    continue
c
         do 110 k = ia1(i), ia1(i+1)-1
            jak = ja1(k)
c
c ...       Produto da linha i pelo vetor x (retangulo a direita):
c
            t   = t + ar(k)*x(jak)
  110    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  200 continue
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
      return
      end
      subroutine matvec_csrcr(neq,ia,ja,ia1,ja1,ad,al,au,ar,x,y,
     .                        neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                        i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCR: produto matriz-vetor y = Ax  (A nao-simetrica),    *
c *                 coef. de A no formato CSRC e grafo simetrico.      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor al                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A no formatp CSR           *
c *   au(*)  - parte triangular superior de A no formatp CSC           *
c *            al(k) = Aij,  au(k) = Aji, i < j                        *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),ia1(*),ja1(*),i,k,jak
      real*8  ad(*),al(*),au(*),ar(*),x(*),y(*),t,xi
      real*8  dum4  
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
c
c ... Comunicacao do vetor x no sistema overlapping:
      call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
      do 200 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + al(k)*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + au(k)*xi
  100    continue
c 
         do 110 k = ia1(i), ia1(i+1)-1
            jak = ja1(k)
c
c ...       Produto da linha i pelo vetor x (retangulo a direita):
c
            t   = t + ar(k)*x(jak)
  110    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t         
  200 continue
c ......................................................................
      matvectime = matvectime + MPI_Wtime() - time0
      return
      end              
      subroutine matvec_csrc1(neq,ia,ja,dum0,dum1,ad,al,au,dum2,x,y,
     .                        neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                        i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRC1: produto matriz-vetor y = Ax  (A nao-simetrica),    *
c *                 coef. de A no formato CSRC e grafo simetrico.      *
c *                 (versao com loop interno desenrolado)              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor al                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A no formatp CSR           *
c *   au(*)  - parte triangular superior de A no formatp CSC           *
c *            al(k) = Aij,  au(k) = Aji, i < j                        *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),dum0,dum1,dum2,i,k,jak
      real*8  ad(*),al(*),au(*),x(*),y(*),t,xi
      integer k1,k2,jak1
      real*8 dum4
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
c ----------------------------------------------------------------------
      do 300 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
         k1 = ia(i)
         k2 = ia(i+1)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         if (k2 .eq. k1) goto 120
         if (mod(k2-k1,2) .eq. 0) goto 100
         jak1 = ja(k1)
         t = t + al(k1)*x(jak1)
         y(jak1) = y(jak1) + au(k1)*xi
         k1 = k1 + 1
  100    do 110 k = k1, k2-1, 2
            jak  = ja(k)
            jak1 = ja(k+1)
            t = t + al(k)*x(jak) + al(k+1)*x(jak1)
            y(jak)  = y(jak)  +   au(k)*xi
            y(jak1) = y(jak1) + au(k+1)*xi
  110    continue
c
c ...    Armazena o resultado em y(i):
c
  120    y(i) = t
c ----------------------------------------------------------------------
  300 continue
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
c
c ... Comunicacao do vetor y no sistema non-overlapping:
c
      if (novlp) call communicate(y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli)
c ......................................................................
      return
      end
      subroutine matvec_csrcsym1(neq,ia,ja,dum0,dum1,ad,al,au,dum2,x,y,
     .                           neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                           i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCSYM1: produto matriz-vetor y = Ax  (A simetrica),     *
c *                    coef. de A no formato CSRC.                     *
c *                    (versao com loop interno desenrolado)           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),dum0,dum1,dum2,i,k,jak
      real*8  ad(*),al(*),au(*),x(*),y(*),t,xi,s,s1
      integer k1,k2,n,jak1
      real*8 dum4
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
c ----------------------------------------------------------------------
      do 300 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
         k1 = ia(i)
         k2 = ia(i+1)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         if (k2 .eq. k1) goto 120
         n = mod(k2-k1,2)
         if (n .eq. 0) goto 100
         jak1 = ja(k1)
         t = t + al(k1)*x(jak1)
         y(jak1) = y(jak1) + au(k1)*xi
         k1 = k1 + 1
  100    do 110 k = k1, k2-1, 2
            jak  = ja(k)
            jak1 = ja(k+1)
            s  = al(k)
            s1 = al(k+1)
            t = t + s*x(jak) + s1*x(jak1)
            y(jak)  = y(jak)  +  s*xi
            y(jak1) = y(jak1) + s1*xi
  110    continue
c
c ...    Armazena o resultado em y(i):
c  
  120    y(i) = t
c ......................................................................
  300 continue
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
c
c ... Comunicacao do vetor y no sistema non-overlapping:
c
      if (novlp) call communicate(y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli)
c ......................................................................
      return
      end
      subroutine matvec_csrc2(neq,ia,ja,ad,al,au,x,y,neqf1i,neqf2i,
     .                        i_fmapi,i_xfi,i_rcvsi,i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRC2: produto matriz-vetor y = Ax  (A nao-simetrica),    *
c *                 coef. de A no formato CSRC e grafo simetrico       *
c *                 (versao com loops desenrolados).                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor al                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A no formatp CSR           *
c *   au(*)  - parte triangular superior de A no formatp CSC           *
c *            al(k) = Aij,  au(k) = Aji, i < j                        * 
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************            
      implicit none
      integer neq,ia(*),ja(*),i,k,jak
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c .....................................................................      
      real*8  ad(*),al(*),au(*),x(*),y(*),t,xi,t1,xi1
      integer k1,k2,n,jak1,jak2,i0,k3
      real*8 dum4
c ......................................................................
c
c ... Loop nas linhas:
c
      i0 = 1
      n = mod(neq,2)
      if (n .eq. 0) goto 10
      y(1) = ad(1)*x(1)
      i0 = 2
c -------------------------------------------------------------      
   10 do 300 i = i0, neq, 2
c
c ...    Produto da diagonal de A por x:
c
         xi  = x(i)
         xi1 = x(i+1)
         t   =   ad(i)*xi
         t1  = ad(i+1)*xi1                    
         k1 = ia(i)
         k2 = ia(i+1)
         k3 = ia(i+2)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         if (k2 .eq. k1) goto 120
         n = mod(k2-k1,2)
         if (n .eq. 0) goto 100
         jak = ja(k1)         
         t = t + al(k1)*x(jak)
         y(jak) = y(jak) + au(k1)*xi
         k1 = k1 + 1
  100    do 110 k = k1, k2-1, 2
            jak1 = ja(k)
            jak2 = ja(k+1)
            t = t + al(k)*x(jak1) + al(k+1)*x(jak2)
            y(jak1) = y(jak1) + au(k)*xi
            y(jak2) = y(jak2) + au(k+1)*xi                        
  110    continue
c
c ...    Armazena o resultado em y(i):
c  
  120    y(i) = t
c -------------------------------------------------------------
c
c ...    Loop nos coeficientes nao nulos da linha i+1:
c
         if (k3 .eq. k2) goto 220
         n = mod(k3-k2,2)
         if (n .eq. 0) goto 200
         jak = ja(k2)         
         t1 = t1 + al(k2)*x(jak)
         y(jak) = y(jak) + au(k2)*xi1
         k2 = k2 + 1
  200    do 210 k = k2, k3-1, 2
            jak1 = ja(k)
            jak2 = ja(k+1)
            t1 = t1 + al(k)*x(jak1) + al(k+1)*x(jak2)
            y(jak1) = y(jak1) +   au(k)*xi1
            y(jak2) = y(jak2) + au(k+1)*xi1                        
  210    continue  
c
c ...    Armazena o resultado em y(i):
c
  220 y(i+1) = t1
c ----------------------------------------------------------------------  
  300 continue
c ......................................................................
      return
      end
      subroutine matvec_csrcsym2(neq,ia,ja,ad,al,au,x,y,neqf1i,neqf2i,
     .                           i_fmapi,i_xfi,i_rcvsi,i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCSYM2: produto matriz-vetor y = Ax  (A simetrica),     *
c *                    coef. de A no formato CSRC.                     *
c *                    (versao com loops desenrolados)                 *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer neq,ia(*),ja(*),i,k,jak
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      real*8  ad(*),al(*),au(*),x(*),y(*),t,xi,t1,xi1,s,s1
      integer k1,k2,n,jak1,jak2,i0,k3
      real*8 dum4
c ......................................................................
c
c ... Loop nas linhas:
c
      i0 = 1
      n = mod(neq,2)
      if (n .eq. 0) goto 10
      y(1) = ad(1)*x(1)
      i0 = 2
c -------------------------------------------------------------      
   10 do 300 i = i0, neq, 2
c
c ...    Produto da diagonal de A por x:
c
         xi  = x(i)
         xi1 = x(i+1)
         t   =   ad(i)*xi
         t1  = ad(i+1)*xi1                    
         k1 = ia(i)
         k2 = ia(i+1)
         k3 = ia(i+2)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         if (k2 .eq. k1) goto 120
         n = mod(k2-k1,2)
         if (n .eq. 0) goto 100
         jak = ja(k1)         
         t = t + al(k1)*x(jak)
         y(jak) = y(jak) + au(k1)*xi
         k1 = k1 + 1
         
  100    do 110 k = k1, k2-1, 2
            jak1 = ja(k)
            jak2 = ja(k+1)
            s =  al(k)
            s1 = al(k+1)
            t = t + s*x(jak1) + s1*x(jak2)
            y(jak1) = y(jak1) + s*xi
            y(jak2) = y(jak2) + s1*xi                        
  110    continue
c
c ...    Armazena o resultado em y(i):
c  
  120    y(i) = t
c -------------------------------------------------------------
c
c ...    Loop nos coeficientes nao nulos da linha i+1:
c
         if (k3 .eq. k2) goto 220
         n = mod(k3-k2,2)
         if (n .eq. 0) goto 200
         jak = ja(k2)         
         t1 = t1 + al(k2)*x(jak)
         y(jak) = y(jak) + au(k2)*xi1
         k2 = k2 + 1
  200    do 210 k = k2, k3-1, 2
            jak1 = ja(k)
            jak2 = ja(k+1)
            s = al(k)
            s1 = al(k+1)            
            t1 = t1 + s*x(jak1) + s1*x(jak2)
            y(jak1) = y(jak1) +  s*xi1
            y(jak2) = y(jak2) + s1*xi1                        
  210    continue  
c
c ...    Armazena o resultado em y(i):
c
  220 y(i+1) = t1
c ----------------------------------------------------------------------  
  300 continue
c ......................................................................
      return
      end
      subroutine matvec_csrsym(neq,ia,ja,dum0,dum1,ad,a,dum2,dum3,x,y,
     .                         neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                         i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRSYM: produto matriz-vetor y = Ax  (A simetrica),       *
c *   -------------                   coef. de A no formato CSR        *
c *                               e armazenamento da parte inferior.   *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro       *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor a                               *
c *   ad(neq) - diagonal da matriz A                                   *
c *   a(nad)  - coef. de A, sem a diagonal                             *
c *   al(*)   - nao utilizado                                          *
c *   x(neq+1)- vetor a ser multiplicado                               *
c *   y(neq+1)- nao definido                                           *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y=Ax              *
c *                                                                    *
c **********************************************************************      
      implicit none
      integer neq,ia(*),ja(*),i,k,dum0,dum1,dum2,dum3,jak
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      real*8  ad(*),a(*),x(*),y(*),t,xi
      real*8 dum4
c ......................................................................
      do 110 i = 1, neq
         xi = x(i)
c
c ...    Produto da diagonal de A por x:
c
         t = ad(i)*xi
c
c ...    Parte superior de A, linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            t = t + a(k)*x(ja(k))
  100    continue
         y(i) = y(i) + t
c
c ...    Parte inferior de A, coluna i:
c
         do 105 k = ia(i), ia(i+1)-1
            jak = ja(k)
            y(jak) = y(jak) + a(k)*xi
  105    continue
  110 continue
c ......................................................................
      return
      end
      subroutine matvec_csrsym1(neq,ia,ja,dum0,dum1,ad,a,dum2,dum3,x,y,
     .                      neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .                      dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRSYM1: produto matriz-vetor y = Ax  (A simetrica),      *
c *   -------------   coef. de A no formato CSR e armazenamento        *
c *                   da parte superior, com loops aninhados.          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro       *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor a                               *
c *   ad(neq) - diagonal da matriz A                                   *
c *   a(nad)  - coef. de A, sem a diagonal                             *
c *   al(*)   - nao utilizado                                          *
c *   x(neq+1)- vetor a ser multiplicado                               *
c *   y(neq+1)- nao definido                                           *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y=Ax              *
c *                                                                    *
c **********************************************************************      
      implicit none
      integer neq,ia(*),ja(*),i,k,dum0,dum1,dum2,dum3,jak
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c .....................................................................      
      real*8  ad(*),a(*),x(*),y(*),t,s,xi
      real*8  dum4
c ......................................................................
      y(1:neq) = 0.d0
      do 110 i = 1, neq
         xi = x(i)
c
c ...    Produto da diagonal de A por x:
c
         t  = ad(i)*xi
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s = a(k)
c
c ...    Parte superior de A, linha i:
c
            t = t + s*x(jak)
c
c ...    Parte inferior de A, coluna i:
c
            y(jak) = y(jak) + s*xi
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = y(i) + t
  110 continue
c ......................................................................
      return
      end
      subroutine matvec_csrcr1(neq,ia,ja,ia1,ja1,ad,al,au,ar,x,y,
     .                         neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                         i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCR1: produto matriz-vetor y = Ax  (A nao-simetrica),   *
c *                  coef. de A no formato CSRCR e grafo simetrico.    *
c *                  (versao com loop interno desenrolado)             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor al                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A no formatp CSR           *
c *   au(*)  - parte triangular superior de A no formatp CSC           *
c *            al(k) = Aij,  au(k) = Aji, i < j                        *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),ia1(*),ja1(*),i,k,jak,jak1,k1,k2
      real*8  ad(*),al(*),au(*),ar(*),x(*),y(*),t,xi
      real*8 dum4   
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
c
c ... Comunicacao do vetor x no sistema overlapping:
c
      call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
      do 200 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Produto da linha i pelo vetor x (parte quadrada de A):
c
         k1 = ia(i)
         k2 = ia(i+1)
         if (k2 .eq. k1) goto 120
         if (mod(k2-k1,2) .eq. 0) goto 100
         jak1 = ja(k1)         
         t = t + al(k1)*x(jak1)
         y(jak1) = y(jak1) + au(k1)*xi
         k1 = k1 + 1
  100    do 110 k = k1, k2-1, 2
            jak  = ja(k)
            jak1 = ja(k+1)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior de A):
c
            t = t + al(k)*x(jak) + al(k+1)*x(jak1)
c
c ...       Produto dos coef. da parte triangular superior de A por x(i):
c
            y(jak)  = y(jak)  +   au(k)*xi
            y(jak1) = y(jak1) + au(k+1)*xi                        
  110    continue
c
c ...       Produto da linha i pelo vetor x (retangulo a direita):
c
  120    k1 = ia1(i)
         k2 = ia1(i+1)
         if (k2 .eq. k1) goto 150
         if (mod(k2-k1,2) .eq. 0) goto 130
         t = t + ar(k1)*x(ja1(k1))
         k1 = k1 + 1
  130    do 140 k = k1, k2-1, 2
            jak = k+1
            t = t + ar(k)*x(ja1(k)) + ar(jak)*x(ja1(jak))
  140    continue
c
c ...    Armazena o resultado em y(i):
c
  150    y(i) = t
  200 continue
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
      return
      end
      subroutine matvec_csrcrsym1(neq,ia,ja,ia1,ja1,ad,al,au,ar,x,y,
     .                            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCRSYM1: produto matriz-vetor y = Ax  (A simetrica),    *
c *                     coef. de A no formato CSRCR e grafo simetrico. *
c *                     (versao com loop interno desenrolado)          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor al                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A no formatp CSR           *
c *   au(*)  - parte triangular superior de A no formatp CSC           *
c *            al(k) = Aij,  au(k) = Aji, i < j                        *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),ia1(*),ja1(*),i,k,jak,jak1,k1,k2
      real*8  ad(*),al(*),au(*),ar(*),x(*),y(*),t,xi,s,s1
      real*8 dum4
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
c
c ... Comunicacao do vetor x no sistema overlapping:
c
      call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
      do 200 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
         k1 = ia(i)
         k2 = ia(i+1)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         if (k2 .eq. k1) goto 120
         if (mod(k2-k1,2) .eq. 0) goto 100
         jak1 = ja(k1)
         t = t + al(k1)*x(jak1)
         y(jak1) = y(jak1) + au(k1)*xi
         k1 = k1 + 1
  100    do 110 k = k1, k2-1, 2
            jak  = ja(k)
            jak1 = ja(k+1)
            s  = al(k)
            s1 = al(k+1)
            t = t + s*x(jak) + s1*x(jak1)
            y(jak)  = y(jak)  +  s*xi
            y(jak1) = y(jak1) + s1*xi
  110    continue
c
c ...       Produto da linha i pelo vetor x (retangulo a direita):
c
  120    k1 = ia1(i)
         k2 = ia1(i+1)
         if (k2 .eq. k1) goto 150
         if (mod(k2-k1,2) .eq. 0) goto 130
         t = t + ar(k1)*x(ja1(k1))
         k1 = k1 + 1
  130    do 140 k = k1, k2-1, 2
            jak = k+1
            t = t + ar(k)*x(ja1(k)) + ar(jak)*x(ja1(jak))
  140    continue
c
c ...    Armazena o resultado em y(i):
c
  150    y(i) = t
  200 continue
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csrc_pm(neq,nequ,ia,ja,iapu,japu,ad,al,apul,x,y,
     .                        neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                        i_dspli)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRC_PM:produto matriz-vetor y = Ax  (A Kuu, Kpp e Kpu )  *
c *                   coef. de A no formato CSRC.                      *
c *       | Kuu  -Kpu |                                                *
c *   A = |           |                                                *
c *       | Kpu   Kpp |                                                *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   nequ  - numero de equacoes no bloco Kuu                          *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(nad  ) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor al                              *
c *   iapu(neqp+1)-ia(i) informa a posicao no vetor au do primeiro     *
c *                     coeficiente nao-nulo da linha   i              *
c *   japu(nadpu) - ja(k) informa a coluna do coeficiente que ocupa    *
c *               a posicao k no vetor apul                            *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de Kuu e Kpp no formato CSR   *
c *   apul(*)- Kpu no formato CSR                                      *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer neq,nequ,neqp,ia(*),ja(*),iapu(*),japu(*),i,ii,k,jak
      real*8  ad(*),al(*),apul(*),x(*),y(*),s,t,xi
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas: Kuu e Kpp
c
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
c .....................................................................
c
c ... loop nas linha Kpu
      neqp = neq - nequ
      do 120 i = 1, neqp
        ii = nequ+i
        xi = x(ii)
        do 130 k = iapu(i), iapu(i+1)-1
          jak   = japu(k)
          s     = apul(k)
          
c ... Kpu
          y(ii)  =  y(ii)  + s*x(jak)
c ... Kup
          y(jak) =  y(jak) - s*xi
  130   continue
  120 continue
c .....................................................................
c
c ...
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
c
c ... Comunicacao do vetor y no sistema non-overlapping:
c
      if (novlp) call communicate(y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli)
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csrcr_sym_pm(neq    ,dum0  ,ia
     1                          ,ja         ,iar   ,jar
     2                          ,ad         ,al    ,ar   
     3                          ,x          ,y    
     4                          ,neqf1i     ,neqf2i
     5                          ,i_fmapi    ,i_xfi ,i_rcvsi
     6                          ,i_dspli    ,dum4)
c **********************************************************************
c * Data de criacao    : 27/10/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ * 
c * MATVEC_CSRC_csr_sym_pm: produto matriz-vetor y = Ax  (A simetrica),*
c *                    coef. de A no formato CSRC e grafo simetrico.   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq        - numero de equacoes                                    *
c * ia(neq+1)  - ia(i) informa a posicao no vetor au do primeiro       *
c *                   coeficiente nao-nulo da equacao i                *
c * ja(neq+1)  - ja(k) informa a coluna do coeficiente que ocupa       *
c *             a posicao k no vetor al                                *
c * iar(neq+1) - ia(i) informa a posicao no vetor ar do primeiro       *
c *             coeficiente nao-nulo da equacao i da parte retangular  *
c * jar(nadr)  - ja(k) informa a coluna do coeficiente que ocupa       *
c *             a posicao k no vetor ar                                *
c * ad(neq)    - diagonal da matriz A                                  *
c * al(nad)    - parte triangular inferior de A no formato CSR         *
c * ar(nad)    - parte retangular de A no formatp CSR                  *
c * x(neqovlp) - vetor a ser multiplicado                              *
c * y(neq)     - nao definido                                          *
c * neqf1i     - numero de equacoes no buffer de recebimento (MPI)     *
c * neqf2i     - numero de equacoes no buffer de envio (MPI)           *
c * i_fmapi    - ponteiro para o mapa de comunicacao  (MPI)            *
c * i_xfi      - ponteiro para o buffer de valores    (MPI)            *
c * i_rcvsi    - ponteiro extrutura da comunicacao    (MPI)            *
c * i_dspli    - ponteiro extrutura da comunicacao    (MPI)            *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * y(neq) - vetor contendo o resultado do produto y = Ax              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer*8 ia(*),iar(*),k
      integer neq,ja(*),jar(*),i,jak,dum0
      real*8  ad(*),al(*),ar(*),x(*),y(*),t,xi,s
      real*8 dum4   
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
c
c ... Comunicacao do vetor x no sistema overlapping:
c
      call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
      do 200 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            s = al(k)
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
  100    continue
c
         do 110 k = iar(i), iar(i+1)-1
            jak = jar(k)
c
c ...       Produto da linha i pelo vetor x (retangulo a direita):
c
            t   = t + ar(k)*x(jak)
  110    continue
 
c ...    Armazena o resultado em y(i):
 
         y(i) = t
  200 continue
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csrc_sym_pm(neq      ,dum0  ,ia
     1                             ,ja       ,dum1  ,dum2
     2                             ,ad       ,al    ,dum3
     3                             ,x        ,y   
     4                             ,neqf1i   ,neqf2i
     5                             ,i_fmapi  ,i_xfi ,i_rcvsi
     6                             ,i_dspli  ,dum4)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 15/12/2016                                    *
c * ------------------------------------------------------------------ * 
c * MATVEC_CSRC_SYM_PM: produto matriz-vetor y = Ax  (A simetrica),    *
c *                   coef. de A no formato CSRC.                      *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq   - numero de equacoes                                         *
c * ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro        *
c *                   coeficiente nao-nulo da linha   i                *
c * ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor au                                *
c * ad(neq)- diagonal da matriz A                                      *
c * al(nad)- parte triangular inferior de A, no formato CSR, ou        *
c *          parte triangular superior de A, no formato CSC            *
c * x(neq) - vetor a ser multiplicado                                  *
c * y(neq) - nao definido                                              *
c * neqf1i - numero de equacoes no buffer de recebimento (MPI)         *
c * neqf2i - numero de equacoes no buffer de envio (MPI)               *
c * i_fmapi- ponteiro para o mapa de comunicacao  (MPI)                *
c * i_xfi  - ponteiro para o buffer de valores    (MPI)                *
c * i_rcvsi- ponteiro extrutura da comunicacao    (MPI)                *
c * i_dspli- ponteiro extrutura da comunicacao    (MPI)                *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * y(neq) - vetor contendo o resultado do produto y = Ax              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer*8 ia(*),k
      integer neq,ja(*),dum0,dum1,dum2,dum3,i,jak
      real*8  ad(*),al(*),x(*),y(*),s,t,xi
      real*8 dum4
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Loop nas linhas:
c
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
c
c ... Comunicacao do vetor y no sistema non-overlapping:
c
      if (novlp) call communicate(y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli)
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csr_pm(ni,nj,ia  ,ja,al  ,x   ,y,flag)
c **********************************************************************
c * Data de criacao    : 27/03/2016                                    *
c * Data de modificaco : 30/04/2016                                    *
c * ------------------------------------------------------------------ * 
c * MATVEC_CSR_PM: produto matriz-vetor retangular y = Ax              *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neqi      - numero de linhas da matriz A                           *
c * neqj      - numero de colunas da matriz A                          *
c * ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro        *
c *                   coeficiente nao-nulo da linha   i                *
c * ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor au                                *
c * al(nad)- coeficientes                                              *
c * x(neq) - vetor a ser multiplicado                                  *
c * y(neq) - nao definido                                              *
c * flag   - protudo transposto                                        *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *                                                                    *
c * ------------------------------------------------------------------ *
c * y(*)   - vetor contendo o resultado do produto                     * 
c *        .true.  - y(nj) = (At)x                                     *
c *        .false. - y(ni) = Ax                                        *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'time.fi'
      integer ni,nj,ia(*),ja(*),i,k,jak
      real*8  al(*),x(*),y(*),xi
c ... ponteiros      
      logical flag
c ......................................................................
c
c ... loop nas linha Kup
      if(flag) then
        y(1:nj) = 0.d0
        do i = 1, ni
          xi = x(i)
          do k = ia(i), ia(i+1)-1
            jak   = ja(k)
c ... Kup
            y(jak) =  y(jak) - al(k)*xi
          enddo      
        enddo   
c ......................................................................
c
c ... loop nas linha Kpu
      else 
        do i = 1, ni
          y(i) = 0.0d0
c ......................................................................
c
          do k = ia(i), ia(i+1)-1
            jak   = ja(k)
c ... Kpu
            y(i)  =  y(i)  + al(k)*x(jak)
          enddo
        enddo
      endif
c .....................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine saxpb(a,b,x,n,c)
c **********************************************************************
c *                                                                    *
c *   SAXPB: escalar . vetor + vetor                                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    x    - escalar                                                  *
c *    n    - dimensao                                                 *
c *    c(n) - nao definido                                             *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   c - resultado: a.x + b                                           *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),c(*),x
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i) * x + b(i)
  100 continue
      return
      end
      real*8 function dot(a,b,n)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar a.b                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8  a(*),b(*)
c ......................................................................
      dot = 0.d0
      do 100 i = 1, n
         dot = dot + a(i)*b(i)
  100 continue
c ......................................................................    
      return
      end
      real*8 function dot_par(a,b,neq_doti)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar a.b                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'            
      include 'parallel.fi'
      include 'time.fi'
      integer neq_doti
      integer n,i,k
      real*8  a(*),b(*),tmp
c ......................................................................
      time0 = MPI_Wtime()
      tmp = 0.d0
      do 100 i = 1, neq_doti
         tmp = tmp + a(i)*b(i)
  100 continue
      dottime = dottime + MPI_Wtime() - time0
c ......................................................................
      if (nprcs .gt. 1) then
         call MPI_ALLREDUCE(tmp,dot_par,1,MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_WORLD,ierr)
      else
         dot_par = tmp
      endif
c ......................................................................    
      return
      end      
      real*8 function ddot(x,y,n)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar x.y                                         *
c *                                                                    *
c **********************************************************************
      implicit none    
      real*8 x(*),y(*),temp
      integer i,m,mp1,n,inc
c ......................................................................
      ddot = 0.d0
      temp = 0.d0
      inc  = 10
c ......................................................................
      m = mod(n,inc)
      if(m .eq. 0) goto 40
      do 30 i = 1,m
         temp = temp + x(i)*y(i)
   30 continue
      if(n .lt. inc) goto 60
   40 mp1 = m + 1
      do 50 i = mp1, n, inc
         temp = temp + x(i  )*y(  i) + x(i+1)*y(i+1) + x(i+2)*y(i+2) +
     .                 x(i+3)*y(i+3) + x(i+4)*y(i+4) + x(i+5)*y(i+5) +
     .                 x(i+6)*y(i+6) + x(i+7)*y(i+7) + x(i+8)*y(i+8) +
     .                 x(i+9)*y(i+9)
   50 continue
   60 ddot = temp
      return
      end
c **********************************************************************
      subroutine matvec_sym(a,x,y,nl,nc)
c **********************************************************************
c * Data de criacao    : 18/01/2017                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *
c * MATVEC_SYM: produro matriz vetor simetrico                         *  
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * a    - matriz                                                      *
c * x    -                                                             *
c * y    -                                                             *
c * nl   - numero de linhas                                            *
c * nc   - numero de colunas                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * y = ax                                                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Matriz a com armazenamento da parte triangular inferior            *     
c **********************************************************************
      implicit none
      integer nc,nl,i,j
      real*8  a(nl,*),x(*),y(*),xj,aij
c ......................................................................
      y(1:nl) = 0.d0
      do 200 j = 1, nc 
         xj   = x(j) 
         y(j) = y(j) + a(j,j)*xj  
         do 100 i = j+1, nl 
            aij  = a(i,j) 
            y(i) = y(i) + aij*xj 
            y(j) = y(j) + aij*x(i)
  100    continue
  200 continue
c ......................................................................  
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec(a,x,y,nl,nc)
c **********************************************************************
c * Data de criacao    : 18/01/2017                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *
c * MATVEC_SYM: produro matriz vetor simetrico                         *  
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * a    - matriz                                                      *
c * x    -                                                             *
c * y    -                                                             *
c * nl   - numero de linhas                                            *
c * nc   - numero de colunas                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * y = ax                                                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      integer nc,nl,i,j
      real*8  a(nl,*),x(*),y(*),xj,aij
c ......................................................................
      y(1:nl) = 0.d0
      do 200 j = 1, nc 
         xj   = x(j) 
         do 100 i = 1, nl  
           y(i) = y(i) + a(i,j)*xj
  100    continue
  200 continue
c ......................................................................  
      return
      end
      subroutine aequalb(a,b,n)
c **********************************************************************
c *                                                                    *
c *   AEQUALB:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   a = b                                                            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*)
c ......................................................................
      do 100 i = 1, n
         a(i) = b(i)
  100 continue
      return
      end
      subroutine aminusb(a,b,c,n)
c **********************************************************************
c *                                                                    *
c *   AMINUSB:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = a - b                                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i) - b(i)
  100 continue
      return
      end
      subroutine vsum(a,b,n,c)
c **********************************************************************
c *                                                                    *
c *   VSUM:                                                            *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = a + b                                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i) + b(i)
  100 continue
      return
      end      
      subroutine vsmul(a,x,n,c)
c **********************************************************************
c *                                                                    *
c *   VSUM:                                                            *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    x    - escalar                                                  *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = x*a                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),x,c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i)*x
  100 continue
      return
      end 
c **********************************************************************
      subroutine vet(a,b,c)
c **********************************************************************
c *                                                                    *
c *   VET: Produto vetorial axb                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  c(3),a(3),b(3)
c ......................................................................
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
c ......................................................................    
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 23/10/2016                                    *
c * Data de modificaco : 31/10/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_DOT : calcula o numero de operacoes de pontos flutuantes      *    
c * do produto interno                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n        - dimensao dos vetores                                    *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_dot(n)
      implicit none
      integer n
      flop_dot = 2.d0*n
      return
      end
c ********************************************************************** 
c
c **********************************************************************
c * Data de criacao    : 23/10/2016                                    *
c * Data de modificaco : 31/10/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_CSRC: calcula o numero de operacoes de pontos flutuantes      *    
c * da op matvec CSRC                                                  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * nl       - numero de linhas                                        *
c * nad      - numero de termos fora da diagonal                       *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_csrc(nl,nad)
      implicit none
      integer nl
      integer*8 nad
      flop_csrc = nl + 4.d0*nad
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 23/10/2016                                    *
c * Data de modificaco : 29/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_CG : calcula o numero de operacoes de pontos flutuantes       *    
c * do CG                                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacores                                     *
c * nad      - numero de termos fora da diagonal                       *
c * it       - numero de iteracoes                                     *
c * icod     - 1 - CG                                                  *
c *            2 - PCG                                                 *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_cg(neq,nad,it,icod,mpi)
      implicit none
      include 'mpif.h'      
      integer neq,it,icod,ierr
      integer*8 nad
      real*8 flop_csrc,flop_dot,fmatvec,fdot,flops,gflops
      logical mpi
c
      fmatvec = flop_csrc(neq,nad)
      fdot    = flop_dot(neq)
c ... CG
      if(icod .eq. 1) then
        flops = (fmatvec + 2.d0*fdot + 6.d0*neq + 2.d0)*it     
c ... PCG
      elseif(icod .eq. 2) then
        flops = (fmatvec + 2.d0*fdot + 7.d0*neq + 2.d0)*it 
      endif
c .....................................................................
c
c ...
      if(mpi) then
        call MPI_ALLREDUCE(flops,gflops,1,MPI_REAL8
     .                    ,MPI_SUM,MPI_COMM_WORLD,ierr)
        flops = gflops
      endif
c .....................................................................
c
c ...
      flop_cg = flops
c .....................................................................
c
c ...      
      return
      end
c ********************************************************************** 
c
c **********************************************************************
c * Data de criacao    : 31/10/2016                                    *
c * Data de modificaco : 29/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_SQRM: calcula o numero de operacoes de pontos flutuantes      *   
c * do SQRM                                                            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacores                                     *
c * nad      - numero de termos fora da diagonal                       *
c * it       - numero de iteracoes                                     *
c * icod     - 1 - SQRM                                                *
c *            2 - RSQRM                                               *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_sqrm(neq,nad,it,icod,mpi)
      implicit none
      include 'mpif.h'      
      integer neq,it,icod,ierr
      integer*8 nad
      real*8 flop_csrc,flop_dot,fmatvec,fdot,flops,gflops
      logical mpi
c
      fmatvec = flop_csrc(neq,nad)
      fdot    = flop_dot(neq)
c ... SQRM
      if(icod .eq. 1) then
        flops = (fmatvec + 5.d0*fdot + 6.d0*neq + 2.d0)*it     
c ... RSQRM
      elseif(icod .eq. 2) then
        flops = (fmatvec + 5.d0*fdot + 6.d0*neq + 14.d0)*it 
      endif
c .....................................................................
c
c ...
      if(mpi) then
        call MPI_ALLREDUCE(flops,gflops,1,MPI_REAL8
     .                    ,MPI_SUM,MPI_COMM_WORLD,ierr)
        flops = gflops
      endif
c .....................................................................
c
c ...
      flop_sqrm = flops
c .....................................................................
c
c ...      
      return
      end
c ********************************************************************** 
c
c **********************************************************************
c * Data de criacao    : 15/11/2016                                    *
c * Data de modificaco : 29/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_BICGSTAB(2) : calcula o numero de operacoes de pontos         *
c * flutuantes do BICGSTAB(2)                                          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacores                                     *
c * nad      - numero de termos fora da diagonal                       *
c * it       - numero de iteracoes                                     *
c * icod     - 1 - BICGSTAB2                                           *
c *            2 - PBICGSTAB2                                          *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_bicgstab2(neq,nad,it,icod,mpi)
      implicit none
      include 'mpif.h'      
      integer neq,it,icod,ierr
      integer*8 nad
      real*8 flop_csrc,flop_dot,fmatvec,fdot,flops,gflops
      logical mpi
c
      fmatvec = flop_csrc(neq,nad)
      fdot    = flop_dot(neq)
c ... BICGSTAB2
      if(icod .eq. 1) then
        flops = 0.d0    
c ... PBICGSTAB2
      elseif(icod .eq. 2) then
        flops = (4.d0*fmatvec + 10.d0*fdot + 23.d0*neq + 17.d0)*it 
      endif
c .....................................................................
c
c ...
      if(mpi) then
        call MPI_ALLREDUCE(flops,gflops,1,MPI_REAL8
     .                    ,MPI_SUM,MPI_COMM_WORLD,ierr)
        flops = gflops
      endif
c .....................................................................
c
c ...
      flop_bicgstab2 = flops
c .....................................................................
c
c ...      
      return
      end
c ********************************************************************** 
 
