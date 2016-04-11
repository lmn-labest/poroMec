c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c *   PRE_DIAG: precondicionador diagonal                              *
c *   --------                                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * m   - indefinido                                                   *
c * ad  - coeficientes da diagonal principal                           *
c * neq - numero de equacoes                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * m   - precondicionador diagonal ( M-1)                             *
c *                                                                    *
c **********************************************************************
      subroutine pre_diag(m,ad,neq)
      implicit none
      real*8 m(*),ad(*)
      integer i,neq
      do i = 1, neq
        m(i) = 1.d0/ad(i)
      enddo
      return
      end
c **********************************************************************
c
c **********************************************************************
       subroutine ildlt(n,ia,ja,al,ad,ldlt,v,shift,fcheck)  
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * ILDLT: fatora imcompleta LDLt com matriz simetrica A no formato    *
c * CSRC.                                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * n - numero de equacoes                                             *
c * ia(1:n+1) - ponteiro das linhas, para o formato CSR                *
c * ja(1:nad) - ponteiro das colunas no formato CSR                    *
c * al(1:nad) - parte triangular inferior de A                         *
c * ad(1:n)   - diagonal da matriz A                                   *
c * ldlt(1:n+nad) - nao definido                                       *
c * v(1:n)    - arranjo auxiliar (nao inicializado)                    *
c * shift     - parametro de deslocamento para manter a fotaracao      * 
c *             estavel ( shift >= 0)                                  *   
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * ldlt(n+1:n+nad) - fator L                                          *
c * ldlt(1:n)       - diagonal fatorada                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * LDLT: fatora a matriz simetrica A na forma LDL. Os coeficientes    *
c * Aij (j<i) sao armazenados em a(1:nad), no formato CSR,             *
c * e a diagonal eh armazenada em ad(1:n). Os fatores L e D sao        *
c * armazenados em ldlt(n+1:n+1+nad) e ldlt(1:n), no mesmo formato CSR *
c * e com a mesma estrututa ia, ja. Pela propria definicao do CSR,     *
c * esta fatoracao eh incompleta. No entanto, pode tornar-se           *
c * completa utilizando-se a tecnica de fill-in.                       *
c * Esta implementacao eh uma adaptacao para o formato CSR do          *
c * algoritmo 5.1.2 em Matriz Computations, p. 84, Golub e Van Loan.   *
c **********************************************************************
      implicit none      
      integer n,ia(*),ja(*),i,j,k,kk,i1,i2,kj,nad
      real*8 al(*),ad(*),ldlt(*),v(*),t
      real*8 dmax,mod,shift
      logical fcheck  
      real*8 tol
      parameter (tol=1.d-14)
c ......................................................................
      nad      = ia(n+1)-1
      do i = 1, nad
        ldlt(i+n) = al(i)
      enddo
c ... A = A + shift*diag(A)
      ldlt(1:n)   = ad(1:n) + shift*ad(1:n)
c ......................................................................
      v(1:n)   = 0.d0
c ......................................................................
c
c ...
      dmax = 0.d0
      do 50 i = 1, n
        mod  = dabs(ad(i))
        dmax = max(mod,dmax)
  50  continue
c ......................................................................
c
c ......................................................................      
      do 1000 j = 1, n   
        i1 = ia(j)
        i2 = ia(j+1)-1
        do 100 i = i1, i2
          k    = ja(i)
          v(k) = ldlt(i+n)*ldlt(k)
  100   continue
c ......................................................................                    
        t = 0.d0
        do 200 i = i1, i2
          t = t + ldlt(i+n)*v(ja(i))
  200   continue
        v(j) = ldlt(j) - t
c ... pivo negativo ou muito pequeno
        if(v(j) .lt. tol*dmax .and. fcheck) then
          print*,'Valor do Pivo invalido !!!',v(j)
          stop
        endif
        ldlt(j) = v(j)
c .....................................................................                  
        do 310 k = j+1, n
          kj = 0
          t  = 0.d0
          do 300 i = ia(k), ia(k+1)-1
            kk = ja(i)
            if (kk .lt. j) then
              t = t + ldlt(i+n)*v(kk)
            else if (kk .eq. j) then
              kj = i + n
            endif
  300     continue
          if (kj .gt. 0 ) then 
            ldlt(kj) = (ldlt(kj) - t)/v(j)
          endif
  310   continue
c .....................................................................
        v(j) = 0.d0  
        do 400 i = i1, i2             
          k    = ja(i) 
          v(k) = 0.d0   
  400   continue
c .....................................................................  
 1000 continue
c .....................................................................
c
c ...
      ldlt(1:n)   = 1.d0/ldlt(1:n)
c .....................................................................
c     open(15,file='ildlt.csr')
c     do i = 1, 1500
c       write(15,'(i9,1d15.2)')i,ldlt(i)
c     enddo 
c     write(15,*)'l'
c     do i = 1, nad
c       write(15,'(i9,2d15.2)')i,al(i),ldlt(i+n)
c     enddo 
c      close(15)
      return
      end
c **********************************************************************
c
c **********************************************************************
c **********************************************************************
      subroutine ildlt_solv(n,ia,ja,jat,iat,kat,ad,a,b,x)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *  
c * ildlt_solv: resolve o sistema com a fatoracao imcompleta LDLt      *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n - numero de equacoes                                             *
c * ia(1:n+1)  - ponteiro das linhas, para o formato CSR               *
c * ja(1:nad)  - ponteiro das colunas no formato CSR                   *
c * jat(1:n+1) - ponteiro das colunas, para o formato CSC              *
c * iat(1:nad) - ponteiro das linhas no formato CSC                    *
c * kat(1:nad  - ponteiro que relaciona o CSC com a matriz a do CSR    * 
c * ad(1:nad)  - fatores D da matriz A ( D: 1/D )                      *
c * a(1:nad)   - fatores L da matriz A                                 *
c * b(1:n)     - vetor independente                                    *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * x(1:n) - vetor solucao.                                            *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Todas as informação da matriz são da parte inferior                * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'time.fi'
      integer n,ia(*),ja(*),jat(*),iat(*),kat(*)
      integer i,k,kk
      real*8  a(*),ad(*),b(*),x(*),tmp
c ......................................................................
c
c ...
      ifatsolvtime = Mpi_Wtime() - ifatsolvtime
      x(1:n) = b(1:n)
c ......................................................................
c
c ... Forward substitution:
       do 110 i = 2, n
        tmp = x(i)
        do 100 k = ia(i), ia(i+1) - 1
           tmp = tmp - a(k)*x(ja(k))
  100   continue
        x(i) = tmp 
  110 continue
c ......................................................................
c
c ...
      do 115 i = 1, n
        x(i) = x(i) * ad(i)
  115 continue
c ......................................................................
c
c ... Backward substitution:
      do 210 i = n-1, 1, -1
        tmp = x(i)
        do 200 k = jat(i), jat(i+1)-1
           kk   = kat(k)
           tmp  = tmp - a(kk)*x(iat(k))
  200   continue
        x(i) = tmp 
  210 continue
c ......................................................................
      ifatsolvtime = Mpi_Wtime() - ifatsolvtime    
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine ildlt_csrc_aux(neq,ia,ja,jat,iat,at)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *  
c * ILDLT_CSRC_AUX: resolve o sistema com a fatoracao imcompleta LDLt  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq        - numero de equacoes                                    *
c * ia(1:n+1)  - ponteiro das linhas, para o formato CSR               *
c * ja(1:nad)  - ponteiro das colunas no formato CSR                   *
c * jat(1:n+1) - nao definido                                          *
c * iat(1:nad) - nao definido                                          *
c * kat(1:nad  - nao definido                                          * 
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * jat(1:n+1) - ponteiro das colunas, para o formato CSC              *
c * iat(1:nad) - ponteiro das linhas no formato CSC                    *
c * kat(1:nad  - ponteiro que relaciona o CSC com a matriz a do CSR    * 
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Todas as informação da matriz são da parte inferior                * 
c **********************************************************************
 
      implicit none
      integer ia(*),ja(*),jat(*),iat(*),at(*)
      integer neq,nl,nc,nnc,k
c ...
      nnc          = 0
      jat(1:neq+1) = 0
c ... loop na colunas
      do nc = 1, neq
c ... loop nas linhas abaixo da diagonal
        do nl = nc + 1, neq
c ...
          do k = ia(nl), ia(nl+1) - 1
            if( nc .eq. ja(k)) then
              jat(nc)  = jat(nc) + 1 
              nnc      = nnc     + 1
              iat(nnc) = nl
              at(nnc)  = k
            endif
          enddo
c ....................................................................
        enddo
c ....................................................................
      enddo
c .....................................................................
c
      do k = neq, 1, -1
        jat(k+1) =  jat(k)    
      enddo
c     
      jat(1) = 1
      do k = 1, neq 
        jat(k+1) =  jat(k+1) + jat(k)   
      enddo
c .....................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine ldltc(n,a,v,w)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * LDLTC : fatoracao LDLt completa com a matriz cheia                 *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n      - numero de equacoes                                        *
c * a(n,n) - matriz cheia                                              *
c * v(n)   - vetor auxiliar de trabalho                                *
c * w(n)   - vetor auxiliar de trabalho                                *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * a(n,n) - matriz fatorada LDLt                                      *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * parte triangular superio inalterada                                *  
c **********************************************************************      
      implicit none
      integer n,i,j,k
      real*8 a(n,*),v(*),w(*),dot
c ......................................................................
      do j = 1, n
         print*,'neq',j
         do i = 1, j-1
            v(i) = a(j,i)*a(i,i)
         enddo
         w(1:j-1) = a(j,1:j-1)
         v(j) = a(j,j) - dot(w,v,j-1)
         if (v(j) .eq. 0.d0) then
           print*, 'Subrotina LDLTC: Pivot nulo !'
           stop
         endif 
         a(j,j) = v(j)
         do k = j+1, n
            w(1:j-1) = a(k,1:j-1)
            a(k,j) = (a(k,j) - dot(w,v,j-1))/v(j) 
         enddo
      enddo
c ......................................................................
c     open(15,file='ldtl.matrix')
c     do i = 1, n
c       write(15,'(100d15.1)')(a(i,j),j=1,i)
c     enddo 
c     do i = 1, n
c        write(15,'(100d15.2)')(a(i,j),j=1,i)
c        write(15,'(100d15.2)')a(i,i)
c     enddo 
c     close(15)
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine ildlt_full_matrix(n,a,v,w)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * ILDLTC_FULL_MATRIX : fatoracao LDLt incompleta com a matriz cheia  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n      - numero de equacoes                                        *
c * a(n,n) - matriz cheia                                              *
c * v(n)   - vetor auxiliar de trabalho                                *
c * w(n)   - vetor auxiliar de trabalho                                *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * a(n,n) - matriz fatorada LDLt incompleta                           *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * parte triangular superio inalterada                                *  
c **********************************************************************      
      implicit none
      integer n,i,j,k
      real*8 a(n,*),v(*),w(*),dot
c ......................................................................
      do j = 1, n   
c       print*,'neq',j
        do i = 1, j-1
          v(i) = a(j,i)*a(i,i)
        enddo
        w(1:j-1) = a(j,1:j-1)
        v(j) = a(j,j) - dot(w,v,j-1)
c
        if (v(j) .eq. 0.d0) then
          print*, 'Subrotina LDLTC: Pivot nulo !'
          stop
        endif 
c
        a(j,j) = v(j)
c
        do k = j+1, n
          w(1:j-1) = a(k,1:j-1)
          if( a(k,j) .ne. 0.0d0) then
            a(k,j) = (a(k,j) - dot(w,v,j-1))/v(j)
          endif 
        enddo
      enddo
c ......................................................................
c     open(15,file='ildtl.matrix')
c     do i = 1, 1500
c       write(15,'(i9,d15.2)')i,a(i,i)
c     enddo 
c     write(15,*)'l'
c     k = 1
c     do i = 1, n
c       do j = 1, i - 1
c         if( a(i,j) .ne. 0.0d0 ) then 
c           write(15,'(i9,d15.2)')k,a(i,j)
c           k = k + 1
c         endif
c       enddo
c     enddo 
      close(15)
c ...
      return
      end
c **********************************************************************