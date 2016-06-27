c **********************************************************************
c *                                                                    *
c *   ASSBLY.F                                            31/08/2005   *
c *   ========                                                         *
c *   Este arquivo contem subrotinas para montagem de matrizes e       *
c *   vetores a partir das matrizes e vetores de elemento:             *
c *                                                                    *
c *   assbly: chama as rotinas de montagem, de acordo com o formato    *
c *           utilizado.                                               *
c *   csr: monta arranjos globais no formato CSR.                      *
c *   skyline: monta arranjos globais no formato SKYLINE.              *
c *   arestas : monta arranjos globais no formato de ARESTAS.          *
c *   ebe: armazena matrizes de elemento.                              *
c *   blockdg: monta a diagonal em blocos.                             *
c *                                                                    *
c **********************************************************************
c **********************************************************************
c *                                                                    *
c *   ASSBLY.F                                            31/08/2005   *
c *   ========                                                         *
c *   Este arquivo contem subrotinas para montagem de matrizes e       *
c *   vetores a partir das matrizes e vetores de elemento:             *
c *                                                                    *
c *   assbly: chama as rotinas de montagem, de acordo com o formato    *
c *           utilizado.                                               *
c *   csr: monta arranjos globais no formato CSR.                      *
c *   skyline: monta arranjos globais no formato SKYLINE.              *
c *   arestas : monta arranjos globais no formato de ARESTAS.          *
c *   ebe: armazena matrizes de elemento.                              *
c *   blockdg: monta a diagonal em blocos.                             *
c *                                                                    *
c **********************************************************************
      subroutine assbly(s ,p ,ld,ia,ja,au,al,ad,b,nst,neq,lhs,rhs,unsym,
     .                  stge)
c **********************************************************************
c *                                                                    *
c *   ASSBLY: montagem da matriz A e do vetor b.                       *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    s(nst,nst) - matriz do elemento nel                             *
c *    p(nst)     - vetor de elemento                                  *
c *    ld(nst)    - numeracao das equacoes do elemento                 *
c *    ia(*)                                                           *
c *    ja(*)                                                           *
c *    au(*)                                                           *
c *    al(*)                                                           *
c *    ad(*)                                                           *
c *    b(*)                                                            *
c *    nst         - numero de graus de liberdade por elemento         *
c *    lhs         - flag para a montagem da matriz A                  *
c *    rhs         - flag para a montagem de b                         *
c *    unsym       - flag para matriz nao-simetrica                    *
c *    stge        - tipo de armazenamento                             *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c **********************************************************************     
      implicit none
      real*8  s(*),p(*),au(*),al(*),ad(*),b(*)
      integer ld(*),ia(*),ja(*),nst,neq,stge
      logical lhs,rhs,unsym
c .....................................................................
c
c ... csr (ia,ja,ad,au,al))     
      if     (stge .eq. 1) then
         call csr(s,p,ld,ia,ja,au,al,ad,b,nst,neq,lhs,rhs,unsym)
c .....................................................................
c
c ...
      elseif (stge .eq. 2) then
          print*,'(ASSBLY) Formato nao implementado !'
          stop
c         call kedg(ix,ip,edge,s,au,nel,nen,nst,ndf,nc,unsym)
c .....................................................................
c
c ...
      elseif (stge .eq. 3) then
          print*,'(ASSBLY) Formato nao implementado !'
          stop      
c         call kel(s,au,nen,ndf,nst,nc,nel,unsym)
c .....................................................................
c
c ...
      elseif (stge .eq. 4) then
          call skyline(s,p,ld,ja,au,al,ad,b,nst,lhs,rhs,unsym)
c .....................................................................
c
c ...
      elseif (stge .eq. 5) then
          call diagonal(s,p,ld,ad,b,nst,lhs,rhs)
c .....................................................................
c
c ... csr (ia,ja,a) 
      elseif (stge .eq. 6) then
       call csrv3(s,p,ld,ia,ja,ad,b,nst,neq,lhs,rhs)
      endif
c .....................................................................
      return
      end
c **********************************************************************
      subroutine assbly_pm(s       ,p          ,ld
     .                    ,ia      ,ja         ,au
     .                    ,al      ,ad         ,b    ,nst
     .                    ,neq     ,nequ       ,neqp
     .                    ,nad     ,nadu       ,nadp ,nadpu 
     .                    ,lhs     ,rhs        ,unsym,stge
     .                    ,dualCsr ,n_blocks_up)
c **********************************************************************
c *                                                                    *
c *   ASSBLY_PM: montagem da matriz A e do vetor b.                    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    s(nst,nst) - matriz do elemento nel                             *
c *    p(nst)     - vetor de elemento                                  *
c *    ld(nst)    - numeracao das equacoes do elemento                 *
c *    ia(*)                                                           *
c *    ja(*)                                                           *
c *    au(*)                                                           *
c *    al(*)                                                           *
c *    ad(*)                                                           *
c *    b(*)                                                            *
c *    nst         - numero de graus de liberdade por elemento         *
c *    neq         - numero de equacoes                                *
c *    nequ        - numero de equacoes em kuu                         *
c *    nad         - numero de posicoes no CSR dos blocos Kuu e Kpp    *
c *    nadu        - numero de posicoes no CSR dos blocos Kuu          *
c *    nadp        - numero de posicoes no CSR dos blocos Kpp          *
c *    lhs         - flag para a montagem da matriz A                  *
c *    rhs         - flag para a montagem de b                         *
c *    unsym       - flag para matriz nao-simetrica                    *
c *    stge        - tipo de armazenamento                             *
c *    n_blocks_up = numero de blocos                                  *
c *                  1 - ( [ Kuu Kpp]  )                               *
c *                  2 - ( [Kuu, Kpp] e [kpu] )                        *
c *                  3 - ( [Kuu], [Kpp] e [kpu])                       *     
c *    dualCsr     - true - armazenamento em blocos Kuu,Kpp e kpu      *
c *                  false- aramzenamento em unico bloco               *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c **********************************************************************     
      implicit none
      real*8  s(*),p(*),au(*),al(*),ad(*),b(*)
      integer ld(*),ia(*),ja(*),nst
      integer neq,nequ,neqp
      integer nad,nadu,nadp,nadpu
      integer stge,n_blocks_up
      logical lhs,rhs,unsym,dualCsr
c ...................................................................... 
c
c ... csrc     
      if (stge .eq. 1) then
c ... bloco Kuu, Kpp e Kpu
         if(dualCsr) then
c ... | Kuu   0  |
c     |  0   Kpp |         
           if(n_blocks_up .eq. 1) then    
c .....................................................................
c
c ... | Kuu   0  |
c     |  0   Kpp | e [kpu]
           else if(n_blocks_up .eq. 2) then    
             call csr2(s ,p        ,ld
     .                ,ia,ia(neq+2)
     .                ,ja,ja(nad+1)
     .                ,al,al(nad+1)
     .                ,ad,b
     .                ,nst,neq     ,nequ
     .                ,lhs,rhs)
c .....................................................................
c
c ... [Kuu], [Kpp] e [kpu]
           else if(n_blocks_up .eq. 3) then
c            print*,neq,nequ,neqp,nad,nadu,nadp
             call csr3(s   ,p         ,ld
     .                ,ia  ,ia(nequ+2),ia(neq+3)
     .                ,ja  ,ja(nadu+1),ja(nadu+nadp+1)
     .                ,al  ,al(nadu+1),al(nadu+nadp+1)
     .                ,ad  ,ad(nequ+1),b
     .                ,nst ,neq ,nequ 
     .                ,lhs ,rhs)
           endif 
c .....................................................................
c
c ... | Kuu Kpu |
c     | kpu Kpp |
         else
           call csr(s,p,ld,ia,ja,au,al,ad,b,nst,neq,lhs,rhs,unsym)
         endif     
c .....................................................................
c
c ...
      elseif (stge .eq. 2) then
          print*,'(ASSBLY) Formato nao implementado !'
          stop
c         call kedg(ix,ip,edge,s,au,nel,nen,nst,ndf,nc,unsym)
c .....................................................................
c
c ...
      elseif (stge .eq. 3) then
          print*,'(ASSBLY) Formato nao implementado !'
          stop      
c         call kel(s,au,nen,ndf,nst,nc,nel,unsym)
c .....................................................................
c
c ...
      elseif (stge .eq. 4) then
          call skyline(s,p,ld,ja,au,al,ad,b,nst,lhs,rhs,unsym)
c .....................................................................
c
c ...
      elseif (stge .eq. 5) then
          call diagonal(s,p,ld,ad,b,nst,lhs,rhs)
c ...................................................................... 
c
c ... csr (ia,ja,a) 
      elseif (stge .eq. 6) then
       call csrv3(s,p,ld,ia,ja,ad,b,nst,neq,lhs,rhs)
      endif
c ......................................................................       
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csr(s,p,ld,ia,ja,au,al,ad,b,nst,neq,lhs,rhs,unsym)
c **********************************************************************
c *                                                                    *
c *   CSR : armazenamento da matriz no formato CSR.                    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    s(nst,nst) - matriz do elemento nel                             *
c *    p(nst)     - vetor de elemento                                  *
c *    ld(nst)    - numeracao das equacoes do elemento                 *
c *    ia(neq+1)  - ia(i) informa a posicao no vetor a do primeiro     *
c *                 coeficiente nao-nulo da equacao i                  *
c *    ia(neq+2:2*neq+2) - apontador dos coeficientes da parte         *
c *                        retangular da matriz                        *
c *    ja(nad)    - ja(k) informa a coluna do coeficiente que ocupa    *
c *                 a posicao k no vetor a                             *
c *    ja(nad+1:nad+nad1) - colunas da parte retangular da matriz      *
c *                         em particoes overlapping                   *
c *    au(nad)      - parte triangular superior da matriz global       *
c *                 armazenada no formato CSR                          *
c *    al(nad)      - parte triangular inferior da matriz global       *
c *                 armazenada no formato CSR                          *
c *    al(nad+1:nad+nad1) - parte retangular da matriz para particoes  *
c *                         overlapping                                *
c *    ad(*)      - diagonal da matriz de coeficientes                 *
c *    b(neq)     - vetor de forcas                                    *
c *    nst        - numero de graus de liberdade por elemento          *
c *    lhs        - flag para a montagem da matriz A                   *
c *    rhs        - flag para a montagem de b                          *
c *    unsym      - flag para matriz nao-simetrica                     *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    au   -  coeficientes acima da diagonal                          *
c *    al   -  coeficientes abaixo da diagonal                         *
c *    ad   -  coeficientes da diagonal                                *
c *    b    -  vetor de forcas corrigido                               *
c *                                                                    *
c **********************************************************************
      implicit none
      integer ld(*),ia(*),ja(*),nst,neq
      integer i,j,k,l,kk,nad,ovlp
      real*8  s(nst,*),p(*),au(*),al(*),ad(*),b(*)
      logical lhs,rhs,unsym
c ......................................................................
c
c ... Verifica se a particao eh overlapping (versao em paralelo):
      nad  = ia(neq+1)-1
c     ovlp = ja(nad+1)
c ...................................................................... 
      do 200 i = 1, nst
         k = ld(i)
         if (k .le. 0 .or. k .gt. neq) goto 200
         if(rhs) b(k) = b(k) - p(i) 
c ......................................................................
         if(lhs) then
            ad(k) = ad(k) + s(i,i)         
            do 110 j = 1, nst
               l = ld(j)
               if (l .le. 0) goto 110
               do 100 kk = ia(k), ia(k+1)-1
                  if (l .eq. ja(kk)) then
                     al(kk) = al(kk) + s(i,j)
                     if(unsym) au(kk) = au(kk) + s(j,i)
                  endif
  100          continue
c ......................................................................
c
c ...     Monta a parte retangular da matriz para particoes overlapping:
c              if (ovlp .gt. 0) then 
c                 do 101 kk = ia(neq+k+1), ia(neq+k+2)-1
c                    if (l.eq.ja(nad+kk)) al(nad+kk) = al(nad+kk)+s(i,j)
c 101             continue
c ......................................................................         
c              endif
  110       continue
         endif
  200 continue
c ...................................................................... 
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csr2(s  ,p   ,ld
     .               ,ia ,iapu
     .               ,ja ,japu
     .               ,al ,apul
     .               ,ad ,b
     .               ,nst,neq ,nequ
     .               ,lhs,rhs)
c **********************************************************************
c * Data de criacao    : 12/12/2015                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * CSR2: armazenamento da matriz no formato CSRC para os blocos uu    *
c * e pp e csr para o bloco up                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * s(nst,nst) - matriz do elemento nel                                *
c * p(nst)   - vetor de elemento                                       *
c * ld(nst)  - numeracao das equacoes do elemento                      *
c * ia(neq+1)- ia(i) informa a posicao no vetor a do primeiro          *
c *            coeficiente nao-nulo da equacao i do bloco Kuu e Kpp    *
c * iaup(neq+1)- iaup(i) informa a posicao no vetor a do primeiro      *
c *            coeficiente nao-nulo da equacao i do bloco Kpu          *      
c * ja(nad)    - ja(k) informa a coluna do coeficiente que ocupa       *
c *              a posicao k no vetor a do bloco Kuu e Kpp             *
c * jaup(nadup)- ja(k) informa a coluna do coeficiente que ocupa       *
c *              a posicao k no vetor a do bloco Kpu                   *
c * al(nad)      - parte triangular inferior da matriz global          *
c *              armazenada no formato CSR dp bloco Kuu e Kpp          *
c * apul(nadup)  - parte triangular inferior da matriz global          *
c *              armazenada no formato CSR dp bloco Kpu                *                                
c * ad(*)      - diagonal da matriz de coeficientes                    *
c * b(*)       - vetor de forcas                                       *
c * nst        - numero de graus de liberdade por elemento             *
c * neq        - numero de equacoes totais                             *
c * nequ       - numero de equacoes de do bloco Kuu                    *
c * lhs        - flag para a montagem da matriz A                      *
c * rhs        - flag para a montagem de b                             *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * al   -  coeficientes abaixo da diagonal dos blocos Kuu e Kpp       *
c * apul -  coeficientes abaixo da diagonal dos blocos Kuu e Kpp       *
c * ad   -  coeficientes da diagonal dos blocos Kuu e Kpp              *
c * b    -  vetor de forcas corrigido                                  *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Kuu e Kpp juntos e Kpu separado                                    *
c **********************************************************************
      implicit none
      integer ld(*),ia(*),ja(*),iapu(*),japu(*),nst,neq,nequ
      integer i,j,k,ki,l,kk
      real*8  s(nst,*),p(*),al(*),ad(*),b(*),apul(*)
      logical lhs,rhs
c ......................................................................
c
c ...................................................................... 
      do 200 i = 1, nst
         k = ld(i)
         if (k .le. 0 .or. k .gt. neq) goto 200
         if(rhs) b(k) = b(k) - p(i)
c ......................................................................
         if(lhs) then
            ad(k) = ad(k) + s(i,i)
c ... bloco Kuu
            if (k .le. nequ) then 
               do 110 j = 1, nst
                  l = ld(j)
                  if( l .le. 0 ) go to 110
                  if (l .le. nequ)  then
                     do 100 kk = ia(k), ia(k+1)-1
                        if (l .eq. ja(kk)) then
                          al(kk) = al(kk) + s(i,j)
                        endif
  100                continue
                  endif   
  110          continue
c .....................................................................
c
c ... bloco Kpp e kpu
            else if (k .gt. nequ) then 
               do 140 j = 1, nst
                  l = ld(j)
                  if( l .le. 0 ) go to 140
c ... bloco Kpp                  
                  if (l .gt. nequ )  then
                     do 130 kk = ia(k), ia(k+1)-1
                        if (l .eq. ja(kk)) then
                          al(kk) = al(kk) + s(i,j)
                        endif
  130                continue
c ..................................................................
c                 
c ... bloco Kpu    
                  else if(l .le. nequ )  then
                     ki = k - nequ     
                     do 150 kk = iapu(ki), iapu(ki+1)-1
                        if (l .eq. japu(kk)) then
                          apul(kk) = apul(kk) + s(i,j)
                        endif
  150                continue
                  endif
  140          continue
            endif     
c ......................................................................
         endif
  200 continue
c ...................................................................... 
      return
      end      
c **********************************************************************
c
c **********************************************************************
      subroutine csr3(s   ,p    ,ld
     .               ,iau ,iap  ,iapu
     .               ,jau ,jap  ,japu
     .               ,alu ,alp  ,alpu
     .               ,adu ,adp  ,b
     .               ,nst ,neq  ,nequ
     .               ,lhs ,rhs)
c **********************************************************************
c * Data de criacao    : 12/12/2015                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * CSR3: armazenamento da matriz no formato CSRC para os blocos uu    *
c * e pp e csr para o bloco up                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * s(nst,nst) - matriz do elemento nel                                *
c * p(nst)   - vetor de elemento                                       *
c * ld(nst)  - numeracao das equacoes do elemento                      *
c * iau(nequ+1)- ia(i) informa a posicao no vetor a do primeiro        *
c *            coeficiente nao-nulo da equacao i do bloco Kuu          *
c * iap(neqp+1)- ia(i) informa a posicao no vetor a do primeiro        *
c *            coeficiente nao-nulo da equacao i do bloco Kpp          *
c * iaup(neq+1)- iaup(i) informa a posicao no vetor a do primeiro      *
c *            coeficiente nao-nulo da equacao i do bloco Kpu          *      
c * jau(nadu)  - ja(k) informa a coluna do coeficiente que ocupa       *
c *              a posicao k no vetor a do bloco Kuu                   *
c * jau(nadp)  - ja(k) informa a coluna do coeficiente que ocupa       *
c *              a posicao k no vetor a do bloco Kpp                   *
c * jaup(nadup)- ja(k) informa a coluna do coeficiente que ocupa       *
c *              a posicao k no vetor a do bloco Kpu                   *
c * al(nad)      - parte triangular inferior da matriz global          *
c *              armazenada no formato CSR dp bloco Kuu e Kpp          *
c * apul(nadup)  - parte triangular inferior da matriz global          *
c *              armazenada no formato CSR dp bloco Kpu                *                                
c * adu(*)     - diagonal da matriz de coeficientes Kuu                *
c * adp(*)     - diagonal da matriz de coeficientes Kpp                *
c * b(*)       - vetor de forcas                                       *
c * nst        - numero de graus de liberdade por elemento             *
c * neq        - numero de equacoes totais                             *
c * nequ       - numero de equacoes de do bloco Kuu                    *
c * lhs        - flag para a montagem da matriz A                      *
c * rhs        - flag para a montagem de b                             *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * al   -  coeficientes abaixo da diagonal dos blocos Kuu e Kpp       *
c * apul -  coeficientes abaixo da diagonal dos blocos Kuu e Kpp       *
c * ad   -  coeficientes da diagonal dos blocos Kuu e Kpp              *
c * b    -  vetor de forcas corrigido                                  *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Kuu ,Kpp e Kpu separados                                           *
c **********************************************************************
      implicit none
      integer ld(*)
      integer iau(*) ,jau(*) ,iap(*),jap(*),iapu(*),japu(*)
      integer nst,neq,nequ
      integer i,j,k,ki,l,kk
      real*8  s(nst,*),p(*)
      real*8  alu(*),alp(*),alpu(*),adu(*),adp(*),b(*)
      logical lhs,rhs
c ......................................................................
c
c ...................................................................... 
      do 200 i = 1, nst
         k = ld(i)
         if (k .le. 0 .or. k .gt. neq) goto 200
         if(rhs) b(k) = b(k) - p(i)
c ......................................................................
         if(lhs) then
c ... bloco Kuu
            if (k .le. nequ) then 
               adu(k) = adu(k) + s(i,i)
               do 110 j = 1, nst
                  l = ld(j)
                  if( l .le. 0 ) go to 110
                  if (l .le. nequ)  then
                     do 100 kk = iau(k), iau(k+1)-1
                        if (l .eq. jau(kk)) then
                          alu(kk) = alu(kk) + s(i,j)
                        endif
  100                continue
                  endif   
  110          continue
c .....................................................................
c
c ... bloco Kpp e kpu
            else if (k .gt. nequ) then 
               ki = k - nequ     
               adp(ki) = adp(ki) + s(i,i)
               do 140 j = 1, nst
                  l = ld(j)
                  if( l .le. 0 ) go to 140
c ... bloco Kpp                  
                  if (l .gt. nequ )  then
                     do 130 kk = iap(ki), iap(ki+1)-1
                        if ( (l - nequ) .eq. jap(kk)) then
                          alp(kk) = alp(kk) + s(i,j)
                        endif
  130                continue
c ..................................................................
c                 
c ... bloco Kpu    
                  else if(l .le. nequ )  then
                     ki = k - nequ     
                     do 150 kk = iapu(ki), iapu(ki+1)-1
                        if (l .eq. japu(kk)) then
                          alpu(kk) = alpu(kk) + s(i,j)
                        endif
  150                continue
                  endif
  140          continue
            endif     
c ......................................................................
         endif
  200 continue
c ...................................................................... 
      return
      end      
c **********************************************************************
c
c **********************************************************************
      subroutine csrv3(s,p,ld,ia,ja,a,b,nst,neq,lhs,rhs)
c **********************************************************************
c * Data de criacao    : 30/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *      
c * CSRV3 : armazenamento da matriz no formato CSR na versao           *
c * de 3 arranjos padrao (ia,ja,a)                                     *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                  *
c * s(nst,nst) - matriz do elemento nel                                *
c * p(nst)     - vetor de elemento                                     *
c * ld(nst)    - numeracao das equacoes do elemento                    *
c * ia(neq+1)  - ia(i) informa a posicao no vetor a do primeiro        *
c *              coeficiente nao-nulo da equacao i                     *
c * ja(nad)    - ja(k) informa a coluna do coeficiente que ocupa       *
c *              a posicao k no vetor a                                *
c * a(nad))    -  matriz global armazenada no formato CSR              *
c * b(neq)     - vetor de forcas                                       *
c * nst        - numero de graus de liberdade por elemento             *
c * lhs        - flag para a montagem da matriz A                      *
c * rhs        - flag para a montagem de b                             *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * a(nad))    -  matriz global armazenada no formato CSR              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * CSR padrao com tres arranjos, ia(neq+1) ,ja(nnz) e a(nnz) onde nnz *
c * e o numero total de nao zeros                                      *  
c **********************************************************************
      implicit none
      integer ld(*),ia(*),ja(*),nst,neq
      integer i,j,k,l,kk
      real*8  s(nst,*),p(*),a(*),b(*)
      logical lhs,rhs
c ...................................................................... 
      do 200 i = 1, nst
         k = ld(i)
         if (k .le. 0 ) goto 200
         if(rhs) b(k) = b(k) - p(i) 
c ......................................................................
         if(lhs) then
           do 110 j = 1, nst
             l = ld(j)
             if (l .le. 0) goto 110
             do 100 kk = ia(k), ia(k+1)-1
               if (l .eq. ja(kk)) then
                 a(kk) = a(kk) + s(i,j)
               endif
  100        continue
c ......................................................................
  110     continue
         endif
  200 continue
c ...................................................................... 
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine skyline(s,p,ld,jd,au,al,ad,b,nst,lhs,rhs,unsym)
c **********************************************************************
c *                                                                    *
c *   SKYLINE: monta arranjos globais no formato SKYLINE.              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    s(nst,nst) - matriz de elemento                                 *
c *    p(nst)     - vetor de elemento                                  *
c *    ld(nst)    - numeracao local das equacoes                       *
c *    jd(neq)    - ponteiro da diagonal                               *
c *    au(nad) - coeficientes acima da diagonal                        *
c *    al(nad) - coeficientes abaixo da diagonal                       *
c *    ad(neq) - diagonal da matriz de coeficientes                    *
c *    b(neq)  - vetor de forcas                                       *
c *    nst   - numero de graus de liberdade por elemento               *
c *    lhs   - flag para a montagem da matriz A                        *
c *    rhs   - flag para a montagem de b                               *
c *    unsym - flag para matriz nao-simetrica                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ad   -  coeficientes da diagonal                                *
c *    au   -  coeficientes acima da diagonal                          *
c *    al   -  coeficientes abaixo da diagonal                         *
c *    b    -  vetor de forcas corrigido                               *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nst,ld(nst),jd(*),i,j,ii,jj,jc
      real*8  ad(*),au(*),al(*),b(*),s(nst,nst),p(nst)
      logical lhs,rhs,unsym
c ......................................................................
c
c.... loop nas linhas da matriz de elemento:
c
      do 200 i = 1, nst
         ii = ld(i)
         if (ii .gt. 0) then
            if (rhs) b(ii) = b(ii) - p(i)
            if (lhs) then
               ad(ii) = ad(ii) + s(i,i)          
c ...................................................................
c
c ............ loop nas colunas:
c
               do 100 j = 1, nst
                  if (ld(j) .gt. ii) then
                     jc = ld(j)
                     jj = ii + jd(jc) - jc + 1
                     au(jj) = au(jj) + s(i,j)
                     if(unsym) al(jj) = al(jj) + s(j,i)
                  endif
100            continue
            endif
         endif
200   continue
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine arestas(ix,ip,edge,s,se,nel,nen,nst,ndf,nc,unsym)
c **********************************************************************
c *                                                                    *
c *   ARESTAS : armazenamento da matriz global em arestas.             *
c *             Diagonal em blocos armazenada separadamente.           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    ip(nnode+1)     - apontador de arestas                          *
c *    edge(2,nedge)   - conetividades das arestas                     *
c *    s(nst,nst)      - matriz do elemento nel                        *
c *    se(nc,nedge)    - matriz global armazenada no formato CSR       *
c *    nel   - numero do elemento                                      *
c *    nen   - numero de nos por elemento                              *
c *    nst   - numero de graus de liberdade por elemento               *
c *    ndf   - numero de graus de liberdade por no                     *
c *    nc    - numero de coeficientes armazenados por aresta           *
c *    unsym - flag para matrizes nao simetricas                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    se(nc,nedge)    - matriz global armazenada no formato ARESTAS   *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,ix(nen+1,*),ip(*),edge(2,*),nel,nst,ndf,nc
      integer i,j,k,no,no1,no2,noi,noj,n,m,ll,kk,l
      real*8  s(nst,nst),se(nc,*)
      logical unsym      
c.......................................................................
      do 400 i = 1, nen-1
         do 300 j = i+1, nen 
            noi = i
            noj = j
            no1 = ix(i,nel)
            no2 = ix(j,nel)
            if (no2 .lt. no1) then
               noi = j
               noj = i
               no1 = ix(j,nel)
               no2 = ix(i,nel)
            endif
            do 200 n = ip(no1), ip(no1+1)-1
               no = edge(2,n)
               if (no2 .eq. no) then
                  m = 1
                  do 110 ll = 1, ndf
                     l = (noi-1)*ndf+ll 
                     do 100 kk = 1, ndf
                        k = (noj-1)*ndf+kk
                        se(m,n) = se(m,n) + s(l,k)
                        m = m + 1
  100                continue
  110             continue
                  if (unsym) then
                     do 130 ll = 1, ndf
                        l = (noj-1)*ndf+ll 
                        do 120 kk = 1, ndf
                           k = (noi-1)*ndf+kk
                           se(m,n) = se(m,n) + s(l,k)
                           m = m + 1
  120                   continue
  130                continue
                  endif
               endif
  200       continue
  300    continue
  400 continue
c ......................................................................  
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine ebe(s,se,nen,ndf,nst,nc,nel,unsym)
c **********************************************************************
c *                                                                    *
c *   EBE : armazenamento das matrizes de elemento (termos fora da     *
c *                                             diagonal em blocos).   *
c *                                                                    *
c *    s(nst,nst) - matriz do elemento nel                             *
c *    se(nc,numel)- matrizes de elemento                              *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero de graus de liberdade por no                     *
c *    nst   - graus de liberdade por elemento                         *
c *    nc    - numero de coeficientes armazenados por elemento         *
c *    nel   - elemento                                                *
c *    unsym - flag para matrizes nao simetricas                       *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,ndf,nst,nc,nel,i,j,k,l,m,ii,jj
      real*8  s(nst,*),se(nc,*)
      logical unsym
c ......................................................................
      k = 0
      if(unsym) then
         do 130 i  = 1, nen
         do 120 ii = 1,ndf
            l = (i-1)*ndf+ii
            do 110 j = 1, nen
               if(j .eq. i) goto 110
               do 100 jj = 1, ndf
                  m = (j-1)*ndf+jj
                  k = k + 1
                        se(k,nel) = s(l,m) 
  100          continue
  110       continue
  120    continue
  130    continue
      else
         do 230  i = 1, nen-1
         do 220 ii = 1, ndf
            l = (i-1)*ndf+ii
            do 210  j = i+1, nen
            do 200 jj = 1, ndf
               m = (j-1)*ndf+jj
               k = k + 1
               se(k,nel) = s(l,m)
  200       continue
  210       continue
  220    continue
  230    continue
      endif
c ......................................................................      
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine blockdg(ix,s,bd,nen,nst,ndf,nc,nel,unsym)
c **********************************************************************
c *                                                                    *
c *   BLOCKDG - diagonal em blocos                                     *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    s(nst,nst)      - matriz do elemento nel                        *
c *    bd(nc,nnode)    - diagonal em blocos                            *
c *    nen   - numero de nos por elemento                              *
c *    nst   - numero de graus de liberdade por elemento               *
c *    ndf   - numero de graus de liberdade por no                     *
c *    nc    - numero de coeficientes por no                           *
c *    nel   - elemento                                                *
c *    unsym - flag para matrizes nao simetricas                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    bd - diagonal em blocos                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,nst,nc,nel,ndf,ix(nen+1,*)
      integer n,i,k,l,kk,ll,no
      real*8  s(nst,nst),bd(nc,*)
      logical unsym
c ======================================================================
c
c ... Diagonal em blocos:
c
      if(unsym) then
c
c ...    matriz nao-simetrica:
c
         do 110 i = 1, nen
            no = ix(i,nel)
            k  = (i-1)*ndf+1
            l  = k
            n  = 1
            do 100 kk = k, k+ndf-1
            do 100 ll = l, l+ndf-1
               bd(n,no) = bd(n,no) + s(kk,ll)
               n = n + 1
  100       continue
  110    continue
      else
c
c ...    matriz simetrica:
c
         do 210 i = 1, nen
            no = ix(i,nel)
            k  = (i-1)*ndf+1
            l  = k
            n  = 1
            do 200 kk = k, k+ndf-1
            do 200 ll = kk, l+ndf-1
               bd(n,no) = bd(n,no) + s(kk,ll)
               n = n + 1
  200       continue
  210    continue
      endif
c ......................................................................      
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine diagonal(s,p,ld,ad,b,nst,lhs,rhs)
c **********************************************************************
c *                                                                    *
c *   SKYLINE: monta arranjos globais no formato SKYLINE.              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    s(nst,nst) - matriz de elemento                                 *
c *    p(nst)     - vetor de elemento                                  *
c *    ld(nst)    - numeracao local das equacoes                       *
c *    jd(neq)    - ponteiro da diagonal                               *
c *    au(nad) - coeficientes acima da diagonal                        *
c *    al(nad) - coeficientes abaixo da diagonal                       *
c *    ad(neq) - diagonal da matriz de coeficientes                    *
c *    b(neq)  - vetor de forcas                                       *
c *    nst   - numero de graus de liberdade por elemento               *
c *    lhs   - flag para a montagem da matriz A                        *
c *    rhs   - flag para a montagem de b                               *
c *    unsym - flag para matriz nao-simetrica                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ad   -  coeficientes da diagonal                                *
c *    au   -  coeficientes acima da diagonal                          *
c *    al   -  coeficientes abaixo da diagonal                         *
c *    b    -  vetor de forcas corrigido                               *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nst,ld(nst),i,j,ii,jj,jc
      real*8  ad(*),b(*),s(nst,nst),p(nst)
      logical lhs,rhs
c ......................................................................
c
c.... loop nas linhas da matriz de elemento:
c
      do 200 i = 1, nst
         ii = ld(i)
         if (ii .gt. 0) then
            if (rhs) b(ii) = b(ii) - p(i)
            if (lhs) then
               ad(ii) = ad(ii) + s(i,i)          
            endif
         endif
200   continue
c ......................................................................
      return
      end
c **********************************************************************
