c **********************************************************************
c *                                                                    *
c *   CSR.F                                               31/08/2005   *
c *                                                                    *
c *   Este arquivo contem subrotinas para montagem da estrutura de     *
c *   dados CSR:                                                       *
c *                                                                    *
c *   csrstruct                                                        *
c *   csrstruct_pm                                                     *
c *   csria                                                            *
c *   csrja                                                            *
c *   csriaup                                                          *
c *   csrjaup                                                          *
c *   csriaup2                                                         *
c *   csrjaup2                                                         *
c *                                                                    *
c **********************************************************************
      subroutine csrstruct_pm(id ,ix  ,num ,nnode,nnodev
     1                     ,numel,nen  ,ndf
     2                     ,neq  ,nequ,neqp
     3                     ,i2   ,i3   
     4                     ,nad  ,naduu,nadpp,nadpu,nadr
     5                     ,lower,diag,upper,right,ija,ja
     6                     ,n_blocks_up,block_pu ,block_pu_sym)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 15/12/2016                                    * 
c * ------------------------------------------------------------------ * 
c * CSRSTRUCT_PM: monta os arranjos ia e ja do formato CSR.            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * id(ndf,nnode)   - numeracao global das equacoes                    *
c * ix(nen+1,numel) - conetividades nodais                             *
c * num(nnode)      - numeracao RCM                                    *
c * nnode - numero de nos total da particao                            *
c * nnodev- numero de nos vertices da particao                         *
c * numel - numero de elementos                                        *
c * nen   - numero de nos por elemento                                 *
c * ndf   - numero max. de graus de liberdade por no                   *
c * neq   - numero de equacoes                                         *
c * nequ  - numero de equacoes de deslocamentos                        *
c * neqp  - numero de equacoes de  pressao                             *
c * i2    - nao definido                                               *
c * i3    - nao definido                                               *
c * nad   - nao definido                                               *
c * naduu - nao definido                                               *
c * nadpp - nao definido                                               *
c * nadpu - nao definido                                               *
c * nadr  - nao definido                                               *
c * lower = .true.  -> inclui a parte triangular inferior no csr       *
c * diag  = .true.  -> inclui a diagonal no csr                        *
c * upper = .true.  -> inclui a parte triangular superior no csr       *
c * ija   = string do nome do vetor ija                                *
c * ja    = string do nome do vetor ja                                 *
c * n_blocks_up = numero de blocos                                     *
c *               1 - ( [ Kuu Kpp]  )                                  *
c *               2 - ( [Kuu, Kpp] e [kpu] )                           *
c *               3 - ( [Kuu], [Kpp] e [kpu])                          *     
c * block_pu     - flag para matrizes Kuu,kpp e Kup nao simetricos     *
c * block_pu_sym - flag para matrizes Kuu, Kpp e Kup simetricos        *
c * ------------------------------------------------------------------ * 
c *   Parametros de saida:                                             *
c * ------------------------------------------------------------------ * 
c * i2    - ponteiro para o arranjo ia(neq+1)                          *
c * i3    - ponteiro para o arranjo ja(nad)                            *
c * block_pu = true                                                    *
c *    nad   - numero de coeficientes nao nulos dos blocos uu e pp     *
c *    naduu - numero de coeficientes nao nulos do bloco uu            *
c *    nadpp - numero de coeficientes nao nulos do bloco pp            *
c * block_pu_sym = true | false; block_pu = false                      *
c * nad   - numero de termos nao nulos fora diagonal principal         *
c * naduu - numero de termos nao nulos fora diagonal principal bloco   *
c *         kuu                                                        *
c * nadpp - numero de termos nao nulos fora diagonal principal bloco   *
c *         kpp                                                        *
c * nadpu - numero de termos nao nulos fora diagonal principal bloco   *
c *         kuu                                                        *
c * nadr  - numero de termos nao nulos na parte retangular             *
c *         ( MPI em overllaping )                                     *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * block_pu = true e n_blocks_up = 1                                  *
c * monta o bloco kuu e kpp separados e nao monta o bloco kup          *
c *                                                                    *
c * block_pu = true e n_blocks_up = 2                                  *
c * monta o bloco kuu e kpp juntos e o bloco kup separado              *
c *                                                                    *
c * block_pu = true e n_blocks_up = 3                                  *
c * monta o bloco kuu, kpp e kup separados                             *
c *                                                                    *
c * block_pu_sym = true                                                *
c * monta o bloco kuu, kpp e kup juntos ( primeiras equacoes u e       *
c * depois as esquacoes de pressao                                     *
c *                                                                    *
c * block_pu_sym = false e block_pu = false                            *                   
c * monta o bloco kuu, kpp e kup juntos sem considera a estrutura      * 
c * blocada, i.e., graus de liberda juntos                             *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer id(ndf,*),ix(nen+1,*),num(*)
      integer nnode,nnodev,numel,nen,ndf,neq,nequ,neqp
      integer nadr,naduu,nadpp,nadpu,n_blocks_up
      logical block_pu,block_pu_sym
      integer*8 nad,nl,nc
c ... ponteiros      
      integer*8 i0,i1,i2,i3
c .....................................................................      
      integer n
      logical lower,diag,upper,right
      character*8 ija,ja
c ......................................................................
c
c ... Grafo da malha:
c
c     ia(i0) => ip(nnode+1) - ip(i) indica a posicao em ips do primeiro
c                             no vizinho ao no i.
c     ia(i1) => ips(ipos) - contem as conectividades nodais de cada no
      call graph(ix,nnode,numel,nen,i0,i1)
c ......................................................................
c      
c ... csrc(uu+pp) + csr(up) nao simetricos
      if(block_pu) then
c ... [Kuu] e [Kpp] apenas. Bloco Kpu a nivel de elemento
        if( n_blocks_up .eq. 1) then
          print*,'CSRSTRUCT: Nao implementado n_blocks 1 !!' 
          stop
c ...................................................................... 
c      
c ... | Kuu   0  |
c     |  0   Kpp | e [kpu]
        else if(n_blocks_up .eq. 2) then
c
c ... Montagem do arranjo ia(neq + 1 + neqp +1):
c
c ... ia(i2)=>ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro
c                               coeficiente nao-nulo da equacao i dos 
c                               blocos Kuu e Kpp
c ... ia(i2+neq+1)=>ia(neq+1) - ia(neq+1+i) informa a posicao no vetor a 
c                               do primeiro coeficiente nao-nulo da 
c                                equacao i do bloco Kpu 
          n = (neq + 1) + (neqp + 1) 
          i2 = alloc_4(ija,1,n)
c ... blocos [ Kuu, Kpp]  e [Kup]
          call csriaup(id,num,ia(i0),ia(i1),ia(i2),ia(i2+neq+1) 
     .                ,nnode,ndf,neq 
     .                ,nequ,neqp,nad,nadpu,nadr,lower,diag,upper 
     .                ,right)
c ......................................................................     
c
c ... Montagem do arranjo ja(nad):
c
c ... ia(i3)=>ja(nad) - ja(k) informa a coluna do coeficiente que ocupa
c                       a posicao k no vetor a do bloco Kuu e Kpp
c
c ... ia(i3+nad)=>ja(nadpu) - ja(k) informa a coluna do coeficiente que 
c                             ocupa a posicao k no vetor a do bloco Kpu
c
          i3 = alloc_4(ja,1,nad+nadpu)
          call csrjaup(id,num,ia(i0),ia(i1),ia(i3),ia(i3+nad)
     .                ,nnode,nnodev,ndf,neq
     .                ,nequ,nad,nadpu,lower,diag,upper,right)
          call sortgraph(ia(i2),ia(i3),neq)
          call sortgraph(ia(i2+neq+1),ia(i3+nad),neqp)
c ......................................................................
c
c ...
          call get_nads(ia(i2),naduu,nadpp,neq,nequ)
c ......................................................................      
c
c ... blocos [Kuu] [Kpp] e [Kpu]
        else if(n_blocks_up .eq. 3) then
c
c ... Montagem do arranjo ia(nequ + 1 + neqp +1 + neqp + 1):
c
c ... ia(i2)=>ia(nequ+1)        - ia(i) informa a posicao no vetor a 
c                                 do primeiro coeficiente nao-nulo da 
c                                 equacao i do bloco Kuu
c ... ia(i2+nequ+1)=>ia(neqp+1) - ia(i) informa a posicao no vetor a 
c                                 do primeiro coeficiente nao-nulo da 
c                                 equacao i dos bloco Kpp
c ... ia(i2+neq+2)=>ia(neqp+1)  - informa a posicao no vetor a 
c                                 do primeiro coeficiente nao-nulo da 
c                                 equacao i do bloco Kpu 
          n = (nequ + 1) + 2*(neqp + 1)
          i2 = alloc_4(ija,1,n)
c ... blocos uu, pp e up
          call csriaup2(id    ,num          ,ia(i0)       ,ia(i1)
     .                 ,ia(i2),ia(i2+nequ+1),ia(i2+neq+2)
     .                 ,nnode ,ndf          
     .                 ,neq   ,nequ  ,neqp         
     .                 ,nad   ,naduu ,nadpp ,nadpu 
     .                 ,lower ,diag  ,upper )
c .......................................................................     
c
c ... Montagem do arranjo ja(nadu+nadp+nadpu):
c
c ... ia(i3)=>ja(nadu)      - ja(k) informa a coluna do coeficiente que ocupa
c                          a posicao k no vetor a do bloco Kuu
c
c ... ia(i3+nadu)=>ja(nadp) - ja(k) informa a coluna do coeficiente que ocupa
c                           a posicao k no vetor a do bloco Kpp
c
c ... ia(i3+nadu+nadp)=>ja(nadpu) - ja(k) informa a coluna do coeficiente que 
c                                   ocupa a posicao k no vetor a do bloco Kpu
c
          i3  = alloc_4(ja,1,naduu+nadpp+nadpu)
          call csrjaup2(id,num,ia(i0),ia(i1)
     .                 ,ia(i3),ia(i3+naduu  ),ia(i3+naduu+nadpp)
     .                 ,nnode ,ndf
     .                 ,neq   ,nequ       
     .                 ,lower ,diag          ,upper,right)
          call sortgraph(ia(i2)       ,ia(i3)            ,nequ)
          call sortgraph(ia(i2+nequ+1),ia(i3+naduu)      ,neqp)
          call sortgraph(ia(i2+neq +2),ia(i3+naduu+nadpp),neqp)
c ......................................................................
      endif
c ......................................................................
c 
c ... csrc(uu+pp+up) simetrico
      else if(block_pu_sym) then
        nadpu = 0
        naduu = 0
        nadpp = 0
c
c ... Montagem do arranjo ia(neq+1):
c
c ... ia(i2)=>ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro
c                               coeficiente nao-nulo da equacao i   
c
        n = neq + 1
        if (right) n = 2*n
        i2 = alloc_4(ija,1,n)
        call csria(id,num,ia(i0),ia(i1),ia(i2),nnode,ndf,neq,nad,nadr,
     .             lower,diag,upper,right)
c ......................................................................     
c
c ... Montagem do arranjo ja(nad):
c
c ... ia(i3)=>ja(nad) - ja(k) informa a coluna do coeficiente que ocupa
c                       a posicao k no vetor a  
c
        i3 = alloc_4(ja,1,nad+nadr)
        call csrja2(id,num,ia(i0),ia(i1),ia(i3),nnode,ndf
     .             ,neq,nequ,nad,lower,diag,upper,right)
        call sortgraph(ia(i2),ia(i3),neq)
        if (right) then
          call sortgraph(ia(i2+neq+1),ia(i3+nad),neq)      
        endif
c ......................................................................
c 
c ... csrc(uu+pp+up) ou csrc para problema desacoplados  
c     | Kuu Kpu |
c     | kpu Kpp |
      else
        nadpu = 0
        naduu = 0
        nadpp = 0
c
c ... Montagem do arranjo ia(neq+1):
c
c ... ia(i2)=>ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro
c                               coeficiente nao-nulo da equacao i   
c
        n = neq + 1
        if (right) n = 2*n
        i2 = alloc_8(ija,1,n)
        call csria(id,num,ia(i0),ia(i1),ia(i2),nnode,ndf,neq,nad,nadr,
     .             lower,diag,upper,right)
c ......................................................................     
c
c ... Montagem do arranjo ja(nad):
c
c ... ia(i3)=>ja(nad) - ja(k) informa a coluna do coeficiente que ocupa
c                       a posicao k no vetor a  
c
        nl = 1
        nc = nad+nadr
        i3 = dalloc_4(ja,nl,nc)
        call csrja(id,num,ia(i0),ia(i1),ia(i3),nnode,ndf,neq,nad,lower,
     .             diag,upper,right)
        call sort_graph_csr(ia(i2),ia(i3),neq,nad,right) 
c ......................................................................
c
c ...      
c       call printIaJa(ia(i2),ia(i3),ia(i2+neq+1),ia(i3+nad),neq,neqp,
c    .                 nad,nadup)
c ......................................................................      
      endif
c ......................................................................
c
c ...
      i1 = dealloc('iaux1   ')
      i0 = dealloc('iaux0   ')
      i2 = locate (ija)
      i3 = locate (ja)
c ......................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csrstruct(id,ix,num,nnode,numel,nen,ndf,neq,i2,i3,nad,
     .                     nadr,lower,diag,upper,right,ija,ja)
c **********************************************************************
c *                                                                    *
c *   CSRSTRUCT: monta os arranjos ia e ja do formato CSR.             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    id(ndf,nnode)   - numeracao global das equacoes                 *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    num(nnode)      - numeracao RCM                                 *
c *    nnode - numero de nos total da particao                         *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    neq   - numero de equacoes                                      *
c *    lower = .true.  -> inclui a parte triangular inferior no csr    *
c *    diag  = .true.  -> inclui a diagonal no csr                     *
c *    upper = .true.  -> inclui a parte triangular superior no csr    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    i2    - ponteiro para o arranjo ia(neq+1)                       *
c *    i3    - ponteiro para o arranjo ja(nad)                         *
c *    nad   - numero de coeficientes nao nulos                        *
c *    nadr  - numero de coeficientes nao nulos da parte retangular    *                                              *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer id(ndf,*),ix(nen+1,*),num(*)
      integer nnode,numel,nen,ndf,neq,nad,nadr
c ... ponteiros      
      integer*8 i0,i1,i2,i3
c .....................................................................      
      integer n
      logical lower,diag,upper,right
      character*8 ija,ja
c ......................................................................
c
c ... Grafo da malha:
c
c     ia(i0) => ip(nnode+1) - ip(i) indica a posicao em ips do primeiro
c                             no vizinho ao no i.
c     ia(i1) => ips(ipos) - contem as conectividades nodais de cada no
      call graph(ix,nnode,numel,nen,i0,i1)
c ......................................................................      
c
c ... Montagem do arranjo ia(neq+1):
c
c ... ia(i2)=>ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro
c                               coeficiente nao-nulo da equacao i   
c
      n = neq + 1
      if (right) n = 2*n
      i2 = alloc_4(ija,1,n)
      call csria(id,num,ia(i0),ia(i1),ia(i2),nnode,ndf,neq,nad,nadr,
     .           lower,diag,upper,right)
c ......................................................................     
c
c ... Montagem do arranjo ja(nad):
c
c ... ia(i3)=>ja(nad) - ja(k) informa a coluna do coeficiente que ocupa
c                       a posicao k no vetor a  
c
      i3 = alloc_4(ja,1,nad+nadr)
      call csrja(id,num,ia(i0),ia(i1),ia(i3),nnode,ndf,neq,nad,lower,
     .           diag,upper,right)
      call sortgraph(ia(i2),ia(i3),neq)
      if (right) then
         call sortgraph(ia(i2+neq+1),ia(i3+nad),neq)      
      endif
c ......................................................................
      i1 = dealloc('iaux1   ')
      i0 = dealloc('iaux0   ')
      i2 = locate (ija)
      i3 = locate (ja)
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csria(id,num,ip,ips,ia,nnode,ndf,neq,nad,nadr,lower,
     .                 diag,upper,right)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 15/12/2016                                    * 
c * ------------------------------------------------------------------ * 
c * CSRIA: monta o arranjo ia do formato CSR                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * id(ndf,nnode)- numeracao global das equacoes                       *
c * num(nnode)   - renumeracao nodal                                   *
c * ip(nnode+1)  - ip(i) indica a posicao em ips do primeiro no        *
c *                   conectado ao no i                                *
c * ips(ipos)    - contem as conetividades nodais de cada no           *
c * nnode - numero de nos                                              *
c * ndf   - numero max. de graus de liberdade por no                   *
c * neq   - numero de equacoes                                         *
c * lower = .true.  -> inclui a parte triangular inferior no csr       *
c * diag  = .true.  -> inclui a diagonal no csr                        *
c * upper = .true.  -> inclui a parte triangular superior no csr       *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro         *
c *                   coeficiente nao-nulo da equacao i                *
c * nad   - numero de coeficientes nao nulos                           *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      integer nnode,ndf,neq,nadr
      integer*8 nad,ia(*)
      integer id(ndf,nnode),ip(nnode+1),ips(*),num(nnode)
      integer i,j,k,ii,jj,kk,neqi,neqj,no
      logical lower,diag,upper,right
c ......................................................................
c
c ... Inicializa o arranjo ia:
c
      do 50 i = 1, neq
         ia(i) = 0
         if (right) ia(i+neq+1) = 0
   50 continue
c ----------------------------------------------
c
c ... Loop nos vertices:
c
      do 140 no = 1, nnode
         i = num(no)
c
c ...    Loop nas equacoes do vertice i:
c
         do 130 ii = 1, ndf
            neqi = id(ii,i)
            if (neqi .gt. 0 .and. neqi .le. neq) then
c ----------------------------------------------
c
c ...          Loop nas equacoes do vertice i:
c
               do 100 kk = 1, ndf
                  neqj = id(kk,i)
                  if (neqj .le. 0) goto 100
                  if (neqj .lt. neqi) then
                     if (lower) ia(neqi) = ia(neqi) + 1
                  elseif (neqj .eq. neqi) then
                     if (diag) ia(neqi)  = ia(neqi) + 1
                  elseif (neqj .gt. neqi) then
                     if (upper) ia(neqi) = ia(neqi) + 1
                  endif
  100          continue
c
c ...          Loop nos vertices conectados ao vertice i:
c
               do 120 k = ip(i), ip(i+1)-1
                  j = ips(k)
c
c ...             Loop nas equacoes do vertice j:
c
                  do 110 jj = 1, ndf
                     neqj = id(jj,j)
                     if (neqj .le. 0) goto 110
                     if (neqj .lt. neqi) then
                        if (lower) ia(neqi) = ia(neqi) + 1
c ... parte retangular                        
                     elseif (neqj .gt. neq) then
                        if (right) ia(neq+1+neqi) = ia(neq+1+neqi) + 1
c ...                        
                     elseif (neqj .gt. neqi) then
                        if (upper) ia(neqi) = ia(neqi) + 1
                     endif
  110             continue
  120          continue
c ----------------------------------------------
            endif
  130    continue
  140 continue
c ----------------------------------------------
      do 200 i = neq, 1, -1
         ia(i+1) = ia(i)
         if (right) ia(i+neq+2) = ia(i+neq+1)         
  200 continue
      ia(1) = 1
      if (right) ia(neq+2) = 1
      do 210 i = 1, neq
         ia(i+1) = ia(i+1) + ia(i)
         if (right) ia(i+neq+2) = ia(i+neq+2) + ia(i+neq+1)         
  210 continue
      nad  = ia(neq+1)-1
      nadr = 0
      if (right) nadr = ia(2*neq+2)-1      
c ......................................................................      
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csrja(id,num,ip,ips,ja,nnode,ndf,neq,nad,lower,diag,
     .                 upper,right)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 15/12/2016                                    * 
c * ------------------------------------------------------------------ * 
c * CSRJA: monta o arranjo ja do formato CSR (matrizes simetricas      *
c *                                             e nao simetricas)      *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * id(ndf,nnode)- numeracao global das equacoes                       *
c * num(nnode)   - renumeracao nodal                                   *
c * ip(nnode+1)  - ip(i) indica a posicao em ips do primeiro           *
c * ips(ipos)    - contem as conetividades nodais de cada no           *
c * nnode - numero de nos                                              *
c * ndf   - numero max. de graus de liberdade por no                   *
c * neq   - numero de equacoes                                         *
c * nad   - numero de coeficientes nao nulos                           *
c * lower = .true.  -> inclui a parte triangular inferior no csr       *
c * diag  = .true.  -> inclui a diagonal no csr                        *
c * upper = .true.  -> inclui a parte triangular superior no csr       *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * ja(nad) - ja(k) informa a coluna do coeficiente que ocupa          *
c *           a posicao k no vetor a                                   *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      integer*8 nad
      integer nnode,ndf,neq
      integer id(ndf,*),num(*),ip(*),ips(*),ja(*)
      integer i,j,k,ii,jj,kk,neqi,neqj,no,n,m
      logical lower,diag,upper,right
c ......................................................................
      m = 0
      n = 0
c
c ... Loop nos vertices:
c
      do 140 no = 1, nnode
         i = num(no)
c
c ...    Loop nas equacoes do vertice i:
c
         do 130 ii = 1, ndf
            neqi = id(ii,i)
            if (neqi .gt. 0 .and. neqi .le. neq) then
c ----------------------------------------------------------
c
c ...          Loop nas equacoes do vertice i:
c
               do 100 kk = 1, ndf
                  neqj  = id(kk,i)
                  if (neqj .le. 0) goto 100
                  if (neqj .lt. neqi) then
                     if (lower) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  elseif (neqj .eq. neqi) then
                     if (diag) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  elseif (neqj .gt. neqi) then
                     if (upper) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  endif                  
  100          continue
c
c ...          Loop nos vertices conectados ao vertice i:
c
               do 120 k = ip(i), ip(i+1)-1
                  j = ips(k)
c
c ...             Loop nas equacoes do vertice j:
c
                  do 110 jj = 1, ndf
                     neqj  = id(jj,j)
                     if (neqj .le. 0) goto 110
                     if (neqj .lt. neqi) then
                        if (lower) then
                           n = n + 1
                           ja(n) = neqj
                        endif
c ... parte retangular
                     elseif (neqj .gt. neq) then
                        if (right) then
                           m = m + 1
                           ja(nad+m) = neqj
                        endif
                     elseif (neqj .gt. neqi) then
                        if (upper) then
                           n = n + 1
                           ja(n) = neqj
                        endif
                     endif
  110             continue
  120          continue
c ----------------------------------------------------------
            endif
  130    continue
  140 continue
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csrja2(id,num,ip,ips,ja,nnode,ndf
     .                 ,neq,nequ,nad,lower,diag,upper,right)
c **********************************************************************
c * Data de criacao    : 26/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * CSRJA2: monta o arranjo ja do formato CSR para os blocos           *
c * (uu, pp e up) todos na mesmo CSRC                                  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * id(ndf,nnode)- numeracao global das equacoes                       *
c * num(nnode)   - renumeracao nodal                                   *
c * ip(nnode+1)  - ip(i) indica a posicao em ips do primeiro           *
c * ips(ipos)    - contem as conetividades nodais de cada no           *
c * nnode - numero de nos                                              *
c * ndf   - numero max. de graus de liberdade por no                   *
c * neq   - numero de equacoes                                         *
c * nequ  - numero de equacoes de deslocamentos                        *
c * nad   - numero de coeficientes nao nulos                           *
c * lower = .true.  -> inclui a parte triangular inferior no csr       *
c * diag  = .true.  -> inclui a diagonal no csr                        *
c * upper = .true.  -> inclui a parte triangular superior no csr       *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ja(nad) - ja(k) informa a coluna do coeficiente que ocupa          *
c *           a posicao k no vetor a                                   *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Kuu e Kpp e kpu juntos                                             *
c **********************************************************************
      implicit none
      integer nnode,ndf,neq,nequ,nad
      integer id(ndf,*),num(*),ip(*),ips(*),ja(*)
      integer i,j,k,ii,jj,kk,neqi,neqj,no,n,m
      logical lower,diag,upper,right
c ......................................................................
      m = 0
      n = 0
c
c ... Loop nos vertices:
c
      do 140 no = 1, nnode
         i = num(no)
c
c ...    Loop nas equacoes do vertice i:
c
         do 130 ii = 1, ndf - 1
            neqi = id(ii,i)
            if (neqi .gt. 0 .and. neqi .le. neq) then
c
c ...          Loop nas equacoes do vertice i:
c
               do 100 kk = 1, ndf - 1
                  neqj  = id(kk,i)
                  if (neqj .le. 0) goto 100
                  if (neqj .lt. neqi) then
                     if (lower) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  elseif (neqj .eq. neqi) then
                     if (diag) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  elseif (neqj .gt. neqi) then
                     if (upper) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  endif                  
  100          continue
c
c ...          Loop nos vertices conectados ao vertice i:
c
               do 120 k = ip(i), ip(i+1)-1
                  j = ips(k)
c
c ...             Loop nas equacoes do vertice j:
c
                  do 110 jj = 1, ndf - 1
                     neqj  = id(jj,j)
                     if (neqj .le. 0) goto 110
                     if (neqj .lt. neqi) then
                        if (lower) then
                           n = n + 1
                           ja(n) = neqj
                        endif
                     elseif (neqj .gt. neq) then
                        if (right) then
                           m = m + 1
                           ja(nad+m) = neqj
                        endif
                     elseif (neqj .gt. neqi) then
                        if (upper) then
                           n = n + 1
                           ja(n) = neqj
                        endif
                     endif
  110             continue
  120          continue
c ......................................................................
            endif
  130    continue
  140 continue
c ......................................................................
c
c ... Equacoes pp
c
c ... Loop nos nos:
c
      do 190 no = 1, nnode
         i = num(no)
c
c ...    equacao da pressao do no i:
c
         neqi = id(ndf,i)
         if(neqi .gt. 0 .and. neqi .le. neq) then
c
c ...          Loop nas equacoes do vertice i:
c
           do 193 kk = 1, ndf    
             neqj  = id(kk,i)
             if (neqj .le. 0) goto 193
             if (neqj .lt. neqi) then
               if (lower) then
                 n = n + 1
                 ja(n) = neqj
               endif
             elseif (neqj .eq. neqi) then
               if (diag) then
                 n = n + 1
                 ja(n) = neqj
               endif
             elseif (neqj .gt. neqi) then
               if (upper) then
                 n = n + 1
                 ja(n) = neqj
               endif
             endif                  
  193      continue
c
c ...          Loop nos vertices conectados ao vertice i:
c
           do 195 k = ip(i), ip(i+1)-1
             j = ips(k)
c
c ...             Loop nas equacoes do vertice j:
c
             do 200 jj = 1, ndf 
               neqj  = id(jj,j)
               
               if (neqj .le. 0) goto 200
               if (neqj .lt. neqi) then
                 if (lower) then
                   n = n + 1
                   ja(n) = neqj
                 endif
               elseif (neqj .gt. neqi) then
                 if (upper) then
                   n = n + 1
                   ja(n) = neqj
                 endif
               endif
  200        continue
  195      continue
c ......................................................................
         endif
  190 continue
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csriaup(id,num,ip,ips,ia,iaup,nnode,ndf,neq,nequ,
     .                   neqp,nad,nadup,nadr,lower,diag,upper,right)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 12/12/2015                                    * 
c * ------------------------------------------------------------------ *  
c * CSRIAUP: monta o arranjo ia do formato CSRC para os blocos (uu pp) *
c * e csr para o bloco up                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * id(ndf,nnode)- numeracao global das equacoes                       *
c * num(nnode)   - renumeracao nodal                                   *
c * ip(nnode+1)  - ip(i) indica a posicao em ips do primeiro no        *
c *                conectado ao no i                                   *
c * ips(ipos)    - contem as conetividades nodais de cada no           *
c * nnode - numero de nos                                              *
c * ndf   - numero max. de graus de liberdade por no                   *
c * neq   - numero de equacoes                                         *
c * nequ  - numero de equacoes de deslocamentos                        *
c * neqp  - numero de equacoes de  pressao                             *
c * lower = .true.  -> inclui a parte triangular inferior no csr       *
c * diag  = .true.  -> inclui a diagonal no csr                        *
c * upper = .true.  -> inclui a parte triangular superior no csr       *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro         *
c *                   coeficiente nao-nulo da equacao i                *
c * iaup(neqp+1) - ia(i) informa a posicao no vetor a do primeiro      *
c *                   coeficiente nao-nulo da equacao i                *
c * nad   - numero de coeficientes nao nulos dos blocos uu e pp        *
c * nadup - numero de coeficientes nao nulos do bloco up               *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Kuu e Kpp juntos e Kpu separado                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,neq,nequ,neqp,nad,nadup,nadr
      integer id(ndf,*),ip(*),ips(*),ia(*),iaup(*),num(*)
      integer i,j,k,ii,jj,kk,neqi,neqj,no
      logical lower,diag,upper,right
c ......................................................................
c
c ... Inicializa o arranjo ia:
c
      do 50 i = 1, neq+1
         ia(i) = 0
c ... overllaping
c        if (right) ia(i+neq+1) = 0
   50 continue
c ----------------------------------------------
c
c ... Loop nos nos:
c
      do 150 no = 1, nnode
         i = num(no)
         do 140 ii = 1, ndf
c
c ...    Loop nas equacoes do no i:
c
            neqi = id(ii,i)
c ... equacoes uu
            if (neqi .gt. 0 .and. neqi .le. nequ) then
c ----------------------------------------------
c
c ...          Loop nas equacoes de deslocamantos do no i:
c
               do 100 kk = 1, ndf-1
                  neqj = id(kk,i)
                  if (neqj .le. 0) goto 100
                  if (neqj .lt. neqi) then
                     if (lower) ia(neqi) = ia(neqi) + 1
                  elseif (neqj .eq. neqi) then
                     if (diag) ia(neqi)  = ia(neqi) + 1
                  elseif (neqj .gt. neqi) then
                     if (upper) ia(neqi) = ia(neqi) + 1
                  endif
  100          continue
c
c ...          Loop nos nos conectados ao no i:
c
               do 120 k = ip(i), ip(i+1)-1
                  j = ips(k)
c
c ...             Loop nas equacoes de deslocamentos do no j:
c
                  do 110 jj = 1, ndf-1
                     neqj = id(jj,j)
                     if (neqj .le. 0) goto 110
                     if (neqj .lt. neqi) then
                        if (lower) ia(neqi) = ia(neqi) + 1
c ... overllaping
c                    elseif (neqj .gt. neq) then
c                       if (right) ia(neq+1+neqi) = ia(neq+1+neqi) + 1
                     elseif (neqj .gt. neqi) then
                        if (upper) ia(neqi) = ia(neqi) + 1
                     endif
  110             continue
  120          continue
c ----------------------------------------------
c
c ... equacoes pp
            else if(neqi .gt. nequ .and. neqi .le. neq) then 
c
c ...          equacao do no i (diagonal):
c
              if (diag) ia(neqi)  = ia(neqi) + 1
c             
c ...          Loop nos nos conectados ao no i:
c
              do 130 k = ip(i), ip(i+1)-1
                j = ips(k)
c
c ...         Loop nas equacoes de pressao do no j:
c
                neqj = id(ndf,j)
                if (neqj .le. 0) goto 130
                if (neqj .lt. neqi) then
                  if (lower) ia(neqi) = ia(neqi) + 1
c ... overllaping
c                elseif (neqj .gt. neq) then
c                  if (right) ia(neq+1+neqi) = ia(neq+1+neqi) + 1
                elseif (neqj .gt. neqi) then
                  if (upper) ia(neqi) = ia(neqi) + 1
                endif
  130         continue
            endif
  140    continue
  150 continue
c ----------------------------------------------
      
      do 200 i = neq, 1, -1
         ia(i+1) = ia(i)
c ... overllaping
c        if (right) ia(i+neq+2) = ia(i+neq+1)         
  200 continue
      ia(1) = 1
c ... overllaping
c     if (right) ia(neq+2) = 1
      do 210 i = 1, neq
         ia(i+1) = ia(i+1) + ia(i)
c ... overllaping
c        if (right) ia(i+neq+2) = ia(i+neq+2) + ia(i+neq+1)         
  210 continue
      nad  = ia(neq+1)-1
c ... overllaping
c     nadr = 0
c     if (right) nadr = ia(2*neq+2)-1
c ...................................................................... 
c
c ... Inicializa o arranjo iaup:
c
      do 250 i = 1, neqp+1
         iaup(i) = 0
c ... overllaping
c        if (right) ia(i+neq+1) = 0
  250 continue
c ----------------------------------------------
c
c ... Loop nos nos:
c
      do 300 no = 1, nnode
        i = num(no)
c ... equacoes up
        neqi = id(ndf,i)
        if(neqi .gt. 0) then
c ----------------------------------------------
c
c ...  Loop nas equacoes de deslocamento do no i:
c
            do 310 jj = 1, ndf-1
              neqj = id(jj,i)
              if (neqj .gt. 0) then     
                iaup(neqi-nequ) = iaup(neqi-nequ) + 1
              endif
  310     continue
c .....................................................................
c
c ...  Loop nos nos conectados ao no i:
c
          do 340 k = ip(i), ip(i+1)-1
            j = ips(k)
c
c ... Loop nas equacoes de deslocamento do no j:
c
            do 330 jj = 1, ndf-1
              neqj = id(jj,j)
              if (neqj .gt. 0) then 
                iaup(neqi-nequ) = iaup(neqi-nequ) + 1
              endif
  330       continue
  340     continue
        endif     
c ----------------------------------------------
  300 continue
c ----------------------------------------------
      do 400 i = neqp, 1, -1
        iaup(i+1) = iaup(i)
  400 continue
      iaup(1) = + 1
      do 410 i = 1, neqp
        iaup(i+1) = iaup(i+1) + iaup(i)
  410 continue
      nadup = iaup(neqp+1)-1
c ......................................................................      
c **********************************************************************
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csrjaup(id,num,ip,ips,ja,jaup,nnode,nnodev,ndf,
     .                   neq,nequ,
     .                   nad,nadup,lower,diag,upper,right)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 12/12/2015                                    * 
c * ------------------------------------------------------------------ * 
c * CSRJAUP: monta o arranjo ja do formato CSRC para os blocos         *
c * uu e pp e csr para o bloco up                                      * 
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * id(ndf,nnode)- numeracao global das equacoes                       *
c * num(nnode)   - renumeracao nodal                                   *
c * ip(nnode+1)  - ip(i) indica a posicao em ips do primeiro           *
c * ips(ipos)    - contem as conetividades nodais de cada no           *
c * nnode - numero de nos                                              *
c * ndf   - numero max. de graus de liberdade por no                   *
c * neq   - numero de equacoes                                         *
c * nequ  - numero de equacoes de deslocamentos                        *
c * nad   - numero de coeficientes nao nulos                           *
c * nadup - numero de coeficientes nao nulos do bloco up               *
c * lower = .true.  -> inclui a parte triangular inferior no csr       *
c * diag  = .true.  -> inclui a diagonal no csr                        *
c * upper = .true.  -> inclui a parte triangular superior no csr       *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ja(nad) - ja(k) informa a coluna do coeficiente que ocupa          *
c *           a posicao k no vetor a                                   *
c * jaup(nadup) - ja(k) informa a coluna do coeficiente que ocupa      *
c *           a posicao k no vetor a                                   *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Kuu e Kpp juntos e Kpu separado                                    *
c **********************************************************************
      implicit none
      integer nnode,nnodev,ndf,neq,nequ,nad,nadup
      integer id(ndf,*),num(*),ip(*),ips(*),ja(*),jaup(*)
      integer i,j,k,ii,jj,kk,neqi,neqj,no,n
      logical lower,diag,upper,right
c ......................................................................
      n = 0
c ... Equacoes uu
c
c ... Loop nos nos:
c
      do 160 no = 1, nnode
         i = num(no)
c
c ...    Loop nas equacoes de deslocamento do no i:
c
         do 150 ii = 1, ndf - 1
            neqi = id(ii,i)
            if (neqi .gt. 0 .and. neqi .le. nequ) then
c ----------------------------------------------------------
c
c ...          Loop nas equacoes de deslocamento do no i:
c
               do 100 kk = 1, ndf-1
                  neqj  = id(kk,i)
                  if (neqj .le. 0) goto 100
                  if (neqj .lt. neqi) then
                     if (lower) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  elseif (neqj .eq. neqi) then
                     if (diag) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  elseif (neqj .gt. neqi) then
                     if (upper) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  endif                  
  100          continue
c
c ...          Loop nos conectados ao no i:
c
               do 120 k = ip(i), ip(i+1)-1
                  j = ips(k)
c
c ...             Loop nas equacoes de deslocamento do no j:
c
                  do 110 jj = 1, ndf - 1
                     neqj  = id(jj,j)
                     if (neqj .le. 0) goto 110
                     if (neqj .lt. neqi) then
                        if (lower) then
                           n = n + 1
                           ja(n) = neqj
                        endif
c ... overllaping
c                    elseif (neqj .gt. neq) then
c                       if (right) then
c                          m = m + 1
c                          ja(nad+m) = neqj
c                       endif
                     elseif (neqj .gt. neqi) then
                        if (upper) then
                           n = n + 1
                           ja(n) = neqj
                        endif
                     endif
  110             continue
  120          continue
            endif   
c ----------------------------------------------------------
c
  150    continue
  160 continue
c ......................................................................
c
c ... Equacoes pp
c
c ... Loop nos nos:
c
      do 190 no = 1, nnode
         i = num(no)
c
c ...    equacao da pressao do no i:
c
         neqi = id(ndf,i)
c ---------------------------------------------------------      
         if(neqi .gt. nequ .and. neqi .le. neq) then  
c
c ...    equacao da pressao do no i (diagonal):
c
           if (diag) then
             n     = n + 1
             ja(n) = neqi
           endif
c
c ...          Loop nos nos conectados ao no i:
c
           do 170 k = ip(i), ip(i+1)-1
             j = ips(k)
c
c ...             Loop nas equacoes de pressao do no j:
c
             neqj  = id(ndf,j)
             if (neqj .le. 0) goto 170
             if (neqj .lt. neqi) then
             if (lower) then
               n = n + 1
               ja(n) = neqj
             endif
c ... overllaping
c            elseif (neqj .gt. neq) then
c              if (right) then
c                m = m + 1
c                ja(nad+m) = neqj
c              endif
             elseif (neqj .gt. neqi) then
               if (upper) then
                 n = n + 1
                 ja(n) = neqj
               endif
             endif
  170      continue
         endif
  190 continue
c .....................................................................
c
c ...
      n = 0
c ... loop nas equacoes pu
c
c ... Loop nos nos:
c
      do 290 no = 1, nnode
         i = num(no)
         neqi = id(ndf,i)
c ---------------------------------------------------------      
         if(neqi .gt. 0) then 
c
c ... Loop nas equacoes de deslocamento do no i:
c
             do 310 jj = 1, ndf-1
               neqj  = id(jj,i)
               if (neqj .le. 0 ) goto 310
               n = n + 1
               jaup(n) = neqj
  310        continue
c .....................................................................
c
c
c ... Loop nos nos conectados ao no i:
c
           do 270 k = ip(i), ip(i+1)-1
             j = ips(k)
c
c ... Loop nas equacoes de deslocamento do no j:
c
             do 320 jj = 1, ndf-1
               neqj  = id(jj,j)
               if (neqj .le. 0) goto 320
               n = n + 1
               jaup(n) = neqj
  320        continue
  270      continue
         endif
  290 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csriaup2(id   ,num  ,ip   ,ips
     .                   ,iau  ,iap  ,iapu
     .                   ,nnode,ndf
     .                   ,neq  ,nequ ,neqp
     .                   ,nad  ,naduu,nadpp,nadup
     .                   ,lower,diag ,upper)
c **********************************************************************
c * Data de criacao    : 01/01/2015                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *                                                                   *
c * CSRIAUP2: monta o arranjo ia do formato CSRC para os blocos        *
c * uu e pp e csr para o bloco up separados                            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * id(ndf,nnode)- numeracao global das equacoes                       *
c * num(nnode)   - renumeracao nodal                                   *
c * ip(nnode+1)  - ip(i) indica a posicao em ips do primeiro no        *
c *                conectado ao no i                                   *
c * ips(ipos)    - contem as conetividades nodais de cada no           *
c * nnode - numero de nos                                              *
c * iau(nequ+1) - nao definido                                         *
c * iap(neqp+1) - nao definido                                         *
c * iaup(neq+1) - nao definido                                         *
c * ndf   - numero max. de graus de liberdade por no                   *
c * neq   - numero de equacoes                                         *
c * nequ  - numero de equacoes de deslocamentos                        *
c * neqp  - numero de equacoes de  pressao                             *
c * ja(nad) - ja(k) informa a coluna do coeficiente que ocupa          *
c *           a posicao k no vetor a do bloco Kuu                      *
c * lower = .true.  -> inclui a parte triangular inferior no csr       *
c * diag  = .true.  -> inclui a diagonal no csr                        *
c * upper = .true.  -> inclui a parte triangular superior no csr       *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * iau(nequ+1) - ia(i) informa a posicao no vetor a do primeiro       *
c *               coeficiente nao-nulo da equacao i do bloco Kuu       *
c * iap(neqp+1) - ia(i) informa a posicao no vetor a do primeiro       *
c *               coeficiente nao-nulo da equacao i do bloco Kpp       *
c * iaup(neq+1) - ia(i) informa a posicao no vetor a do primeiro       *
c *               coeficiente nao-nulo da equacao i  do bloco Kup      *
c * nad   - numero de coeficientes nao nulos dos blocos Kuu e Kpp      *
c * naduu - numero de coeficientes nao nulos do bloco Kuu              *
c * nadpp - numero de coeficientes nao nulos do bloco Kpp              *
c * nadup - numero de coeficientes nao nulos do bloco Kup              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Kuu ,Kpp e Kpu separados                                           *
c **********************************************************************
      implicit none
      integer nnode,ndf,neq,nequ,neqp,nad,naduu,nadpp,nadup
      integer id(ndf,*),ip(*),ips(*),iau(*),iap(*),iapu(*),num(*)
      integer i,j,k,ii,jj,kk,neqi,neqii,neqj,no
      logical lower,diag,upper
c ......................................................................
c
c ... Inicializa o arranjo ia:
c
      iau(1:nequ+1) = 0
c 
      iap(1:neqp+1) = 0
c ----------------------------------------------
c
c ... Loop nos nos:
c
      do 150 no = 1, nnode
         i = num(no)
         do 140 ii = 1, ndf
c
c ...    Loop nas equacoes do no i:
c
            neqi = id(ii,i)
c ... equacoes uu
            if (neqi .gt. 0 .and. neqi .le. nequ) then
c ----------------------------------------------
c
c ...          Loop nas equacoes de deslocamantos do no i:
c
               do 100 kk = 1, ndf-1
                  neqj = id(kk,i)
                  if (neqj .le. 0) goto 100
                  if (neqj .lt. neqi) then
                     if (lower) iau(neqi) = iau(neqi) + 1
                  elseif (neqj .eq. neqi) then
                     if (diag) iau(neqi)  = iau(neqi) + 1
                  elseif (neqj .gt. neqi) then
                     if (upper) iau(neqi) = iau(neqi) + 1
                  endif
  100          continue
c
c ...          Loop nos nos conectados ao no i:
c
               do 120 k = ip(i), ip(i+1)-1
                  j = ips(k)
c
c ...             Loop nas equacoes de deslocamentos do no j:
c
                  do 110 jj = 1, ndf-1
                     neqj = id(jj,j)
                     if (neqj .le. 0) goto 110
                     if (neqj .lt. neqi) then
                        if (lower) iau(neqi) = iau(neqi) + 1
                     elseif (neqj .gt. neqi) then
                        if (upper) iau(neqi) = iau(neqi) + 1
                     endif
  110             continue
  120          continue
c ----------------------------------------------
c
c ... equacoes pp
            else if(neqi .gt. nequ .and. neqi .le. neq) then 
              neqii = neqi-nequ
c
c ...          equacao do no i (diagonal):
c
              if (diag) iap(neqii)  = iap(neqii) + 1
c             
c ...          Loop nos nos conectados ao no i:
c
              do 130 k = ip(i), ip(i+1)-1
                j = ips(k)
c
c ...         Loop nas equacoes de pressao do no j:
c
                neqj = id(ndf,j)
                if (neqj .le. 0) goto 130
                if (neqj .lt. neqi) then
                  if (lower) iap(neqii) = iap(neqii) + 1
                elseif (neqj .gt. neqi) then
                  if (upper) iap(neqii) = iap(neqii) + 1
                endif
  130         continue
            endif
  140    continue
  150 continue
c .....................................................................
c
c ...
      do 200 i = nequ, 1, -1
         iau(i+1) = iau(i)
  200 continue
      iau(1) = 1
c .....................................................................
c
c ...
      do 201 i = neqp, 1, -1
         iap(i+1) = iap(i)
  201 continue
      iap(1) = 1
c .....................................................................
c
c ...
      do 210 i = 1, nequ
         iau(i+1) = iau(i+1) + iau(i)
  210 continue
      naduu= iau(nequ+1)-1
c .....................................................................
c
c ...
      do 211 i = 1, neqp
         iap(i+1) = iap(i+1) + iap(i)
  211 continue
      nadpp= iap(neqp+1)-1
c ...................................................................... 
c
c ... Inicializa o arranjo iaup:
c
      do 250 i = 1, neqp+1
         iapu(i) = 0
  250 continue
c ----------------------------------------------
c
c ... Loop nos nos:
c
      do 300 no = 1, nnode
        i = num(no)
c ... equacoes up
        neqi = id(ndf,i)
        if(neqi .gt. 0) then
            neqii = neqi-nequ
c ----------------------------------------------
c
c ...  Loop nas equacoes de deslocamento do no i:
c
            do 310 jj = 1, ndf-1
              neqj = id(jj,i)
              if (neqj .gt. 0) then     
                iapu(neqii) = iapu(neqii) + 1
              endif
  310     continue
c .....................................................................
c
c ...  Loop nos nos conectados ao no i:
c
          do 340 k = ip(i), ip(i+1)-1
            j = ips(k)
c
c ... Loop nas equacoes de deslocamento do no j:
c
            do 330 jj = 1, ndf-1
              neqj = id(jj,j)
              if (neqj .gt. 0) then 
                iapu(neqii) = iapu(neqii) + 1
              endif
  330       continue
  340     continue
        endif     
c ----------------------------------------------
  300 continue
c ----------------------------------------------
      do 400 i = neqp, 1, -1
        iapu(i+1) = iapu(i)
  400 continue
      iapu(1) = 1
      do 410 i = 1, neqp
        iapu(i+1) = iapu(i+1) + iapu(i)
  410 continue
      nadup = iapu(neqp+1)-1
c ...
      nad   = naduu + nadpp + nadup
c ......................................................................      
c **********************************************************************
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine csrjaup2(id    ,num,ip  ,ips
     .                   ,jau   ,jap,jaup
     .                   ,nnode ,ndf
     .                   ,neq   ,nequ
     .                   ,lower ,diag,upper,right)
c **********************************************************************
c * Data de criacao    : 01/01/2015                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *                                                                      *
c * CSRJAUP2: monta o arranjo ja do formato CSRC para os blocos        *
c * uu e pp e csr para o bloco up                                      * 
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * id(ndf,nnode)- numeracao global das equacoes                       *
c * num(nnode)   - renumeracao nodal                                   *
c * ip(nnode+1)  - ip(i) indica a posicao em ips do primeiro           *
c * ips(ipos)    - contem as conetividades nodais de cada no           *
c * jau          - nao definido                                        *
c * jap          - nao definido                                        *
c * japu         - nao definido                                        *
c * nnode        - numero de nos                                       *
c * ndf          - numero max. de graus de liberdade por no            *
c * neq          - numero de equacoes                                  *
c * nequ         - numero de equacoes de deslocamentos                 *
c * lower = .true.  -> inclui a parte triangular inferior no csr       *
c * diag  = .true.  -> inclui a diagonal no csr                        *
c * upper = .true.  -> inclui a parte triangular superior no csr       *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * jau(nadu) - ja(k) informa a coluna do coeficiente que ocupa        *
c *           a posicao k no vetor a do bloco Kuu                      *
c * jap(nadp) - ja(k) informa a coluna do coeficiente que ocupa        *
c *           a posicao k no vetor a do bloco Kpp                      *
c * jaup(nadup) - ja(k) informa a coluna do coeficiente que ocupa      *
c *           a posicao k no vetor a do bloco Kpu                      *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Kuu ,Kpp e Kpu separados                                           *
c **********************************************************************
      implicit none
      integer nnode,nnodev,ndf,neq,nequ
      integer id(ndf,*),num(*),ip(*),ips(*),jau(*),jap(*),jaup(*)
      integer i,j,k,ii,jj,kk,neqi,neqj,no,n
      logical lower,diag,upper,right
c ......................................................................
      n = 0
c ... Equacoes uu
c
c ... Loop nos nos:
c
      do 160 no = 1, nnode
         i = num(no)
c
c ...    Loop nas equacoes de deslocamento do no i:
c
         do 150 ii = 1, ndf - 1
            neqi = id(ii,i)
            if (neqi .gt. 0 .and. neqi .le. nequ) then
c ----------------------------------------------------------
c
c ...          Loop nas equacoes de deslocamento do no i:
c
               do 100 kk = 1, ndf-1
                  neqj  = id(kk,i)
                  if (neqj .le. 0) goto 100
                  if (neqj .lt. neqi) then
                     if (lower) then
                        n = n + 1
                        jau(n) = neqj
                     endif
                  elseif (neqj .eq. neqi) then
                     if (diag) then
                        n = n + 1
                        jau(n) = neqj
                     endif
                  elseif (neqj .gt. neqi) then
                     if (upper) then
                        n = n + 1
                        jau(n) = neqj
                     endif
                  endif                  
  100          continue
c
c ...          Loop nos conectados ao no i:
c
               do 120 k = ip(i), ip(i+1)-1
                  j = ips(k)
c
c ...             Loop nas equacoes de deslocamento do no j:
c
                  do 110 jj = 1, ndf - 1
                     neqj  = id(jj,j)
                     if (neqj .le. 0) goto 110
                     if (neqj .lt. neqi) then
                        if (lower) then
                           n = n + 1
                           jau(n) = neqj
                        endif
                     elseif (neqj .gt. neqi) then
                        if (upper) then
                           n = n + 1
                           jau(n) = neqj
                        endif
                     endif
  110             continue
  120          continue
            endif   
c ----------------------------------------------------------
c
  150    continue
  160 continue
c ......................................................................
c
c ......................................................................
      n = 0
c ... Equacoes pp
c
c ... Loop nos nos:
c
      do 190 no = 1, nnode
         i = num(no)
c
c ...    equacao da pressao do no i:
c
         neqi = id(ndf,i)
c ---------------------------------------------------------      
         if(neqi .gt. nequ .and. neqi .le. neq) then  
c
c ...    equacao da pressao do no i (diagonal):
c
           if (diag) then
             n     = n + 1
             jap(n) = neqi - nequ
           endif
c
c ...          Loop nos nos conectados ao no i:
c
           do 170 k = ip(i), ip(i+1)-1
             j = ips(k)
c
c ...             Loop nas equacoes de pressao do no j:
c
             neqj  = id(ndf,j)
             if (neqj .le. 0) goto 170
             if (neqj .lt. neqi) then
             if (lower) then
               n = n + 1
               jap(n) = neqj - nequ
             endif
             elseif (neqj .gt. neqi) then
               if (upper) then
                 n = n + 1
                 jap(n) = neqj - nequ
               endif
             endif
  170      continue
         endif
  190 continue
c .....................................................................
c
c ...
      n = 0
c ... loop nas equacoes pu
c
c ... Loop nos nos:
c
      do 290 no = 1, nnode
         i = num(no)
         neqi = id(ndf,i)
c ---------------------------------------------------------      
         if(neqi .gt. 0) then 
c
c ... Loop nas equacoes de deslocamento do no i:
c
             do 310 jj = 1, ndf-1
               neqj  = id(jj,i)
               if (neqj .le. 0 ) goto 310
               n = n + 1
               jaup(n) = neqj
  310        continue
c .....................................................................
c
c
c ... Loop nos nos conectados ao no i:
c
           do 270 k = ip(i), ip(i+1)-1
             j = ips(k)
c
c ... Loop nas equacoes de deslocamento do no j:
c
             do 320 jj = 1, ndf-1
               neqj  = id(jj,j)
               if (neqj .le. 0) goto 320
               n = n + 1
               jaup(n) = neqj
  320        continue
  270      continue
         endif
  290 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c * GET_NADS :                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ia    - vetor ia dos blocos Kuu e Kpp                           *
c *    neq   - numero de equacoes                                      *
c *    nequ  - numero de equacoes de deslocamentos                     *
c *    naduu - nao definido                                            *
c *    nadpp - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    naduu - numero de coeficientes nao nulos do bloco uu            *
c *    nadpp - numero de coeficientes nao nulos do bloco pp            *
c *                                                                    *
c **********************************************************************
      subroutine get_nads(ia,naduu,nadpp,neq,nequ)
      implicit none
      integer ia(*),naduu,nadpp,neq,nequ
      naduu = ia(nequ+1) - ia(1)
      nadpp = ia(neq+1)  - ia(nequ+1)
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 15/12/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * SORT_GRAPH_CSR: reordena as coeficientes das linhas em ordem       *
c * crescente de coluna                                                *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ia(*) - informa a posicao do primeiro coeficiente nao-nulo         * 
c *                da equacao i                                        *
c * ja(*) - informa a coluna do coeficiente                            *
c * neq   - numero de equacoes                                         *
c * nad   - numero total de termo nao nulos fora da diagonal da parte  *
c *         quadrada da matriz de coeficientes                         *
c * ovlp  - overlaping                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * ia (integer*8),ja (integer), nad(interger*8) e neq (integer)       *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine sort_graph_csr(ia,ja,neq,nad,ovlp)
      implicit none
      integer*8 ia(*),nad
      integer ja(*),neq
      logical ovlp
      call sortgraphv2(ia,ja,neq)
      if(ovlp) call sortgraphv2(ia(neq+2),ja(nad+1),neq)  
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine printIaJa(ia,ja,iar,jar,neq,neqr,nad,nadr)
      implicit none
      integer neq,neqr,nad,nadr,i
      integer ia(*),ja(*),iar(*),jar(*)
c ...      
      open(14,file='ia.dat',action='write')
      open(15,file='ja.dat',action='write')
c ......................................................................
c
c ...
      do i = 1, neq
        write(14,*)i,ia(i),ia(i+1) - ia(i)    
      enddo
      write(14,*)i,ia(neq+1)    
      if( neqr .ne. 0 ) then
        do i = 1, neqr
          write(14,*)i,iar(i),iar(i+1) - iar(i)    
        enddo
        write(14,*)i,iar(neqr+1)    
      endif
      do i = 1, nad
        write(15,*)i,ja(i)    
      enddo
      if( nadr .ne. 0) then
        do i = 1, nadr
          write(15,*)i,jar(i)    
        enddo
      endif
      close(14)
      close(15)
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine csrc_to_full(ia,ja,al,ad,a,neq)
      implicit none
      integer ia(*),ja(*)
      real*8 ad(*),al(*),a(neq,*)
      integer neq,i,j,jak,kk
c ...
      do j = 1, neq
        do i = 1, neq
          a(i,j) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ...
      do i = 1, neq
        a(i,i) = ad(i)
      enddo
c .....................................................................
c
c ...
      do i = 1, neq
        do j = ia(i), ia(i+1) - 1
          jak      = ja(j)
c ... caso onde o csr tem termos nulos 
          if( al(j) .eq. 0.d0) then
            a(i,jak) = ad(i)*1.d-16
          else 
            a(i,jak) = al(j)
          endif
        enddo  
      enddo
c .....................................................................
c
c ...
      do j = 1, neq
        do i = j+1 , neq
          a(j,i) = a(i,j)
        enddo
      enddo
c .....................................................................
      open(15,file='afullf.matrix')
      do i = 1, neq
        write(15,'(100d10.2)')(a(i,j),j=1,neq)
      enddo 
      close(15)
      return
      end
c *********************************************************************
