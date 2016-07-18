      subroutine datastruct_pm(ix,id,num,nnode,nnodev,numel,nen,ndf,nst
     .               ,neq,nequ,neqp,stge
     .               ,unsym,nad,naduu,nadpp,nadpu
     .               ,i_ia,i_ja,i_au,i_al,i_ad,ija,ja,au
     .               ,al,ad,ovlp,n_blocks_up,block_pu ,block_pu_sym)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 26/06/2016                                    * 
c * ------------------------------------------------------------------ * 
c * DATASTRUCT: monta a estrutura de dados para a matriz de            *
c *               coeficientes do sistema de equacoes de acordo com    *
c *               o formato especificado.                              *
c * ------------------------------------------------------------------ *                                                             *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ix(nen+1,numel) - conetividades nodais                             *
c * id(ndf,nnode)   - numeracao global das equacoes                    *
c * num(nnode)      - arranjo auxiliar temporario                      *
c * nnode - numero de nos total da particao                            *
c * nnodev- numero de nos vertices da particao                         *
c * numel - numero de elementos                                        *
c * nen   - numero de nos por elemento                                 *
c * ndf   - numero max. de graus de liberdade por no                   *
c * nst   - nen*ndf                                                    *
c * neq   - numero de equacoes                                         *
c * nequ  - numero de equacoes de deslocamentos                        *
c * neqp  - numero de equacoes de  pressao                             *
c * stge  - estrutura de dados, 1 = CSR, 2 = ARESTAS, 3 = EBE,         *
c *                             4 = skyline                            *
c * unsym - flag para matrizes nao simetricas                          *
c * n_blocks_up = numero de blocos                                     *
c *               1 - ( [Kuu  Kpp]  )                                  *
c *               2 - ( [Kuu, Kpp] e [kpu] )                           *
c *               3 - ( [Kuu], [Kpp] e [kpu])                          *     
c * block_pu     - flag para matrizes Kuu,kpp e Kup nao simetricos     *
c * block_pu_sym - flag para matrizes Kuu, Kpp e Kup simetricos        *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * i_ia  = ponteiro para o arranjo ia(neq+1) do CSR (stge = 1)        *
c *       = ponteiro para o arranjo edge(2,nedge)    (stge = 2)        *
c *       = nao definido (stge = 3)                                    *
c *       = i_ja         (stge = 4)                                    *
c * i_ja  = ponteiro para o arranjo ja do CSR        (stge = 1)        *
c *       = ponteiro para o arranjo ipedg(nnode+1)   (stge = 2)        *
c *       = nao definido (stge = 3)                                    *
c *       = ponteiro da diagonal  (stge = 4)                           *
c * i_au  - ponteiro para o arranjo au(nad)                            *
c * i_al  - ponteiro para o arranjo al(nad)                            *
c * i_ad  - ponteiro para a diagonal                                   *
c * block_pu = true                                                    *     *
c *    nad   - numero de coeficientes nao nulos dos blocos uu e pp     *
c *    naduu - numero de coeficientes nao nulos do bloco uu            *
c *    nadpp - numero de coeficientes nao nulos do bloco pp            *
c *    nadpu - numero de coeficientes nao nulos do bloco pu            *
c * block_pu_sym = true | false; block_pu = false                      *                                    *
c *    nad   - numero de coeficientes nao nulos                        *
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
c * monta o bloco kuu, kpp e kup juntos se considera a estrutura       * 
c * blocada, i.e., graus de liberda juntos                             *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer ix(nen+1,*),id(ndf,*),num(*),nnode,nnodev
      integer numel,nen,ndf,nst,neq, n_blocks_up
      integer stge,nad,nnzr,naduu,nadpp,nadpu,nequ,neqp
c ... ponteiros      
      integer*8 i_ia,i_ja,i_au,i_al,i_ad
      integer*8 i_bd,i_lde
c .....................................................................      
      integer nedge,nste
      logical unsym,bdfl,ovlp,block_pu ,block_pu_sym
      character*8 ija,ja,au,al,ad
c ......................................................................
      i_ia    = 1
      i_ja    = 1
      i_au    = 1
      i_al    = 1
      i_ad    = 1
      i_bd    = 1
      i_lde   = 1
      nad     = 0
      naduu   = 0
      naduu   = 0
      nadpu   = 0
      nedge   = 0
      nste    = 0
      bdfl    = .false.
c ......................................................................
      if(stge .eq. 1) then
c      
c ...    Armazenamento CSR da matriz de coeficientes:
c        --------------------------
c        | ia | ja | al | au | ad |
c        --------------------------
c
c ...    estrutura de dados do csr:
c
         call csrstruct_pm(id,ix,num,nnode,nnodev
     .                    ,numel,nen,ndf,neq,nequ,neqp
     .                    ,i_ia,i_ja,nad,naduu,nadpp,nadpu
     .                    ,.true.,.false.,.false.,ovlp,ija,ja
     .                    ,n_blocks_up,block_pu ,block_pu_sym)
c     
c ...    matriz de coeficientes:
c
         i_al = alloc_8(al,1,nad+nadpu)
         i_au = i_al
         if(unsym) i_au = alloc_8(au,1,nad+nadpu) 
         i_ad = alloc_8(ad,1,neq)         
c ......................................................................
      elseif(stge .eq. 2) then
c
c ...    Armazenamento por ARESTAS:
c        --------------------------------------
c        | edge | ipedge | lde | au | bd | ad |
c        --------------------------------------
         call edgstruct(ix,nnode,numel,nen,nedge,i_ia,i_ja)
c
c ...    Numeracao das eqs. por aresta:
c         
         nste = 2*ndf         
         i_lde = alloc_4('lde     ',nste,nedge)
         call neqedg(ia(i_ia),id,ia(i_lde),nedge,ndf,nste)
c
c ...    Numero total de coeficientes armazenados por arestas:
c         
         nad = ndf*ndf*nedge
         if(unsym) nad = nad*2
         i_au = alloc_8('au      ',1,nad)
         i_al = i_au
c          
c ...    diagonal em blocos, se ndf > 1:
c
         if(ndf .gt. 1) then
            bdfl = .true.
            if (unsym) then
              i_bd = alloc_8('bd      ',ndf*ndf,nnode)
            else
              i_bd = alloc_8('bd      ',(ndf*(ndf+1))/2,nnode)
            endif
         endif
         i_ad   = alloc_8('ad      ',1,neq)         
c ......................................................................         
      elseif(stge .eq. 3) then
c      
c ...    Armazenamento EBE:
c        ---------------
c        | a | bd | ad |
c        ---------------
c
c ...    Numeracao das eqs. por elemento:
c
         i_lde = alloc_4('lde     ',nst,numel)
         call neqel(ix,id,ia(i_lde),numel,nen,ndf,nst)
c
c ...    Numero total de coeficientes armazenados por elementos:
c         
         nad = ((nst*(nst-ndf))/2)*numel
         if(unsym) nad = nad*2
c         
c ...    matrizes de elemento:
c
         i_au = alloc_8('au      ',1,nad)
         i_al = i_au
c                   
c ...    diagonal em blocos, se ndf > 1:
c
         if(ndf .gt. 1) then
            bdfl = .true.
            if (unsym) then
              i_bd = alloc_8('bd      ',ndf*ndf,nnode)
            else
              i_bd = alloc_8('bd      ',(ndf*(ndf+1))/2,nnode)
            endif
         endif
         i_ad   = alloc_8('ad      ',1,neq)         
c ......................................................................
c
c ... Armazenamento SKYLINE da matriz de coeficientes:
c
c     ---------------------
c     | ja | au | al | ad | 
c     ---------------------
c
      elseif(stge .eq. 4) then
c      
c ...    determina o perfil da matriz:
c
         i_ja = alloc_4('ja      ',1,neq)
         i_ia = i_ja
         call profil(ix,id,ia(i_ja),nnode,numel,nen,ndf,neq,nad,.true.)
c
c ...    matriz de coeficientes:
c
         i_au = alloc_8('au      ',1,nad)
         i_al = i_au
         if(unsym) i_al = alloc_8('al      ',1,nad)
         i_ad   = alloc_8('ad      ',1,neq)
      elseif(stge .eq. 5) then
c
c ...    Matriz diagonal de coeficientes:
c
         i_ad   = alloc_8('ad      ',1,neq)
c
c ......................................................................
      else if(stge .eq. 6) then
c      
c ...    Armazenamento CSR3 da matriz de coeficientes:
c        --------------------------
c        | ia | ja | a |
c        --------------------------
c
c ...    estrutura de dados do csr:
c
c ... 
        if(unsym) then
          call csrstruct(id,ix,num,nnode,numel,nen,ndf,neq,i_ia,i_ja
     .                  ,nad,nnzr,.true.,.true.,.true.,ovlp,ija,ja)
c .....................................................................
c
c ...
        else
          call csrstruct(id,ix,num,nnode,numel,nen,ndf,neq,i_ia,i_ja
     .                  ,nad,nnzr ,.false.,.true.,.true.,ovlp,ija,ja)
        endif
c .....................................................................
c     
c ...    matriz de coeficientes:
c
        i_ad = alloc_8(ad,1,nad)
        i_au = i_ad
        i_al = i_ad                        
      endif
c ......................................................................      
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine datastruct(ix   ,id  ,num  ,nnode
     .                     ,numel,nen ,ndf  ,nst
     .                     ,neq  ,stge,unsym,nad   ,nad1
     .                     ,i_ia ,i_ja,i_au ,i_al  ,i_ad
     .                     ,ija  ,ja  
     .                     ,au   ,al    ,ad
     .                     ,ovlp)
c **********************************************************************
c *                                                                    *
c *   DATASTRUCT: monta a estrutura de dados para a matriz de          *
c *               coeficientes do sistema de equacoes de acordo com    *
c *               o formato especificado.                              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    id(ndf,nnode)   - numeracao global das equacoes                 *
c *    num(nnode)      - arranjo auxiliar temporario                   *
c *    nnode - numero de nos total da particao                         *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - nen*ndf                                                 *
c *    neq   - numero de equacoes                                      *
c *    stge  - estrutura de dados, 1 = CSR, 2 = ARESTAS, 3 = EBE,      *
c *                                4 = skyline                         *
c *    unsym - flag para matrizes nao simetricas                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    i_ia  - ponteiro para o arranjo ia(neq+1) do CSR (stge = 1)     *
c *            ponteiro para o arranjo edge(2,nedge)    (stge = 2)     *
c *            nao definido (stge = 3)                                 *
c *            i_ja         (stge = 4)                                 *
c *    i_ja  - ponteiro para o arranjo ja do CSR        (stge = 1)     *
c *            ponteiro para o arranjo ipedg(nnode+1)   (stge = 2)     *
c *            nao definido (stge = 3)                                 *
c *            ponteiro da diagonal  (stge = 4)                        *
c *    i_au  - ponteiro para o arranjo au(nad)                         *
c *    i_al  - ponteiro para o arranjo al(nad)                         *
c *    i_ad  - ponteiro para a diagonal                                *
c *    nad   - numero de posicoes da matriz de coeficientes            *
c *    nad1  - numero de posicoes da matriz de coeficientes da parte   *
c *            retangular                                              *
c **********************************************************************
      use Malloc
      implicit none
      integer ix(nen+1,*),id(ndf,*),num(*),nnode,numel,nen,ndf,nst,neq
      integer stge,nad,nad1
c ... ponteiros      
      integer*8 i_ia,i_ja,i_au,i_al,i_ad
      integer*8 i_bd,i_lde
c .....................................................................      
      integer nedge,nste
      logical unsym,bdfl,ovlp
      character*8 ija,ja,au,al,ad
c ......................................................................
      i_ia    = 1
      i_ja    = 1
      i_au    = 1
      i_al    = 1
      i_ad    = 1
      i_bd    = 1
      i_lde   = 1
      nad     = 0
      nad1    = 1
      nedge   = 0
      nste    = 0
      bdfl    = .false.
c ......................................................................
      if(stge .eq. 1) then
c      
c ...    Armazenamento CSR da matriz de coeficientes:
c        --------------------------
c        | ia | ja | al | au | ad |
c        --------------------------
c
c ...    estrutura de dados do csr:
c
         call csrstruct(id,ix,num,nnode,numel,nen,ndf,neq,i_ia,i_ja,nad,
c     .                 nad1, lower ,  diag , upper,right)
c ... matvec novo:
     .                  nad1,.true.,.false.,.false.,ovlp,ija,ja)
c ... matvec antigo:
c     .                  nad1,.false.,.false.,.true.,.false.)
c     
c ...    matriz de coeficientes:
c
         i_al = alloc_8(al,1,nad+nad1)
         i_au = i_al
         if(unsym) i_au = alloc_8(au,1,nad) 
         i_ad = alloc_8(ad,1,neq)         
c ......................................................................
      elseif(stge .eq. 2) then
c
c ...    Armazenamento por ARESTAS:
c        --------------------------------------
c        | edge | ipedge | lde | au | bd | ad |
c        --------------------------------------
         call edgstruct(ix,nnode,numel,nen,nedge,i_ia,i_ja)
c
c ...    Numeracao das eqs. por aresta:
c         
         nste = 2*ndf         
         i_lde = alloc_4('lde     ',nste,nedge)
         call neqedg(ia(i_ia),id,ia(i_lde),nedge,ndf,nste)
c
c ...    Numero total de coeficientes armazenados por arestas:
c         
         nad = ndf*ndf*nedge
         if(unsym) nad = nad*2
         i_au = alloc_8('au      ',1,nad)
         i_al = i_au
c          
c ...    diagonal em blocos, se ndf > 1:
c
         if(ndf .gt. 1) then
            bdfl = .true.
            if (unsym) then
              i_bd = alloc_8('bd      ',ndf*ndf,nnode)
            else
              i_bd = alloc_8('bd      ',(ndf*(ndf+1))/2,nnode)
            endif
         endif
         i_ad   = alloc_8('ad      ',1,neq)         
c ......................................................................         
      elseif(stge .eq. 3) then
c      
c ...    Armazenamento EBE:
c        ---------------
c        | a | bd | ad |
c        ---------------
c
c ...    Numeracao das eqs. por elemento:
c
         i_lde = alloc_4('lde     ',nst,numel)
         call neqel(ix,id,ia(i_lde),numel,nen,ndf,nst)
c
c ...    Numero total de coeficientes armazenados por elementos:
c         
         nad = ((nst*(nst-ndf))/2)*numel
         if(unsym) nad = nad*2
c         
c ...    matrizes de elemento:
c
         i_au = alloc_8('au      ',1,nad)
         i_al = i_au
c                   
c ...    diagonal em blocos, se ndf > 1:
c
         if(ndf .gt. 1) then
            bdfl = .true.
            if (unsym) then
              i_bd = alloc_8('bd      ',ndf*ndf,nnode)
            else
              i_bd = alloc_8('bd      ',(ndf*(ndf+1))/2,nnode)
            endif
         endif
         i_ad   = alloc_8('ad      ',1,neq)         
c ......................................................................
c
c ... Armazenamento SKYLINE da matriz de coeficientes:
c
c     ---------------------
c     | ja | au | al | ad | 
c     ---------------------
c
      elseif(stge .eq. 4) then
c      
c ...    determina o perfil da matriz:
c
         i_ja = alloc_4('ja      ',1,neq)
         i_ia = i_ja
         call profil(ix,id,ia(i_ja),nnode,numel,nen,ndf,neq,nad,.true.)
c
c ...    matriz de coeficientes:
c
         i_au = alloc_8('au      ',1,nad)
         i_al = i_au
         if(unsym) i_al = alloc_8('al      ',1,nad)
         i_ad   = alloc_8('ad      ',1,neq)
      elseif(stge .eq. 5) then
c
c ...    Matriz diagonal de coeficientes:
c
         i_ad   = alloc_8('ad      ',1,neq)                        
c ......................................................................
      else if(stge .eq. 6) then
c      
c ...    Armazenamento CSR3 da matriz de coeficientes:
c        --------------------------
c        | ia | ja | a |
c        --------------------------
c
c ...    estrutura de dados do csr:
c
c ... 
        if(unsym) then
          call csrstruct(id,ix,num,nnode,numel,nen,ndf,neq,i_ia,i_ja
     .                  ,nad,nad1,.true.,.true.,.true.,ovlp,ija,ja)
c .....................................................................
c
c ...
        else
          call csrstruct(id,ix,num,nnode,numel,nen,ndf,neq,i_ia,i_ja
     .                  ,nad,nad1 ,.false.,.true.,.true.,ovlp,ija,ja)
        endif
c .....................................................................
c     
c ...    matriz de coeficientes:
c
        i_ad = alloc_8(ad,1,nad)
        i_au = i_ad
        i_al = i_ad                        
      endif
c ......................................................................  
      return
      end
