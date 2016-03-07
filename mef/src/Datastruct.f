      subroutine datastruct(ix,id,num,nnode,nnodev,numel,nen,ndf,nst
     .               ,neq,nequ,neqp,stge
     .               ,unsym,nad,naduu,nadpp,nadpu
     .               ,i_ia,i_ja,i_au,i_al,i_ad,ija,ja,au
     .               ,al,ad,ovlp,n_blocks_up,dualCsr)
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
c *    nnodev- numero de nos vertices da particao                      *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - nen*ndf                                                 *
c *    neq   - numero de equacoes                                      *
c *    nequ  - numero de equacoes de deslocamentos                     *
c *    neqp  - numero de equacoes de  pressao                          *
c *    stge  - estrutura de dados, 1 = CSR, 2 = ARESTAS, 3 = EBE,      *
c *                                4 = skyline                         *
c *    unsym - flag para matrizes nao simetricas                       *
c *    n_blocks_up = numero de blocos                                  *
c *                  1 - ( [Kuu  Kpp]  )                               *
c *                  2 - ( [Kuu, Kpp] e [kpu] )                        *
c *                  3 - ( [Kuu], [Kpp] e [kpu])                       *     
c *    dualCsr - flag para matrizes (Kuu+Kpp) e Kup separadas          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    i_ia  = ponteiro para o arranjo ia(neq+1) do CSR (stge = 1)     *
c *          = ponteiro para o arranjo edge(2,nedge)    (stge = 2)     *
c *          = nao definido (stge = 3)                                 *
c *          = i_ja         (stge = 4)                                 *
c *    i_ja  = ponteiro para o arranjo ja do CSR        (stge = 1)     *
c *          = ponteiro para o arranjo ipedg(nnode+1)   (stge = 2)     *
c *          = nao definido (stge = 3)                                 *
c *          = ponteiro da diagonal  (stge = 4)                        *
c *    i_au  - ponteiro para o arranjo au(nad)                         *
c *    i_al  - ponteiro para o arranjo al(nad)                         *
c *    i_ad  - ponteiro para a diagonal                                *
c * dualCsr = true                                                     *
c *    nad   - numero de coeficientes nao nulos dos blocos uu e pp     *
c *    naduu - numero de coeficientes nao nulos do bloco uu            *
c *    nadpp - numero de coeficientes nao nulos do bloco pp            *
c *    nadpu - numero de coeficientes nao nulos do bloco pu            *
c * dualCsr = false                                                    *
c *    nad   - numero de coeficientes nao nulos                        *
c *                                                                    *
c *    Obs: esta rotina deve ser utilizada com a estrutura de dados    *
c *         definida no arquivo 'common.h'                             *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer ix(nen+1,*),id(ndf,*),num(*),nnode,nnodev
      integer numel,nen,ndf,nst,neq, n_blocks_up
      integer stge,nad,naduu,nadpp,nadpu,nequ,neqp
c ... ponteiros      
      integer*8 i_ia,i_ja,i_au,i_al,i_ad
      integer*8 i_bd,i_lde
c .....................................................................      
      integer nedge,nste
      logical unsym,bdfl,ovlp,dualCsr
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
         call csrstruct(id,ix,num,nnode,nnodev
     .                 ,numel,nen,ndf,neq,nequ,neqp
     .                 ,i_ia,i_ja,nad,naduu,nadpp,nadpu
     .                 ,.true.,.false.,.false.,ovlp,ija,ja
     .                 ,n_blocks_up,dualCsr)
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
      endif
c ......................................................................      
      return
      end
