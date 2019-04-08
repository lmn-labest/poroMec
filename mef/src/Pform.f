      subroutine pform_pm(ix      ,eloads      ,eloadsp
     1                   ,ie      ,e
     2                   ,x       ,id         ,ia       ,ja
     3                   ,au      ,al         ,ad       ,b
     4                   ,u0      ,u          ,tx1      ,tx2
     5                   ,deps    ,dp         ,plastic  ,elplastic
     6                   ,vpropel 
     7                   ,xl      ,ul         ,p0l      ,dpl   
     8                   ,pl      ,sl         ,ld       ,txnl
     9                   ,tx1l    ,tx2l       ,depsl    ,plasticl 
     1                   ,vpropell  
     2                   ,numel   ,nen        ,nenv     ,ndf 
     3                   ,ndm     ,nst        ,npi      ,ntn
     4                   ,neq     ,nequ       ,neqp 
     5                   ,nad     ,nadu       ,nadp     ,nadpu,nadr
     6                   ,lhs     ,rhs        ,unsym 
     7                   ,stge    ,isw        ,ilib     ,nlit
     8                   ,i_colorg,i_elcolor  ,numcolors
     9                   ,block_pu,n_blocks_pu,fplastic ,vprop)
c **********************************************************************
c * Data de criacao    : 12/12/2015                                    *
c * Data de modificaco : 11/04/2017                                    * 
c * ------------------------------------------------------------------ * 
c * PFORM_PM:                                                          *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ix(nen+1,numel) - conetividades nodais                             *
c * eloads(7,numel)   - cargas nos elementos mecanico                  *
c * eloadsp(7,numel)  - cargas nos elementos hidraulico                *
c * ie(numat)       - tipo de elemento                                 *
c * e(10,numat)     - constantes fisicas dos materiais                 *
c * x(ndm,nnode)    - coordenadas nodais                               *
c * id(ndf,nnode)   - numeracao global das equacoes                    *
c * ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro         *
c *             coeficiente nao-nulo da equacao i, no formato CSR      *
c * ja(nad)   - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor a, no formato CSR (stge = 1)      *
c *           - ponteiro da diagonal do skyline (stge = 4)             *
c * ad(neq)   - nao definido                                           *
c * au(nad)   - nao definido                                           *
c * al(nad)   - nao definido                                           *
c * b(neq)    - vetor de forcas nodais equivalentes                    *
c * u0(ndf,nnode) - solucao do passo de tempo anterior                 *
c * u(ndf,nnode) - solucao (com valores prescritos)                    *
c * tx1(7,npi,*) - tensoes                                             *
c *                  elastico : nao definido                           *
c *                  plastico : tensao nos pontos de integracao no     * 
c *                  passo de tempo anterior                           *
c * tx2(7,npi,*) - tensoes                                             *
c *                  elastico : nao definido                           *
c *                  plastico : tensao nos pontos de integracao        *
c * deps(7,npi,*) - elastico : nao definido                            *
c *               plastico :                                           *
c *                delta deformacao entre iteracoes nao lineares       *
c *                delta pressoes entre as iteracoes                   *
c * dp(*)      - delta p ( p(n  ,0  ) - p(0) )                         *  
c * plastic(3,*) - deformacao volumetricas plasticas no passo de tempo *
c *                anterior                                            *
c *                deformacao volumetricas plasticas                   *
c *                paramentro de endurecimento nos pontos de integracao*
c * elplastic(*) - identificao se o elemento plastificou ou nao(0 ou 1)*
c * vpropel(7,npi,*) - propriedades variaveis por pontos de itegracao  *
c *     1 - porosideade nos pontos de integracao                       *
c *     2 - permeabilidade konzey-Carman                               *
c *     3 - massa especifica                                           *
c *     4 - modulo volumetrico                                         *
c *     5 - modulo de cisalhamento                                     *
c *     6 - inverso do modulo de biot                                  *
c *     7 - coeficiente de biot                                        *
c * xl(ndm,nen)     - nao definido                                     *
c * ul(nst)         - nao definido                                     *
c * dpl(nst)        - nao definido                                     *
c * pl(nst)         - nao definido                                     *
c * sl(nst,nst)     - nao definido                                     *
c * ld(nst)         - nao definido                                     *
c * txnl(6    ,nen) - nao definido                                     *
c * tx1l(6,npi)     - nao definido                                     *
c * tx2l(6,npi)     - nao definido                                     *
c * despl(6,npi)    - nao definido                                     *
c * plastic(3,npi)  - nao definido                                     *
c * vpropell(7,npi) - nao definido                                     *
c * numel - numero de elementos                                        *
c * nen   - numero de nos por elemento                                 *
c * nenv  - numero de nos de vertice por elemento                      *
c * ndf   - numero max. de graus de liberdade por no                   *
c * ndm   - dimensao                                                   *
c * nst   - nen*ndf                                                    *
c * npi   - numero de pontos de integracao por elemento                *
c * ntn   - numero max. de derivadas por no                            *
c * neq   - numero de equacoes                                         *
c * nequ  - numero de equacoes em kuu                                  *   
c * nad   - numero de posicoes no CSR (storage = 1)                    *
c * nadu  - numero de posicoes no CSR (storage = 1)                    *
c * nadp  - numero de posicoes no CSR (storage = 1)                    *
c * nadpu - numero de posicoes no CSR (storage = 1,block_pu = true)    *
c * nadr  - numero de termos nao nulos na parte retangular             *
c *         ( MPI em overllaping )                                     *      
c * lhs   - flag para montagem de ad (letf hand side)                  *
c * rhs   - flag para correcao de b  (right hand side)                 *
c * unsym - flag para matrizes nao simetricas                          *
c * stge  - = 1, armazenamento CSR                                     *
c *         = 2, armazenamento por arestas                             *
c *         = 3, armazenamento EBE                                     *
c *         = 4, armazenamento SKYLINE                                 *
c *         = 0, nao monta matriz                                      *
c * isw   - codigo de instrucao para as rotinas de elemento            *
c * ilib  - determina a biblioteca do elemento                         *
c * nlit     - iteracao nao linear                                     * 
c * i_colorg(2,*)- inicio e final de cada cor no arranjo i_elcolor     *
c * i_elcolor(i) - elemento  agrupados por corda cor                   *
c * numcolors - numero de cores                                        *
c * block_pu    - true - armazenamento em blocos Kuu,Kpp e kpu         *
c *               false- aramzenamento em unico bloco                  *      
c * n_block_pu  - numeros de blocos                                    *
c * plastic  - (true/false)                                            * 
c * vprop(*) -                                                         * 
c *           1 - prop variavel                 (true|false)           *
c *           2 - konzey-Caraman                (true|false)           *
c *           3 - massa especifica homogenizada (true|false)           *
c *           4 - mecanico                      (true|false)           *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c *  b(neq) - vetor de forcas corrigido     (rhs = true)               *
c * ad(neq) - coeficientes da diagonal      (lhs = true)               *
c * au(nad) - parte triangular sup. de A (lhs = true)                  *
c * al(nad) - parte triangular inf. de A (lhs = true, unsym = true)    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'openmp.fi'
      include 'omp_lib.h'
      include 'transiente.fi'
      include 'termprop.fi'
      integer*8 nad,ia(*),aux1
      integer numel,nen,nenv,ndf,ndm,nst,nadu,nadp,nadpu,nadr
      integer stge,isw,numat,nlit,ntn,npi
      integer neq,nequ,neqp,n_blocks_pu
      integer ix(nen+1,*),eloads(7,*),eloadsp(7,*)
      integer ie(*),id(ndf,*),ld(nst)
      integer ja(*),elplastic(*)
      integer iel,ma,nel,no,i,j,k,kk,ilib
      integer i_colorg(2,*),i_elcolor(*),numcolors
      integer ic,jc
      real*8  e(prop,*),x(ndm,*),ad(*),au(*),al(*),b(*)
      real*8  u(ndf,*),u0(ndf,*),tx1(ntn,npi,*),tx2(ntn,npi,*),dp(*)
      real*8  deps(ntn+1,npi,*),plastic(3,npi,*),vpropel(nvprop,npi,*)
c ... varaiaveis locais por elemento
      real*8  xl(ndm,nenv),ul(nst),p0l(nenv),el(prop)
      real*8  pl(nst),sl(nst,nst),plasticl(3,npi),txnl(ntn,nen)
      real*8  tx1l(ntn,npi),tx2l(ntn,npi),depsl(ntn+1,npi),dpl(nen)
      real*8  vpropell(nvprop,npi)
      integer ep
c ... 
      logical lhs,rhs,unsym,block_pu,fplastic,vprop(*)
c .....................................................................
c
c ... Zera a matriz de coeficientes:
c
      if(lhs) then
        call dazero(au,nad)
        aux1 = nad+nadpu+nadr 
        call dazero(al,aux1)
        call azero(ad,neq)      
      endif
c ..................................................................... 
c
c ... openmp
      if(omp_elmt)then
c ... loop nos cores
c$omp parallel num_threads(nth_elmt)
c$omp.default(none)
c$omp.private(nel,ic,jc,no,k,ma,iel,i,j,kk,xl,ld,ul,pl,dpl,el,sl)
c$omp.private(p0l,tx1l,tx2l,vpropell,depsl,plasticl,txnl,ep)
c$omp.shared(numcolors,i_colorg,i_elcolor,nen,nenv,ndf,ndm)
c$omp.shared(eloads,eloadsp)
c$omp.shared(id,u,u0,ix,dp,tx1,tx2,ie,e,vpropel,block_pu,n_blocks_pu)
c$omp.shared(stge,unsym,rhs,lhs,ilib,nlit,ntn,npi,isw,nst,dt,elplastic)
c$omp.shared(neq,neqp,nequ,nad,nadu,nadp,nadpu,plastic,fplastic,vprop)
       do ic = 1, numcolors
c$omp do
c ... Loop nos elementos:
c     ------------------
         do jc = i_colorg(1,ic), i_colorg(2,ic)
           nel = i_elcolor(jc)
           kk = 0
c ... loop nos nos de deslocamentos
           do i = 1, nen
             no = ix(i,nel)
             do j = 1, ndf - 1
               kk     = kk + 1
               ld(kk) = id(j,no)
               pl(kk) = 0.d0
               ul(kk) = u(j,no)
             enddo
c ......................................................................              
           enddo   
c ......................................................................               
c
c...... loop nos nos do pressao
           do i = 1, nenv
             no       = ix(i,nel)
             kk       = kk + 1
             ld(kk)   = id(ndf,no)
             pl(kk)   = 0.d0
             ul(kk)   = u(ndf,no)
             p0l(i)   = u0(ndf,no)
             dpl(i)   = dp(no)
             do j = 1, ndm
               xl(j,i) = x(j,no)
             enddo   
           enddo    
c ......................................................................
c
c ... plasticidade
           if(fplastic) then
c ...
            do j = 1, npi              
c ... tensoes no passo de tempo anterior
              tx1l(1:ntn,j) = tx1(1:ntn,j,nel)
            enddo
c .....................................................................
c
c ...
            if(isw .eq. 2) then
c ...
              do j = 1, npi
c ... tensoes
                tx2l(1:ntn,j) = tx2(1:ntn,j,nel)
c ... delta deformacoes
                depsl(1:6,j)  = deps(1:6,j,nel)
c ... delta pressao                        
                depsl(7,j)    = deps(7,j,nel)
c ... dilatacao volumetrica plastica do passo de tempo anterior
                plasticl(1,j) = plastic(1,j,nel)
c ... dilatacao volumetrica plastica
                plasticl(2,j) = plastic(2,j,nel)
c ... paramentro de endurecimento 
                plasticl(3,j) = plastic(3,j,nel)
              enddo
            endif  
           endif  
c ......................................................................   
c
c ... propriedade variavel 
          if(vprop(1)) then
            do j = 1, npi
              vpropell(1:nvprop,j) = vpropel(1:nvprop,j,nel)  
            enddo
          endif
c ......................................................................     
c
c ...... Arranjos de elemento:
           ma  = ix(nen+1,nel)
           iel = ie(ma)
           el(1:prop) = e(1:prop,ma)
c ......................................................................
c
c ...... Chama biblioteca de elementos:
           call elmlib_pm(el   ,eloads(1,nel),eloadsp(1,nel)
     1                  ,xl    ,ul           ,p0l
     2                  ,dpl   ,pl           ,sl  ,tx1l    ,tx2l
     3                  ,depsl ,plasticl     ,vpropell,ep
     4                  ,ndm   ,nst          ,nel ,iel,isw
     5                  ,ma    ,nlit         ,ilib,block_pu) 
c ......................................................................
c
c ... plasticidade armazena nos arranjos globais
           if(fplastic .and. isw .eq. 2) then      
             if( ep .ne. 0 )  elplastic(nel) = ep 
c ...
             do j = 1, npi  
c ... tensoes
               tx2(1:ntn,j,nel) = tx2l(1:ntn,j) 
c ... delta deformacoes
               deps(1:6,j,nel)  = depsl(1:6,j)   
c ... delta pressao                        
               deps(7,j,nel)    = depsl(7,j)  
c ... dilatacao volumetrica plastica
               plastic(2,j,nel) = plasticl(2,j) 
c ... paramentro de endurecimento 
               plastic(3,j,nel) = plasticl(3,j)
             enddo
c ..................................................................... 
           endif  
c ...................................................................... 
c
c ...... Monta arrranjos locais em arranjos globais:
           call assbly_pm(sl   ,pl         ,ld
     1                ,ia      ,ja         ,au
     2                ,al      ,ad         ,b    ,nst
     3                ,neq     ,nequ       ,neqp
     4                ,nad     ,nadu       ,nadp ,nadpu
     5                ,lhs     ,rhs        ,unsym,stge
     6                ,block_pu,n_blocks_pu)
          enddo   
c$omp end do
c .....................................................................
        enddo
c$omp end parallel
c .....................................................................
c
c ... sequencial
      else  
c ... Loop nos elementos:
c     ------------------
        do nel = 1, numel
          kk = 0
c ... loop nos nos de deslocamentos
          do i = 1, nen
            no = ix(i,nel)
            do j = 1, ndf - 1
              kk     = kk + 1
              ld(kk) = id(j,no)
              pl(kk) = 0.d0
              ul(kk) = u(j,no)
            enddo
c ......................................................................
         enddo
c ......................................................................      
c
c...... loop nos nos do pressao
          do i = 1, nenv
            no       = ix(i,nel)
            kk       = kk + 1
            ld(kk)   = id(ndf,no)
            pl(kk)   = 0.d0
            ul(kk)   = u(ndf,no)
            p0l(i)   = u0(ndf,no)
            dpl(i)   = dp(no)
            do j = 1, ndm
              xl(j,i) = x(j,no)
            enddo   
          enddo
c .....................................................................
c
c ... plasticidade
          if(fplastic) then
c ...
            do j = 1, npi              
c ... tensoes no passo de tempo anterior
              tx1l(1:ntn,j) = tx1(1:ntn,j,nel)
            enddo
c .....................................................................
c
c ...
            if(isw .eq. 2) then
c ...
              do j = 1, npi
c ... tensoes
                tx2l(1:ntn,j) = tx2(1:ntn,j,nel)
c ... delta deformacoes
                depsl(1:6,j)  = deps(1:6,j,nel)
c ... delta pressao                        
                depsl(7,j)    = deps(7,j,nel)
c ... dilatacao volumetrica plastica do passo de tempo anterior
                plasticl(1,j) = plastic(1,j,nel)
c ... dilatacao volumetrica plastica
                plasticl(2,j) = plastic(2,j,nel)
c ... paramentro de endurecimento 
                plasticl(3,j) = plastic(3,j,nel)
              enddo
            endif  
          endif  
c ......................................................................
c
c ... propriedade variavel 
          if(vprop(1)) then
            do j = 1, npi
              vpropell(1:nvprop,j) = vpropel(1:nvprop,j,nel)  
            enddo
          endif
c ......................................................................     
c
c ...... Arranjos de elemento:
          ma         = ix(nen+1,nel)
          iel        = ie(ma)
          el(1:prop) = e(1:prop,ma)
c ......................................................................
c
c ...... Chama biblioteca de elementos:
          call elmlib_pm(el    ,eloads(1,nel)  ,eloadsp(1,nel)
     1                  ,xl    ,ul       ,p0l
     2                  ,dpl   ,pl       ,sl   ,tx1l    ,tx2l
     3                  ,depsl ,plasticl ,vpropell,ep
     4                  ,ndm   ,nst      ,nel  ,iel,isw
     5                  ,ma    ,nlit     ,ilib ,block_pu)        
c ......................................................................
c
c ... plasticidade armazena nos arranjos globais
          if(fplastic .and. isw .eq. 2) then  
            if( ep .ne. 0 )  elplastic(nel) = ep 
c ...
            do j = 1, npi  
c ... tensoes
              tx2(1:ntn,j,nel) = tx2l(1:ntn,j) 
c ... delta deformacoes
              deps(1:6,j,nel)  = depsl(1:6,j)   
c ... delta pressao                        
              deps(7,j,nel)    = depsl(7,j)  
c ... dilatacao volumetrica plastica
              plastic(2,j,nel) = plasticl(2,j) 
c ... paramentro de endurecimento 
              plastic(3,j,nel) = plasticl(3,j)
            enddo
c ..................................................................... 
          endif  
c ...................................................................... 
c
c ...... Monta arrranjos locais em arranjos globais:
          call assbly_pm(sl   ,pl         ,ld
     .               ,ia      ,ja         ,au
     .               ,al      ,ad         ,b    ,nst
     .               ,neq     ,nequ       ,neqp
     .               ,nad     ,nadu       ,nadp ,nadpu
     .               ,lhs     ,rhs        ,unsym,stge
     .               ,block_pu,n_blocks_pu)
c ...................................................................... 
        enddo 
c ......................................................................
      endif
c ......................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine tform_pm(ix    ,x    ,e    ,ie
     1                   ,ic    ,xl   ,ul   ,dpl     ,tx1l   
     2                   ,vpropell 
     3                   ,pl    ,u    ,dp   ,tx1     ,vpropel
     4                   ,t     ,tb   ,flux ,fnno
     5                   ,nnode ,numel,nen  ,nenv
     6                   ,ndm   ,ndf  ,nst  ,ntn  ,npi
     7                   ,isw   ,ilib ,i_xfi,novlp,plastic,vprop)
c **********************************************************************
c * Data de criacao    : 12/12/2015                                    *
c * Data de modificaco : 05/04/2017                                    * 
c * ------------------------------------------------------------------ * 
c * TFORM: caluclo das tensoes e dos fluxo nodal nos vertices          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ix(nen+1,numel)  - conetividades nodais                            *
c * x(ndm,nnode)     - coordenadas nodais                              *
c * e(10,numat)      - constantes fisicas dos materiais                *
c * ie(numat)        - tipo de elemento                                *
c * ic(nnode)        - nao definido                                    *
c * xl(ndm,nen)      - nao definido                                    *
c * ul(ndf,nen)      - nao definido                                    *
c * pl(ntn*nen)      - nao definido                                    *
c * dpl(nst)         - nao definido                                    *
c * tx1l(ntn,npi)    - nao definido                                    *
c * vpropell         - nao definido                                    *
c * u(ndf,nnode)     - solucao corrente                                *
c * dp(nnodev)       - delta p ( p(n  ,0  ) - p(0) )                   *
c * tx1              - elastico : nao definido                         *
c *                    plastico : tensao nos pontos de integracao      *
c * vpropel(7,npi,*) - propriedades variaveis por pontos de itegracao  *
c *     1 - porosideade nos pontos de integracao                       *
c *     2 - permeabilidade konzey-Carman                               *
c *     3 - massa especifica                                           *
c *     4 - modulo volumetrico                                         *
c *     5 - modulo de cisalhamento                                     *
c *     6 - inverso do modulo de biot                                  *
c *     7 - coeficiente de biot                                        *
c * t(ntn,nnodev)    - nao definido                                    *
c * tb(ntn,nnodev)   - nao definido                                    *
c * te(ntn,nnodev)   - nao definido                                    *
c * flux(ndm,nnodev) - nao definido                                    *
c * fnno             -  identifica dos nos de vertices                 *
c *                     ( 1 - vertice | 0 )                            *
c * nnodev           - numero de nos de vertices                       *
c * numel            - numero de elementos                             *
c * nen              - numero max. de nos por elemento                 *
c * nenv             - numero de nos de vertice por elemento           *
c * ndf              - numero max. de graus de liberdade por no        *
c * nst              - nst = nen*ndf                                   *
c * ndm              - dimensao                                        *
c * ndf              - numero max. de graus de liberdade por no        *
c * ntn              - numero max. de derivadas por no                 *
c * npi              - numero de pontos de integracao por elemento     *
c * nprcs            - numero de processos                             *
c * isw              - codigo de instrucao para as rotinas             *
c *                    de elemento                                     *
c * ilib             - determina a biblioteca do elemento              *
c * i_xf             - buffer de comunicao para o numero de elementos  *
c *                    conectidos aos nos                              *
c * novlp            - chave indicativa de non-overlapping             *
c * plastic  - (true/false)                                            *  
c * vprop(*) -                                                         * 
c *           1 - prop variavel                 (true|false)           *
c *           2 - konzey-Caraman                (true|false)           *
c *           3 - massa especifica homogenizada (true|false)           *
c *           4 - mecanico                      (true|false)           *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * t(ntn   ,nnodev) - tensoes medias nodais                           *
c * tb(ntn  ,nnodev) - tensoes efetivas biot medias nodais             *
c * te(ntn  ,nnodev) - tensoes efetivas medias nodais                  *
c * flux(ndm,nnodev) - fluxo medias nodais                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Calcula as tensoes apenas nos vertices                             *     
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'termprop.fi'
c ... mpi
      integer*8 i_xfi
      logical novlp
c ...
      integer nnode,numel,nen,nenv,ndf,nst,ndm,ntn,npi
c ......................................................................      
      integer ix(nen+1,*),ie(*),ic(*),fnno(*)
      integer nel,ma,iel,i,j,k,k1,no,kk
      integer ilib,isw,desloc1,desloc2
      real*8  xl(ndm,nenv),ul(nst),dpl(nenv),pl(nenv*(2*ntn+ndm))
      real*8  tx1l(ntn,npi),vpropell(nvprop,npi)
      real*8  x(ndm,*),e(prop,*),tx1(ntn,npi,*),vpropel(nvprop,npi,*)
      real*8  u(ndf,*),el(prop) ,dp(*) 
      real*8  t(ntn,*),tb(ntn,*),flux(ndm,*)
      logical plastic,vprop(*)
c ...
      logical ldum
      integer idum
      real*8 ddum
c ......................................................................
c
c ... deslocamentos do vetor local p
c     p(1:nenv*ntn)                    - tensao total 
c     p(nenv*ntn+1:2*nenv*ntn)         - tensao efetiva de biot 
c     p(2*nenv*ntn+1:nenv*(2*ntn+ndm)) - flux de darcy             
      desloc1 = nenv*ntn
      desloc2 = 2*nenv*ntn
c ......................................................................
c
c ...
      do 30 i = 1, nnode        
        ic(i) = 0
c ... tensao
        do 10 j = 1, ntn
           t(j,i) = 0.d0
          tb(j,i) = 0.d0
   10   continue 
c ......................................................................
c
c ... flux
        do 20 j = 1, ndm
          flux(j,i) = 0.d0
   20   continue 
c ......................................................................
   30 continue 
c ......................................................................
c
c ... Loop nos elementos:
      do 900 nel = 1, numel
        kk = 0
c ...
        do 300 i = 1, nenv*(2*ntn+ndm)
          pl(i) = 0.d0
  300   continue
c .......................................................................
c
c ... loop nos arranjos locais 
        do 400 i = 1, nen
          no = ix(i,nel)
c ... loop nos deslocamentos
          do 410 j = 1, ndf - 1
            kk     = kk + 1
            ul(kk) = u(j,no)
  410     continue
c ......................................................................
  400   continue
c ......................................................................
c
c ... loop nas pressoes
        do 510 i = 1, nenv
          no     = ix(i,nel)
          kk     = kk + 1
          ul(kk) = u(ndf,no)
          dpl(i) = dp(no)
          do 500 j = 1, ndm
            xl(j,i) = x(j,no)
  500     continue
  510   continue
c ......................................................................
c
c ... plasticidade
          if(plastic) then
            do j = 1, npi  
c ... tensoes
              tx1l(1:ntn,j)  = tx1(1:ntn,j,nel)
            enddo  
          endif  
c ......................................................................  
c
c ... propriedade variavel 
          if(vprop(1)) then
            do j = 1, npi
              vpropell(1:nvprop,j) = vpropel(1:nvprop,j,nel)
            enddo
          endif
c ......................................................................     
c
c ...... form element array
        ma  = ix(nen+1,nel)
        iel = ie(ma)      
        el(1:prop) = e(1:prop,ma)
c ......................................................................
c
c ...... Chama biblioteca de elementos:
        call elmlib_pm(el  ,idum,idum
     1                ,xl  ,ul  ,ddum
     2                ,dpl ,pl  ,ddum,tx1l,tx1l 
     3                ,ddum,ddum,vpropel  ,idum
     4                ,ndm ,nst ,nel,iel,isw 
     5                ,ma  ,idum,ilib,ldum)
c ......................................................................
c
c ...... media do vetor global
        do 800 i = 1, nenv
           no = ix(i,nel)
           if (no .le. 0) go to 800
           ic(no) = ic(no) + 1  
           do 700 j = 1, ntn
c ... tensao total
              k       = (i-1)*ntn + j
              t(j,no) = t(j,no)  + pl(k)
c ... tensao efetiva de biot 
              k1      = (i-1)*ntn + j + desloc1
              tb(j,no)= tb(j,no) + pl(k1)
  700      continue
           do 710 j = 1, ndm
              k = (i-1)*ndm + j + desloc2
              flux(j,no) = flux(j,no) + pl(k)
  710      continue
  800   continue
  900 continue
c .......................................................................
c
c ... Comunica vetor de contagem de elementos por no'
      if (novlp) call allgatheri(ic,i_xfi)
c .......................................................................
      do 1000 i = 1, nnode
c ... no de vertice
        if(fnno(i) .eq. 1) then
c ... tensao
          do 1010 j = 1, ntn
c ... tensao total            
            t(j,i)  = t(j,i)/ic(i) 
c ... tensao biot           
            tb(j,i) = tb(j,i)/ic(i)
 1010     continue 
c ......................................................................
c
c ... fluxo
          do 1030 j = 1, ndm
            flux(j,i) = flux(j,i)/ic(i)
 1030     continue 
c ......................................................................
        endif
 1000 continue
c ......................................................................
c
c ......................................................................
      return
      end      
c **********************************************************************
c
c **********************************************************************
      subroutine initial_stress(ix ,ie       ,e
     1                   ,x        ,id       ,b
     2                   ,u        ,tx       ,tx0      ,dp
     3                   ,xl       ,ul       ,dpl      ,pl            
     4                   ,ld       ,txnl     ,tx1l
     5                   ,numel    ,nen      ,nenv     ,ndf 
     6                   ,ndm      ,nst      ,npi      ,ntn       
     7                   ,neq      ,stge     ,ilib     
     8                   ,block_pu ,plastic)
c **********************************************************************
c * Data de criacao    : 26/12/2015                                    *
c * Data de modificaco : 11/05/2017                                    * 
c * ------------------------------------------------------------------ * 
c * INITIAL_STRESS: calculo da influencia das tensoes iniciais         *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ix(nen+1,numel) - conetividades nodais                             *
c * ie(numat)       - tipo de elemento                                 *
c * e(10,numat)     - constantes fisicas dos materiais                 *
c * x(ndm,nnode)    - coordenadas nodais                               *
c * id(ndf,nnode)   - numeracao global das equacoes                    *
c * b(neq)    - vetor de forcas nodais equivalentes                    *
c * u(ndf,nnode) - solucao (com valores prescritos)                    *
c * tx           - elastico : nao definido                             *
c *                plastico : tensao nos pontos de integracao          *
c * tx0         - tensao inicial nodal                                 *
c * dp(*)       - delta p ( p(n  ,0  ) - p(0) )                        *  
c * xl(ndm,nen)     - nao definido                                     *
c * ul(nst)         - nao definido                                     *
c * dpl(nst)        - nao definido                                     *
c * pl(nst)         - nao definido                                     *
c * ld(nst)         - nao definido                                     *
c * txnl(ntn  ,nen) - nao definido                                     *
c * tx1l(6,npi)     - nao definido                                     *
c * numel - numero de elementos                                        *
c * nen   - numero de nos por elemento                                 *
c * nenv  - numero de nos de vertice por elemento                      *
c * ndf   - numero max. de graus de liberdade por no                   *
c * ndm   - dimensao                                                   *
c * nst   - nen*ndf                                                    *
c * npi   - numero de pontos de integracao por elemento                *
c * stge  - = 1, armazenamento CSR                                     *
c *         = 2, armazenamento por arestas                             *
c *         = 3, armazenamento EBE                                     *
c *         = 4, armazenamento SKYLINE                                 *
c *         = 0, nao monta matriz                                      *
c * isw   - codigo de instrucao para as rotinas de elemento            *
c * ilib  - determina a biblioteca do elemento                         *
c *               false- aramzenamento em unico bloco                  *      
c * block_pu    - true - armazenamento em blocos Kuu,Kpp e kpu         *
c *               false- aramzenamento em unico bloco                  *      
c * n_block_pu  - numeros de blocos                                    *
c * plastic  - (true/false)                                            * 
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * b(neq) - vetor de forcas devido as tensoes inicias (elastic)       *
c * tx      - tensao nos pontos de integracao com as tensoes iniciais  *
c *           (plastic)                                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'openmp.fi'
      include 'omp_lib.h'
      include 'transiente.fi'
      include 'termprop.fi'
      integer idum,ddum
      integer numel,nen,nenv,ntn,ndf,ndm,nst
      integer stge,numat,npi
      integer neq
      integer ix(nen+1,*),ie(*),id(ndf,*),ld(nst)
      integer iel,ma,nel,no,i,j,k,kk,ilib
      integer ic,jc
      real*8  e(prop,*),x(ndm,*),b(*)
      real*8  u(ndf,*),tx0(ntn,*),tx(ntn,npi,*),dp(*)
      real*8  xl(ndm,nenv),ul(nst),el(prop)
      real*8  pl(nst),txnl(ntn,nen),tx1l(ntn,npi),dpl(nen)
      logical block_pu,plastic
c .....................................................................
c
c ... Loop nos elementos:
c     ------------------
      do nel = 1, numel
        kk = 0
c ... loop nos nos de deslocamentos
        do i = 1, nen
          no = ix(i,nel)
          do j = 1, ndf - 1
            kk     = kk + 1
            ld(kk) = id(j,no)
            pl(kk) = 0.d0
            ul(kk) = u(j,no)
          enddo
c ......................................................................
        enddo
c ......................................................................      
c
c...... loop nos nos do pressao
        do i = 1, nenv
          no       = ix(i,nel)
          kk       = kk + 1
          ld(kk)   = id(ndf,no)
          pl(kk)   = 0.d0
          ul(kk)   = u(ndf,no)
          dpl(i)   = dp(no)
          do j = 1, ndm
            xl(j,i) = x(j,no)
          enddo   
        enddo
c .....................................................................
c            
c ... tensao inicial
        do i = 1, nen  
          no          = ix(i,nel)  
          txnl(1:ntn,i) = tx0(1:ntn,no)
        enddo  
c ......................................................................     
c
c ...... Arranjos de elemento:
        ma  = ix(nen+1,nel)
        iel = ie(ma)
        el(1:prop) = e(1:prop,ma)
c ......................................................................
c
c ...... Chama biblioteca de elementos:
        call elmlib_pm(el  ,idum     ,idum  
     1                ,xl  ,ul       ,ddum
     2                ,dpl ,pl       ,ddum ,txnl,tx1l  
     3                ,ddum,ddum     ,ddum ,idum
     4                ,ndm ,nst ,nel ,iel ,5
     5                ,ma  ,idum     ,ilib,block_pu)        
c ......................................................................

c ... atualiza as tensoes no pontos de integracao: plastico
        if(plastic) then
c ...
          do j = 1, npi  
c ... tensoes
            tx(1:ntn,j,nel) = tx1l(1:ntn,j) 
          enddo
c .....................................................................
c
c ...... Monta o vetor b global: elastico
        else
          do j = 1, nst
            k = ld(j)
            if (k .gt. 0 .and. k .le. neq) b(k) = b(k) - pl(j)
          enddo
        endif
c ...................................................................... 
      enddo 
c ......................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine initial_pc0(ix    ,ie       ,e
     1                   ,x        ,tx0      ,u        ,plastic
     3                   ,xl       ,ul       ,txnl     ,plasticl   
     5                   ,numel    ,nen      ,nenv     ,ndf 
     6                   ,ndm      ,nst      ,npi      ,ntn       
     7                   ,ilib)
c **********************************************************************
c * Data de criacao    : 24/02/2017                                    *
c * Data de modificaco : 11/05/2017                                    * 
c * ------------------------------------------------------------------ * 
c * INITIAL_PC0: calculo da influencia das tensoes iniciais            *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ix(nen+1,numel) - conetividades nodais                             *
c * ie(numat)       - tipo de elemento                                 *
c * e(10,numat)     - constantes fisicas dos materiais                 *
c * x(ndm,nnode)    - coordenadas nodais                               *
c * u(ndf,nnode) - solucao (com valores prescritos)                    *
c * tx0         - tensao inicial nodal                                 *
c * plastic(3,*) - deformacao volumetricas plasticas no passo de tempo *
c *                anterior                                            *
c *                deformacao volumetricas plasticas                   *
c *                paramentro de endurecimento nos pontos de integracao*
c * xl(ndm,nen)     - nao definido                                     *
c * ul(nst)         - nao definido                                     *
c * txnl(ntn  ,nen) - nao definido                                     *
c * plastic(3,npi)  - nao definido                                     *
c * numel - numero de elementos                                        *
c * nen   - numero de nos por elemento                                 *
c * nenv  - numero de nos de vertice por elemento                      *
c * ndf   - numero max. de graus de liberdade por no                   *
c * ndm   - dimensao                                                   *
c * nst   - nen*ndf                                                    *
c * npi   - numero de pontos de integracao por elemento                *
c * npi   - numero de pontos de integracao por elemento                *
c * ntn   - numero max. de tensoes por no                              *
c * ilib  - determina a biblioteca do elemento                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * plastic(3,j,nel) - inicializado                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'openmp.fi'
      include 'omp_lib.h'
      include 'transiente.fi'
      include 'termprop.fi'
      integer idum,ddum
      integer numel,nen,nenv,ntn,ndf,ndm,nst
      integer stge,numat,npi
      integer neq
      integer ix(nen+1,*),ie(*)
      integer iel,ma,nel,no,i,j,k,kk,ilib
      integer ic,jc
      real*8  e(prop,*),x(ndm,*)
      real*8  u(ndf,*),tx0(ntn,*),plastic(3,npi,*)
      real*8  xl(ndm,nenv),ul(nst),el(prop)
      real*8  pl(nst),txnl(ntn,nen),plasticl(3,npi)
      logical ldum
c .....................................................................
c
c ... Loop nos elementos:
c     ------------------
      do nel = 1, numel
        kk = 0
c ... loop nos nos de deslocamentos
        do i = 1, nen
          no = ix(i,nel)
          do j = 1, ndf - 1
            kk     = kk + 1
            ul(kk) = u(j,no)
          enddo
c ......................................................................
        enddo
c ......................................................................      
c
c...... loop nos nos do pressao
        do i = 1, nenv
          no       = ix(i,nel)
          kk       = kk + 1
          ul(kk)   = u(ndf,no)
          do j = 1, ndm
            xl(j,i) = x(j,no)
          enddo   
        enddo
c .....................................................................
c            
c ... tensao inicial
        do i = 1, nen  
          no            = ix(i,nel)  
          txnl(1:ntn,i) = tx0(1:ntn,no)
        enddo  
c ......................................................................     
c
c ...... Arranjos de elemento:
        ma  = ix(nen+1,nel)
        iel = ie(ma)
        el(1:prop) = e(1:prop,ma)
c ......................................................................
c
c ...... Chama biblioteca de elementos:
        call elmlib_pm(el  ,idum     ,idum 
     1                ,xl  ,ul       ,ddum
     2                ,ddum,ddum     ,ddum ,txnl,ddum  
     3                ,ddum,plasticl ,ddum ,idum
     4                ,ndm ,nst ,nel ,iel  ,9
     5                ,ma  ,idum     ,ilib ,ldum)        
c ......................................................................

c ... atualiza as tensoes no pontos de integracao: plastico
        do j = 1, npi  
c ... paramentro de endurecimento 
          plastic(3,j,nel) = plasticl(3,j)
        enddo
c ...................................................................... 
      enddo 
c ......................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine porosity_form(ix ,x       ,e      ,ie   
     1                   ,u       ,porosity,plastic,dporo 
     2                   ,vpropel ,ic      ,fnno 
     3                   ,xl      ,ul      ,pl,plasticl,vpropell 
     4                   ,nnode   ,numel   ,nen    ,nenv
     5                   ,ndm     ,ndf     ,nst    ,npi 
     6                   ,isw     ,ilib    ,i_xfi  ,novlp
     7                   ,fplastic,vprop)
c **********************************************************************
c * Data de criacao    : 26/09/2016                                    *
c * Data de modificaco : 11/05/2017                                    * 
c * ------------------------------------------------------------------ * 
c * POROSITY : calculo da porosidade                                   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ix(nen+1,numel)  - conetividades nodais                            *
c * x(ndm,nnode)     - coordenadas nodais                              *
c * e(10,numat)      - constantes fisicas dos materiais                *
c * ie(numat)        - tipo de elemento                                *
c * u(ndf,nnode)     - solucao corrente                                *
c * porosity(nnode)  - nao definido                                    *
c * plastic(3,*) - deformacao volumetricas plasticas no passo de tempo *
c *                anterior                                            *
c *                deformacao volumetricas plasticas                   *
c *                paramentro de endurecimento nos pontos de integracao*
c * dporo(nnode)     - nao definido                                    *
c * fnno             -  identifica dos nos de vertices                 *
c *                     ( 1 - vertice | 0 )                            *
c * ic(nnode)        - nao definido                                    *
c * vpropel(7,npi,*) - propriedades variaveis por pontos de itegracao  *
c *     1 - porosideade nos pontos de integracao                       *
c *     2 - permeabilidade konzey-Carman                               *
c *     3 - massa especifica                                           *
c *     4 - modulo volumetrico                                         *
c *     5 - modulo de cisalhamento                                     *
c *     6 - inverso do modulo de biot                                  *
c *     7 - coeficiente de biot                                        *
c * xl(ndm,nen)      - nao definido                                    *
c * ul(ndf,nen)      - nao definido                                    *
c * pl(ntn*nen)      - nao definido                                    *
c * plastic(3,npi)  - nao definido                                     *
c * vpropell(5,npi) - nao definido                                     *
c * nnodev           - numero de nos de vertices                       *
c * numel            - numero de elementos                             *
c * nen              - numero max. de nos por elemento                 *
c * nenv             - numero de nos de vertice por elemento           *
c * ndf              - numero max. de graus de liberdade por no        *
c * nst              - nst = nen*ndf                                   *
c * npi              - numero de pontos de integracao por elemento     *
c * ndm              - dimensao                                        *
c * ndf              - numero max. de graus de liberdade por no        *
c * ntn              - numero max. de derivadas por no                 *
c * nprcs            - numero de processos                             *
c * isw              - codigo de instrucao para as rotinas             *
c *                    de elemento                                     *
c * ilib             - determina a biblioteca do elemento              *
c * i_xf             - buffer de comunicao para o numero de elementos  *
c *                    conectidos aos nos                              *
c * novlp            - chave indicativa de non-overlapping             * 
c * vprop(*) -                                                         * 
c *           1 - prop variavel                 (true|false)           *
c *           2 - konzey-Caraman                (true|false)           *
c *           3 - massa especifica homogenizada (true|false)           *
c *           4 - mecanico                      (true|false)           *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * pososity(nnode)  - porosidade atualizada                           *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Calcula as porosidades apenas nos vertices                         *     
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'termprop.fi'
c ... mpi
      integer*8 i_xfi
      logical novlp
c ...
      integer nnode,numel,nen,nenv,ndf,nst,ndm,ntn,npi
c ......................................................................      
      integer ix(nen+1,*),ie(*),ic(*),fnno(*)
      integer nel,ma,iel,i,j,k,k1,no,kk
      integer ilib,isw
      real*8 xl(ndm,nenv),ul(nst),dpl(nenv),pl(nenv),plasticl(3,npi)
      real*8 vpropell(nvprop,*)
      real*8 x(ndm,*),e(prop,*),plastic(3,npi,*),vpropel(nvprop,npi,*)
      real*8 u(ndf,*),el(prop),porosity(*),dporo(*)
c ...
      logical ldum
      integer idum
      real*8 ddum
c ...
      logical fplastic,vprop(*)
c ......................................................................
c
c ...
      do 30 i = 1, nnode
        ic(i)    = 0
        dporo(i) = 0.d0
   30 continue
c ......................................................................
c
c ... Loop nos elementos:
      do 900 nel = 1, numel
        kk = 0
c ... 
        pl(1:nenv) = 0.d0
c .......................................................................
c
c ... loop nos arranjos locais 
        do 400 i = 1, nen
          no = ix(i,nel)
c ... loop nos deslocamentos
          do 410 j = 1, ndf - 1
            kk     = kk + 1
            ul(kk) = u(j,no)
  410     continue
c ......................................................................
  400   continue
c ......................................................................
c
c ... loop nas pressoes
        do 510 i = 1, nenv
          no     = ix(i,nel)
          kk     = kk + 1
          ul(kk) = u(ndf,no)
          do 500 j = 1, ndm
            xl(j,i) = x(j,no)
  500     continue
  510   continue
c ......................................................................
c
c ... plasticidade
        if(fplastic) then
          do j = 1, npi  
c ... dilatacao volumetrica plastica do passo de tempo anterior
            plasticl(1,j) = plastic(1,j,nel)
c ... dilatacao volumetrica plastica
            plasticl(2,j) = plastic(2,j,nel)
          enddo  
        endif  
c ...................................................................... 
c
c ...
        if(vprop(1)) then
          do j = 1, npi  
c ... porosidade nos pontos de inegracao
            vpropell(1,j) = vpropel(1,j,nel)
          enddo  
        endif
c ...................................................................... 
c
c ...... Arranjos de elemento:
        ma  = ix(nen+1,nel)
        iel = ie(ma)
        el(1:prop) = e(1:prop,ma)
c ......................................................................
c
c
c ...... Chama biblioteca de elementos:
        call elmlib_pm(el  ,idum    ,idum 
     1                ,xl  ,ul      ,ddum
     2                ,ddum,pl      ,ddum     ,ddum,ddum
     3                ,ddum,plasticl,vpropell ,idum
     4                ,ndm ,nst     ,nel      ,iel ,isw
     5                ,ma  ,idum    ,ilib     ,ldum)
c ......................................................................
c
c ...... media do vetor global
        do 800 i = 1, nenv
           no = ix(i,nel)
           if (no .le. 0) go to 800
           ic(no)    = ic(no) + 1
           dporo(no) = dporo(no) + pl(i)
  800   continue
  900 continue
c .......................................................................
c
c ... Comunica vetor de contagem de elementos por no'
      if (novlp) call allgatheri(ic,i_xfi)
c .......................................................................
c
c ...
      if (vprop(1)) then
        do i = 1, nnode
          if(fnno(i) .eq. 1 ) then
            porosity(i) = dporo(i)/ic(i)
          endif
        enddo
      else
        do i = 1, nnode
          if(fnno(i) .eq. 1 ) then
            porosity(i) = porosity(i) + dporo(i)/ic(i)
          endif
        enddo
      endif   
c ......................................................................
c
c ......................................................................
      return
      end      
c **********************************************************************
c
c **********************************************************************
      subroutine update_prop(ix   ,x     ,e   ,ie      ,vpropel
     1                      ,u    ,plastic
     2                      ,xl    ,ul   ,plasticl,vpropell  
     3                      ,numel,nen   ,nenv
     4                      ,ndm   ,ndf  ,nst     ,npi 
     5                      ,isw   ,ilib ,fplastic,vprop ,up_porosity)
c **********************************************************************
c * Data de criacao    : 26/09/2016                                    *
c * Data de modificaco : 11/05/2017                                    * 
c * ------------------------------------------------------------------ * 
c * POROSITY : calculo da porosidade                                   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ix(nen+1,numel)  - conetividades nodais                            *
c * x(ndm,nnode)     - coordenadas nodais                              *
c * e(10,numat)      - constantes fisicas dos materiais                *
c * ie(numat)        - tipo de elemento                                *
c * vpropel(7,npi,*) - propriedades variaveis por pontos de itegracao  *
c *     1 - porosideade nos pontos de integracao                       *
c *     2 - permeabilidade konzey-Carman                               *
c *     3 - massa especifica                                           *
c *     4 - modulo volumetrico                                         *
c *     5 - modulo de cisalhamento                                     *
c *     6 - inverso do modulo de biot                                  *
c *     7 - coeficiente de biot                                        *
c * u(ndf,nnode)     - delta u e p                                     *
c * plastic(3,*) - deformacao volumetricas plasticas no passo de tempo *
c *                anterior                                            *
c *                deformacao volumetricas plasticas                   *
c *                paramentro de endurecimento nos pontos de integracao*      
c * xl(ndm,nen)      - nao definido                                    *
c * ul(ndf,nen)      - nao definido                                    *
c * plastic(3,npi)  - nao definido                                     *
c * vpropell(7,npi)  - nao definido                                    *
c * numel            - numero de elementos                             *
c * nen              - numero max. de nos por elemento                 *
c * nenv             - numero de nos de vertice por elemento           *
c * ndm              - dimensao                                        *
c * ndf              - numero max. de graus de liberdade por no        *
c * isw              - codigo de instrucao para as rotinas             *
c * ilib             - determina a biblioteca do elemento              *
c * fplatic          - (true|false)                                    *            
c * vprop(*) -                                                         * 
c *           1 - prop variavel                 (true|false)           *
c *           2 - konzey-Caraman                (true|false)           *
c *           3 - massa especifica homogenizada (true|false)           *
c *           4 - mecanico                      (true|false)           *
c * up_porosity - atualizada o porosidade para o proximo passo de tempo*
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * vpropel(7,npi,*) - atualizados                                     *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'termprop.fi'
c ... mpi
      integer*8 i_xfi
      logical novlp
c ...
      integer nnode,numel,nen,nenv,ndf,nst,ndm,npi
c ......................................................................      
      integer ix(nen+1,*),ie(*)
      integer nel,ma,iel,i,j,k,k1,no,kk
      integer ilib,isw
      real*8  xl(ndm,nenv),ul(nst),vpropell(nvprop,*),plasticl(3,npi)
      real*8  x(ndm,*),e(prop,*),vpropel(nvprop,npi,*),plastic(3,npi,*)
      real*8  u(ndf,*),el(prop)
c ...
      logical ldum
      integer idum
      real*8 ddum
c ...
      logical fplastic,vprop(*),up_porosity
c ......................................................................
c
c ... Loop nos elementos:
      do 900 nel = 1, numel
        kk = 0
c ... loop nos arranjos locais 
        do 400 i = 1, nen
          no = ix(i,nel)
c ... loop nos deslocamentos
          do 410 j = 1, ndf - 1
            kk     = kk + 1
            ul(kk) = u(j,no)
  410     continue
c ......................................................................
  400   continue
c ......................................................................
c
c ... loop nas pressoes
        do 510 i = 1, nenv
          no     = ix(i,nel)
          kk     = kk + 1
          ul(kk) = u(ndf,no)
          do 500 j = 1, ndm
            xl(j,i) = x(j,no)
  500     continue
  510   continue
c ......................................................................
c
c
c ... plasticidade
        if(fplastic) then
          do j = 1, npi  
c ... dilatacao volumetrica plastica do passo de tempo anterior
            plasticl(1,j) = plastic(1,j,nel)
c ... dilatacao volumetrica plastica
            plasticl(2,j) = plastic(2,j,nel)
          enddo  
        endif  
c ...................................................................... 
c
c ... propriedade variavel 
        do j = 1, npi
          vpropell(1:nvprop,j) = vpropel(1:nvprop,j,nel)             
        enddo
c ......................................................................     
c
c ...... Arranjos de elemento:
        ma  = ix(nen+1,nel)
        iel = ie(ma)
        el(1:prop) = e(1:prop,ma)
c ......................................................................
c
c
c ...... Chama biblioteca de elementos:
        call elmlib_pm(el  ,idum    ,idum 
     1                ,xl  ,ul      ,ddum
     2                ,ddum,ddum    ,ddum,ddum,ddum
     3                ,ddum,plasticl,vpropell,idum
     4                ,ndm ,nst     ,nel ,iel ,isw
     5                ,ma  ,idum    ,ilib,ldum)
c ......................................................................
c
c ... propriedade variavel 
        do j = 1, npi
c ... porosidade
          if(up_porosity)vpropel(1,j,nel) = vpropell(1,j) 
c ... permeabilidade
          if(vprop(2)) vpropel(2,j,nel) = vpropell(2,j)
c ... massa especifica
          if(vprop(3)) vpropel(3,j,nel) = vpropell(3,j) 
c ... propriedades mecanicas
          if(vprop(4)) then
            vpropel(4,j,nel) = vpropell(4,j)
            vpropel(5,j,nel) = vpropell(5,j)
            vpropel(6,j,nel) = vpropell(6,j)
            vpropel(7,j,nel) = vpropell(7,j)
          endif 
c .....................................................................           
        enddo
c ......................................................................     
  900 continue
c .......................................................................
c
c ......................................................................
      return
      end      
c **********************************************************************
c
c **********************************************************************
      subroutine initial_porosity(ix    ,e    ,ie   ,ic
     1                         ,porosity,fnno
     2                         ,nnode   ,numel,nen     ,nenv
     3                         ,ndm     ,ndf  ,i_xfi   ,novlp)
c **********************************************************************
c * Data de criacao    : 26/09/2016                                    *
c * Data de modificaco : 11/05/2017                                    * 
c * ------------------------------------------------------------------ * 
c * DELTA_PORISITY: calculo da porosidade                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ix(nen+1,numel)  - conetividades nodais                            *
c * e(10,numat)      - constantes fisicas dos materiais                *
c * ie(numat)        - tipo de elemento                                *
c * ic(nnode)        - nao definido                                    *
c * pososity(nnodev) - nao definido                                    *
c * fnno             -  identifica dos nos de vertices                 *
c *                     ( 1 - vertice | 0 )                            *
c * nnodev           - numero de nos de vertices                       *
c * numel            - numero de elementos                             *
c * nen              - numero max. de nos por elemento                 *
c * nenv             - numero de nos de vertice por elemento           *
c * ndf              - numero max. de graus de liberdade por no        *
c * ndm              - dimensao                                        *
c * ndf              - numero max. de graus de liberdade por no        *
c * ilib             - determina a biblioteca do elemento              *
c * i_xf             - buffer de comunicao para o numero de elementos  *
c *                    conectidos aos nos                              *
c * novlp            - chave indicativa de non-overlapping             * 
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * pososity(nnodev) - porosidade                                      *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Calcula as tensoes apenas nos vertices                             *     
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'termprop.fi'
c ... mpi
      integer*8 i_xfi
      logical novlp
c ...
      integer nnode,numel,nen,nenv,ndf,ndm
c ......................................................................      
      integer ix(nen+1,*),ie(*),ic(*),fnno(*)
      integer nel,ma,i,no
      real*8  e(prop,*),porosity(*),poro
c ......................................................................
c
c ...
      do 30 i = 1, nnode
        ic(i)        = 0
        porosity(i)  = 0.d0
   30 continue
c ......................................................................
c
c ... Loop nos elementos:
      do 900 nel = 1, numel
c ...... Arranjos de elemento:
        ma         = ix(nen+1,nel)
        poro       = e(8,ma)
c ......................................................................
c
c ...... media do vetor global
        do 800 i = 1, nenv
           no = ix(i,nel)
           if (no .le. 0) go to 800
           ic(no)        = ic(no) + 1
           porosity(no)  = porosity(no) + poro
  800   continue
  900 continue
c .......................................................................
c
c ... Comunica vetor de contagem de elementos por no'
      if (novlp) call allgatheri(ic,i_xfi)
c .......................................................................
      do 1000 i = 1, nnode
        if(fnno(i) .eq. 1 ) then
          porosity(i) = porosity(i)/ic(i)
        endif
 1000 continue
c ......................................................................
      return
      end      
c **********************************************************************
c
c **********************************************************************
      subroutine initial_prop(ix    ,e      ,ie   ,vpropel   
     1                       ,numel,nen     ,numat,npi
     2                       ,vprop)
c **********************************************************************
c * Data de criacao    : 28/02/2017                                    *
c * Data de modificaco : 29/05/2017                                    * 
c * ------------------------------------------------------------------ * 
c * INITIAL_PROP: calculo da prorpriedaes iniciais                     *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ix(nen+1,numel)  - conetividades nodais                            *
c * e(10,numat)      - constantes fisicas dos materiais                *
c * ie(numat)        - tipo de elemento                                *
c * vpropel          - nao definido                                    *
c * numel            - numero de elementos                             *
c * nen              - numero max. de nos por elemento                 *
c * numat            - numero de elementos                             * 
c * npi              - numero de pontos de integracao por elemento     *
c * vprop(*) -                                                         * 
c *           1 - prop variavel                 (true|false)           *
c *           2 - konzey-Caraman                (true|false)           *
c *           3 - massa especifica homogenizada (true|false)           *
c *           4 - mecanico                      (true|false)           *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * vprop            - vetor inicializado                              *
c *     1 - porosideade nos pontos de integracao                       *
c *     2 - permeabilidade konzey-Carman                               *
c *     3 - massa especifica homogenizada                              *
c *     4 - modulo volumetrico                                         *
c *     5 - modulo de cisalhamento                                     *
c *     6 - modulo de biot                                             *
c *     7 - coeficiente de biot                                        *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'termprop.fi'
c ... mpi
      integer*8 i_xfi
      logical novlp
c ...
      integer numel,nen,npi
c ......................................................................      
      integer ix(nen+1,*),ie(*)
      integer nel,ma,i,numat
      real*8  e(prop,*),vpropel(nvprop,npi,*)
      real*8  perm,poro,ym,ps,ibiot,cbiot,ro,rof,lam,mu,k
      real*8 scale
      parameter (scale = 1.d-06)
      logical vprop(*)
c ......................................................................
c
c ... calculo das propriedades da fase solida
      do 100 ma = 1, numat
c ...
        ym   = e(1,ma)
        ps   = e(2,ma)
        ro   = e(6,ma)
        rof  = e(7,ma)
        poro = e(8,ma)
c .....................................................................
c
c ...
        lam  = ym*ps/( (1.0+ps)*( 1.0 - 2.0*ps ) )
        mu   = (0.5d0*ym)/(1.0+ps)
        k    = lam + (2.d0/3.d0)*mu 
c ......................................................................
c
c ... massa especifica da parte solida
        if( vprop(3) ) e(6,ma) = (ro - poro*rof)/(1.d0 - poro)
c .....................................................................
c
c ... calculo modulo volumetrico e de cisalhamento da parte solida
        if( vprop(4) ) call ks_and_mus(ym,ps,poro,k,mu)
        e(1,ma) = k
        e(2,ma) = mu
c .....................................................................
  100 continue
c ......................................................................
c
c ... Loop nos elementos:
      do 300 nel = 1, numel
c ...... Arranjos de elemento:
        ma         = ix(nen+1,nel)
c ...
        k          = e(1,ma)
        mu         = e(2,ma)
        perm       = e(3,ma)
        ro         = e(6,ma)
        rof        = e(7,ma)
        poro       = e(8,ma)
c ... modulo de volumetrico         
        ibiot      = 1/e(4,ma)
c ......................................................................
c
c ... modulo de volumetrico         
        cbiot       = e(5,ma)
c ......................................................................
c
c ...
        do 200 i = 1, npi
c ... porosidade
          vpropel(1,i,nel) = poro
c ......................................................................
c
c ... permebilidade
          vpropel(2,i,nel) = perm
c ......................................................................
c
c ... massa especifica homogenizada
          if( vprop(3) ) then
            vpropel(3,i,nel) = (1.0 - poro)*ro +  poro*rof
          else
            vpropel(3,i,nel) = ro
          endif
          vpropel(3,i,nel) = vpropel(3,i,nel)*scale 
c ......................................................................
          
c
c ... hashin-shtrikman
          if( vprop(4) ) then
            call hashin_shtrikman(poro,k,mu
     .                         ,vpropel(4,i,nel),vpropel(5,i,nel)
     .                         ,vpropel(6,i,nel),vpropel(7,i,nel))  
          else
            vpropel(4,i,nel)  = k 
            vpropel(5,i,nel)  = mu
            vpropel(6,i,nel)  = ibiot 
            vpropel(7,i,nel)  = cbiot
          endif
c ......................................................................
c
  200   continue
c ......................................................................
  300 continue
c ......................................................................
      return
      end      
c **********************************************************************
c
c **********************************************************************
      subroutine consolidation_pressure(ix  ,ie        ,ic   
     1                         ,plasticl,pl   
     2                         ,pc      ,plastic ,fnno
     3                         ,nnode   ,numel,nen     ,nenv
     4                         ,ndm     ,ndf  ,nst     ,npi
     5                         ,isw     ,ilib ,i_xfi   ,novlp)
c **********************************************************************
c * Data de criacao    : 24/01/2017                                    *
c * Data de modificaco : 11/05/2017                                    * 
c * ------------------------------------------------------------------ * 
c * CONSILADATION_PRESURE:  da porosidade                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ix(nen+1,numel)  - conetividades nodais                            *
c * x(ndm,nnode)     - coordenadas nodais                              *
c * e(10,numat)      - constantes fisicas dos materiais                *
c * ie(numat)        - tipo de elemento                                *
c * ic(nnode)        - nao definido                                    *
c * xl(ndm,nen)      - nao definido                                    *
c * ul(ndf,nen)      - nao definido                                    *
c * pl(ntn*nen)      - nao definido                                    *
c * dpl(nst)         - nao definido                                    *
c * u(ndf,nnode)     - solucao corrente                                *
c * dp(nnodev)       - delta p ( p(n  ,0  ) - p(0) )                   *
c * dpososity(nnodev)- nao definido                                    *
c * fnno             -  identifica dos nos de vertices                 *
c *                     ( 1 - vertice | 0 )                            *
c * nnode            - numero de nos                                   *
c * numel            - numero de elementos                             *
c * nen              - numero max. de nos por elemento                 *
c * nenv             - numero de nos de vertice por elemento           *
c * ndf              - numero max. de graus de liberdade por no        *
c * nst              - nst = nen*ndf                                   *
c * ndm              - dimensao                                        *
c * ndf              - numero max. de graus de liberdade por no        *
c * ntn              - numero max. de derivadas por no                 *
c * nprcs            - numero de processos                             *
c * isw              - codigo de instrucao para as rotinas             *
c *                    de elemento                                     *
c * ilib             - determina a biblioteca do elemento              *
c * i_xf             - buffer de comunicao para o numero de elementos  *
c *                    conectidos aos nos                              *
c * novlp            - chave indicativa de non-overlapping             * 
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * dpososity(nnodev)- delta phi ( phi(n  ,0  ) - phi(0) )             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Calcula as tensoes apenas nos vertices                             *     
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'termprop.fi'
c ... mpi
      integer*8 i_xfi
      logical novlp
c ...
      integer nnode,numel,nen,nenv,ndf,nst,ndm,ntn,npi
c ......................................................................      
      integer ix(nen+1,*),ic(*),fnno(*),ie(*)
      integer nel,ma,iel,i,j,k,k1,no,kk
      integer ilib,isw
      real*8  pl(nenv),plasticl(3,npi)
      real*8  plastic(3,npi,*),pc(*)
c ...
      logical ldum
      integer idum
      real*8 ddum
c ......................................................................
c
c ...
      do 30 i = 1, nnode
        ic(i)        = 0
        pc(i)        = 0.d0
   30 continue
c ......................................................................
c
c ... Loop nos elementos:
      do 900 nel = 1, numel
        kk = 0
c ... 
        do 300 i = 1, nenv
          pl(i) = 0.d0
  300   continue
c ......................................................................
c
c ... plasticidade
        do j = 1, npi  
c ... paramentro de endurecimento 
          plasticl(3,j) = plastic(3,j,nel)
        enddo  
c ......................................................................     
c
c ...... Arranjos de elemento:
        ma         = ix(nen+1,nel)
        iel        = ie(ma)
c ......................................................................
c
c ...... Chama biblioteca de elementos:
        call elmlib_pm(ddum,idum    ,idum
     1                ,ddum,ddum    ,ddum
     2                ,ddum,pl      ,ddum,ddum,ddum
     3                ,ddum,plasticl,ddum,idum
     4                ,ndm ,nst     ,nel  ,iel,isw
     5                ,ma  ,idum    ,ilib,ldum)
c ......................................................................
c
c ...... media do vetor global
        do 800 i = 1, nenv
           no = ix(i,nel)
           if (no .le. 0) go to 800
           ic(no) = ic(no) + 1
           pc(no) = pc(no) + pl(i)
  800   continue
  900 continue
c .......................................................................
c
c ... Comunica vetor de contagem de elementos por no'
      if (novlp) call allgatheri(ic,i_xfi)
c .......................................................................
      do 1000 i = 1, nnode
        if(fnno(i) .eq. 1 ) then
          pc(i) = pc(i)/ic(i)
        endif
 1000 continue
c ......................................................................
c
c ......................................................................
      return
      end      
c **********************************************************************
c
c **********************************************************************
      subroutine deltat_critico(ix   ,iq ,ie  ,e
     1                         ,x    ,xl
     2                         ,numel,nen,nenv,ndf
     3                         ,ndm  ,nst,isw ,ilib
     4                         ,my_id,mpi)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 31/08/2018                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_SQRM: calcula o numero de operacoes de pontos flutuantes      *   
c * do SQRM                                                            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ix(nen+1,numel) - conetividades nodais                             *
c * iq(7,numel)     - cargas nos elementos                             *
c * ie(numat)       - tipo de elemento                                 *
c * e(10,numat)     - constantes fisicas dos materiais                 *
c * x(ndm,nnode)    - coordenadas nodais                               *
c * xl(ndm,nen)     - nao definido                                     *
c * numel - numero de elementos                                        *
c * nen   - numero de nos por elemento                                 *
c * nenv  - numero de nos de vertice por elemento                      *
c * ndf   - numero max. de graus de liberdade por no                   *
c * ndm   - dimensao                                                   *
c * nst   - nen*ndf                                                    *
c * isw   - codigo de instrucao para as rotinas de elemento            *
c * ilib  - determina a biblioteca do elemento                         *   
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
      include 'mpif.h' 
      include 'omp_lib.h'
      include 'transiente.fi'
      include 'termprop.fi'
c ... mpi
      integer ierr,my_id
      logical mpi
c ...
      integer numel,nen,nenv,ndf,ndm,nst,nad,nadpu,stge,isw,numat,nlit
      integer neq,nequ
      integer ix(nen+1,*),iq(7,*),ie(*)
      integer iel,ma,nel,no,i,j,k,kk,nadr,ilib
      real*8  e(prop,*),x(ndm,*)
      real*8  xl(ndm,nenv),el(prop),dtc,dtc_min,g_dtc_min
c ... 
      logical ldum
      integer idum
      real*8 ddum
c ----------------------------------------------------------------------
c
c ... Loop nos elementos:
c     ------------------
      dtc_min = 3.1536d+10
      do 900 nel = 1, numel
        kk = 0
c...... loop nos nos do pressao
        do 400 i = 1, nenv
          no       = ix(i,nel)
          kk       = kk + 1
          do 300 j = 1, ndm
            xl(j,i) = x(j,no)
  300     continue
  400   continue 
c ......................................................................
c
c ...... Arranjos de elemento:
        ma  = ix(nen+1,nel)
        iel = ie(ma)
        do 610 i = 1, prop
          el(i) = e(i,ma)
  610   continue
c ......................................................................
c
c ...... Chama biblioteca de elementos:
        call elmlib_pm(el  ,iq(1,nel),iq(1,nel)
     1                ,xl  ,ddum     ,ddum
     2                ,ddum,dtc      ,ddum,ddum,ddum
     3                ,ddum,ddum     ,ddum,idum
     4                ,ndm,nst       ,nel ,iel,isw
     5                ,ma ,idum      ,ilib,ldum)
c ......................................................................
c 
c ...
        dtc_min = min(dtc, dtc_min)      
c ......................................................................
  900 continue
c ......................................................................
c
c ...    
      if(mpi) then
        call MPI_Allreduce(dtc_min,g_dtc_min,1,MPI_REAL8
     .                    ,MPI_MIN,MPI_COMM_WORLD,ierr)
      else
        g_dtc_min = dtc_min
      endif
c ......................................................................
c
c ...
      if( dtc_min .lt. dt ) then
        if(my_id .eq. 0 ) then
          print*,'Delta critico:',g_dtc_min
          print*,'Delta t      :',dt 
        endif
        call stop_mef()
      endif 
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine pform_mec(ix      ,iq         ,ie       ,e
     1                    ,x       ,id         ,ia       ,ja
     2                    ,au      ,al         ,ad       ,b
     3                    ,u       ,tx0
     4                    ,xl      ,ul         ,pl   
     5                    ,sl      ,ld         ,txnl
     6                    ,numel   ,nen        ,nenv     ,ndf 
     7                    ,ndm     ,nst        ,neq      ,nad  ,nadr 
     8                    ,lhs     ,rhs        ,unsym 
     9                    ,stge    ,isw        ,ilib     ,nlit
     1                    ,i_colorg,i_elcolor  ,numcolors,fstress0)
c **********************************************************************
c * Data de criacao    : 09/04/2016                                    *
c * Data de modificaco : 10/11/2016                                    * 
c * ------------------------------------------------------------------ * 
c * PFORM_MEC:                                                         *
c * ------------------------------------------------------------------ *      
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ix(nen+1,numel) - conetividades nodais                             *
c * iq(7,numel)     - cargas nos elementos                             *
c * ie(numat)       - tipo de elemento                                 *
c * e(10,numat)     - constantes fisicas dos materiais                 *
c * x(ndm,nnode)    - coordenadas nodais                               *
c * id(ndf,nnode)   - numeracao global das equacoes                    *
c * ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro         *
c *             coeficiente nao-nulo da equacao i, no formato CSR      *
c * ja(nad)   - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor a, no formato CSR (stge = 1)      *
c *           - ponteiro da diagonal do skyline (stge = 4)             *
c * ad(neq)   - nao definido                                           *
c * au(nad)   - nao definido                                           *
c * al(nad)   - nao definido                                           *
c * b(neq)    - vetor de forcas nodais equivalentes                    *
c * u(ndf,nnode) - solucao (com valores prescritos)                    *
c * dp(nnodev)   - delta p ( p(n  ,0  ) - p(0) )                       *
c * xl(ndm,nen)  - nao definido                                        *
c * ul(nst)      - nao definido                                        *
c * pl(nst)      - nao definido                                        *
c * sl(nst,nst)  - nao definido                                        *
c * ld(nst)      - nao definido                                        *
c * numel - numero de elementos                                        *
c * nen   - numero de nos por elemento                                 *
c * nenv  - numero de nos de vertice por elemento                      *
c * ndf   - numero max. de graus de liberdade por no                   *
c * ndm   - dimensao                                                   *
c * nst   - nen*ndf                                                    *
c * neq   - numero de equacoes                                         *
c * nequ  - numero de equacoes em kuu                                  *
c * nad   - numero de posicoes no CSR (storage = 1)                    *
c * nadr  - numero de termos nao nulos na parte retangular             *
c *         ( MPI em overllaping )                                     *
c * lhs   - flag para montagem de ad (letf hand side)                  *
c * rhs   - flag para correcao de b  (right hand side)                 *
c * unsym - flag para matrizes nao simetricas                          *
c * stge  - = 1, armazenamento CSR                                     *
c *         = 2, armazenamento por arestas                             *
c *         = 3, armazenamento EBE                                     *
c *         = 4, armazenamento SKYLINE                                 *
c *         = 0, nao monta matriz                                      *
c * isw   - codigo de instrucao para as rotinas de elemento            *
c * ilib  - determina a biblioteca do elemento                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c *     b(neq) - vetor de forcas corrigido     (rhs = true)            *
c *    ad(neq) - coeficientes da diagonal      (lhs = true)            *
c *    au(nad) - parte triangular sup. de A (lhs = true)               *
c *    al(nad) - parte triangular inf. de A (lhs = true, unsym = true) *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'openmp.fi'
      include 'omp_lib.h'
      include 'transiente.fi'
      include 'termprop.fi'
      integer numel,nen,nenv,ndf,ndm,nst,nad,nadr
      integer stge,isw,numat,nlit
      integer neq
      integer ix(nen+1,*),iq(7,*),ie(*),id(ndf,*),ld(nst)
      integer ia(*),ja(*)
      integer iel,ma,nel,no,i,j,k,kk,ilib
      integer i_colorg(2,*),i_elcolor(*),numcolors
      integer ic,jc
      real*8  e(prop,*),x(ndm,*),ad(*),au(*),al(*),b(*)
      real*8  u(ndf,*),tx0(6,*)
      real*8  xl(ndm,nenv),ul(nst),el(prop)
      real*8  pl(nst),sl(nst,nst),txnl(6,nen)
      logical lhs,rhs,unsym,fstress0
c .....................................................................
c
c ... Zera a matriz de coeficientes:
c
      if(lhs) then
        call azero(au,nad)
        call azero(al,nad+nadr)
        call azero(ad,neq)      
      endif
c ..................................................................... 
c
c ... openmp
      if(omp_elmt)then
c ... loop nos cores
c$omp parallel num_threads(nth_elmt)
c$omp.default(none)
c$omp.private(nel,ic,jc,no,k,ma,iel,i,j,kk,xl,ld,ul,pl,el,sl,txnl)
c$omp.shared(numcolors,i_colorg,i_elcolor,nen,nenv,ndf,ndm)
c$omp.shared(id,u,ix,tx0,ie,e)
c$omp.shared(stge,unsym,rhs,lhs,ilib,nlit,isw,nst,dt)
c$omp.shared(neq,nad,fstress0)
        do 150 ic = 1, numcolors
c$omp do
c ... Loop nos elementos:
c     ------------------
          do  100 jc = i_colorg(1,ic), i_colorg(2,ic)
            nel = i_elcolor(jc)
c           print*,omp_get_thread_num(),jc,nel,ic
            kk = 0
c ... loop nos nos de deslocamentos
            do i = 1, nen
              no = ix(i,nel)
              do j = 1, ndf
                kk     = kk + 1
                ld(kk) = id(j,no)
                ul(kk) = u(j,no)
                pl(kk) = 0.d0
              enddo
c ......................................................................              
           enddo   
c ......................................................................               
c
c ... loop nos nos vertices
            do i = 1, nenv
              no       = ix(i,nel)
              do j = 1, ndm
                xl(j,i) = x(j,no)
              enddo   
            enddo    
c ......................................................................
c            
c ... tensao inicial
          if(fstress0) then
            do i = 1, nen  
              no        = ix(i,nel)  
              txnl(1,i) = tx0(1,no)
              txnl(2,i) = tx0(2,no) 
              txnl(3,i) = tx0(3,no) 
              txnl(4,i) = tx0(4,no) 
              txnl(5,i) = tx0(5,no) 
              txnl(6,i) = tx0(6,no)
            enddo  
          endif  
c ......................................................................                 
c
c ...... Arranjos de elemento:
           ma  = ix(nen+1,nel)
           iel = ie(ma)
           do i = 1, prop
             el(i) = e(i,ma)
           enddo   
c ......................................................................
c
c ...... Chama biblioteca de elementos:
           call elmlib_mec(el,iq(1,nel),xl ,ul,pl ,sl ,txnl
     1                    ,dt,ndm     ,nst ,nel     ,iel,isw
     2                    ,ma,nlit    ,ilib)
c ......................................................................
c
c ...... Monta arrranjos locais em arranjos globais:
           call assbly(sl ,pl   ,ld  ,ia ,ja ,au
     1                ,al ,ad   ,b   ,nst,neq ,lhs
     2                ,rhs,unsym,stge)
c ......................................................................
  100     continue   
c$omp end do
c .....................................................................
  150   continue
c$omp end parallel
c .....................................................................
c
c ... sequencial
      else  
c ... Loop nos elementos:
c     ------------------
        do 250 nel = 1, numel
          kk = 0
c ... loop em todos os nos
          do i = 1, nen
            no = ix(i,nel)
            do j = 1, ndf
              kk     = kk + 1
              ld(kk) = id(j,no)
              pl(kk) = 0.d0
              ul(kk) = u(j,no)
            enddo
c ......................................................................
         enddo
c ......................................................................      
c
c ... loop nos nos vertices
          do i = 1, nenv
            no       = ix(i,nel)
            do j = 1, ndm
              xl(j,i) = x(j,no)
            enddo   
          enddo
c .....................................................................
c            
c ... tensao inicial
          if(fstress0) then
            do i = 1, nen  
              no        = ix(i,nel)  
              txnl(1,i) = tx0(1,no)
              txnl(2,i) = tx0(2,no) 
              txnl(3,i) = tx0(3,no) 
              txnl(4,i) = tx0(4,no) 
              txnl(5,i) = tx0(5,no) 
              txnl(6,i) = tx0(6,no)
            enddo  
          endif  
c ......................................................................      
c
c ...... Arranjos de elemento:
          ma  = ix(nen+1,nel)
          iel = ie(ma)
          do i = 1, prop
            el(i) = e(i,ma)
          enddo   
c ......................................................................
c
c ...... Chama biblioteca de elementos:
          call elmlib_mec(el,iq(1,nel),xl ,ul  ,pl ,sl ,txnl
     1                   ,dt,ndm     ,nst ,nel ,iel,isw
     2                   ,ma,nlit    ,ilib)
c ......................................................................
c
c ...... Monta arrranjos locais em arranjos globais:
          call assbly(sl ,pl   ,ld  ,ia ,ja ,au
     1               ,al ,ad   ,b   ,nst,neq ,lhs
     2               ,rhs,unsym,stge)
c ......................................................................
  250   continue 
c ......................................................................
      endif
c ......................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine tform_mec(ix    ,x    ,e   ,ie
     1                    ,ic    ,xl   ,ul  
     2                    ,pl    ,u    ,tx0 ,t   
     3                    ,nnodev,numel,nenv,nen
     4                    ,ndm   ,ndf  ,nst ,ntn  
     5                    ,isw   ,ilib)
c **********************************************************************
c * Data de criacao    : 09/04/2016                                    *
c * Data de modificaco : 31/10/2016                                    * 
c * ------------------------------------------------------------------ * 
c * TFORM: caluclo das tensoes e dos fluxo nodal nos vertices          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ix(nen+1,numel)  - conetividades nodais                            *
c * x(ndm,nnode)     - coordenadas nodais                              *
c * e(10,numat)      - constantes fisicas dos materiais                *
c * ie(numat)        - tipo de elemento                                *
c * ic(nnode)        - nao definido                                    *
c * xl(ndm,nen)      - nao definido                                    *
c * ul(ndf,nen)      - nao definido                                    *
c * pl(ntn*nen)      - nao definido                                    *
c * u(ndf,nnode)     - solucao corrente                                *
c * tx0(ntn,nnodev)  - tensao inicial                                  *
c * t(ntn,nnodev)    - nao definido                                    *
c * nnodev           - numero de nos de vertices                       *
c * numel            - numero de elementos                             *
c * nen              - numero max. de nos por elemento                 *
c * nenv             - numero de nos de vertice por elemento           *
c * ndf              - numero max. de graus de liberdade por no        *
c * nst              - nst = nen*ndf                                   *
c * ndm              - dimensao                                        *
c * ndf              - numero max. de graus de liberdade por no        *
c * ntn              - numero max. de derivadas por no                 *
c * ovlp             - chave indicativa de overlapping                 *
c * nprcs            - numero de processos                             *
c * isw              - codigo de instrucao para as rotinas             *
c *                    de elemento                                     *
c * ilib             - determina a biblioteca do elemento              *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * t(ntn   ,nnodev) - tensoes medias nodais                           *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Calcula as tensoes apenas nos vertices                             *     
c **********************************************************************
      use Malloc
      implicit none
c     include 'mpif.h'
      include 'parallel.fi'
      include 'termprop.fi'
      integer nnodev,numel,nen,nenv,ndf,nst,ndm,ntn
c ......................................................................      
      integer ix(nen+1,*),ie(*),ic(nnodev)
      integer nel,ma,iel,i,j,k,k1,no,kk
      integer ilib,isw
      real*8  xl(ndm,nenv),ul(nst),pl(nenv*ntn)
      real*8  x(ndm,*),e(prop,*)
      real*8  u(ndf,*),el(prop),tx0(ntn,*)
      real*8  t(ntn,*)
c ...
      logical ldum
      integer idum
      real*8 ddum
c ......................................................................
c
c ... deslocamentos do vetor local p
c ......................................................................
c
c ...
      do 30 i = 1, nnodev
        ic(i) = 0
c ... tensao
        do 10 j = 1, ntn
           t(j,i) = 0.d0
   10   continue 
c ......................................................................
c
   30 continue 
c ......................................................................
c
c ... Loop nos elementos:
      do 900 nel = 1, numel
        kk = 0
c ...
        do 300 i = 1, nenv*ntn
          pl(i) = 0.d0
  300   continue
c .......................................................................
c
c ... loop em todos os nos
        do 400 i = 1, nen
          no = ix(i,nel)
c ... loop nos deslocamentos
          do 410 j = 1, ndf
            kk     = kk + 1
            ul(kk) = u(j,no)
  410     continue
c ......................................................................
  400   continue
c ......................................................................
c
c ... loop nos nos vertices
        do 510 i = 1, nenv
          no     = ix(i,nel)
          kk     = kk + 1
          do 500 j = 1, ndm
            xl(j,i) = x(j,no)
  500     continue
  510   continue
c ......................................................................
c
c ...... form element array
        ma  = ix(nen+1,nel)
        iel = ie(ma)      
        do 610 i = 1, prop
          el(i) = e(i,ma)
  610   continue
c ......................................................................
c
c ...... Chama biblioteca de elementos:
        call elmlib_mec(el,idum,xl,ul,pl,ddum,ddum,ddum,ndm,nst,nel
     .                 ,iel,isw,ma,idum,ilib)
c ...... media do vetor global
        do 800 i = 1, nenv
           no = ix(i,nel)
           if (no .le. 0) go to 800
           ic(no) = ic(no) + 1  
           do 700 j = 1, ntn
c ... tensao total
              k       = (i-1)*ntn + j
              t(j,no) = t(j,no)  + pl(k)
  700      continue
  800   continue
  900 continue
c .......................................................................
c
c ... Comunica vetor de contagem de elementos por no'
c     if (novlp) call allgatheri(ic,i_xfi)
c .......................................................................
      do 1000 i = 1, nnodev
c ... tensao
        do 1010 j = 1, ntn
c ... tensao total            
         t(j,i)  = t(j,i)/ic(i)  + tx0(j,i)
 1010   continue 
c ......................................................................
 1000 continue
c ......................................................................
c
c ......................................................................
      return
      end      
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 24/11/2016                                    *
c * Data de modificaco : 00/00/0006                                    * 
c * ------------------------------------------------------------------ * 
c * TFORM: caluclo das tensoes e dos fluxo nodal nos vertices          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * txe              - nao definido                                    *
c * t                - tensao total                                    *
c * u                - pressao u(4,i)                                  *
c * nnode            - numero de nos                                   *
c * ntn              - numero de tensoes                               *
c * ndf              - graus de liberade                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * te(ntn   ,nnode ) - tensoes efetiva                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine effective_stress(txe  ,tx ,u
     .                           ,nnode,ntn,ndf)
      implicit none
      integer ntn,ndf,i,j,nnode
      real*8 txe(ntn,*),tx(ntn,*),u(ndf,*)
c ...
      if(ntn .eq. 6 ) then
        do i = 1, nnode
          txe(1,i) =  tx(1,i) + u(ndf,i) 
          txe(2,i) =  tx(2,i) + u(ndf,i) 
          txe(3,i) =  tx(3,i) + u(ndf,i) 
          txe(4,i) =  tx(4,i)  
          txe(5,i) =  tx(5,i)  
          txe(6,i) =  tx(6,i)  
        enddo
      endif  
c ......................................................................
      return
      end  
c **********************************************************************
c     subroutine initial_pc(ix,ie,e,txp,u,numel,nen,npi,ndf,ntn)
c     include 'termprop.fi'
c     implicit none
c     integer ix(nen+1,*),ie(*),numel,nen,npi,ndf
c     real*8 e(prop,*),txp(ntn,npi,*),u(ndf,*)
c     real*8 p,q,mcs
c 
c
c     do i = 1, numel        
c         ma         = ix(nen+1,nel)
c         iel        = ie(ma)
c         pc         = e(12)
c         do j = 1, npi
c           sm = txp(1,j,i) + txp(2,j,i) + txp(3,j,i) + 3*u(ndf, 
c           p  = (1.0d0/3.d0)*(txp(1,j,i)
c         enddo
c **********************************************************************
c **********************************************************************
c     subroutine pform(ix,iq,ie,e,x,id,ia,ja,au,al,ad,b,u,u0,ut,v,a,du,
c    .         hi,tx,txp,eps,w,xl,ul,vl,zl,wl,pl,sl,tl,ld,numel,nen,ndf,
c    .         ndm,npi,nst,neq,nad,lhs,rhs,unsym,stge,isw,ilib,
c    .         i_colorg,i_elcolor,numcolors,nlit,flaghidr)
c **********************************************************************
c *                                                                    *
c *   PFORM:                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    iq(7,numel)     - cargas nos elementos                          *
c *    ie(numat)       - tipo de elemento                              *
c *    e(10,numat)     - constantes fisicas dos materiais              *
c *    x(ndm,nnode)    - coordenadas nodais                            *
c *    id(ndf,nnode)   - numeracao global das equacoes                 *
c *    ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro      *
c *                coeficiente nao-nulo da equacao i, no formato CSR   *
c *    ja(nad)   - ja(k) informa a coluna do coeficiente que ocupa     *
c *                a posicao k no vetor a, no formato CSR (stge = 1)   *
c *              - ponteiro da diagonal do skyline (stge = 4)          *
c *    ad(neq)   - nao definido                                        *
c *    au(nad)   - nao definido                                        *
c *    al(nad)   - nao definido                                        *
c *    b(neq)    - vetor de forcas nodais equivalentes                 *
c *    u(ndf,nnode) - solucao (com valores prescritos)                 *
c *    v(ndf,nnode) - derivada no tempo da solucao                     *
c *    a(ndf,nnode) - derivada segunda no tempo da solucao             *
c *    w(ndm,nnode) - campo de velocidades de conveccao                *
c *    xl(ndm,nen)  - nao definido                                     *
c *    ul(nst)      - nao definido                                     *
c *    vl(nst)      - nao definido                                     *
c *    zl(nst)      - nao definido                                     *
c *    wl(ndm,nen)  - nao definido                                     *
c *    pl(nst)      - nao definido                                     *
c *    sl(nst,nst)  - nao definido                                     *
c *    ld(nst)      - nao definido                                     *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    ndm   - dimensao                                                *
c *    nst   - nen*ndf                                                 *
c *    neq   - numero de equacoes                                      *
c *    nad   - numero de posicoes no CSR          (storage = 1)        *
c *    adfl  - flag para montagem de ad                                *
c *    bdfl  - flag para armazenamento da diagonal em blocos           *
c *    bfl   - flag para correcao de b                                 *
c *    unsym - flag para matrizes nao simetricas                       *
c *    stge  - = 1, armazenamento CSR                                  *
c *            = 2, armazenamento por arestas                          *
c *            = 3, armazenamento EBE                                  *
c *            = 4, armazenamento SKYLINE                              *
c *            = 0, nao monta matriz                                   *
c *    isw   - codigo de instrucao para as rotinas de elemento         *
c *    ilib  - determina a biblioteca do elemento                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     b(neq) - vetor de forcas corrigido     (rhs = true)            *
c *    ad(neq) - coeficientes da diagonal      (lhs = true)            *
c *    au(nad) - parte triangular sup. de A (lhs = true)               *
c *    al(nad) - parte triangular inf. de A (lhs = true, unsym = true) *
c *                                                                    *
c **********************************************************************
c     implicit none
c      include 'openmp.fi'
c     include 'omp_lib.h'
c     include 'transiente.fi'
c     include 'termprop.fi'
c     integer numel,nen,ndf,ndm,nst,neq,nad,stge,isw,numat,nlit,npi
c     integer ix(nen+1,*),iq(7,*),ie(*),id(ndf,*),ld(ndf,nen)
c     integer ia(*),ja(*)
c     integer iel,ma,nel,no,i,j,k,nadr,ilib
c     integer i_colorg(2,*),i_elcolor(*),numcolors
c     integer ic,jc
c     real*8  e(prop,*),x(ndm,*),ad(neq),au(*),al(*),b(neq),du(*)
c     real*8  u(ndf,*),v(ndf,*),a(ndm,*),w(ndm,*),ut(1,*),hi(3,*)
c     real*8  xl(ndm,nen),ul(ndf,nst),vl(ndf,nen),zl(ndf,nen),el(prop)
c     real*8  wl(ndm,nen),pl(ndf,nen),sl(nst,nst),tl(nen),utl(nen)
c     real*8  tx(npi,*),eps(npi,*),u0(ndf*nen,*),txp(npi,*)
c     logical lhs,rhs,unsym,flaghidr
c ----------------------------------------------------------------------
c ... Verifica se a matriz eh retangular e calc. o no. de coeficientes:
c     nadr = 0
c     if(stge .eq. 1) nadr = ja(nad+1)
c     if(nadr .gt. 0) nadr = ia(2*neq+2)-1
c ......................................................................      
c
c.... Zera a matriz de coeficientes:
c
c     if(lhs) then
c        call azero(au,nad)
c        call azero(al,nad+nadr)
c        call azero(ad,neq)      
c     endif
c ----------------------------------------------------------------------
c
c.... Loop nos elementos:
c     ------------------
cc$omp parallel 
cc$omp.private(nel,jc,no,k,ma,iel,i,j,xl,wl,ld,ul,vl,pl,sl,zl,tl,utl,el)
c     do 1000 ic = 1, numcolors
cc$omp do        
c       do 900 jc = i_colorg(1,ic), i_colorg(2,ic)
c          nel = i_elcolor(jc)      
c
c...... Arranjos locais:
c
c       do 600 i = 1, nen
c         no = ix(i,nel)
c         if (no .gt. 0) go to 300
c            do 100 j = 1, ndm
c              xl(j,i) = 0.d0
c              wl(j,i) = 0.d0
c 100        continue
c            do 200 j = 1, ndf
c              ld(j,i) = 0
c              ul(j,i) = 0.d0
c              vl(j,i) = 0.d0
c              zl(j,i) = 0.d0               
c              pl(j,i) = 0.d0
c              tl(i)   = 0.d0
c 200        continue
c            go to 600
c 300     continue
c         do 400 j = 1, ndm
c            xl(j,i) = x(j,no)
c            wl(j,i) = w(j,no)             
c 400     continue
c         do 500 j = 1, ndf
c            k = id(j,no)
c            ld(j,i) = k
c            pl(j,i) = 0.d0
c            ul(j,i) = u(j,no)
c            vl(j,i) = v(j,no)
c            zl(j,i) = a(j,no)
c 500     continue
c         tl(i) = du(no)
c         utl(i)= ut(1,no)
c 600   continue
c ......................................................................
c
c ...... Arranjos de elemento:
c
c       ma  = ix(nen+1,nel)
c       iel = ie(ma)
c       do 610 i = 1, prop
c          el(i) = e(i,ma)
c 610   continue
c ...... Chama biblioteca de elementos:
c       call elmlib(el,iq(1,nel),xl,ul,vl,zl,wl,pl,sl,tl,utl,hi(1,nel),
c    .              u0(1,nel),tx(1,nel),txp(1,nel),eps(1,nel),ndm,nst,
c    .              nel,iel,isw,ma,nlit,ilib,flaghidr,1.d0)
c ----------------------------------------------------------------------
c ...... Monta arrranjos locais em arranjos globais:
c       call assbly(sl,pl,ld,ia,ja,au,al,ad,b,nst,neq,lhs,rhs,unsym,
c    .              stge)
c 900 continue
cc$omp end do
c1000 continue
cc$omp end parallel
c ......................................................................
c     return
c     end
c     subroutine tform(ix,x,e,ie,ic,xl,ul,pl,tl,u,du,hi,tx,eps,t,nnode,
c    .             numel,nen,ndm,npi,ndf,nst,ntn,i_xfi,ilib,flaghidr)
c **********************************************************************
c *                                                                    *
c *   Subroutine TFORM                                                 *
c *   ----------------                                                 *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    x(ndm,nnode)    - coordenadas nodais                            *
c *    e(10,numat)     - constantes fisicas dos materiais              *
c *    ie(numat)       - tipo de elemento                              *
c *    ic(nnode)   - nao definido                                      *
c *    xl(ndm,nen) - nao definido                                      *
c *    ul(ndf,nen) - nao definido                                      *
c *    pl(ntn*nen) - nao definido                                      *
c *    u(ndf,nnode)- solucao corrente                                  *
c *    stres(nte,numel) - tensoes por elemento                         *
c *    t(ntn,nnode)- nao definido                                      *
c *    nnode - numero de nos                                           *
c *    numel - numero de elementos                                     *
c *    nen   - numero max. de nos por elemento                         *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - nst = nen*ndf                                           *
c *    ndm   - dimensao                                                *
c *    ntn   - numero max. de derivadas por no                         *
c *    ovlp  - chave indicativa de overlapping                         *
c *    nprcs - numero de processos                                     *
c *    ilib - determina a biblioteca do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *    t(ntn,nnode)- tensoes medias nodais                             *
c *                                                                    *
c **********************************************************************
c     use Malloc
c     implicit none
c     include 'mpif.h'
c     include 'parallel.fi'
c     include 'termprop.fi'
c     integer nnode,numel,nen,ndf,nst,ndm,ntn,npi
c ... ponteiro      
c     integer*8 i_xfi
c     logical flaghidr
c ......................................................................      
c     integer ix(nen+1,*),ie(*),ic(*),dum
c     integer nel,ma,iel,i,j,k,no
c     integer ilib
c     real*8  x(ndm,*),e(prop,*),xl(ndm,*),ul(ndf,*),pl(*),tl(*)
c     real*8  utl(nen),u(ndf,*),t(ntn,*),du(*),hi(3,*),tx(npi,*)
c     real*8  eps(npi,*)
c ......................................................................
c     do 20 i = 1, nnode
c        ic(i) = 0
c        do 10 j = 1, ntn
c           t(j,i) = 0.d0
c  10    continue
c  20 continue
c ......................................................................
c .... loop on elements:
c     do 900 nel = 1, numel
c ...... set up local arrays
c       do 600 i = 1, nen
c         no = ix(i,nel)
c         if (no .gt. 0) go to 300
c         do 100 j = 1, ndm
c           xl(j,i) = 0.d0
c 100     continue
c         do 200 j = 1, ndf
c           ul(j,i)  = 0.d0
c 200     continue
c         go to 600
c 300     continue
c         do 400 j = 1, ndm
c           xl(j,i) = x(j,no)
c 400     continue
c         do 500 j = 1, ndf
c           ul(j,i) = u(j,no)
c 500     continue
c         tl(i) = du(no)
c 600   continue
c ...... form element array
c       ma  = ix(nen+1,nel)
c        ma  = ix(21,nel)
c       iel = ie(ma)      
c ...... Biblioteca de elementos:
c       if (iel .eq. 10 .or. iel .eq. 11 .or. iel .eq. 12)  goto 900
c       call elmlib(e(1,ma),dum,xl,ul,pl,pl,pl,pl,pl,tl,utl,hi(1,nel),
c    .    pl,tx(1,nel),tx(1,nel),eps(1,nel),ndm,nst,nel,iel,3,ma,1,ilib,
c    .    flaghidr)
c ...... add to global array
c 650   do 800 i = 1, nen
c          no = ix(i,nel)
c          if (no .le. 0) go to 800
c          ic(no) = ic(no) + 1
c          do 700 j = 1, ntn
c             k = (i-1)*ntn + j
c             t(j,no) = t(j,no) + pl(k)
c 700      continue
c 800   continue
c 900 continue
c .......................................................................
c ... Comunica vetor de contagem de elementos por no'
c     if (novlp) call allgatheri(ic,i_xfi)
c .......................................................................
c     do 1000 i = 1, nnode
c     do 1000 j = 1, ntn
c        t(j,i) = t(j,i)/ic(i)
c1000 continue
c ......................................................................
c     return
c     end      
c     subroutine t_efetive(ix,x,e,ie,ic,xl,ul,pl,tl,u,du,hi,tx,eps,t
c    .         ,nnode,numel,nen,ndm,npi,ndf,nst,ntn,i_xfi,ilib,flaghidr,
c    .                                         plastic)
c **********************************************************************
c *                                                                    *
c *   Subroutine T_EFETIVE                                             *
c *   ----------------                                                 *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    x(ndm,nnode)    - coordenadas nodais                            *
c *    e(10,numat)     - constantes fisicas dos materiais              *
c *    ie(numat)       - tipo de elemento                              *
c *    ic(nnode)   - nao definido                                      *
c *    xl(ndm,nen) - nao definido                                      *
c *    ul(ndf,nen) - nao definido                                      *
c *    pl(ntn*nen) - nao definido                                      *
c *    u(ndf,nnode)- solucao corrente                                  *
c *    stres(nte,numel) - tensoes por elemento                         *
c *    t(ntn,nnode)- nao definido                                      *
c *    nnode - numero de nos                                           *
c *    numel - numero de elementos                                     *
c *    nen   - numero max. de nos por elemento                         *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - nst = nen*ndf                                           *
c *    ndm   - dimensao                                                *
c *    ntn   - numero max. de derivadas por no                         *
c *    ovlp  - chave indicativa de overlapping                         *
c *    nprcs - numero de processos                                     *
c *    ilib - determina a biblioteca do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *    t(ntn,nnode)- tensoes medias nodais                             *
c *                                                                    *
c **********************************************************************
c     use Malloc
c     implicit none
c     include 'mpif.h'
c     include 'parallel.fi'
c     include 'termprop.fi'
c     integer nnode,numel,nen,ndf,nst,ndm,ntn,npi
c ... ponteiro      
c     integer*8 i_xfi
c     logical flaghidr
c ......................................................................      
c     integer ix(nen+1,*),ie(*),ic(*),dum
c     integer nel,ma,iel,i,j,k,no
c     integer ilib
c     real*8 plastic(2,*)
c     real*8  x(ndm,*),e(prop,*),xl(ndm,*),ul(ndf,*),pl(*),tl(*)
c     real*8  utl(nen),u(ndf,*),du(*),hi(3,*),tx(npi,*)
c     real*8 t(ntn,*)
c     real*8  eps(npi,*)
c ......................................................................
c      do 20 i = 1, nnode
c         ic(i) = 0
c         do 10 j = 1, ntn
c            t(j,i) = 0.d0
c   10    continue
c   20 continue
c ......................................................................
c .... loop on elements:
c     do 900 nel = 1, numel
c ...... set up local arrays
c       do 600 i = 1, nen
c         no = ix(i,nel)
c         if (no .gt. 0) go to 300
c         do 100 j = 1, ndm
c           xl(j,i) = 0.d0
c 100     continue
c         do 200 j = 1, ndf
c           ul(j,i)  = 0.d0
c 200     continue
c         go to 600
c 300     continue
c         do 400 j = 1, ndm
c           xl(j,i) = x(j,no)
c 400     continue
c         do 500 j = 1, ndf
c           ul(j,i) = u(j,no)
c 500     continue
c         tl(i) = du(no)
c 600   continue
c ...... form element array
c       ma  = ix(nen+1,nel)
c        ma  = ix(21,nel)
c       iel = ie(ma)      
c ...... Biblioteca de elementos:
c        if (iel .eq. 11 .or. iel .eq. 12)  then
c         call elmlib(e(1,ma),dum,xl,ul,pl,pl,pl,pl,pl,tl,utl,hi(1,nel),
c    .    pl,tx(1,nel),tx(1,nel),eps(1,nel),ndm,nst,nel,iel,5,ma,1,ilib,
c    .    flaghidr,plastic(1,nel))
c        else
c           plastic(1,nel) = 0
c           plastic(2,nel) = 0
c        endif
c 900 continue
c .......................................................................
c     return
c     end 
c ***********************************************************************
cr
