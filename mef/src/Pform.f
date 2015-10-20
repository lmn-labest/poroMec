c*****************************Svn***************************************      
c*$Date: 2015-06-09 16:23:50 -0300 (Tue, 09 Jun 2015) $                 
c*$Rev: 971 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************    
      subroutine pformpmec(ix,iq,ie,e,x,id,ia,ja,au,al,ad,b,u,
     .                     xl,ul,pl,sl,ld,numel,nen,nenv,ndf,
     .                     ndm,nst,neq,nequ,nad,nadpu,lhs,rhs,unsym,
     .                     stge,isw,ilib,nlit,dualCsr)
c **********************************************************************
c *                                                                    *
c *   PFORMPMEC:                                                       *
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
c *    xl(ndm,nen)  - nao definido                                     *
c *    ul(nst)      - nao definido                                     *
c *    pl(nst)      - nao definido                                     *
c *    sl(nst,nst)  - nao definido                                     *
c *    ld(nst)      - nao definido                                     *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    nenv  - numero de nos de vertive por elemento                   *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    ndm   - dimensao                                                *
c *    nst   - nen*ndf                                                 *
c *    neq   - numero de equacoes                                      *
c *    nequ  - numero de equacoes em kuu                               *   
c *    nad   - numero de posicoes no CSR (storage = 1)                 *
c *    nadpu - numero de posicoes no CSR (storage = 1,dualCsr = true)  *      
c *    lhs   - flag para montagem de ad                                *
c *    rhs   - flag para correcao de b                                 *
c *    unsym - flag para matrizes nao simetricas                       *
c *    stge  - = 1, armazenamento CSR                                  *
c *            = 2, armazenamento por arestas                          *
c *            = 3, armazenamento EBE                                  *
c *            = 4, armazenamento SKYLINE                              *
c *            = 0, nao monta matriz                                   *
c *    isw   - codigo de instrucao para as rotinas de elemento         *
c *    ilib  - determina a biblioteca do elemento                      *
c *    dualCsr     - true - armazenamento em blocos Kuu,Kpp e kpu      *
c *                  false- aramzenamento em unico bloco               *      
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     b(neq) - vetor de forcas corrigido     (rhs = true)            *
c *    ad(neq) - coeficientes da diagonal      (lhs = true)            *
c *    au(nad) - parte triangular sup. de A (lhs = true)               *
c *    al(nad) - parte triangular inf. de A (lhs = true, unsym = true) *
c *                                                                    *
c **********************************************************************
      implicit none
c      include 'openmp.fi'
      include 'omp_lib.h'
      include 'transiente.fi'
      include 'termprop.fi'
      integer numel,nen,nenv,ndf,ndm,nst,nad,nadpu,stge,isw,numat,nlit
      integer neq,nequ
      integer ix(nen+1,*),iq(7,*),ie(*),id(ndf,*),ld(nst)
      integer ia(*),ja(*)
      integer iel,ma,nel,no,i,j,k,kk,nad1,ilib
      integer ic,jc
      real*8  e(prop,*),x(ndm,*),ad(*),au(*),al(*),b(*)
      real*8  u(ndf,*)
      real*8  xl(ndm,nenv),ul(nst),el(prop)
      real*8  pl(nst),sl(nst,nst)
      logical lhs,rhs,unsym,dualCsr
c ----------------------------------------------------------------------
c
c.... Zera a matriz de coeficientes:
c
      if(lhs) then
         call azero(au,nad)
         call azero(al,nad+nadpu)
         call azero(ad,neq)      
      endif
c ----------------------------------------------------------------------
c
c.... Loop nos elementos:
c     ------------------
      kk = 0
      do 900 nel = 1, numel
c ... loop nos nos
         do 600 i = 1, nen
            no = ix(i,nel)
            do 500 j = 1, ndf - 1
               kk     = kk + 1
               ld(kk) = id(j,no)
               pl(kk) = 0.d0
               ul(kk) = u(j,no)
  500       continue
  600    continue
c
c...... loop nos nos dos vertices
         do 400 i = 1, nenv
            no = ix(i,nel)
            kk = kk + 1
            ld(kk)    = id(ndf,no)
            pl(kk)    = 0.d0
            ul(kk)    = u(ndf,no)
            do 300 j = 1, ndm
               xl(j,i) = x(j,no)
  300       continue
  400    continue 
c ......................................................................
c
c ...... Arranjos de elemento:
c
         ma  = ix(nen+1,nel)
         iel = ie(ma)
         do 610 i = 1, prop
            el(i) = e(i,ma)
  610    continue
c ...... Chama biblioteca de elementos:
         call elmlibpmec(el,iq(1,nel),xl,ul,pl,sl,dt,ndm,nst,nel,iel,
     .                   isw,ma,nlit,ilib)
c ......................................................................
c
c ...... Monta arrranjos locais em arranjos globais:
         call assbly(sl,pl,ld,ia,ja,au,al,ad,b,nst,neq,nequ,nad,
     .              lhs,rhs,unsym,stge,dualCsr)
  900 continue
c ......................................................................
      if (lhs) then 
        call printadal(ad,al,al(nad+1),ja,ja(nad+1),neq,nad,nadpu)
        call printAx(sl,b,neq)
      endif
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine pform(ix,iq,ie,e,x,id,ia,ja,au,al,ad,b,u,u0,ut,v,a,du,
     .         hi,tx,txp,eps,w,xl,ul,vl,zl,wl,pl,sl,tl,ld,numel,nen,ndf,
     .         ndm,npi,nst,neq,nad,lhs,rhs,unsym,stge,isw,ilib,
     .         i_colorg,i_elcolor,numcolors,nlit,flaghidr)
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
      implicit none
c      include 'openmp.fi'
      include 'omp_lib.h'
      include 'transiente.fi'
      include 'termprop.fi'
      integer numel,nen,ndf,ndm,nst,neq,nad,stge,isw,numat,nlit,npi
      integer ix(nen+1,*),iq(7,*),ie(*),id(ndf,*),ld(ndf,nen)
      integer ia(*),ja(*)
      integer iel,ma,nel,no,i,j,k,nad1,ilib
      integer i_colorg(2,*),i_elcolor(*),numcolors
      integer ic,jc
      real*8  e(prop,*),x(ndm,*),ad(neq),au(*),al(*),b(neq),du(*)
      real*8  u(ndf,*),v(ndf,*),a(ndm,*),w(ndm,*),ut(1,*),hi(3,*)
      real*8  xl(ndm,nen),ul(ndf,nst),vl(ndf,nen),zl(ndf,nen),el(prop)
      real*8  wl(ndm,nen),pl(ndf,nen),sl(nst,nst),tl(nen),utl(nen)
      real*8  tx(npi,*),eps(npi,*),u0(ndf*nen,*),txp(npi,*)
      logical lhs,rhs,unsym,flaghidr
c ----------------------------------------------------------------------
c ... Verifica se a matriz eh retangular e calc. o no. de coeficientes:
      nad1 = 0
      if(stge .eq. 1) nad1 = ja(nad+1)
      if(nad1 .gt. 0) nad1 = ia(2*neq+2)-1
c ......................................................................      
c
c.... Zera a matriz de coeficientes:
c
      if(lhs) then
         call azero(au,nad)
         call azero(al,nad+nad1)
         call azero(ad,neq)      
      endif
c ----------------------------------------------------------------------
c
c.... Loop nos elementos:
c     ------------------
c$omp parallel 
c$omp.private(nel,jc,no,k,ma,iel,i,j,xl,wl,ld,ul,vl,pl,sl,zl,tl,utl,el)
      do 1000 ic = 1, numcolors
c$omp do        
        do 900 jc = i_colorg(1,ic), i_colorg(2,ic)
           nel = i_elcolor(jc)      
c
c...... Arranjos locais:
c
        do 600 i = 1, nen
          no = ix(i,nel)
          if (no .gt. 0) go to 300
             do 100 j = 1, ndm
               xl(j,i) = 0.d0
               wl(j,i) = 0.d0
  100        continue
             do 200 j = 1, ndf
               ld(j,i) = 0
               ul(j,i) = 0.d0
               vl(j,i) = 0.d0
               zl(j,i) = 0.d0               
               pl(j,i) = 0.d0
               tl(i)   = 0.d0
  200        continue
             go to 600
  300     continue
          do 400 j = 1, ndm
             xl(j,i) = x(j,no)
             wl(j,i) = w(j,no)             
  400     continue
          do 500 j = 1, ndf
             k = id(j,no)
             ld(j,i) = k
             pl(j,i) = 0.d0
             ul(j,i) = u(j,no)
             vl(j,i) = v(j,no)
             zl(j,i) = a(j,no)
  500     continue
          tl(i) = du(no)
          utl(i)= ut(1,no)
  600   continue
c ......................................................................
c
c ...... Arranjos de elemento:
c
        ma  = ix(nen+1,nel)
        iel = ie(ma)
        do 610 i = 1, prop
           el(i) = e(i,ma)
  610   continue
c ...... Chama biblioteca de elementos:
        call elmlib(el,iq(1,nel),xl,ul,vl,zl,wl,pl,sl,tl,utl,hi(1,nel),
     .              u0(1,nel),tx(1,nel),txp(1,nel),eps(1,nel),ndm,nst,
     .              nel,iel,isw,ma,nlit,ilib,flaghidr,1.d0)
c ----------------------------------------------------------------------
c ...... Monta arrranjos locais em arranjos globais:
        call assbly(sl,pl,ld,ia,ja,au,al,ad,b,nst,neq,lhs,rhs,unsym,
     .              stge)
  900 continue
c$omp end do
 1000 continue
c$omp end parallel
c ......................................................................
      return
      end
      subroutine tform(ix,x,e,ie,ic,xl,ul,pl,tl,u,du,hi,tx,eps,t,nnode,
     .             numel,nen,ndm,npi,ndf,nst,ntn,i_xfi,ilib,flaghidr)
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
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'termprop.fi'
      integer nnode,numel,nen,ndf,nst,ndm,ntn,npi
c ... ponteiro      
      integer*8 i_xfi
      logical flaghidr
c ......................................................................      
      integer ix(nen+1,*),ie(*),ic(*),dum
      integer nel,ma,iel,i,j,k,no
      integer ilib
      real*8  x(ndm,*),e(prop,*),xl(ndm,*),ul(ndf,*),pl(*),tl(*)
      real*8  utl(nen),u(ndf,*),t(ntn,*),du(*),hi(3,*),tx(npi,*)
      real*8  eps(npi,*)
c ......................................................................
      do 20 i = 1, nnode
         ic(i) = 0
         do 10 j = 1, ntn
            t(j,i) = 0.d0
   10    continue
   20 continue
c ......................................................................
c .... loop on elements:
      do 900 nel = 1, numel
c ...... set up local arrays
        do 600 i = 1, nen
          no = ix(i,nel)
          if (no .gt. 0) go to 300
          do 100 j = 1, ndm
            xl(j,i) = 0.d0
  100     continue
          do 200 j = 1, ndf
            ul(j,i)  = 0.d0
  200     continue
          go to 600
  300     continue
          do 400 j = 1, ndm
            xl(j,i) = x(j,no)
  400     continue
          do 500 j = 1, ndf
            ul(j,i) = u(j,no)
  500     continue
          tl(i) = du(no)
  600   continue
c ...... form element array
        ma  = ix(nen+1,nel)
c        ma  = ix(21,nel)
        iel = ie(ma)      
c ...... Biblioteca de elementos:
        if (iel .eq. 10 .or. iel .eq. 11 .or. iel .eq. 12)  goto 900
        call elmlib(e(1,ma),dum,xl,ul,pl,pl,pl,pl,pl,tl,utl,hi(1,nel),
     .    pl,tx(1,nel),tx(1,nel),eps(1,nel),ndm,nst,nel,iel,3,ma,1,ilib,
     .    flaghidr)
c ...... add to global array
  650   do 800 i = 1, nen
           no = ix(i,nel)
           if (no .le. 0) go to 800
           ic(no) = ic(no) + 1
           do 700 j = 1, ntn
              k = (i-1)*ntn + j
              t(j,no) = t(j,no) + pl(k)
  700      continue
  800   continue
  900 continue
c .......................................................................
c ... Comunica vetor de contagem de elementos por no'
      if (novlp) call allgatheri(ic,i_xfi)
c .......................................................................
      do 1000 i = 1, nnode
      do 1000 j = 1, ntn
         t(j,i) = t(j,i)/ic(i)
 1000 continue
c ......................................................................
      return
      end      
      subroutine t_efetive(ix,x,e,ie,ic,xl,ul,pl,tl,u,du,hi,tx,eps,t
     .         ,nnode,numel,nen,ndm,npi,ndf,nst,ntn,i_xfi,ilib,flaghidr,
     .                                         plastic)
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
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'termprop.fi'
      integer nnode,numel,nen,ndf,nst,ndm,ntn,npi
c ... ponteiro      
      integer*8 i_xfi
      logical flaghidr
c ......................................................................      
      integer ix(nen+1,*),ie(*),ic(*),dum
      integer nel,ma,iel,i,j,k,no
      integer ilib
      real*8 plastic(2,*)
      real*8  x(ndm,*),e(prop,*),xl(ndm,*),ul(ndf,*),pl(*),tl(*)
      real*8  utl(nen),u(ndf,*),du(*),hi(3,*),tx(npi,*)
      real*8 t(ntn,*)
      real*8  eps(npi,*)
c ......................................................................
c      do 20 i = 1, nnode
c         ic(i) = 0
c         do 10 j = 1, ntn
c            t(j,i) = 0.d0
c   10    continue
c   20 continue
c ......................................................................
c .... loop on elements:
      do 900 nel = 1, numel
c ...... set up local arrays
        do 600 i = 1, nen
          no = ix(i,nel)
          if (no .gt. 0) go to 300
          do 100 j = 1, ndm
            xl(j,i) = 0.d0
  100     continue
          do 200 j = 1, ndf
            ul(j,i)  = 0.d0
  200     continue
          go to 600
  300     continue
          do 400 j = 1, ndm
            xl(j,i) = x(j,no)
  400     continue
          do 500 j = 1, ndf
            ul(j,i) = u(j,no)
  500     continue
          tl(i) = du(no)
  600   continue
c ...... form element array
        ma  = ix(nen+1,nel)
c        ma  = ix(21,nel)
        iel = ie(ma)      
c ...... Biblioteca de elementos:
c        if (iel .eq. 11 .or. iel .eq. 12)  then
          call elmlib(e(1,ma),dum,xl,ul,pl,pl,pl,pl,pl,tl,utl,hi(1,nel),
     .    pl,tx(1,nel),tx(1,nel),eps(1,nel),ndm,nst,nel,iel,5,ma,1,ilib,
     .    flaghidr,plastic(1,nel))
c        else
c           plastic(1,nel) = 0
c           plastic(2,nel) = 0
c        endif
  900 continue
c .......................................................................
      return
      end 
c ***********************************************************************
      subroutine printadal(ad,al,apul,ja,japu,neq,nad,nadpu)
      implicit none
      integer i,ja(*),japu(*),neq,nad,nadpu
      real*8 ad(*),al(*),apul(*)
c ...  
      open(14,file='adCsr.dat',action='write')
c      
      write(14,*)'ad'
      do i = 1, neq
        write(14,*),i,ad(i)          
       enddo
c       
      write(14,*)'al'
      do i = 1, nad
        write(14,*),i,ja(i),al(i)   
      enddo
c      
      write(14,*)'alpu'
      do i = 1, nadpu 
        write(14,*),i,japu(i),apul(i)   
      enddo
      close(14)
c .....................................................................      
      return
      end
      subroutine printAx(sl,b,neq)
      implicit none
      integer i,j,neq
      real*8 sl(neq,*),x(100),b(*)
      do i = 1 , 68 
        x(i) = 0.d0
        do j = 1 , 68 
          if(j .ge. 61) then
            x(i) = x(i) - sl(i,j)*b(j)
          else
            x(i) = x(i) + sl(j,i)*b(j)
          endif
        enddo
      enddo
      do i = 1 , 68 
        print*,i,x(i)
      enddo 
      return
      end
