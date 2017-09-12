      subroutine partdmesh(i_np,i_ep,i_ix
     .                    ,i_elGL,i_elLG,i_noGL,i_noLG
     .                    ,i_fmap,i_dspl,i_rcvs,i_nnof,i_ranks,i_sizes  
     .                    ,i_elmdist,numel,nnode,maxno
     .                    ,my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
     .                    ,my_nnode,my_numel,my_numel_nov,my_numel_ov
     .                    ,nprcs,npart,mype,npes,MPIW,ovlp,partmesh,k)
c *********************************************************************
c *                                                                   *
c *  PARTDMESH : Gera o mapa da particao npart                        *
c *  ---------                                                        *
c *                                                                   *
c *  Parametros de entrada:                                           *
c *  ----------------------                                           *
c *                                                                   *
c *  i_np                                                             *
c *  i_ep                                                             *
c *  i_ix                                                             *
c *  numel                                                            *
c *  nnode                                                            *
c *  maxno                                                            *
c *  nnode                                                            *
c *  npart                                                            *
c *  partmesh                                                         *
c *  Parametros de saida:                                             *
c *  --------------------                                             *
c *********************************************************************
      use Malloc
      implicit none
c ... Do particionamento
      integer*8 i_nincid,i_incid,i_np,i_ep,i_Vi
      integer*8 i_elGL,i_elLG,i_noGL,i_noLG
      integer*8 i_fmap,i_dspl,i_rcvs,i_nnof,i_ranks,i_sizes
      integer*8 i_g,i_l
      integer my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
      integer my_nnode,my_numel,my_numel_nov,my_numel_ov
      integer nprcs
      integer maxgrade
      logical ovlp
c ... variaveis da malha
      integer numel,nnode,maxno
      integer*8 i_ix
c ... variaveis do Mpi   
      integer npart,npes,mype,MPIW,i,k
      integer*8 i_elmdist
      logical partmesh
c ... variaveis auxiliares
c      real*8 use_mem      
c =====================================================================
c
c ===
      if(partmesh)then
        i_sizes  = dealloc('sizes    ')
        i_ranks  = dealloc('ranks    ')
        i_dspl   = dealloc('dspl     ')
        i_rcvs   = dealloc('rcvs     ')
        i_nnof   = dealloc('nnof     ')
        i_fmap   = dealloc('fmap     ')
        i_noLG   = dealloc('noLG     ')
        i_noGl   = dealloc('noGL     ')
        i_elLG   = dealloc('elLG     ')
        i_elGL   = dealloc('elGL     ')
      endif
c =====================================================================
c
c === 
      if(npes .gt. 1)then
        i_g   = alloc_4('g        ',maxno+1,numel)
        call mzero(ia(i_g),numel*(maxno+1))
        call global_int(ia(i_g),ia(i_ix),ia(i_elmdist),numel,maxno+1
     .                 ,npes   ,mype    ,MPIW)  
        if( k .lt. 1 ) then
          i_g      = dealloc('g       ')
          return
        endif
        i_l  = i_ix
        i_ix = i_g
      endif
c =====================================================================
c
c ====  contando o numero de elemenos do meu processo
c     print*,"contando numero de elementos por particao ..."
      call nelmpart(ia(i_np),ia(i_ep),ia(i_ix),numel,maxno
     .             ,my_numel,my_numel_nov,my_numel_ov,ovlp,npart)
c =====================================================================
c ... Correspondencia entre numeracao de elementos local e global:
c     elGL(i) = numero local do elemento global i
c
      i_elGL = alloc_4('elGL    ',  1, numel)
      call mzero (ia(i_elGL),numel)
c
c ... Correspondencia entre numeracao de elementos local e global:
c     elLG(i) = numero global do elemento local i
      i_elLG = alloc_4('elLG    ', 1, my_numel)
      call mzero (ia(i_elLG),my_numel)
c ===  gerando as propriedades nodais
c
c === propnode
      call propnode(ia(i_ix),i_nincid,i_incid,nnode,numel,maxno,maxgrade
     .             ,.true.)
c =====================================================================
c
c =====================================================================
c ...
      call getmapelm(ia(i_np),ia(i_ep),ia(i_ix),ia(i_elLG),ia(i_elGL)
     .              ,numel,my_numel,my_numel_ov,my_numel_nov,maxno
     .              ,ovlp,npart) 
c .....................................................................
c
c =====================================================================
c =====================================================================
c
c ====  nos
c      
c ... Conjunto de nos locais:
      i_Vi = alloc_4('Vi      ', 1,   nnode)
      call mzero(ia(i_Vi),  nnode)
c .....................................................................
c
c ... contando o numero de elemenos do meu processo
      call nnodepart(ia(i_np),ia(i_ep),ia(i_ix)
     .              ,ia(i_incid),ia(i_nincid),ia(i_Vi),ia(i_elLG)
     .              ,nnode,maxno,maxgrade,my_numel,my_numel_nov,my_nno1
     .              ,my_nno2,my_nno3,my_nno4,my_nnode,ovlp,npart)
c     
c ... Correspondencia entre numeracao nodal global -> local:
c     noGL(i) = numero local do no global i
      i_noGl = alloc_4('noGL    ',  1,   nnode)
      call mzero(ia(i_noGL),nnode)
c
c ... Correspondencia enter numeracao nodal local e global:
c     noLG(i) = numero global do no local i
c
      i_noLG = alloc_4('noLG    ',  1,my_nnode)
      call mzero(ia(i_noLG),my_nnode)
c === mapa dos nos
c     print*,"mapa de nos..."
      call getmapnod(ia(i_np),ia(i_ep),ia(i_ix),ia(i_noLG),ia(i_noGL)
     .              ,ia(i_incid),ia(i_nincid),ia(i_Vi),ia(i_elLG)
     .              ,nnode,maxno,maxgrade,my_numel,my_numel_nov,my_nno1
     .              ,my_nno2,my_nno3,my_nno4,my_nno1a,my_nnode
     .              ,ovlp,npart)
c ... libera memoria     
      i_Vi     = dealloc('Vi      ')
c ... recupera os ponteiros      
      i_noLG   = locate('noLG     ')
      i_noGl   = locate('noGL     ')
c      i_incid  = locate('incid   ')
c      i_nincid = locate('nincid  ')
c .....................................................................
c 
c ... Correspondencia entre numeracao nodal de interface e global:
      i_fmap  = alloc_4('fmap    ', 1, 2*my_nnode )
      call mzero(ia(i_fmap),2*my_nnode)
      i_nnof  = alloc_4('nnof    ', 1,       2 )
      call mzero(ia(i_nnof),2)
c .....................................................................
c
c ... Quantidade e lista de vizinhos:
c
      i_rcvs = alloc_4('rcvs    ', 1, 2*nprcs)
      call tes_mem(i_rcvs,'i_rcvs  ')
      call mzero(ia(i_rcvs),2*nprcs)
      i_dspl = alloc_4('dspl    ', 1, 2*nprcs)
      call tes_mem(i_dspl,'dspl    ')
      call mzero(ia(i_dspl),2*nprcs)
c ......................................................................
c
c ... Correspondencia entre numeracao nodal de interface e global:
c
      i_ranks  = alloc_4('ranks   ',       1, 2*nprcs)
      call tes_mem(i_ranks,'ranks   ')
      call mzero(ia(i_ranks),2*nprcs)
      i_sizes  = alloc_4('sizes   ',       1,       2)
      call mzero(ia(i_sizes),2)
      call tes_mem(i_sizes,'sizes   ')
c .....................................................................
c
c
c =====================================================================
c
c ===
       call getmapviz(ia(i_np),ia(i_ep),ia(i_ix),ia(i_noLG)
     .              ,ia(i_elLG)
     .              ,ia(i_incid),ia(i_nincid),ia(i_ranks)
     .              ,ia(i_sizes)
     .              ,my_numel,my_numel_ov
     .              ,my_nnode,my_nno1,my_nno2,my_nno3,my_nno4
     .              ,maxno
     .              ,maxgrade
     .              ,ovlp,nprcs,npart) 
c .....................................................................      
c     print*,'getfmap'
      call getfmap(ia(i_np),ia(i_ep),ia(i_ix),ia(i_noLG)
     .            ,ia(i_incid),ia(i_nincid),ia(i_ranks)
     .            ,ia(i_sizes),ia(i_fmap),ia(i_nnof)
     .            ,ia(i_rcvs),ia(i_dspl),nnode  
     .            ,my_nnode,my_nno1,my_nno2,my_nno4
     .            ,maxgrade,maxno
     .            ,ovlp,nprcs,npart) 
c .....................................................................
c
c =====================================================================
c
c === libera a memoria
      i_incid =  dealloc('incid   ')
      i_nincid = dealloc('nincid  ')
      if(npes .gt. 1)then
        i_g      = dealloc('g       ')
        i_ix     = i_l 
      endif
c ... recalcula os ponteiros      
      i_np     = locate('np       ') 
      i_ep     = locate('ep       ')
      i_sizes  = locate('sizes    ')
      i_ranks  = locate('ranks    ')
      i_dspl   = locate('dspl     ')
      i_rcvs   = locate('rcvs     ')
      i_nnof   = locate('nnof     ')
      i_fmap   = locate('fmap     ')
      i_noLG   = locate('noLG     ')
      i_noGl   = locate('noGL     ')
      i_elLG   = locate('elLG     ')
      i_elGL   = locate('elGL     ')
c =====================================================================
c
c ===
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine nelmpart(np,ep,ix,numel,maxno,my_numel
     .                   ,my_numel_nov,my_numel_ov,ovlp,npart)
c *********************************************************************
c *                                                                   *
c *  NELPART :                                                        *
c *  -------                                                          *
c *                                                                   *
c *  Parametros de entrada:                                           *
c *  ----------------------                                           *
c *    np           - particinamento dos nos                          *
c *    ep           - particionamento dos elementos                   *
c *    ix           - conectiva dos elementos                         *
c *    numel        - numero de elementos                             *
c *    maxno        - numero maximo de nos por elementos              *
c *    ovlp         - indicador de overlapping                        *
c *                                                                   *
c *  Parametros de saida:                                             *
c *  --------------------                                             *
c *                                                                   *
c *    my_numel     - numero de elementos  locais                     *
c *    my_numel_nov - numero de elementos  locais non-overlaping      *
c *    my_numel_ov  - numero de elementos  locais overlaping          *
c *                                                                   *
c *********************************************************************
      implicit none
      integer numel,maxno,npart
      integer np(*),ep(*),ix(maxno+1,*)
      integer my_numel,my_numel_nov,my_numel_ov
      integer nel,nno,ii
      logical ovlp
      integer aux,cont
c =====================================================================
c
c === non-overlapping
      my_numel_nov = 0
      do nel = 1,numel
        if( ( npart + 1 ) .eq. ep( nel ) )then
          my_numel_nov = my_numel_nov + 1
        endif
      enddo 
c =====================================================================
c
c === overlapping
      my_numel_ov = 0
      if( ovlp ) then
        cont = 0
        do nel = 1, numel
          if( ( npart + 1 ) .ne. ep(nel) ) then
            aux = 0
              do nno = 1, maxno
                ii   = ix( nno , nel )
                if( ii .gt. 0 ) then
                  if( ( npart + 1 ) .eq. np(ii) ) then
                    aux = 1
                  endif
                endif
              enddo
              if( aux .eq. 1 ) then 
                cont = cont + 1 
              endif
          endif
        enddo
        my_numel_ov = cont
      endif
c =====================================================================
c
c ===
      my_numel = my_numel_nov + my_numel_ov
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c *                                                                   *
c *  GETMAPELM : Gera o mapa local dos elementos                      *
c *  ---------                                                        *
c *                                                                   *
c *  Parametros de entrada:                                           *
c *  ----------------------                                           *
c *    np           - particinamento dos nos                          *
c *    ep           - particionamento dos elementos                   *
c *    ix           - conectivadades  dos elementos                   *
c *    elLG         -                                                 *
c *    elGL         -                                                 *
c *    nnode        - numero de nos                                   *
c *    numel        - numero de elementos                             *
c *    my_numel     - numero de elementos locais                      *
c *    my_numel_nov - numero de elementos locais non-overlaping       *
c *    my_numel_no  - numero de elementos locais overlaping           *
c *    maxno        - numero maximo de nos por elementos              *
c *    ovlp         - indicador de overlapping                        *
c *    npart        - meu rank                                        *
c *                                                                   *
c *  Parametros de saida:                                             *
c *  --------------------                                             *
c *                                                                   *
c *    elLG         - mapa de elementos local->global                 *
c *    elGL         - mapa de elementos global0>local                 *
c *                                                                   *
c *********************************************************************
      subroutine getmapelm(np,ep,ix,elLG,elGL,numel,my_numel,my_numel_ov
     .                    ,my_numel_nov,maxno,ovlp,npart) 
c === Variaveis
      implicit none
c ... grahp metis
      integer np(*),ep(*)
c ... malha
      integer numel,maxno,ii
      integer ix(maxno+1,*)
c ... particionamento
      logical ovlp
      integer my_numel,my_numel_ov,my_numel_nov
      integer elGL(*),elLG(*)
      logical incluir
c ... MPi
      integer npart
c ... auxiliares
      integer cont,nel,part,no
c =====================================================================
c
c _____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                NUMERAÇÃO LOCAL DE ELEMENTOS NON-OVERLAPPING         |
c|                                                                     |
c|_____________________________________________________________________|
c
c === mapa de elementos local -> global
      cont = 1
      do nel = 1, numel
        if( (npart + 1) .eq. ep(nel) ) then
          elLG(cont) = nel
          cont = cont + 1
        endif
      enddo
c =====================================================================      
c _____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                NUMERAÇÃO LOCAL DE ELEMENTOS OVERLAPPING             |
c|                                                                     |
c|_____________________________________________________________________|
      if (ovlp) then
c ... Lopp global nos elementos        
        do nel = 1, numel
c ... Teste se o elemento é overlapping
          part = ep(nel)
          if( part .ne. (npart + 1))then
            incluir = .false.
            do no = 1, maxno
              ii = ix(no,nel)
              if( ii .gt. 0 ) then
                if(np(ii) .eq. (npart +1)) incluir =.true.
              endif
            enddo
            if(incluir)then
              elLG(cont) = nel
              cont = cont + 1
            endif
          endif
        enddo
      endif  
c
c === reordena elLG
      call ord(            elLG(1),my_numel_nov)
      call ord(  elLG(my_numel_nov+1), my_numel_ov)
c    com overlapping      
c =====================================================================
c _____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                    MAPA GLOBAL-LOCAL DE ELEMENTOS                   |
c|                                                                     |
c|_____________________________________________________________________|
c
c ===
c     cont =  0
c     do nel = 1, my_numel
c       cont = cont + 1
c       elGL( elLG(nel) ) = cont
c     enddo  
c =====================================================================
c
c ===
      return
      end
c =====================================================================      
c *********************************************************************
c
c *********************************************************************
c *                                                                   *
c *  NNODEPART : conta o numero de nos V1,V2,V3,V4                    *
c *  ---------                                                        *
c *                                                                   *
c *  Parametros de entrada:                                           *
c *  ----------------------                                           *
c *    np           - particinamento dos nos                          *
c *    ep           - particionamento dos elementos                   *
c *    el           - conectiva dos elementos                         *
c *    inicid       - elemento por                                    *
c *    nincid       - numero de aresta por no                         *
c *    elLG         - map de elemento local - global                  *
c *    Vi           - grupo de nos (V1,V2,V3,V4)                      *
c *    nnode        - numero de nos da malha                          *
c *    my_nnode     -                                                 *
c *    my_nno1      -                                                 *
c *    my_nno2      -                                                 *
c *    my_nno3      -                                                 *
c *    my_nno4      -                                                 *
c *    my_numel     - numero de elementos locais                      *
c *    my_numel_nov - numero de elementos locais non-overllaping      *
c *    maxgrade     - numero maximo de conectividade da malha         *
c *    maxno        - numeo maximo de nos por elemento                *
c *    ovlp         - indicador de overlapping                        *
c *    nprcs        - numero de processo                              *
c *    npart        - meu rank                                        *
c *                                                                   *
c *  Parametros de saida:                                             *
c *  --------------------                                             *
c *    my_nnode     - numero da malha local                           *
c *    my_nno1      - numero de nos V1                                *
c *    my_nno2      - numero de nos V2                                *
c *    my_nno3      - numero de nos V3                                *
c *    my_nno4      - numero de nos V4                                *
c *                                                                   *
c *********************************************************************
      subroutine nnodepart(np,ep,el,incid,nincid,Vi,elLG
     .                    ,nnode,maxno,maxgrade
     .                    ,my_numel,my_numel_nov,my_nno1,my_nno2
     .                    ,my_nno3,my_nno4,my_nnode,ovlp,npart) 
c === Variaveis
      implicit none
c ... grahp metis
      integer np(*),ep(*)
c ... malha
      integer nnode,maxno
      integer el(maxno+1,*)
c ... particionamento
      logical ovlp
      integer  my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
      integer my_numel,my_numel_nov,my_nnode,maxgrade
      integer elLG(*),incid(maxgrade,*)
      integer nincid(*),Vi(*)
c ... MPi
      integer npart
c ... auxiliares
      integer nno,ninc,incel,k,i,j,kk,cont 
c =====================================================================
c _____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                         IDENTIFICA NÓS V1                           |
c|                                                                     |
c|_____________________________________________________________________|
c
c ===
      my_nno1 = 0
      do nno = 1, nnode
        k = 0
        ninc = nincid(nno)
        do j = 1, ninc
          incel = incid(j,nno)
          if ( ep( incel ) .eq. (npart + 1) ) k = k + 1
        enddo
        if( k .eq. ninc ) then
          Vi(nno) = 1
          my_nno1 = my_nno1 + 1
        endif
      enddo
c =====================================================================      
c _____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                         IDENTIFICA NÓS V2 e V3                      |
c|                                                                     |
c|_____________________________________________________________________|
c
c ===
      my_nno2 = 0
      my_nno3 = 0
      do i = 1, my_numel_nov
        do j = 1, maxno
          nno  = el(j,elLG(i))
          if( nno .gt. 0 ) then
            kk   = np(nno)
c ... Se o no não pertence a pratição, então e V3
            if( kk .ne. (npart + 1) )then
              if(Vi(nno) .ne. 3)then
                Vi( nno ) = 3
                my_nno3 = my_nno3 + 1
              endif  
c ... se o no pertence a particao e não é v1 e v2          
            else
              if( (Vi(nno) .ne. 1) .and. (Vi(nno) .ne. 2))then
                Vi( nno ) = 2
                my_nno2 = my_nno2 + 1 
              endif  
            endif
          endif  
        enddo
      enddo
c =====================================================================      
c _____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                          IDENTIFICA NÓS V4                          |
c|                                                                     |
c|_____________________________________________________________________|
c
c === 
      my_nno4  = 0
      if(ovlp)then
c ... Loop local nos elementos overlapping da minha partição
        do i = my_numel_nov + 1,my_numel
c
c ... Loop nos nós do elemento
          do j = 1,maxno
            nno = el(j,elLG(i))
            if( nno .gt. 0 ) then
c ... Se o nó não for V2 nem V3, ele é V4
              kk= Vi(nno)
              if( (kk.ne.2) .and. (kk.ne.3) .and. (kk.ne.4))then
                Vi(nno) = 4
                my_nno4 = my_nno4 + 1
              endif
            endif
          enddo
        enddo
      endif
c =====================================================================
      my_nnode = my_nno1 + my_nno2 + my_nno3 + my_nno4
c === 
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c *                                                                   *
c *  GETMAPNODE                                                       *
c *  -------                                                          *
c *                                                                   *
c *  Parametros de entrada:                                           *
c *  ----------------------                                           *
c *    np           - particinamento dos nos                          *
c *    ep           - particionamento dos elementos                   *
c *    el           - conectiva dos elementos                         *
c *    inicid       - elemento por                                    *
c *    nincid       - numero de aresta por no                         *
c *    noLG         -                                                 *
c *    noGL         -                                                 *
c *    elLG         - map de elemento local - global                  *
c *    Vi           - grupo de nos (V1,V2,V3,V4)                      *
c *    nnode        - numero de nos da malha                          *
c *    my_nnode     - numero de nos da malha local                    *
c *    my_nno1      - numero de nos V1                                *
c *    my_nno2      - numero de nos V2                                *
c *    my_nno3      - numero de nos V3                                *
c *    my_nno4      - numero de nos V4                                *
c *    my_nno1a     -                                                 *
c *    my_numel     - numero de elementos locais                      *
c *    my_numel_nov - numero de elementos locais non-overllaping      *
c *    maxgrade     - numero maximo de conectividade da malha         *
c *    maxno        - numero maximo de nos por elemento               *
c *    ovlp         - indicador de overlapping                        *
c *    nprcs        - numero de processo                              *
c *    npart        - meu rank                                        *
c *                                                                   *
c *  Parametros de saida:                                             *
c *  --------------------                                             *
c *    noLG         - mapa local>global de nos                        *
c *    noGL         - mapa global>local de nos                        *
c *    my_nno1a     - numero de nos V1a                               *
c *                                                                   *
c *********************************************************************
      subroutine getmapnod(np,ep,el,noLG,noGL,incid,nincid,Vi,elLG
     .              ,nnode,maxno,maxgrade
     .              ,my_numel,my_numel_nov,my_nno1,my_nno2
     .              ,my_nno3,my_nno4,my_nno1a,my_nnode
     .              ,ovlp,npart) 
c === Variaveis
      implicit none
c ... grahp metis
      integer np(*),ep(*)
c ... malha
      integer nnode,maxno
      integer el(maxno+1,*)
c ... particionamento
      logical ovlp
      integer  my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
      integer my_numel,my_numel_nov,my_nnode,maxgrade
      integer elLG(*),incid(maxgrade,*)
      integer nincid(*),Vi(*),nolG(*),noGl(*)
c ... MPi
      integer npart
c ... auxiliares
      integer nno,ninc,incel,k,i,j,kk,cont 
c _____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                NUMERAÇÃO LOCAL DE NÓS NON-OVERLAPPING               |
c|                                                                     |
c|_____________________________________________________________________|
c
c ... V1
c
      cont = 0
      do nno = 1, nnode
        if(Vi(nno) .eq. 1) then
          cont = cont + 1
          noLG(cont) = nno
        endif
      enddo
c 
c ... V2
c
      do nno = 1, nnode
        if(Vi(nno) .eq. 2) then
          cont = cont + 1
          noLG(cont) = nno
        endif
      enddo
c
c ... V3
c
      do nno = 1, nnode
        if(Vi(nno) .eq. 3) then
          cont = cont + 1
          noLG(cont) = nno
        endif
      enddo
c _____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                  NUMERAÇÃO LOCAL DE NÓS OVERLAPPING                 |
c|                                                                     |
c|_____________________________________________________________________|
      my_nno1a = 0
      if(ovlp) then
c
c ... V4
c
        do nno = 1, nnode
          if(Vi(nno) .eq. 4)then
            cont = cont + 1
            noLG(cont) = nno
          endif
        enddo
c .....................................................................        
c
c ... Agrupa V1a ao final de V1
c
c        
c .. identifica os nos nno1a como V5        
c
        do i= my_nno1 + my_nno2 + 1, my_nno1+my_nno2+my_nno3
          nno = noLG(i)
          ninc  = nincid(nno)
          do j = 1, ninc
            incel = incid(j,nno)
            if (ep(incel) .eq. (npart+1))then
              do k = 1, maxno
                kk   = el(k,incel)
                if( kk .gt. 0 ) then
c ... Se no for V1 ele é nno1a              
                  if( Vi(kk) .eq. 1 )then
                    Vi(kk)   = 5
                    my_nno1a = my_nno1a + 1
                  endif
                endif
              enddo  
            endif
          enddo
        enddo
c .....................................................................        
        kk = 0
         do i = my_nno1,1,-1
           nno = noLG(i)
           k   = Vi(nno)
           if( k .eq. 5) then
             j                  = noLG(my_nno1-kk)
             noLG(my_nno1 - kk) = noLG(i)
             noLG(i)            = j 
             kk                 = kk + 1
           endif
         enddo
c .....................................................................           
      endif   
c ____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                          REORDENA noLG                              |
c|                                                                     |
c|_____________________________________________________________________|
c
c ===
      call ord(noLG( 1 )                              ,my_nno1-my_nno1a)
      call ord(noLG( my_nno1 - my_nno1a + 1)          ,        my_nno1a)
      call ord(noLG( my_nno1 + 1)                     ,         my_nno2)
      call ord(noLG( my_nno1 + my_nno2 + 1 )          ,         my_nno3)
      call ord(noLG( my_nno1 + my_nno2 + my_nno3 + 1 ),         my_nno4)
c ======================================================================
c _____________________________________________________________________
c|                                                                     |
c|                                                                     |
c|                       MAPA GLOBAL-LOCAL DE NÓS                      |
c|                                                                     |
c|_____________________________________________________________________|
c
c ===
       call mzero(noGl,nnode)
       do nno = 1 , my_nnode
         noGL(noLG(nno)) = nno
       enddo
c =====================================================================
c
c === 
      return
      end
c =====================================================================
c *********************************************************************
c
c
c *********************************************************************
c *                                                                   *
c *  GETMAPVIZ                                                        *
c *  -------                                                          *
c *                                                                   *
c *  Parametros de entrada:                                           *
c *  ----------------------                                           *
c *    np           - particinamento dos nos                          *
c *    ep           - particionamento dos elementos                   *
c *    el           - conectiva dos elementos                         *
c *    inicid       - elemento por                                    *
c *    nincid       - numero de aresta por no                         *
c *    ranks        -                                                 *
c *    sizes        -                                                 *
c *    nnode        - numero de nos da malha                          *
c *    my_nnode     - numero da malha local                           *
c *    my_nno1      - numero de nos V1                                *
c *    my_nno2      - numero de nos V2                                *
c *    my_nno3      - numero de nos V3                                *
c *    my_nno4      - numero de nos V4                                *
c *    my_numel     - numero de elementos locais                      *
c *    my_numel_nov - numero de elementos locais non-overllaping      *
c *    my_numel_ov  - numero de elementos locais overlaping           *
c *    maxgrade     - numero maximo de conectividade da malha         *
c *    maxno        - numero maximo de nos por elemento               *
c *    ovlp         - indicador de overlapping                        *
c *    nprcs        - numero de processo                              *
c *    npart        - meu rank                                        *
c *                                                                   *
c *  Parametros de saida:                                             *
c *  --------------------                                             *
c *    ranks(2*nprcs)- vizinhos                                       *
c *    sizes(2)      - numero de vizinhos                             *
c *                                                                   *
c *********************************************************************
       subroutine getmapviz(np,ep,el,noLG
     .              ,elLG
     .              ,incid,nincid,ranks
     .              ,sizes
     .              ,my_numel,my_numel_ov
     .              ,my_nnode,my_nno1,my_nno2,my_nno3,my_nno4
     .              ,maxno
     .              ,maxgrade
     .              ,ovlp,nprcs,npart) 
c === Variaveis
c ... Particionamento
      implicit none
      integer np(*),ep(*),noLG(*),elLG(*)
      integer incid(maxgrade,*),nincid(*)
      integer ranks(*)
      integer sizes(*)
      integer maxgrade
      integer nprcs
      logical ovlp
c ... malha local      
      integer my_nnode,my_nno1,my_nno2,my_nno3,my_nno4
      integer my_numel_ov,my_numel
c ... malha
      integer el(maxno+1,*),maxno
c ... Mpi
      integer npart
c ... auxiliar
      logical contviz(128)
      integer i,l,part,nel,j,nno,kk,k,incel
c =====================================================================      
c
c ===
c ===
c ... inicializando contviz
      do i = 1, nprcs
        contviz(i) = .false.
      enddo
c .....................................................................      
c
c ... overlaping
      if(ovlp) then
c ... recv     
c ... loop nos V3 e V4
        do i = my_nno1 + my_nno2 + 1, my_nnode 
          nno = noLG(i)
          part = np(nno)
c ... a particao q tiver esse no e minha vizinha(recv)  
          if(contviz(part) .eqv. .false.) contviz(part) = .true.
        enddo
c ... gerando os vizinhos
        i = 1
        do part = 1,nprcs
          if(contviz(part))then
            ranks(i) = part - 1 
            i = i + 1
          endif
        enddo
        sizes(1) = i - 1
c ... send      
        do i = 1, nprcs
          contviz(i) = .false.
        enddo
c ... loop nos elementos em overlaping
        do i = my_numel_ov + 1, my_numel
          nel = elLG(i)
c ... loop nos nos desse elemento
          do j = 1 ,maxno
            kk = el(j,nel)
            if( kk .gt. 0 ) then
c ... se tiver um no dessa particao teremos um elemnto overleping
c     desse no
              if(ep(kk).ne.np(nel).and.np(nel).ne.(npart+1))then
c ... para na contar vizinhos repetidos          
                if(contviz(part) .eqv. .false.) contviz(part) = .true.
              endif
            endif
          enddo
        enddo
c ... loop nos nos V2 normal
        do i = my_nno1 + 1, my_nno1 + my_nno2
          nno = noLG(i)
          l   = nincid(nno)
          do j = 1 , l
            nel   = incid(j,nno)
            part = ep(nel)
            if( part .ne. (npart+1))then
c ... para na contar vizinhos repetidos          
              if(contviz(part) .eqv. .false.) contviz(part) = .true.
            endif
          enddo
        enddo
c ... loop nos nos V2 em overlaping com outra particao 
        do i = my_nno1 + 1, my_nno1 + my_nno2
          nno = noLG(i)
          l   = nincid(nno)
          do j = 1 , l
            nel   = incid(j,nno)
            part = ep(nel)
            if( part .ne. (npart+1))then
              do k = 1, maxno
                kk = el(k,nel)
                if ( kk .gt. 0) then
                  if((part.ne.np(kk)).and.(np(kk).ne.(npart+1)))then
                    part = np(kk)
c... para na contar vizinhos repetidos 
                    if(contviz(part) .eqv. .false.)contviz(part)=.true.
                  endif
                endif
              enddo
            endif
          enddo
        enddo
c ... loop nos nos V3, pois eles geram os nos V1a
        do i = my_nno1 + my_nno2 + 1 , my_nno1 + my_nno2 + my_nno3
          nno = noLG(i)
          l = nincid(nno)
          do j = 1, l
            incel = incid(j,nno)
c ... ha possibilidade de ser um elemento de overlapping  
            if(ep(incel) .eq. (npart+1)) then
              do k = 1, maxno
                kk = el(k,incel)
                if( kk .gt. 0 ) then
                  if(np(kk) .eq.(npart+1)) then
                    part = np(nno)
c ... para na contar vizinhos repetidos
                    if(contviz(part) .eqv. .false.)contviz(part)=.true.
                  endif
                endif
              enddo
            endif
          enddo  
        enddo
c ... gerando os vizinhos
      i = 1
      do part = 1,nprcs
        if(contviz(part))then
          ranks(nprcs+i) = part - 1 
          i = i + 1
        endif
      enddo
      sizes(2) = i - 1
c .....................................................................
c
c ... non-overlaping      
      else
c ... loop nos nos V2 , V3
        do i = my_nno1 + 1, my_nnode - my_nno4
          nno = noLG(i)
          l   = nincid(nno)
          do j = 1 , l
            nel   = incid(j,nno)
            part = ep(nel)
            if( part .ne. (npart+1))then
c ... para na contar vizinhos repetidos          
              if(contviz(part) .eqv. .false.) contviz(part) = .true.
            endif
          enddo
        enddo
c ... gerando os vizinhos
        i = 1
        do part = 1,nprcs
          if(contviz(part))then
            ranks(i) = part - 1 
            i = i + 1
          endif
        enddo
        sizes(1) = i - 1
c .....................................................................
c
      endif
c ...
      call ord(ranks,sizes(1))
c .....................................................................
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c *                                                                   *
c *  GETFMAP                                                          *
c *  -------                                                          *
c *                                                                   *
c *  Parametros de entrada:                                           *
c *  ----------------------                                           *
c *    np           - particinamento dos nos                          *
c *    ep           - particionamento dos elementos                   *
c *    el           - conectiva dos elementos                         *
c *    inicid       - elemento por                                    *
c *    nincid       - numero de aresta por no                         *
c *    ranks        - vizinhos                                        *
c *    sizes        - numeros vizinhos                                *
c *    fmap         -                                                 *
c *    nnof         -                                                 *
c *    rcvs         -                                                 *
c *    dspl         -                                                 *
c *    nnode        - numero de nos da malha                          *
c *    my_nnode     - numero da malha local                           *
c *    my_nno1      - numero de nos V1                                *
c *    my_nno2      - numero de nos V2                                *
c *    my_nno3      - numero de nos V3                                *
c *    my_nno4      - numero de nos V4                                *
c *    my_nno1a     - numero de nos V1a                               *
c *    my_numel     - numero de elementos locais                      *
c *    maxgrade     - numero maximo de conectividade da malha         *
c *    maxno        - numero maximo de nos por elemento               *
c *    ovlp         - indicador de overlapping                        *
c *    nprcs        - numero de processo                              *
c *    npart        - meu rank                                        *
c *                                                                   *
c *  Parametros de saida:                                             *
c *  --------------------                                             *
c *    fmap         - nos de comunica por tipo ordenado por vizinho   *
c *    nnof         - numero de nos de comunicacao                    *
c *    rcvs         -                                                 *
c *    dspl         -                                                 *
c *                                                                   *
c *********************************************************************
      subroutine getfmap(np,ep,el,noLG
     .                  ,incid,nincid,ranks,sizes
     .                  ,fmap,nnof,rcvs,dspl,nnode
     .                  ,my_nnode,my_nno1,my_nno2,my_nno4
     .                  ,maxgrade,maxno
     .                  ,ovlp,nprcs,npart) 
c === Variaveis
      implicit none
c ...particionamento      
      integer np(*),ep(*),noLG(*)
      integer incid(maxgrade,*),nincid(*)
      integer ranks(*)
      integer sizes(*)
      integer maxgrade
      integer nprcs
      logical ovlp
      integer fmap(*),nnof(*)
      integer rcvs(*),dspl(*)
c ... malha local      
      integer my_nnode,my_nno1,my_nno2,my_nno4
c ... malha      
      integer nnode,maxno
      integer el(maxno+1,*)
c ... Mpi
      integer npart
c ... auxiliar      
      integer i,l,part,nel,j,nno,kk,k,m,viz,ii,jj
      logical nc
c =====================================================================
c
c ===
c      do i = 1, 2*nnode
c        print*,"V",i,Vi(i)
c      enddo
      if(ovlp)then
c ... rcvs
        kk = 0
        do i = 1,sizes(1)
            m = 0
          viz = ranks(i) + 1
c ... loop nos V3 e V4
          do j = my_nno1 + my_nno2 + 1, my_nnode
            nno = noLG(j)
c ... V3 e V4 so ha envio do viz q contem akele no
            nc = .true.
            if( np(nno) .eq. viz .and. nc) then
              m = m + 1
              fmap(kk+m) = nno
              nc = .false.
            endif
          enddo
          call ord(fmap(kk+1),m)
          kk   = kk + m
          rcvs(i) = m
        enddo
        nnof(1) = kk  
c ... Send
c ... loop vizinhos
        kk = 0
        do i = 1, sizes(2)
          m   = 0
          viz = ranks(nprcs+i) + 1
c ... loop nos V1a
          do j = my_nno1 - my_nno1 + 1, my_nno1
            nno = noLG(j)
            l   = nincid(nno)
            nc = .true.
            do k = 1, l
              nel  = incid(k,nno)
              part = ep(nel)
c ... se o elemento penterce a minha pid, procura os nos V3 dele  
              if(part .eq. (npart + 1)) then 
                do jj = 1, maxno
                  ii = el(jj,nel)
                  if( ii .gt. 0 ) then
                    if( np(ii) .eq. viz .and. nc)then
                      m  = m +1
                      fmap(my_nnode+kk+m) = nno
                      nc = .false.
                    endif
                  endif  
                enddo
              endif
            enddo  
          enddo 
c ... loop nos V2
          do j = my_nno1 + 1, my_nno1+my_nno2
            nno = noLG(j)
            l   = nincid(nno)
            nc  = .true.
            do k = 1, l
              nel  = incid(k,nno)
              part = ep(nel)
c ... no V2  em overlaping com nenhum dos nos do elemento 
c     pertecendo a part viz  
              if( (part .eq. viz) .and. nc )then
                m  = m + 1
                fmap(my_nnode+kk+m) = nno
                nc = .false.
              endif
c ... no V2  s/overlaping de outras particoes
              if( part .ne. (npart + 1) )then
                do jj = 1, maxno
                  ii = el(jj,nel)
                  if( ii .gt. 0 ) then
                    if( np(ii) .eq. viz .and. nc)then
                      m  = m +1
                      fmap(my_nnode+kk+m) = nno
                      nc = .false.
                    endif
                  endif
                enddo
              endif  
c ... no V2  c/overlaping de outras particoes
              if( part .eq. (npart + 1) )then
                do jj = 1, maxno
                  ii = el(jj,nel)
                  if( ii .gt. 0 ) then
                    if( np(ii) .eq. viz .and. nc)then
                      m  = m +1
                      fmap(my_nnode+kk+m) = nno
                      nc = .false.
                    endif
                  endif
                enddo
              endif
c .....................................................................        
            enddo
          enddo
          call ord(fmap(my_nnode+kk+1),m)
          kk   = kk + m
          rcvs(nprcs+i) = m
        enddo  
c .....................................................................        
c
c ...
        nnof(2) = kk
c .....................................................................        
c
c ... rcvs
        do i = 2, sizes(1)
          dspl(i) = dspl(i-1) + rcvs(i-1)
        enddo
c .....................................................................        
c
c ... send
        do i = 2, sizes(2)
          dspl(i+nprcs) = dspl(nprcs+i-1) + rcvs(nprcs+i-1)
        enddo
c .....................................................................      
c .....................................................................
c
c ... non-overlapping
      else 
        kk   = 0
c ... loop nos vizinhos      
        do i = 1, sizes(1)
          m    = 0
          viz  = ranks(i) + 1
c ... loop nos V2  e V3       
         do j = my_nno1 + 1,my_nnode-my_nno4
            nno = noLG(j)
c ... flag para nao contarmos o mesmo noh 2 vezes para um memos vizinho
c nc = .true.  no ainda nao contabilizado
c nc = .flase. no jah contabilizado
            nc = .true.
            l   = nincid(nno)
            do k = 1, l
              nel   = incid(k,nno)
              part = ep(nel)
              if(part .eq. viz .and. nc) then
                m       = m + 1
                fmap(kk+m) = nno
                nc = .false.
              endif
            enddo  
          enddo
c .....................................................................
c
          call ord(fmap(kk+1),m)
          kk      = kk   + m
          rcvs(i) = m
        enddo
c .....................................................................      
        nnof(1) = kk
        nnof(2) = kk
c ...
        do i = 2, sizes(1)
          dspl(i) = dspl(i-1) + rcvs(i-1)
        enddo
      endif
c =====================================================================
      return
      end
c=======================================================================     
c *********************************************************************
c
c *********************************************************************
      subroutine ord(a,n)
c *********************************************************************
c *                                                                    
c *   ORD: Re-ordena um vetor de inteiros com algoritmo 'bubble sort'  
c *                                                                    
c *   Parâmetros de entrada:                                           
c *   ----------------------                                           
c *      a(n) - vetor de inteiros                                      
c *      n    - cardinalidade de a                                     
c *                                                                    
c *   Parâmetros de saida:                                             
c *   -------------------                                              
c *      a(n) - vetor de inteiros reordenado                           
c *                                                                    
c *********************************************************************
      integer a(*),n,i,j,k,m,p
c ... bubble sort
   10 continue
      p=0
      do i=1,n-1
         j=i+1
         k=a(i)
         m=a(j)
         if(k.gt.m)then
            a(i)=m
            a(j)=k
            p=1
         endif
      enddo
      if(p.eq.1)goto 10
c ... dumb sort
c      do i=n,2,-1
c        m=a(i)
c        k=i
c        do j=1,i-1
c          p=a(j)
c          if(p.gt.m)then
c            k=j
c            m=p
c          endif
c        enddo
c        j=a(i)
c        a(i)=m
c        a(k)=j
c      enddo
      return
      end
