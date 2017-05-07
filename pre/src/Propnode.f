c *********************************************************************
c * PRONODE  : Gera  as conectividades por no                         *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ix          - conecticidades dos elementos                        *
c * i_nincid    -  ponteiro                                           *
c * i_incid     -  ponteiro                                           *
c * nnode       - numero de nos                                       *
c * numerl      - numero de elementos                                 *
c * nen         - numero de nos por elemento                          *
c * maxgrade    -                                                     *
c * flag        - .true. conectividade com o material                 *
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * i_nincid    -  numero de aresta por no                            *
c * i_incid     -  conectividade por no                               *
c * maxgrade    -  numero de aresta maximo por no                     *
c * ----------------------------------------------------------------- *  
c *********************************************************************
      subroutine propnode(ix,i_nincid,i_incid,nnode,numel,nen,maxgrade
     .                   ,flag)
c === Variaveis
      Use Malloc
      implicit none
c ... malha
      integer ix(*)  
      integer nnode,numel,nen 
c ... particionamento
      integer maxgrade
      integer*8 i_nincid,i_incid
      logical flag
c =====================================================================
c
c ===
      i_nincid = alloc_4('nincid  ',      1,nnode)
      call mynodegrade(ix,ix,nnode,numel,nen,nen
     .                ,ia(i_nincid),maxgrade,flag)
c =====================================================================
      
c ===
      i_incid =  alloc_4('incid   ',maxgrade,nnode)
      call myelmincid(ix,ix,ia(i_incid),ia(i_nincid),nnode,numel
     .               ,nen,nen,maxgrade,flag)
c =====================================================================
c
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * NODEGRADE: Funcao do arquivo colormesh do George. Determina       * 
c * o numero maximo de aparicoes dos nos na malha.                    *  
c * ----------------------------------------------------------------- *  
c * Parametros de entrada :                                           *  
c * ----------------------------------------------------------------- *  
c * ix1    - conectivida da malha(com material)                       *  
c * ix2    - conectivida da malha(sem material)                       *  
c * nnode  - numero de nos                                            *  
c * numel  - numero de elementos                                      *  
c * nenv   - numero de nos de vertices por elementos                  *  
c * nen    - numero de nos por elementos                              *  
c * flag   -.true.com o material                                      *  
c * ----------------------------------------------------------------- *  
c * Parametros de saida :                                             *  
c * ----------------------------------------------------------------- *  
c * nincid   - numero de incidencias por nos                          *  
c * maxgrade - numero de incidencias maximo na malha                  *  
c *********************************************************************
      subroutine mynodegrade(ix1,ix2,nnode,numel,nenv,nen,nincid
     .                      ,maxgrade,flag)
c === Variaveis
      implicit none
c ... Parametros da malha
      integer nnode,numel,nen,nenv
      integer ix1(nen+1,*),ix2(nen,*)
      integer nincid(*)
      integer maxgrade
      logical flag
c ... variaveis auxiliares
      integer i,j,node,grade 
c =====================================================================
c
c ===
      maxgrade = 0
      grade    = 0
c ... zera o vetor     
      do i = 1, nnode
        nincid(i) = 0
      enddo
c .....................................................................      
c
c ...
      if(flag) then
        do i = 1, numel
          do j = 1, nenv
            node = ix1(j,i)
            if( node .gt. 0 ) then
              nincid(node) = nincid(node) + 1
              grade        = nincid(node)
              if ( grade .gt . maxgrade ) maxgrade = grade
            endif
          enddo
        enddo
      else
        do i = 1, numel
          do j = 1, nenv
            node = ix2(j,i)
            if( node .gt. 0 ) then
              nincid(node) = nincid(node) + 1
              grade        = nincid(node)
              if ( grade .gt . maxgrade ) maxgrade = grade
            endif
          enddo
        enddo
      endif  
c .....................................................................
c      print*,"Maxgrade",maxgrade
c      print*,"nincid"
c      do i = 1, nnode
c        print *,i, nincid(i)
c      enddo  
c =====================================================================
      return
      end
c *********************************************************************      
c
c *********************************************************************
c * ELMINID  : Funcao do arquivo Colormesh do George. Determina o     * 
c * numero maximo de aparicoes dos nos na malha.                      *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ix1    - conectivida da malha(sem material)                       *
c * ix2    - conectivida da malha(sem material)                       *
c * nnode  - numero de nos                                            *
c * numel  - numero de elementos                                      *
c * nenv   - numero de nos de vertices por elementos                  *  
c * men    - numero de nos por elementos                              *
c * nincid   - numero de incidencias por nos                          *
c * maxgrade - numero de incidencias maximo na malha                  *
c * flag   -.true.com o material                                      *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * incid    - incidencia dos lementos por no                         *
c *********************************************************************
      subroutine myelmincid(ix1,ix2,incid,nincid,nnode,numel,nenv,nen
     .                     ,maxgrade,flag)
c === Variaveis
      implicit none
c ... Parametros da malha
      integer nnode,numel,nen,nenv
      integer ix1(nen+1,*),ix2(nen,*)
      integer maxgrade
      integer incid(maxgrade,*),nincid(*)
c ... variaveis auxiliares
      integer i,j,node,ipos
      logical flag
c =====================================================================
c
c ===
c ... zera o vetor     
      do i = 1, nnode
        nincid(i) = 0
      enddo
      do i = 1, nnode
        do j = 1 ,maxgrade
          incid(j,i) = 0
        enddo
      enddo  
c .....................................................................      
c
c ...
      if(flag) then
        do i = 1, numel
          do j = 1, nenv
            node = ix1(j,i)
            if( node .gt. 0 ) then
              nincid(node) = nincid(node) + 1
              ipos = nincid(node)
              incid(ipos,node) = i
            endif
          enddo
        enddo
      else  
        do i = 1, numel
          do j = 1, nenv
            node = ix2(j,i)
            if( node .gt. 0 ) then
              nincid(node) = nincid(node) + 1
              ipos = nincid(node)
              incid(ipos,node) = i
            endif
          enddo
        enddo
      endif  
c =====================================================================
c
c === 
c      do i = 1,nnode
c         print*,i,(incid(j,i),j = 1, maxgrade)
c      enddo   
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *****************************************************************
