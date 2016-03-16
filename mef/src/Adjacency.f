      subroutine adjtria3(ix,nodcon,nelcon,nnode,numel,nen,nedge)
c **********************************************************************
c *                                                                    *
c *   ADJTRIA3                                                         *
c *   --------                                                         *
c *                                                                    *
c *   Determina as adjacencias dos elementos (triangulos)              *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   ix(nen+1,numel) - conetividades nodais dos elementos             *
c *   nodcon(nnode)   - nao definido (usado como arranjo auxiliar)     *
c *   nelcon(3,numel) - nao definido                                   *
c *   nnode           - numero de nos                                  *
c *   numel           - numero de elementos                            *
c *   nedge           - nao definido                                   *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   nelcon - elementos adjacentes ao elemento j sao dados por        *
c *            nelcon(i,j), i = 1,2,3.                                 *
c *            nelcon(i,j) = -1 (face pertence ao contorno)            *
c *   nedge  - numero de arestas                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,numel,nen,nedge,nel,is,no,imiss,no1,no2,nel2,is2
      integer ix(nen+1,*),nodcon(*),nelcon(3,*),isnod(2,3),no21,no22
      data isnod /1,2,2,3,3,1/
c.......................................................................
      nedge = 0
      do 100 nel = 1, numel
      do 100 is = 1, 3
         nelcon(is,nel) = 0
  100 continue
      do 200 no = 1, nnode
         nodcon(no) = 0
  200 continue
c ...........................
  300 continue
      imiss = 0
      do 400 nel = 1, numel
      do 400 is = 1, 3
         if (nelcon(is,nel) .ne. 0) go to 400
         no1 = ix(isnod(1,is),nel)
         no2 = ix(isnod(2,is),nel)
         if (nodcon(no1) .eq. 0 .or. nodcon(no2) .eq. 0) then
            nodcon(no1) = nel
            nodcon(no2) = nel
            imiss = 1
         endif
  400 continue
      do 550 nel = 1, numel
      do 550 is = 1, 3
         if (nelcon(is,nel) .ne. 0) go to 550
         no1  = ix(isnod(1,is),nel)
         no2  = ix(isnod(2,is),nel)
         nel2 = nodcon(no1)
         if(nel2 .gt. 0) then
         if(nel2.eq.nodcon(no2).and. nel2.ne.nel) then
            do 510 is2 = 1, 3
               no21 = ix(isnod(1,is2),nel2)
               no22 = ix(isnod(2,is2),nel2)
               if (no21 .eq. no2 .and. no22 .eq. no1) then
                  nelcon(is,nel)  = nel2
                  nelcon(is2,nel2)= nel
                  nodcon(no1) = 0
                  nodcon(no2) = 0
                  imiss = 1
                  nedge = nedge + 1
               endif
  510       continue
         endif
         endif
  550 continue
      do 600 nel = 1, numel
      do 600 is = 1, 3
         if (nelcon(is,nel) .ne. 0) go to 600
         no1 = ix(isnod(1,is),nel)
         no2 = ix(isnod(2,is),nel)
         if (nodcon(no1).eq.nodcon(no2) .and. nodcon(no1).eq.nel) then
            nelcon(is,nel)   = -1
            nodcon(no1) = 0
            nodcon(no2) = 0
            imiss = 1
            nedge = nedge + 1
         endif
  600 continue
      if (imiss .eq. 1) go to 300
c ......................................................................      
      return
      end
      subroutine adjquad4(ix,nodcon,nelcon,nnode,numel,nen,nedge)
c **********************************************************************
c *                                                                    *
c *   ADJQUAD4                                                         *
c *   --------                                                         *
c *                                                                    *
c *   Determina as adjacencias dos elementos (quadrilateros)           *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   ix(nen+1,numel) - conetividades nodais dos elementos             *
c *   nodcon(nnode)   - nao definido (usado como arranjo auxiliar)     *
c *   nelcon(3,numel) - nao definido                                   *
c *   nnode           - numero de nos                                  *
c *   numel           - numero de elementos                            *
c *   nen             - numero de nos por elemento                     *
c *   nedge           - nao definido                                   *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   nelcon - elementos adjacentes ao elemento j sao dados por        *
c *            nelcon(i,j), i = 1,2,3,4.                               *
c *   nedge  - numero de arestas                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,numel,nen,nedge,nel,is,no,imiss,no1,no2,nel2,is2
      integer ix(nen+1,*),nodcon(*),nelcon(4,*),isnod(2,4),no21,no22
      data isnod /1,2,2,3,3,4,4,1/
c.......................................................................
      nedge = 0
      do 100 nel = 1, numel
      do 100 is = 1, 4
         nelcon(is,nel) = 0
  100 continue
      do 200 no = 1, nnode
         nodcon(no) = 0
  200 continue
c ...........................
  300 continue
      imiss = 0
      do 400 nel = 1, numel
      do 400 is = 1, 4
         if (nelcon(is,nel) .ne. 0) go to 400
         no1 = ix(isnod(1,is),nel)
         no2 = ix(isnod(2,is),nel)
         if (nodcon(no1) .eq. 0 .or. nodcon(no2) .eq. 0) then
            nodcon(no1) = nel
            nodcon(no2) = nel
            imiss = 1
         endif
  400 continue
      do 550 nel = 1, numel
      do 550 is = 1, 4
         if (nelcon(is,nel) .ne. 0) go to 550
         no1  = ix(isnod(1,is),nel)
         no2  = ix(isnod(2,is),nel)
         nel2 = nodcon(no1)
         if(nel2 .gt. 0) then
         if(nel2.eq.nodcon(no2).and. nel2.ne.nel) then
            do 510 is2 = 1, 4
               no21 = ix(isnod(1,is2),nel2)
               no22 = ix(isnod(2,is2),nel2)
               if (no21 .eq. no2 .and. no22 .eq. no1) then
                  nelcon(is,nel)  = nel2
                  nelcon(is2,nel2)= nel
                  nodcon(no1) = 0
                  nodcon(no2) = 0
                  imiss = 1
                  nedge = nedge + 1
               endif
  510       continue
         endif
         endif
  550 continue
      do 600 nel = 1, numel
      do 600 is = 1, 4
         if (nelcon(is,nel) .ne. 0) go to 600
         no1 = ix(isnod(1,is),nel)
         no2 = ix(isnod(2,is),nel)
         if (nodcon(no1).eq.nodcon(no2) .and. nodcon(no1).eq.nel) then
            nelcon(is,nel)   = -1
            nodcon(no1) = 0
            nodcon(no2) = 0
            imiss = 1
            nedge = nedge + 1
         endif
  600 continue
      if (imiss .eq. 1) go to 300
      nedge = nedge + 2*numel
c ......................................................................      
      return
      end
      subroutine adjhexa8(ix,nodcon,nelcon,nnode,numel,noMax)
c **********************************************************************
c *                                                                    *
c *   HEXA8FACE - determina as adjacencias do hexaedro                 *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   ix(9,numel)- conetividades nodais dos elementos                  *
c *   nodcon(nnode)   - nao definido (usado como arranjo auxiliar)     *
c *   nelcon(n,numel) - nao definido                                   *
c *   nnode           - numero de nos                                  *
c *   numel           - numero de elementos                            *
c *   nVert           - numero de vertives por elemento                *
c *   noMax           - numero maximo de nos por elemento              *
c *   nedge           - nao definido                                   *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   nelcon - elementos adjacentes                                    *
c *   need   - numero de arestas por elemento                          *
c *   nedge  - numero de arestas                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      integer numel,nnode,ipass,noMax
      integer ix(noMax+1,numel),nodcon(nnode),nelcon(6,numel)
      integer node(4),i,j,k,el,miss,hexa8face
c ============================================================================
      do 55 i = 1, numel
      do 50 j = 1, 6
         nelcon(j,i) = -2
   50 continue
   55 continue
      do 60 i = 1, nnode
         nodcon(i) = -2
   60 continue  
c .........................................................................
      ipass = 0
  100 continue
      miss = 0
      do 210 i = 1, numel
      do 200 j = 1, 6
         if(nelcon(j,i) .eq. -2) then
            call hexa8fnod(i,j,ix,node,noMax)
            if ((nodcon(node(1)).eq.-2).and.(nodcon(node(2)).eq.-2).and.
     .          (nodcon(node(3)).eq.-2).and.(nodcon(node(4)).eq.-2))then
                 nodcon(node(1)) = i
                 nodcon(node(2)) = i
                 nodcon(node(3)) = i
                 nodcon(node(4)) = i
                 miss = 1
                 goto 210
            endif
         endif
  200 continue
  210 continue
c ......................................................................
      do 320 i = 1, numel
      do 310 j = 1, 6
         if(nelcon(j,i).eq.-2) then
            call hexa8fnod(i,j,ix,node,noMax)
            el = nodcon(node(1))
            if(el .ne. -2) then
             if((nodcon(node(2)).eq.el).and.(nodcon(node(3)).eq.el).and.
     .          (nodcon(node(4)).eq.el).and.(el.ne.i)) then
                 k = hexa8face(el,ix,node,noMax)
                 if (k.eq.-1) then
                    print*,'*** <ADJHEXA8> Erro na vizinhanca ***',j
                    stop
                 endif
                 nelcon(k,el) = i
                 nelcon(j,i)  = el
                 nodcon(node(1)) = -2
                 nodcon(node(2)) = -2
                 nodcon(node(3)) = -2
                 nodcon(node(4)) = -2
                 miss = 1             
             endif
            endif
         endif
  310 continue
  320 continue
c ......................................................................
      do 410 i = 1, numel
      do 400 j = 1, 6
         if(nelcon(j,i).eq.-2) then
            call hexa8fnod(i,j,ix,node,noMax)
            el = nodcon(node(1))
            if((nodcon(node(2)).eq.el).and.(nodcon(node(3)).eq.el).and.
     .         (nodcon(node(4)).eq.el).and.(el.eq.i)) then
                 nelcon(j,i) = -1
                 nodcon(node(1)) = -2
                 nodcon(node(2)) = -2
                 nodcon(node(3)) = -2
                 nodcon(node(4)) = -2
                 miss = 1        
            endif
         endif
  400 continue
  410 continue
      ipass = ipass + 1
      if (miss .eq. 1) goto 100
c ......................................................................
      return
      end
      subroutine Hexa8fnod(k,j,ix,node,noMax)
c **********************************************************************
c *                                                                    *
c *   HEXA8FNOD - determina os nos da face j do elemento k             *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   k      - numero do elemento                                      *
c *   j      - numero da face do elemento                              *
c *   ix(9,numel)- conetividades nodais dos elementos                  *
c *   node(4)- nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   node - nos da face j (numeracao local)                           *
c *                                                                    *
c **********************************************************************
      implicit none
      integer noMax
      integer k,j,ix(noMax+1,*),node(*),fnode(4,6),i
      data fnode /1,2,3,4,5,8,7,6,1,5,6,2,4,3,7,8,1,4,8,5,2,6,7,3/
c ======================================================================
      do 100 i = 1, 4
         node(i) = ix(fnode(i,j),k)
  100 continue
      return
      end
      integer function hexa8face(k,ix,node,noMax)
c **********************************************************************
c *                                                                    *
c *   HEXA8FACE - determina a face do hexaedro k adjacente a face j    *
c *   ---------               cujos nos estao armazenados em node(4)   *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   k          - numero do elemento adjacente                        *
c *   node - nos da face j (numeracao local)                           *
c *   ix(9,numel - conetividades nodais dos elementos                  *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer noMax
      integer k,ix(noMax+1,*),node(*),ind(3,4),no(4),i,j
      data ind /2,3,4,3,4,1,4,1,2,1,2,3/
c ======================================================================
      hexa8face = -1
      do 200 i = 1, 6
         call hexa8fnod(k,i,ix,no,noMax)
         do 100 j = 1, 4
            if(no(1) .eq. node(j)) then
               if((no(2) .eq. node(ind(3,j))).and.
     .            (no(3) .eq. node(ind(2,j))).and.
     .            (no(4) .eq. node(ind(1,j)))) then
                   hexa8face = i
                   return
               endif
            endif
  100    continue
  200 continue
c ......................................................................  
      return
      end
      integer function edg(l,k,e,n)
c **********************************************************************
c *                                                                    *
c *   EDG - determina qual lado do triangulo l e adjacente             *
c *   ---   ao triangulo k.                                            *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   l          - numero do elemento adjacente                        *
c *   k          - numero do elemento                                  *
c *   e(n,numel) - adjacencias dos elementos                           *
c *   n          - numero de arestas do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   edg - lado do elemento l adjacente ao elemento k                 *
c *                                                                    *
c **********************************************************************
      implicit none
      integer l,k,n,e(n,*),i
c ......................................................................
      do 10 i = 1, n
         if (e(i,l) .eq. k) then
            edg = i
            return
         endif
   10 continue
c ......................................................................   
      print*, '*** Erro na funcao EDG: elementos nao adjacentes ***'
      print*, l,e(1,l),e(2,l),e(3,l)
      print*, k,e(1,k),e(2,k),e(3,k)
      stop
      end
      subroutine adjtetra4(numel,nnode,nen,nodcon,neighbors,ix)
c **********************************************************************
c *                                                                    *
c *   TETRA4NEIGHBORS                                     07/08/2014   *
c *                                                                    *
c *   Rotina para identificar os vizinhos dos elementos                *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   numel - numero de elementos                                      *
c *   nnode - numero de nos                                            *
c *   nen - numero de nos por elemento                                 *      
c *   nodcon - variavel auxiliar do tamanho de nnode                   *
c *   neighbors - armazena os vizinhos de cada elemento                *
c *   ix(nen+1,*) - conetividades nodais dos elementos                 *
c *   node(3) - nos da face k do elemento nel                          *      
c *                                                                    *
c *   Parametros de saida:                                             *
c *   n - numero da face do elemento vizinho                           * 
c *       ou -1 caso nao tenha vizinho                                 *
c *                                                                    *
c *                                                                    *  
c **********************************************************************
      implicit none
      integer numel,nnode,nen,nodcon(*),neighbors(4,*)
      integer ix(nen+1,*),node(3)
      integer i,j,miss,el,k,tetra4face
c ----------------------------------------------------------------------
      do 200 i = 1, numel
      do 100 j = 1, 4
          neighbors(j,i) = -2
  100 continue
  200 continue
      do 300 i = 1, nnode
         nodcon(i) = -2   
  300 continue
c ----------------------------------------------------------------------
      miss = 1
      do while (miss .eq. 1)
          miss = 0
          do 500 i = 1, numel
          do 400 j = 1, 4
             if (neighbors(j,i) .eq. -2) then
                call tetra4fnod(i,j,ix,node,nen)
                if (nodcon(node(1)).eq.-2 .and. nodcon(node(2)).eq.-2
     .                              .and. nodcon(node(3)).eq.-2) then
                   nodcon(node(1)) = i
                   nodcon(node(2)) = i
                   nodcon(node(3)) = i 
                   miss = 1
                   goto 500
                endif
             endif
  400     continue
  500     continue
c ----------------------------------------------------------------------          
          do 700 i = 1, numel
          do 600 j = 1, 4
             if (neighbors(j,i) .eq. -2) then
                call tetra4fnod(i,j,ix,node,nen)
                el = nodcon(node(1))
                if(el .ne. -2) then
                   if(nodcon(node(2)) .eq. el .and. nodcon(node(3))
     .                                .eq. el .and. el .ne. i) then
                      k = tetra4face(el,ix,node,nen)
                      if(k .eq. -1) then
                          print*, 'Erro !'
                          stop
                      endif
                      neighbors(k,el) = i
                      neighbors(j,i)  = el
                      nodcon(node(1)) = -2
                      nodcon(node(2)) = -2
                      nodcon(node(3)) = -2
                      miss = 1
                   endif
                endif
             endif
  600     continue
  700     continue
c ---------------------------------------------------------------------- 
          do 900 i = 1, numel
          do 800 j = 1, 4
             if (neighbors(j,i) .eq. -2) then
                call tetra4fnod(i,j,ix,node,nen)
                el = nodcon(node(1))
                if (nodcon(node(2)) .eq. el .and. nodcon(node(3)) 
     .                              .eq. el .and. el .eq. i) then
                   neighbors(j,i) =  -1
                   nodcon(node(1)) = -2
                   nodcon(node(2)) = -2
                   nodcon(node(3)) = -2
                   miss = 1
                endif
             endif
  800     continue
  900     continue
c ----------------------------------------------------------------------
      enddo
c ----------------------------------------------------------------------
      do 1100 i = 1, numel
      do 1000 j = 1, 4
         el = neighbors(j,i)
         if (el .eq. -1) neighbors(j,i) = 0
1000  continue
1100  continue 
c ----------------------------------------------------------------------      
      return
      end
c **********************************************************************
      subroutine mk_elconn_hex_quad(ix    ,nelcon ,numel
     .                             ,nnode ,nnodev ,nen
     .                             ,nenv  ,nMaxViz)
c **********************************************************************
c *                                                                    *
c *   MK_ELCON_HEX_QUAD - gera a connectividade dos elementos          *
c *   -----------------   quadraticos a partir dos verticeis dos       *
c *                       elementos lineares                           *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   ix(*,numel) - conetividades nodais dos elementos(Vertices apenas)*
c *   nelcon(j,i) - elementos vizinhos ao elemento ao i                *
c *   numel       - numero de elementos                                *
c *   nnodev      - numero de nos dos vertices                         *
c *   nnode       - numero total de nos                                *
c *   nen         - numero maximo de nos por elemento                  *
c *   nenv        - numero maximo de vertices por elemento             *
c *   maxViz      - numero maximo de vizinhos                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *   ix(*,numel) - conetividades nodais dos elementos atualizadas     *
c **********************************************************************
      implicit none
      integer i,j,k,nElViz,numFace
      integer nenv,nen,nMaxViz,numel,nnodev,nnode,nno,noFace,noF1,noF2
      integer ix(nen+1,*),nelcon(nMaxViz,*)
      integer node1(4),node2(4)
      integer hexa20(4,6),hexa20r(4,6)
      integer hexa8face
c ... hexaedro
      noFace = 4
c .....................................................................
      nno = nnodev + 1
      call Hexa20faceNodi(hexa20,.false.)
      call Hexa20faceNodi(hexa20r,.true.)
c ... loop nos elementos
      do 100 i = 1, numel
        do 110 j = nenv, nen
          if( ix(j,i) .eq. 0 )  then
            ix(j,i) = nno
            nno = nno + 1
          endif             
 110    continue
c ... loop nas faces
c        do j = 1, nMaxViz
c          do k = 1, noFace
c            noF1 = hexa20(4+k,j)
c            if( ix(noF1,i) .eq. 0 )  then
c              ix(noF1,i) = nno
c              nno = nno + 1
c            endif         
c          enddo
c        enddo
c ... vizinhos
        do 120 j = 1, nMaxViz
          numFace = 0
          call Hexa8fnod(i,j,ix,node1,nen)
c ... elemento vizinho a face
          nElViz  = nelcon(j,i)
          if(nElViz .gt. 0) then
            numFace = hexa8face(nElViz,ix,node1,nen)
            call Hexa8fnod(nElViz,numFace ,ix,node2,nen)
            do 130 k = 1, noFace
              noF2 = hexa20r(k,numFace)
              if( ix(noF2,nElViz) .eq. 0 ) then
                noF1            = hexa20(k,j)
                ix(noF2,nElViz) = ix(noF1,i)
              endif
  130       continue
          endif
  120   continue
  100 continue
c .....................................................................
c
c ...
      nnode = nno - 1
c .....................................................................
c
c ...
c     do i = 1, numel
c       print*,i,ix(1:nen+1,i)
c     enddo
c ..................................................................... 
      return
      end
c **********************************************************************
      integer function tetra4face(nel,ix,node,nen)
c **********************************************************************
c *                                                                    *
c *   TETRA4FACE                                          07/08/2014   *
c *                                                                    *
c *   Rotina para comparar as faces do tetraedro                       *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   nel - numero do elemento                                         *
c *   nen - numero de nos por elemento                                 *      
c *   ix(nen+1,*) - conetividades nodais dos elementos                 *
c *   node(3) - nos da face k do elemento nel                          *      
c *                                                                    *
c *   Parametros de saida:                                             *
c *   n - numero da face do elemento vizinho                           * 
c *       ou -1 caso nao tenha vizinho                                 *
c *                                                                    *
c *                                                                    *  
c **********************************************************************
      integer nel,nen,ix(nen+1,*),node(*)
      integer i,j,n,no(3),edge(2,3)
      data edge/2,3,3,1,1,2/
c ----------------------------------------------------------------------
      tetra4face = -1
      do 200 i = 1, 4
         call tetra4fnod(nel,i,ix,no,nen)
         do 100 j = 1, 3
            if(no(1) .eq. node(j)) then
               if(no(2) .eq. node(edge(2,j)) .and. no(3) .eq.
     .                       node(edge(1,j))) then
                  tetra4face = i
                  return
               endif
            endif
 100     continue
 200  continue
      return
      end
      subroutine tetra4fnod(nel,k,ix,node,nen)
c **********************************************************************
c *                                                                    *
c *   tetra4fnod                                     07/08/2014        *
c *                                                                    *
c *   Rotina para identificar os nos das faces do tetraedro            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   nel - numero do elemento                                         *
c *   k   - numero da face                                             *
c *   nen - numero de nos por elemento                                 *      
c *   ix(nen+1,*) - conetividades nodais dos elementos                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   node(3) - nos da face k do elemento nel                          *
c *                                                                    *
c *                                                                    *  
c **********************************************************************
      integer nel,k,nen,ix(nen+1,*),node(*)
      integer fnode(3,4),i,j
c ... normal externa
c     data fnode/4,3,2,4,2,1,4,1,3,1,2,3/
c ... normal interna
      data fnode/2,3,4,1,4,3,1,2,4,1,3,2/
c ----------------------------------------------------------------------
      do 100 j = 1, 3
         i       = fnode(j,k)
         node(j) = ix(i,nel)
  100 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine mk_elconn_tetra_quad_v1(ix    ,incid  ,nincid
     .                                  ,numel ,nnode  ,nnodev 
     .                                  ,nen   ,nenv   ,maxgrade)
c **********************************************************************
c *                                                                    *
c *   MK_ELCON_TETRA_QUAD - gera a connectividade dos elementos        *
c *   -----------------   quadraticos a partir dos verticeis dos       *
c *                       elementos lineares                           *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   ix(*,numel) - conetividades nodais dos elementos(Vertices apenas)*
c *   incid(j,i)  - elementos ligados ao no i                          *
c *   nincid(j,i) - numero de elementos ligados ao no i                *
c *   numel       - numero de elementos                                *
c *   nnode       - numero total de nos                                *
c *   nnodev      - numero de nos dos vertices                         *
c *   nen         - numero maximo de nos por elemento                  *
c *   nenv        - numero maximo de vertices por elemento             *
c *   maxgrade    - numero maximo de elementos que compartilham um no  *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *   ix(*,numel) - conetividades nodais dos elementos atualizadas     *
c **********************************************************************
      implicit none
      integer i,j,k,l,nElViz,numFace
      integer nenv,nen,maxgrade,numel,nnodev,nnode,nno,nel
      integer no1,no2,no3,no1v,no2v,no3v
      integer ix(nen+1,*),incid(maxgrade,*),nincid(*)
      integer iEdge(3,6),nedge
c ... tetraedro
      nedge =  6
      call tetra10edgeNod(iEdge)
c ... 
      nno = nnodev
c ... loop nos elementos
      do 100 i = 1, numel
        do 110 j = nenv, nen
          if( ix(j,i) .eq. 0 )  then
            nno = nno + 1
            ix(j,i) = nno
          endif             
 110    continue
c .....................................................................
c
c ... loop nas arestas
        do 120 j = 1, nedge
c ... no vertices
          no1     = ix(iEdge(1,j),i)
          no2     = ix(iEdge(2,j),i)
c ... no central
          no3     = ix(iEdge(3,j),i)
c ... loop nos elementos que compatilham esse no
          do 130 k = 1, nincid(no1)          
            nel = incid(k,no1)
            if( nel .ne. i) then 
c ... loop nas arestas
              do 140 l = 1, nedge
                no1v     = ix(iEdge(1,l),nel)
                no2v     = ix(iEdge(2,l),nel)
c ...
                if((no1 .eq. no1v) .and. (no2 .eq. no2v) .or.
     .             (no1 .eq. no2v) .and. (no2 .eq. no1v) ) then
c ... no central
                  no3v     = ix(iEdge(3,l),nel)
                  if( no3v .eq. 0 ) then
c ...
                    ix(iEdge(3,l),nel) = no3
                  endif
c .....................................................................
                endif
c .....................................................................
 140          continue
c .....................................................................
            endif
c .....................................................................
 130      continue
c .....................................................................
 120    continue
c .....................................................................
 100  continue
c .....................................................................
c
c ..................................................................... 
c
c ...
      nnode = nno
c .....................................................................
c
c ...
c     do i = 1, numel
c       print*,i,ix(1:nen,i)
c     enddo
      return
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   hexa20faceNodi numeracao dos nos medios das faces                *
c *   -------------                                                    *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   hexa        - nao definido                                       *
c *   rev         - true - numeracao reversa; false numera direta      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *   hexa - numera dos nos medios das faces                           *
c **********************************************************************
      subroutine Hexa20faceNodi(hexa,rev)
      implicit none
      integer hexa(*)
      integer hexa20(24),hexa20Rev(24)
      logical rev
      data hexa20 / 9,10,11,12
     .            , 16,15,14,13
     .            , 17,13,18, 9
     .            , 11,19,15,20
     .            , 12,20,16,17
     .            , 18,14,19,10/
      data hexa20Rev / 12,11,10, 9
     .               , 13,14,15,16
     .               ,  9,18,13,17
     .               , 20,15,19,11
     .               , 17,16,20,12
     .               , 10,19,14,18/
      if(rev) then
        hexa(1:24) = hexa20Rev(1:24) 
      else 
        hexa(1:24) = hexa20(1:24) 
      endif
      return
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   tetra10faceNodi numeracao dos nos medios das faces               *
c *   -------------                                                    *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   hexa        - nao definido                                       *
c *   rev         - true - numeracao reversa; false numera direta      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *   hexa - numera dos nos medios das faces                           *
c **********************************************************************
      subroutine tetra10faceNodi(tetra,rev)
      implicit none
      integer tetra(*)
      integer tetra10(12),tetra10Rev(12)
      logical rev
      data tetra10 / 8, 9, 10  
     .            ,  6, 7, 9   
     .            ,  5,10, 7   
     .            ,  5, 6, 8/  
      data tetra10Rev / 10, 9, 8   
     .                ,  6, 7, 9   
     .                ,  7,10, 5   
     .                ,  8, 6, 5/  
      if(rev) then
        tetra(1:12) = tetra10Rev(1:12) 
      else 
        tetra(1:12) = tetra10(1:12) 
      endif
      return
      end
c *********************************************************************
c
c ********************************************************************* 
      subroutine mkCoorQuad(x       ,xq
     .                     ,el
     .                     ,numel   ,nen             
     .                     ,tnnode  ,nnode  ,ndm)      
      implicit none
      integer maxEdge
      parameter (maxEdge = 12) 
c ...
      integer i,j,k
      integer numel,tnnode,nnode,ndm
      integer nen,nedge,no1,no2,no3
      integer el(nen+1,*)
      integer iEdge(3,maxEdge)
      real*8 x(ndm,*),xq(ndm,*)
c
      nedge = 0  
c ... tetraedros de 10 nos 
      if( nen .eq. 10 ) then
        nedge =  6
        call tetra10edgeNod(iEdge) 
c ... hexaedros de 20 nos 
      else if( nen .eq. 20 ) then
        nedge = 12
        call hexa20edgeNod(iEdge) 
      endif
c ...
      do i = 1, nnode
        do j = 1, ndm
          xq(j,i) = x(j,i)
        enddo
      enddo
c .....................................................................
c
c ...
      do i = 1, numel
        do j = 1, nedge
c ... no vertices
          no1     = el(iEdge(1,j),i)
          no2     = el(iEdge(2,j),i)
c ... no central
          no3     = el(iEdge(3,j),i)
          do k = 1, ndm
            xq(k,no3) = 0.5d0*(x(k,no1) + x(k,no2)) 
          enddo
        enddo
      enddo
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************                
c * HEXA20EDGENOD : numera de nos por aresta do hexaedro de 20 nos    *
c *********************************************************************                
      subroutine Hexa20edgeNod(iEdge)  
      integer iEdge(36),edge(36)
      data edge/  1, 2, 9
     .          , 2, 3,10
     .          , 3, 4,11
     .          , 4, 1,12
     .          , 5, 6,13
     .          , 6, 7,14
     .          , 7, 8,15
     .          , 8, 5,16
     .          , 1, 5,17
     .          , 2, 6,18
     .          , 3, 7,19
     .          , 4, 8,20/
c ...
      iEdge(1:36) = edge(1:36) 
c .....................................................................
      return
      end
c *********************************************************************                
c
c *********************************************************************       
c * TETRA10EDGENOD : numera de nos por aresta do tetraedro de 10 nos  *
c *********************************************************************                
      subroutine tetra10edgeNod(iEdge)  
      integer iEdge(18),edge(18)
      data edge/  1, 2, 5
     .          , 1, 3, 6
     .          , 1, 4, 7
     .          , 2, 3, 8
     .          , 3, 4, 9
     .          , 4, 2,10/
c ...
      iEdge(1:18) = edge(1:18) 
c .....................................................................
      return
      end
c *********************************************************************                
