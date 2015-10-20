c *********************************************************************
c * WRITEMESH: escreve a malha particionada no formato do vtk         *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el     - conectividade com material                               *
c *  x     - coordenadas                                              *
c * nnode  - numero de nos                                            *
c * numel  - numero de elementos                                      *
c * nen    - numero de nos por elementos                              *
c * ndm    - numero de dimensoes                                      *
c * filein - prefix do arquivo de saida                               *
c * bvtk   - true BINARY vtk false ASCII vtk                          *
c * legacy - true (formato padrão .vtk) false (formato xlm .vtu)      *
c * nout   - arquivo de saida                                         *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine writeMeshGeo(el     ,x  ,nnode ,numel
     .                       ,nen    ,ndm,filein,bvtk
     .                       ,legacy ,nout)
c ===
      use Malloc 
      implicit none
c ... variaveis da malha      
      integer nnode,numel,nen,ndm
      integer el(nen+1,numel)
      real*8  x(ndm,nnode)
      integer nel,nno
c ... locais     
      integer*8 i_p
      data i_p/1/
      character*15 aux1
      character*30 aux
c ... variaveis dums
      real*8 ddum
      real*4 fdum
c ... arquivo      
      integer nout
      character*80 fileout,name,filein
      logical bvtk,legacy
      integer cod,cod2,gdl
c =====================================================================
c
c ===
      if(legacy) then
        fileout = name(filein,0,0)
      else  
        fileout = name(filein,0,1)
      endif  
      if(bvtk)then
        open(unit=nout,file=fileout,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nout,file=fileout)
      endif  
c     print*,fileout,nout
c =====================================================================
c
c === cabecalho
      if(legacy) then
        write(aux,'(30a)')"Malha part metis" 
        call head_vtk(aux,bvtk,nout) 
      else  
        call head_vtu(nnode,numel,bvtk,nout) 
      endif  
c =====================================================================
c
c === Coordenadas
      if(legacy) then
        call coor_vtk(x,nnode,ndm,bvtk,nout)
      else  
        call coor_vtu(x,nnode,ndm,bvtk,nout)
      endif
c =====================================================================
c
c === Elementos
      if(legacy) then
        call elm_vtk(el,numel,nen,bvtk,nout)
      else  
        call elm_vtu(el,numel,nen,bvtk,nout)
      endif  
c =====================================================================
c
c === cell 
      if(legacy) then
        call cell_data_vtk(numel,bvtk,nout)
      else  
        call cell_data_vtu(bvtk,nout)
      endif  
c ... materiais      
      i_p = alloc_4('p       ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = el(nen+1,nel)
      enddo
      write(aux1,'(15a)')"mat" 
c ... cod = 1 variaveis interias
      cod = 1
      gdl = 1
      if(legacy) then
        call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif 
      i_p = dealloc('p       ')
c .....................................................................
c
c ... elIdG
      i_p = alloc_4('p       ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = nel
      enddo
      write(aux1,'(15a)')"elIdG" 
c ... cod = 1 variaveis interias
      cod = 1
      gdl = 1
      if(legacy) then
        call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif 
      i_p = dealloc('p       ')
c .....................................................................
c =====================================================================
c
c === nos  
c ...       
      if(legacy) then
        call point_data_vtk(nnode,bvtk,nout)
      else  
        call point_data_vtu(bvtk,nout)
      endif  
      i_p = alloc_4('p       ', 1,nnode)
      do nno = 1, nnode
        ia(i_p+nno-1) = nno    
      enddo
      write(aux1,'(15a)')"noIdG"
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 int(4bytes) 
      gdl =  1
      cod =  1
      cod2 = 1
      if(legacy) then
        call pont_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call pont_prop_vtu(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      i_p = dealloc('p       ')

      if(legacy .eqv. .false.) then
        call point_data_finalize_vtu(bvtk,nout)
        call finalize_vtu(bvtk,nout)
      endif
c .....................................................................
c =====================================================================
      close(nout)
      return
      end
c =====================================================================
c *********************************************************************
c 


c *********************************************************************
c * WRITEMESH: escreve a malha particionada no formato do vtk         *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * np     - nos particionados                                        *
c * ep     - elementos particionado                                   *
c * el     - conectividade com material                               *
c *  x     - coordenadas                                              *
c * nnode  - numero de nos                                            *
c * numel  - numero de elementos                                      *
c * nen    - numero de nos por elementos                              *
c * ndm    - numero de dimensoes                                      *
c * filein - prefix do arquivo de saida                               *
c * bvtk   - true BINARY vtk false ASCII vtk                          *
c * legacy - true (formato padrão .vtk) false (formato xlm .vtu)      *
c * nout   - arquivo de saida                                         *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine writemesh(np,ep,el,x,nnode,numel,nen,ndm,filein,bvtk
     .                    ,legacy,nout)
c ===
      use Malloc 
      implicit none
c ... variaveis da malha      
      integer nnode,numel,nen,ndm
      integer np(nnode),ep(numel)
      integer el(nen+1,numel)
      real*8  x(ndm,nnode)
      integer nel,nno
c ... locais     
      integer*8 i_p
      data i_p/1/
      character*15 aux1
      character*30 aux
c ... variaveis dums
      real*8 ddum
      real*4 fdum
c ... arquivo      
      integer nout
      character*80 fileout,name,filein
      logical bvtk,legacy
      integer cod,cod2,gdl
c =====================================================================
c
c ===
      if(legacy) then
        fileout = name(filein,0,106)
      else  
        fileout = name(filein,0,116)
      endif  
      if(bvtk)then
        open(unit=nout,file=fileout,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nout,file=fileout)
      endif  
c     print*,fileout,nout
c =====================================================================
c
c === cabecalho
      if(legacy) then
        write(aux,'(30a)')"Malha part metis" 
        call head_vtk(aux,bvtk,nout) 
      else  
        call head_vtu(nnode,numel,bvtk,nout) 
      endif  
c =====================================================================
c
c === Coordenadas
      if(legacy) then
        call coor_vtk(x,nnode,ndm,bvtk,nout)
      else  
        call coor_vtu(x,nnode,ndm,bvtk,nout)
      endif
c =====================================================================
c
c === Elementos
      if(legacy) then
        call elm_vtk(el,numel,nen,bvtk,nout)
      else  
        call elm_vtu(el,numel,nen,bvtk,nout)
      endif  
c =====================================================================
c
c === cell 
      if(legacy) then
        call cell_data_vtk(numel,bvtk,nout)
      else  
        call cell_data_vtu(bvtk,nout)
      endif  
c ... materiais      
      i_p = alloc_4("p     ", 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = el(nen+1,nel)
      enddo
      write(aux1,'(15a)')"mat" 
c ... cod = 1 variaveis interias
      cod = 1
      gdl = 1
      if(legacy) then
        call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif 
      i_p = dealloc("p     ")
c .....................................................................
c
c ... particionamento
      write(aux1,'(15a)')"part_metis" 
      cod = 1
      gdl = 1
      if(legacy) then
        call cell_prop_vtk(ep,fdum,ddum,numel,aux1,cod,gdl,bvtk,nout)
      else  
        call cell_prop_vtu(ep,fdum,ddum,numel,aux1,cod,gdl,bvtk,nout)
        call cell_data_finalize_vtu(bvtk,nout)
      endif  
c =====================================================================
c
c === nos  
c ...       
      if(legacy) then
        call point_data_vtk(nnode,bvtk,nout)
      else  
        call point_data_vtu(bvtk,nout)
      endif  
      i_p = alloc_4("p     ", 1,nnode)
      do nno = 1, nnode
        ia(i_p+nno-1) = np(nno)
      enddo
      write(aux1,'(15a)')"part_por_no"
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 int(4bytes) 
      gdl =  1
      cod =  1
      cod2 = 1
      if(legacy) then
        call pont_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call pont_prop_vtu(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      i_p = dealloc("p     ")

      if(legacy .eqv. .false.) then
        call point_data_finalize_vtu(bvtk,nout)
        call finalize_vtu(bvtk,nout)
      endif
c .....................................................................
c =====================================================================
      close(nout)
      return
      end
c =====================================================================
c *********************************************************************
c 
c *********************************************************************
c * MAKE_FACE: gera a conectividades dasfaces                         *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el         - conectividade                                        *
c * ie         - tipo do elemento                                     *
c * eload      - cargas nos elementos                                 *
c * face       - indefinido                                           *
c * carga      - indefinido                                           *
c * tipoface   - indefinido                                           *
c * numel      - numero de elementos                                  *
c * maxno      - numero maximo de no por elemento                     *
c * maxface    - numero maximo de faces                               *
c * lineface   - indefinido                                           *
c * triaface   - indefinido                                           *
c * quadface   - indefinido                                           *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * face       - face onde as carga sao aplicadas                     *
c * carga      - carga na face                                        *
c * tipoface   - tipo da face ( 1 linha; 2 tria; 3 quad)              *
c * lineface   - numero de linha com carga                            *
c * triaface   - numero de face triangulares com carga                *
c * quadface   - numero de faca quandrigulares com carga              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine make_face(el      ,ie      ,eload   ,face    ,carga  
     .                    ,tipoface,numel   ,maxno   ,maxface ,maxnoface
     .                    ,lineface,triaface,quadface)
      implicit none
      integer numel,maxno
      integer i,j,k,c,maxface,nface,maxnoface,ty
      integer lineface,triaface,quadface
      integer carga(*),tipoface(*)
      integer eload(maxface+1,*),face(maxnoface,*),el(maxno+1,*),ie(*)
      integer tetraFace(3,4)
      data tetraface/1,3,2  
     .              ,1,2,4
     .              ,2,3,4
     .              ,3,1,4/
      lineface = 0 
      triaface = 0 
      quadface = 0 
      nface    = 0
c ...
      do i = 1, numel
        ty = ie(el(maxno+1,i))
c ... triangulo
        if( ty .eq. 3 ) then
c ....................................................................
c
c ... quadrilatero
        else if(ty .eq. 4 ) then
c ....................................................................
c
c ... Tetraedro
        else if( ty .eq. 6 ) then
c ... verifica se ha carga nas faces do elemento
          do j = 1, 4
            c = eload(j,i)
            if( c .ne. 0) then
              triaface          = triaface + 1
              nface             = nface + 1
              do k = 1, 3
                face(k,nface) = el(tetraFace(k,j),i)
              enddo
              carga(nface)      = c
              tipoface(nface)   = 2
            endif
          enddo
c ....................................................................
c 
c ...  hexaedros      
        else if(ty .eq. 7) then
        endif
c ..................................................................... 
      enddo
c ..................................................................... 
      return
      end
c *********************************************************************
