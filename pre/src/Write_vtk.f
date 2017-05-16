c *********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 04/11/2016                                    * 
c * ------------------------------------------------------------------ *  
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
c ...
      character*8 malloc_name
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
c =====================================================================
c
c === cabecalho
      if(legacy) then
        write(aux,'(30a)')"Malha part metis" 
        call head_vtk(aux,bvtk,0.0,0,.false.,nout) 
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
      malloc_name = 'p'
      i_p = alloc_4(malloc_name, 1,numel)
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
      i_p = dealloc(malloc_name)
c .....................................................................
c
c ... elementos com numeracao global     
      write(aux1,'(15a)')"gElm"
      gdl =  1   
      cod =  1
      cod2 = 1
      malloc_name = 'p'
      i_p = alloc_4(malloc_name, 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = nel
      enddo
      if(legacy) then
        call point_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,ndm,gdl,cod
     .                     ,cod2   ,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,ndm,gdl,cod
     .                     ,cod2   ,bvtk,nout)
      endif
      i_p = dealloc(malloc_name)
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
      malloc_name = 'p'
      i_p = alloc_4(malloc_name, 1,nnode)
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
        call point_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      endif
      i_p = dealloc(malloc_name)
c .....................................................................
c
c ... nos com numeracao global     
      write(aux1,'(15a)')"gNode"
      gdl =  1   
      cod =  1
      cod2 = 1
      malloc_name = 'p'
      i_p = alloc_4(malloc_name,1,nnode)
      do nno = 1, nnode
        ia(i_p+nno-1) = nno
      enddo
      if(legacy) then
        call point_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      endif
      i_p = dealloc(malloc_name)
c .....................................................................
c
c ...
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
c * WRITEMESHPART:escreve a parte da malha do processo no formato vtk *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * np           - nos particionados                                  *
c * ep           - elementos particionado                             *
c *  x           - coordenadas                                        *
c * ie           - tipo do elemento                                   *
c * elLG         - mapa local -> global de elementos                  *
c * noLG         - mapa local -> global de nos                        *
c * nodedist     - distribuicao das propriedades nodais               *
c * elmdist      - distribuicao das propriedades dos elmentos         *
c * nnode        - numero de nos                                      *
c * numel        - numero de elementos                                *
c * my_nnode     - numero de nos locais                               *
c * my_numel     - numero de elementos locais                         *
c * my_numel_nov - numero de elementos locais non-overllaping         *
c * my_numel_ov  - numero de elementos locais overllaping             *
c * maxno    - numero de maximo de nos por elementos                  *
c * ndm      - numero de dimensoes                                    *
c * filein - prefix do arquivo de saida                               *
c * bvtk     - true BINARY vtk false ASCII vtk                        *
c * legacy   - true (formato padrão .vtk) false (formato xlm .vtu)    *
c * nout     - arquivo de saida                                       *
c * npart    - numero da particao a ser escrita                       *
c * div      - numero de particoes totais deste processo              *
c * my_rank  - id do processo (MPI)                                   *
c * npes     - numero de processos (MPI)                              *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine writemeshpart(np,ep,el,x,ie,elLG,noLG,noGL,nodedist
     .                        ,elmdist,nnode,numel,my_nnode,my_numel
     .                        ,my_numel_nov ,my_numel_ov 
     .                        ,maxno,ndm,filein,bvtk,legacy,nout,npart
     .                        ,div,mype,npes,MPIW)
c ===
      use Malloc 
      implicit none
      include "elementos.fi"
      integer nnode,numel,my_nnode,my_numel,my_numel_nov,my_numel_ov
      integer maxno,ndm
      integer np(*),ep(*),nodedist(*),elmdist(*)
      integer elLG(*),noLG(*),noGL(*),ty
      integer el(maxno+1,*),ie(*)
      integer i,lel,tolel
      integer mype,npes,MPIW,npart,div
      integer nout
      integer cod
      real*4 fdum
      real*8 ddum
      integer nt1,nq1,nh1,nt2,nq2,nh2,ntr1,ntr2,np1,np2,ntrq1,ntrq2
      integer nhq1,nhq2
      integer*8 i_lel,i_lno,i_p,i_g
      real*8 x(ndm,*)
      character*30 aux
      character*15 aux1
      character*80 fileout,name,filein
      logical bvtk,legacy,flagTria3,flagQuad4,flagTetra4,flagHexa8
     .        ,flagPrism6,flagTetra10,flagHexa20
c =====================================================================
c
c === casa particular que ha processos  ociosos  na etapa   
c da geracao do mapa     
      if( div .lt. 1 ) then
c ... coor      
        i_g   = alloc_8('g        ',  ndm,nnode)
        call global_double(ia(i_g),x,nodedist,nnode,ndm,npes,mype
     .                    ,MPIW)   
        i_g  = dealloc('g        ')
c ... elementos        
        i_g   = alloc_4('g        ',maxno+1,numel)
        call global_int(ia(i_g),el,elmdist,numel,maxno+1,npes,mype
     .                 ,MPIW)   
        i_g  = dealloc('g        ')
        return
      endif
c ===
      if(legacy) then
        fileout = name(filein,npart,107)
      else
        fileout = name(filein,npart,111)
      endif
      if(bvtk)then
        open(unit=nout,file=fileout,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nout,file=fileout)
      endif
c =====================================================================
c
c === Caracteriscas da funcao vtk
      nt1          = ntria3(1)
      nt2          = ntria3(2)
      nq1          = nquad4(1)
      nq2          = nquad4(2)
      ntr1         = ntetra4(1)
      ntr2         = ntetra4(2)
      nh1          = nhexa8(1)
      nh2          = nhexa8(2)
      np1          = nprism6(1)
      np2          = nprism6(2)
      ntrq1        = ntetra10(1)
      ntrq2        = ntetra10(2)
      nhq1         = nhexa20(1)
      nhq2         = nhexa20(2)
      ntria3(1:4)  = 0
      nquad4(1:4)  = 0
      ntetra4(1:4) = 0
      nhexa8(1:4)  = 0
      nprism6(1:4) = 0
      ntetra10(1:4)= 0
      nhexa20(1:4) = 0
      tolel        = 1
      flagTria3    = .true.
      flagQuad4    = .true.
      flagTetra4   = .true.
      flagHexa8    = .true.
      flagPrism6   = .true.
      flagTetra10  = .true.
      flagHexa20   = .true.
      do i = 1, my_numel_nov
        lel = elLG(i)
        ty = geom_el_type(ie(el(maxno+1,lel)))
c ... tria3
        if( ty .eq. type_tria3) then
          if(flagTria3) ntria3(2) = tolel    
          ntria3(1) = ntria3(1) + 1
          tolel     = tolel + 1
          flagTria3 = .false.
c .....................................................................
c
c ... quad4
        else if( ty .eq. type_quad4) then
          if(flagQuad4) nquad4(2) = tolel
          nquad4(1) = nquad4(1) + 1
          tolel     = tolel + 1
          flagQuad4 = .false.     
c .....................................................................
c
c ... tetra4
        else if( ty .eq. type_tetra4) then
          if(flagTetra4) ntetra4(2) = tolel
          ntetra4(1) = ntetra4(1) + 1
          tolel     = tolel + 1
          flagTetra4 = .false.     
c .....................................................................
c
c ... hexa8
        else if( ty .eq. type_hexa8) then
          if(flagHexa8) nhexa8(2) = tolel
          nhexa8(1) = nhexa8(1) + 1
          tolel     = tolel + 1
          flagHexa8 = .false.     
c .....................................................................
c
c ... prism6
        else if( ty .eq. type_prism6) then
          if(flagprism6) nprism6(2) = tolel
          nprism6(1) = nprism6(1) + 1
          tolel      = tolel + 1
          flagPrism6 = .false. 
c .....................................................................
c
c ... tetra10
        else if( ty .eq. type_tetra10) then
          if(flagtetra10) ntetra10(2) = tolel
          ntetra10(1) = ntetra10(1) + 1
          tolel      = tolel + 1
          flagTetra10= .false.   
c .....................................................................
c
c ... hexa20
        else if( ty .eq. type_hexa20) then
          if(flaghexa20) nhexa20(2) = tolel
          nhexa20(1)  = nhexa20(1) + 1
          tolel       = tolel + 1
          flaghexa20 = .false.       
        endif  
c .....................................................................
      enddo
      flagTria3    = .true.
      flagQuad4    = .true.
      flagTetra4   = .true.
      flagHexa8    = .true.
      flagPrism6   = .true.
      flagTetra10  = .true.
      flagHexa20   = .true.
      do i = my_numel_nov + 1, my_numel_nov+my_numel_ov
        lel = elLG(i)
        ty = geom_el_type(ie(el(maxno+1,lel)))
c ... tria3
        if( ty .eq. type_tria3) then
          if(flagTria3) ntria3(4) = tolel    
          ntria3(3) = ntria3(3) + 1
          tolel     = tolel + 1
          flagTria3 = .false.
c .....................................................................
c
c ... quad4
        else if( ty .eq. type_quad4) then
          if(flagQuad4) nquad4(4) = tolel
          nquad4(3) = nquad4(3) + 1
          tolel     = tolel + 1
          flagQuad4 = .false.     
c .....................................................................
c
c ... tetra4
        else if( ty .eq. type_tetra4) then
          if(flagTetra4) ntetra4(4) = tolel
          ntetra4(3) = ntetra4(3) + 1
          tolel      = tolel + 1
          flagTetra4 = .false.     
c .....................................................................
c
c ... hexa8
        else if( ty .eq. type_hexa8) then
          if(flagHexa8) nhexa8(4) = tolel
          nhexa8(3) = nhexa8(3) + 1
          tolel     = tolel + 1
          flagHexa8 = .false.     
c .....................................................................
c
c ... prism6
        else if( ty .eq. type_prism6) then
          if(flagPrism6) nPrism6(4) = tolel
          nPrism6(3) = nPrism6(3) + 1
          tolel      = tolel + 1
          flagPrism6 = .false. 
c .....................................................................
c
c ... tetra10
        else if( ty .eq. type_tetra10) then
          if(flagtetra10) ntetra10(4) = tolel
          ntetra10(3) = ntetra10(3) + 1
          tolel      = tolel + 1
          flagTetra10= .false.   
c .....................................................................
c
c ... hexa20
        else if( ty .eq. type_hexa20) then
          if(flagHexa20) nhexa20(4) = tolel
          nhexa20(3)  = nhexa20(3) + 1
          tolel       = tolel + 1
          flaghexa20 = .false.       
c .....................................................................    
        endif  
c .....................................................................
      enddo
c =====================================================================
c
c === gerando as coordenados locais
      if( npes .gt. 1 )then
        i_lno = alloc_8('lno      ',  ndm,my_nnode)
        i_g   = alloc_8('g        ',  ndm,nnode)
        call global_double(ia(i_g),x,nodedist,nnode,ndm,npes,mype
     .                    ,MPIW)   
        call localcoor(ia(i_g),ia(i_lno),noLG,my_nnode,ndm)
        i_g  = dealloc('g        ')
      else
        i_lno = alloc_8('lno      ',  ndm,my_nnode)
        call localcoor(x,ia(i_lno),noLG,my_nnode,ndm)
      endif
c =====================================================================
c
c === gerando a malha local
      if( npes .gt. 1 )then
        i_lel = alloc_4('lep     ',maxno+1,my_numel)
        i_g   = alloc_4('g       ',maxno+1,numel)
        call global_int(ia(i_g),el,elmdist,numel,maxno+1,npes,mype
     .                 ,MPIW)   
        call localel(ia(i_g),ia(i_lel),elLG,noGL,my_numel,maxno)
        i_g  = dealloc('g        ')
      else
        i_lel = alloc_4('lep      ',maxno+1,my_numel)
        call localel(el,ia(i_lel),elLG,noGL,my_numel,maxno)
      endif
c =====================================================================
c
c === gerando as propriedades dos materias
      i_p = alloc_4('p        ',    1,my_numel)
      call copy_mat(ia(i_p),ia(i_lel),my_numel,maxno)
c =====================================================================
c
c === 
c
c ... head
      write(aux,'(30a)')'Part mesh'
      if(legacy) then
        call head_vtk(aux,bvtk,0.0,0,.false.,nout)
      else
        call head_vtu(my_nnode,my_numel,bvtk,nout)
      endif
c ... escrevendo coordenadas
c     print*,'coor' 
      if(legacy) then
        call coor_vtk(ia(i_lno),my_nnode,ndm,bvtk,nout)
      else
        call coor_vtu(ia(i_lno),my_nnode,ndm,bvtk,nout)
      endif  
c ... escrevendo elementos      
c     print*,'elmt' 
      if(legacy) then
        call  elm_part_vtk(ia(i_lel),my_numel,maxno,bvtk,nout)
      else
        call elm_vtu(ia(i_lel),my_numel,maxno,bvtk,nout)
      endif  
c ... escrevendo propriedades dos elementos    
      if(legacy) then
        call cell_data_vtk(my_numel,bvtk,nout)
      else
        call cell_data_vtu(bvtk,nout)
      endif
c ... materias      
      write(aux1,'(15a)')'mat'
      cod = 1
      if(legacy) then
        call cell_prop_vtk(ia(i_p),fdum,ddum,my_numel,aux1,cod,1,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,my_numel,aux1,cod,1,bvtk
     .                    ,nout)
      endif 
c ... partcionamento      
      write(aux1,'(15a)')'part'
      do i = 1, my_numel
        ia(i_p+i-1) = ep(elLG(i)) 
      enddo
      if(legacy)then
        call cell_prop_vtk(ia(i_p),fdum,ddum,my_numel,aux1,cod,1,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,my_numel,aux1,cod,1,bvtk
     .                    ,nout)
      endif
c .....................................................................
c
c ...
      if(legacy .eqv. .false.) then
        call cell_data_finalize_vtu(bvtk,nout)
        call finalize_vtu(bvtk,nout)
      endif  
c .....................................................................
c =====================================================================
c
c === 
      ntria3(1)   = nt1
      ntria3(2)   = nt2
      nquad4(1)   = nq1
      nquad4(2)   = nq2
      ntetra4(1)  = ntr1
      ntetra4(2)  = ntr2
      nhexa8(1)   = nh1
      nhexa8(2)   = nh2
      nprism6(1)  = np1
      nprism6(2)  = np2
      ntetra10(1) = ntrq1
      ntetra10(2) = ntrq2
      nhexa20(1)  = nhq1
      nhexa20(2)  = nhq2
      i_p   = dealloc('p       ')
      i_lel = dealloc('lep     ')
      i_lno = dealloc('lno     ')
      close(nout)
      return
      end
c =====================================================================
c
c *********************************************************************
c
c *********************************************************************
c * Data de criacao    : 00/00/0000                                   *
c * Data de modificaco : 14/05/2017                                   *
c * ----------------------------------------------------------------- *      
c * WRITEMESHGEO: escreve a malha global com carregamento do formato  *
c * vtk.                                                              *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el         - conectividade com material                           *
c *  x         - coordenadas                                          *
c * ie         - tipo do elemento                                     *
c * id         - restricoes mecanico                                  *
c * nload      - cargas nos nos (mecanico)                            *
c * eload      - cargas nos elementos (poro-mecanico-mecanico)        *
c * eloadp      - cargas nos elementos (poro-mecanico-hidraulico)     *
c * f          - carregamento mecanicos                               *
c * u0         - cargas mecanicas iniciais                            *
c * idt        - restricoes termico                                   *
c * nloadt     - cargas nos nos (termico)                             *
c * eloadt     - cargas nos elementos (termico)                       *
c * ft         - carregamento termicos                                *
c * ut0        - cargas termicas iniciais                             *
c * nnode      - numero de nos                                        *
c * numel      - numero de elementos                                  *
c * ndf        - graus de liberdade mecanicas                         *
c * ndft       - graus de liberdade termicas                          *
c * nen        - numero de nos por elementos                          *
c * ndm        - numero de dimensoes                                  *
c * filein     - prefix do arquivo de saida                           *
c * bvtk       - true BINARY vtk false ASCII vtk                      *
c * macros     - macros lidas pela rdata                              *
c * legacy     - true (formato padrão .vtk) false (formato xlm .vtu) *
c * nout       - arquivo de saida principal                           *
c * nelemtload - arquivo de saida para face com cargas(elloads)       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine write_mesh_geo(el  ,x   ,ie 
     1                       ,id  ,nload ,eload  ,eloadp
     2                       ,f   ,u0    ,tx0
     3                       ,idt ,nloadt,eloadt
     4                       ,ft  ,ut0
     5                       ,nnode,numel,ndf,ndft
     6                       ,nen,ndm,filein,bvtk
     7                       ,macros,legacy,nout,nelemtload)
c ===
      use Malloc 
      implicit none
      include 'elementos.fi'
c ... variaveis da malha
      integer nnode,numel,nen,ndm,ntn,ntn1
      integer el(nen+1,numel),ie(*)
      real*8  x(ndm,nnode)
      integer nel,nno
c ... variaveis do problema
c     mecanico
      integer ndf
      integer id(ndf,*),nload(ndf,*),eload(*),eloadp(*)
      real*8  f(ndf,*),u0(ndf,*),tx0(*)
c     termico      
      integer ndft
      integer idt(ndft,*),nloadt(ndf,*),eloadt(*)
      real*8  ft(ndft,*),ut0(ndft,*)
c ... locais     
      integer*8 i_p,i_b1,i_b2
      data i_p/1/ i_b1/1/ i_b2/1/
      character*15 aux1
      character*30 aux
c ... faces
      integer lineface,triaface,quadface,nface
c ... macrp lida pela rdat      
      character*15 macros(*)
c ... variaveis dums
      real*8 ddum
      real*4 fdum
      integer idum
c ...
      character*8 malloc_name
c ... arquivo      
      integer nout,nelemtload
      character*80 fileout,fileelmtload,name,filein
      logical bvtk,legacy
      integer cod,cod2,gdl
c =====================================================================
c ... auxiliares
      integer nmacros,nmc,j,i
      character*15 macro1(18),macro2(18),rc
c ......................................................................
      data macro1/'constraindisp  ','nodalforces    ','nodalloads     ',
     1           '               ','               ','               ',
     2           'constraintemp  ','nodalsources   ','nodalthermloads',
     3           'constrainpmec  ','initialpres    ','initialstress  ',
     4           '               ','               ','               ',
     5           'initialtemp    ','               ','end            '/
c      
      data macro2/'elmtthermloads ','elmtloads      ','elmtpresloads ',
     1            '               ','               ','               ',
     2            '               ','               ','               ',
     3            '               ','               ','               ',
     4            '               ','               ','               ',
     5            '               ','               ','end            '/
      data nmc /18/
c ......................................................................
c
c ...
      if( ndm .eq. 2) then
        ntn = 4
      else
        ntn = 6
      endif
c ......................................................................
c
c ===
      if(legacy) then
        fileout = name(filein,0,110)
      else
        fileout = name(filein,0,117)
      endif
      if(bvtk)then
        open(unit=nout,file=fileout,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nout,file=fileout)
      endif  
c =====================================================================
c
c === cabecalho
      if(legacy) then
        write(aux,'(30a)')"Geometria lida" 
        call head_vtk(aux,bvtk,0.0,0,.false.,nout)
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
c === point data
      if(legacy) then
        call point_data_vtk(nnode,bvtk,nout)
      else  
        call point_data_vtu(bvtk,nout)
      endif  
c .....................................................................
c
c ... nos com numeracao global     
      write(aux1,'(15a)')"gNode"
      gdl =  1   
      cod =  1
      cod2 = 1
      malloc_name = 'p'
      i_p = alloc_4(malloc_name,1,nnode)
      do nno = 1, nnode
        ia(i_p+nno-1) = nno
      enddo
      if(legacy) then
        call point_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      endif
      i_p = dealloc(malloc_name)
c .....................................................................
c
c ...
      nmacros = 0 
  100 continue
      nmacros = nmacros + 1 
      write(rc,'(15a)')macros(nmacros)
      print*,rc
      do 200 j = 1, nmc
        if(rc .eq. macro1(j)) go to 300
  200 continue
      goto 100
c .....................................................................
  300 continue
      goto( 400, 450, 500   !constraindisp,nodalforces ,nodalloads
     1    , 550, 600, 650   !             ,            ,   
     2    , 700, 750, 800   !constraintemp,nodalsources,nodalthermloads
     3    , 850, 900, 950   !constrainpmec,initialpres ,initialstress
     4    ,1000,1050,1100   !             ,            ,   
     5    ,1150,1200,1250)j !initialtemp  ,            ,end            
c .....................................................................
c
c ... constriandisp
  400 continue
c ... gdb graus de liberdade
c     cod  1 vetorial
c     cod2 1 int(4bytes)
      write(aux1,'(15a)')"constraindisp"
      gdl =  ndf 
      cod =  2
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(id,fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call point_prop_vtu(id,fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      goto 100
c .....................................................................
c
c ... nodalforces  
  450 continue
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 3 real(8bytes) 
      write(aux1,'(15a)')"forces"
      gdl =  ndf
      cod =  2
      cod2 = 3
      if(legacy) then
        call point_prop_vtk(idum,fdum,f,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call point_prop_vtu(idum,fdum,f,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      goto 100
c .....................................................................      
c 
c ... nodalloads
  500 continue
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 interio(4bytes) 
      write(aux1,'(15a)')"nodalloads"
      gdl =  ndf
      cod =  2
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(nload,fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call point_prop_vtu(nload,fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      goto 100
c ......................................................................
c
c ......................................................................
  550 continue
  600 continue
  650 continue
      goto 100
c ......................................................................
c          
c ... constraintemp
  700 continue
      write(aux1,'(15a)')"constraintemp"
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 int(4bytes) 
      gdl =  ndft
      cod =  1
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(idt,fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call point_prop_vtu(idt,fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      goto 100
c ......................................................................
c
c ...     
  750 continue
c ... nodalsources          
      write(aux1,'(15a)')"sources"
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 3 real(8bytes) 
      gdl =  ndft
      cod =  1
      cod2 = 3
      if(legacy)then
        call point_prop_vtk(idum,fdum,ft,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call point_prop_vtu(idum,fdum,ft,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      goto 100
c .....................................................................
c
c ... nodalthermloads
  800 continue
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 interio(4bytes) 
      write(aux1,'(15a)')"nodalthermloads"
      gdl =  ndft
      cod =  1
      cod2 = 1
      if(legacy)then
        call point_prop_vtk(nloadt,fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      else
        call point_prop_vtu(nloadt,fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                     ,cod2  ,bvtk,nout)
      endif
      goto 100
c .....................................................................
c
c ...
  850 continue
      write(aux1,'(15a)')"constrainpmec"
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 int(4bytes) 
      gdl =  ndf 
      cod =  1
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(id,fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      else
        call point_prop_vtu(id,fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      goto 100
c .....................................................................
c
c ...
  900 continue
      write(aux1,'(15a)')"initialpres"
      malloc_name = 'p'
      i_p = alloc_8(malloc_name,1,nnode)
      call get_pres(u0,ia(i_p),el,nnode,numel,nen,ndf,.true.)
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 3 real(8bytes) 
      gdl =  1   
      cod =  1
      cod2 = 3
      if(legacy) then
        call point_prop_vtk(idum,fdum,ia(i_p),nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      else
        call point_prop_vtu(idum,fdum,ia(i_p),nnode,aux1,ndm,gdl,cod
     .                     ,cod2,bvtk,nout)
      endif
      i_p = dealloc(malloc_name)
      goto 100
c .....................................................................
c
c ...
  950 continue
c ... gerando o tensor completo 
      if( ntn .eq. 6 ) then
        malloc_name ='tensor'
        ntn1        = 9
        i_p         = alloc_8(malloc_name,ntn1 ,nnode)
        call make_full_tensor(tx0,ia(i_p),nnode,ntn, ntn1 )
      endif
c .....................................................................
c
c ... initialstress          
      write(aux1,'(15a)')"initialstress"
c ... gdb graus de liberdade
c     cod  3 tensor     
c     cod2 3 real(8bytes) 
      gdl  =  ntn1              
      cod  =  3
      cod2 =  3
      if(legacy)then
        call point_prop_vtk(idum,fdum,ia(i_p),nnode,aux1,ndm,gdl
     .                    ,cod ,cod2,bvtk   ,nout)
      else
        call point_prop_vtu(idum,fdum,ia(i_p),nnode,aux1,ndm,gdl
     .                    ,cod ,cod2,bvtk   ,nout)
      endif
c .....................................................................
c
c ... dealloc
      i_p = dealloc(malloc_name)
c .....................................................................
      goto 100
c .....................................................................
c
c ...
 1000 continue
 1050 continue
 1100 continue
      goto 100
c .....................................................................
c
c ...
 1150 continue
c ... intialtemp            
      write(aux1,'(15a)')"initialtemp"
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 real(8bytes) 
      gdl =  ndft
      cod =  1
      cod2 = 3
      if(legacy)then
        call point_prop_vtk(idum,fdum,ut0,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call point_prop_vtu(idum,fdum,ut0,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      goto 100
c .....................................................................
c
c ...
 1200 continue
      goto 100
c ... end      
 1250 continue
      if(legacy .eqv. .false.) call point_data_finalize_vtu(bvtk,nout)
c .....................................................................
c =====================================================================
c
c === cell data
      if(legacy) then
        call cell_data_vtk(numel,bvtk,nout)
      else
        call cell_data_vtu(bvtk,nout)
      endif
c ... materiais   
      malloc_name = 'p'   
      i_p = alloc_4(malloc_name, 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = el(nen+1,nel)
      enddo
      write(aux1,'(15a)')"mat" 
c ... cod = 1 variaveis interias
      cod = 1
      if(legacy)then
        call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,1,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,cod,1,bvtk
     .                    ,nout)
      endif
      i_p = dealloc(malloc_name)
c ... elementos com numeracao global     
      write(aux1,'(15a)')"gElm"
      gdl =  1   
      cod =  1
      cod2 = 1
      malloc_name = 'p'
      i_p = alloc_4(malloc_name, 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = nel
      enddo
      if(legacy) then
        call point_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,ndm,gdl,cod
     .                     ,cod2   ,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,ndm,gdl,cod
     .                     ,cod2   ,bvtk,nout)
      endif
      i_p = dealloc(malloc_name)
c .....................................................................
c =====================================================================
c ...
      nmacros = 0 
  101 continue
      nmacros = nmacros + 1 
      write(rc,'(15a)')macros(nmacros)
      print*,rc
      do 201 j = 1, nmc
        if(rc .eq. macro2(j)) go to 301
  201 continue
      goto 101
c .....................................................................
  301 continue
      goto( 401, 451, 501
     .    , 551, 601, 651
     .    , 701, 751, 801
     .    , 851, 901, 951  
     .    ,1001,1051,1101  
     .    ,1151,1201,1251)j
c .....................................................................
c
c ... elmtthermloads
  401 continue
c ... gdb graus de liberdade
c     cod  1 interio(4bytes)
      write(aux1,'(15a)')"elmtthermloads"
      gdl =  7   
      cod =  1
      if(legacy)then
        call cell_prop_vtk(eloadt,fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(eloadt,fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif
c .....................................................................
c
c ... arquivo auxiliar com as faces
      if(legacy) then
        fileelmtload = name(filein,0,118)
      else
        fileelmtload = name(filein,0,119)
      endif
      if(bvtk)then
        open(unit=nelemtload,file=fileelmtload,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nelemtload,file=fileelmtload )
      endif
c .....................................................................
c
c ... criando a estrutura para as faces 
      malloc_name = 'face'
      i_p  = alloc_4(malloc_name, maxface*maxnoface,numel)
      malloc_name = 'carga'
      i_b1 = alloc_4(malloc_name,       1,numel)
      malloc_name = 'faceType'
      i_b2 = alloc_4(malloc_name,       1,numel)
      call make_face(el      ,ie      ,eloadt  ,ia(i_p) ,ia(i_b1) 
     .              ,ia(i_b2),numel   ,nen     ,maxface ,maxnoface
     .              ,lineface,triaface,quadface)
c .....................................................................
c
c ... cabecalho
      nface = lineface + triaface + quadface 
      write(aux,'(30a)')"face com carregamento"
      if(legacy) then
        call head_vtk(aux,bvtk,0.0,0,.false.,nelemtload)
      else
        call head_vtu(nnode,nface,bvtk,nelemtload)
      endif
c .....................................................................
c
c ... Coordenadas
      if(legacy) then
        call coor_vtk(x,nnode,ndm,bvtk,nelemtload)
      else  
        call coor_vtu(x,nnode,ndm,bvtk,nelemtload)
      endif  
c .....................................................................
c
c ... faces
      if(legacy) then
        call face_vtk(ia(i_p)  ,ia(i_b1),ia(i_b2)
     .               ,maxnoface,lineface,triaface
     .               ,quadface ,bvtk    ,nelemtload)
      else
        call face_vtu(ia(i_p)  ,ia(i_b1),ia(i_b2)
     .               ,maxnoface,lineface,triaface
     .               ,quadface ,bvtk    ,nelemtload)
      endif
c .....................................................................
c
c ... cargas 
      if(legacy) then
        call cell_data_vtk(nface,bvtk,nelemtload)
      else
        call cell_data_vtu(bvtk,nelemtload)
      endif
c
      write(aux1,'(15a)')"cargas"
      gdl =  1   
      cod =  1
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(ia(i_b1),fdum,ddum,nface,aux1,ndm,gdl,cod
     .                    ,cod2    ,bvtk,nelemtload)
      else
        call point_prop_vtu(ia(i_b1),fdum,ddum,nface,aux1,ndm,gdl,cod
     .                    ,cod2    ,bvtk,nelemtload)
      endif
c .....................................................................
c
c ...
      if(legacy .eqv. .false.) then
        call cell_data_finalize_vtu(bvtk,nelemtload)
        call finalize_vtu(bvtk,nelemtload)
      endif  
c .....................................................................
c
c ...  
      malloc_name = 'face'
      i_p  = dealloc(malloc_name)
      malloc_name = 'carga'
      i_b1 = dealloc(malloc_name)
      malloc_name = 'faceType'
      i_b2 = dealloc(malloc_name)
      close(nelemtload)
c .....................................................................
      goto 101
c .....................................................................
c
c ...
  451 continue
c ... gdb graus de liberdade
c     cod  1 interio(4bytes)
      write(aux1,'(15a)')"elmtloads"
      gdl =  7   
      cod =  1
      if(legacy)then
        call cell_prop_vtk(eload,fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(eload,fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif
      goto 101
c .....................................................................      
c
c ...
  501 continue
c ... gdb graus de liberdade
c     cod  1 interio(4bytes)
      write(aux1,'(15a)')"elmtpresloads"
      gdl =  7   
      cod =  1
      if(legacy)then
        call cell_prop_vtk(eloadp,fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(eloadp,fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif
      goto 101
c .....................................................................      
c
c ...
  551 continue
  601 continue
  651 continue
  701 continue
  751 continue
  801 continue
  851 continue
  901 continue
  951 continue
 1001 continue
 1051 continue
 1101 continue
 1151 continue
 1201 continue
      goto 101 
c .....................................................................
c
c ... end      
 1251 continue
      if(legacy .eqv. .false.) then
        call cell_data_finalize_vtu(bvtk,nout)
        call finalize_vtu(bvtk,nout)
      endif  
c ..................................................................... 
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
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 27/03/2016                                    *
c * Data de modificaco : 25/10/2016                                    * 
c * ------------------------------------------------------------------ *   
c * GET_PRES : obtem as presoes                                        *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * u (ndf,*)    - deslocamento e pressoes                             *
c * pres(*)      - indefinido                                          *
c * el(maxno+1,*)- conectividade                                       *
c * nnode        - nos de nos ( verticeis ou total)                    *
c * numel        - numero de elmentos                                  *
c * maxno        - numero de nos maximo por elemento                   *
c * ndf          - grau de liberdade                                   *
c * quad         - pressoes nos nos intermediarios por interpolacao    *
c *                linear (true| false)                                *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * pres(*)    - pressoes                                              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine get_pres(u,pres,el,nnode,numel,maxno,ndf,quad)
      implicit none
      integer maxEdge
      parameter (maxEdge = 12) 
c ...
      integer i,j
      integer nnode,numel,ndf,nedge,no1,no2,no3,maxno
      integer el(maxno+1,*)
      integer iEdge(3,maxEdge)
      real*8 u(ndf,*),pres(*)
      logical quad
c .....................................................................  
c
c ... pressa no nos lineres
      do i = 1, nnode
        pres(i) = u(ndf,i)
      enddo
c .....................................................................  
c
c ... pressao no nos quadraticos
      if(quad) then
c ...
        nedge = 0  
c ... tetraedros de 10 nos 
        if( maxno .eq. 10 ) then
          nedge =  6
          call tetra10edgeNod(iEdge) 
c ... hexaedros de 20 nos 
        else if( maxno .eq. 20 ) then
          nedge = 12
          call hexa20edgeNod(iEdge) 
        endif
c ...
        do i = 1, numel
          do j = 1, nedge
c ... no vertices
            no1     = el(iEdge(1,j),i)
            no2     = el(iEdge(2,j),i)
c ... no central
            no3     = el(iEdge(3,j),i)
            pres(no3) = 0.5d0*(u(ndf,no1) +  u(ndf,no2))
          enddo
        enddo
      endif  
c .....................................................................

      return
      end
c *********************************************************************
c       
c *********************************************************************
c * make_full_tensor : transforma um tensor simetrico em um tensor    *
c * geral                                                             *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * tx     - tensor simetrico                                         *
c * tensor - indefinido                                               *
c * nnode  - numero de nos de vertices                                *
c * ntn    - numero total de termos do tensor simetrico               *
c * n      - numero total de termos do tensor geral                   *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * tensor - tensor geral                                             *
c * ----------------------------------------------------------------- *
c * OBS:                                                              *
c * tx    (Sxx,Syy,Szz, Sxy, Syz, Sxz)                                *
c * tensor(Sxx,Sxy,Sxz, Syx, Syy, Syz, Szx, Syz, Szz)                 *
c *********************************************************************      
      subroutine make_full_tensor(tx,tensor,nnode,ntn,n)
      implicit none
      real*8 tx(ntn,*),tensor(n,*)
      integer nnode,ntn,n,i
c ...
      if(ntn .eq. 4) then
        do i = 1, nnode
c ... sgima xx
          tensor(1,i) = tx(1,i)
c ... sgima xy
          tensor(2,i) = tx(4,i)
c ... sgima xz
          tensor(3,i) = 0.d0   
c ... sgima yx
          tensor(4,i) = tx(4,i)
c ... sgima yy
          tensor(5,i) = tx(2,i)
c ... sgima yz
          tensor(6,i) = 0.0d0  
c ... sgima zx
          tensor(7,i) = 0.0d0
c ... sgima zy
          tensor(8,i) = 0.0d0  
c ... sgima zz
          tensor(9,i) = tx(3,i)  
        enddo
c ...
      else if(ntn .eq. 6) then
        do i = 1, nnode
c ... sgima xx
          tensor(1,i) = tx(1,i)
c ... sgima xy
          tensor(2,i) = tx(4,i)
c ... sgima xz
          tensor(3,i) = tx(6,i)
c ... sgima yx
          tensor(4,i) = tx(4,i)
c ... sgima yy
          tensor(5,i) = tx(2,i)
c ... sgima yz
          tensor(6,i) = tx(5,i)
c ... sgima zx
          tensor(7,i) = tx(6,i)
c ... sgima zy
          tensor(8,i) = tx(5,i)
c ... sgima zz
          tensor(9,i) = tx(3,i)
        enddo
      endif
c .....................................................................
c
c ...
      return
      end
c **********************************************************************
