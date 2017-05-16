c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 14/05/2017                                    * 
c * ------------------------------------------------------------------ * 
c * WRITEMYDATPART : escreve o arquivo no formato do mefpar            *
c * ------------------------------------------------------------------ *
c * Parametros de entrada :                                            *
c * ------------------------------------------------------------------ *
c * el(maxno+1,numel) - conectividade global                           *
c * x(ndm,nnode)      - coordenadas globais                            *
c * elGL(numel)       - map global -> local de elementos               *
c * elLG(my_numel)    - map local -> global de elementos               *
c * noLG(my_nnode)    - map local -> global de nos                     *
c * noGL(nnode)       - map global -> local de nos                     *
c * ranks             -                                                *
c * sizes             -                                                *
c * fmap              -                                                *
c * nnof              -                                                *
c * rvcs              -                                                *
c * dspl              -                                                *
c * ie                -                                                *
c * e                 -                                                *
c * id                -                                                *
c * nload             -                                                *
c * eload             -                                                *
c * eloadp            -                                                *
c * f                 -                                                *
c * u0                -                                                *
c * tx0               -                                                *
c * idt               -                                                *
c * nloadt            -                                                *
c * eloadt            -                                                *
c * ft                -                                                *
c * ut0               -                                                *
c * nodedist          - distribuicao das propriedades nodais           *
c * elmdist           - distribuicao das propriedades por elementos    *
c * nnodev            - numero de nos dos vertices                     *   
c * nnode             - numero de nos global                           *
c * numel             - numero de elementos global                     *
c * maxnov            - numero maximo de nos de vertices por elem nto  *
c * maxno             - numero maximo de nos por elementos             *
c * ndm               - numero de dimensoes                            *
c * numat             - numero de materiais                            *
c * ndf               - graus de liberdade                             *
c * npi               - numero de pontos de integracao                 *
c * file_prop         - arquivos com propriedades do material          *
c * lines             - linhas apos a macro end mesh                   *
c * nlines            - numero de linha apos a macro end mesh          *
c * plines            - linhas pre macro mesh                          *
c * nplines           - numero de linha pre macro mesh                 *
c * filein            - nome do arquivo                                *
c * nout              -                                                *
c * npart             - numero da particao                             *
c * div               - numero da particao totais deste processo       *
c * my_rank           - id do processo (MPI)                           *
c * npes              - numero de processos (MPI)                      *
c * ------------------------------------------------------------------ *
c * Parametros de saida :                                              *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * Obs:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine  writemydatpart(el      ,x     ,elLG
     1                          ,noLG    ,noGL  ,ranks
     2                          ,sizes   ,fmap  ,nnof
     3                          ,rcvs    ,dspl  ,ie    ,e
     4                          ,id      ,nload ,eload ,eloadp
     5                          ,f       ,u0    ,tx0
     6                          ,idt     ,nloadt,eloadt
     7                          ,ft      ,ut0
     8                          ,nodedist,elmdist
     9                          ,nnodev  ,nnode,numel,maxnov,maxno,ndm
     1                          ,numat   ,ndf,ndft,npi,file_prop
     2                          ,my_nno1 ,my_nno2,my_nno3,my_nno4
     3                          ,my_nno1a,my_nnode
     4                          ,my_numel,my_numel_nov,my_numel_ov
     5                          ,lines   ,nlines,plines,nplines
     6                          ,macros  ,filein,nprcs,ovlp
     7                          ,npart   ,div,my_rank,npes,MPIW,nout)
c ===
      use Malloc
      implicit none
      include 'load.fi'
      include 'elementos.fi'
c ... malha      
      integer el(maxno+1,*)
      real*8 x(*)
      real*8 e(10,*)
      integer ie(*)
      integer maxnov,maxno,ndm,numel,nnodev,nnode,numat,npi,ntn
c ... problema termico      
      integer ndft
      integer idt(ndft,*),nloadt(*),eloadt(7,*)
      real*8 ft(ndft,*),ut0(ndft,*)
c ... problema mecanico
      integer ndf 
      real*8 f(ndf,*),u0(ndf,*),tx0(*)
      integer id(*),nload(*),eload(7,*),eloadp(7,*)
c ... particionamento
      integer nprcs
      integer my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a,my_nnode
      integer my_numel,my_numel_ov,my_numel_nov
      integer elLG(*),noLG(*),noGL(*)
      integer ranks(*),sizes(*)
      integer fmap(*),nnof(*)
      integer rcvs(*),dspl(*)
      logical ovlp
      integer*8 i_lel,i_xl,i_g1,i_g2
      data i_lel/1/, i_xl/1/, i_g1/1/, i_g2/1/ 
c ... Mpi
      integer my_rank,npes,MPIW,nodedist(*),elmdist(*)
      integer npart,div
c ... Arquivo
      integer nout
      character*80 filein,name,fileout
      character*200 lines(*),plines(*)
      character*15 macros(*)
      integer nlines,nplines
c ...
      character*80 file_prop(*)
c ... auxiliares
      integer nmacros,nmc,j,i,lel,ty,tolel
      logical flagtria3,flagquad4,flagtetra4,flaghexa8,flagprism6,
     .        flagtetra10,flaghexa20
      integer nt1,nq1,nh1,nt2,nq2,nh2,ntr1,ntr2,np1,np2,ntq1,ntq2,
     .        nhq1,nhq2
      character*15 macro(42),rc
      character*8 malloc_name
c ......................................................................
      data macro/'materials      ','bar2           ','tria3          ',
     1           'quad4          ','tetra4         ','hexa8          ',
     2           '               ','termprop       ','adiabat        ',
     3           'coordinates    ','constraindisp  ','constraintemp  ',
     4           'nodalforces    ','nodalsources   ','nodalloads     ',
     5           'nodalthermloads','elmtloads      ','elmtthermloads ',
     6           'loads          ','hidrprop       ','prism6         ',
     7           'initialdisp    ','initialtemp    ','velocityfield  ',
     8           'parallel       ','               ','               ',
     9           'coordinatesbin ','               ','               ',
     1           'quad4bin       ','tetra4bin      ','hexa8bin       ',
     2           'quad8bin       ','hydrostatic    ','hydrostress    ',
     3           'constrainpmec  ','initialpres    ','initialstress  ',
     4           'elmtpresloads  ','fmaterials     ','end            '/
      data nmc /42/
c ......................................................................

c =====================================================================
c
c ...
      if( ndm .eq. 2) then
        ntn = 4
      else
        ntn = 6
      endif
c ......................................................................
c
c === caso especial de processos ociosos
c     print*,div
      if (div .lt. 1) goto 15
      fileout = name(filein,npart,101)
      open(nout,file=fileout)
c =====================================================================
c
c === macros comandos pre a "mesh"
      do i = 1, nplines
        write(nout,'(a)')trim(plines(i))
      enddo  
c =====================================================================
c
c === escrevendo no arquivo
c
c ...
      call writeparameters(my_nnode,my_numel,numat,maxno,ndf,ndft,ndm
     .                    ,npi,nout)
c ....................................................................
c
c ... alocando memoria para as conectividades de elementos locais
      i_lel = alloc_4('lep      ',maxno+1,my_numel)
      call localel(el,ia(i_lel),elLG,noGL,my_numel,maxno)
c ....................................................................
c
c ... calulo da distribuicao de elementos por particao
      nt1           = ntria3(1)
      nt2           = ntria3(2)
      nq1           = nquad4(1)
      nq2           = nquad4(2)
      ntr1          = ntetra4(1)
      ntr2          = ntetra4(2)
      nh1           = nhexa8(1)
      nh2           = nhexa8(2)
      np1           = nprism6(1)
      np2           = nprism6(2)
      ntq1          = ntetra10(1)
      ntq2          = ntetra10(2)
      nhq1          = nhexa20(1)
      nhq2          = nhexa20(2)
      ntria3(1:4)   = 0
      nquad4(1:4)   = 0
      ntetra4(1:4)  = 0
      nhexa8(1:4)   = 0
      nprism6(1:4)  = 0
      ntetra10(1:4) = 0
      nhexa20(1:4)  = 0
      tolel         = 1
      flagTria3     = .true.
      flagQuad4     = .true.
      flagTetra4    = .true.
      flagHexa8     = .true.
      flagprism6    = .true.
      flagtetra10   = .true.
      flaghexa20    = .true.
c .....................................................................
c
c ...
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
          if(flagPrism6) nprism6(2) = tolel
          nprism6(1) = nprism6(1) + 1
          tolel      = tolel + 1
          flagPrism6 = .false.
c .....................................................................
c
c ... tetra10
        else if( ty .eq. type_tetra10) then
          if(flagtetra10) ntetra10(2) = tolel
          ntetra10(1) = ntetra10(1) + 1
          tolel       = tolel + 1
          flagtetra10 = .false.      
c .....................................................................
c
c ... hexa20
        else if( ty .eq. type_hexa20) then
          if(flaghexa20) nhexa20(2) = tolel
          nhexa20(1) = nhexa20(1) + 1
          tolel      = tolel + 1
          flaghexa20 = .false.      
c .....................................................................
c
c ...
        else
          print*,'Non-Existing Element Type:',ty
          call finalize() 
        endif 
c .....................................................................
      enddo
c ... elementos em overllaping
      flagtria3    = .true.
      flagquad4    = .true.
      flagtetra4   = .true.
      flaghexa8    = .true.
      flagprism6   = .true.
      flagtetra10  = .true.
      flaghexa20   = .true.
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
          if(flagPrism6) nprism6(4) = tolel
          nprism6(3) = nprism6(3) + 1
          tolel      = tolel + 1
          flagPrism6 = .false.     
c .....................................................................
c
c ... tetra10
        else if( ty .eq. type_tetra10) then
          if(flagtetra10) ntetra10(4) = tolel
          ntetra10(3) = ntetra10(3) + 1
          tolel       = tolel + 1
          flagtetra10 = .false.      
c .....................................................................
c
c ... hexa20
        else if( ty .eq. type_hexa20) then
          if(flaghexa20) nhexa20(4) = tolel
          nhexa20(3) = nhexa20(3) + 1
          tolel      = tolel + 1
          flaghexa20 = .false.      
c .....................................................................
c
c ...
        else
          print*,'Non-Existing Element Type:',ty
          call finalize() 
        endif 
c .....................................................................
      enddo
c .....................................................................  
      nmacros = 0 
  100 continue
      nmacros = nmacros + 1 
      write(rc,'(15a)')macros(nmacros)
      do 200 j = 1, nmc
        if(rc .eq. macro(j)) go to 300
  200 continue
      goto 100
c .....................................................................  
c
c .....................................................................
  300 continue
c .....................................................................
c
c .....................................................................
      goto (400 ,450 , 500,  !materials      ,bar2         ,tria3
     1      550 ,600 , 650,  !quad4          ,tetra4       ,hexa8
     2      700 ,750 , 800,  !quad8          ,termprop     ,adiabt
     3      850 ,900 , 950,  !coordinates    ,constraindisp,constraintemp
     4      1000,1050,1100,  !nodalforces    ,nodalsources ,nodalloads
     5      1150,1200,1250,  !nodalthermloads,elmtloads    ,elmtthermdloads
     6      1300,1350,1400,  !loads          ,hidprop      ,prism6
     7      1450,1500,1550,  !initialdisp    ,initialtemp  ,velocityfield
     8      1600,1650,1700,  !parallel       ,insert       ,return
     9      1750,1800,1850,  !coordinatesbin ,             ,
     1      1900,1950,2000,  !quad4bin       ,tetra4bin    ,hexa8bin
     1      2050,2100,2150,  !quad8bin       ,hydrostatic  ,hydrostress
     2      2200,2250,2300,  !constrainpmec  ,initialpres  ,initialstress
     3      2350,2400,2450)j !elmtpresloads  ,fmaterials   ,end
c .....................................................................
c
c ... Materiais
  400 continue
        call writemate(ie,e,numat,nout)
      goto 100
c .....................................................................
c
c ... conetividade bar2
  450 continue
        print*, "Escrita file.dat bar2 nao implementado."
      goto 100
c .....................................................................
c
c ... conectividade tria3
  500 continue
        if( npes .gt. 1 )then
          i_lel = alloc_4('lep      ',maxno+1,my_numel)
          i_g1  = alloc_4('g        ',maxno+1,numel)
          call global_int(ia(i_g1),el,elmdist,numel,maxno+1,npes,my_rank
     .                   ,MPIW)
          call localel(ia(i_g1),ia(i_lel),elLG,noGL,my_numel,maxno)
          call writeelmef(ia(i_lel),my_numel,maxno,nout)
          i_g1  = dealloc('g        ')
          i_lel = dealloc('lep      ')
        else
          if(ntria3(1) .gt. 0 ) then 
            write(nout,'(a)')"tria3"
            tolel = ntria3(2)+ntria3(1) - 1
            call writeelmef(ia(i_lel),ie,ntria3(2),tolel
     .                     ,maxno,3,type_tria3,nout)
            write(nout,'(a)')"end tria3"
          endif
c ... elementos em overlaping
          if(ntria3(3) .gt. 0) then
            write(nout,'(a)')"tria3ov"
            tolel = ntria3(3)+ntria3(4) - 1
            call writeelmef(ia(i_lel),ie,ntria3(4)
     .                     ,tolel,maxno,3,type_tria3,nout)
            write(nout,'(a)')"end tria3ov"
c .....................................................................
          endif
        endif
      goto 100
c .....................................................................
c
c ... conectividade quad4
  550 continue
        if(npes.gt.1)then
          i_lel = alloc_4('lep      ',maxno+1,my_numel)
          i_g1  = alloc_4('g        ',maxno+1,numel)
          call global_int(ia(i_g1),el,elmdist,numel,maxno+1,npes,my_rank
     .                   ,MPIW)
          call localel(ia(i_g1),ia(i_lel),elLG,noGL,my_numel,maxno)
          call writeelmef(ia(i_lel),my_numel,maxno,nout)
          i_g1  = dealloc('g        ')            
          i_lel = dealloc('lep      ')
        else
          if(nquad4(1) .gt. 0) then
            write(nout,'(a)')"quad4"
            tolel = nquad4(2)+nquad4(1) - 1
            call writeelmef(ia(i_lel),ie,nquad4(2),tolel
     .                     ,maxno,4,type_quad4,nout)
            write(nout,'(a)')"end quad4"
          endif
c ... elementos em overlaping
          if(nquad4(3) .gt. 0 ) then
            write(nout,'(a)')"quad4ov"
            tolel = nquad4(4)+nquad4(3) - 1
            call writeelmef(ia(i_lel),ie,nquad4(4)
     .                     ,tolel,maxno,4,type_quad4,nout)
            write(nout,'(a)')"end quad4ov"
c .....................................................................
          endif
        endif  
      goto 100
c .....................................................................
c
c ... conectividae tetra4
  600 continue
        if(npes.gt.1)then
          i_lel = alloc_4('lep      ',maxno+1,my_numel)
          i_g1  = alloc_4('g        ',maxno+1,numel)
          call global_int(ia(i_g1),el,elmdist,numel,maxno+1,npes,my_rank
     .                   ,MPIW)
          call localel(ia(i_g1),ia(i_lel),elLG,noGL,my_numel,maxno)
          call writeelmef(ia(i_lel),my_numel,maxno,nout)
          i_g1   = dealloc('g        ')            
          i_lel = dealloc('lep      ')
        else
          if(ntetra4(1) .gt. 0) then
            write(nout,'(a)')"tetra4"
            tolel = ntetra4(2)+ntetra4(1) - 1
            call writeelmef(ia(i_lel),ie,ntetra4(2)
     .                     ,tolel,maxno,4,type_tetra4,nout)
            write(nout,'(a)')"end tetra4"
          endif
c ... elementos em overlaping
          if(ntetra4(3) .gt. 0) then
            write(nout,'(a)')"tetra4ov"
            tolel = ntetra4(4)+ntetra4(3) - 1
            call writeelmef(ia(i_lel),ie,ntetra4(4)
     .                     ,tolel,maxno,4,type_tetra4,nout)
            write(nout,'(a)')"end tetra4ov"
          endif
c .....................................................................
c
c ...
          if(ntetra10(1) .gt. 0) then
            write(nout,'(a)')"tetra10"
            tolel = ntetra10(2)+ntetra10(1) - 1
            call writeelmef(ia(i_lel),ie,ntetra10(2)
     .                     ,tolel,maxno,10,type_tetra10,nout)
            write(nout,'(a)')"end tetra10"
          endif
c ........................................................................
c
c ... elementos em overlaping
          if(ntetra10(3) .gt. 0) then
            write(nout,'(a)')"tetra10ov"
            tolel = ntetra10(4)+ntetra10(3) - 1
            call writeelmef(ia(i_lel),ie,ntetra10(4)
     .                     ,tolel,maxno,10,type_tetra10,nout)
            write(nout,'(a)')"end tetra10ov"
          endif
c .....................................................................
        endif
      goto 100
c .....................................................................
c
c ... conectividade hexa8
  650 continue
c ...
        if(npes.gt.1)then
          i_lel = alloc_4('lep      ',maxno+1,my_numel)
          i_g1  = alloc_4('g        ',maxno+1,numel)
          call global_int(ia(i_g1),el,elmdist,numel,maxno+1,npes,my_rank
     .                   ,MPIW)
          call localel(ia(i_g1),ia(i_lel),elLG,noGL,my_numel,maxno)
          call writeelmef(ia(i_lel),my_numel,maxno,nout)
          i_g1   = dealloc('g        ')            
          i_lel = dealloc('lep      ')
c .....................................................................
c
c ...
        else
          if(nhexa8(1) .gt. 0) then
            write(nout,'(a)')"hexa8"
            tolel = nhexa8(2)+nhexa8(1) - 1
            call writeelmef(ia(i_lel),ie,nhexa8(2)
     .                     ,tolel,maxno,8,type_hexa8,nout)
            write(nout,'(a)')"end hexa8"
          endif
c ........................................................................
c
c ... elementos em overlaping
          if(nhexa8(3) .gt. 0) then
            write(nout,'(a)')"hexa8ov"
            tolel = nhexa8(4)+nhexa8(3) - 1
            call writeelmef(ia(i_lel),ie,nhexa8(4)
     .                     ,tolel,maxno,8,type_hexa8,nout)
            write(nout,'(a)')"end hexa8ov"
          endif
c .....................................................................
c
c ...
          if(nhexa20(1) .gt. 0) then
            write(nout,'(a)')"hexa20"
            tolel = nhexa20(2)+nhexa20(1) - 1
            call writeelmef(ia(i_lel),ie,nhexa20(2)
     .                     ,tolel,maxno,20,type_hexa20,nout)
            write(nout,'(a)')"end hexa20"
          endif
c ........................................................................
c
c ... elementos em overlaping
          if(nhexa20(3) .gt. 0) then
            write(nout,'(a)')"hexa20ov"
            tolel = nhexa20(4)+nhexa20(3) - 1
            call writeelmef(ia(i_lel),ie,nhexa20(4)
     .                     ,tolel,maxno,20,type_hexa20,nout)
            write(nout,'(a)')"end hexa20ov"
          endif
c .....................................................................
        endif  
      goto 100
c .....................................................................
c
c ...
  700 continue
      goto 100
c ... Propriedades variaveis com a temperatura:termprop
  750 continue
        write(nout,'(a)')'termprop'
        call writetermprop(numat,nout)
        write(nout,'(a)')'end termprop'
      goto 100
c .....................................................................
  800 continue
c ... Curva adiabatica        
        write(nout,'(a)')'adiabat'
        call writeadiabat(numat,nout)
        write(nout,'(a)')'end adiabat'
      goto 100
c .....................................................................
c
c ... Coordenadas
  850 continue
        write(nout,'(a)')'coordinates'
c ... escrevendo coordenandas no mpi (coordenadas espalhada pelos
c     processos
        if( npes .gt. 1 )then
          i_xl   = alloc_8('xl       ',ndm,my_nnode)
          i_g1   = alloc_8('Gcoor    ',  ndm,nnode)
          call azero(ia(i_g1),ndm*nnode)
c ... recupera os valores gobais da coordenadas 
          call global_double(ia(i_g1),x,nodedist,nnode,ndm,npes
     .                             ,my_rank,MPIW)   
          call localcoor(ia(i_g1),ia(i_xl),noLG,my_nnode,ndm)
          call writecoormef(ia(i_xl),my_nnode,ndm,nout)
          i_g1  = dealloc('Gcoor    ')
          i_xl  = dealloc('xl       ')
c ... escrevendo sequencial                  
        else
          i_xl = alloc_8('xl       ',ndm,my_nnode)
          call localcoor(x,ia(i_xl),noLG,my_nnode,ndm)
          call writecoormef(ia(i_xl),my_nnode,ndm,nout)
          i_xl = dealloc('xl       ')
        endif
        write(nout,'(a)')'end coordinates'
      goto 100
c .....................................................................
c
c ... Restricoes constraindisp(deslocamentos)
  900 continue
        write(nout,'(a)')'constraindisp'
        if( npes .gt. 1)then
          i_g1   = alloc_4('Gid      ',  ndf,nnode)
          call mzero(ia(i_g1),ndf*nnode)
c ... recupera os valores gobais de id          
          call global_int(ia(i_g1),id,nodedist,nnode,ndf,npes
     .                   ,my_rank,MPIW)   
          call wpbound(ia(i_g1),noLG,my_nnode,ndf,nout)   
          i_g1  = dealloc('Gid      ')
        else  
          call wpbound(id,noLG,my_nnode,ndf,nout)   
        endif
        write(nout,'(a)')'end contraindisp'
      goto 100
c .....................................................................
c
c ... Resticoes termicas(temperaturas)
  950 continue
        write(nout,'(a)')'constraintemp'
        if( npes .gt. 1)then
          i_g1   = alloc_4('Gidt     ',  ndft,nnode)
          call mzero(ia(i_g1),ndft*nnode)
c ... recupera os valores gobais de id          
          call global_int(ia(i_g1),idt,nodedist,nnode,ndft,npes
     .                   ,my_rank,MPIW)   
          call wpbound(ia(i_g1),noLG,my_nnode,ndft,nout)   
          i_g1  = dealloc('Gidt     ')
        else  
          call wpbound(idt,noLG,my_nnode,ndft,nout)
        endif  
        write(nout,'(a)')'end constraintemp'
      goto 100
c .....................................................................
c
c ... nodalforces - forcas nodais
 1000 continue
        write(nout,'(a)')'nodalforces'
        if( npes .gt. 1)then
          i_g1  = alloc_4('Gid      ',  ndf,nnode)
          i_g2  = alloc_8('Gf       ',  ndf,nnode)
          call mzero(ia(i_g1),ndf*nnode)
          call azero(ia(i_g2),ndf*nnode)
c ... recupera os valores gobais de id          
          call global_int(ia(i_g1),id,nodedist,nnode,ndf,npes
     .                   ,my_rank,MPIW)   
c ... recupera os valores gobais de f          
          call global_double(ia(i_g2),f,nodedist,nnode,ndf,npes
     .                      ,my_rank,MPIW)   
          call wpforces(ia(i_g2),ia(i_g1),noLG
     .                 ,my_nnode,my_nno1,my_nno2,ndf,nout)
          i_g2  = dealloc('Gf       ')
          i_g1  = dealloc('Gid      ')
        else  
          call wpforces(f,id,noLG
     .                 ,my_nnode,my_nno1,my_nno2,ndf,nout)
        endif
        write(nout,'(a)')'end nodalforces'
      goto 100
c .....................................................................
c
c ... nodalsources - fontes termicas nodais:
 1050 continue
        write(nout,'(a)')'nodalsources'
        if( npes .gt. 1)then
          i_g1  = alloc_4('Gidt     ',  ndft,nnode)
          i_g2  = alloc_8('Gft      ',  ndft,nnode)
          call mzero(ia(i_g1),ndft*nnode)
          call azero(ia(i_g2),ndft*nnode)
c ... recupera os valores gobais de id          
          call global_int(ia(i_g1),idt,nodedist,nnode,ndft,npes
     .                   ,my_rank,MPIW)   
c ... recupera os valores gobais de f          
          call global_double(ia(i_g2),ft,nodedist,nnode,ndft,npes
     .                      ,my_rank,MPIW)   
          call wpforces(ia(i_g2),ia(i_g1),noLG
     .                 ,my_nnode,my_nno1,my_nno2,ndft,nout)
          i_g2  = dealloc('Gft      ')
          i_g1  = dealloc('Gidt     ')
        else  
          call wpforces(ft,idt,noLG
     .                 ,my_nnode,my_nno1,my_nno2,ndft,nout)
        endif
        write(nout,'(a)')'end nodalsources'
      goto 100
c .....................................................................
c
c ...nodalloads - nos com cargas cariaveis no tempo:
 1100 continue
        write(nout,'(a)')'nodalloads'
        if( npes .gt. 1)then
          i_g1  = alloc_4('Gid      ',  ndf,nnode)
          i_g2  = alloc_4('Gnload   ',  ndf,nnode)
          call mzero(ia(i_g1),ndf*nnode)
          call mzero(ia(i_g2),ndf*nnode)
c ... recupera os valores gobais de id          
          call global_int(ia(i_g1),id,nodedist,nnode,ndf,npes
     .                   ,my_rank,MPIW)   
c ... recupera os valores gobais de nload       
          call global_int(ia(i_g2),nload,nodedist,nnode,ndf,npes
     .                   ,my_rank,MPIW)   
          call wploads(ia(i_g2),ia(i_g1),noLG,my_nnode,my_nno1,my_nno2
     .                ,ndf,nout)
          i_g2 = dealloc('Gnload   ')
          i_g1 = dealloc('Gid      ')
        else  
          call wploads(nload,id,noLG,my_nnode,my_nno1,my_nno2,ndf,nout)
        endif  
        write(nout,'(a)')'end nodalloads'
      goto 100
c .....................................................................
c
c ...nodalthermloads - nos com cargas termicas variaveis no tempo:
 1150 continue
        write(nout,'(a)')'nodalthermloads'
        if( npes .gt. 1)then
          i_g1  = alloc_4('Gidt     ', ndft,nnode)
          i_g2  = alloc_4('Gnloadt  ', ndft,nnode)
          call mzero(ia(i_g1),ndft*nnode)
          call mzero(ia(i_g2),ndft*nnode)
c ... recupera os valores gobais de id          
          call global_int(ia(i_g1),idt,nodedist,nnode,ndft,npes
     .                   ,my_rank,MPIW)   
c ... recupera os valores gobais de nload       
          call global_int(ia(i_g2),nloadt,nodedist,nnode,ndft,npes
     .                   ,my_rank,MPIW)   
          call wploads(ia(i_g2),ia(i_g1),noLG,my_nnode,my_nno1,my_nno2
     .                ,ndft,nout)
          i_g2 = dealloc('Gnloadt  ')
          i_g1 = dealloc('Gidt     ')
        else  
          call wploads(nloadt,idt,noLG,my_nnode,my_nno1,my_nno2,ndft
     .                ,nout)
        endif
        write(nout,'(a)')'end nodalthermloads'
      goto 100
c .....................................................................
c
c ... elmtloads - cargas nos elementos:
 1200 continue
        write(nout,'(a)')'elmtloads'
        if( npes .gt. 1)then
          i_g1  = alloc_4('Geload   ',    7,numel)
          call mzero(ia(i_g1),7*numel)
c ... recupera os valores gobais de eload          
          call global_int(ia(i_g1),eload,elmdist,numel,7,npes
     .                   ,my_rank,MPIW)   
          call wpelloads(ia(i_g1),elLG,my_numel,my_numel_ov,nout)
          i_g1 = dealloc('Geload   ')
        else  
          call wpelloads(eload,elLG,my_numel,my_numel_ov,nout)
        endif  
        write(nout,'(a)')'end elmtloads'
      goto 100
c .....................................................................
c
c ... elmtthermloads - cargas termicas nos elementos
 1250 continue
        write(nout,'(a)')'elmtthermloads'
        if( npes .gt. 1)then
          i_g1  = alloc_4('Geloadt  ',    7,numel)
          call mzero(ia(i_g1),7*numel)
c ... recupera os valores gobais de eload          
          call global_int(ia(i_g1),eloadt,elmdist,numel,7,npes
     .                   ,my_rank,MPIW)   
          call wpelloads(ia(i_g1),elLG,my_numel,my_numel_ov,nout)
          i_g1 = dealloc('Geloadt  ')
        else  
          call wpelloads(eloadt,elLG,my_numel,my_numel_ov,nout)
        endif  
        write(nout,'(a)')'end elmtthermloads'
      goto 100
c .....................................................................
c
c ... Definicao das cargas variaveis no tempo: (loads)
 1300 continue
        write(nout,'(a)')'loads'
        call wload(nout)
        write(nout,'(a)')'end loads'
      goto 100
c .....................................................................
c
c ... hidrprop
 1350 continue
        write(nout,'(a)')'hidrprop'
        call writetermprop(numat,nout)
        write(nout,'(a)')'end hidrprop'
      goto 100
c .....................................................................
c
c ... prism6
 1400 continue
        if( npes .gt. 1 )then
          i_lel = alloc_4('lep      ',maxno+1,my_numel)
          i_g1  = alloc_4('g        ',maxno+1,numel)
          call global_int(ia(i_g1),el,elmdist,numel,maxno+1,npes,my_rank
     .                   ,MPIW)
          call localel(ia(i_g1),ia(i_lel),elLG,noGL,my_numel,maxno)
          call writeelmef(ia(i_lel),my_numel,maxno,nout)
          i_g1  = dealloc('g        ')            
          i_lel = dealloc('lep      ')
        else
          if(nprism6(1) .gt. 0 ) then 
            write(nout,'(a)')"prism6"
            tolel = nprism6(2)+nprism6(1) - 1
            call writeelmef(ia(i_lel),ie,nprism6(2),tolel
     .                     ,maxno,6,type_prism6,nout)
            write(nout,'(a)')"end prism6"
          endif  
c ... elementos em overlaping
          if(nprism6(3) .gt. 0) then
            write(nout,'(a)')"prism6ov"
            tolel = nprism6(3)+nprism6(4) - 1
            call writeelmef(ia(i_lel),ie,nprism6(4)
     .                     ,tolel,maxno,6,type_prism6,nout)
            write(nout,'(a)')"end prism6ov"
c .....................................................................
          endif
        endif  
      goto 100
c .....................................................................
c
c ... initialdisp - deslocamentos inicias:       
 1450 continue
        print*,"(1450) initialdisp deslocamentos inicias."
      goto 100
c .....................................................................
c
c ... initialtemp - temperaturas inicias:       
 1500 continue
        write(nout,'(a)')'initialtemp'
c ... escrevendo tempratura inciais  (temperatura espalhada pelos
c     processos)
        if( npes .gt. 1 )then
          i_xl   = alloc_8('xl       ',ndft,my_nnode)
          i_g1   = alloc_8('Gut0     ',  ndft,nnode)
          call azero(ia(i_g1),ndft*nnode)
c ... recupera os valores gobais das temperaturas
          call global_double(ia(i_g1),ut0,nodedist,nnode,ndft,npes
     .                               ,my_rank,MPIW)   
          call localcoor(ia(i_g1),ia(i_xl),noLG,my_nnode,ndft)
          call writecoormef(ia(i_xl),my_nnode,ndft,nout)
          i_g1  = dealloc('Gut0     ')
          i_xl = dealloc('xl       ')
c ... escrevendo sequencial                  
        else
          i_xl = alloc_8('xl      ',ndft,my_nnode)
          call localcoor(ut0,ia(i_xl),noLG,my_nnode,ndft)
          call writecoormef(ia(i_xl),my_nnode,ndft,nout)
          i_xl = dealloc('xl      ')
        endif
        write(nout,'(a)')'end initialtemp'
      goto 100
c .....................................................................
c
c ...
 1550 continue
        print*,"(1450) velocityfield nao implementado."
      goto 100
c .....................................................................
c
c ... Paralelo:
 1600 continue
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
        ntetra10(1) = ntq1
        ntetra10(2) = ntq2
        nhexa20(1)  = nhq1
        nhexa20(2)  = nhq2
        call wparallel(elLG,noLG,ranks,sizes,fmap,nnof
     .                ,rcvs,dspl,id,idt
     .                ,my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
     .                ,my_nnode,my_numel,my_numel_ov
     .                ,nnodev,nnode,numel,ndf,ndft,nprcs,ovlp,nout)
c .....................................................................
      goto 100
c .....................................................................
c ...
 1650 continue
        print*, "1650"
      goto 100
c .....................................................................
c ...
 1700 continue
        print*, "1700"
      goto 100
c .....................................................................
c
c ... coordenada arquivo binario auxiliar
 1750 continue
        goto 850
      goto 100
c ....................................................................
c
c ...
 1800 continue
        print*, "1800"
      goto 100
c .....................................................................
c
c ...
 1850 continue
        print*, "1850"
      goto 100
c .....................................................................
c
c ... quad4 arquivo binario auxiliar
 1900 continue
        goto 550
      goto 100
c .....................................................................
c
c ... tetra4 arquivo binario auxiliar
 1950 continue
        goto 600
      goto 100
c .....................................................................
c
c ... hexa8 arquivo binario auxiliar
 2000 continue
        goto 650
      goto 100
c .....................................................................
c
c ... 
 2050 continue
         goto 700
      goto 100
c .....................................................................
c
c ... 
 2100 continue
        write(nout,'(a,1x,a)')'hydrostatic',name_hidrostatic(1)
      goto 100
c .....................................................................
c
c ... 
 2150 continue
        write(nout,'(a,1x,a)')'hydrostress',name_hidrostatic(2)
      goto 100
c .....................................................................
c
c ... Restricoes constrainpmec(deslocamentos-pressao)
 2200 continue
        write(nout,'(a)')'constrainpmec'
        if( npes .gt. 1)then
c         i_g1   = alloc_4('Gid      ',  ndf,nnode)
c         call mzero(ia(i_g1),ndf*nnode)
c ... recupera os valores gobais de id          
c         call global_int(ia(i_g1),id,nodedist,nnode,ndf,npes
c    .                   ,my_rank,MPIW)   
c         call wpbound(ia(i_g1),noLG,my_nnode,ndf,nout)   
c         i_g1  = dealloc('Gid      ')
        else  
          call wpbound(id,noLG,my_nnode,ndf,nout)   
        endif
        write(nout,'(a)')'end constrainpmec'
      goto 100
c .....................................................................
c
c ... 
 2250 continue
        write(nout,'(a)')'initialpres'
c ... escrevendo pressao inciais  (pressao espalhada pelos
c     processos)
        if( npes .gt. 1 )then
c         i_xl   = alloc_8('xl       ',   1,my_nnode)
c         i_g1   = alloc_8('Gut0     ',   1,nnode)
c         call azero(ia(i_g1),1*nnode)
c ... recupera os valores gobais das temperaturas
c         call global_double(ia(i_g1),ut0,nodedist,nnode,ndft,npes
c    .                               ,my_rank,MPIW)   
c         call localcoor(ia(i_g1),ia(i_xl),noLG,my_nnode,ndft)
c         call writecoormef(ia(i_xl),my_nnode,ndft,nout)
c         i_g1  = dealloc('Gut0     ')
c         i_xl = dealloc('xl       ')
c ... escrevendo sequencial                  
        else
          malloc_name = 'xl'
          i_xl = alloc_8(malloc_name,1,my_nnode)
          call azero(ia(i_xl),my_nnode)
          call local_no_pmec(el,u0,ia(i_xl),elLG,noGL,my_numel
     .                      ,maxnov,maxno,ndf)
          call write_init_mef(ia(i_xl),my_nnode,1,nout)
          i_xl = dealloc(malloc_name)
        endif
        write(nout,'(a)')'end initialpres'
      goto 100
c .....................................................................
c
c ... 
 2300 continue
        write(nout,'(a)')'initialstress'
c ... escrevendo pressao inciais  (pressao espalhada pelos
c     processos)
        if( npes .gt. 1 )then
c         i_xl   = alloc_8('xl       ',   1,my_nnode)
c         i_g1   = alloc_8('Gut0     ',   1,nnode)
c         call azero(ia(i_g1),1*nnode)
c ... recupera os valores gobais das temperaturas
c         call global_double(ia(i_g1),ut0,nodedist,nnode,ndft,npes
c    .                               ,my_rank,MPIW)   
c         call localcoor(ia(i_g1),ia(i_xl),noLG,my_nnode,ndft)
c         call writecoormef(ia(i_xl),my_nnode,ndft,nout)
c         i_g1  = dealloc('Gut0     ')
c         i_xl = dealloc('xl       ')
c ... escrevendo sequencial                  
        else
          malloc_name = 'xl'
          i_xl = alloc_8(malloc_name,ntn,my_nnode)
          call azero(ia(i_xl),my_nnode)
          call localcoor(tx0,ia(i_xl),noLG,my_nnode,ntn)
          call write_init_mef(ia(i_xl),my_nnode,ntn,nout)
          i_xl = dealloc(malloc_name)
        endif
        write(nout,'(a)')'end initialstress'
      goto 100

      goto 100
c .....................................................................
c
c ... 
 2350 continue
c ... elmtpresloads - cargas nos elementos:
        write(nout,'(a)')'elmtpresloads'
        if( npes .gt. 1)then
          i_g1  = alloc_4('Geload   ',    7,numel)
          call mzero(ia(i_g1),7*numel)
c ... recupera os valores gobais de eload          
          call global_int(ia(i_g1),eload,elmdist,numel,7,npes
     .                   ,my_rank,MPIW)   
          call wpelloads(ia(i_g1),elLG,my_numel,my_numel_ov,nout)
          i_g1 = dealloc('Geload   ')
        else  
          call wpelloads(eloadp,elLG,my_numel,my_numel_ov,nout)
        endif  
        write(nout,'(a)')'end elmtpresloads'
      goto 100
c .....................................................................
      goto 100
c .....................................................................
c
c ... 
 2400 continue
        call writefmate(file_prop,numat,nout)
      goto 100
c .....................................................................
c
c ...
 2450 continue
        write(nout,'(a)')"end mesh"
        i_lel = dealloc('lep      ')
c .....................................................................
c =====================================================================
c
c === macros comandos apos a "end mesh"
      do i = 1, nlines
        write(nout,'(a)')trim(lines(i))
      enddo  
c =====================================================================
c
c ===
      close(nout)
      return
c ===
   15 continue
c ... caso especial de processos ociosos
      nmacros = 0
  105 continue
      nmacros = nmacros + 1 
      write(rc,'(15a)')macros(nmacros)
      do 205 j = 1, nmc
        if(rc .eq. macro(j)) go to 305
  205 continue
      goto 105
  305 continue
      goto( 405, 455, 505
     .    , 555, 605, 655
     .    , 705, 755, 805
     .    , 855, 905, 955
     .    ,1005,1055,1105
     .    ,1155,1205,1255
     .    ,1305,1355,1405
     .    ,1455,1505,1555
     .    ,1605,1655,1705
     .    ,1755,1805,1855
     .    ,1905,1955,2005
     .    ,2055,2105,2155) j
c ...
  405 continue
  705 continue
  755 continue
  805 continue
 1305 continue
 1355 continue
 1405 continue
 1455 continue
 1555 continue
 1605 continue
 1655 continue
 1705 continue
 1805 continue
 1855 continue
 2055 continue
 2105 continue
      goto 105
c .....................................................................
c
c ... conctividades
  455 continue
  505 continue
  555 continue
  605 continue
  655 continue
 1905 continue
 1955 continue
 2005 continue
c ... conectividades no mpi (conectividade espalhada pelos processos)
        i_g1   = alloc_4('Gel     ',maxno+1,numel)
        call global_int(ia(i_g1),el,elmdist,numel,maxno+1,npes,my_rank
     .                 ,MPIW)
        i_g1  = dealloc('Gel     ')
      goto 105    
c .....................................................................
c
c ... Coordenadas
  855 continue
 1755 continue
c ... coordenandas no mpi (coordenadas espalhada pelos processos)
          i_g1   = alloc_8('Gcoorf  ',  ndm,nnode)
c         call azero(ia(i_g1),ndm*nnode)
          call global_double(ia(i_g1),x,nodedist,nnode,ndm,npes,my_rank
     .                      ,MPIW)   
          i_g1  = dealloc('Gcoorf  ')
      goto 105    
c .....................................................................
c
c ... Restricoes constraindisp(deslocamentos)
  905 continue
c ... restricoes no mpi (restricoes espalhada pelos processos)
          i_g1   = alloc_4('Gidf    ',  ndf,nnode)
          call global_int(ia(i_g1),id,nodedist,nnode,ndf,npes,my_rank
     .                    ,MPIW)   
          i_g1  = dealloc('Gidf    ')
      goto 105    
c .....................................................................
c
c ... Resticoes termicas(temperaturas)
  955 continue
c ... restricoes no mpi (restricoes espalhada pelos processos)
          i_g1   = alloc_4('Gidtf   ',  ndft,nnode)
          call global_int(ia(i_g1),idt,nodedist,nnode,ndft,npes,my_rank
     .                    ,MPIW)   
          i_g1  = dealloc('Gidtf   ')
      goto 105    
c .....................................................................
c
c ... nodalforces - forcas nodais
 1005 continue
c ... nodalforces no mpi (nodalforces espalhada pelos processos)
          i_g1  = alloc_4('Gidf    ',  ndf,nnode)
          i_g2  = alloc_8('Gff     ',  ndf,nnode)
          call global_int(ia(i_g1),id,nodedist,nnode,ndf,npes
     .                   ,my_rank,MPIW)   
          call global_double(ia(i_g2),f,nodedist,nnode,ndf,npes
     .                   ,my_rank,MPIW)   
          i_g2  = dealloc('Gidf    ')
          i_g1  = dealloc('Gff     ')
      goto 105    
c .....................................................................
c
c ... nodalsources - fontes termicas
 1055 continue
c ... nodalsources no mpi (nodalsources espalhada pelos processos)
          i_g1  = alloc_4('Gidtf   ',  ndft,nnode)
          i_g2  = alloc_8('Gftf    ',  ndft,nnode)
          call global_int(ia(i_g1),idt,nodedist,nnode,ndft,npes
     .                   ,my_rank,MPIW)   
          call global_double(ia(i_g2),ft,nodedist,nnode,ndft,npes
     .                   ,my_rank,MPIW)   
          i_g2  = dealloc('Gidtf   ')
          i_g1  = dealloc('Gftf    ')
      goto 105    
c .....................................................................
c
c ...nodalloads - nos com cargas variaveis no tempo:
 1105 continue
c ... nodalloads no mpi (nodalloads espalhada pelos processos)
          i_g1  = alloc_4('Gidf    ',  ndf,nnode)
          i_g2  = alloc_4('Gnloadf ',  ndf,nnode)
          call global_int(ia(i_g1),id,nodedist,nnode,ndf,npes
     .                   ,my_rank,MPIW)   
          call global_int(ia(i_g2),nload,nodedist,nnode,ndf,npes
     .                   ,my_rank,MPIW)   
          i_g2  = dealloc('Gidf    ')
          i_g1  = dealloc('Gnloadf ')
      goto 105    
c .....................................................................
c
c ...nodalthermloads - nos com cargas variaveis no tempo:
 1155 continue
c ... nodalthermloads no mpi (nodalthermloads espalhada pelos processos)
          i_g1  = alloc_4('Gidtf   ', ndft,nnode)
          i_g2  = alloc_4('Gnloadtf', ndft,nnode)
          call global_int(ia(i_g1),idt,nodedist,nnode,ndft,npes
     .                   ,my_rank,MPIW)   
          call global_int(ia(i_g2),nloadt,nodedist,nnode,ndft,npes
     .                   ,my_rank,MPIW)   
          i_g2  = dealloc('Gidtf   ')
          i_g1  = dealloc('Gnloadtf')
      goto 105    
c .....................................................................
c
c ...eload - elementos com cargas mecanicas:
 1205 continue
c ... eload no mpi (eload espalhada pelos processos)
          i_g1  = alloc_4('Geloadf ',    7,numel)
          call global_int(ia(i_g1),eload,elmdist,numel,7,npes
     .                   ,my_rank,MPIW)   
          i_g1  = dealloc('Geloadf ')
      goto 105    
c .....................................................................
c
c ...eloadt- elementos com cargas termicas:
 1255 continue
c ... eloadt no mpi (eloadt espalhada pelos processos)
          i_g1  = alloc_4('Geloadtf',    7,numel)
          call global_int(ia(i_g1),eloadt,elmdist,numel,7,npes
     .                   ,my_rank,MPIW)   
          i_g1  = dealloc('Geloadtf')
      goto 105    
c .....................................................................
c
c ... initialtemp - temperaturas inicias:       
 1505 continue
c ... intialtemp   no mpi (inttialtemp espalhada pelos processos)
          i_g1   = alloc_8('Gut0f   ',  ndft,nnode)
          call global_double(ia(i_g1),ut0,nodedist,nnode,ndft,npes
     .                      ,my_rank,MPIW)   
          i_g1  = dealloc('Gut0f   ')
      goto 105    
c .....................................................................
c =====================================================================
c ... end mesh   
 2155 continue
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * WRITECOORMEF: escreve as coordenadas no formato do mefpar         *
c *********************************************************************
      subroutine writecoormef(x,nnode,ndm,nout)
      implicit none
c === variaveis
c ... malha
      integer nnode,ndm
      real*8 x(ndm,*)
c *********************************************************************
c ... arquivo
      integer nout
c ... auxiliar
      integer i,j
c =====================================================================
c
c === 
      do i = 1, nnode
        if( ndm .eq. 1)then
          write(nout,'(i8,1x,es20.10)')i,(x(j,i),j=1,ndm)
        else if( ndm .eq. 2)then
          write(nout,'(i8,1x,2es20.10)')i,(x(j,i),j=1,ndm)
        else if(ndm .eq. 3)then 
          write(nout,'(i8,1x,3es20.10)')i,(x(j,i),j=1,ndm)
        endif  
      enddo
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * WRITE_INIT_MEF: escreve os valores inicias no formato do mefpar   *
c *********************************************************************
      subroutine write_init_mef(u,nnode,ndf,nout)
      implicit none
c === variaveis
c ... malha
      integer nnode,ndf
      real*8 u(ndf,*),add
c *********************************************************************
c ... arquivo
      integer nout
c ... auxiliar
      integer i,j
c =====================================================================
c
c === 
      do i = 1, nnode
c ... verifica se todas as forcas sao nulas            
        add = 0.d0
        do j = 1, ndf
          add = add + dabs(u(j,i))
        enddo
c .....................................................................
c
c ...
        if( add .ne. 0.d0) then
          write(nout,'(i8)',advance='no') i
          do j = 1, ndf
            write(nout,'(1x,es20.10)',advance='no') u(j,i)
          enddo
          write(nout,'(a1)')' '
        endif
c .....................................................................
      enddo
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * WRITEELMEF : escreve os elementos no formato do mefpar            *
c *********************************************************************
      subroutine writeelmef(el,ie,it,nit,maxno,nen,typeElm,nout)
c ===
      implicit none
      include 'elementos.fi'
c ... malha
      integer maxno,typeElm,mat,nen
      integer el(maxno+1,*),ie(*)
c ... particionamento
      integer it,nit
c ... arquivo
      integer nout
c ... auxiliar
      integer i,j
c =====================================================================
c
c ===
      do i = it, nit
        mat = el(maxno+1,i)
        if( geom_el_type(ie(mat)) .eq. typeElm) then
          write(nout,'(22i10)' )i,(el(j,i),j=1,nen),mat
        endif
      enddo
c =====================================================================
c
c ===
      return
      end
c *********************************************************************
c
c *********************************************************************
c * WRITETERMPROP:                                                    *
c *********************************************************************
      subroutine writetermprop(numat,nout)
c ===
      implicit none
      include 'termprop.fi'
c ... auxiliar
      integer i,j,k 
      integer pol
c ... arquivo
      integer nout
c ... malha      
      integer numat
c =====================================================================
c
c ===
      do i = 1, numat
        do j = 1 , maxmate
          pol =eprop(j,i)
c ... verifica se ha propriedade variavel no tempo       
          if(pol .gt. 0 .and. j .ne. 9) then
            write(nout,'(3i8)')i,j,pol
            do k = 1, pol
              write(nout,'(2es18.9)')nprop(k,j,i),nprop(k+pol,j,i)
            enddo
          endif
        enddo
      enddo
c =====================================================================
c
c ===
      return
      end
c *********************************************************************
c
c *********************************************************************
c * WRITEADIABT:                                                      *
c *********************************************************************
      subroutine writeadiabat(numat,nout)
c ===
      implicit none
      include 'termprop.fi'
c ... auxiliar
      integer i,j,k 
      integer pol
c ... arquivo
      integer nout
c ... malha      
      integer numat
c =====================================================================
c
c ===
      do i = 1, numat
          pol =eprop(9,i)
c ... verifica se ha propriedade variavel no tempo       
          if(pol .gt. 0) then
            write(nout,'(3i8)')i,pol
            do k = 1, pol
              write(nout,'(2es18.9)')tad(k,i),tempad(k,i)
            enddo
          endif
      enddo
c =====================================================================
c
c ===
      return
      end
c *********************************************************************
c
c *********************************************************************
c * WRITEPARAMETERS :                                                 *
c *********************************************************************
      subroutine writeparameters(nnode,numel,numat,maxno,ndf,ndft
     .                          ,ndm,npi,nout)
      implicit none
      include 'readfile.fi'
c === variaveis
c ... malha
      integer nnode,numel,numat,maxno,ndf,ndft,ndm,npi
c ... arquivo
      integer nout
c =====================================================================
c
c ===
        write(nout,'(a4)')'mesh'
        write(nout,'(2(a5,1x,i9,1x),2(a5,1x,i2,1x))',advance='no')
     .     'nnode',nnode,'numel',numel,'numat',numat,'maxno',maxno
c ...
        if(read_parameter(5)) then 
          write(nout,'(a8,1x,i2,1x)',advance='no')'ndf',ndf
        endif
        if(read_parameter(6)) then 
          write(nout,'(a8,1x,i2,1x)',advance='no')'therm',ndft
        endif
        if(read_parameter(7)) then 
          write(nout,'(a8,1x,i2,1x)',advance='no')'dim',ndm
        endif
        if(read_parameter(8)) then 
          write(nout,'(a8,1x,i2,1x)',advance='no')'npi',npi
        endif
        write(nout,*)
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
c * WRITEMATE :                                                       *
c *********************************************************************
      subroutine writemate(ie,e,numat,nout)
      implicit none
      include 'termprop.fi'
c === variaveis
c ... malha
      integer numat
      integer ie(*)
      real*8 e(maxprop,*)
c ... arquivo
      integer nout
c ... auxiliar
      integer i,j
c =====================================================================
c
c === 
      write(nout,'(a)')'materials'
      do i = 1, numat
        write(nout,'(2(i2,1x),30(1x,es18.6))')i,ie(i)
     .       ,(e(j,i),j=1,read_prop(i))
      enddo
      write(nout,'(a)')'end materials'
c ... Pre antigo(yuri)
c     write(nout,'(a)')'materials'
c     do i = 1, numat
c       write(nout,'(2(i2,1x),4es18.9)')
c    .              i,ie(i),(e(j,i),j=1,4)
c     enddo
c     write(nout,'(a)')'end materials'
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * Data de criacao    : 21/03/2017                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * WRITEFMATE: escreve a macro de fmaterias                           *
c * ------------------------------------------------------------------ *
c * Par�metros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * fprop(*)   - nome dos arquivos com as propriedades                 *
c * numat      - numero de materias                                    *
c * nout       - arquivo de escrita                                    *
c * ------------------------------------------------------------------ *
c * Par�metros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c *********************************************************************
      subroutine writefmate(fprop,numat,nout)
      implicit none
c === variaveis
c ... malha
      integer numat
      character*80 fprop(*)
c ... arquivo
      integer nout
c ... auxiliar
      integer i
c =====================================================================
c
c === 
      write(nout,'(a)')'fmaterials'
      do i = 1, numat
        write(nout,'(i2,1x,a80)')i,fprop(i)
      enddo
      write(nout,'(a)')'end fmaterials'
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * Restricoes nodais                                                 *
c *********************************************************************
      subroutine wpbound(id,noLG,my_nnode,ndf,nout)
c === variaveis
      implicit none
c ... malha
      integer my_nnode,ndf
      integer id(ndf,*)
c ... particionamento
      integer noLG(*)
c ... arquivo
      integer nout
c ... auxiliares
      integer j,k,noL,noG
c =====================================================================
c
c ===
c ... loop nos nos locais
      do noL = 1, my_nnode
        noG  = noLG(noL)
c ... Existe restricao neste no?        
        k = 0
        do j =1 ,ndf
          k = k + id(j,noG)
        enddo
c ... impressao
        if(k .ne. 0) write(nout,'(i8,1x,6i2)')nol,(id(k,noG),k=1,ndf)
      enddo
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c 
c *********************************************************************
c * Forcas nodais e valores pescritos.                                *
c *********************************************************************
      subroutine wpforces(f,id,noLG,my_nnode,nno1,nno2,ndf,nout)
c === variaveis
      implicit none
c ... malha
      integer ndf
      integer id(ndf,*)
      real*8 f(ndf,*),add,aux(6)
c ... particionamento
      integer noLG(*),nno1,nno2,my_nnode
c ... arquivo
      integer nout
c ... auxiliares
      integer j,noL,noG,nnoL
c =====================================================================
c
c ===
c ... loop nos nos locais
      nnoL = nno1 + nno2
      do noL = 1, my_nnode
        noG  = noLG(noL)
c ... Existe restricao neste no?        
        add = 0.d0
        do j =1 ,ndf
          add = add + dabs(f(j,noG))
          if(id(j,noG) .ne. 0) add = add + 1.d0
        enddo
        if( add .ne. 0.d0 )then
          aux(1:ndf) = f(1:ndf,noG)
c ... o no for V3 ou V4 a contribuicao da forcas 
c     sera zero nesta particao caso seja deslocamento nada se altera
          if(nol .gt. nnoL)then
            do j = 1, ndf
              if(id(j,noG) .eq. 0) aux(j) =0.d0
            enddo
          endif
c ... verifica se todas as forcas sao nulas            
          add = 0.d0
          do j = 1, ndf
            add = add + dabs(aux(j))
          enddo
c ... so imprime se ao menos houver uma forca nao nula          
          if( add .ne. 0.d0) then 
            write(nout,'(i8,1x,6es18.9)')noL,(aux(j),j=1,ndf)
          endif  
        endif
      enddo  
c =====================================================================
c
c ===
      return
      end
c ===================================================================== 
c *********************************************************************
c
c *********************************************************************
c * WPLOADS : escrita das cargas por no                               *
c *********************************************************************
      subroutine wploads(nloads,id,noLG,my_nnode,nno1,nno2,ndf,nout)
c === variaveis      
      implicit none
c ... malha
      integer id(ndf,*),nloads(ndf,*),ndf
c ... particionamento
      integer nno1,nno2,my_nnode
      integer noLG(*)
c ... arquivo
      integer nout 
c ... auxiliares
      integer nnoL,noL,j,noG
      integer add,aux(6)
c =====================================================================
c
c ===
c ... loop nos nos locais
      nnoL = nno1 + nno2
      do noL = 1, my_nnode
        noG  = noLG(noL)
c ... Existe restricao neste no?        
        add = 0
        do j =1 ,ndf
          add = add + nloads(j,noG)
        enddo
        if( add .ne. 0 )then
          aux(1:ndf) = nloads(1:ndf,noG)
c ... se o no for V3 ou V4 a contribuicao da forcas 
c     sera zerada nesta particao, caso seja deslocamento, nada se altera   
          if(nol .gt. nnoL)then
            do j = 1, ndf
              if(id(j,noG) .eq. 0) aux(j) = 0
            enddo
          endif
c ... verifica se todos os loads desse no sao nulos            
          add = 0
          do j = 1, ndf
            add = add + aux(j)
          enddo
c ... so imprime se ao menos houver um load nesse no nao nulo          
          if( add .ne. 0) then 
            write(nout,'(i8,1x,9i9)')noL,(aux(j),j=1,ndf)
          endif  
        endif
      enddo  
c =====================================================================
c
c ===
      return
      end
c ===================================================================== 
c **********************************************************************
c
c *********************************************************************
c * WPELLOADS : escrita das cargas por elemento                       *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * eloads     -> carga por elemento                                  *
c * elLG       -> mapa local-global de elementos                      *
c * my_numel   -> numero de elementos da particao                     *
c * my_numel_ov-> numero de elementos da particao em overllaping      *
c * ndf        -> grau de liberdade                                   *
c * nout       -> arquivo de saida                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine wpelloads(eloads,elLG,my_numel,my_numel_ov,nout)
c === variaveis      
      implicit none
c ... malha
      integer eloads(7,*)
c ... particionamento
      integer my_numel,my_numel_ov
      integer elLG(*)
c ... arquivo
      integer nout 
c ... auxiliares
      integer numel_nov,nelL,nelG,j
      integer add,aux(7)
c =====================================================================
c
c ===
c ... loop nos elementos locais over e non-overllaping
c     numel_nov = my_numel - my_numel_ov
      do nelL = 1, my_numel
        nelG  = elLG(nelL)
        aux(1:7) = eloads(1:7,nelG)
c ... verifica se ha loads nesse elementos            
        add = 0
        do j = 1, 7
          add = add + aux(j)
        enddo
c ... so imprime se ao menos houver um load nesse elemento nao nulo          
        if( add .ne. 0) then 
          write(nout,'(i8,1x,7i9)')nelL,(aux(j),j=1,7)
        endif  
      enddo  
c =====================================================================
c
c ===
      return
      end
c ===================================================================== 
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   Escrita das cargas variaveis no tempo                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    idl(6) - variavel auxuliar                                      *
c *    nin - arquivo de entrada                                        *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    load(1,n) - tipo da carga n (n <= numload)                      *
c *    load(2,n) - numero de termos da carga n                         *
c *    fload(i,j,k) - coeficiente do termo j da carga k                *
c *                                                                    *
c **********************************************************************
      subroutine wload(nout)
c ...  
      implicit none
      include 'string.fi'
      include 'load.fi'
      integer  nout
      integer  i,j,k,itype,nparc
c .....................................................................
C
c ...
         do i = 1, numload
            itype = load(1,i)
            if (itype .gt. 0) then
               nparc = load(2,i)
               if     (itype .eq. 1) then
c ...             valor da carga:            
                  write(nout,10) i,itype,fload(1,1,i)
c .....................................................................
c
c ...
               elseif (itype .eq. 2) then
c ...             u(t) = a.t + b 
c ...             valores  b,a:
                  write(nout,10) i,itype,(fload(j,1,i),j=1,2)
c .....................................................................
c
c ...
               elseif (itype .eq. 3) then
c ...             -kdu/dx = H.(u-uext) 
c ...             valores  uext,H:
                  write(nout,10) i,itype,(fload(j,1,i),j=1,2)
c .....................................................................
c
c ...
               elseif (itype .eq. 4) then
c ...             u(t) = Sum(k)[ak.sin(bk.t + ck)]
c ...             valores  ak,bk,ck:
                  write(nout,10) i,itype,
     .                               ((fload(j,k,i),j=1,3),k=1,nparc)
c .....................................................................
c
c ...
               elseif (itype .eq. 5) then
c ...             u(t) = Sum(k)[ak.cos(bk.t + ck)]
c ...             valores  ak,bk,ck:
                  write(nout,10) i,itype,
     .                               ((fload(j,k,i),j=1,3),k=1,nparc)
c .....................................................................
c
c ...
               elseif (itype .eq. 6) then
                  write(nout,11) i,itype,nparc
                  do k = 1, nparc
                     write(nout,12) (fload(j,k,i),j=1,3)
                  enddo
c .....................................................................
c
c ...
               elseif (itype .eq. 7) then
c ...       kdu/dx = emiss * const(Stef-Boltz) *(u4-uext4)+H(uext-u)
                  write(nout,14) i,itype,(fload(j,1,i),j=1,3),load(2,i)
c ...      valores (tempo,temperatura)            
                  do k = 1, nparc
                     write(nout,15) (fload(k,j,i),j=2,3)
                  enddo
c .....................................................................
c
c ...
               elseif (itype .eq. 8) then
c ...       carga constante a partir de um dado t0
                  write(nout,11) i,itype,load(2,i)
                  do k = 1, nparc
                     write(nout,15) (fload(k,j,i),j=1,2)
                  enddo
c .....................................................................
c
c ...
               elseif (itype .eq. 9) then
c ...       interpolacao de pontos
                  write(nout,11) i,itype,load(2,i)
                  do k = 1, nparc
                     write(nout,15) (fload(k,j,i),j=1,2)
                  enddo
c .....................................................................
c
c ...    troca de calor de Newton(com uExt e h variavel):    
               elseif (itype .eq. 10) then
c ...       interpolacao de pontos
                  write(nout,11) i,itype,load(2,i)
                  do k = 1, nparc
                     write(nout,16)
     .               ,fload(k,1,i),fload(k,2,i),fload(k,3,i)
                  enddo
c .....................................................................
c
c ...
               elseif (itype .eq. 17) then
c ...             carga do problema geofisico                    
                  write(nout,13) i,itype,nparc,(fload(j,1,i),j=1,nparc)
c .....................................................................
c
c ...
               elseif (itype .eq. 18) then
c ...             variacao senoidal                  
                  write(nout,13) i,itype,nparc,fload(1,1,i)
     .                          ,fload(2,1,i),fload(3,1,i),fload(1,2,i)
c .....................................................................
c                 
c ... 
c 39 - valor de uma carga ou deslocamento variavel no tempo
c 40 - forca distribuida constante no contorno
c 41 - carga normal constante no contorno 
c      (F=carga*normal,normal - calculado nivel elemento)
c 42 - valor por carga hidrostatica      
c      (f = alfa*(f0 + density*g*(h-x)) )
c 43 - fluxo normal de massa
          elseif (   (itype .eq. 39)    
     .          .or. (itype .eq. 40) .or. (itype .eq. 41) 
     .          .or. (itype .eq. 42) .or. (itype .eq. 43) ) then
c ...       interpolacao de pontos
                  write(nout,11,advance = "no") i,itype
                  write(nout,'(80a)')name_loads(i)
                endif
            endif
         enddo
c ....................................................................
c
c ...
      return
   10 format (2(i4,1x),99es18.9)
   11 format (3(i4,1x))
   12 format (3(es18.9,1x))
   13 format (3(i4,1x),99es18.9)
   14 format (2(i4,1x),3(es18.9,1x),i9)
   15 format (2(es18.9,1x))
   16 format (3(es18.9,1x))
   17 format (99(es18.9,1x))
      end                        
c ********************************************************************* 
c
c ********************************************************************* 
      subroutine wparallel(elLG,noLG,ranks,sizes,fmap,nnof,rvcs,dspl
     .                    ,id,idt
     .                    ,my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
     .                    ,my_nnode,my_numel,my_numel_ov
     .                    ,nnodev,nnode,numel,ndf,ndft,nprcs,ovlp,nout)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 09/10/2016                                    *
c * ------------------------------------------------------------------ * 
c * WPARALLEL : Escreve informacoes da paralelizacao geradas pelo      *
c * pre-mpi atraves do macro-comando parallel em wdat                  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * Obs:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'elementos.fi'
      integer nprcs
      integer my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
      integer my_nnode,my_numel,my_numel_ov
      integer elLG(*),noLG(*)
      integer ranks(*),sizes(*)
      integer fmap(*),nnof(*)
      integer rvcs(*),dspl(*)
      integer id(ndf,*),idt(ndft,*)
      integer nout,nnodev,nnode,numel,ndf,ndft
      logical ovlp
c ... auxiliar      
      integer i
c-----------------------------------------------------------------------
       write(nout,'(a)')  'parallel'
       write(nout,'(a4,1x,i9)') 'novG',nnodev 
       write(nout,'(a4,1x,i9)') 'noG' ,nnode 
       write(nout,'(a4,1x,i9)') 'nelG',numel               
c ......................................................................
c
c ..... Metodo de sub-divisao de dominio
c
        if(ovlp)then
          write(nout,'(a)') 'overlapping'
        else
          write(nout,'(a)') 'non-overlapping'
        endif
c ......................................................................
c
c ..... Numeros de no's pelo tipo de nó
c
        write(nout,'(5(a4,1x,i9,1x))')
     .     'nno1',my_nno1
     .    ,'nno2',my_nno2
     .    ,'nno3',my_nno3
     .    ,'nno4',my_nno4
     .    ,'no1a',my_nno1a
        if (ovlp)write(nout,'(a,1x,i9)') 'numel_ov',my_numel_ov
c ......................................................................
c
c ..... Mapa Local -> Global de no's
c
        write(nout,'(a)') 'noLG'
        call Wrtvet_4(noLG,my_nnode,nout)
c       do i = 1, my_nnode
c         write(nout,'(i10)')noLG(i) 
c       enddo  
c ......................................................................
c
c ..... Numeros de elemento's pelo tipo
c
        write(nout,'(8(a8,1x,2(i9,1x),1x))')
     .     'nbar2  ' , nbar2(1)   , nbar2(2)
     .    ,'ntria3 ' , ntria3(1)  , ntria3(2)
     .    ,'nquad4 ' , nquad4(1)  , nquad4(2)
     .    ,'ntetra4' , ntetra4(1) , ntetra4(2)
     .    ,'nhexa8 ' , nhexa8(1)  , nhexa8(2)
     .    ,'nprism6' , nprism6(1) , nprism6(2)
     .    ,'ntetra10', ntetra10(1), ntetra10(2)
     .    ,'nhexa20' , nhexa20(1) , nhexa20(2)
c ......................................................................
c
c ..... Mapa Local -> Global de elementos
c
        write(nout,'(a4,1x,i9)') 'elLG'        
        call Wrtvet_4(elLG,my_numel,nout)
c        do i = 1, my_numel
c          write(nout,'(i10)')elLG(i) 
c        enddo  
c ......................................................................
c
c ..... Mapa Interface -> Global de no's
       write(nout,'(a4,1x,2i9)') 'fmap',nnof(1),nnof(2)
       call Wrtvet_4(fmap,nnof(1),nout)
       if(ovlp)call Wrtvet_4(fmap(my_nnode+1),nnof(2),nout)
c      do i = 1, nnof(1) 
c        write(nout,'(i10)') fmap(i)
c      enddo  
c ......................................................................
c
c ..... Restricoes Globais
c
c       if ( ndf .gt. 0) then
c         write(nout,'(a)') 'idG'
c         call Wrtmtx_4(id,ndf,nnode,nout)  
c         write(nout,'(a)')  'end idG'
c       endif  
c       if ( ndft .gt. 0) then
c         write(nout,'(a)') 'idGt'
c         call Wrtmtx_4(idt,ndft,nnode,nout)  
c         write(nout,'(a)')  'end idGt'
c       endif  
c ......................................................................
c
c ..... Estruturas de comunicacao send/recv
c
      write(nout,'(a)') 'nviz'
      write(nout,'(128(i9,1x))')sizes(1) 
      if(ovlp)write(nout,'(128(i9,1x))') sizes(2)
c .......
      write(nout,'(a)') 'ranks'
      write(nout,'(129(i9,1x))') (ranks(i),i=1,sizes(1))
      if(ovlp) then
        write(nout,'(129(i9,1x))')(ranks(i),i=nprcs+1,nprcs+sizes(2))
      endif  
c .......
      write(nout,'(a)') 'rcvs'
      write(nout,'(128(i9,1x))')(rvcs(i),i=1,sizes(1)) 
      if(ovlp)then
        write(nout,'(257(i9,1x))')
     .       (rvcs(i),i=nprcs+1,nprcs+sizes(2)) 
      endif  
c .......
      write(nout,'(a)') 'dspl'
      write(nout,'(257(i9,1x))')(dspl(i),i=1,sizes(1))
      if(ovlp)then
        write(nout,'(257(i9,1x))')
     .       (dspl(i),i=nprcs+1,nprcs+sizes(2)) 
      endif
c antigo:
c         write(nout,'(a)') 'nviz'
c         write(nout,'(128(i3,1x))') ia(i_sizes-1+pid)
c         write(nout,'(a)') 'ranks'
c         write(nout,'(129(i3,1x))') (ia(i_ranks+j+(pid-1)*nprcs),
c     .   j=0,ia(i_sizes-1+pid)-1)
c.......................................................................
      write(nout,'(a)')  'end parallel'
c.......................................................................
c ... loop nos processos
c     enddo
c ......................................................................
      return
      end
c =====================================================================
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine Wrtvet_4(array,lin,nin)
c **********************************************************************
c *                                                                    *
c *   Wrtvet_4                                                         *
c *   --------                                                         *
c *                                                                    *
c *   Escreve arranjo 1D inteiro                                       *
c *                                                                    *
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *   Parâmetros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer array(*),lin,nin,i
c ......................................................................
      do i=1,lin
         write(nin,'(i9)')array(i)
c       write(nin,'(i,i9)')i,array(i)
      enddo
c ......................................................................
      return
      end
c **********************************************************************
c
      subroutine Wrtmtx_4(array,lin,col,nin)
c **********************************************************************
c *                                                                    *
c *   Wrtmtx_4                                                         *
c *   --------                                                         *
c *                                                                    *
c *   escreve arranjo 2D inteiro                                       *
c *                                                                    *
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *   Parâmetros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer array(lin,*),lin,col,nin,i,j,k
c ......................................................................
      do j=1,col
        k = 0
        do i = 1, lin
          k = k + array(i,j)
        enddo
        if(k.ne.0)write(nin,'(i7,1x,6i9)') j,(array(i,j),i=1,lin)
      enddo
c ......................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine controlepre(noLG,elLG,ranks,sizes,rcvs,dspl,fmap,nnof
     .                    ,my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
     .                    ,my_nnode,my_numel_nov,my_numel_ov
     .                    ,nnode,numel,ovlp,my_rank,nprcs,filein,nout)
c ====
      implicit none
      include 'elementos.fi'
      integer noLG(*),elLG(*),ranks(*),sizes(*)
      integer rcvs(*),dspl(*),fmap(*),nnof(*)
      integer my_nno1,my_nno2,my_nno3,my_nno4,my_nno1a
      integer my_nnode,my_numel,my_numel_nov,my_numel_ov
      integer nnode,numel
      integer my_rank,nprcs
      integer i
      logical ovlp
      character filein*80,fileout*100,name*80
      integer nout
c =====================================================================
c
c ===
      fileout = name(filein,my_rank,108)
      open(nout,file=fileout)
c ...
       write(nout,'(a)')  'parallel'
       write(nout,'(a4,1x,i9)') 'noGL',nnode 
       write(nout,'(a4,1x,i9)') 'nelG',numel               
c ......................................................................
c
c ..... Metodo de sub-divisao de dominio
c
        if(ovlp)then
          write(nout,'(a)') 'overlapping'
        else
          write(nout,'(a)') 'non-overlapping'
        endif
c ......................................................................
c
c ..... Numeros de no's pelo tipo de nó
c
        write(nout,'(5(a4,1x,i9,1x))')
     .     'nno1',my_nno1
     .    ,'nno2',my_nno2
     .    ,'nno3',my_nno3
     .    ,'nno4',my_nno4
     .    ,'no1a',my_nno1a
        if (ovlp)write(nout,'(a,1x,i9)') 'numel_ov',my_numel_ov
c ......................................................................
c
c ..... Mapa Local -> Global de no's
c
        write(nout,'(a)') 'noLG'
        call Wrtvet_4(noLG,my_nnode,nout)
c ......................................................................
c
c       write(nout,'(5(a7,1x,2(i9,1x),1x))')
c    .     'nbar2  ', nbar2(1)  , nbar2(2)
c    .    ,'ntria3 ', ntria3(1) , ntria3(2)
c    .    ,'nquad4 ', nquad4(1) , nquad4(2)
c    .    ,'ntetra4', ntetra4(1), ntetra4(2)
c    .    ,'nhexa8 ', nhexa8(1) , nhexa8(2)
c ......................................................................
c
c ..... Mapa Local -> Global de elementos
c
        my_numel = my_numel_ov + my_numel_nov
        write(nout,'(a4,1x,i9)') 'elLG'        
        call Wrtvet_4(elLG,my_numel,nout)
c ......................................................................
c
       write(nout,'(a4,1x,2i9)') 'fmap',nnof(1),nnof(2)
       call Wrtvet_4(fmap,nnof(1),nout)
       if(ovlp)call Wrtvet_4(fmap(nnode+1),nnof(2),nout)
c
c ..... Estruturas de comunicacao send/recv
c
      write(nout,'(a)') 'nviz'
      write(nout,'(128(i9,1x))')sizes(1) 
      if(ovlp)write(nout,'(128(i9,1x))') sizes(2)
c .......
      write(nout,'(a)') 'ranks'
      write(nout,'(129(i9,1x))') (ranks(i),i=1,sizes(1))
      if(ovlp) then
        write(nout,'(129(i9,1x))')(ranks(i),i=nprcs+1,nprcs+sizes(2))
      endif  
c .......
      write(nout,'(a)') 'rcvs'
      write(nout,'(128(i9,1x))')(rcvs(i),i=1,sizes(1)) 
      if(ovlp)then
        write(nout,'(257(i9,1x))')
     .       (rcvs(i),i=nprcs+1,nprcs+sizes(2)) 
      endif  
c .......
      write(nout,'(a)') 'dspl'
      write(nout,'(257(i9,1x))')(dspl(i),i=1,sizes(1))
      if(ovlp)then
        write(nout,'(257(i9,1x))')
     .       (dspl(i),i=nprcs+1,nprcs+sizes(2)) 
      endif
c ......................................................................
      write(nout,'(a)')  'end parallel'
      write(nout,'(a)')  'stop'
c =====================================================================
c
c === 
      return
      end
c =====================================================================
