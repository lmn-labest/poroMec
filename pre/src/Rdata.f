      subroutine rdatm(nnodev,nnode,numel,numat,maxnov,maxno
     1                ,ndf  ,ndft,ndm,npi,nst  ,nstt
     2                ,i_ix ,i_ie,i_e,i_x
     3                ,i_id ,i_nload,i_eload 
     4                ,i_f  ,i_u,i_v,i_vt     ,i_a ,i_tx0
     5                ,i_idt,i_nloadt,i_eloadt
     6                ,i_ft ,i_ut    ,i_ut0   ,i_du,i_w
     7                ,nin,verbose,mrload     ,file_prop,ncont)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 03/05/2017                                    * 
c * ------------------------------------------------------------------ * 
c * RDAT: leitura de dados.                                            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * nin     - arquivo de entrada                                       *  
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * nnodev  - numero de nos dos vertices                               * 
c * nnode   - numero de nos                                            *
c * numel   - numero de elementos                                      *
c * numat   - numero de materiais                                      *
c * maxnov  - numero maximo de nos de vertices por elem nto            *
c * maxno   - numero maximo de nos por elementos                       *
c * ndf     - numero max. de graus de liberdade por no                 *
c * ndm     - dimensao (1, 2 ou 3)                                     *
c * npi     - numero de pontos de integracao                           *
c * plastic - paramento plastico                                       *
c * nst     - numero de graus de liberdade por elemento                *
c * nstt    - numero de graus de liberdade por elemento                *
c * i_ix    - ponteiro para conetividades                              *
c * i_id    - ponteiro para restricoes nodais (mecanico)               *
c * i_ie    - ponteiro para materiais                                  *
c * i_nload - ponteiro para o arranjo nload (mecanico)                 *
c * i_eload - ponteiro para o arranjo eload (mecanico)                 *
c * i_e     - ponteiro para o arranjo e                                *
c * i_x     - ponteiro para o arranjo x                                *
c * i_f     - ponteiro para o arranjo f (mecanico)                     *
c * i_u     - ponteiro para o arranjo u (mecanico)                     *
c * i_tx0   - ponteiro para o arranjo tx0 (mecanico-poromec)           *
c * i_v     - ponteiro para o arranjo v (mecanico)                     *
c * i_a     - ponteiro para o arranjo a (mecanico)                     *
c * i_idt    - ponteiro para restricoes nodais  (termico)              *
c * i_nloadt - ponteiro para o arranjo nloadt   (termico)              *
c * i_eloadt - ponteiro para o arranjo eloadt   (termico)              *
c * i_ft     - ponteiro para o arranjo ft       (termico)              *
c * i_ut     - ponteiro para o arranjo ut       (termico)              *
c * i_vt     - ponteiro para o arranjo vt       (termico)              *
c * i_w     - ponteiro para o arranjo w         (termico)              *
c * verbose - modo com saida para o prompt                             *
c * file_prop - nome dos arquivos de propriedades                      *
c * ------------------------------------------------------------------ * 
c * Obs:                                                               *
c * ------------------------------------------------------------------ * 
c * ix    - conetividades nodais dos elementos                         *
c * id    - restricoes nodais                                          *
c * ie    - tipo de elemento                                           *
c * e     - constantes fisicas                                         *
c * x     - coordenadas nodais                                         *
c * f     - forcas e prescricoes nodais                                *
c * nload(i,j) - numero identificador da carga na direcao i do no j    *
c * load(1,n)  - tipo da carga n                                       *
c * load(2,n)  - numero de termos da carga n                           *
c * fload(i,j,k) - coeficiente i do termo j da carga k                 *
c **********************************************************************
      use Malloc
      implicit none
c ......................................................................
      include 'elementos.fi'
      include 'load.fi'
      include 'string.fi'
      include 'parallel.fi'
      include 'termprop.fi'
c ......................................................................
      integer nnodev,nnode,numel,numat,maxnov,maxno,maxgrade,dum
      integer ndf,ndft,ntn,ndm,npi,nst,nstt
c ... saida para o prompt
      logical verbose
c ... ponteiros      
      integer*8 i_e,i_x,i_f,i_nload,i_eload,i_u,i_v,i_a,i_tx0
      integer*8 i_idt,i_nloadt,i_eloadt,i_ft,i_ut,i_vt,i_w,i_ut0,i_du
      integer*8 i_ix,i_id,i_ie,i_geomel,i_nincid,i_incid
c ... leitura das macros
      integer ncont 
      logical mrload(ncont)
c ...
      character*80 file_prop(*)
c ... elementos quadraticos
      logical f_read_el,el_quad
c ...  
      logical init_poro_pmec,fmec,fpmec      
c ......................................................................
      integer nin
      integer i,j,nmc,totnel,nmacros
      character*15 macro(42),rc
      character*80 fname
      character*30 string
      integer naux,nincl /7/
c ......................................................................
      data macro/'materials      ','bar2           ','tria3          ',
     1           'quad4          ','tetra4         ','hexa8          ',
     2           'quad8          ','termprop       ','adiabat        ',
     3           'coordinates    ','constraindisp  ','constraintemp  ',
     4           'nodalforces    ','nodalsources   ','nodalloads     ',
     5           'nodalthermloads','elmtloads      ','elmtthermloads ',
     6           'loads          ','hidrprop       ','prism6         ',
     7           'initialdisp    ','initialtemp    ','velocityfield  ',
     8           'parallel       ','insert         ','return         ',
     9           'coordinatesbin ','               ','               ',
     1           'quad4bin       ','tetra4bin      ','hexa8bin       ',
     2           'quad8bin       ','hydrostatic    ','hydrostress    ',
     3           'constrainpmec  ','initialpres    ','initialstress  ',
     4           '               ','fmaterials     ','end            '/
      data nmc /42/      
c ......................................................................
c
c ...
      init_poro_pmec = .true.
c ......................................................................
c
c ... Leitura dos parametros da malha: nnode,numel,numat,nen,ndf,ndm
c 
c     print*,'Reading mesh'
      ndf  = 0
      ndft = 0
      ndm  = 0
      npi  = 0
      call parameters(dum,numel,numat,dum,ndf,ndft,ndm,npi,nin)
      nst    = maxno*ndf
      nstt   = maxno*ndft
c ......................................................................
c
c ...
      if( ndm .eq. 2) then
        ntn = 4
c ... mec
        if(ndf .eq. 2) then
          fmec = .true.
c ... poro mec
        else if(ndf .eq. 3) then          
          fpmec = .true.
        endif
c .....................................................................
c
c ...
      else
        ntn = 6
c ... mec
        if(ndf .eq. 3) then
          fmec = .true.
c ... poro mec
        else if(ndf .eq. 4) then
          fpmec = .true.
        endif
      endif
c ......................................................................
c
c ... 
      nbar2(1:4)    = 0
      ntria3(1:4)   = 0
      nquad4(1:4)   = 0
      ntetra4(1:4)  = 0 
      nhexa8(1:4)   = 0
      nquad8(1:4)   = 0
      ntetra10(1:4) = 0
      nhexa20(1:4)  = 0
c ......................................................................
c
c ...
c     Alocacao de arranjos na memoria:
c     ---------------------------------------------------------------
c     | ix | id | ie | nload | eload | inum | e | x | f | u | v | a |
c     ---------------------------------------------------------------
c
      i_ix    = alloc_4('ix      ',maxno+1,numel)
      i_ie    = alloc_4('ie      ',    1  ,numat)
      i_e     = alloc_8('e       ',maxprop,numat)
      i_x     = alloc_8('x       ',  ndm  ,nnode)
c ... para o caso onde soh ha condicao mista em termos fluxo 
      if(ndft .ne. 0) then
        mrload(8) = .true.
        i_idt    = alloc_4('idt     ',  ndft,nnode)
        call mzero(ia(i_idt),nnode*ndft)
      endif
c ......................................................................
      call mzero(ia(i_ix),numel*(maxno+1))
      call mzero(ia(i_ie),numat) 
      call azero(ia(i_e),numat*maxprop)
      call azero(ia(i_x),nnode*ndm)
c ......................................................................
c
c ... macro comandos lidas
      do j = 1, ncont
        mrload(j) = .false.
      enddo
c ......................................................................
      totnel  = 0            
      nmacros = 0
c ......................................................................
  100 continue
      call readmacro(nin,.true.)
      write(rc,'(15a)') (word(j),j=1,15)
      do 200 j = 1, nmc
         if (rc .eq. macro(j)) go to 300
  200 continue
      go to 100
  300 continue
c ......................................................................
c
c ... Controle de Macros (PRE-processador):
c
      nmacros = nmacros + 1
      write(macros(nmacros),'(15a)') rc
c ......................................................................
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
     3      2350,2400,2450)j !               ,fmaterials   ,end
c ......................................................................
c
c ... Propriedades dos materiais:
c
  400 continue
      call mate(ia(i_ie),ia(i_e),numat,nin)
      go to 100
c ......................................................................
c
c ... Conetividades bar2:    
c
  450 continue
      nbar2(1) = 0
      if(verbose) print*, 'loading bar2 ...'
      call elconn(ia(i_ix),3,maxno,nbar2(1),numel,nin)
      if(verbose) print*, 'load.'
      nbar2(2) = totnel+1
      totnel   = totnel + nbar2(1)      
      go to 100
c ......................................................................
c
c ... Conetividades tria3:          
c
  500 continue
      f_read_el = .true.
      ntria3(1) = 0
      if(verbose) print*, 'loading tria3 ...'
      call elconn(ia(i_ix),4,maxno,ntria3(1),numel,nin)
      if(verbose) print*, 'load.'
      ntria3(2) = totnel+1
      totnel    = totnel + ntria3(1)
      go to 100
c ......................................................................
c
c ... Conetividades quad4:          
c
  550 continue
      f_read_el = .true.
      nquad4(1) = 0
      if(verbose) print*, 'loading quad4 ...'
      call elconn(ia(i_ix),5,maxno,nquad4(1),numel,nin)
      if(verbose) print*, 'load.'
      nquad4(2) = totnel+1
      totnel    = totnel + nquad4(1)
      go to 100 
c ......................................................................
c
c ... Conetividades tetra4:          
c
  600 continue
      f_read_el = .true.
      ntetra4(1) = 0
      if(verbose) print*, 'loading tetra4 ...'
      call elconn(ia(i_ix),maxno+1,4,ntetra4(1),numel,nin)
      if(verbose) print*, 'load.'
      ntetra4(2) = totnel+1
      totnel    = totnel + ntetra4(1)
c ... transforma os elementos lineares em quadraticos (20 nos)
      if(maxno .eq. 10) then
        el_quad = .true.
c ...  
        i_nincid = alloc_4('nincid  ',1,nnodev) 
        call mynodegrade(ia(i_ix),ia(i_ix),nnodev,numel,maxnov,maxno
     .                 ,ia(i_nincid)  ,maxgrade,.true.) 
        i_incid  = alloc_4('incid   ',maxgrade,nnodev)
        call myelmincid(ia(i_ix),ia(i_ix),ia(i_incid),ia(i_nincid)
     .                 ,nnodev
     .                 ,numel    ,maxnov     ,maxno     
     .                 ,maxgrade ,.true.)
c ... gera a conectividade dos elementos quadraticos
        call mk_elconn_quad_v1(ia(i_ix),ia(i_incid),ia(i_nincid)
     .                        ,numel   ,nnodev     ,nnode 
     .                        ,4       ,maxno      ,maxno+1     
     .                        ,maxgrade)
c .....................................................................
c
c ...
        ntetra10(1:4) = ntetra4(1:4)   
        ntetra4(1:4)  = 0 
c .....................................................................
c
c ...
        i_incid     = dealloc('incid   ')
        i_nincid    = dealloc('nincid  ')
      endif
c .....................................................................
      go to 100
c ......................................................................
c
c ... Conetividades hexa8:          
c
  650 continue
      f_read_el = .true.
      nhexa8(1) = 0
      if(verbose) print*, 'loading hexa8 ...'
      call elconn(ia(i_ix),maxno+1,8,nhexa8(1),numel,nin)
      if(verbose) print*, 'load.'
      nhexa8(2) = totnel+1
      totnel    = totnel + nhexa8(1)
c ... transforma os elementos lineares em quadraticos (20 nos)
      if(maxno .eq. 20) then
        el_quad = .true.
c ...  
        i_nincid = alloc_4('nincid  ',1,nnodev) 
        call mynodegrade(ia(i_ix),ia(i_ix),nnodev,numel,maxnov,maxno
     .                 ,ia(i_nincid)  ,maxgrade,.true.) 
        i_incid  = alloc_4('incid   ',maxgrade,nnodev)
        call myelmincid(ia(i_ix),ia(i_ix),ia(i_incid),ia(i_nincid)
     .                 ,nnodev
     .                 ,numel    ,maxnov     ,maxno     
     .                 ,maxgrade ,.true.)
c ... gera a conectividade dos elementos quadraticos
        call mk_elconn_quad_v1(ia(i_ix),ia(i_incid),ia(i_nincid)
     .                        ,numel   ,nnodev     ,nnode 
     .                        ,8       ,maxno      ,maxno+1     
     .                        ,maxgrade)
c .....................................................................
c
c ...
        nhexa20(1:4) = nhexa8(1:4)   
        nhexa8(1:4)  = 0 
c .....................................................................
c
c ...
        i_incid     = dealloc('incid   ')
        i_nincid    = dealloc('nincid  ')
      endif
c .....................................................................
      go to 100
c ......................................................................
c
c ... Conetividades quad8:          
c
  700 continue
      nquad8(1) = 0
      if(verbose) print*, 'loading quad8 ...'
      call elconn(ia(i_ix),9,8,nquad8(1),numel,nin)
      if(verbose) print*, 'load.'
      nquad8(2) = totnel+1
      totnel    = totnel + nquad8(1)
      go to 100
c ......................................................................
c
c ... Propriedades variaveis com a temperatura:termprop
c
  750 continue
      if(verbose) print*, 'loading termprop  ...'
      call rtermprop(numat,nin)
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
  800 continue
      if(verbose) print*, 'loading adiabat ...'
      call read_adiabat_pre(numat,nin)
      if(verbose) print*, 'load.'
      go to 100                                      
c ......................................................................
c
c ... Coordenadas:
c
  850 continue
      if(verbose) print*, 'loading coordinates ...'
      call coord(ia(i_x),nnodev,ndm,nin)
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
c
c ... constraindisp - restricoes nodais (deslocamentos)
c
  900 continue
      if(verbose) print*, 'loading constraindisp ...'
      mrload(1) = .true.
      if(f_read_el) then
c
        i_id    = alloc_4('id      ',  ndf,nnode)
        call mzero(ia(i_id),nnode*ndf)
c      
        call bound(ia(i_id),nnode,ndf,nin,1)
c ... gerando os carregamentos nos nohs intermediarios quadraticos
        if(el_quad) then
          call mk_bound_quad(ia(i_id),ia(i_ix),numel,ndf,ndf,maxno)
        endif       
       else
        print*,'MACRO: constraindisp !! elementos nao lidos'
      endif  
      if(verbose) print*, 'load.'   
      go to 100
c ......................................................................
c
c ... constraintemp -restricoes nodais (temperaturas):
c
  950 continue
      if(verbose) print*, 'loading constraintemp ...'
c      
      call bound(ia(i_idt),nnode,ndft,nin,1)
      if(verbose) print*, 'load.'
      go to 100      
      
c ......................................................................
c
c ... nodalforces - forcas nodais:
c
 1000 continue
      if(verbose) print*, 'loading nodalforces ...'
      mrload(4) = .true.
      if(verbose) print*, 'load.'
c
      if(f_read_el) then
        i_f     = alloc_8('f       ',  ndf,nnode)
        call azero(ia(i_f),nnode*ndf)
c
        call forces(ia(i_f),nnodev,ndf,nin)
c ...
        if(el_quad) then
          call mk_forces_quad(ia(i_id),ia(i_f),ia(i_ix)
     .                       ,numel,ndf,ndf,maxno)
        endif
c ......................................................................
      else
        print*,'MACRO: nodalforces !! elementos nao lidos'
      endif
      if(my_id .eq. 0) print*,'load.'
      go to 100
c ......................................................................
c
c ... nodalsources - fontes termicas nodais:
c
 1050 continue
      if(verbose) print*, 'loading nodalsources ...'
      mrload(11) = .true.
      i_ft     = alloc_8('ft      ',  ndft,nnode)
      call azero(ia(i_ft),nnode*ndft)
c      
      call forces(ia(i_ft),nnode,ndft,nin)
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
c
c ... nodalloads - nos com cargas variaveis no tempo:
c
 1100 continue
      if(verbose) print*, 'loading nodalloads ...'
      mrload(2) = .true.
      i_nload = alloc_4('nload   ',  ndf,nnode)
      call mzero(ia(i_nload),nnode*ndf)
c      
      call bound(ia(i_nload),nnode,ndf,nin,2) 
      if(verbose) print*, 'load.'
      goto 100
c ......................................................................
c
c ... nodalthermloads - nos com cargas termicas variaveis no tempo:
c
 1150 continue
      if(verbose) print*, 'loading nodalthermloads ...'
      mrload(9) = .true.
      i_nloadt = alloc_4('nloadt  ',  ndft,nnode)
      call mzero(ia(i_nloadt),nnode*ndft)
c      
      call bound(ia(i_nloadt),nnode,ndft,nin,2) 
      if(verbose) print*, 'load.'
      goto 100
c ......................................................................
c
c ... elmtloads - cargas nos elementos:
c      
 1200 continue
      if(verbose) print*, 'loading elmtloads ...'
      mrload(3) = .true.
      i_eload = alloc_4('eload   ',maxface+1,numel)
      call mzero(ia(i_eload),numel*(maxface+1))
c      
      call bound(ia(i_eload),numel,maxface+1,nin,3) 
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
c
c ... elmtthermloads - cargas termicas nos elementos:                                         
c      
 1250 continue
      if(verbose) print*, 'loading elmtthermloads ...'
      mrload(10) = .true.
      i_eloadt = alloc_4('eloadt  ',maxface+1,numel)
      call mzero(ia(i_eloadt),numel*(maxface+1)) 
c      
      call bound(ia(i_eloadt),numel,(maxface+1),nin,3) 
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
c
c ... Definicao das cargas variaveis no tempo:
c
 1300 continue
      if(verbose) print*, 'loading loads ...'
      call rload(nin)
      if(verbose) print*, 'load.'
      goto 100
c ......................................................................
c
c ... hidrprop                                
c
 1350 continue
      if(verbose) print*, 'loading hidrprop ...'
      call rtermprop(numat,nin)
      if(verbose) print*, 'load.'
      goto 100
c ......................................................................
c
c ... Conetividades prism6:          
c 
 1400 continue
      nprism6(1) = 0
      if(verbose) print*, 'loading prism6 ...'
      call elconn(ia(i_ix),7,6,nprism6(1),numel,nin)
      if(verbose) print*, 'load.'
      nprism6(2) = totnel+1
      totnel     = totnel + nprism6(1)
      goto 100
c ......................................................................
c
c ... initialdisp - deslocamentos iniciais:
c
 1450 continue
      if(verbose) print*, 'loading initialdisp ...'
      mrload(5) = .true.
      i_u     = alloc_8('u       ',  ndf,nnode)
      call azero(ia(i_u),nnode*ndf)
c      
      call forces(ia(i_u),nnode,ndf,nin)  
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
c
c ... initialtemp - temperaturas iniciais:
c
 1500 continue
      if(verbose) print*, 'loading initialtemp ...'
      mrload(12) = .true.
      i_ut0    = alloc_8('ut0     ',  ndft,nnode)         
      call azero(ia(i_ut0),nnode*ndft)  
c      
      call forces(ia(i_ut0),nnode,ndft,nin)  
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
c
c ... velocityfield - campo de velocidades de conveccao:                                         
c
 1550 continue
      if(verbose) print*, 'loading velocityfield ...'
      mrload(14) = .true.
      i_w      = alloc_8('w       ',   ndm,nnode)
      call azero(ia(i_w),  nnode*ndm)      
c
      call coord(ia(i_w),nnode,ndm,nin)
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
c
c ... Paralelo:                                         
c      
 1600 continue
      call read_par(nin,nnode,numel,ndf,ndft) 
      go to 100 
c ......................................................................
c
c ... (insert) Desvia leitura para arquivo auxiliar:
c      
 1650 continue
      nmacros = nmacros - 1
      naux = nin
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,80)
      open(nincl, file=fname,status='old',err=1655)
      nin = nincl
      go to 100
 1655 continue  
      print*,'File ',trim(fname),' not found !'
      call finalize(ierr)
c ......................................................................
c
c ... (return) Retorna leitura para arquivo de dados basico:
c      
 1700 continue
      nmacros = nmacros - 1
      close(nincl)
      nin = naux
      go to 100
c ......................................................................
c      
c ... (coordenadas no formato binario arquivo auxiliar)      
c      
 1750 continue
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(unit=nincl,file=fname,access='stream'
     .    ,form='unformatted',convert='big_endian')
      call coordbin(ia(i_x),nnode,ndm,nincl)
      close(nincl)
      go to 100                  
c ......................................................................
c
 1800 continue
      go to 100                  
c ......................................................................
c
 1850 continue
      go to 100                  
c ......................................................................
c
c ... quad4 bin (elementos no formato binario arquivo auxiliar)
 1900 continue
      nquad4(1) = 0
c ... numero de elementos quad4      
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1910,end =1910) nquad4(1)
c ... nome do arquivo bin      
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(unit=nincl,file=fname,access='stream'
     .    ,form='unformatted',convert='big_endian')
c .....................................................................
c
c ...     
      call elmbin(ia(i_ix),nquad4(1),4,nincl)
c .....................................................................
c
c ...
      nquad4(2) = totnel + 1
      totnel    = totnel + nquad4(1)
c .....................................................................
c
c ...
      close(nincl)
      go to 100
c .....................................................................
c
c ...      
 1910 continue 
      print*,"Erro na leitura do numeros de elementos quad4(QUAD4BIN)!"
      call finalize()
c ......................................................................
c
c ... tetra4bin (elementos no formato binario arquivo auxiliar)
 1950 continue
      ntetra4(1) = 0
c ... numero de elementos quad4      
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1960,end =1960) ntetra4(1)
c ... nome do arquivo bin      
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(unit=nincl,file=fname,access='stream'
     .    ,form='unformatted',convert='big_endian')
c .....................................................................
c
c ...     
      call elmbin(ia(i_ix),numel,4,nincl)
c .....................................................................
c
c ...
      ntetra4(2) = totnel + 1
      totnel     = totnel + ntetra4(1)
c .....................................................................
c
c ...
      close(nincl)
      go to 100                  
c .....................................................................
c
c ...      
 1960 continue 
      print*,
     .      "Erro na leitura do numeros de elementos tetra4(TETRA4BIN)!"
      call finalize()
c ......................................................................
c
c ......................................................................
c
c ... hexa8 bin (elementos no formato binario arquivo auxiliar)
 2000 continue
      nhexa8(1) = 0
c ... numero de elementos quad8      
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2010,end =2010) nhexa8(1)
c ... nome do arquivo bin      
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(unit=nincl,file=fname,access='stream'
     .    ,form='unformatted',convert='big_endian')
c .....................................................................
c
c ...
      call elmbin(ia(i_ix),nhexa8(1),8,nincl)
c .....................................................................
c
c ...
      nhexa8(2) = totnel + 1
      totnel    = totnel + nhexa8(1)
c .....................................................................
c
c ...
      close(nincl)
      go to 100                  
c ...      
 2010 continue 
      print*,"Erro na leitura do numeros de elementos hexa8(HEXA8BIN)!"
      call finalize()
c ......................................................................
c
c ... quad8 bin (elementos no formato binario arquivo auxiliar)
 2050 continue
      nquad8(1) = 0
c ... numero de elementos quad8      
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2060,end =2060) nquad8(1)
c ... nome do arquivo bin      
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(unit=nincl,file=fname,access='stream'
     .    ,form='unformatted',convert='big_endian')
c .....................................................................
c
c ...     
      call elmbin(ia(i_ix),numel,8,nincl)
c .....................................................................
c
c ...
      nquad8(2) = totnel+1
      totnel    = totnel + nquad8(1)
c ......................................................................
c
c ...
      close(nincl)
      go to 100                  
c ......................................................................
c
c ...      
 2060 continue 
      print*,
     .      "Erro na leitura do numeros de elementos quad8(quad8BIN)!"
      call finalize()
c ......................................................................
c
c ...
 2100 continue
      if(verbose) print*, 'hydrostatic ...'
      call readmacro(nin,.false.)
      write(name_hidrostatic(1),'(80a)') (word(j),j=1,strl)
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
c
c ...
 2150 continue
      if(verbose) print*, 'hydrostress ...'
      call readmacro(nin,.false.)
      write(name_hidrostatic(2),'(80a)') (word(j),j=1,strl)
      if(verbose) print*, 'load.'
      go to 100
c ......................................................................
c
c ... constraintemp:
 2200 continue
      mrload(15) = .true.
      if(verbose) print*, 'loading constrainpmec ...'
c     mrload(1) = .true.
      if(f_read_el) then
c
        i_id    = alloc_4('id      ',  ndf,nnode)
        call mzero(ia(i_id),nnode*ndf)
c      
        call bound(ia(i_id),nnode,ndf,nin,1)
c ... gerando os carregamentos nos nohs intermediarios quadraticos
        if(el_quad) then
          call mk_bound_quad(ia(i_id),ia(i_ix),numel,ndf-1,ndf,maxno)
        endif       
       else
        print*,'MACRO: constrainpmec !! elementos nao lidos'
      endif  
      if(verbose) print*, 'load.'   
      go to 100 
c ......................................................................
c
c ... 
 2250 continue
      mrload(16) = .true.
      if(verbose) print*,'load intialpres ...'
      if(f_read_el) then
        if(init_poro_pmec) i_u     = alloc_8('u       ',  ndf,nnode)
        call azero(ia(i_u),nnode*ndf)
        call init_poro_mec(ia(i_u),nnode,ndf,ndf,ndf,nin)
      else
        print*,'MACRO: initialpres !! elementos nao lidos'
      endif
      if(my_id .eq. 0) print*,'load.'
      go to 100
c ......................................................................
c
c ... 
 2300 continue
      mrload(17) = .true.
      if(verbose) print*,'load intialstress ...'
      if(f_read_el) then
        i_tx0   = alloc_8('tx0     ',  ntn,nnode)
        call azero(ia(i_tx0),nnode*ntn)
        call forces(ia(i_tx0),nnode,ntn,nin)  
c ... gerando as tensoes nos nohs intermediarios quadraticos
        if(el_quad) then
          call mk_initial_quad(ia(i_tx0),ia(i_ix),numel,ntn,maxno)
        endif   
      else
        print*,'MACRO: initialstress !! elementos nao lidos'
      endif
      if(my_id .eq. 0) print*,'load.'
      go to 100
c ......................................................................
c
c ...
 2350 continue
      go to 100  
c ......................................................................
c
c ...
 2400 continue
      call rfmate(file_prop,ia(i_ie),ia(i_e),numat,nin)
      go to 100                  
c ......................................................................
c
c ... End:
c
 2450 continue
c
c ... Inclui macro de paralelo (PRE-processador):
c      
      if (nprcs .gt. 1) then
         write(macros(nmacros),'(15a)') 'parallel'
         nmacros = nmacros + 1
         write(macros(nmacros),'(15a)') 'end'
      endif
c ......................................................................
c
c ...
      call mk_coor_quad(ia(i_x),ia(i_ix),numel,maxno,ndm)  
c ......................................................................
      return
      end
      subroutine elconn(ix,nen1,nen,nel,numel,nin)
c **********************************************************************
c *                                                                    *
c *   Conetividades nodais.                                            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'      
      integer ix(nen1,*),nen1,nen,nel,numel,nin,j,k,m
      character*12 string
c ......................................................................
      do 100 j = 1, numel
         read(nin,*) k,(ix(m,k),m=1,nen),ix(nen1,k)
  100 continue
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)  
      if (string .ne. 'end') goto 200
      nel=numel
      return
c ......................................................................
c     nel = 0
c     call readmacro(nin,.true.)
c     write(string,'(12a)') (word(j),j=1,12)
c     do 100 while(string .ne. 'end')
c       read(string,*,err = 200,end = 200) k
c       if(k .lt. 1 .or. k .gt. numel) goto 200      
c       do 10 j = 1, nen
c         call readmacro(nin,.false.)
c         write(string,'(12a)') (word(m),m=1,12)
c         read(string,*,err = 200,end = 200) ix(j,k)
c  10   continue
c       call readmacro(nin,.false.)
c       write(string,'(12a)') (word(m),m=1,12)
c       read(string,*,err = 200,end = 200) ix(nen1,k)   
c       nel = nel + 1
c       call readmacro(nin,.true.)
c       write(string,'(12a)') (word(j),j=1,12)
c 100 continue
c     return
c ......................................................................
  200 continue
      print*,'*** Error reading the elements!'
      call finalize()                   
      end
      subroutine coord(x,nnode,ndm,nin)
c **********************************************************************
c *                                                                    *
c *   Coordenadas.                                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer nnode,ndm,nin,j,k,m
      real*8 x(ndm,nnode)
      character*30 string
c ......................................................................
      do 100 j = 1, nnode
         read(nin,*) k,(x(m,k),m=1,ndm)
  100 continue
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
      if (string .ne. 'end') goto 200
      return
c ......................................................................
c      call readmacro(nin,.true.)
c      write(string,'(12a)') (word(j),j=1,12)
c      do 100 while(string .ne. 'end')
c         read(string,*,err = 200,end = 200) k
c         if(k .lt. 1 .or. k .gt. nnode) goto 200
c         do 10 j = 1, ndm
c            call readmacro(nin,.false.)
c            write(string,'(30a)') (word(m),m=1,30)
c            read(string,*,err = 200,end = 200) x(j,k)         
c   10    continue
c         call readmacro(nin,.true.)
c         write(string,'(12a)') (word(j),j=1,12)   
c  100 continue
c      return
c ......................................................................                        
  200 continue
      print*,'*** Error in reading the coordinates!'
      call finalize()             
      end
c **********************************************************************
c      
      subroutine coordbin(x,nnode,ndm,nin)
c **********************************************************************
c *                                                                    *
c *   Coordenadas em arquivo auxiliar binario                          *
c *                                                                    *
c **********************************************************************
      implicit none
c ===      
      real*8 x(*)
      integer ndm,nnode,nin
      integer i,j,no,k
c =====================================================================      
c
c ===      
      do i = 1,nnode
        read(nin)no
c       print*,no
        k  = (no-1)*ndm
        do j = 1 , ndm
          read(nin)x(k+j)
c         print*,x(k+j)
        enddo
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
c      
      subroutine elmbin(ix,numel,nen,nin)
c **********************************************************************
c *                                                                    *
c *   Coordenadas em arquivo auxiliar binario                          *
c *                                                                    *
c **********************************************************************
      implicit none
c ===      
      integer ix(nen+1,*)
      integer nen,numel,nin
      integer i,j,no,k
c =====================================================================      
c
c ===      
      do i = 1,numel
        read(nin)no
        do j = 1 , nen
          read(nin)ix(j,no)
c         print*,no,ix(j,no)
        enddo
c ... material        
        read(nin)ix(nen+1,no)
c       print*,ix(nen+1,no)
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
      subroutine mate(ie,e,numat,nin)
c **********************************************************************
c *                                                                    *
c *   Materiais.                                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      include 'termprop.fi'
      integer ie(*),numat,nin,i,j,m,ma
      real*8  e(maxprop,*)
      character*30 string
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
c .....................................................         
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) ma
         if(ma .lt. 1 .or. ma .gt. numat) goto 200
c .....................................................            
         call readmacro(nin,.false.)
         write(string,'(12a)') (word(m),m=1,12)
         read(string,*,err = 200,end = 200) ie(ma)
c .....................................................     
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(m),m=1,30)
         i = 0
c
c ... linha de comando original:
c         do 100 while(string .ne. ' ')
c
c ... formato necessario para funcionamento em linux:
         do 100 while ( (string .ne.   ' '    ) .and.
     .                  (string .ne. CHAR(13) ) )
            i = i + 1
            if(i .gt. maxprop) goto 200
            read(string,*,err = 200,end = 200) e(i,ma)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(m),m=1,30)
  100    continue
         read_prop(ma) = i
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................                        
  200 continue
      print*,'*** Erro na leitura dos materiais !'
      call finalize()             
      end
c **********************************************************************  
c
c **********************************************************************  
      subroutine rfmate(file,ie,e,numat,nin)
c **********************************************************************
c * Data de criacao    : 20/03/2017                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * RFMATE :leitura dos nomes dos arquivos de matetias                 * 
c * ------------------------------------------------------------------ *
c * Paramtros de entrada:                                              *
c * ------------------------------------------------------------------ *
c * file   -> nao definido                                             *
c * numat  -> numero de materias                                       *
c * nin    -> arquivo de leitura                                       *
c * ------------------------------------------------------------------ *
c * Paramtros de  Saida:                                               *
c * ------------------------------------------------------------------ *
c * file   -> nomer dos arquivos de materiais                          *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * Obs:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'string.fi'
      include 'termprop.fi'
      integer ie(*),numat,nin,j,m,ma
      real*8 e(maxprop,*)
      character*80 file(*),string
      integer nincl
      data nincl /12/ 
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(80a)') (word(j),j=1,80)
c ......................................................................        
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) ma
         if(ma .lt. 1 .or. ma .gt. numat) goto 200
c ......................................................................            
         call readmacro(nin,.false.)
         write(file(ma),'(80a)') (word(m),m=1,80)
c ......................................................................
         open(nincl, file= file(ma),status= 'old',err=900,action='read')
         call file_prop(ie(ma),e(1,ma),nincl)
         close(nincl)
c ......................................................................
         call readmacro(nin,.true.)
         write(string,'(80a)') (word(j),j=1,80)
  110 continue
      return
c ......................................................................                        
  200 continue
      print*,'*** Erro na leitura dos materiais !'
      call finalize()       
c ...
 900  continue
      print*,'File ',trim(file(ma)),' not found !'
      call finalize() 
c .....................................................................      
      end      
c **********************************************************************
c 
c **********************************************************************
      subroutine forces(f,nnode,ndf,nin)
c **********************************************************************
c *                                                                    *
c *   Forcas nodais e valores prescritos.                              *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'      
      integer nnode,ndf,nin,i,j,k
      real*8 f(ndf,*)
      character*30 string            
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) k
         if(k .lt. 1 .or. k .gt. nnode) goto 200      
         do 10 j = 1, ndf
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(i),i=1,30)
            read(string,*,err = 200,end = 200) f(j,k)
   10    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(i),i=1,12)
  100 continue
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura das forcas nodais !'
      call finalize()             
      end
      subroutine bound(id,nnode,ndf,nin,code)
c **********************************************************************
c *                                                                    *
c *   Coordenadas.                                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer id(ndf,*),nnode,ndf,nin,code,i,j,k
      character*30 string
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
c .....................................................         
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) j
         if(j .lt. 1 .or. j .gt. nnode) goto 200
c .....................................................            
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(k),k=1,30)
         i = 0
c
c ... linha de comando original:
c         do 100 while(string .ne. ' ')
c
c ... formato necessário ao funcionamento em linux:
         do 100 while ( (string .ne.   ' '    ) .and.
     .                  (string .ne. CHAR(13) ) )
            i = i + 1
            if(i .gt. ndf) goto 200
            read(string,*,err = 200,end = 200) id(i,j)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(k),k=1,30)               
  100    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................                        
  200 continue
      if    (code .eq. 1) then
         print*,'*** Erro na leitura das restricoes !'
      elseif(code .eq. 2) then
         print*,'*** Erro na leitura das cargas nodais !'
      elseif(code .eq. 3) then
         print*,'*** Erro na leitura das cargas nos elementos !'
      endif
      call finalize()             
      end            
      subroutine bound0(id,nnode,ndf,nin)
c **********************************************************************
c *                                                                    *
c *   Restricoes nodais.                                               *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer id(ndf,*),nnode,ndf,nin,i,j,k
      character*12 string      
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) k
         if(k .lt. 1 .or. k .gt. nnode) goto 200      
         do 10 j = 1, ndf
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(i),i=1,12)
            read(string,*,err = 200,end = 200) id(j,k)
   10    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(i),i=1,12)
  100 continue
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura das restricoes nodais !'
      call finalize()             
      end
      subroutine rload(nin)
c **********************************************************************
c *                                                                    *
c *   Leitura das cargas variaveis no tempo                            *
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
      implicit none
      include 'string.fi'
      include 'load.fi'
      integer  nin,i,j,k
      character*30 string            
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
c ...    numero da carga:
         read(string,*,err = 200,end = 200) i
         if(i .lt. 1 .or. i .gt. numload) goto 200 
         call readmacro(nin,.false.)
         write(string,'(12a)') (word(j),j=1,12)
c ...    tipo da carga:
         read(string,*,err = 200,end = 200) load(1,i)
         load(2,i) = 1
         if     (load(1,i) .eq. 1) then
c ...       valor da carga:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(1,1,i)
c .....................................................................
c
c ...
         elseif (load(1,i) .eq. 2) then
c ...       u(t) = a.t + b 
c ...       valor de b:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(1,1,i)
c ...       valor de a:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(2,1,i)
c .....................................................................
c
c ...
         elseif (load(1,i) .eq. 3) then
c ...       -kdu/dx = H.(u-uext) 
c ...       valor de uext:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(1,1,i)
c ...       valor de H:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(2,1,i)
c .....................................................................
c
c ...
         elseif (load(1,i) .eq. 6) then
c ...       numero de parcelas:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
            do k = 1, load(2,i)
               call readmacro(nin,.true.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(1,k,i)
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(2,k,i)
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(3,k,i)
            enddo
c .....................................................................
c
c ...
         elseif (load(1,i) .eq. 7) then
c ...       kdu/dx = emiss * const(Stef-Boltz) *(u4-uext4)+H(uext-u)
c ...       Emissividade:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(1,1,i)
c ...       Constante de Stefan-Boltzmann:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(2,1,i)
c ...       constante de conveccao (h):            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(3,1,i)
c ...       numero de parcelas (uext no tempo):            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
            do k = 1, load(2,i)
c ...       valor de tempo: 
               call readmacro(nin,.true.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,2,i)
c ...       valor de uext:
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,3,i)
            enddo   
c .....................................................................
c
c ...
         elseif (load(1,i) .eq. 8) then
c ...       carga constante a partir de um dado t0
c ...       valor da carga: 
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
            do k = 1, load(2,i) 
              call readmacro(nin,.true.)
              write(string,'(30a)') (word(j),j=1,30)   
              read(string,*,err = 200,end = 200) fload(k,1,i)  
              call readmacro(nin,.false.)
              write(string,'(30a)') (word(j),j=1,30)
              read(string,*,err = 200,end = 200) fload(k,2,i) 
            enddo
c .....................................................................
c
c ...
         elseif (load(1,i) .eq. 9) then
c ...       interpolacao de pontos
c           numero de parcelas:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
            do k = 1, load(2,i)
               call readmacro(nin,.true.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,1,i)
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,2,i)
            enddo
c .....................................................................
c
c ...  
         elseif (load(1,i) .eq. 10) then
c ...    troca de calor de Newton(com uExt e h variavel):            
c ...       numero de parcelas (uext no tempo):            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
            do k = 1, load(2,i)
c ...       valor de tempo: 
               call readmacro(nin,.true.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,1,i)
c ...       valor de h: 
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,2,i)
c ...       valor de uext:
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,3,i)
            enddo   
c .....................................................................
c
c ...
         elseif (load(1,i) .eq. 17) then
c ...   fonte de geofisica
c ...       numero de parcelas:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
c ... amplitude            
           call readmacro(nin,.false.)
           write(string,'(30a)') (word(j),j=1,30)
           read(string,*,err = 200,end = 200) fload(1,1,i)          
c ... valor da frequencia de corte            
           call readmacro(nin,.false.)
           write(string,'(30a)') (word(j),j=1,30)
           read(string,*,err = 200,end = 200) fload(2,1,i)          
         elseif (load(1,i) .eq. 18) then
c ...   fonte de senoidal 
c ...       numero de parcelas:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
c ... amplitude            
           call readmacro(nin,.false.)
           write(string,'(30a)') (word(j),j=1,30)
           read(string,*,err = 200,end = 200) fload(1,1,i)
c ... valor da frequencia           
           call readmacro(nin,.false.)
           write(string,'(30a)') (word(j),j=1,30)
           read(string,*,err = 200,end = 200) fload(2,1,i)
c ... fase                          
           call readmacro(nin,.false.)
           write(string,'(30a)') (word(j),j=1,30)
           read(string,*,err = 200,end = 200) fload(3,1,i)
c ... numero de peridos             
           call readmacro(nin,.false.)
           write(string,'(30a)') (word(j),j=1,30)
           read(string,*,err = 200,end = 200) fload(1,2,i)
c .....................................................................
c
c ... forca distribuida no contorno
         elseif (load(1,i) .eq. 40) then
c ... nome do arquivo de cargas
            call readmacro(nin,.false.)
            write(name_loads(i),'(80a)') (word(j),j=1,strl)
c .....................................................................
c
c ... forca distribuida normal ao contorno 
         elseif (load(1,i) .eq. 41) then
c ... nome do arquivo de cargas
            call readmacro(nin,.false.)
            write(name_loads(i),'(80a)') (word(j),j=1,strl)
c .....................................................................
c
c ... pressao hidrostatica               o 
         elseif (load(1,i) .eq. 42) then
c ... nome do arquivo de cargas
            call readmacro(nin,.false.)
            write(name_loads(i),'(80a)') (word(j),j=1,strl)
c .....................................................................
c
c ...       
         else     
           print*,'***tipo de carregamento nao existente!!',load(1,i) 
           call finalize()
c .....................................................................
         endif
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  100 continue
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura das cargas !'
      call finalize()
      end
      subroutine rtermprop(numat,nin)
c **********************************************************************
c *                                                                    *
c *   Leitura das propriedades variaveis com a temperatura             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    idl(6) - variavel auxuliar                                      *
c *    nin - arquivo de entrada                                        *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    nprop (i,j,k) - valor i da propriedade j do material tipo k     *
c *    eprop (j,k) - numero de termos da propriedade j do material tipo*
c *                 k                                                  *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      include 'termprop.fi'
      integer  nin,i,j,k,ma, numat,m
      character*30 string            
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
c ...    numero do material:
         read(string,*,err = 200,end = 200) ma
         if(ma .lt. 1 .or. ma .gt. numat) goto 200 
         call readmacro(nin,.false.)
         write(string,'(12a)') (word(j),j=1,12)
c ...    numero da propriedade variavel no material ma:
         read(string,*,err = 200,end = 200) i
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(j),j=1,30)
c ...    numero de pontos da poligonal:
         read(string,*,err = 200,end = 200) eprop(i,ma)
         k = eprop(i,ma)
         do m = 1, k
c ...       valor da temperatura no ponto j:            
             call readmacro(nin,.true.)
             write(string,'(30a)') (word(j),j=1,30)
             read(string,*,err = 200,end = 200) nprop(m,i,ma)
             call readmacro(nin,.false.)
             write(string,'(30a)') (word(j),j=1,30)
c ...       valor da propriedade no ponto j:  
             read(string,*,err = 200,end = 200) nprop(m+k,i,ma)
         enddo
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
  100 continue
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura das propriedades variaveis !'
      call finalize()
      end
      subroutine read_adiabat_pre(numat,nin)
c **********************************************************************
c *                                                                    *
c *   Leitura da curva de elevacao adiabatica e calculo da afinidade   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    idl(6) - variavel auxuliar                                      *
c *    nin - arquivo de entrada                                        *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    nprop (i,9,k) - valor i da da afinidade do material tipo k      *
c *    eprop (9,k) - numero de termos da afinidade do material tipo    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      include 'termprop.fi'
      integer  nin,i,j,k,ma, numat,m
      character*30 string     
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
c ...    numero do material:
         read(string,*,err = 200,end = 200) ma
         if(ma .lt. 1 .or. ma .gt. numat) goto 200 
         call readmacro(nin,.false.)
         write(string,'(12a)') (word(j),j=1,12)
c ...    numero de pontos da poligonal:
         read(string,*,err = 200,end = 200) eprop(9,ma)
         k = eprop(9,ma)
         do m = 1, k
c ...       Leitura do tempo (dias) e da Temperatura(oC):
           read(nin,*) tad(m,ma),tempad(m,ma)
         enddo
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
  100 continue
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura da curva de elevacao adiabatica!'
      call finalize()
      end    
      subroutine boundc(nnode,ndf,id,f,u)
c **********************************************************************
c *                                                                    *
c *   Valores prescritos.                                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 f(ndf,*),u(ndf,*)
c ......................................................................
      do 200 i = 1, nnode
         do 100 j = 1, ndf
            k = id(j,i)
            if(k .gt. 0) then
               u(j,i) = f(j,i)
            endif
  100    continue
  200 continue
c ......................................................................
      return
      end      
      subroutine initc(nnode,ndf,id,u0,u)
c **********************************************************************
c *                                                                    *
c *   Valores prescritos.                                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 u0(ndf,*),u(ndf,*)
c ......................................................................
      do 200 i = 1, nnode
         do 100 j = 1, ndf
            k = id(j,i)
            if(k .le. 0) then
               u(j,i) = u0(j,i)
            endif
  100    continue
  200 continue
c ......................................................................  
      return
      end   
      subroutine parameters(nnode,numel,numat,maxno,ndf,ndft,ndm,npi
     .                     ,nin)
c **********************************************************************
c *                                                                    *
c *   Parameters                                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      include 'readfile.fi'
      character*12 string
      integer nnode,numel,numat,maxno,ndf,ndft,ndm,npi,nin,n,j
c ...
      read_parameter(1:num_var) = .false.
c ......................................................................
      n   = 0
      call readmacro(nin,.true.)
      write(string,'(8a)') (word(j),j=1,8)
      do while (strl .ne. 0)
         if     (string .eq. 'nnode') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) nnode
            read_parameter(1) = .true.
         elseif (string .eq. 'numel') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) numel
            read_parameter(2) = .true.
         elseif (string .eq. 'numat') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) numat
            read_parameter(3) = .true.
         elseif (string .eq. 'maxno') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) maxno
            read_parameter(4) = .true.
         elseif (string .eq. 'ndf') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) ndf
            read_parameter(5) = .true.
         elseif (string .eq. 'therm') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) ndft
            read_parameter(6) = .true.
         elseif (string .eq. 'dim') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) ndm
            read_parameter(7) = .true.
         elseif (string .eq. 'npi') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) npi
            read_parameter(8) = .true.
         endif
         call readmacro(nin,.false.)
         write(string,'(8a)') (word(j),j=1,8)
      end do
c ......................................................................
      return
c ......................................................................
  100 continue
      call finalize()
      end
      subroutine readmacro(nin,newline)
c **********************************************************************
c *                                                                    *
c *   Subroutine READMACRO: le uma macro numa linha nova ou a partir   *
c *                         da coluna line_col de uma linha lida       *
c *                         anteriormente e armazenada em line(*)      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   nin - numero do arquivo de entrada                               *
c *   newline = .true.  - nova linha deve ser lida                     *
c *           = .false. - macro deve ser lida a partir da coluna       *
c *                       line_col em line(*)                          *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer j,k,nin
      logical newline
c ......................................................................
      if(newline) then
         line_col = 1
         read(nin,'(500a1)',err = 200,end = 200) (line(j),j=1,maxstrl)
      endif
c ......................................................................      
      do j = 1, maxstrl
         word(j) = ' '
      enddo
      strl = 0
      if(line_col .gt. maxstrl) return
c ......................................................................      
      do while (line(line_col) .eq. ' ')
         line_col = line_col + 1
         if (line_col .gt. maxstrl) return         
      end do
c ......................................................................      
      do while ( (line(line_col) .ne. ' ') .and.
     .           (line(line_col) .ne. CHAR(13)) )
         strl = strl + 1
         line_col = line_col + 1
         if (line_col .gt. maxstrl) goto 100
      end do
c ......................................................................      
  100 continue      
      k = line_col-strl
      do j = 1, strl
         write(word(j),'(a)') line(k)
         k = k + 1
      end do      
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura do arquivo de dados !'
      call finalize()             
c ......................................................................      
      end 
c ======================================================================      
c ======================================================================
c
c     Leitura da estrutura de dados do paralelo
c
c ======================================================================      
      subroutine read_par(nin,nnode,numel,ndf,ndft)
c **********************************************************************
c *                                                                    *
c *   READ_PAR                                                         *
c *   --------                                                         *
c *                                                                    *
c *   Le informacoes da paralelizacao geradas pelo pre-processador     *
c *   atraves do macro-comando read_par em rdat e cria os arranjos:    *
c *                                                                    *
c *      noLG(nnode) - mapa Local->Global de nos                       *
c *      noGL(nnoG)  - mapa Global-<Local de nos                       *
c *      elLG(numel) - mapa Local->Global de elementos                 *
c *      fmap(nnof)  - mapa Interface->Global de nos                   *
c *      idG(ndf,nnoG) - restricoes nodais globais                     *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      include 'parallel.fi'
      include 'elementos.fi'
      integer nin,nnode,numel,ndf,ndft,i,j,k
      character*160 comando,string
c ......................................................................     
      read(nin,*) string,nnoG      
      read(nin,*) string,nelG      
c ......................................................................
c
c ... Metodo de sub-divisao de dominio
c
      ovlp  = .true.
      novlp = .false.
      read(nin,*)  comando
      if(comando(1:3) .eq. 'non') then
         ovlp  = .false.
         novlp = .true.
      endif
c ......................................................................
c
c ... Numeros de no's do processo, parte non-overlapping
c
      read(nin,'(a80)')  comando
      read(comando,*) string,nno1,
     .                string,nno2,
     .                string,nno3,
     .                string,nno4,
     .                string,nno1a
      if (ovlp)read(nin,*) string,numel_ov
c ......................................................................
c
c ... Numeros de elementos do processo, parte non-overlapping
c
      numel_nov = numel - numel_ov
c ......................................................................
c
c ... {noLG} = Mapa Local -> Global de no's
c
      read(nin,*) string
      i_noLG = alloc_4('noLG    ', 1, nnode)
      call Levet_4(ia(i_noLG),nnode,nin)
      i_noGL = alloc_4('noGL    ', 1, nnoG)
      call mzero(ia(i_noGL),nnoG)      
      do i = 1, nnode
         j = ia(i_noLG+i-1)
         ia(i_noGL+j-1) = i
      enddo        
c ......................................................................
c
c ... {noGL} = Mapa Global -> Local de no's
c
c      read(nin,*) string,nnoG
c      i_noGL = alloc_4('noGL    ', 1, nnoG)
c      call Lemtx_4(ia(i_noGL),1,nin)
c ......................................................................
c
c ..... Numeros de elementos pelo tipo
c
      read(nin,'(a)')  comando
      read(comando,*) string, nbar2(1)  , nbar2(2)  ,
     .                string, ntria3(1) , ntria3(2) ,
     .                string, nquad4(1) , nquad4(2) ,
     .                string, ntetra4(1), ntetra4(2),
     .                string, nhexa8(1) , nhexa8(2)
c ......................................................................
c
c ... {elLG} = Mapa Local -> Global de elementos
c
c      read(nin,*) string,nelG
      read(nin,*) string
      i_elLG = alloc_4('elLG    ', 1, numel)
      call Levet_4(ia(i_elLG),numel,nin)
c ......................................................................
c
c ... {rcvs0} = lista dos tamanhos dos blocos de nós {Vfi}
c
c      i_rcvs0 = alloc_4('rcvs0   ', 1, nprcs)
c      if (ovlp) then
c        neqfi=nno1a+nno2
c      else
c        neqfi=nno2+nno3
c      endif
c      call MPI_allgather(neqfi,1,MPI_INTEGER,ia(i_rcvs0),1,MPI_INTEGER,
c     .                   MPI_COMM_WORLD,ierr)
c ......................................................................
c
c ... {dspl0} = Ponteiro p/ início de cada bloco em {Vf}
c
c      i_dspl0 = alloc_4('dspl0   ', 1, nprcs)
c      ia(i_dspl0) = 0
c      do i = 2, nprcs
c        ia(i_dspl0+i-1) = ia(i_dspl0+i-2) + ia(i_rcvs0+i-2)
c      enddo
c ......................................................................
c
c ... {fmap0} = Mapa Interface -> Global de no's
c
c      read(nin,*) string,nnof
c      i_fmap0 = alloc_4('fmap0   ', 1, nnof)
c      call Levet_4(ia(i_fmap0),nnof,nin)
      read(nin,*) string,nnof1,nnof2
      i_fmap0 = alloc_4('fmap0   ', 1, nnof1+nnof2)
      call Levet_4(ia(i_fmap0),nnof1,nin)
      if(ovlp)call Levet_4(ia(i_fmap0+nnof1),nnof2,nin)
c ......................................................................
c
c ...   Restricoes Globais
c
c      if(ndf .gt. 0) then
c         read(nin,*) string
c         i_idG = alloc_4('idG     ', ndf, nnoG)
c         call Lemtx_4(ia(i_idG),ndf,nin)
c      endif
c      if (ndft .gt. 0) then
c         read(nin,*) string
c         i_idGt = alloc_4('idGt    ', ndft, nnoG)
c         call Lemtx_4(ia(i_idGt),ndft,nin)    
c      endif
c ......................................................................
c
c ... Estruturas de comunicacao sendr
c
c      read(nin,*) string !nviz
c      read(nin,*) nviz
c      i_ranks  = alloc_4('ranks   ', 1, nprcs)
c      read(nin,*) string !ranks
c      read(nin,*) (ia(i_ranks+j), j = 0,nviz-1)
      nviz1 = 0
      nviz2 = 0
      read(nin,*) string !nviz
      read(nin,*) nviz1
      nviz2 = nviz1
      if(ovlp)read(nin,*) nviz2
c ....
      i_ranks  = alloc_4('ranks   ', 1, nviz1+nviz2)
      read(nin,*) string !ranks
      read(nin,*) (ia(i_ranks+j), j = 0,nviz1-1)
      if(ovlp)read(nin,*) (ia(i_ranks+nviz1+j), j = 0,nviz2-1)
c ....
      i_rcvs0 = alloc_4('rcvs0   ', 1, nviz1+nviz2)
      read(nin,*) string !rcvs
      read(nin,*) (ia(i_rcvs0+j), j = 0,nviz1-1)
      if(ovlp)read(nin,*) (ia(i_rcvs0+nviz1+j), j = 0,nviz2-1)
c ....
      i_dspl0 = alloc_4('dspl0   ', 1, nviz1+nviz2)
      read(nin,*) string !dspl
      read(nin,*) (ia(i_dspl0+j), j = 0,nviz1-1)
      if(ovlp)read(nin,*) (ia(i_dspl0+nviz1+j), j = 0,nviz2-1)
c ... end parallel
      read(nin,*) string
c ......................................................................
      return
      end
c ......................................................................
      subroutine Levet_4(array,lin,nin)
c **********************************************************************
c *                                                                    *
c *   Levet_4                                                          *
c *   -------                                                          *
c *                                                                    *
c *   Le arranjo inteiro                                               *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de entrada:                                          *
c *   ----------------------                                           *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de saida:                                            *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer array(*),lin,nin,i
c ......................................................................
      do i=1,lin
        read(nin,*)array(i)
      enddo
c ......................................................................
      return
      end
c ......................................................................
      subroutine Lemtx_4(array,lin,nin)
c **********************************************************************
c *                                                                    *
c *   Lemtx_4                                                          *
c *   -------                                                          *
c *                                                                    *
c *   Le arranjo inteiro                                               *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de entrada:                                          *
c *   ----------------------                                           *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de saida:                                            *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer array(lin,*),lin,nin,i,k
      character*80 comando
c ......................................................................
      read(nin,'(a80)') comando
      do while(comando(1:3) .ne. 'end')
        read(comando,*) k,(array(i,k),i=1,lin)
        read(nin,'(a80)') comando
      enddo
c ......................................................................
      return
      end
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco : 09/04/2016                                    * 
c * ------------------------------------------------------------------ *
c * MK_BOUND_QUAD: gera as restricoes nos deslocamente nos pontos      *  
c *              intermediarios                                        *
c * ------------------------------------------------------------------ *
c * Par�metros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * id(ndf,*)  - restricoes nos vertices                               *
c * numel      - numero de elementos                                   *
c * ndf1       - grau de liberdade liberdade de deslocamento           *
c * ndf2       - grau de liberdade liberdade total do problema         *
c * nen        - numero de nos elementos quadraticos                   *
c * ------------------------------------------------------------------ *
c * Par�metros de saida:                                               *
c * ------------------------------------------------------------------ *
c * id(ndf1,*)  - restricoes atualizadas                               *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * mecanico     : ndf1 = x, y e z  |  ndf2 = x, y e z                 *  
c * poromecanico : ndf1 = x, y,e z  |  ndf2 = x, y, z e p              *
c **********************************************************************
      subroutine mk_bound_quad(id,el,numel,ndf1,ndf2,nen)
      implicit none
      integer maxEdge
      parameter (maxEdge = 12) 
c ...
      integer i,j,k
      integer el(nen+1,*)
      integer iEdge(3,maxEdge)
      integer numel,ndf1,ndf2,nedge,no1,no2,no3,nen
      integer id(ndf2,*)
c ...
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
      do i = 1, numel
        do j = 1, nedge
c ... no vertices
          no1     = el(iEdge(1,j),i)
          no2     = el(iEdge(2,j),i)
c ... no central
          no3     = el(iEdge(3,j),i)
          do k = 1, ndf1
            if( id(k,no1) .eq. 1 .and. id(k,no2) .eq. 1) then
              id(k,no3) = 1
            endif
          enddo
        enddo
      enddo
c .....................................................................
c
c ...
c     do i = 1, 27
c       print*,i,id(1:4,i)
c     enddo
c .....................................................................  
      return
      end
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco : 13/10/2016                                    * 
c * ------------------------------------------------------------------ *
c * MK_FORCES_QUAD: gera as valores das cargas  nos pontos             *
c * intermediarios                                                     *
c * ------------------------------------------------------------------ *
c * Par�metros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * id(ndf,*)  - restricoes atualizadas                                *
c * f (ndf,*)  - valor das cargas nos vertices                         *
c * el(nen+1,*)- conectividade nodal                                   *
c * numel      - numero de elementos                                   *
c * ndf1       - grau de liberdade liberdade de deslocamento           *
c * ndf2       - grau de liberdade liberdade total do problema         *
c * nen        - numero de nos elementos quadraticos                   *
c * ------------------------------------------------------------------ *
c * Par�metros de saida:                                               *
c * ------------------------------------------------------------------ *
c * id(ndf,*)  - restricoes atualizadas                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine mk_forces_quad(id,f,el,numel,ndf1,ndf2,nen)
      implicit none
      integer maxEdge
      parameter (maxEdge = 12) 
c ...
      integer i,j,k
      integer el(nen+1,*)
      integer iEdge(3,maxEdge)
      integer id(ndf2,*)
      integer numel,ndf1,ndf2,nedge,no1,no2,no3,nen
      real*8  f(ndf2,*)
c ...
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
      do i = 1, numel
        do j = 1, nedge
c ... no vertices
          no1     = el(iEdge(1,j),i)
          no2     = el(iEdge(2,j),i)
c ... no central
          no3     = el(iEdge(3,j),i)
          do k = 1, ndf1
            if( id(k,no1) .eq. 1 .and. id(k,no2) .eq. 1) then
              f(k,no3) = 0.5d0*(f(k,no1) +  f(k,no2))
            endif 
          enddo
        enddo
      enddo
c .....................................................................
c
c ...
c     do i = 1, 70
c       print*,i,f(1:4,i)
c     enddo
c .....................................................................  
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *
c * MK_INITIAL_QUAD:gera os valores iniciais nos pontos intermediarios *
c * ------------------------------------------------------------------ *
c * Par�metros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * f(ndf,*)   - valor das cargas nos vertices                         *
c * el(nen+1,*)- conectividade nodal                                   *
c * numel      - numero de elementos                                   *
c * ndf        - grau de liberdade                                     *
c * nen        - numero de nos elementos quadraticos                   *
c * ------------------------------------------------------------------ *
c * Par�metros de saida:                                               *
c * ------------------------------------------------------------------ *
c * f(ndf,*)   - valor das cargas atualizadas                          *
c * ------------------------------------------------------------------ * 
c *  OBS:                                                              *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine mk_initial_quad(f,el,numel,ndf,nen)
      implicit none
      integer maxEdge
      parameter (maxEdge = 12) 
c ...
      integer i,j,k
      integer el(nen+1,*)
      integer iEdge(3,maxEdge)
      integer numel,ndf,nedge,no1,no2,no3,nen
      real*8  f(ndf,*)
c ...
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
      do i = 1, numel
        do j = 1, nedge
c ... no vertices
          no1     = el(iEdge(1,j),i)
          no2     = el(iEdge(2,j),i)
c ... no central
          no3     = el(iEdge(3,j),i)
          do k = 1, ndf
            f(k,no3) = 0.5d0*(f(k,no1) +  f(k,no2))
          enddo
        enddo
      enddo
c .....................................................................
c
c ...
c     do i = 1, 70
c       print*,i,f(1:4,i)
c     enddo
c .....................................................................  
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 27/09/2016                                    *
c * Data de modificaco : 14/10/2016                                    * 
c * ------------------------------------------------------------------ *   
c * SET_PRINT_VTK : leitualeitura das configuracoes basicas de excucao *
c * ------------------------------------------------------------------ *
c * Parametros de entrada :                                            *
c * -----------------------------------------------------------------  *
c * fprint    - nao definido                                           *
c * nin       - arquivo de entrada                                     *
c * -----------------------------------------------------------------  *
c * Parametros de saida :                                              *
c * -----------------------------------------------------------------  *
c * fprint - deslocamento    (1)                                       *   
c *          stress          (2)                                       *
c *          quadratico      (3)                                       *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine set_print_vtk(fprint,nin)
      implicit none
      include 'string.fi'
      character*15 string,macro(3)
      logical fprint(*),fexit
      integer j,nmacro
      integer nin
      data nmacro /3/
      data macro/'quadratic      ','desloc        ','stress         '/
c ......................................................................
c
c ...
      do j = 1, nmacro
        fprint(j) = .false.
      enddo
c ......................................................................
c
c ...
      fexit = .true.
      call readmacro(nin,.false.)
      write(string,'(15a)') (word(j),j=1,12)
      do while (strl .ne. 0 )
c ... quadratic 
        if     (string .eq. macro(1)) then
          fprint(1) = .true.
c .....................................................................
c
c ... desloc   
        elseif (string .eq. macro(2)) then
          fprint(2) = .true. 
c .....................................................................
c
c ... stress   
        elseif (string .eq. macro(3)) then
          fprint(3) = .true.
        endif
c .....................................................................
        call readmacro(nin,.false.)
        write(string,'(15a)') (word(j),j=1,15)                 
      end do
c .....................................................................
c
c ...
      do j = 1, nmacro
        print*,macro(j),fprint(j)
      enddo
c ......................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine init_poro_mec(f,nnode,ndf,n1,n2,nin)
c **********************************************************************
c *                                                                    *
c *   Valores iniciais para o problema poro mecanico                   *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'      
      integer n1,n2,nnode,ndf,nin,i,j,k
      real*8 f(ndf,*)
      character*30 string            
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) k
         if(k .lt. 1 .or. k .gt. nnode) goto 200      
         do 10 j = n1, n2
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(i),i=1,30)
            read(string,*,err = 200,end = 200) f(j,k)
   10    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(i),i=1,12)
  100 continue
      return
c ......................................................................
  200 continue
      print*,'*** Error in reading the nodal forces !'
      stop             
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 13/11/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *
c * READ_GRAVITY : leitura do campo de gravidade                       *
c * ------------------------------------------------------------------ *
c * Parametros de entrada :                                            *
c * -----------------------------------------------------------------  *
c * lines     - nao definido                                           *
c * nin       - arquivo de entrada                                     *
c * -----------------------------------------------------------------  *
c * Parametros de saida :                                              *
c * -----------------------------------------------------------------  *
c * str    - macro gravity                                             *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * gravity gx gy gz                                                   * 
c **********************************************************************
      subroutine read_gravity_v2(str,nin)
      include 'string.fi'
      character*15 gx,gy,gz
      character*200 str,aux
      integer nin
c ... gx
      call readmacro(nin,.false.)
      write(gx,'(15a)') (word(i),i=1,15)
      write(gx,'(15a)')  adjustl(gx)
c .....................................................................
c
c ... gy
      call readmacro(nin,.false.)
      write(gy,'(15a)') (word(i),i=1,15)
      write(gy,'(15a)')  adjustl(gy)
c .....................................................................
c  
c ... gz
      call readmacro(nin,.false.)
      write(gz,'(15a)') (word(i),i=1,15)
      write(gz,'(15a)')  adjustl(gz)
c .....................................................................
c
c ...
      write( aux, '(1x,A,1x,A,1x,A)' ) ,trim(gx),trim(gy),trim(gz)
      write( str, '( 2A )' ) 'gravity',trim(aux)
      return
c .....................................................................
c
c ...
      end
c **********************************************************************
c
c *********************************************************************
c * Data de criacao    : 20/03/2017                                   *
c * Data de modificaco : 00/00/0000                                   *
c * ------------------------------------------------------------------*
c * file_prop : : leitura do arquivo com as propriedades do material  *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * etype     - tipo do elemento                                      *
c * e(*)      - propriedades fisicas do material                      *
c * nin       - arquivo de entrada                                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * OBS:                                                              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine file_prop(etype,e,nin)
      implicit none
      include 'string.fi'
      character*15 string,macro(15)
      character*80 fname
      logical etyp
      real*8 e(*)
      integer i,j,nmacro,etype
      integer nin
      data nmacro /15/
      data macro/'modE           ','poisson        ','permeability   ',
     1           'mbiot          ','cbiot          ','density        ',
     2           'fdensity       ','ivoid          ','l_plastic      ',
     3           'k_plastic      ','mcs            ','pc0            ',
     4           'elType         ','               ','               '/
c .....................................................................
c
c ...
      call readmacro(nin,.true.)
      write(string,'(15a)') (word(j),j=1,12)
      do while (string .ne. 'end')
c ... modulo de elasticidade
         if (string .eq. macro(1)) then
            i = 1
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(1)
c .....................................................................
c
c ... poisson
         else if (string .eq. macro(2)) then
            i = 2
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(2)
c .....................................................................
c
c ... permeability
         else if (string .eq. macro(3)) then
            i = 3
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(3)
c .....................................................................
c
c ... mbiot        
         else if (string .eq. macro(4)) then
            i = 4
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(4)
c .....................................................................
c
c ... cbiot        
         else if (string .eq. macro(5)) then
            i = 5
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(5)
c .....................................................................
c
c ... density      
         else if (string .eq. macro(6)) then
            i = 6
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(6)
c .....................................................................
c
c ... fdensity      
         else if (string .eq. macro(7)) then
            i = 7
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(7)
c .....................................................................
c
c ... ivoid      
         else if (string .eq. macro(8)) then
            i = 8
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(8)
c .....................................................................
c
c ... l_plastic      
         else if (string .eq. macro(9)) then
            i = 9
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(9)
c .....................................................................
c
c ... k_plastic      
         else if (string .eq. macro(10)) then
            i = 10
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(10)
c .....................................................................
c
c ... mcs      
         else if (string .eq. macro(11)) then
            i = 10
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(11)
c .....................................................................
c
c ... pc0      
         else if (string .eq. macro(12)) then
            i = 12
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(12)
c .....................................................................
c
c ... elType   
         else if (string .eq. macro(13)) then
            i = 13
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) etype
c .....................................................................
         endif 
c .....................................................................
c
c ... 
         call readmacro(nin,.true.)
         write(string,'(15a)') (word(j),j=1,15)
      end do
c ......................................................................
c
c ...
      return  
c ......................................................................
 100  continue
      print*,'*** Erro in reading the matprop file ! ',macro(i)
      call finalize()                         
 200  continue
      print*,'File ',trim(fname),' not found !'
      call finalize()      
c ......................................................................
      end
c *********************************************************************
