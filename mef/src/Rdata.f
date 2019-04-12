      subroutine rdat_pm(nnode     ,nnodev    ,numel      ,numat
     1                  ,nen       ,nenv      ,ntn        ,ndf
     2                  ,ndm       ,nst       ,i_ix       ,i_ie 
     3                  ,i_inum    ,i_e       ,i_x        ,i_id 
     4                  ,i_nload   ,i_eload   ,i_eloadp   ,i_f
     5                  ,i_u       ,i_u0      ,i_tx0      ,i_dp
     6                  ,i_v
     6                  ,i_tx1p    ,i_tx2p    ,i_epsp     ,i_plastic
     7                  ,i_porosity,i_fnno    ,i_elplastic,i_vpropel
     8                  ,fstress0  ,fporomec  ,fmec      ,fterm
     9                  ,print_quad,plastic   ,vprop     ,nin     )
c **********************************************************************
c * Data de criacao    : 10/01/2016                                    *
c * Data de modificaco : 12/04/2019                                    *
c * ------------------------------------------------------------------ *
c * RDAT: leitura de dados do problema poromecanico.                   *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * nin     - arquivo de entrada                                       *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * nnode - numero total de nos                                        *
c * nnodev- numero de nos dos vertices                                 *
c * numel - numero de elementos                                        *
c * numat - numero de materiais                                        *
c * nen   - numero max. de nos por elemento                            *
c * nenv  - numero max. de nos geometicos por elemento                 *
c * ntn   - numero tensoes no tensor de tensao (ntn = 4 2D; ntn = 6 3D)*
c * ndf   - numero max. de graus de liberdade por no                   *
c * ndm   - dimensao (1, 2 ou 3)                                       *
c * nst   - numero de graus de liberdade por elemento                  *
c * i_ix    - ponteiro para conetividades                              *
c * i_id    - ponteiro para restricoes nodais (poro_mecanico)          *
c * i_ie    - ponteiro para materiais                                  *
c * i_nload - ponteiro para o arranjo nload (poro_mecanico)            *
c * i_eload - ponteiro para o arranjo eload (poro_mecanico-mecanico)   *
c * i_eloadp- ponteiro para o arranjo eload (poro_mecanico-hidraulico) *
c * i_inum  - ponteiro para o arranjo inum                             *
c * i_e     - ponteiro para o arranjo e                                *
c * i_x     - ponteiro para o arranjo x                                *
c * i_f     - ponteiro para o arranjo f (poro_mecanico)                *
c * i_u     - ponteiro para o arranjo u (poro_mecanico)                *
c * i_u0    - ponteiro para o arranjo u0(poro_mecanico)                *
c * i_v     - ponteiro para o arranjo v (termico)                      *
c * i_tx0   - ponteiro para o arranjo tx0(poro_mecanico)               *
c * i_dp    - ponteiro para o arranjo deltaP(poro_mecanico)            *
c * i_txp1  - ponteiro para o arranjo txp(poro_mecanico)               *
c * i_txp2  - ponteiro para o arranjo txp(poro_mecanico)               *
c * i_epsp  - ponteiro para o arranjo espp(poro_mecanico)              *
c * i_plastic   - ponteiro para o arranjo espp(poro_mecanico)          *
c * i_porosity  - ponteiro para o arranjo espp(poro_mecanico)          *
c * i_fnno      - ponteiro para o arranjo fnno                         *
c * i_elplastic - ponteiro para o arranjo elplastic                    *
c * i_vporel    - ponteiro para o arranjo vpropel                      *
c * fstress0- leitura de tensoes iniciais (true/false)                 *
c * fpropmec- (true/false)                                             *
c * fmec    - (true/false)                                             *
c * fterm   - (true/false)                                             *
c * fstress0- leitura de tensoes iniciais (true/false)                 *
c * print_quad - escrita da malha com elmentos quadraticos(true|false) *
c * plastic  - (true/false)                                            * 
c * vprop(*) - propriedades variaveis                                  *
c *           1 - prop variavel  (true|false)                          *
c *           2 - konzey-Caraman (true|false)                          *
c *           3 - mecanico       (true|false)                          *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * ix    - conetividades nodais dos elementos                         *
c * id    - restricoes nodais                                          *
c * ie    - tipo de elemento                                           *
c * e     - constantes fisicas                                         *
c * x     - coordenadas nodais                                         *
c * f     - forcas e prescricoes nodais                                *
c * u0    - condicoes de contorno e inicial                            *
c * v     - deriva primeria da solucao
c * tx0   - tensoes inicias                                            *
c * dp    - delta P ( apenas alocado)                                  *
c * txp   - tensoes nos pontos de integracao                           *
c * epsp  - delta deformacao e pressao entre iteracoes nao lineares    *
c *         deformacao volumetrica plastica total e parametro          *
c *         de encruamento nos pontos de integracao                    *
c * fnno  - identifica dos nos de vertices ( 1 - vertice | 0 )         *
c * elplastic - identificao se o elemento plastificou ou nao (0 ou 1)  *
c * vpropel   - propriedades variaveis por ponto de integracao         *
c *                                                                    *
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
      integer nnodev,nnode,numel,numat,nen,nenv,ndf,ndft,ndm,nst,ntn
      integer maxgrade,iplastic
c ... ponteiros      
      integer*8 i_e,i_x,i_f,i_nload,i_eload,i_eloadp,i_inum
      integer*8 i_u,i_u0,i_tx0,i_tx1p,i_tx2p,i_epsp,i_plastic,i_dp
      integer*8 i_ix,i_id,i_ie,i_porosity,i_elplastic,i_vpropel
      integer*8 i_nelcon,i_nodcon,i_nincid,i_incid,i_fnno,i_aux,i_v
c ......................................................................
      integer nin
      integer j,nmc,totnel,nmacros,itmp
      character*15 macro(39),rc
      character*80 fname
      integer naux
      integer nincl /7/
      logical fstress0,fporomec,fmec,fterm,plastic,print_quad,vprop(*)
      logical f_read_el /.false./
      logical el_quad  /.false./
      logical mk_el_quad  /.false./
      logical fMixedMesh  /.false./
c ......................................................................
      data macro/'materials      ','bar2           ','tria3          ',
     1           'quad4          ','tetra4         ','hexa8          ',
     2           'hexa20         ','tetra10        ','               ',
     3           'coordinates    ','constrainpmec  ','constraindisp  ',
     4           'nodalforces    ','elmtloads      ','nodalloads     ',
     5           'constraintemp'  ,'nodalthermloads','               ',
     6           'hydrostatic    ','hydrostress    ','elmtpresloads  ',
     7           'loads          ','               ','               ',
     8           'initialdisp    ','initialpres    ','initialstress  ',
     9           'parallel       ','insert         ','return         ',
     1           'tria3ov        ','quad4ov        ','tetra4ov       ',
     2           'hexa8ov        ','tetra10ov      ','hexa20ov       ',
     3           '               ','fmaterials     ','end            '/
      data nmc /39/      
c ......................................................................
c
c ... Leitura dos parametros da malha: nnode,numel,numat,nen,ndf,ndm
      if(my_id .eq. 0) print*,'loading parameters ...'
      call parameters(nnodev,numel,numat,nen,ndf,ndm,nin)
      if(my_id .eq. 0) print*,'load.'
      nnode  = nnodev
      nst    = 0
c ......................................................................
c
c ... tipo do problema
      fmec     = .false.
      fporomec = .false.
      fterm    = .false.
c ......................................................................
c
c ...
      i_tx1p      = 1 
      i_tx2p      = 1 
      i_epsp      = 1
      i_plastic   = 1
      i_elplastic = 1
c ... temporario
      if( nen .eq. 10 )then
        npi = 4
      else if( nen .eq. 20 )then
        npi = 64
      endif       
c .....................................................................
c
c ... problema termico
      if( ndf .eq. 1) fterm = .true.
c .....................................................................
c
c ... 
      if(ndm .eq. 2) then
c ... mecancico 
        if(ndf .eq. 2) then
          fmec = .true.
c ... poro mecanico
        else if(ndf .eq. 3) then
          fporomec = .true.
        endif
c .....................................................................
c
c ...
      elseif(ndm .eq. 3) then
c ... mecancico 
        if( ndf .eq. 3) then
          fmec = .true.
c ... poro mecanico
        else if(ndf .eq. 4) then
          fporomec = .true.
        endif
      endif
c ......................................................................
c
c ...
      nbar2(1:4)    = 0
      ntria3(1:4)   = 0
      ntetra4(1:4)  = 0
      nhexa8(1:4)   = 0
      ntetra10(1:4) = 0
      nhexa20(1:4)  = 0
c ......................................................................
c
c ... numero do tensor de tensoes
c ... | sxx syy szz sxy  0 0 0|
      if( ndm .eq. 2) then
        ntn = 4
c ... | sxx syy szz  sxy syz sxz |
      else if(ndm .eq. 3) then
        ntn = 6
      endif
c ......................................................................
c   
c     Alocacao de arranjos na memoria: mecancio
c     ---------------------------------------------------------------
c     | ix | ie | e | x | eload |
c     ---------------------------------------------------------------
c
      if(fmec) then
        i_ix    = alloc_4('ix      ',nen+1,numel)
        i_ie    = alloc_4('ie      ',    1,numat)
        i_eload = alloc_4('eload   ',    7,numel)
        i_e     = alloc_8('e       ', prop,numat)
        i_x     = alloc_8('x       ',  ndm,nnode)  
        call mzero(ia(i_ix),numel*(nen+1))
        call mzero(ia(i_ie),numat) 
        call mzero(ia(i_eload),numel*7)
        call azero(ia(i_e),numat*prop)
        call azero(ia(i_x),nnodev*ndm)        
      endif 
c ......................................................................
c   
c     Alocacao de arranjos na memoria: poromecanico
c     elastic
c     ---------------------------------------------------------------
c     | ix | ie | e | x | eload |
c     ---------------------------------------------------------------
c     plastic
c     ---------------------------------------------------------------
c     | ix | ie | e | x | eload | txp1 | txp2 | espsp | plastcic |
c     ---------------------------------------------------------------
c     proppriedades variaveis com a porosidade
c     ---------------------------------------------------------------
c     | vpropel |                                                    
c     ---------------------------------------------------------------
      if(fporomec) then
        i_ix    = alloc_4('ix      ',nen+1,numel)
        i_ie    = alloc_4('ie      ',    1,numat)
        i_eload = alloc_4('eload   ',    7,numel)
        i_eloadp= alloc_4('eloadp  ',    7,numel)
        i_e     = alloc_8('e       ', prop,numat)
        i_x     = alloc_8('x       ',  ndm,nnodev)
        call mzero(ia(i_ix),numel*(nen+1))
        call mzero(ia(i_ie),numat) 
        call azero(ia(i_e),numat*prop)
        call mzero(ia(i_eload),numel*7)
        call mzero(ia(i_eloadp),numel*7)
        call azero(ia(i_x),nnodev*ndm)        
        if(plastic) then  
          i_elplastic= alloc_4('elplast ',      1,numel)
          i_tx1p     = alloc_8('tx1p    ',ntn*npi,numel)
          i_tx2p     = alloc_8('tx2p    ',ntn*npi,numel)
          i_epsp     = alloc_8('epsp    ',7*npi,numel)
          i_plastic  = alloc_8('vplastic',3*npi,numel)
          call mzero(ia(i_elplastic),       numel)
          call azero(ia(i_tx1p)     ,ntn*npi*numel)
          call azero(ia(i_tx2p)     ,ntn*npi*numel)
          call azero(ia(i_epsp)     ,7*npi*numel)
          call azero(ia(i_plastic)  ,3*npi*numel)
        endif
        i_vpropel = 1
        if(vprop(1))then
          i_vpropel  = alloc_8('vpropel ',nvprop*npi,numel)
          call azero(ia(i_vpropel)       ,nvprop*npi*numel)
        endif
      endif 
c ......................................................................
c
c ......................................................................
c   
c     Alocacao de arranjos na memoria: termico
c     ---------------------------------------------------------------
c     | ix | ie | e | x | eload |
c     ---------------------------------------------------------------
c
      if(fmec .or. fterm) then
        i_ix    = alloc_4('ix      ',nen+1,numel)
        i_ie    = alloc_4('ie      ',    1,numat)
        i_eload = alloc_4('eload   ',    7,numel)
        i_e     = alloc_8('e       ', prop,numat)
        i_x     = alloc_8('x       ',  ndm,nnode)  
        call mzero(ia(i_ix),numel*(nen+1))
        call mzero(ia(i_ie),numat) 
        call mzero(ia(i_eload),numel*7)
        call azero(ia(i_e),numat*prop)
        call azero(ia(i_x),nnodev*ndm)        
      endif 
c ......................................................................
      totnel  = 0            
      nmacros = 0
c ......................................................................
c ...
c     Alocacao de arranjos na memoria:
c     ---------------------------------------------------------------
c     | id  nload | inum | e | x | f | u | u0 | tx0 | | v |
c     ---------------------------------------------------------------
   50 continue
      if (totnel .eq. numel) then
        i_inum  = alloc_4('inum    ',    1,nnode)  
        i_id    = alloc_4('id      ',  ndf,nnode)
        i_nload = alloc_4('nload   ',  ndf,nnode)
        i_f     = alloc_8('f       ',  ndf,nnode)
        i_u     = alloc_8('u       ',  ndf,nnode)
        i_u0    = alloc_8('u0      ',  ndf,nnode)
        i_tx0   = alloc_8('tx0     ',  ntn,nnode)
        call mzero(ia(i_inum) ,nnode)  
        call mzero(ia(i_id)   ,nnode*ndf)
        call mzero(ia(i_nload),nnode*ndf)
        call azero(ia(i_f)    ,nnode*ndf)
        call azero(ia(i_u)    ,nnode*ndf)
        call azero(ia(i_u0)   ,nnode*ndf)
        call azero(ia(i_tx0)  ,nnode*ntn)
c ...
        if(fporomec) then
          i_dp       = alloc_8('dpres   ',    1,nnode)
          i_porosity = alloc_8('poro    ',    1,nnode)
          call azero(ia(i_dp)      ,nnode)
          call azero(ia(i_porosity),nnode)
        endif    
c .....................................................................
c
c ... term
        if(fterm) then
          i_v   = alloc_8('v       ',  ndf,nnode)
          call mzero(ia(i_v) ,nnode*ndf)      
        endif
c .....................................................................
      endif
c .....................................................................
c
c ...               
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
      go to (400, 450, 500,    ! eaterials   ,bar2         ,tria3
     1       550, 600, 650,    !quad4        ,tetra4       ,hexa8
     2       700, 750, 800,    !hexa20       ,tetra10      ,
     3       850, 900, 950,    !coordinates  ,             ,constraindisp
     4      1000,1050,1100,    !nodalforces  ,elmtloads    ,nodalloads
     4      1150,1200,1250,    !constraintemp,nodalthermload
     5      1300,1350,1400,    !hstaticpres  ,hstaticstress,elmtpresloads
     6      1450,1500,1550,    !loads        ,             ,
     7      1600,1650,1700,    !initialdisp  ,             ,initialstress
     8      1750,1800,1850,    !parallel     ,insert       ,return
     9      1900,1950,2000,    !tria3ov      ,quad4ov      ,tetra4ov
     1      2050,2100,2150,    !hexa8ov      ,tetra10ov    ,hexa20ov
     2      2200,2250,2300) j  !             ,fmaterials   ,end
c ......................................................................
c
c ... Propriedades dos materiais:
c
  400 continue
c     if(my_id .eq. 0) print*,'loading materials ...'
      call mate(ia(i_ie),ia(i_e),numat,nin)
      call check_element(ia(i_ie),numat,plastic,vprop)
c     if(my_id .eq. 0) print*,'done.'
      go to 100
c ......................................................................
c
c ... Conetividades bar2:    
c
  450 continue
      nbar2(1) = 0
      call elconn(ia(i_ix),nen+1,2,nbar2(1),numel,fMixedMesh,nin)
      nbar2(2) = totnel+1
      totnel   = totnel + nbar2(1)      
      go to 100
c ......................................................................
c
c ... Conetividades tria3:          
c
  500 continue
      if(my_id .eq. 0) print*,'loading tria3 ...'
      f_read_el = .true.
      ntria3(1) = 0
      nenv      = 3
      nst       = nenv*ndf
      call elconn(ia(i_ix),nen+1,3,ntria3(1),numel,fMixedMesh,nin)
      ntria3(2) = totnel+1
      totnel    = totnel + ntria3(1)
c ... transforma os elementos lineares em quadraticos (3 nos)
c     if( nen .eq. 3) then
c     endif
c ...
c .....................................................................
      if(my_id .eq. 0) print*,'done.'
      go to 50
c ......................................................................
c
c ... Conetividades quad4:
c
  550 continue
      if(my_id .eq. 0) print*,'loading quad4 ...'
      f_read_el = .true.
      nquad4(1) = 0
      nenv      = 4
      nst       = nenv*ndf
      call elconn(ia(i_ix),nen+1,4,nquad4(1),numel,fMixedMesh,nin)
      nquad4(2) = totnel+1
      totnel    = totnel + nquad4(1)
c ... transforma os elementos lineares em quadraticos (8 nos)
c     if( nen .eq. 8) then
c     endif
c .....................................................................
      if(my_id .eq. 0) print*,'done.'
      go to 50  
c ......................................................................
c
c ... Conetividades tetra4:
c
  600 continue
      if(my_id .eq. 0) print*,'loading tetra4 ...'
      f_read_el  = .true.
      ntetra4(1) = 0
      nenv       = 4
      call elconn(ia(i_ix),nen+1,nenv,ntetra4(1),numel,fMixedMesh,nin)
      ntetra4(2) = totnel+1
      totnel     = totnel + ntetra4(1)
      nst        = nenv*ndf
c ... transforma os elementos lineares em quadraticos (10 nos)
      if( nen .eq. 10) then
        if(fmec) then
          itmp =  nen*ndf 
        else if(fporomec) then
          itmp = nen*(ndf-1) + nenv   
        endif
        nst   = max(nst,itmp)
        el_quad      = .true.
        mk_el_quad   = .true.
c .....................................................................
c
c ...  Multicore finite element assembling:
        i_nincid = alloc_4('nincid  ',1,nnodev) 
c ... Compute the maxgrade of the mesh and element incidences:
        call nodegrade(ia(i_ix),nnodev,numel,nenv,nen,ia(i_nincid)
     .                ,maxgrade) 
        i_incid  = alloc_4('incid   ',maxgrade,nnode)
        call elmincid(ia(i_ix),ia(i_incid),ia(i_nincid),nnodev,numel
     .               ,nenv    ,nen        ,maxgrade)
c ... gera a conectividade dos elementos quadraticos
        call mk_elconn_quad_v1(ia(i_ix),ia(i_incid),ia(i_nincid)
     .                        ,numel     ,nnode      ,nnodev
     .                        ,nen       ,nenv       ,maxgrade)
c .....................................................................
c
c ...
        i_incid     = dealloc('incid   ')
        i_nincid    = dealloc('nincid  ')
c .....................................................................
      endif
c .....................................................................
      if(my_id .eq. 0) print*,'done.'
      go to 50
c ......................................................................
c
c ... Conetividades hexa8:
c
  650 continue
      if(my_id .eq. 0) print*,'loading hexa8 ...'
      f_read_el = .true.
      nhexa8(1) = 0
      nenv      = 8
      call elconn(ia(i_ix),nen+1,nenv,nhexa8(1),numel,fMixedMesh,nin)
      nhexa8(2) = totnel + 1
      totnel    = totnel + nhexa8(1)
      nst       = nenv*ndf
c ... transforma os elementos lineares em quadraticos (20 nos)
      if(nen .eq. 20) then
        if(fmec) then
          itmp =  nen*ndf 
        else if(fporomec) then
          itmp = nen*(ndf-1) + nenv   
        endif
        nst   = max(nst,itmp) 
        el_quad      = .true.
        mk_el_quad   = .true.
c .....................................................................
c
c ...  Multicore finite element assembling:
        i_nincid = alloc_4('nincid  ',1,nnodev) 
c ... Compute the maxgrade of the mesh and element incidences:
        call nodegrade(ia(i_ix),nnodev,numel,nenv,nen,ia(i_nincid)
     .                ,maxgrade) 
        i_incid  = alloc_4('incid   ',maxgrade,nnode)
        call elmincid(ia(i_ix),ia(i_incid),ia(i_nincid),nnodev,numel
     .               ,nenv    ,nen        ,maxgrade)
c ... gera a conectividade dos elementos quadraticos
        call mk_elconn_quad_v1(ia(i_ix),ia(i_incid),ia(i_nincid)
     .                        ,numel   ,nnode      ,nnodev
     .                        ,nen     ,nenv       ,maxgrade)
c .....................................................................
c
c ...
        i_incid     = dealloc('incid   ')
        i_nincid    = dealloc('nincid  ')
c .....................................................................
      endif
c .....................................................................
      if(my_id .eq. 0) print*,'done.'
      go to 50
c ......................................................................
c
c ... Conetividades hexa20:
c
  700 continue
      if(my_id .eq. 0) print*,'loading hexa20 ...'
      f_read_el  = .true.
      el_quad    = .true.
      nhexa20(1) = 0
      nenv       = 8
      if(fmec) then
        itmp =  nen*ndf 
      else if(fporomec) then
        itmp = nen*(ndf-1) + nenv   
      endif
      nst        = max(nst,itmp) 
      call elconn(ia(i_ix),nen+1,20,nhexa20(1),numel,fMixedMesh,nin)
      nhexa20(2)  = totnel+1
      totnel      = totnel + nhexa20(1)
c ......................................................................
      if(my_id .eq. 0) print*,'done.'
      go to 50
c ......................................................................
c
c ... Conetividades tetra10:
c
  750 continue
      if(my_id .eq. 0) print*,'loading tetra10 ...'
      f_read_el   = .true.
      el_quad     = .true.
      ntetra10(1) = 0
      nenv        = 4
      if(fmec) then
        itmp =  nen*ndf 
      else if(fporomec) then
        itmp = nen*(ndf-1) + nenv   
      endif
      nst         = max(nst,itmp) 
      call elconn(ia(i_ix),nen+1,10,ntetra10(1),numel,fMixedMesh,nin)
      ntetra10(2)  = totnel+1
      totnel       = totnel + ntetra10(1)
c ......................................................................
      if(my_id .eq. 0) print*,'done.'
      go to 50
c ......................................................................
c
c ...
c      
  800 continue
      go to 100
c ......................................................................
c
c ... Coordenadas:
c
  850 continue
      if(my_id .eq. 0) print*,'loading coordinates ...'
      call coord(ia(i_x),nnodev,ndm,nin)
      if(my_id .eq. 0) print*,'done.'
      go to 100
c ......................................................................
c
c ... constrainpmec - restricoes nodais (deslocamentos + pressao)
c
  900 continue
      if(my_id .eq. 0) print*,'loading constrainpmec ...'
      if(f_read_el) then
        call bound(ia(i_id),nnodev,ndf,nin,1)
c ... malha quadratica gerada internamente
        if(mk_el_quad) then
          call mk_bound_quad(ia(i_id),ia(i_ix),numel,ndf-1,ndf,nen)
        endif
c ......................................................................
      else
        print*,'MACRO: constrainpmec !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      go to 100
c ......................................................................
c
c ... constraindisp - restricoes nodais (deslocamentos)
c
  950 continue
      if(my_id .eq. 0) print*,'loading constraindisp ...'
      if(f_read_el) then
        call bound(ia(i_id),nnodev,ndf,nin,1)
c ... malha quadratica gerada internamente
        if(mk_el_quad) then
          call mk_bound_quad(ia(i_id),ia(i_ix),numel,ndf,ndf,nen)
        endif
c ......................................................................
      else
        print*,'MACRO: constraindisp !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      go to 100
c ......................................................................
c
c ... nodalforces - forcas nodais:
c
 1000 continue
      if(f_read_el) then
        if(my_id .eq. 0) print*,'loading nodalforces ...'
        call forces(ia(i_f),nnodev,ndf,nin)
c ... malha quadratica gerada internamente
        if(mk_el_quad) then
          call mk_forces_quad(ia(i_id),ia(i_f),ia(i_ix),numel,ndf,nen)
        endif
c ......................................................................
      else
        print*,'MACRO: nodalforces !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      go to 100
c ......................................................................
c
c ... elmtloads - cargas nos elementos
c
 1050 continue
      if(my_id .eq. 0) print*,'loading elmtloads ...'
      if(f_read_el) then
        call bound(ia(i_eload),numel,7,nin,3) 
      else
        print*,'MACRO: elmtloads !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      go to 100
c ......................................................................
c
c ... nodalloads - nos com cargas variaveis no tempo:
c
 1100 continue
      if(my_id .eq. 0) print*,'loading nodalloads ...'
      if(f_read_el) then
        call bound(ia(i_nload),nnodev,ndf,nin,2) 
      else
        print*,'MACRO: nodalloads !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      goto 100
c ......................................................................
c
c ... nodalthermloads - nos com cargas variaveis no tempo:
c
 1150 continue
      if(my_id .eq. 0) print*,'loading constraintemp ...'
      if(f_read_el) then
        call bound(ia(i_id),nnodev,ndf,nin,1)
c ... malha quadratica gerada internamente
        if(mk_el_quad) then
          call mk_bound_quad(ia(i_id),ia(i_ix),numel,ndf,ndf,nen)
        endif
c ......................................................................
      else
        print*,'MACRO: constraindisp !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      goto 100
c ......................................................................
c
c ... nodalthermloads - nos com cargas variaveis no tempo:
c
 1200 continue
      if(my_id .eq. 0) print*,'loading nodalthermloads ...'
      if(f_read_el) then
        call bound(ia(i_nload),nnodev,ndf,nin,2) 
      else
        print*,'MACRO: nodalloads !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      goto 100
c ......................................................................
c
c ...
c
 1250 continue
      goto 100
c ......................................................................
c
c ...
c
 1300 continue
      if(my_id .eq. 0) print*,'loading hstaticpres ...'
      call init_hydrostatic_pres(ia(i_x),ia(i_u0)
     .                          ,nnodev,ndf,ndm,nin)
      if(my_id .eq. 0) print*,'done.'
      goto 100
c ......................................................................
c
c ...
c
 1350 continue
      if(my_id .eq. 0) print*,'loading hstaticstress ...'
      fstress0 = .true. 
      call init_hydrostatic_stress(ia(i_tx0),ia(i_x)
     .                             ,nnodev,ntn,ndm,nin)
      if(my_id .eq. 0) print*,'done.'
      go to 100
c ......................................................................
c
c ...
c
 1400 continue
      if(my_id .eq. 0) print*,'loading elmtpresloads ...'
      if(f_read_el) then
        call bound(ia(i_eloadp),numel,7,nin,3) 
      else
        print*,'MACRO: elmtpresloads !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      go to 100
c ......................................................................
c
c ... Definicao das cargas variaveis no tempo:
c
 1450 continue
      if(my_id .eq. 0) print*,'loading loads ...'
      call rload(nin)
      if(my_id .eq. 0) print*,'done.'
      goto 100
c ......................................................................
c
c ...
c
 1500 continue
      goto 100
c ......................................................................
c
c ...
c
 1550 continue
      goto 100
c ......................................................................
c
c ... initialdisp - deslocamentos iniciais:
c
 1600 continue
c     if(fReadEl) then
c       call init_poro_mec(ia(i_u0),nnodev,ndf,1,ndf-1,nin)  
c     else
c       print*,'MACRO: initialdisp !! elementos nao lidos'
c     endif
      go to 100
c ......................................................................
c
c ... intialpres
c
 1650 continue
      if(my_id .eq. 0) print*,'loading initialpres ...'
      if(f_read_el) then
        call init_poro_mec(ia(i_u0),nnodev,ndf,ndf,ndf,nin)
      else
        print*,'MACRO: initialpres !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      go to 100 
c ......................................................................
c
c ...
c      
 1700 continue
      if(my_id .eq. 0) print*,'loading initialstress ...'
      if(f_read_el) then
        fstress0 = .true.  
        call init_poro_mec(ia(i_tx0),nnodev,ntn,1,ntn,nin)
c ... malha quadratica gerada internamente
        if(mk_el_quad) then
          call mk_initial_quad(ia(i_tx0),ia(i_ix),numel,ntn,nen)
        endif
c ......................................................................
c
c ......................................................................
      else
        print*,'MACRO: initialstress !! Unread Elements'
        call stop_mef()
      endif
      if(my_id .eq. 0) print*,'done.'
      go to 100
c ......................................................................
c
c ... Paralelo:                                         
c      
 1750 continue
      if(my_id .eq. 0) print*,'loading read_par...'
      call read_par(nin,nnode,numel) 
      if(my_id .eq. 0) print*,'done.'
      go to 100 
c ......................................................................
c
c ... (insert) Desvia leitura para arquivo auxiliar:
c
 1800 continue
      nmacros = nmacros - 1
      naux = nin
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(nincl, file= fname,status= 'old',err=1851,action='read')
      nin = nincl
      go to 100
 1851 continue
      print*,'File ',trim(fname),' not found !'
      call stop_mef()
c ......................................................................
c
c ... (return) Retorna leitura para arquivo de dados basico:
c
 1850 continue
      nmacros = nmacros - 1
      close(nincl)
      nin = naux
      go to 100
c ......................................................................
c
c ... tria3ov
c
 1900 continue
      if(my_id .eq. 0) print*,'loading tria3ov ...'
      ntria3(3) = 0
      call elconn(ia(i_ix),nen+1,3,ntria3(3),numel,.true.,nin)
      ntria3(4) = totnel+1
      totnel    = totnel + ntria3(3)
      if(my_id .eq. 0) print*,'done.'
      go to 50
c ......................................................................  
c
c ... quad4ov
c
 1950 continue
      if(my_id .eq. 0) print*,'loading quad4ov ...'
      nquad4(3) = 0
      call elconn(ia(i_ix),nen+1,4,nquad4(3),numel,.true.,nin)
      nquad4(4) = totnel+1
      totnel     = totnel + nquad4(3)
      if(my_id .eq. 0) print*,'done.'
      go to 50
c ......................................................................
c
c ... tetra4ov                              
c      
 2000 continue
      if(my_id .eq. 0) print*,'loading tetra4ov ...'
      ntetra4(3) = 0
      call elconn(ia(i_ix),nen+1,4,ntetra4(3),numel,.true.,nin)
      ntetra4(4) = totnel+1
      totnel     = totnel + ntetra4(3)
      if(my_id .eq. 0) print*,'done.'
      go to 50                   
c ......................................................................
c
c ... hexa8ov                              
c      
 2050 continue
      if(my_id .eq. 0) print*,'loading hexa8ov ...'
      nhexa8(3) = 0
      call elconn(ia(i_ix),nen+1,8,nhexa8(3),numel,.true.,nin)
      nhexa8(4) = totnel+1
      totnel    = totnel + nhexa8(3)
      if(my_id .eq. 0) print*,'done.'
      go to 50                  
c ......................................................................
c
c ... tetra10ov                            
c      
 2100 continue
      if(my_id .eq. 0) print*,'loading tetra10ov ...'
      ntetra10(3) = 0
      call elconn(ia(i_ix),nen+1,10,ntetra10(3),numel,.true.,nin)
      ntetra10(4) = totnel+1
      totnel     = totnel + ntetra10(3)
      if(my_id .eq. 0) print*,'done.'
      go to 50                  
c ......................................................................
c
c ... hexa20ov                             
c      
 2150 continue
      if(my_id .eq. 0) print*,'loading hexa20ov ...'
      nhexa20(3) = 0
      call elconn(ia(i_ix),nen+1,20,nhexa20(3),numel,.true.,nin)
      nhexa20(4) = totnel+1
      totnel    = totnel + nhexa20(3)
      if(my_id .eq. 0) print*,'done.'
      go to 50                                
c ......................................................................
c
c ...
c      
 2200 continue
      go to 100
c ...................................................................... 
c
c ... fmaterials - arquivos com propriedades do materias
c      
 2250 continue
      if(my_id .eq. 0) print*,'loading fmaterials ...'
      call fmate(ia(i_ie),ia(i_e),numat,my_id,nin)
      call check_element(ia(i_ie),numat,plastic,vprop)
      if(my_id .eq. 0) print*,'done.'
      go to 100                  
c ......................................................................
c
c ... End:
c
 2300 continue
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
      i_fnno  = alloc_4('ffno    ',    1,nnode)
c ... identifica os nos de vertices( Necessario para o mpi)
      if(mk_el_quad) then
        call count_node_vertice(ia(i_ix),ia(i_fnno)
     .                         ,nnodev,nnode,numel,nenv,nen,.false.)
c ... identifica e conta os nos de vertices( Necessario para o mpi)
      else
        call count_node_vertice(ia(i_ix),ia(i_fnno)
     .                         ,nnodev,nnode,numel,nenv,nen,.true.)   
      endif
c ......................................................................
c
c
c ... Inicializa as condicoes de contorno no vetor u:
c
      if(ndf.gt.0) then 
c ... cc macro nodalsources
        call boundc(nnode,ndf,ia(i_id),ia(i_f),ia(i_u0))
c ... cc macro nodal lodas
        call bound_nodalloads(ia(i_x),ia(i_id),ia(i_u0)
     2                       ,ia(i_nload),nnode,ndf,ndm)
      endif  
      call aequalb(ia(i_u),ia(i_u0),nnode*ndf)
c .....................................................................
c
c ... escrita de todes os nohs
      if(el_quad) then
        if(print_quad) then
          if(mk_el_quad) then
            i_aux   = i_x
            i_x     = alloc_8('xq      ',  ndm,nnode)
c ... gerando as coordenada quadraticas 
            call mk_coor_quad(ia(i_aux) ,ia(i_x),ia(i_ix)
     .                     ,numel       ,nen  
     .                     ,nnode       ,nnodev  ,ndm)
c .....................................................................
c
c ...
            i_aux   = dealloc('x       ')
            i_x     =  locate('xq      ')
            i_inum  =  locate('inum    ')  
            i_id    =  locate('id      ')
            i_nload =  locate('nload   ')
            i_f     =  locate('f       ')
            i_u     =  locate('u       ')
            i_u0    =  locate('u0      ')
            i_tx0   =  locate('tx0     ')
            if(fporomec) then
              i_dp       = locate('dpres   ')
              i_porosity = locate('poros   ')
            endif
c ...
            i_ix    = locate('ix      ')
            i_ie    = locate('ie      ')
            i_eload = locate('eload   ')
            i_e     = locate('e       ')
            i_fnno  = locate('ffno    ')
            if(plastic) then  
              i_elplastic = locate('elplast ')
              i_tx1p      = locate('tx1p    ')
              i_tx2p      = locate('tx2p    ')
              i_epsp      = locate('epsp    ')
              i_plastic   = locate('plastic ')
            endif
            if(vprop(1)) then
              i_vpropel  = locate('vpropel ')
            endif
c .....................................................................
c
c ...
            ntetra10(1:4) = ntetra4(1:4)
            nhexa20(1:4)  = nhexa8(1:4)
            ntetra4(1:4)  = 0
            nhexa8(1:4)   = 0
          endif
c ......................................................................
c
c ... escrita apenas dos vertices
        else
          if(mk_el_quad .eqv. .false.) then
            ntetra4(1:4)  = ntetra10(1:4)
            nhexa8(1:4)   = nhexa20(1:4)
            ntetra10(1:4) = 0
            nhexa20(1:4)  = 0
          endif 
c ......................................................................
        endif
c ......................................................................
      endif  
c ......................................................................
      return
c ......................................................................
      end
      subroutine elconn(ix,nen1,nen,nel,numel,flag,nin)
c **********************************************************************
c *                                                                    *
c *   Conetividades nodais.                                            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'      
      integer ix(nen1,*),nen1,nen,nel,numel,nin,j,k,m
      character*12 string
      logical flag
c ......................................................................
      if(.not. flag) then
        do  50 j = 1, numel
          read(nin,*) k,(ix(m,k),m=1,nen),ix(nen1,k)
   50   continue
        call readmacro(nin,.true.)
        write(string,'(12a)') (word(j),j=1,12)  
        if (string .ne. 'end') goto 200
        nel=numel
        return
c ......................................................................

      else
        nel = 0
        call readmacro(nin,.true.)
        write(string,'(12a)') (word(j),j=1,12)
        do 100 while(string .ne. 'end')
          read(string,*,err = 200,end = 200) k
          if(k .lt. 1 .or. k .gt. numel) goto 200      
          do 10 j = 1, nen
             call readmacro(nin,.false.)
             write(string,'(12a)') (word(m),m=1,12)
             read(string,*,err = 200,end = 200) ix(j,k)
   10     continue
          call readmacro(nin,.false.)
          write(string,'(12a)') (word(m),m=1,12)
          read(string,*,err = 200,end = 200) ix(nen1,k)   
          nel = nel + 1
          call readmacro(nin,.true.)
          write(string,'(12a)') (word(j),j=1,12)
  100   continue
        return
      endif
c ......................................................................
  200 continue
      print*,'*** Erro na leitura dos elementos !',k,j
      stop                   
      end
c **********************************************************************
c
c **********************************************************************
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
      print*,'*** Erro na leitura das coordenadas !'
      stop             
      end
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
      real*8  e(prop,*)
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
c ... formato necessário ao funcionamento em linux:
         do 100 while ( (string .ne.   ' '    ) .and.
     .                  (string .ne. CHAR(13) ) )
            i = i + 1
            if(i .gt. prop) goto 200
            read(string,*,err = 200,end = 200) e(i,ma)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(m),m=1,30)
  100    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura dos materiais !'
      stop             
      end 
c **********************************************************************
c
c **********************************************************************
      subroutine fmate(ie,e,numat,my_id,nin)
c **********************************************************************
c *                                                                    *
c *   Materiais.                                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      include 'termprop.fi'
      character*80 fname
      integer ie(*),numat,nin,i,j,m,ma,my_id
      integer nincl /12/
      real*8  e(prop,*)
      character*30 string
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
c .....................................................
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) ma
         if(ma .lt. 1 .or. ma .gt. numat) goto 200
c .....................................................
c
c ... abrindo arquivo
         call readmacro(nin,.false.)
         write(fname,'(80a)') (word(j),j=1,strl)
         open(nincl, file= fname,status= 'old',err=900,action='read')
         call file_prop(ie(ma),e(1,ma),my_id,nincl)
         close(nincl)
c .....................................................................
c
c ...
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura dos materiais !'
      call stop_mef()
c ...
 900  continue
      if(my_id.eq.0) then
        print*,'File ',trim(fname),' not found !'
      endif
      call stop_mef()
c .....................................................................
      stop             
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
      stop             
      end
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
      print*,'*** Erro na leitura das forcas nodais !'
      stop             
      end
c **********************************************************************
c *
c **********************************************************************
      subroutine init_hydrostatic_pres(x,u,nnode,ndf,ndm,nin)
c **********************************************************************
c *                                                                    *
c *   Inicias nos com pressao hidroestatica                            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'gravity.fi'
      include 'string.fi'      
      integer ndm,nincl,nnode,ndf,nin,i
      parameter (nincl = 8)
      real*8 u(ndf,*),x(ndm,*),density,p0,h,scale
      character*12 string 
      character*80 fname
c ... nome do arquivo
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(i),i=1,strl)
      open(nincl, file= fname,status= 'old',err=201,action='read')
c ......................................................................
      call readmacro(nincl,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 200,end = 200) density
c ... pressao inicial      
      call readmacro(nincl,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 200,end = 200) p0
c ... datum      
      call readmacro(nincl,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 200,end = 200) h
c ... fator de escala      
      call readmacro(nincl,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 200,end = 200) scale
      do i = 1, nnode
        u(ndf,i) = ( p0 + density*gravity_mod*(h-x(3,i) ) )*scale
      enddo
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura da funcao init_hydrostatic_pres !'
      call stop_mef()
  201 continue
      print*,'File ',trim(fname),' not found !'
      call stop_mef()
      end
c **********************************************************************
c
c **********************************************************************
      subroutine init_hydrostatic_stress(tx,x,nnode,ntn,ndm,nin)
c **********************************************************************
c *                                                                    *
c *   Inicias nos com tensoes hidroestatica                            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'gravity.fi'
      include 'string.fi'      
      integer ndm,nincl,nnode,ntn,nin,i
      parameter (nincl = 8)
      real*8 tx(ntn,*),x(ndm,*),density,f0,k,h,scale
      character*12 string 
      character*80 fname 
c ......................................................................
c
c ... nome do arquivo
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(i),i=1,strl)
      open(nincl, file= fname,status= 'old',err=201,action='read')
c ... densidade
      call readmacro(nincl,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 200,end = 200) density
c ... tensao inicial      
      call readmacro(nincl,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 200,end = 200) f0
c ... datum      
      call readmacro(nincl,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 200,end = 200) h
c ...fator multiplicativo entre a tensao vertical e horizontal      
      call readmacro(nincl,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 200,end = 200) k
c ... fator de escala      
      call readmacro(nincl,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 200,end = 200) scale
      do i = 1, nnode
        tx(3,i)   = ( f0 + density*gravity_mod*(h-x(3,i))) *scale
        tx(1,i)   = k*tx(3,i)
        tx(2,i)   = k*tx(3,i)
        tx(4:6,i) = 0.d0
      enddo
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura da funcao init_hydrostatic_stress !'
      call stop_mef()
  201 continue
      print*,'File ',trim(fname),' not found !'
      call stop_mef()
      end
c **********************************************************************
c *
c **********************************************************************
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
         print*,'*** Erro na leitura das restricoes !',j,i
      elseif(code .eq. 2) then
         print*,'*** Erro na leitura das cargas nodais !',j,i
      elseif(code .eq. 3) then
         print*,'*** Erro na leitura das cargas nos elementos !',j,i
      endif
      stop
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
      stop             
      end
      subroutine rload(nin)
c **********************************************************************
c *                                                                    *
c *   Leitura das cargas variaveis no tempo                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
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
      integer  nin,nincl,i,j,k,l,nc,ty
      parameter (nincl = 8)
      character*30 string
      character*80 fname  
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
         ty =  load(1,i)
c .....................................................................
c
c ...   
         if     (ty .eq. 1) then
c ...       valor da carga:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(1,1,i)
c .....................................................................
c
c ...       u(t) = a.t + b 
         elseif (ty .eq. 2) then
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
c ...       -kdu/dx = -H.(u-uext) 
         elseif (ty .eq. 3) then
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
         elseif (ty .eq. 6) then
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
c ...       kdu/dx = emiss * const(Stef-Boltz) *(u4-uext4)+H(uext-u)
         elseif (ty .eq. 7) then
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
         elseif (ty .eq. 8) then
c ...       numero de parcelas:            
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
c 39 - valor de uma carga ou deslocamento variavel no tempo
c 40 - forca distribuida constante no contorno
c 41 - carga normal constante no contorno 
c      (F=carga*normal,normal - calculado nivel elemento)
c 42 - valor por carga hidrostatica      
c      (f = alfa*(f0 + density*g*(h-x)) )
c 43 - fluxo normal de massa
         elseif ( (ty .eq. 39) 
     .       .or. (ty .eq. 40) .or. (ty .eq. 41) 
     .       .or. (ty .eq. 42) .or. (ty .eq. 43) ) then
c ... nome do arquivo            
            call readmacro(nin,.false.)
            write(fname,'(80a)') (word(j),j=1,strl)
            open(nincl, file= fname,status= 'old',err=201,action='read')
c ... numero de parcelas
            call readmacro(nincl,.true.) 
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)    
c ... numero de parcelas temporais da carga
            call readmacro(nincl,.false.) 
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(3,i)
c ...
            do l = 1, load(3,i)
              call readmacro(nincl,.true.)
              do k = 1, load(2,i)                 
                 write(string,'(30a)') (word(j),j=1,30)
                 read(string,*,err = 200,end = 200) fload(l,k,i)
                 call readmacro(nincl,.false.)
              enddo
            enddo 
c .....................................................................
c
c ... 
            close(nincl)   
c .....................................................................
         endif
c .....................................................................
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  100 continue
      return
c ......................................................................
  200 continue
      print*,'*** Load reading error !',i
      call stop_mef()
  201 continue
      print*,'File ',trim(fname),' not found !'
      call stop_mef()
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
      stop
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
      subroutine parameters_pm(nnode,numel,numat,nen,ndf,ndm,plastic
     .                         ,nin)
c **********************************************************************
c *                                                                    *
c *   Parameters                                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      character*12 string
      integer nnode,numel,numat,nen,ndf,ndft,ndm,plastic,nin,n,j
      logical flag(7)
      character*8 macro(7)
      integer i,nmc         
      data macro/'nnode   ' ,'numel   ','numat   '
     .          ,'maxno   ' ,'ndf     ','dim     '
     .          ,'plastic '/
      data nmc /7/
c ......................................................................
      flag(1:nmc) = .false.
c ...
      plastic = 0
      flag(7) = .true.
c ......................................................................
      n       = 0
      call readmacro(nin,.true.)
      write(string,'(8a)') (word(j),j=1,8)
      do while (strl .ne. 0)
         if     (string .eq. 'nnode') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) nnode
            flag(1) = .true.
            n = n + 1
         elseif (string .eq. 'numel') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) numel
            flag(2) = .true.
            n = n + 1
         elseif (string .eq. 'numat') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) numat
            flag(3) = .true.
            n = n + 1
         elseif (string .eq. 'maxno') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) nen
            flag(4) = .true.
            n = n + 1
         elseif (string .eq. 'ndf') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) ndf
            flag(5) = .true.
            n = n + 1
         elseif (string .eq. 'dim') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) ndm
            flag(6) = .true.
            n = n + 1
         elseif (string .eq. 'plastic') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) plastic
            n = n + 1
         endif
         call readmacro(nin,.false.)
         write(string,'(8a)') (word(j),j=1,8)
      end do
c .....................................................................
c
c ...
      do i = 1, nmc
        if(.not. flag(i)) then 
          print *,"Missing macro: ",macro(i)
          call stop_mef()
        endif
      enddo
      return
c ......................................................................
  100 continue
      print*,'*** Erro in reading the control variables !'
      call stop_mef()    
c ......................................................................
      end
c **********************************************************************
c
c **********************************************************************
      subroutine parameters(nnode,numel,numat,nen,ndf,ndm,nin)
c **********************************************************************
c *                                                                    *
c *   Parameters                                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      character*12 string
      integer nnode,numel,numat,nen,ndf,ndft,ndm,nin,n,j
      logical flag(6)
      character*8 macro(6)
      integer i,nmc         
      data macro/'nnode   ' ,'numel   ','numat   '
     .          ,'maxno   ' ,'ndf     ','dim     '/
      data nmc /6/
c ......................................................................
      flag(1:nmc) = .false.
c ......................................................................
      n       = 0
      call readmacro(nin,.true.)
      write(string,'(8a)') (word(j),j=1,8)
      do while (strl .ne. 0)
         if     (string .eq. 'nnode') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) nnode
            flag(1) = .true.
            n = n + 1
         elseif (string .eq. 'numel') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) numel
            flag(2) = .true.
            n = n + 1
         elseif (string .eq. 'numat') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) numat
            flag(3) = .true.
            n = n + 1
         elseif (string .eq. 'maxno') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) nen
            flag(4) = .true.
            n = n + 1
         elseif (string .eq. 'ndf') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) ndf
            flag(5) = .true.
            n = n + 1
         elseif (string .eq. 'dim') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) ndm
            flag(6) = .true.
            n = n + 1
         endif
         call readmacro(nin,.false.)
         write(string,'(8a)') (word(j),j=1,8)
      end do
c .....................................................................
c
c ...
      do i = 1, nmc
        if(.not. flag(i)) then 
          print *,"Missing macro: ",macro(i)
          call stop_mef()
        endif
      enddo
      return
c ......................................................................
  100 continue
      print*,'*** Erro in reading the control variables !'
      call stop_mef()    
c ......................................................................
      end
c **********************************************************************
c
c **********************************************************************
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
      print*,'*** Error reading the data file! '
      stop             
c ......................................................................
      end 
c ======================================================================
c **********************************************************************
      subroutine readmacrov2(nin,newline,ier)
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
c *   ier     = codigo de erro    (  0 - sem erro                      *
c *                                  1 - com erro)                     *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer j,k,ier,nin
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
      ier = 0     
      return
c ......................................................................
  200 continue
      ier = 1
      return
c ......................................................................
      end 
c ======================================================================
c
c ======================================================================
c
c     Leitura da estrutura de dados do paralelo
c
c ======================================================================
      subroutine read_par(nin,nnode,numel)
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
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
c      include 'mpif.h'
      include 'parallel.fi'
      include 'elementos.fi'
      integer nin,nnode,numel,i,j,k
      character*500 comando,string
c ...................................................................... 
      read(nin,*) string,nnovG      
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
      read(comando,*) string, nbar2(1)   ,nbar2(2)   
     .               ,string, ntria3(1)  ,ntria3(2)  
     .               ,string, nquad4(1)  ,nquad4(2)  
     .               ,string, ntetra4(1) ,ntetra4(2) 
     .               ,string, nhexa8(1)  ,nhexa8(2)
     .               ,string, nprism6(1) ,nprism6(2)
     .               ,string, ntetra10(1),ntetra10(2)
     .               ,string, nhexa20(1) ,nhexa20(2)
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
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de saida:                                             *
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
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de saida:                                             *
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
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco : 09/04/2016                                    * 
c * ------------------------------------------------------------------ *
c * MK_BOUND_QUAD: gera as restricoes nos deslocamente nos pontos      *
c *              intermediarios                                        *
c * ------------------------------------------------------------------ *
c * Parâmetros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * id(ndf,*)  - restricoes nos vertices                               *
c * numel      - numero de elementos                                   *
c * ndf1       - grau de liberdade liberdade de deslocamento           *
c * ndf2       - grau de liberdade liberdade total do problema         *
c * nen        - numero de nos elementos quadraticos                   *
c * ------------------------------------------------------------------ *
c * Parâmetros de saida:                                               *
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
c * Data de modificaco :                                               *
c * ------------------------------------------------------------------ *
c * MK_FORCES_QUAD: gera as valores das cargas  nos pontos             *
c * intermediarios                                                     *
c * ------------------------------------------------------------------ *
c * Parâmetros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * id(ndf,*)  - restricoes atualizadas                                *
c * f (ndf,*)  - valor das cargas nos vertices                         *
c * el(nen+1,*)- conectividade nodal                                   *
c * numel      - numero de elementos                                   *
c * ndf        - grau de liberdade                                     *
c * nen        - numero de nos elementos quadraticos                   *
c * ------------------------------------------------------------------ *
c * Parâmetros de saida:                                               *
c * ------------------------------------------------------------------ *
c * id(ndf,*)  - restricoes atualizadas                                *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine mk_forces_quad(id,f,el,numel,ndf,nen)
      implicit none
      integer maxEdge
      parameter (maxEdge = 12) 
c ...
      integer i,j,k
      integer el(nen+1,*)
      integer iEdge(3,maxEdge)
      integer id(ndf,*)
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
          do k = 1, ndf-1
            if( id(k,no1) .eq. 1 .and. id(k,no2) .eq. 1) then
              f(k,no3) = 0.5d0*(f(k,no1) +  f(k,no2))
            endif 
          enddo
        enddo
      enddo
c .....................................................................
c
c ...
c     do i = 1, 44
c       print*,i,f(1:4,i)
c     enddo
c .....................................................................
      return
      end
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco :                                               *
c * ------------------------------------------------------------------ *
c * MK_INITIAL_QUAD:gera os valores iniciais nos pontos intermediarios *
c * ------------------------------------------------------------------ *
c * Parâmetros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * f(ndf,*)   - valor das cargas nos vertices                         *
c * el(nen+1,*)- conectividade nodal                                   *
c * numel      - numero de elementos                                   *
c * ndf        - grau de liberdade                                     *
c * nen        - numero de nos elementos quadraticos                   *
c * ------------------------------------------------------------------ *
c * Parâmetros de saida:                                               *
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
c *********************************************************************
c
c *********************************************************************
c * Data de criacao    : 20/03/2017                                   *
c * Data de modificaco : 02/04/2017                                   *
c * ------------------------------------------------------------------*
c * read_constitutive_equation : leitura do tipo regime das relacoes  *
c * constitutivas                                                     *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * iplastic  - plasticidade (true|false)                             *
c * ivprop    - propriedades variaveis                                *
c *           1 - prop por pontos de integracao (true|false)          * 
c *           2 - kozeny-Caraman                (true|false)          *
c *           3 - massa especifica              (true|false)          *                       
c *           4 - mecanico                      (true|false)          *
c * nin       - arquivo de entrada                                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * OBS:                                                              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine read_constitutive_equation(iplastic,ivprop,my_id,nin)
      implicit none
      include 'string.fi'
      character(len=16) string,macro(8)
      character(len=80) fname
      logical iplastic,ivprop(*)
      integer i,j,nmacro
      integer nin,nprcs
      integer my_id
      integer nincl /7/
      data nmacro /9/
      data macro/'plastic         ','vprop           ',
     1           'kozeny_carman   ','density         ',
     2           'hashin_shtrikman','                ',
     3           '                ','help            '/
c .....................................................................
c
c ...
      iplastic   = .false.
      ivprop(1:3)= .false.
c .....................................................................
c      
c ... arquivo de config
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(nincl, file= fname,status= 'old',err=200,action='read')
c .....................................................................
c
c ...
      call readmacro(nincl,.true.)
      write(string,'(16a)') (word(j),j=1,16)
      do while (string .ne. 'end')
c ... plastic
         if (string .eq. macro(1)) then
      
           iplastic = .true. 
c .....................................................................
c
c ... vprop
         else if (string .eq. macro(2)) then
           ivprop(1) = .true. 
c .....................................................................
c
c ... konzeycarman
         else if (string .eq. macro(3)) then
           ivprop(2) = .true. 
c .....................................................................
c
c ... density         
         else if (string .eq. macro(4)) then
           ivprop(3) = .true. 
c .....................................................................
c
c ... Hashin-Shtrikman
         else if (string .eq. macro(5)) then
           ivprop(4) = .true. 
c .....................................................................
c
c ... help            
         else if (string .eq. macro(8)) then
           goto 300            
c .....................................................................
         endif 
c .....................................................................
c
c ... 
         call readmacro(nincl,.true.)
         write(string,'(16a)') (word(j),j=1,16)
      end do
c ......................................................................
c
c ...
      close(nincl)
c ......................................................................
c
c ... caso variacao da permeabilidade, prop mecanicas ou massa
c     especifica habilitadas 
      if( ivprop(2) .or. ivprop(3) .or. ivprop(4) ) then
        if ( .not.  ivprop(1)) ivprop(1) = .true.
      endif
c ......................................................................
c
c ... 
      if(my_id .eq. 0 ) then
        print*,'Plastic         :',iplastic
        print*,'vprop           :',ivprop(1)
        print*,'Kozeny-Carman   :',ivprop(2)
        print*,'Density         :',ivprop(3)
        print*,'Hashin-Shtrikman:',ivprop(4)
      endif
c .....................................................................
c
c .....................................................................
      return  
c ......................................................................
 100  continue
      print*,'*** Erro in reading the matprop file !',macro(i)
      call stop_mef()                          
 200  continue
      print*,'File ',trim(fname),' not found !'
      call stop_mef()      
 300  continue
      if(my_id.eq.0) then
        write(*,'(2x,a)') '******************************'
        write(*,'(2x,a)') 'Example usage of macro conseq:'
        write(*,'(2x,a)') '------------------------------'
        do i = 1, 5
          write(*,'(2x,a)') macro(i)
        enddo  
        write(*,'(2x,a)') 'end'
        write(*,'(2x,a)') '------------------------------'
        write(*,'(2x,a)') '******************************'
      endif
      call stop_mef()
c ......................................................................
      end
c *********************************************************************
c
c *********************************************************************
c * Data de criacao    : 20/03/2017                                   *
c * Data de modificaco : 26/03/2017                                   *
c * ------------------------------------------------------------------*
c * file_prop : : leitura do arquivo com as propriedades do material  *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * etype     - tipo do elemento                                      *
c * e(*)      - propriedades fisicas do material                      *
c * my_id     - id do mpi                                             *
c * nin       - arquivo de entrada                                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * OBS:                                                              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine file_prop(etype,e,my_id,nin)
      implicit none
      include 'string.fi'
      character*15 string,macro(15)
      character*16 ex(14)
      character*80 fname
      logical etyp,fread(15)
      real*8 e(*)
      integer i,j,nmacro,etype
      integer nin
      integer my_id
      data nmacro /15/
      data macro/'modE           ','poisson        ','permeability   ',
     1           'mbiot          ','cbiot          ','density        ',
     2           'fdensity       ','porosity       ','l_plastic      ',
     3           'k_plastic      ','mcs            ','pc0            ',
     4           'elType         ','               ','               '/
c ... exemplo
      data ex/'elType        37',
     1        'modE         1.0',
     2        'poisson      0.3',
     3        'permeability 1.0',
     4        'mbiot        1.0',
     5        'cbiot        1.0',
     6        'density      1.0', 
     7        'fdensity     1.0',
     8        'porosity     0.72',
     9        'l_plastic    0.2',
     1        'k_plastic    0.2',
     2        'mcs          1.2',
     3        'pc0          1.0',
     4        'end             '/
c .....................................................................
c
c ...
      fread(1:15) = .false.
      call readmacro(nin,.true.)
      write(string,'(15a)') (word(j),j=1,12)
      do while (string .ne. 'end')
c ... modulo de elasticidade
         if (string .eq. macro(1)) then
            i = 1
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(1)
c .....................................................................
c
c ... poisson
         else if (string .eq. macro(2)) then
            i = 2
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(2)
c .....................................................................
c
c ... permeability
         else if (string .eq. macro(3)) then
            i = 3
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(3)
c .....................................................................
c
c ... mbiot        
         else if (string .eq. macro(4)) then
            i = 4
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(4)
c .....................................................................
c
c ... cbiot        
         else if (string .eq. macro(5)) then
            i = 5
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(5)
c .....................................................................
c
c ... density      
         else if (string .eq. macro(6)) then
            i = 6
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(6)
c .....................................................................
c
c ... fdensity      
         else if (string .eq. macro(7)) then
            i = 7
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(7)
c .....................................................................
c
c ... ivoid      
         else if (string .eq. macro(8)) then
            i = 8
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(8)
c .....................................................................
c
c ... l_plastic      
         else if (string .eq. macro(9)) then
            i = 9
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(9)
c .....................................................................
c
c ... k_plastic      
         else if (string .eq. macro(10)) then
            i = 10
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(10)
c .....................................................................
c
c ... mcs      
         else if (string .eq. macro(11)) then
            i = 11
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(11)
c .....................................................................
c
c ... pc0      
         else if (string .eq. macro(12)) then
            i = 12
            fread(i) = .true.
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) e(12)
c .....................................................................
c
c ... elType   
         else if (string .eq. macro(13)) then
            i = 13
            fread(i) = .true.
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
c ... elemento elasticos
      if( (etype .eq. 16) .or. (etype .eq. 17) ) then
        do i = 1, 7
          if(.not. fread(i)) then
            write(*,'(1x,a,a)') 'Property missing: ',trim(macro(i))
            call stop_mef()                    
          endif
        enddo
c ......................................................................
c
c ... elementos plasticos
      else if( (etype .eq. 36) .or. (etype .eq. 37) ) then
        do i = 1, 12
          if(.not. fread(i)) then
            write(*,'(1x,a,a)') 'Property missing: ',trim(macro(i))
            call stop_mef()      
          endif
        enddo
c ......................................................................
c
c ... 
      else
        if(.not. fread(13)) then
          write(*,'(1x,a,a)') 'Property missing: ',trim(macro(13))
          call stop_mef()      
        endif
      endif
c ......................................................................
c
c ...
      return  
c ......................................................................
 100  continue
      print*,'*** Erro in reading the matprop file ! ',macro(i)
      if(my_id.eq.0) then
        print*,'*** Error reading solver fmaterials macro !'
        write(*,'(2x,a)') '******************************'
        write(*,'(2x,a)') 'Example prop file usage :'
        write(*,'(2x,a)') '------------------------------'
        do i = 1, 14
          write(*,'(2x,a)') ex(i)
        enddo  
        write(*,'(2x,a)') '------------------------------'
        write(*,'(2x,a)') '******************************'
      endif
      call stop_mef()                          
 200  continue
      print*,'File ',trim(fname),' not found !'
      call stop_mef()      
c ......................................................................
      end
c *********************************************************************
c
c *********************************************************************
c * Data de criacao    : 00/00/0000                                   *
c * Data de modificaco : 27/05/2017                                   *
c * ------------------------------------------------------------------*
c * read_config : leitura das configuracoes basicas de excucao        *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * maxmem        - memoria do vetor de trabalho                      *
c * omp_elmt  - flag do openmp na fase de elemento                    *
c * nth_elmt  - numero de threads usado na fase de elemento           *
c * omp_solv  - flag do openmp na fase do solver                      *
c * nth_solv  - numero de threads usado do solver                     *
c * reord     - reordanaco do sistema de equacao                      *
c * newton_raphson  - reordanaco do sistema de equacao                *
c * bvtk      - saida binario para o vtk legacy                       *
c * legacy    - saida legacy do vtk                                   *
c * mpi       - true|fasle                                            *
c * nprcs     - numero de processos do MPI                            *
c * nin       - arquivo de entrada                                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * OBS:                                                              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine read_config(maxmem  
     1                      ,omp_elmt,omp_solver
     2                      ,nth_elmt,nth_solver
     3                      ,reord   ,newton_raphson 
     4                      ,bvtk    ,legacy    
     5                      ,mpi     ,nprcs
     6                      ,my_id   ,nin)
      implicit none
      include 'string.fi'
      character(len=15) string,macro(10)
      character(len=20) ex(10)
      character(len=80) fname
      integer*8 maxmem
      integer nth_elmt,nth_solver
      logical omp_elmt,omp_solver,reord,bvtk,legacy,newton_raphson 
      integer i,j,cont,nmacro
      integer nin,nprcs,my_id
      logical mpi
      integer nincl /7/
      data nmacro /10/
      data macro/'memory         ','omp_elmt       ','omp_solver     ',
     1           'nth_elmt       ','nth_solver     ','reord          ',
     2           'binary_vtk     ','legacy_vtk     ','newton_raphson ',
     3           'help'/
c .....................................................................
c
c ... exemplo
      data ex /'memory          1000','omp_elmt        true',   
     1         'nth_elmt           4','omp_solver      true',
     2         'nth_elmt           2','binary_vtk      true',
     3         'legacy_vtk      true','newton_raphson false',
     4         'reord           true','end                 '/
c .....................................................................
c
c ... arquivo de config
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(nincl, file= fname,status= 'old',err=200,action='read')
c .....................................................................
c
c ...
      cont = 0
      call readmacro(nincl,.true.)
      write(string,'(15a)') (word(j),j=1,15)
      do while ( (string .ne. 'end') .and. (cont .le. 9) )
         cont = cont + 1
c ... memory
         if (string .eq. macro(1)) then
            i = 1
            call readmacro(nincl,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) maxmem
            if(mpi) maxmem = maxmem/nprcs
c ... convertendo de Mbytes para para numero de inteiros e 4 bytes
            maxmem = (maxmem*1024*1024)/4
c .....................................................................
c
c ... 
         elseif (string .eq. macro(2)) then
            i = 2
            call readmacro(nincl,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) omp_elmt
c .....................................................................
c
c ... 
         elseif (string .eq. macro(3)) then
            i = 3
            call readmacro(nincl,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) omp_solver
c .....................................................................
c
c ... 
         elseif (string .eq. macro(4)) then
            i = 4
            call readmacro(nincl,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) nth_elmt
c .....................................................................
c
c ... 
          elseif (string .eq. macro(5)) then
            i = 5
            call readmacro(nincl,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) nth_solver
c .....................................................................
c
c ... 
         elseif (string .eq. macro(6)) then
            i = 6
            call readmacro(nincl,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) reord
            if(my_id.eq.0) then
              write(*,'(2x,a15,1x)',advance='no')'reord         :'
              write(*,*)reord    
            endif
c .....................................................................
c
c ... 
          elseif (string .eq. macro(7)) then
            i = 7
            call readmacro(nincl,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) bvtk
            if(my_id.eq.0) then
              write(*,'(2x,a15,1x)',advance='no')'binary_vtk    :'
              write(*,*)bvtk    
            endif
c .....................................................................
c
c ... 
          elseif (string .eq. macro(8)) then
            i = 8
            call readmacro(nincl,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) legacy
            if(my_id.eq.0) then
              write(*,'(2x,a15,1x)',advance='no')'legacy_vtk    :'
              write(*,*)legacy    
            endif
c .....................................................................
c
c ... 
          elseif (string .eq. macro(9)) then
            i = 9
            call readmacro(nincl,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) newton_raphson
            if(my_id.eq.0) then
              write(*,'(2x,a15,1x)',advance='no')'newton_raphson:'
              write(*,*)newton_raphson    
            endif
c .....................................................................
c
c ... 
          elseif (string .eq. macro(10)) then
            goto 300
          endif 
c .....................................................................
c
c ... 
         call readmacro(nincl,.true.)
         write(string,'(15a)') (word(j),j=1,15)
      end do
c ......................................................................
c
c ...
      if( cont .eq. 0 ) goto 300
      close(nincl)
      return  
c ......................................................................
 100  continue
      print*,'*** Erro in reading the config file !',macro(i)
      call stop_mef()                          
 200  continue
      print*,'File ',trim(fname),' not found !'
      call stop_mef() 
 300  continue      
      if(my_id.eq.0) then
        print*,'*** Error reading solver macro !'
        write(*,'(2x,a)') '******************************'
        write(*,'(2x,a)') 'Example usage of macro config:'
        write(*,'(2x,a)') '------------------------------'
        do i = 1, 10
          write(*,'(2x,a)') ex(i)
        enddo  
        write(*,'(2x,a)') '------------------------------'
        write(*,'(2x,a)') '******************************'
      endif
      call stop_mef()
c ......................................................................
      end
c *********************************************************************
c
c *********************************************************************
c * Data de criacao    : 13/10/2016                                   *
c * Data de modificaco : 27/05/2017                                   *
c * ------------------------------------------------------------------*
c * read_solver_config : leitura das configuracoes basicas do solver  *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * solver    - nao definido                                          *
c * tol       - nao definido                                          *
c * maxit     - nao definido                                          *
c * precond   - nao definido                                          *
c * nkrylov   - nao definido                                          *
c * fhist_log - nao definido                                          *
c * fprint    - nao definido                                          *
c * prename   - prefixo do arquivo de entrada                         *
c * nout      - arquivo de log do solver(it)                          *
c * my_id     - id mpi                                                *
c * nprcs     - numero de processos de mpi                            *
c * nin       - arquivo de entrada                                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * solver    - solver ( 1 - CG)                                      *
c * tol       - tolerancia do solver                                  *
c * maxit     - numero de iteracoes                                   *
c * precond   - precondicionador ( 1 - NONE; 2 - Diag )               *
c * nkrylov   - numero de base de krylov                              *
c * fhist_log - escreve o historico de iteracao do solver             *
c * fprint    - saida de tela do solver                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * OBS:                                                              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine read_solver_config_pm(solver    ,tol
     1                                ,maxit     ,precond
     2                                ,nkrylov   ,fhist_log
     3                                ,fprint
     4                                ,prename   ,nout
     5                                ,alfap     ,alfau
     6                                ,ctol      ,cmaxit 
     7                                ,nprcs     ,my_id
     8                                ,nin)
      implicit none
      include 'string.fi'
      character*15 string,macro(12),ex(12)
      character*80 fname,prename,name
      real*8 tol,ctol,smachn,alfap,alfau
      integer my_id,nprcs,solver,maxit,cmaxit,precond,nkrylov
      integer i,j,nmacro,erro
      integer nin,nout
      logical fhist_log,fprint
      integer nincl /12/
      data nmacro /12/
      data macro/'name           ','it             ','tol            ',
     .           'precond        ','nkrylov        ','histlog        ',
     .           'alfap          ','alfau          ','ctol           ',
     .           'cmaxit         ','fprint         ','help           '/
c ... exemplo
      data ex   /'name         cg','it         1000','tol      1.e-11',
     .           'precond    diag','nKrylov     300','end            ',
     .           '               ','               ','               ',
     .           '               ','               ','               '/
c .....................................................................
c      
c ... arquivo de config
      call readmacrov2(nin,.false.,erro)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(nincl, file= fname,status= 'old',err=200,action='read')
c .....................................................................
c
c ...
      call readmacrov2(nincl,.true.,erro)
      if(erro .eq. 1) goto 100 
      write(string,'(15a)') (word(j),j=1,12)
      do while (string .ne. 'end')
c ... solver
         if (string .eq. macro(1)) then
            i = 1
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(6a)') (word(j),j=1,6)
            call set_solver(solver,nin,my_id)
c .....................................................................
c
c ... maxit
         elseif (string .eq. macro(2)) then
            i = 2
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) maxit
            if(my_id.eq.0) write(*,'(a10,1x,i9)')'It     :',maxit 
c .....................................................................
c
c ... tol
         elseif (string .eq. macro(3)) then
            i = 3
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) tol
            if(tol .eq. 0.d0) tol = smachn()
            if(my_id.eq.0) write(*,'(a10,1x,d13.6)')'Tol    :',tol   
c .....................................................................
c
c ... precond
         elseif (string .eq. macro(4)) then
            i = 4
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(6a)') (word(j),j=1,6)
            call set_precond(word,precond,nin,my_id)
c .....................................................................
c
c ... nkrylov
         elseif (string .eq. macro(5)) then
            i = 5
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) nkrylov
            if(my_id.eq.0) write(*,'(a10,1x,i9)')'nKrylov:',nkrylov
c .....................................................................
c
c ... histlog
         elseif (string .eq. macro(6)) then
            i = 6
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) fhist_log
            if(my_id.eq.0) then
              write(*,'(a10,1x)',advance='no')'Histlog:'
              write(*,*)fhist_log
            endif
c ... alphap 
         elseif (string .eq. macro(7)) then
            i = 7
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) alfap
            if(my_id.eq.0) write(*,'(a10,1x,d13.6)')'alpap  :',alfap
c .....................................................................
c
c ... alphau 
         elseif (string .eq. macro(8)) then
            i = 8
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) alfau
            if(my_id.eq.0) write(*,'(a10,1x,d13.6)')'alpau  :',alfau
c .....................................................................
c
c ... ctol   
         elseif (string .eq. macro(9)) then
            i = 9
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) ctol
            if(my_id.eq.0) write(*,'(a10,1x,d13.6)')'ctol   :',ctol
c .....................................................................
c
c ... cmaxit 
         elseif (string .eq. macro(10)) then
            i = 10
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) cmaxit
            if(my_id.eq.0) write(*,'(a10,1x,i9)')'cmaxit :',cmaxit
c .....................................................................
c
c ... fprint 
         elseif (string .eq. macro(11)) then
            i = 11
            call readmacrov2(nincl,.false.,erro)
            if(erro .eq. 1) goto 300 
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) fprint    
            if(my_id.eq.0) then
              write(*,'(a10,1x)',advance='no')'fprint :'
              write(*,*)fprint   
            endif
c .....................................................................
c
c ... help  
         elseif (string .eq. macro(12)) then
           goto 300
c .....................................................................
         endif
c .....................................................................
c
c ... 
         call readmacrov2(nincl,.true.,erro)
         if(erro .eq. 1) goto 300 
         write(string,'(15a)') (word(j),j=1,15)
      end do
c ......................................................................
c
c ...
      close(nincl)
c ......................................................................
c
c ...
      if(fhist_log) then 
        fname = name(prename,nprcs,17)
        open(nout,file=fname)
        write(nout,'(a)') '#solver it r/r0'
      endif
c ......................................................................
      return  
c ......................................................................
 100  continue
      if(my_id.eq.0) then
        print*,'*** Error reading macro: '
     .      ,trim(macro(i)),' in macro solver!'
      endif
      call stop_mef()                        
 200  continue
      if(my_id.eq.0) then
        print*,'File ',trim(fname),' not found !'
      endif
      call stop_mef()
 300  continue
      if(my_id.eq.0) then
        print*,'*** Error reading solver macro !'
        write(*,'(2x,a)') '******************************'
        write(*,'(2x,a)') 'Example usage of macro solver:'
        write(*,'(2x,a)') '------------------------------'
        do i = 1, 5
          write(*,'(2x,a)') ex(i)
        enddo  
        write(*,'(2x,a)') '------------------------------'
        write(*,'(2x,a)') '******************************'
      endif
      call stop_mef()
c ......................................................................
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 27/09/2016                                    *
c * Data de modificaco : 12/04/2018                                    *
c * ------------------------------------------------------------------ *
c * SET_PRINT_VTK : leitualeitura das configuracoes basicas de excucao *
c * ------------------------------------------------------------------ *
c * Parametros de entrada :                                            *
c * -----------------------------------------------------------------  *
c * fprint    - nao definido                                           *
c * my_id     -id do processo do mpi                                   *
c * nin       - arquivo de entrada                                     *
c * -----------------------------------------------------------------  *
c * Parametros de saida :                                              *
c * -----------------------------------------------------------------  *
c * fprint -                                                           *
c *          quadratic       (1)                                       *
c *          deslocamento    (2)                                       *
c *          pressao         (3)                                       *
c *          delta pressa    (4)                                       *
c *          stress Total    (5)                                       *
c *          stress Biot     (6)                                       *
c *          stress Terzaghi (7)                                       *
c *          fluxo de darcy  (8)                                       *
c *          delta prosidade (9)                                       *
c *          pconsolidation  (10)                                      *
c *          eplastic        (11)                                      *
c *          temperatura     (12)                                      *
c *          fluxo de calor  (13)                                      *      
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine set_print_vtk_pm(fprint,my_id,nin)
      implicit none
      include 'string.fi'
      character*15 string,macro(15)
      character*80 fname
      logical fprint(*),fexit
      integer j,nmacro,my_id
      integer nin
      logical fplastic
      integer nincl /12/
      data nmacro /15/
      data macro/'quadratic     ','end           ','               '
     1          ,'desloc        ','pressure      ','dpressure      '
     2          ,'totalstress   ','biotstress    ','terzaghistress '
     3          ,'darcyflux     ','porosity      ','pconsolidation '
     4          ,'elplastic     ','temperature   ','heatFlux       '/
c ......................................................................
c
c ...
      !  1 - quadratic
      !  2 - end
      !  3 -
      !  4 - desloc
      !  5 - pressure
      !  6 - dpressure
      !  7 - totalstress
      !  8 - biotstress
      !  9 - terzaghistress
      ! 10 - darcyflux
      ! 11 - porosidade
      ! 12 - pconsolidation
      ! 13 - elplastic
      ! 14 - temperature
      ! 15 - heatflux
c ......................................................................
c
c ...
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(nincl, file= fname,status= 'old',err=200,action='read')
c ......................................................................
c
c ...
      do j = 1,nmacro-1
        fprint(j) = .false.
      enddo
c ......................................................................
c
c ...
      call readmacro(nincl,.true.)
      write(string,'(15a)') (word(j),j=1,15)
      do while (string .ne. 'end')
c ... quadratic
        if (string .eq. macro(1)) then
          fprint(1) = .true.
c .....................................................................
c
c ... desloc
        else if (string .eq. macro(4)) then
          fprint(2) = .true.
c .....................................................................
c
c ... pressure 
        elseif (string .eq. macro(5)) then
          fprint(3) = .true. 
c .....................................................................
c
c ... delta pressure 
        elseif (string .eq. macro(6)) then
          fprint(4) = .true.
c .....................................................................
c
c ... stress total  
        elseif (string .eq. macro(7)) then 
          fprint(5) = .true.
c .....................................................................
c
c ... stress biot   
        elseif (string .eq. macro(8)) then
          fprint(6) = .true. 
c .....................................................................
c
c ... stress terzaghi 
        elseif (string .eq. macro(9)) then
          fprint(7) = .true. 
c .....................................................................
c
c ... fdarcy  
        elseif (string .eq. macro(10)) then
          fprint(8) = .true.
c .....................................................................
c
c ... porosidade
        elseif (string .eq. macro(11)) then
          fprint(9) = .true.
c .....................................................................
c
c ... pconsolidation 
        elseif (string .eq. macro(12)) then
          fprint(10) = .true.
c .....................................................................
c
c ... eplastic          
        elseif (string .eq. macro(13)) then
          fprint(11) = .true.
c .....................................................................
c
c ... temperature       
        elseif (string .eq. macro(14)) then
          fprint(12) = .true.
c .....................................................................
c
c ... heatflux          
        elseif (string .eq. macro(15)) then
          fprint(13) = .true.
        endif
c .....................................................................
        call readmacro(nincl,.true.)
        write(string,'(15a)') (word(j),j=1,15)
      end do
c .....................................................................
c
c ...
      do j = 1, nmacro-1
        if(my_id.eq.0)print*,macro(j),':',fprint(j)
      enddo
c ......................................................................
c
c ...
      close(nincl)
      return
c ......................................................................
c
c ...
 200  continue
      if(my_id.eq.0) then
        print*,'File ',trim(fname),' not found !'
      endif
      call stop_mef()
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 23/10/2016                                    *
c * Data de modificaco : 04/03/2017                                    *
c * ------------------------------------------------------------------ *
c * READ_GRAVITY : leitura do campo de gravidade                       *
c * ------------------------------------------------------------------ *
c * Parametros de entrada :                                            *
c * -----------------------------------------------------------------  *
c * g         - nao definido                                           *
c * mg        - nao definido                                           *
c * nin       - arquivo de entrada                                     *
c * -----------------------------------------------------------------  *
c * Parametros de saida :                                              *
c * -----------------------------------------------------------------  *
c * g      - gravidade                                                 *
c * mg     - modulo da gravidade                                       *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine read_gravity(g,mg,nin)
      include 'string.fi'
      character*30 string
      character*80 fname
      integer nin
      integer nincl /12/
      real*8 g(3),mg
c ...
      g(1:3) = 0.d0
c .....................................................................
c
c ...
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(nincl, file= fname,status= 'old',err=900,action='read')
c ......................................................................
c
c ... gx
      call readmacro(nincl,.true.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =910,end =910) g(1)
c .....................................................................
c
c ... gy
      call readmacro(nincl,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =920,end =920) g(2)  
c .....................................................................
c  
c ... gz
      call readmacro(nincl,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =930,end =930) g(3)    
c .....................................................................
c  
c ... mg
      call readmacro(nincl,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =940,end =940) mg    
c .....................................................................
c
      return
c .....................................................................
c
c ...
 900  continue
      if(my_id.eq.0) then
        print*,'File ',trim(fname),' not found !'
      endif
      call stop_mef()
c .....................................................................
c
c ...
  910 continue
      print*,'Erro na leitura da macro (GRAVITY) gx !'
      call stop_mef()
c .....................................................................
c
c ...
  920 continue
      print*,'Erro na leitura da macro (GRAVITY) gy !'
      call stop_mef()
c .....................................................................
c
c ...
  930 continue
      print*,'Erro na leitura da macro (GRAVITY) gz !'
      call stop_mef()
c .....................................................................
c
c ...
  940 continue
      print*,'Erro na leitura da macro (GRAVITY) mg !'
      call stop_mef()
c .....................................................................
c
c ...
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 09/10/2016                                    *
c * Data de modificaco : 08/11/2016                                    *
c * ------------------------------------------------------------------ *
c * COUNT_NODE_VERTICE : conta o numero de nos de vertices             *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ix    - conectividade                                              *
c * fnno  - nao definido                                               *
c * nnodev- nao definido                                               *
c * nnode - numero total de nos                                        *
c * numel - numero de elementos                                        *
c * numat - numero de materiais                                        *
c * maxnov- numero max. de nos geometicos por elemento                 *
c * maxno - numero max. de nos por elemento                            *
c * fcount- true : conta os numeros de vertices                        * 
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * fnno  - (1 - no de vertice : 0 - no intermediario)                 *
c * nnodev- numero de nos dos vertices(fcount = .true.)                *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine count_node_vertice(ix,fnno,nnodev,nnode,numel
     .                             ,maxnov,maxno,fcount)
      implicit none
      integer ix(maxno+1,*),nnodev,nnode,numel,maxnov,maxno
      integer fnno(*),i,j,no
      logical fcount
c ...
      fnno(1:nnode) = 0
c .....................................................................
c
c ...
      do i = 1, numel
        do j = 1, maxnov
          no = ix(j,i)
          if( no .ne. 0) then
            if(fnno(no) .eq. 0) fnno(no) = 1
          endif
        enddo
      enddo
c .....................................................................
c
c ... conta o numero de vertices
      if(fcount) then
        nnodev = 0
        do i = 1, nnode
          if(fnno(i) .eq. 1) then
            nnodev = nnodev + 1
          endif
        enddo
      endif  
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 20/01/2017                                    *
c * Data de modificaco : 02/04/2017                                    *
c * ------------------------------------------------------------------ *
c * CHECK_ELEMENT : verifica a escolha dos elementos                   *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ie      - tipo de elemento por material                            *
c * numat   - numero de materiais                                      *
c * plastic - true|false                                               *
c * vpro    - true|false                                               *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine check_element(ie,numat,plastic,vprop)
      implicit none
      common /el_lib3/ element_library_3
      character(len=1024) str
      character(len=60) element_library_3(8)
      integer ie(*),numat,i,pel(4),eel(4),vel(4),nvel(4),ty
      logical flag,plastic,vprop(*)
c ... elementos de plasticidade
      pel(1) = 35
      pel(2) = 36
      pel(3) = 37
      pel(4) = 38
c ...
      eel(1) = 15
      eel(2) = 16
      eel(2) = 17
      eel(2) = 18
c ... elementos com propriedades nos pontos de integracao
      vel(1) = 16
      vel(2) = 18
      vel(3) = 36
      vel(4) = 38  
c ... elementos com propriedades nos pontos de integracao
      nvel(1) = 15
      nvel(2) = 17
      nvel(3) = 35
      nvel(4) = 37  
c .....................................................................
c
c ... checa a variacao das propriedade nos elementos
      flag = .false.
      do i = 1, numat
        ty = ie(i)
        if( ty .eq. vel(1) .or. ty .eq. vel(2) 
     . .or. ty .eq. vel(3) .or. ty .eq. vel(4) ) then
          flag = .true.
        endif
      enddo
c .....................................................................
c
c ...
      if(vprop(1)) then
        if(.not. flag) then 
          print*,'Erro: Invalid elements!!'
          str ='Elements with properties at integration points are: '
          print*,trim(str),vel(1:4)
          goto 100
        endif
      else
        if(flag) then 
          print*,'Erro: Invalid elements!!'
 
          str ='Elements without properties at integration points are: '
          print*,trim(str),nvel(1:4)
          goto 100
        endif    
      endif  
c .....................................................................
c
c ... checa elementos plasticos
      flag = .false.
      do i = 1, numat
        ty = ie(i)
        if( ty .eq. pel(1) .or. ty .eq. pel(2)
     .   .or. ty .eq. pel(3) .or. ty .eq. pel(4) ) then
          flag = .true.
        endif
      enddo
c .....................................................................
c
c ...
      if(plastic) then
        if(.not. flag) then 
          print*,'Erro: Invalid elements!!'
          print*,'Plastic elements are: ',pel(1:2)
          goto 100
        endif
      else
        if(flag) then 
          print*,'Erro: Invalid elements!!'
          print*,'Elastic elements are: ',eel(1:2)
          goto 100
        endif    
      endif  
c .....................................................................
c
c ...
      return
c ......................................................................
c
c ...
  100 continue
      print*,'Elements library:' 
      do i = 1, 8
        print*,element_library_3(i) 
      enddo
      call stop_mef()
c ......................................................................
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 20/01/2017                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *
c * INIT_PC0 : inicializa do pressao de consolidacao inicial           *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ie      - tipo de elemento por material                            *
c * ix      - conectividade do material                                *
c * e       - propriedade do material                                  *
c * plastic - true|false                                               * 
c * npi     - numero de pontos de integracao                           *
c * numel   - numero de elementos                                      *
c * nen     - numero de nos por elemento                               *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine init_pc0(ie,ix,e,plastic,npi,numel,nen)
      implicit none
      include 'termprop.fi'
      integer ie(*),ix(nen+1,*),numel,npi,ma,nen,i,j
      real*8 e(prop,*),plastic(3,npi,*)
      do i = 1, numel 
        ma             = ix(nen+1,i)        
        do j = 1, npi 
          plastic(3,j,i) = e(12,ma)
        enddo
      enddo
      return
      end
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *
c * INIT_STRESS:                                                       *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine init_stress(ix,ie,e,x,tx0,numel,nenv,nen,ndm,ntn)
      implicit none
      include 'gravity.fi'
      include 'termprop.fi'
      integer ie(*),ix(nen+1,*),numel,npi,ma,nen,nenv,ndm,ntn,i,j,no
      real*8 e(prop,*),tx0(ntn,*),x(ndm,*),ps,ro,x3
c .....................................................................
c
c ...
      do i = 1, numel 
        ma  = ix(nen+1,i)
        ps  = e(2,ma)
        ro  = e(6,ma)*1.d-6
        do j = 1, nenv
          no     = ix(j,i) 
          x3     = x(3,no)
c ... tensao inicial
          tx0(3,no)   = -ro*gravity_mod*(10.0d0-x3)
          tx0(1,no)   = tx0(3,no)
          tx0(2,no)   = tx0(3,no)
c         tx0(1,no)   = (ps/(1.0-ps))*tx0(3,no)
c         tx0(2,no)   = tx0(1,no) 
          tx0(4:6,no) = 0.d0
c ......................................................................
        enddo
      enddo
c     open(15, file= 'stress.txt',action= 'write')
c     do i = 1, 44
c       write(15,'(i9,6es14.6)')i,(tx0(j,i),j=1,6)
c     enddo  
c     stop
      return
      end
c ......................................................................
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 29/01/2017                                    * 
c * ------------------------------------------------------------------ *  
c * READPNODE : le o arquivo auxiliar com os nos que terao alguma      *
c * de suas grandeza impressas no tempo                                *
c * ------------------------------------------------------------------ *
c * parametros de entrada :                                            *
c * ------------------------------------------------------------------ *
c * fname   - nome do arquivo de entrada                               *
c * i_no    -                                                          *
c * i_nfile -                                                          *
c * num_node-                                                          *
c * flag    -                                                          *
c * nout    - arquivo de entrada                                       *
c * ------------------------------------------------------------------ *
c * parametros de saida                                                *
c * ------------------------------------------------------------------ *
c * i_no    -ponteiros para os nos                                     *
c * i_nfile -ponteiros para os arquivos                                *
c * num_node-numero de nos                                             *
c * flag    -verifica sucesso na abertura do arquivo                   *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine readpnode(fname,i_no,i_nfile,num_node,flag,nout)
      use Malloc
      implicit none
      include 'string.fi'
      integer*8 i_no,i_nfile
      integer num_node,i,j,no
      integer nout
      logical flag
      character*80 fname
      character*30 string
c .....................................................................
c
c ... 
      open(nout,file=fname,status= 'old' , err=1000 , action='read')
      flag = .true.
c ... numero total de nos      
      call readmacro(nout,.true.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*) num_node
      print*,num_node,' print set nodes'
c .....................................................................
c
c ...
      if(num_node .gt. 50) then
        print*,'number max print node is 50 !!!' 
        flag = .false.
        return 
      endif
c .....................................................................
c
c ... memoria para os nos
      i_no        = alloc_4('no      ',1,num_node)
c ... memoria para os arquivos      
      i_nfile     = alloc_4('nnodew  ',1,num_node)
c ... lendo nos
      do i = 1, num_node
         call readmacro(nout,.true.)
         write(string,'(30a)') (word(j),j=1,30)
         read(string,*) no
         ia(i_no+i-1) = no
         ia(i_nfile+i-1) = 51 + (i-1)
      enddo
c .....................................................................
      close(nout)
      return
 1000 continue
      print*,'File ',trim(fname),' not found !'
      flag = .false.
      end
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 29/01/2017                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * READPPG : le o arquivo auxiliar com os nos que terao alguma        *
c * de suas grandeza impressas no tempo                                *
c * -----------------------------------------------------------------  *
c * parametros de entrada :                                            *
c * -----------------------------------------------------------------  *
c * fname   - nome do arquivo de entrada                               *
c * i_el    -                                                          *
c * i_nfile -                                                          *
c * num_nel -                                                          *
c * flag    -                                                          *
c * nout    - arquivo de entrada                                       *
c * -----------------------------------------------------------------  *
c * parametros de saida                                                *
c * -----------------------------------------------------------------  *
c * i_no    -ponteiros para os nos                                     *
c * i_nfile -ponteiros para os arquivos                                *
c * num_node-numero de nos                                             *
c * flag    -verifica sucesso na abertura do arquivo                   *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine readppi(fname,i_el,i_nfile_pi,num_pel,flag,nout)
      use Malloc
      implicit none
      include 'string.fi'
      integer*8 i_el,i_nfile_pi
      integer num_pel,i,j,el
      integer nout
      logical flag
      character*80 fname
      character*30 string
c .....................................................................
c
c ... 
      open(nout,file=fname,status= 'old' , err=1000 , action='read')
      flag = .true.
c ... numero total de nos      
      call readmacro(nout,.true.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*) num_pel
      print*,num_pel,' print set elements'
c .....................................................................
c
c ...
      if(num_pel .gt. 50) then
        print*,'number max print elements is 50 !!!' 
        flag = .false.
        return 
      endif
c .....................................................................
c
c ... memoria para os nos
      i_el        = alloc_4('elw     ',1,num_pel)
c ... memoria para os arquivos      
      i_nfile_pi  = alloc_4('pgw     ',1,num_pel)
c ... lendo nos
      do i = 1, num_pel
         call readmacro(nout,.true.)
         write(string,'(30a)') (word(j),j=1,30)
         read(string,*) el
         ia(i_el+i-1)  = el
         ia(i_nfile_pi+i-1) = 101 + (i-1)
      enddo
c .....................................................................
      close(nout)
      return
 1000 continue
      print*,'File ',trim(fname),' not found !'
      flag = .false.
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 26/04/2017                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * BOUND_NODALLOADS :     auxiliar com os nos que terao alguma        *
c * de suas grandeza impressas no tempo                                *
c * -----------------------------------------------------------------  *
c * parametros de entrada :                                            *
c * -----------------------------------------------------------------  *
c * x       - nome do arquivo de entrada                               *
c * id      - valores prescrritos                                      *
c * u0      - valores inicias                                          *
c * nnode   - numero de nois totais                                    *
c * ndf     - graus de liberdade                                       *
c * ndm     - dimensao                                                 *
c * -----------------------------------------------------------------  *
c * parametros de saida                                                *
c * -----------------------------------------------------------------  *
c * u0      -valores inicias atualizados com a macro nodalloads        *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine bound_nodalloads(x,id,u0,nload,nnode,ndf,ndm)
      implicit none
      real*8 x(ndm,*),u0(ndf,*),vc(3),c,t
      integer id(ndf,*),nnode,ndf,ndm,i,j,k,nc,nload(ndf,*)
       
c... Cargas nodais e deslocamentos prescritos no tempo t:
      t = 0.d0 
      do i = 1, nnode
        do j = 1, ndf
          nc = nload(j,i)
          if( nc .gt. 0 ) then
            k = id(j,i)
            if (k .eq. 1) then
              call tload(nc,t,x(1,i),u0(j,i),c,vc)
              u0(j,i) = c
            endif
          endif
        enddo
      enddo
c .....................................................................
      return
      end
c **********************************************************************