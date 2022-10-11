c *********************************************************************
c * READMETIS: leitura apenas das conectividades da malha necessaria  *
c * para o particionador do metis                                     *
c * ----------------------------------------------------------------- *
c * Paramtros de entrada:                                             *
c * ----------------------------------------------------------------- *
c * i_ix0  -> ponteiro conectividades sem o MATERIAL                  *
c * i_nen  -> ponteiro para o numero de nos por elemento              *
c * nnodev -> numero de nos dos vertices                              *
c * nnode  -> nao definido                                            *
c * numel  -> numeor de elementos                                     *
c * maxnov -> numero maximo de nos de vertices por elem nto           *
c * maxno  -> numero maximo de nos por elementos                      *
c * nin    -> arquivo de entrada                                      *
c * ----------------------------------------------------------------- *
c * Paramtros de  Saida:                                              *
c * ----------------------------------------------------------------- *
c * ix0(maxno,numel)-> conectividade                                  *
c * nen(numel)      -> numnero de nos por elemento                    *
c * nnodev          -> numero de nos dos vertices                     *
c * nnode           -> numero de nos dos totais                       *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine read_metis(i_ix0,i_nen ,nnodev,nnode,numel
     .                     ,maxnov,maxno,nin  ,verbose)
c ===
      use Malloc
      implicit none
      include 'string.fi'
      include 'elementos.fi'
      include 'parallel.fi'
c ... ponteiros
      integer*8 i_ix0,i_nen,i_nincid,i_incid
c ... malha
      integer numel,maxno,maxnov,nnode,nnodev,maxgrade
c ... auxiliar
      integer dum
      integer i,j,totnel
c ... arquivo
      integer nincl,naux,nin
      logical bin,fileaux,verbose
      data bin /.false./, fileaux/.false./
c ... macros
      integer nmc
      data nmc/12/
      character*80 fname
      character*30 string
      character*15 macro(12),rc
      data macro/'tria3         ','quad4         ','tetra4        ',
     .           'hexa8         ','insert        ','return        ',
     .           'mesh          ','prism6        ','quad4bin      ',
     .           'tetra4bin     ','hexa8bin      ','              '/
c ...
      data nincl /1000/
c ...
      totnel = 0
c ... lendo macros
100   continue
      call readmacro(nin,.true.)
      write(rc,'(15a)'),(word(j),j=1,15)
      do j = 1,nmc
        if(rc .eq. macro(j)) goto 200
      enddo
      goto 100
c .....................................................................
c
c ...
200   continue
c ... tria3
      if( j .eq. 1)then
        ntria3(1) = 0
        if(verbose) print*, 'loading tria3 ...'
        call readconmetis(ia(i_ix0),ia(i_nen),numel,3,maxno
     .                   ,ntria3(1),nin)
        if(verbose) print*, 'load.'
        ntria3(2) = totnel + 1
        totnel    = totnel + ntria3(1)
        if( totnel .eq. numel ) goto 300
        go to 100
c .......................................................................
c
c ... quad4
      else if(j.eq.2)then
        nquad4(1) = 0
        if(verbose) print*, 'loading quad4 ...'
        call readconmetis(ia(i_ix0),ia(i_nen),numel,4,maxno
     .                   ,nquad4(1),nin)
        if(verbose) print*, 'load.'
        nquad4(2) = totnel + 1
        totnel    = totnel + nquad4(1)
        if( totnel .eq. numel ) goto 300
        go to 100
c .......................................................................
c
c ... tetra4
      else if(j.eq.3)then
        ntetra4(1) = 0
        maxnov = max(4,maxnov)
        if(verbose) print*, 'loading tetra4 ...'
        call readconmetis(ia(i_ix0),ia(i_nen),numel,4,maxno
     .                   ,ntetra4(1),nin)
        if(verbose) print*, 'load.'
        ntetra4(2) = totnel + 1
        totnel     = totnel + ntetra4(1)
        if( totnel .eq. numel ) then
c ... transforma os elementos lineares em quadraticos (10 nos)
          if(maxno .eq. 10) then
c ...
            i_nincid = alloc_4('nincid  ',1,nnodev)
            call mynodegrade(ia(i_ix0)   ,ia(i_ix0)
     .                      ,nnodev      ,numel,maxno,maxno
     .                      ,ia(i_nincid),maxgrade,.false.)
            i_incid  = alloc_4('incid   ',maxgrade,nnodev)
            call myelmincid(ia(i_ix0),ia(i_ix0)
     .                   ,ia(i_incid),ia(i_nincid),nnodev
     .                   ,numel      ,maxno       ,maxno
     .                   ,maxgrade   ,.false.)
c ... gera a conectividade dos elementos quadraticos
            call mk_elconn_quad_v1(ia(i_ix0),ia(i_incid),ia(i_nincid)
     .                            ,numel   ,nnodev     ,nnode
     .                            ,maxnov  ,maxno      ,maxno
     .                            ,maxgrade)
c .......................................................................
c
c ...
            call ml_nen_quad(ntetra4(2),ntetra4(1),ia(i_nen),10)
c .......................................................................
c
c ...
            i_incid     = dealloc('incid   ')
            i_nincid    = dealloc('nincid  ')
          endif
c .......................................................................
          if( totnel .eq. numel ) goto 300
        endif
c .......................................................................
c
c ...
        go to 100
c .......................................................................
c
c ... hexa8
      else if(j.eq.4)then
        nhexa8(1) = 0
        maxnov = max(8,maxnov)
        if(verbose) print*, 'loading hexa8 ...'
        call readconmetis(ia(i_ix0),ia(i_nen),numel,8,maxno
     .                   ,nhexa8(1),nin)
        if(verbose) print*, 'load.'
        nhexa8(2) = totnel + 1
        totnel    = totnel + nhexa8(1)
        if( totnel .eq. numel ) then
c ... transforma os elementos lineares em quadraticos (20 nos)
          if(maxno .eq. 20) then
c ...
            i_nincid = alloc_4('nincid  ',1,nnodev)
            call mynodegrade(ia(i_ix0)   ,ia(i_ix0)
     .                      ,nnodev      ,numel,maxno,maxno
     .                      ,ia(i_nincid),maxgrade,.false.)
            i_incid  = alloc_4('incid   ',maxgrade,nnodev)
            call myelmincid(ia(i_ix0),ia(i_ix0)
     .                   ,ia(i_incid),ia(i_nincid),nnodev
     .                   ,numel      ,maxno       ,maxno
     .                   ,maxgrade   ,.false.)
c ... gera a conectividade dos elementos quadraticos
            call mk_elconn_quad_v1(ia(i_ix0),ia(i_incid),ia(i_nincid)
     .                            ,numel   ,nnodev     ,nnode
     .                            ,maxnov  ,maxno      ,maxno
     .                            ,maxgrade)
c .....................................................................
c
c ...
            call ml_nen_quad(nhexa8(2),nhexa8(1),ia(i_nen),20)
c .....................................................................
c
c ...
            i_incid     = dealloc('incid   ')
            i_nincid    = dealloc('nincid  ')
          endif
c .....................................................................
c
c ...
         go to 300
c .......................................................................
        endif
c .......................................................................
c
c ...
        go to 100
c .......................................................................
c
c ... insert(desvio do fluxo do arquivo de entrada)
      else if(j.eq.5)then
        fileaux =.true.
        naux = nin
        call readmacro(nin,.false.)
        write(fname,'(80a)')(word(j),j=1,80)
        open(nincl,file=fname,status='old',err=225)
        nin=nincl
        go to 100
  225   continue
        print*,'File ',trim(fname),' File not found !'
        call finalize()
c .......................................................................
c
c ... return (retorn ao arquivo de entrada principal
      else if(j.eq.6)then
        close(nincl)
        nin=naux
        go to 100
c .......................................................................
c
c ... mesh
      else if(j.eq.7)then
        call parameters(nnode,numel,dum,maxno,dum,dum,dum,dum,nin)
c ...
        nnodev = nnode
        maxnov = 0
c ......................................................................
c
c ...
        i_nen  = alloc_4('nen     ',    1,numel)
        i_ix0  = alloc_4('ix0     ',maxno,numel)
        call mzero(ia(i_ix0),numel*maxno)
c ......................................................................
        go to 100
c ......................................................................
c
c ... prism6
      else if(j.eq.8)then
        nprism6(1) = 0
        if(verbose) print*, 'loading prism6 ...'
        call readconmetis(ia(i_ix0),ia(i_nen),numel,6,6
     .                   ,nprism6(1),nin)
        if(verbose) print*, 'load.'
        nprism6(2) = totnel + 1
        totnel     = totnel + nprism6(1)
        if( totnel .eq. numel ) goto 300
        go to 100
      endif
c ... quad4bin
c     else if(j.eq.9)then
c        bin = .true.
c ... numero de elementos quad4
c        nquad4(2) = 1
c        call readmacro(nin,.false.)
c        write(string,'(30a)') (word(i),i=1,30)
c        read(string,*) nquad4(1)
c        print*,word,nquad4(1)
c .....................................................................
c
c ... arquivo auxiliar binario de elementos
c        naux = nin
c        call readmacro(nin,.false.)
c        write(fname,'(80a)')(word(j),j=1,strl)
c        open(unit=nincl,file=fname,access='stream'
c     .      ,form='unformatted',convert='big_endian')
c        nin=nincl
c        go to 300
c ... tetra4bin
c      else if(j.eq.10)then
c        bin = .true.
c ... numero de elementos tetra4
c        nquad4(2) = 1
c        read(word,*) ntetra4(1)
c .....................................................................
c        go to 100
c ... hexa8bin
c      else if(j.eq.11)then
c        bin = .true.
c ... numero de elementos hexa8
c        nhexa8(2) = 1
c        call readmacro(nin,.false.)
c        write(string,'(30a)') (word(i),i=1,30)
c        read(string,*) nhexa8(1)
c        print*,word,nhexa8(1)
c .....................................................................
c
c ... arquivo auxiliar binario de elementos
c        call readmacro(nin,.false.)
c        naux = nin
c        write(fname,'(80a)')(word(j),j=1,strl)
c        open(unit=nincl,file=fname,access='stream'
c     .      ,form='unformatted',convert='big_endian')
c        nin=nincl
c        go to 100
c ......................................................................
c
c ...
300   continue
c ...
      if(fileaux .or. bin) then
        close(nin)
        nin = naux
      endif
      close(nin)
c .....................................................................
c
c ...
      return
      end
c *********************************************************************
c
c *********************************************************************
c * READCONMETIS: leitura das conectividades sem o material           *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine readconmetis(el,nen,numel,lnen,maxno,nel,nin)
      implicit none
      include 'string.fi'
      integer numel,lnen,maxno
      integer nin
      integer el(maxno,numel),nen(*)
      integer j,m,k,nel,dum
      character*12 string
c ...
      do j = 1, numel
        read(nin,*)k,(el(m,k),m=1,lnen),dum
        nen(k) = lnen
      enddo
      nel = numel
      call readmacro(nin,.true.)
      return
c ......................................................................
c     nel = 0
c     call readmacro(nin,.true.)
c     write(string,'(12a)') (word(j),j=1,12)
c     do 100 while(string .ne. 'end')
c        read(string,*,err = 200,end = 200) k
c        if(k .lt. 1 .or. k .gt. numel) goto 200
c        do 10 j = 1, lnen
c           call readmacro(nin,.false.)
c           write(string,'(12a)') (word(m),m=1,12)
c           read(string,*,err = 200,end = 200) el(j,k)
c  10    continue
c        nen(k) = lnen
c        nel    = nel + 1
c        call readmacro(nin,.true.)
c        write(string,'(12a)') (word(j),j=1,12)
c 100 continue
c     return
c ......................................................................
 200  continue
      print*,'*** Error reading elements !'
      call finalize()
      return
      end
c *********************************************************************
c
c *********************************************************************
c * READCONMETISBIN: leitura das conectividades sem o material        *
c * arquivo auxiliar binario                                          *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine readconmetisbin(el,numel,nen,nin)
      implicit none
      integer el(nen,numel)
      integer numel,nen,nin,dum
      integer j,m,k
      do j = 1, numel
         read(nin)k,(el(m,k),m=1,nen),dum
c        write(*,*)k,(el(m,k),m=1,nen)
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine ml_nen_quad(ni,nf,nen,lnen)
      implicit none
      integer nen(*),ni,nf,i,lnen
      do i = ni, nf
        nen(i) = lnen
      enddo
      return
      end
c *********************************************************************
