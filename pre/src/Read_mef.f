c *********************************************************************
c * Data de criacao    : 00/00/0000                                   *
c * Data de modificaco : 20/03/2017                                   *
c * ----------------------------------------------------------------- *
c * READ_MEF : leitura de arquivo para o mefpar atraves da rdat       *
c * modifcada                                                         *
c * ----------------------------------------------------------------- *
c * Paramtros de entrada:                                             *
c * ----------------------------------------------------------------- *
c * nnodev -> numero de nos dos vertices                              *
c * nnode  -> nao definido                                            *
c * numel  -> numeor de elementos                                     *
c * maxnov -> numero maximo de nos de vertices por elem nto           *
c * maxno  -> numero maximo de nos por elementos                      *
c * ndf     - numero max. de graus de liberdade por no                *
c * ndft    - numero max. de graus de liberdade por no                *
c * ndm     - dimensao (1, 2 ou 3)                                    *
c * npi     - numero de pontos de integracao                          *
c * i_ix    - ponteiro para conetividades                             *
c * i_id    - ponteiro para restricoes nodais (mecanico)              *
c * i_ie    - ponteiro para materiais                                 *
c * i_nload - ponteiro para o arranjo nload (mecanico)                *
c * i_eload - ponteiro para o arranjo eload (mecanico)                *
c * i_e     - ponteiro para o arranjo e                               *
c * i_x     - ponteiro para o arranjo x                               *
c * i_f     - ponteiro para o arranjo f (mecanico)                    *
c * i_u     - ponteiro para o arranjo u (mecanico)                    *
c * i_tx0   - ponteiro para o arranjo tx0 (mecanico-poromec)          *
c * i_v     - ponteiro para o arranjo v (mecanico)                    *
c * i_a     - ponteiro para o arranjo a (mecanico)                    *
c * i_idt    - ponteiro para restricoes nodais  (termico)             *
c * i_nloadt - ponteiro para o arranjo nloadt   (termico)             *
c * i_eloadt - ponteiro para o arranjo eloadt   (termico)             *
c * i_ft     - ponteiro para o arranjo ft       (termico)             *
c * i_ut     - ponteiro para o arranjo ut       (termico)             *
c * i_vt     - ponteiro para o arranjo vt       (termico)             *
c * i_w     - ponteiro para o arranjo w         (termico)             *
c * lines  ->                                                         *
c * nlines ->                                                         *
c * verbose-> mode verbose da rdat                                    *
c * rload  ->                                                         *
c * fprop  -> nome dos arquivos de propriedades                       *
c * nin    -> arquivo de entrada                                      *
c * ----------------------------------------------------------------- *
c * Paramtros de  Saida:                                              *
c * lines  -> linhas lidas de macro-comandos apos o end mesh          *
c * nlines -> numero de linha em lines                                *
c * plines -> linhas lidas de macro-comandos pre mesh                 *
c * pnlines-> numero de linha em plines                               *
c * rload  -> macro-comandos descritiva do problema lidas ou nao      *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * Obs:                                                              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine read_mef(nnodev,nnode,numel,numat,maxnov,maxno
     1                   ,ndf   ,ndft ,ndm  ,npi
     2                   ,i_ix  ,i_ie ,i_e  ,i_x
     3                   ,i_id  ,i_nload    ,i_eload ,i_eloadp
     4                   ,i_f   ,i_u        ,i_tx0   ,i_v    ,i_a
     5                   ,i_idt ,i_nloadt   ,i_eloadt
     6                   ,i_ft  ,i_ut       ,i_ut0   ,i_du   ,i_vt
     7                   ,i_w
     8                   ,lines,nlines,plines,pnlines
     9                   ,verbose,rload  ,fprop  ,ncont,nin)
c ===
      implicit none
      include 'string.fi'
c ... malha
      integer nnodev,nnode,numel,numat,maxnov,maxno,ndf,ndft,ndm,dum,npi
      integer plastic
      integer*8 i_ix,i_ie,i_e,i_x
c ... mecanico - poromec
      integer*8 i_f,i_u,i_v,i_a,i_tx0,i_id,i_nload,i_eload,i_eloadp
c ... termico
      integer*8 i_ft,i_ut,i_ut0,i_du,i_vt,i_w,i_idt,i_nloadt,i_eloadt
c ...
      character*80 fprop(*)
c ... arquivo
      real*8 g(3)
      integer nin
      character*200 lines(*)
      character*200 plines(*),string2
      integer nlines,ncont,nmc,pnlines
      character*15 macro(6),string,rc
      logical verbose,rload(*)
c ... auxiliar
      integer j
c
c ...
      data macro/'mesh           ','config         ','solver         ',
     .           'setprint       ','gravity        ','conseq         '/
      data nmc /6/
c .....................................................................
c =====================================================================
c
c ===
c ...
      pnlines = 0
  100 continue
      call readmacro(nin,.true.)
      write(rc,'(15a)'),(word(j),j=1,15)
      do 200 j = 1, nmc
        if (rc .eq. macro(j)) go to 300
  200 continue
      go to 100
  300 continue
c ......................................................................
c
c ......................................................................
      go to ( 400, 450, 450,
     .        450, 450, 450) j
c ......................................................................
c
c ...
  400 continue
      call rdatm(nnodev,nnode,numel,numat,maxnov ,maxno
     1          ,ndf   ,ndft ,ndm  ,npi  ,dum  ,dum
     2          ,i_ix  ,i_ie,i_e,i_x
     3          ,i_id  ,i_nload,i_eload  ,i_eloadp
     4          ,i_f   ,i_u,i_v,i_vt     ,i_a ,i_tx0
     5          ,i_idt ,i_nloadt,i_eloadt
     6          ,i_ft  ,i_ut    ,i_ut0   ,i_du ,i_w
     7          ,nin,verbose,rload       ,fprop,ncont)
      go to 2000
c ....................................................................
c
c ...
  450 continue
      pnlines = pnlines + 1
      call readmacro(nin,.false.)
      write(plines(pnlines),*)  (word(j),j=1,strl)
      plines(pnlines) = trim(rc) // ' ' // trim(plines(pnlines))
      go to 100
c .....................................................................
c
c .....................................................................
c
c .....................................................................
c
c ... le os restantes da macro comando apos "end mesh"
 2000 continue
      call rmacros(lines,nlines,nin)
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c *GETINFOTERM : obter informcoes do terminal                         *
c * ----------------------------------------------------------------- *
c * Parametros de entrada                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * Parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine getinfoterm(maxmem,filein,fileout,nprcs ,ovlp
     1                      ,vtkf  ,vtkp  ,mefp   ,vtkbin,vtkgr
     2                      ,nin)
      implicit none
      character*80 filein,filearg,fileout
      integer*8 maxmem
      integer nprcs,j,nin
      integer nargs,ierr
      logical ovlp
      logical vtkf,vtkp,mefp,vtkbin,vtkgr
c ... argumentos por linha de comando
      nargs = iargc()
      if(nargs .eq. 1) then
        call getarg(1,filearg)
        call args(maxmem,filearg,filein,fileout,nprcs,ovlp
     1           ,vtkf,vtkp,mefp,vtkbin
     2           ,vtkgr)
        open(nin,file=filein,status='old',err=25)
        return
   25   continue
        print*,'File ',trim(filein),' not found !'
        call finalize()
      endif
c .....................................................................
c
c ... entrada de dados pelo terminal
   10 continue
      print*,'Input file:'
      read(*,'(a)')filein
      open(nin,file=filein,status='old',err=15)
      goto 20
   15 continue
      print*,'File ',trim(filein),' not found !'
      goto 10
   20 continue
      print*,'Output file prefix:'
      read(*,*)fileout
      print*,'Number of processes:'
      read(*,*)nprcs
c
   30 continue
      print*,'Domain sub-division method:'
      print*,'  1 = non-overlapping'
      print*,'  2 = overlapping'
      read(*,*)j
      if (j .eq. 1) then
        ovlp = .false.
      elseif (j .eq. 2) then
        ovlp = .true.
      else
        goto 30
      endif
   40 continue
      print*,'Output.'
      print*,'  1 = Global mesh (vtk).'
      print*,'  2 = Partitioned mesh (vtk).'
      print*,'  4 = Partitioned mesh (mefpar).'
      print*,' 10 = Global mesh with loading(vtk).'
      print*,' More of an option at the same time.'
      print*,'  3 = option 1 e 2'
      print*,'  5 = option 1 e 4'
      print*,'  6 = option 2 e 4'
      print*,'  7 = option 1, 2 e 4'
      print*,' 13 = option 1, 2 e 10'
      print*,' 15 = option 1, 4 e 10'
      print*,' 16 = option 2, 4 e 10'
      print*,' 17 = option 1, 2 , 4 e 10 '
      read(*,*)j
c ...
      if (j .eq. 1) then
        vtkf = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 2) then
        vtkp = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 3) then
        vtkf = .true.
        vtkp = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 4) then
        mefp = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 5) then
        vtkf = .true.
        mefp = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 6) then
        vtkp = .true.
        mefp = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 7) then
        vtkf = .true.
        vtkp = .true.
        mefp = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 10) then
        vtkgr = .true.
c .....................................................................
c
c ...

      elseif (j .eq. 13) then
        vtkf  = .true.
        vtkp  = .true.
        vtkgr = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 15) then
        vtkf = .true.
        mefp = .true.
        vtkgr = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 16) then
        vtkp = .true.
        mefp = .true.
        vtkgr = .true.
c .....................................................................
c
c ...
      elseif (j .eq. 17) then
        vtkf = .true.
        vtkp = .true.
        vtkgr = .true.
        mefp = .true.
c .....................................................................
c
c ...
      else
        goto 40
      endif
c .....................................................................
c
c ... formato de escrita vtk
      if(vtkp .or. vtkf .or. vtkgr) then
        print*,' Binary vtk 1'
        print*,' Ascii vtk  2'
        read(*,*)j
        if(j .eq. 1) then
          vtkbin = .true.
        elseif(j.eq.2)then
          vtkbin = .false.
        endif
      endif
c .....................................................................
c
c ...
      return
      end
c *********************************************************************
c
c *********************************************************************
c * RMACROS: Le as linhas apos a end mes                              *
c *********************************************************************
       subroutine rmacros(lines,nlines,nin)
c ===
       implicit none
       include 'string.fi'
       include 'parallel.fi'
       integer nin
       integer i,j,k
       character*15 rc
       integer nlines
       character*200 lines(*)
c ======================================================================
c ===
c ... Leitura de macro-comandos
c
      k = 1
  100 continue
      call readmacro(nin,.true.)
      write(rc,'(15a)') (word(j),j=1,15)
c      print*,rc,line
      if (rc .eq. 'stop          ') go to 200
c ......................................................................
c
c ... Comentario ou macro-comandos apos mesh:
c
      j = 0
      i = maxstrl
      do while(j .eq. 0)
         if (line(i) .ne. ' ') j = i
         i = i-1
      enddo
      write(lines(k),'(200a)') line(1:j)
c     print*,lines(k)
      k = k +1
      go to 100
c ....................................................................
  200 continue
c       print *,rc
        write(lines(k),'(200a)')rc
        nlines = k
c        print *,lines(k)
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * ARGS : leitura da opcoes apartir de um arquivo auxiliar formencido*
c * atraves da linha de comando                                       *
c * ----------------------------------------------------------------- *
c * Parametros de entrada                                             *
c * ----------------------------------------------------------------- *
c * filearg - arquivo de comando auxiliar                             *
c * ----------------------------------------------------------------- *
c * Parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * maxmem  - memoria                                                 *
c * filein  - arquivo de entrada                                      *
c * fileout - arquivo de saida                                        *
c * nprcs   - numero de divisoes da malha                             *
c * olvp    - metodo de particionamento (true - ovelaping)            *
c * vtkf    - malha global com particionamento formato vtk            *
c * vtkp    - particionamento vtk                                     *
c * mefp    - particionamento mefpar                                  *
c * vtkbin  - arquivo vtk binario                                     *
c * vtkgr   - malha com os carregamento globais                       *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine args(maxmem,filearg,filein,fileout,nprcs,ovlp
     1               ,vtkf  ,vtkp   ,mefp  ,vtkbin
     2               ,vtkgr)
      implicit none
      include 'string.fi'
      character*80 filearg,filein,fileout
      integer*8 maxmem
      integer nin
      integer i,j
      integer nmacro,ierr
      logical ovlp,vtkf,vtkp,mefp,vtkbin,vtkgr
      integer nprcs
      parameter(nmacro = 12)
      character*15 rc,macro(nmacro)
      character*20 ex(11)
      character*30 string
      integer lnmacro
      logical lmacro(nmacro)
c ...
      data macro/'end            ','output         ','div            ',
     1           'method         ','partVtk        ','partMeshVtk    ',
     2           'partMeshMef    ','meshLoads      ','vtkBin         ',
     3           'input          ','memory         ','help           '/
c .....................................................................
c
c ... exemplo
      data ex/'input  mesh.dat'   ,'output     part','div          12',
     .        'method overllaping','partVtk     yes','partMeshVtk yes',
     .        'partMeshMef yes'   ,'meshLoads   yes','vtkBin       no',
     .        'memory     1000'   ,'end            '/
c .....................................................................
      lmacro(1:nmacro) = .false.
c ... default
      nprcs   = 2
      fileout = 'output_prepar'
      ovlp    = .false.
      vtkf    = .true.
      vtkp    = .false.
      mefp    = .true.
      vtkgr   = .false.
      vtkbin  = .false.
c ...
      nin = 1
      open(nin,file=filearg,status='old',err=15)
      goto 50
   15 continue
      print*,'File ',trim(filearg),' not found !'
      call finalize()
c .....................................................................
c
c ...
   50 continue
      call readmacro(nin,.true.)
      write(rc,'(15a)') (word(i),i=1,15)
      do j = 1, nmacro
         if (rc .eq. macro(j)) go to 100
      enddo
      goto 2500
c .....................................................................
c
c ...
  100 continue
      go to(200 ,300 , 400     !end        ,output     ,div
     1     ,500 ,600 , 700     !method     ,partVtk    ,partMeshVtk
     2     ,800 ,900 ,1000     !partMeshMef,meshLoads  ,vtkBin
     3     ,1100,1200,1300) j  !input      ,memory     ,help
c .....................................................................
c
c ... macro end
  200 continue
      goto 2000
c .....................................................................
c
c ... nome do arquivo de saida
  300 continue
      call readmacro(nin,.false.)
      write(fileout,'(80a)') (word(i),i=1,strl)
      print*,'Output file prefix: ',trim(fileout)
      lmacro(j) = .true.
      lnmacro = lnmacro + 1
      go to 50
c .....................................................................
c
c ... numero de divisoes
  400 continue
      lmacro(j) = .true.
      call readmacro(nin,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*) nprcs
      print*,'Number of partitions ',nprcs
      lnmacro = lnmacro + 1
      go to 50
c .....................................................................
c
c ... metodo de particionamento
  500 continue
      call readmacro(nin,.false.)
      write(rc,'(15a)') (word(i),i=1,15)
      if(rc .eq. 'overllaping    ') then
        print*,'Method: Overllaping.'
        ovlp = .true.
      else
        print*,'Method: Non-overllaping.'
        ovlp = .false.
      endif
c     print*,'metodo ',word,ovlp
      lmacro(j) = .true.
      lnmacro = lnmacro + 1
      go to 50
c .....................................................................
c
c ... particionamento da malha global
  600 continue
      call readmacro(nin,.false.)
      write(rc,'(15a)') (word(i),i=1,15)
      if(rc .eq. 'yes            ') then
        vtkf = .true.
      else
        vtkf = .false.
      endif
c     print*,'partVkt ',word
      lmacro(j) = .true.
      lnmacro = lnmacro + 1
      go to 50
c .....................................................................
c
c ... malha particionada formato vtk
  700 continue
      call readmacro(nin,.false.)
      write(rc,'(15a)') (word(i),i=1,15)
      if(rc .eq. 'yes            ') then
        vtkp = .true.
      else
        vtkp = .false.
      endif
c     print*,'partMeshVtk ',word
      lmacro(j) = .true.
      lnmacro = lnmacro + 1
      go to 50
c .....................................................................
c
c ... malha particionada formato mef
  800 continue
      call readmacro(nin,.false.)
      write(rc,'(15a)') (word(i),i=1,15)
      if(rc .eq. 'yes            ') then
        mefp = .true.
      else
        mefp = .false.
      endif
c     print*,'partMeshMef ',word
      lmacro(j) = .true.
      lnmacro = lnmacro + 1
      go to 50
c .....................................................................
c
c ... malha com os cargamentos
  900 continue
      call readmacro(nin,.false.)
      write(rc,'(15a)') (word(i),i=1,15)
      if(rc .eq. 'yes            ') then
        vtkgr= .true.
      else
        vtkgr= .false.
      endif
c     print*,'meshLoads ',word
      lmacro(j) = .true.
      lnmacro = lnmacro + 1
      go to 50
c .....................................................................
c
c ... arquivo vtk binario
 1000 continue
      call readmacro(nin,.false.)
      write(rc,'(15a)') (word(i),i=1,15)
      if(rc .eq. 'yes            ') then
        vtkbin= .true.
      else
        vtkbin= .false.
      endif
c     print*,'vtkBin ',word
      lmacro(j) = .true.
      lnmacro = lnmacro + 1
      go to 50
c .....................................................................
c
c ... arquivo de entrada de dados
 1100 continue
      call readmacro(nin,.false.)
      write(filein,'(80a)') (word(i),i=1,strl)
      print*,'Input file: ',trim(filein)
      lmacro(j) = .true.
      lnmacro = lnmacro + 1
      go to 50
c .....................................................................
c
c ... arquivo de entrada de dados
 1200 continue
      call readmacro(nin,.false.)
      write(string,'(15a)') (word(j),j=1,15)
      read(string,*,err = 2222,end = 2222) maxmem
      print*,'Memory(MB): ',maxmem
c ... convertendo de Mbytes para para numero de inteiros e 4 bytes
      maxmem = (maxmem*1024*1024)/4
      go to 50
c .....................................................................
c
c ... help
 1300 continue
      go to 2501
c .....................................................................
c
c ...
 2000 continue
      close(nin)
      return
 2222 continue
      print*,'*** Error reading config !',macro(i)
      call finalize()
 2500 continue
      print*,'*** Error reading ',trim(filearg),' file !'
 2501 continue
      write(*,'(2x,a)') '******************************'
      write(*,'(2x,a)') 'Example usage of file:'
      write(*,'(2x,a)') '------------------------------'
      do i = 1, 11
          write(*,'(2x,a)') ex(i)
      enddo
      write(*,'(2x,a)') '------------------------------'
      write(*,'(2x,a)') '******************************'
      call finalize()
      end
c *********************************************************************
