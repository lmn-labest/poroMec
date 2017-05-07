c *********************************************************************
c * CALL_METIS : chamando o metis                                     *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * i_np  - ponteiro para o part do nos                               *
c * i_ep  - ponteiro para o part do elementos                         *
c * i_ix0 - conectividade                                             *
c * i_nen - numero nos por elementos                                  *
c * nnode - numero de nos global                                      *
c * numel - numero de elementos                                       *
c * maxno - numero de maximo de nos por elemento                      *
c * mype  - numero do processo MPI                                    *
c * npes  - numero de processo do MPI                                 *
c * ndiv  - numeros de divisoes da malha                              *
c * MPIW  - comunicador do MPI                                        *
c * ----------------------------------------------------------------- *
c * Parametros de saida:                                              *
c * ----------------------------------------------------------------- *
c * i_np  - ponteiro para o part do nos preenchidos                   *
c * i_ep  - ponteiro para o part do elementos preenchidos             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine call_metis(i_np,i_ep,i_ix0,i_nen,nnode,numel,maxno
     .                     ,mype,npes,ndiv ,MPIW)
      implicit none
c ... metis
      integer*8 i_np,i_ep
c ... malha 
      integer nnode,numel,maxno
      integer*8 i_ix0,i_nen
      integer mype,npes,MPIW,ndiv
c =====================================================================
c
c === com parmetis e o metis
#ifdef PARMETIS
c ... para mais de 2 processos chamar o parmetis
      if(npes .gt. 1)then
        print*,"Parmetis" 
        call parmetis(i_np,i_ep,i_ix0,nnode,numel,maxno,mype,npes,ndiv
     .               ,MPIW)
c ... para mais de 2 processos chamar o metis
      else
        print*,"sequencial metis" 
        call metis(i_np,i_ep,i_ix0,i_nen,nnode,numel,maxno)
      endif  
c =====================================================================
c
c === com o metis
#else
      print*,"sequencial metis" 
      call metis(i_np,i_ep,i_ix0,i_nen,nnode,numel,maxno)
#endif
c =====================================================================
c
c === 
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * METIS : gera o grafo do porticionamento do metis                  *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * i_np  - ponteiro para o part do nos                               *
c * i_ep  - ponteiro para o part do elementos                         *
c * i_ix0 - conectividade                                             *
c * nnode - numero de nos global                                      *
c * numel - numero de elementos                                       *
c * maxno - numero de maximo de nos por elemento                      *
c * nen   - numero de nos por elemento                                *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine metis(i_np,i_ep,i_ix0,i_nen,nnode,numel,maxno)
c ===
      use Malloc
      implicit none
c ... metis
      integer*8 i_np,i_ep
c ... malha 
      integer nnode,numel,maxno
      integer*8 i_ix0,i_nen
c =====================================================================
c
c === alloca
      i_np  = alloc_4('np      ',  1,nnode)  
      i_ep  = alloc_4('ep      ',  1,numel)  
c =====================================================================
c
c ===
      call pdmetis(ia(i_np),ia(i_ep),ia(i_ix0),ia(i_nen)
     .            ,numel,nnode,maxno)
c =====================================================================
c
c ===
      i_ix0 = dealloc('ix0     ')
      i_nen = locate ('nen     ')
      i_ep  = locate ('ep      ')
      i_np  = locate ('np      ')
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * PDMETIS :                                                         *
c *********************************************************************
      subroutine pdmetis(np  ,ep  ,el  ,nen,numel,nnode,maxno)
c ===
      use Malloc
      implicit none
      include 'elementos.fi'
      include 'parallel.fi'
c ... metis
#ifdef METIS5
      integer, pointer, dimension(:) :: vsize    => NULL()
      integer, pointer, dimension(:) :: vwgt     => NULL()
      integer, pointer, dimension(:) :: tpwgts   => NULL()
c      integer, pointer, dimension(:) :: options  => NULL()
      integer ncommon
      integer options(40)
      integer esize
      integer*8 i_eind,i_eptr 
      integer i
#endif
      integer   np  (*),ep  (*)
      integer   edgecut,maxno
      integer   etype,nn,ne,nflag,nparts
c ... malha 
      integer   el  (maxno,*),nen(*)
      integer nnode,numel
c =====================================================================
c
c === metis 5
#ifdef METIS5
      ncommon = 1  
      if(ntria3(1)   .gt. 0) ncommon = 2
      if(nquad4(1)   .gt. 0) ncommon = 2
      if(ntetra4(1)  .gt. 0) ncommon = 3
      if(nhexa8(1)   .gt. 0) ncommon = 4
      if(ntetra10(1) .gt. 0) ncommon = 3
      if(nhexa20(1)  .gt. 0) ncommon = 4
c ... CSR para o grafo da malha
      i_eptr = alloc_4("eptr    ",  1,numel+1)
      call meshToCsrIa(nen,ia(i_eptr),numel,esize)
      i_eind = alloc_4("eind    ",  1,esize)
      call meshToCsrJa(el,nen,ia(i_eptr),ia(i_eind),numel,maxno)
c .....................................................................
c
c ... metis optons
      call METIS_SetDefaultOptions(options)
c ... numera FORTRAN
      options(18) = 1
c .....................................................................
c
c ... metis
      nn    = nnode
      ne    = numel
      nflag = 1
      nparts = nprcs
      print*,"Metis-5"
      call METIS_PartMeshDual(ne      ,nn     ,ia(i_eptr),ia(i_eind) 
     .                       ,vwgt    ,vsize  ,ncommon   ,nparts
     .                       ,tpwgts  ,options,edgecut   ,ep  
     .                       ,np)
c =====================================================================
c
c === metis 4
#else
c === Tipo da malha
      if(ntria3(1)  .gt. 0)etype = 1
      if(ntetra4(1) .gt. 0)etype = 2
      if(nhexa8(1)  .gt. 0)etype = 3
      if(nquad4(1)  .gt. 0)etype = 4
c ====================================================================
c
c === metis
c ... inteiros no metis 64 sao todos de 8 bytes
      nn     = nnode
      ne     = numel
      nflag  = 1
      nparts = nprcs
      print*,"Metis-4"
c     print*,"Metis",nn,ne,npart,nprcs
      call METIS_PartMeshDual(ne,nn,el  ,etype,nflag,nparts,edgecut
     .                       ,ep  ,np  )
      print*,"saindo Metis"
c ====================================================================
#endif
c
c ===
      return
      end
c =====================================================================
c
c *********************************************************************
c
c *********************************************************************
c * STRUTC_MERSH - particionamento para malha estruturadas(GAMBIARRAR)*
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine struct_mesh(el,x,np,ep,nnode,numel,nenv,nen,ndm,nprcs) 
      implicit none
      integer el(nen+1,*),np(*),ep(*)
      integer i,j,nprcs,numel,nenv,nen,ndm,nnode
      integer nx,ny,hx,hy
      real*8  t,dt,x(ndm,*),xmin,xmax,ymin,ymax,lx,ly
c....    Particiona a malha por coordenada
c
       xmin=1.d20
       xmax=-1.d20
       ymin=1.d20
       ymax=-1.d20
       do i=1,nnode
          if(x(1,i).lt.xmin)xmin=x(1,i)
          if(x(1,i).gt.xmax)xmax=x(1,i)
          if(x(2,i).lt.ymin)ymin=x(2,i)
          if(x(2,i).gt.ymax)ymax=x(2,i)
       enddo
       lx=xmax-xmin
       ly=ymax-ymin
 
       if(nprcs.eq.2)then
          nx=2
          ny=1
       elseif(nprcs.eq.4)then
          nx=2
          ny=2
       elseif(nprcs.eq.6)then
          nx=3
          ny=2
       elseif(nprcs.eq.8)then
          nx=4
          ny=2
       elseif(nprcs.eq.10)then
          nx=5
          ny=2
       elseif(nprcs.eq.16)then
          nx=4
          ny=4
       elseif(nprcs.eq.24)then
          nx=6
          ny=4
       elseif(nprcs.eq.32)then
          nx=8
          ny=4
       elseif(nprcs.eq.64)then
          nx=8
          ny=8
       elseif(nprcs.eq.80)then
          nx=10
          ny=8          
       else
         print*,'numero de divisoes invalida'
         stop
       endif
       do i=1,nnode
          hx=min(int((dfloat(nx)*(x(1,i)-xmin)/lx))+1,nx)
          hy=min(int((dfloat(ny)*(x(2,i)-ymin)/ly))+1,ny)
          np(i)=(hy-1)*nx+hx
       enddo
 
       do i=1,numel
          do j=1,nen
             ep(i)=max(ep(i),np(el(j,i)))
          enddo
       enddo
c ......................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
c * STRUTC_MERSH - particionamento para malha estruturadas(GAMBIARRAR)*
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine struct_cubo(el,x,np,ep,nnode,numel,nenv,nen,ndm,nprcs)
      implicit none
      integer el(nen+1,*),np(*),ep(*)
      integer i,j,k,kk,nprcs,numel,nenv,nen,ndm,nnode
      integer nx,ny,nz,nel
      real*8 hx(3),hx0(3),dx(3),lx(3),xm(3)
      real*8  t,dt,x(ndm,*),xmin,xmax,ymin,ymax,zmin,zmax
c ...    Particiona a malha por coordenada
c
      xmin=1.d20
      xmax=-1.d20
      ymin=1.d20
      ymax=-1.d20
      zmin=1.d20
      zmax=-1.d20
      do i=1,nnode
          if(x(1,i).lt.xmin)xmin=x(1,i)
          if(x(1,i).gt.xmax)xmax=x(1,i)
          if(x(2,i).lt.ymin)ymin=x(2,i)
          if(x(2,i).gt.ymax)ymax=x(2,i)
          if(x(3,i).lt.zmin)zmin=x(3,i)
          if(x(3,i).gt.zmax)zmax=x(3,i)
      enddo
      lx(1)=xmax-xmin
      lx(2)=ymax-ymin
      lx(3)=zmax-zmin
c ... 
      if(nprcs.eq.2)then
          nx=2
          ny=1
          nz=1
      elseif(nprcs.eq.4)then
          nx=2
          ny=2
          nz=1
      elseif(nprcs.eq.6)then
          nx=3
          ny=2
          nz=1
      elseif(nprcs.eq.8)then
          nx=2
          ny=2
          nz=2
      elseif(nprcs.eq.10)then
          nx=5
          ny=2
          nz=1
      elseif(nprcs.eq.16)then
          nx=4
          ny=2
          nz=2
      elseif(nprcs.eq.24)then
          nx=3
          ny=4
          nz=2
      elseif(nprcs.eq.32)then
          nx=4
          ny=4
          nz=2
      elseif(nprcs.eq.64)then
          nx=4
          ny=4
          nz=4
      elseif(nprcs.eq.80)then
          nx=4
          ny=4
          nz=5
      else
         print*,'numero de divisoes invalida'
         stop
      endif
c .....................................................................
c 
c ...
      kk      = 0
      dx(1:3) = 0.d0
      dx(1)   = lx(1)/dfloat(nx)
      dx(2)   = lx(2)/dfloat(ny)
      dx(3)   = lx(3)/dfloat(nz)
c .....................................................................
c 
c ...
      hx0(3) = zmin
      do i = 1, nz
        hx(3)  = hx0(3) + dx(3)
        hx0(2) = ymin
        do j = 1, ny 
          hx(2)  = hx0(2) + dx(2)
          hx0(1) = xmin
          do k = 1, nx 
            hx(1)  = hx0(1) + dx(1)
            kk = kk + 1
            do nel = 1, numel
              call ponto_medio(x,el(1,nel),xm,ndm,nenv)
              if( (       xm(1) .le. hx(1) 
     1              .and. xm(2) .le. hx(2) 
     2              .and. xm(3) .le. hx(3) ) 
     3              .and. 
     4             (      xm(1) .ge. hx0(1)
     5              .and. xm(2) .ge. hx0(2)
     6              .and. xm(3) .ge. hx0(3) ) ) then
                ep(nel) = kk          
              endif
            enddo
            hx0(1) = hx(1)
          enddo
          hx0(2) = hx(2)
        enddo
        hx0(3) = hx(3)
      enddo
c .....................................................................
c
c ...
      do nel = 1, numel
        do j = 1, nen
          kk = el(j,nel)
c ...
          if( np(kk) .eq. 0 ) then
            np(kk) = ep(nel)
c ... troca de particao no caso de se par
          else
            if( mod(kk,2) .eq. 0 ) np(kk) = ep(nel) 
          endif
        enddo
      enddo
c .....................................................................      
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MESHTOCSRIA : convertendo o grafo da malha para o formato csr     *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * nen   - no por elementos                                          *
c * eptr  - nao definido                                              *
c * numel - numero de elementos                                       *
c * maxno - numero maximo de nos por elemento                         *
c * esize - nao definido                                              *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * eptr  - malha no formato CSR ( ponteiro )                         *
c * esize - numero total de elementos no CSR                          *
c *********************************************************************
      subroutine meshToCsrIa(nen,eptr,numel,esize)
      implicit none
      integer nen(*)
      integer eptr(*)
      integer numel,esize
      integer i
c ...
      eptr(1) = 1
      do i = 2, numel + 1
        eptr(i) = eptr(i-1) + nen(i-1)
      enddo
      esize = eptr(numel+1) - eptr(1)
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MESHTOCSRIAJA : convertendo o grafo da malha para o formato csr   *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el    - conecitividade dos elementos                              *
c * nen   - no por elementos                                          *
c * eptr  - ponteiro CSR                                              *
c * eind  - nao definido                                              *
c * numel - numero de elementos                                       *
c * maxno - numero maximo de nos por elemento                         *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * eind  - malha no formato CSR ( conectividade)                     *
c *********************************************************************
      subroutine meshToCsrJa(el,nen,eptr,eind,numel,maxno)
      implicit none
      integer numel,maxno
      integer el(maxno,*),nen(*)
      integer eind(*),eptr(*)
      integer i,j
c ...
      do i = 1, numel
        do j = 1, nen(i)
           eind(eptr(i)+j-1) = el(j,i)
        enddo
      enddo
c .....................................................................
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MESHTOCSRIAJA : convertendo o grafo da malha para o formato csr   *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el    - conecitividade do elemento                                *
c * nen   - no por elementos                                          *
c * maxno - numero de nos por elemento                                *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * xm    - ponto medio geometrico                                    *
c ********************************************************************* 
      subroutine ponto_medio(x,el,xm,ndm,nen)
      implicit none
      real*8 x(ndm,*),xm(*)
      integer el(*),nen,j,k,ndm,no
c      
      xm(1:ndm) = 0.0
      if( nen .eq. 8 ) then
        do j = 1, nen
          no    = el(j)
          do k = 1, ndm
            xm(k) = xm(k) + x(k,no)
          enddo
        enddo
        xm(1:ndm) = xm(1:ndm)/nen
      endif
c      
      return
      end
c *********************************************************************
c
c *********************************************************************
#ifdef PARMETIS
c *********************************************************************
c * PARMETIS : prepar para chamar o Parmetis                          *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * i_np  - ponteiro para o part do nos                               *
c * i_ep  - ponteiro para o part do elementos                         *
c * i_ix0 - conectividade                                             *
c * nnode - numero de nos global                                      *
c * numel - numero de elementos                                       *
c * maxno - numero de nos por elemento                                *
c * mype  - numero do processo MPI                                    *
c * npes  - numero de processo do MPI                                 *
c * ndiv  - numeros de divisoes da malha                              *
c * MPIW  - comunicador do MPI                                        *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine parmetis(i_np,i_ep,i_ix0,nnode,numel,maxno,mype,npes
     .                   ,ndiv,MPIW)
      use Malloc
      implicit none
      include 'mpif.h'
c ... metis
      integer*8 i_np64,i_ep64,i_el64
      integer*8 i_np,i_ep
      integer*8 i_elmdist,i_my_elements,i_elmwgt
      integer*8 i_tpwgts,i_part,i_eptr
      integer*8 i_aux
      real ubvec(10)
      integer nel,ndiv,ncon,options(10)
      integer  edgecut,esizes,mgcnum
      integer  etype,nn,ne,nflag,npart,wflag
      data esizes /1/
      data mgcnum/1/
c ... malha 
      integer nnode,numel,maxno
      integer*8 i_ix0
      integer*8 i_nincid,i_incid
      integer maxgrade
c ... Mpi
      integer mype,npes,MPIW,ierr,i,j,k
c =====================================================================
c
c === gerando a distribuicao do elementos por processo (base do IO do
c     parmetis)
      i_elmdist = alloc_4("elmdist ",  1,npes+1)
      call split_elements(ia(i_ix0),numel,ia(i_elmdist),i_my_elements
     .                   ,npes,mype,esizes,mgcnum,nel,MPIW)
c ....................................................................
c
c ... pesos dos elementos
      i_elmwgt  = alloc_4("elmwgt  ",  1,numel)
      call mzero(ia(i_elmwgt),numel)
c ....................................................................
c
c ...
      ncon = 1
      i_tpwgts = alloc_4("tpwgts  ",  1,ndiv)
      call options_flag(ia(i_tpwgts),ubvec,options,ncon,wflag
     .                 ,nflag,ndiv)  
c ....................................................................
c
c ...
      i_part = alloc_4("part    ",  1,nel)
      i_eptr = alloc_4("eptr    ",  1,nel+1)
      do i = 0, nel
        ia(i_eptr+i) = esizes
      enddo
      call makecsr(nel,ia(i_eptr))
      ia(i_eptr+nel) = ia(i_eptr+nel) - 1
c ....................................................................
c      call MPI_Barrier(MPIW,ierr)
c      do j = 0 , npes
c        call MPI_Barrier(MPIW,ierr)
c        if( mype .eq. j ) then
c          do i = 0, nel - 1
c            print*,mype,"co",(ia(i_my_elements+i*esizes+k),k=0,esizes-1)
c          enddo  
c        endif
c      enddo
c      call MPI_Barrier(MPIW,ierr)
c      do j = 0 , npes
c        call MPI_Barrier(MPIW,ierr)
c        if( mype .eq. j ) then
c          print*,mype,"eptr",(ia(i_eptr+i),i=0,nel)
c        endif
c      enddo
c ===
      call parmetis_v3_partmeshkway(ia(i_elmdist),ia(i_eptr)
     .,ia(i_my_elements),ia(i_elmwgt)
     .,wflag,nflag,ncon,mgcnum,ndiv,ia(i_tpwgts),ubvec
     .,options,edgecut,ia(i_part),MPIW)
c =====================================================================
c
c === particionamento por elementos
      i_ep = alloc_4("ep      ",  1,numel)
      call globalep(ia(i_ep),ia(i_part),ia(i_elmdist),numel,nel,npes
     .             ,mype,MPIW)
c =====================================================================
c
c === particionamento dos nos
c ...  gerando as propriedades nodais
      if( mype .eq. 0 ) then
        call propnode(ia(i_ix0),i_nincid,i_incid,nnode,numel,maxno
     .               ,maxgrade,.false.)
c .....................................................................
      
        i_np = alloc_4("np      ",  1,nnode)
        i_aux= alloc_4("aux     ",  1,ndiv)
        call partnode(ia(i_np),ia(i_ep),ia(i_aux),ia(i_nincid)
     .               ,ia(i_incid),nnode,maxgrade,ndiv)  
        i_incid      = dealloc("incid   ")
        i_nincid     = dealloc("nincid  ")
        i_aux        = dealloc("aux     ")
        i_ix0        = dealloc("ix0     ")
        i_np         = locate("np      ")
c ... escravos        
      else  
        i_np = alloc_4("np      ",  1,nnode)
      endif  
c =====================================================================
c
      call MPI_Bcast(ia(i_np),nnode,MPI_INTEGER,0,MPIW,ierr)
c === liberando
      i_eptr       = dealloc("eptr    ")
      i_part       = dealloc("part    ")
      i_tpwgts     = dealloc("tpwgts  ")
      i_elmwgt     = dealloc("elmwgt  ")
      i_elmdist    = dealloc("elmdist ")
      i_my_elements= dealloc("my_el   ")
c =====================================================================
c
c === relocalizando      
      i_np     =  locate("np      ")
      i_ep     =  locate("ep      ")
c =====================================================================
c
c ===            
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * PARNODE : gera o partiocionamento de nos                          *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ep     - particionamento de elementos                             *
c * el      - conectividade                                           *
c * np      - indefinido                                              *
c * ep      - particionamento global de elementos                     *
c * nincid  - numero de incidencia de elemntos por no                 *
c * incid   - incidencia de elementos por no                          *
c * MPIW  - comunicador do MPI                                        *
c * ----------------------------------------------------------------- *
c * Parametros de saida:                                              *
c * ----------------------------------------------------------------- *
c * np     - particionamento global (todos os processo)               *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine partnode(np,ep,aux,nincid,incid,nnode,maxgrade
     .                   ,nprcs)  
      integer ep(*),np(*)
      integer maxgrade,nnode,nprcs
      integer nincid(*),incid(maxgrade,*)
      integer aux(*)
      integer i,j,npart
c ====================================================================
c
c ===
      do i = 1, nnode
        call mzero(aux,nprcs)
        do j = 1 , nincid(i)
              npart  = ep(incid(j,i))
          aux(npart) = aux(npart) +  1
        enddo
c ... no pertence a particao que ele paparece mais vezes
c     em caso de impate, no ira para a primeira particao
        npart = 0
        do j = 1, nprcs
          if(npart .lt. aux(j))then
            npart = aux(j)
            np(i) = j   
          endif
        enddo
      enddo
c ====================================================================
c
c ===
      return
      end
c ====================================================================
c *********************************************************************
c
c *********************************************************************
c * GLOBALEP: gera o partiocionamento de elementos global             *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ep     -                                                          *
c * part   - particionamento de elementos local                       *
c * elmdist- distribuicao de elementos por no                         *
c * numel  - numero elementos global                                  *
c * nel    - numero local de elementos                                *
c * mype  - numero do processo MPI                                    *
c * npes  - numero de processo do MPI                                 *
c * MPIW  - comunicador do MPI                                        *
c * ----------------------------------------------------------------- *
c * Parametros de saida:                                              *
c * ----------------------------------------------------------------- *
c * ep     - particionamento global (todos os processo)               *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine globalep(ep,part,elmdist,numel,nel,npes,mype,comm)
      implicit none
      include 'mpif.h'
      integer ep(*),part(*),elmdist(*)
      integer nel,numel
      integer npes,mype,comm
      integer i,pe,neld,nela
      integer status(MPI_STATUS_SIZE),ierr
c ====================================================================
c
c ===
      if (mype .eq. 0) then
c ... primeiros elementos pertencem ao root      
        do i = 1, nel
          ep(i) = part(i)
        enddo
        do pe = 2, npes
c ... numero de elementos no processo a ser recebido        
          nela = elmdist(pe+1) - elmdist(pe)
c ... numero de elementos ja recebidos
          neld = elmdist(pe)
          call MPI_Recv(ep(neld+1),nela,MPI_INTEGER,pe-1,1,comm
     .                 ,status,ierr);
        enddo
      else
        call MPI_Send(part,nel,MPI_INTEGER,0,1,comm,ierr);
      endif
      call MPI_Bcast(ep,numel,MPI_INTEGER,0,comm,ierr)
      do i = 1, numel
        ep(i) = ep(i) + 1
      enddo
c ====================================================================
c
c ===
      return
      end
c ====================================================================
c
c *********************************************************************
c * MAKECSR :                                                         *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * n      - tamanho do vetor                                         *
c * v      - vetor (vetor inteiro)                                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida:                                              *
c * ----------------------------------------------------------------- *
c * v      - no formato csr                                           *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine makecsr(n,v)
      implicit none
      integer n,v(*)
      integer i
c ====================================================================
c
c ===
      do i = 2, n
        v(i) = v(i) + v(i-1)
      enddo
      do i = n+1, 2 ,-1
        v(i) = v(i-1)
      enddo
      v(1) = 0
c ====================================================================
      return
      end
c ====================================================================
c
c *********************************************************************
c
c *********************************************************************
c * OPTIONS_FLAG: flags do parmetis                                   *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * tpwgts -                                                          *
c * ubvec  -                                                          *
c * options-                                                          *
c * ncon   -                                                          *
c * wgtflag-                                                          *
c * numflag-                                                          *
c * ndivs  - numero de divisoes da malha                              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine options_flag(tpwgts,ubvec,options,ncon,wgtflag
     .                       ,numflag,ndiv)
      implicit none 
      real tpwgts(*),ubvec(*)
      integer options(*),ncon,wgtflag,numflag,ndiv
      integer i
c .....................................................................
c
c ...
      do i = 1, ndiv*ncon
        tpwgts(i) = 1.0/ndiv
      enddo
c      
      do i = 1, ncon
        ubvec(i) = 1.05
      enddo
c
      options(1) = 1
      options(2) = 7
      options(3) = 0
c
      numflag = 0
      wgtflag = 0
c .....................................................................
c
      return
      end
c *********************************************************************
c
c *********************************************************************
c * SPLIT_ELEMENTS: distribui os elementos por processo               *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el     -conectividades sem o material                             *
c * numel  -numero de elementos                                       *
c * elmdist-                                                          *
c * i_my_el-                                                          *
c * npes   - numero de processos(MPI)                                 *
c * mype   - numero do processo(MPI)                                  *
c * esize  -                                                          *
c * nel    - numero do processo MPI                                   *
c * comm   - comunicador do MPI                                       *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine split_elements(el,numel,elmdist,i_my_el,npes,mype
     .                         ,esizes,mgcnums
     .                         ,nel,comm)  
      use Malloc
      implicit none
      include 'mpif.h'
      include 'elementos.fi'
c ... metis
      integer elmdist(*)
      integer*8 i_my_el,i_your_el
      integer*8 i_np,i_ep
      integer nel,etype,your_nel
      integer  esizes,mgcnums
      integer aesizes(5),amgcnums(5)
      data aesizes /-1,3,4,8,4/
      data amgcnums /-1,2,3,4,2/
c ... malha 
      integer numel
      integer  el(*)
c ... Mpi
      integer status(MPI_STATUS_SIZE)
      integer mype,npes,comm
c ...
      integer i,j,k,kk,pe,ierr 
c =====================================================================
c
c ===
c ...
      if(mype .eq. 0) then
        elmdist(1) = 0
        j = numel
        do i = 1, npes
                      k = j/(npes-i+1)
           elmdist(i+1) = elmdist(i) + k
                     j = j - k
        enddo
      endif
c .....................................................................
c
c ... comunicando a distribuicao
c
      call MPI_Bcast(elmdist,npes+1,MPI_INTEGER,0,comm,ierr)
c
c ... numero de elementos locais
c
      nel = elmdist(mype+2) - elmdist(mype+1)
c      
c ... Tipo da malha
      etype = 0
      if(ntria3(1)  .gt. 0)etype = 1
      if(ntetra4(1) .gt. 0)etype = 2
      if(nhexa8(1)  .gt. 0)etype = 3
      if(nquad4(1)  .gt. 0)etype = 4
c      
c ...
c
      esizes = aesizes(etype+1)
      mgcnums = amgcnums(etype+1)
c      
c ... alocando elementos por processos
c
      i_my_el = alloc_4("my_el   ",nel,esizes)
c
c ...
c
      if( mype .eq. 0) then
        i_your_el = alloc_4("your_el   ",nel,esizes)
        do pe = 1, npes
          your_nel = elmdist(pe+1) - elmdist(pe)
          do i = 0, your_nel - 1
            k  = esizes*(elmdist(pe)+i)
            kk = esizes*i
            do j = 0, esizes - 1
               ia(i_your_el+kk+j) = el(k+j+1)
            enddo
          enddo
c ... enviando para os escravos
          if ( pe .gt. 1) then
            call MPI_Send(ia(i_your_el),your_nel*(esizes),MPI_INTEGER
     .                    ,pe-1,0,comm,ierr)
          else
            do i = 0 , (your_nel)*esizes -1
               ia(i_my_el+i) = ia(i_your_el+i)
            enddo
          endif
        enddo
      i_your_el = dealloc("your_el ")
c ... recebendo do processo raiz
      else
        call MPI_Recv(ia(i_my_el),nel*(esizes),MPI_INTEGER,0,0,comm
     .               ,status,ierr);
      endif
c =====================================================================
c
c ===
      return
      end
c =====================================================================
#endif
