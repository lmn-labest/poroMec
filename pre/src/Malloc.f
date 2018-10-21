c **********************************************************************
c *                                                                    *
c *   MALLOC.F                                            31/08/2005   *
c *                                                                    *
c *   Este arquivo contem subrotinas para gerenciamento de memoria:    *
c *                                                                    *
c *   init_malloc                                                      *
c *   alloc_4                                                          *
c *   alloc_8                                                          *
c *   locate                                                           *
c *   dealloc                                                          *
c *   icopy                                                            *
c *   azero                                                            *
c *   mzero                                                            *
c *   maxtest                                                          *
c *                                                                    *
c **********************************************************************
      module Malloc
         integer, allocatable, dimension(:) :: ia
         integer*8, external :: alloc_4,alloc_8,locate,dealloc
c        integer*8, parameter :: maxmem = 10000000
         integer*8 maxmem /10000/        
      end module
      subroutine init_malloc()
c **********************************************************************
c *                                                                    *
c *   INIT_MALLOC: inicializa a estrutura de dados para as rotinas     *
c *   ------------ de gerenciamento de memoria                         *
c *                                                                    *
c *   Variaveis do COMMON/MALLOC/:                                     *
c *   ---------------------------                                      *
c *                                                                    *
c *   ip(maxnpts) - ponteiro do i-esimo arranjo                        *
c *   arname(maxnpts)   - nome do i-esimo arranjo                      *
c *   nalp        - numero de ponteiros alocados                       *
c *   ip(nalp+1)  - proximo ponteiro livre                             *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
      integer*8 ip(maxnpts),nalp
      character*8 arname(maxnpts)
      integer i,ierr
      common /malloc_info/ arname,ip,nalp
c ......................................................................
c
c ... Inicializacao da estrutura de dados de gerenciamento de memoria:
      allocate(ia(maxmem), stat=ierr)
      if (ierr .ne. 0) then
         print *, 'error: cannot allocate memory pool in heap.'
         call finalize(ierr)
      endif
c ......................................................................

c ...
      do i = 1, maxnpts
         ip(i)     = 0
         arname(i) = '        '
      enddo
      nalp = 0
      ip(nalp+1) = 1
      return
      end
      integer*8 function alloc_4(name,nl,nc)
c **********************************************************************
c *                                                                    *
c *   ALLOC_4: aloca memoria para arranjos de 4 bytes                  *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   name - nome do arranjo a ser alocado                             *
c *   nl   - numero de linhas do arranjo                               *
c *   nc   - numero de colunas do arranjo                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   alloc_4 - ponteiro do arranjo                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
      integer*8 ip(maxnpts),nalp,n,ipi,locate,aloc4
      character*8 arname(maxnpts),name
      common /malloc_info/ arname,ip,nalp
      integer nl,nc
c ......................................................................
      n = nl*nc
      if (n .le. 0) then
         aloc4 = ip(nalp+1)
         return
      endif
c ......................................................................      
      ipi = locate(name)
      if(ipi .gt. 0) then
         print*,'*** <ALLOC_4> Nome de arranjo existente: ',
     .      '(',name,') ***'
         call finalize()
      endif
c ......................................................................      
      nalp   = nalp + 1
      if (nalp+1 .gt. maxnpts) then
         print*,'*** <ALLOC_4> Max. numero de ponteiros excedido: ',
     .        '(',name,') ***'
         call finalize()
      endif
c ......................................................................      
      ipi = ip(nalp)
      ip(nalp+1) = ipi + n
      call maxtest(ip(nalp+1),name)
      arname(nalp) = name
      alloc_4 = ipi
c ......................................................................                  
      return
      end
      integer*8 function alloc_8(name,nl,nc)
c **********************************************************************
c *                                                                    *
c *   ALLOC_8: aloca memoria para arranjos de 8 bytes                  *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   name - nome do arranjo a ser alocado                             *
c *   nl   - numero de linhas do arranjo                               *
c *   nc   - numero de colunas                                         *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   alloc_8 - ponteiro do arranjo                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer maxnpts
      parameter  (maxnpts = 200)
      integer*8 ip(maxnpts),n,ipi,locate,aloc8,nalp
      character*8 arname(maxnpts),name
      common /malloc_info/ arname,ip,nalp
      integer nl,nc
c ......................................................................
      n = nl*nc
      if (n .le. 0) then
         aloc8 = ip(nalp+1)
         return
      endif      
c ......................................................................
      ipi = locate(name)
      if(ipi .gt. 0) then
         print*,'*** <ALLOC_8> Nome de arranjo existente: ',
     .        '(',name,') ***'
         call finalize()
      endif
c ......................................................................            
      nalp = nalp + 1
      if (nalp+1 .gt. maxnpts) then
         print*,'*** <ALLOC_8> Max. numero de ponteiros excedido: ',
     .        '(',name,') ***'
         call finalize()
      endif
c ......................................................................                  
      ipi = ip(nalp)
      if(mod(ipi,2) .eq. 0) then
         ipi      = ipi+1
         ip(nalp) = ipi
      endif
      ip(nalp+1) = ipi + 2*n
      call maxtest(ip(nalp+1),name)
      arname(nalp) = name
      alloc_8 = ipi
c ......................................................................                  
      return
      end
      integer*8 function locate(name)
c **********************************************************************
c *                                                                    *
c *   LOCATE: localiza o ponteiro do arranjo 'name'                    *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   name - nome do arranjo a ser localizado                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   locate - ponteiro do arranjo                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
      integer*8 ip(maxnpts),nalp
      character*8 arname(maxnpts),name
      common /malloc_info/ arname,ip,nalp
      integer i
c ......................................................................
      locate = 0      
      do i = 1, nalp
         if(name .eq. arname(i)) then
               locate = ip(i)
             return
         endif
      enddo
      return
      end
      integer*8 function dealloc(name)
c **********************************************************************
c *                                                                    *
c *   This routine allocates memory for a real*8 array                 *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
      integer i,j        
      integer*8 ip(maxnpts),nalp,ip1,npos0,npos
      character*8 arname(maxnpts),name
      common /malloc_info/ arname,ip,nalp
c ......................................................................
c
c ... Localiza a posicao i do ponteiro a ser removido:
c
      dealloc = 0
      do 100 i = 1, nalp
         if(name .eq. arname(i)) then
            goto 200
         endif
  100 continue
      print*,'*** <DEALLOC> Ponteiro nao encontrado: ','(',name,') ***'
      call finalize()
  200 continue
c
c ... Remove o arranjo da memoria:
c
      ip1   = ip(i+1)
      npos0 = ip1 - ip(i)
      npos  = npos0 - mod(npos0,2)
      do 210 j  = i, nalp-1
         ip(j)  = ip(j+1) - npos
         arname(j) = arname(j+1)
  210 continue
      call icopy(ip1,ip(nalp+1),ip(i))
      ip(nalp)   = ip(nalp+1) - npos
      ip(nalp+1) = 0
      arname(nalp)   = arname(nalp+1)
      arname(nalp+1) = '        '
      nalp = nalp - 1
      return
      end
      subroutine icopy(i1,i2,i3)
c **********************************************************************
c *                                                                    *
c *   MCOPY: desloca um vetor de inteiros na memoria.                  *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer*8 i1,i2,i3,i,j
c ...........................................
      j = i3
      do 100 i = i1, i2-1
         ia(j) = ia(i)
         j = j + 1
  100 continue
      return
      end
      subroutine azero(a,j)
c **********************************************************************
c *                                                                    *
c *   AZERO: zera as posicoes de 1 ate j do vetor a.                   *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8 a(*)
      integer j,k
c ......................................................................
      do 100 k = 1, j
         a(k) = 0.d0
  100 continue
      return
      end
      subroutine mzero(m,j)
c **********************************************************************
c *                                                                    *
c *   MZERO: zera as posicoes de 1 ate j do vetor m.                   *
c *                                                                    *
c **********************************************************************
      implicit none
      integer m(*)
      integer k,j
c ......................................................................
      do 100 k = 1, j
         m(k) = 0
  100 continue
      return
      end
      subroutine maxtest(ip,name)
c **********************************************************************
c *                                                                    *
c *   Verifica se ip > maxmem                                          *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer*8 ip
      character*8 name
c ......................................................................
      if (ip .gt. maxmem) then
         print*,'*** <MAXTEST> Memoria insuficiente: (',name,') ***'
c
         ip = ip - maxmem
         print*, ip,' posicoes'
         ip = ip/1024/1024
         print*, ip,' MB'
c
         call finalize()
      endif
      return
      end
c *********************************************************************
