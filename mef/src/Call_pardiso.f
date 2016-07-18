c **********************************************************************
c * Data de criacao    : 30/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_PARDISO : chama o sover pardiso mkl                           *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * a(*)     - matriz A                                                *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   -                                                         *
c * mtype    - tipo da matriz                                          *
c *           -2 -> simetrico indefinido                               *
c *            2 -> simetrico definido positivo                        *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq)      - vetor solucao                                        *
c * a(*),b(neq) - inalterados                                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Arranjos jat,iat e kat são utilizados na retrosubstituizao do      *
c * solver iLDLt                                                       *
c **********************************************************************  
      subroutine call_mkl_pardiso(neq,ia,ja,a,b,x,z,mtype)
      implicit none
      include 'mpif.h'
      integer ia(*),ja(*),neq
      real*8 a(*),b(*),x(*),z(*)
      real*8 norm,xkx,mem
c ... variavel interna do mkl( 64 btis)
      integer*8 pt(64) 
c ...
      integer iparm(64),msglvl,error
c ...
      integer maxfct,mnum,mtype,phase,nrhs
      integer idum
      real*8 ddum
      real*8 time
c ...
      real*8 dot
c ...
      error  = 0
      maxfct = 1
      mnum   = 1
      nrhs   = 1
c ...
      msglvl = 0
c .....................................................................
c
c ...
      pt(1:64)    = 0
      iparm(1:64) = 0
c .....................................................................
c
c ... simetrico indefinido
      if( mtype .eq. -2) then
        iparm(1)  = 1 ! no solver default
        iparm(2)  = 2 ! fill-in reordering from METIS
        iparm(7)  = 2 ! numbers of iterative refinement steps
        iparm(10) = 8 ! perturbe the pivot elements with 1E-08
        iparm(21) = 1 ! Pivoting for symmetric indefinite matrices.                                   
        iparm(24) = 0 ! Parallel factorization control.
c .....................................................................
c
c ... simetrico definido positivo
      else if( mtype .eq. 2) then
        iparm(1)  = 1 ! no solver default
        iparm(2)  = 2 ! fill-in reordering from METIS
        iparm(7)  = 2 ! numbers of iterative refinement steps
        iparm(24) = 0 ! Parallel factorization control.
      endif         
c .....................................................................
c
c ...
      time = Mpi_Wtime()  
c .....................................................................
c
c ... 
      phase = 13        
      msglvl = 0
#if _MKL_
      call pardiso (pt  , maxfct, mnum, mtype, phase, neq, a, ia, ja,
     .              idum, nrhs  , iparm, msglvl, b, x, error)
#endif
      time = Mpi_Wtime() - time 
c .....................................................................
c
c ...
      mem = (max(iparm(15),iparm(16) + iparm(17)))/1024.d0
c .....................................................................
c 
c ... produto:  x*Kx
      call matvec_csr_sym_v3(neq,ia,ja,a,x,z,.true.)
      xkx = dot(x,z,neq)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq))
c ......................................................................
c
c ... Termination and release of memory
      phase = -1 ! release internal memory
      msglvl = 0
#if _MKL_
      call pardiso (pt, maxfct, mnum, mtype, phase,neq, ddum, idum,idum,
     .              idum, nrhs, iparm, msglvl, ddum, ddum, error)
#endif
c .....................................................................
c
c ...
      write(*,1100)neq,mem,xkx,norm,time
c .....................................................................
c
c ...
      return
c ======================================================================
 1100 format(' (PARDISO) solver:'/
     1 5x,'Number of equations  = ',i20/
     2 5x,'Memory (MB)          = ',f20.2/
     4 5x,'x * Kx               = ',d20.10/
     5 5x,'|| x ||              = ',d20.10/
     6 5x,'CPU time (s)         = ',f20.2/)
      end
c ***********************************************************************  
c
c ***********************************************************************
c     subroutine call_mkl_pardiso_d(neq,ia,ja,a,b,x)
c     implicit none
c     include 'mpif.h'
c     integer ia(*),ja(*),neq
c     real*8 a(*),b(*),x(*),norm
c ... variavel interna do mkl( 64 btis)
c     integer*8 pt(64) 
c ...
c     integer iparm(64),msglvl,error
c ...
c     integer maxfct,mnum,mtype,phase,nrhs
c     integer idum
c     real*8 ddum
c     real*8 time
c ...
c     real*8 dot
c ...
c     error  = 0
c     maxfct = 1
c     mnum   = 1
c     nrhs   = 1
c ...
c     msglvl = 0
c .....................................................................
c
c ... simetrico indefinido 
c     mtype = - 2
c ... simetrico definido positivo
c     mtype =   2 
c .....................................................................
c
c ...
c     pt(1:64)    = 0
c     iparm(1:64) = 0
c .....................................................................
c
c ...
c     iparm(1)  = 0 ! no solver default
c     iparm(2)  = 2 ! fill-in reordering from METIS
c     iparm(3)  = 0 ! numbers of processors
c     iparm(4)  = 1 ! no iterative-direct algorithm
c     iparm(8)  = 2 ! numbers of iterative refinement steps
c     iparm(10) = 13 ! perturbe the pivot elements with 1E-13
c     iparm(11) = 0 ! use nonsymmetric permutation and scaling MPS
c     iparm(18) = -1 ! Output: number of nonzeros in the factor LU
c     iparm(19) = -1 ! Output: Mflops for LU factorization
c     iparm(20) = 0 ! Output: Numbers of CG Iterations
c     iparm(21) = 0 !                                    
c     iparm(24) = 1 !                                     
c .....................................................................
c
c ...
c     time = Mpi_Wtime()  
c     phase = 11 ! only reordering and symbolic factorization
c     msglvl = 1
c     call pardiso(pt  , maxfct, mnum, mtype, phase,neq, a , ia , ja ,
c    .            idum , nrhs, iparm, msglvl,ddum  , ddum, error)
c     write(*,*) 'Reordering completed ... '
c     if (error .ne. 0) then
c        write(*,*) 'The following ERROR was detected: ', error
c        stop
c     endif 
c     write(*,*) 'Number of nonzeros in factors = ',iparm(18)
c     write(*,*) 'Number of factorization MFLOPS = ',iparm(19)
c .....................................................................
c
c... Factorization.
c     phase = 22 ! only factorization
c     msglvl = 1
c     call pardiso (pt, maxfct, mnum, mtype, phase, neq, a, ia, ja,
c    .             idum , nrhs, iparm, msglvl, ddum, ddum, error)
c     write(*,*) 'Factorization completed ... '
c     if (error .ne. 0) then
c        write(*,*) 'The following ERROR was detected: ', error
c        stop
c     endif
c     write(*,*) 'Number of perturbed pivots = ',iparm(14)
c     write(*,*) 'Number of factorization MFLOPS = ',iparm(15)
c     write(*,*) 'Number of factorization MFLOPS = ',iparm(16)
c     write(*,*) 'Number of factorization MFLOPS = ',iparm(17)
c .....................................................................
c
c ... Back substitution and iterative refinement
c     phase = 33        
c     msglvl = 1
c     call pardiso (pt  , maxfct, mnum, mtype, phase, neq, a, ia, ja,
c    .              idum, nrhs  , iparm, msglvl, b, x, error)
c     write(*,*) 'Solve completed ... '
c     time = Mpi_Wtime() - time 
c .....................................................................
c
c ... norm-2 = || x ||
c     norm = dsqrt(dot(x,x,neq))
c     print*,norm,time
c ......................................................................
c
c ... Termination and release of memory
c     phase = -1 ! release internal memory
c     msglvl = 0
c     call pardiso (pt, maxfct, mnum, mtype, phase,neq, ddum, idum,idum,
c    .              idum, nrhs, iparm, msglvl, ddum, ddum, error)
c .....................................................................
c     return
c     end
c ***********************************************************************                  