c *********************************************************************
c * GETIME : retorna o tempo                                          * 
c *                                                                   * 
c *********************************************************************
      real*8 function getime()
      implicit none
c === with MPI 
#ifdef MPI 
      include 'mpif.h'
c ===      
      getime = Mpi_Wtime()
c =====================================================================      
c
c === without MPI 
#else
      real etime
      real tarray(2)
c ===      
      getime = etime(tarray)
c =====================================================================      
c
#endif
c
c ===
      return
      end
c =====================================================================      
 
      
