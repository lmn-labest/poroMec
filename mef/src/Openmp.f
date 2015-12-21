* ***************************svn****************************************
* $LastChangedDate:: 2011-09-29 17:14:54 -0300 (Thu, 29 Sep 2011)      $:*
* $Author:: henrique                                                 $:*
* $Rev:: 949                                                         $:*
* $Id:: Openmp.f 949 2011-09-29 20:14:54Z henrique                   $:* 
* **********************************************************************
c *********************************************************************
c * INIT_OPENMP: incia o openmp                                       *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * my_id       - identificacao do processo mpi                       *
c * num_threads - numero de threads                                   *
c * openmp      - openmp habilitado no codigo                         *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * OBS: se force for false o numero de threads maximo e disponivel   * 
c * ou definido pela variavel de ambiente OMP_NUM_THREADS             * 
c * Windows : set OMP_NUM_THREADS=2                                   * 
c * linux   : export OMP_NUM_THREADS=2  ( fedora,Centos)              * 
c * linux   : setenv OMP_NUM_THREADS=2  ( depende da bash do sistema) * 
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine init_openmp(num_threads,openmp,my_id)
      implicit none
      include 'omp_lib.h'
      integer my_id
      integer num_threads,max_num_threads_available
      logical openmp,forced
c ... checa se eh solicitado o uso do openmp sem o codigo ter sido
c     compilado com flag do opnemp    
      max_num_threads_available = 0                      
!$    max_num_threads_available      = omp_get_max_threads()
      if (max_num_threads_available .eq. 0 .and. openmp) then
        print*," Erro: Openmp disable in compile time,", 
     .         " but enable in run-time. "
        stop
      endif  
c .....................................................................
c
c ... numero de threads
      if (openmp) then
c ... deixa o sistema determinar o numero maximo de threads      
!$      num_threads = omp_get_max_threads()
      else
        num_threads = 1
      endif
c .....................................................................
c
c ... openmp in compile time and run time
      if((my_id.eq.0).and.openmp)then
        print*,"**************    PORO_MEC     ***********************"
        print*,"**************     OPENMP      ***********************"
        print*,'Max available threads:',num_threads
        print*,"******************************************************"
c ...  openmp disable in complile time     
      elseif(my_id.eq.0)then
        print*,"**************    PORO_MEC     ***********************"
        print*,"**************     OPENMP      ***********************"
        print*,"Openmp disable in compile time.                      "
        print*,"******************************************************"
      endif
c .....................................................................
c
c ...      
      return
      end
* **********************************************************************
