c **********************************************************************
c * Data de criacao    : 07/02/2017                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *      
c * CAL_RESIDUO: calculo do residus e dos criteiros de parada          *  
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * id     - numeracao nodal das equacoes                              *
c * b      - residuo                                                   *
c * du     - incremento da solucao                                     *
c * nno_pload - numero de nos ( mpi nno1+nno2)                         *
c * ndf       -graus de liberade                                       *
c * neq_dot   - numero de equacoes (mpi neq1+neq2)                     *
c * nli       - iteracao na linear                                     *
c * istop     - nao definido                                           *
c * type_crit - tipo de criterio de paraa                              *
c                  1 - |R|/|R0| < tol                                  *
c                  2 - |Ru|/|Ru0| < tol                                *
c                      |Rp|/|Rp0| < tol                                *
c                  3 - |du|/|du0| < tol                                *
c                      |dp|/|dp0| < tol                                *
c * tol        - tolerancia                                            *
c * my_id      -                                                       *
c * mpi        -                                                       *
c * nout       - arquivo de log do nao linear                          *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * istop - 2 criterios de convergencia satisfeitos                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************      
       subroutine cal_residuo(id        ,b         ,du   
     1                      ,nno_pload ,ndf       ,neq_dot ,nli  
     2                      ,istop     ,type_crit ,tol
     3                      ,my_id     ,mpi       ,nout)
      implicit none
      include 'mpif.h'
      integer nout
      integer i,j,neqi,istop,my_id,nli,type_crit
      integer id(ndf,*),nno_pload,ndf,neq_dot,ierr
      real*8 b(*),du(*),resu,resp,gresu,gresp,resu0,resp0,res0,res
      real*8 tol
      real*8 dot_par
      save resu0,resp0,res0
      logical mpi
      external dot_par
c ...
      if(type_crit .eq. 1) then
        res = dsqrt(dot_par(b,b,neq_dot))
        if(nli .eq. 1) res0 = res
        if(my_id.eq.0) then
          write(*,'(1x,a,i9)')'nonlinear iteration ',nli
          write(*,'(2(1x,a,e13.6))')'resid/resid0',res/res0,'resid',res
          write(nout,'(i7,4e20.10)')nli,res/res0,res   
        endif
        istop = 0
        if( res/res0 .lt. tol) istop = istop + 2
c .....................................................................
c
c ...
      elseif(type_crit .eq. 2) then
c ...
        resu = 0.0d0
        do i = 1, nno_pload
          do j = 1, ndf - 1
            neqi = id(j,i)          
            if( neqi .ne. 0 )then
              resu = resu + b(neqi)*b(neqi) 
            endif
          enddo
        enddo
c .....................................................................
c
c ...
        resp = 0.0d0
        do i = 1, nno_pload
          neqi = id(ndf,i)
          if( neqi .ne. 0 )then
            resp = resp + b(neqi)*b(neqi) 
          endif
        enddo
c .....................................................................
c
c ...
        if(mpi) then
          call MPI_ALLREDUCE(resu,gresu,1,MPI_REAL8
     .                      ,MPI_SUM,MPI_COMM_WORLD,ierr)
c
          call MPI_ALLREDUCE(resp,gresp,1,MPI_REAL8
     .                      ,MPI_SUM,MPI_COMM_WORLD,ierr)
        else
          gresu = resu
          gresp = resp
        endif
c .....................................................................
c
c ...
        gresu = dsqrt(gresu)
        gresp = dsqrt(gresp)
        if(nli .eq. 1) then
          resu0 = gresu
          resp0 = gresp
        else if(nli .eq. 2) then 
          resu0 = max(resu0,gresu)
          resp0 = max(resu0,gresu)
        endif
c
        if(my_id.eq.0) then
          write(*,'(1x,a,i9)')'nonlinear iteration ',nli
          write(*,'(2(1x,a,e13.6))')
     .    'resid(u)/resid0(u)',gresu/resu0,'resid(u)',gresu
          write(*,'(2(1x,a,e13.6))')
     .    'resid(p)/resid0(p)',gresp/resp0,'resid(p)',gresp
          write(nout,'(i7,4e20.10)')nli
     .         ,gresu/resu0,gresu,gresp/resp0,gresp   

        endif
c .....................................................................
c 
c ...
        istop = 0
        if( gresu/resu0 .lt. tol) istop = istop + 1
        if( gresp/resp0 .lt. tol) istop = istop + 1
c .....................................................................
c
c ...
      elseif(type_crit .eq. 3) then
c ...
        resu = 0.0d0
        do i = 1, nno_pload
          do j = 1, ndf - 1
            neqi = id(j,i)          
            if( neqi .ne. 0 )then
              resu = resu + du(neqi)*du(neqi) 
            endif
          enddo
        enddo
c .....................................................................
c
c ...
        resp = 0.0d0
        do i = 1, nno_pload
          neqi = id(ndf,i)
          if( neqi .ne. 0 )then
            resp = resp + du(neqi)*du(neqi) 
          endif
        enddo
c ...............................................................
        if(mpi) then
          call MPI_ALLREDUCE(resu,gresu,1,MPI_REAL8
     .                      ,MPI_SUM,MPI_COMM_WORLD,ierr)
 
          call MPI_ALLREDUCE(resp,gresp,1,MPI_REAL8
     .                      ,MPI_SUM,MPI_COMM_WORLD,ierr)
        else
          gresu = resu
          gresp = resp
        endif
c .....................................................................
c
c ...
        gresu = dsqrt(gresu)
        gresp = dsqrt(gresp)
        if(nli .eq. 1) then
          resu0 = gresu
          resp0 = gresp
        endif
c
        if(my_id.eq.0) then
          write(*,'(1x,a,i9)')'nonlinear iteration ',nli
          write(*,'(2(1x,a,e13.6))')
     .    'du/du0',gresu/resu0,'du',gresu
          write(*,'(2(1x,a,e13.6))')
     .    'rp/rp0',gresp/resp0,'rp',gresp
          write(nout,'(i7,4e20.10)')nli
     .         ,gresu/resu0,gresu,gresp/resp0,gresp   
        endif
c .....................................................................
c 
c ...
        istop = 0
        if( gresu/resu0 .lt. tol) istop = istop + 1
        if( gresp/resp0 .lt. tol) istop = istop + 1
c .....................................................................
      endif
c .....................................................................
c 
c ...        
      return
      end
