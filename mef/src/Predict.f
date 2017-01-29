c ======================================================================
c      Preditor-multicorretor para derivada de ordem 1.                |
c ====================================================================== 
      subroutine predict(nnode,ndf,id,u,v)
c **********************************************************************
c *                                                                    *
c *   PREDICT                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 u(ndf,*),v(ndf,*),b
c ......................................................................      
      b = dt*(1.d0-alfa)
      do 110 i = 1, nnode
      do 100 j = 1, ndf
         k = id(j,i)
         if(k .gt. 0) then         
            u(j,i) = u(j,i) + v(j,i)*b
            v(j,i) = 0.d0
         endif
  100 continue
  110 continue
c ......................................................................    
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine update(nnode,ndf,id,u,v,x)
c **********************************************************************
c *                                                                    *
c *   UPDATE                                                           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'      
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 u(ndf,*),v(ndf,*),x(*),b
c ......................................................................            
      b = alfa*dt
      do 110 i = 1, nnode
         do 100 j = 1, ndf
            k = id(j,i)
            if(k .gt. 0) then
                u(j,i) =  u(j,i) + x(k)*b
                v(j,i) =  v(j,i) + x(k)
            endif
  100    continue
  110 continue
c ......................................................................  
      return
      end
c *********************************************************************
c
c *********************************************************************
c ======================================================================
c      Preditor-multicorretor para derivada de ordem 2.                |
c ====================================================================== 
      subroutine predict2(nnode,ndf,id,u,v,a,dt,alfa,beta)
c **********************************************************************
c *                                                                    *
c *   PREDICT2                                                         *
c *                                                                    *
c **********************************************************************
      integer nnode,ndf,k
      integer id(ndf,*)
      real*8  u(ndf,*),v(ndf,*),a(ndf,*)
      real*8  dt,alfa,beta,c,d
c ......................................................................
      c = dt*dt*.5d0*(1.d0-2.d0*beta)
      d = dt*(1.d0-alfa)
      do 110 i = 1, nnode
        do 100 j = 1, ndf
          k = id(j,i)
          if (k .gt. 0) then
            u(j,i) = u(j,i) + dt*v(j,i) + c*a(j,i)
            v(j,i) = v(j,i) + d*a(j,i)
            a(j,i) = 0.d0
          endif
  100   continue
  110 continue
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine update2(nnode,ndf,id,u,v,a,x,dt,alfa,beta)
c **********************************************************************
c *                                                                    *
c *   UPDATE2                                                          *
c *                                                                    *
c **********************************************************************
      integer nnode,ndf,k
      integer id(ndf,*)
      real*8  u(ndf,*),v(ndf,*),a(ndf,*),x(*)
      real*8  dt,alfa,beta,c,d
c ......................................................................
      c = beta*dt*dt
      d = alfa*dt
      do 110 i = 1, nnode
        do 100 j = 1, ndf
          k = id(j,i)
          if (k .gt. 0) then
            a(j,i) = a(j,i) + x(k)
            u(j,i) = u(j,i) + c*x(k)
            v(j,i) = v(j,i) + d*x(k)
          endif
  100   continue
  110 continue
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine delta_predict_pm(nnode,ndf,u,u0,fnno)
c **********************************************************************
c * Data de criacao    : 70/00/0000                                    *
c * Data de modificaco : 08/11/2016                                    *
c * ------------------------------------------------------------------ *
c * DELTA_PREDIT_PM :  calculo o incremento inicial                    *                                                            *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * nnode - numero total de nos                                        *
c * ndf   - numero max. de graus de liberdade por no                   *
c * u     - deslocamento e pressoes do tempo n+1                       *
c * u0    - deslocamento e pressoes do tempo n                         *
c * fnno  - identifica dos nos de vertices ( 1 - vertice | 0 )         *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * u     - delta de deslocamentos e pressao                           *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Note que caso os deslocamentos ou pressao forem prescritos         *
c * com variacao temporal du inicial sera diferente de zero            * 
c **********************************************************************
      implicit none
      integer nnode,ndf,i,j,fnno(*)
      real*8 u(ndf,*),u0(ndf,*),b
c ......................................................................      
c
c ...
      do 110 i = 1, nnode
        do 100 j = 1, ndf - 1
          u(j,i) = u(j,i) - u0(j,i)
  100   continue
        if(fnno(i) .eq. 1 ) u(ndf,i) = u(ndf,i) - u0(ndf,i)  
  110 continue
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine delta_update_pm(nnode,ndf,id,du,x,fnno)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 08/11/2016                                    *
c * ------------------------------------------------------------------ *
c * DELTA_UPDATE_PM :  atualizacao os incrementos do graus de liberdade*
c * do poro                                                            *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * nnode - numero total de nos                                        *
c * ndf   - numero max. de graus de liberdade por no                   *
c * id    - numeracao nodal das equacoes                               *
c * du    - incremento deslocamentos e pressoes i                      *
c * x     - solucao do sistemas de equacoes                            *
c * fnno  - identifica dos nos de vertices ( 1 - vertice | 0 )         *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c  * du    - incremento deslocamentos e pressoes da iteracao i + 1     *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      integer nnode,nnodev,ndf,k
      integer id(ndf,*),fnno(*)
      real*8  du(ndf,*),x(*)
c ......................................................................
c
c ... atualizacao dos incrementos dos deslocamentos
      do 110 i = 1, nnode
        do 100 j = 1, ndf-1
          k = id(j,i)
          if (k .gt. 0) then
            du(j,i) =  du(j,i) + x(k)
          endif
  100   continue
  110 continue
c ......................................................................
c
c ... atualizao da deltap
      do 120 i = 1, nnode
        if(fnno(i) .eq. 1 ) then
          k = id(ndf,i)
          if (k .gt. 0) then
            du(ndf,i) = du(ndf,i) + x(k)
          endif
        endif
  120 continue
c ......................................................................
c
c ...
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine update_mec(nnode,ndf,id,u,x)
c **********************************************************************
c *                                                                    *
c *   UPDATE_MEC : atualizacao dos graus de liberdade do mecanico      *
c *   estatico                                                         *
c *                                                                    *
c **********************************************************************
      integer nnode,nnodev,ndf,k
      integer id(ndf,*)
      real*8  u(ndf,*),x(*)
c ......................................................................
c
c ... atualizacao dos incrementos dos deslocamentos
      do 110 i = 1, nnode
        do 100 j = 1, ndf
          k = id(j,i)
          if (k .gt. 0) then
            u(j,i) = x(k)
          endif
  100   continue
  110 continue
c ......................................................................
c
c ...
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine update_plastic(ix,e,plastic,nen,npi,numel)
c **********************************************************************
c * Data de criacao    : 11/01/2016                                    *
c * Data de modificaco : 20/01/2017                                    *
c * ------------------------------------------------------------------ *
c * UPDATE_PLASTIC : atualizacao da dilatacao volumetrica              *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * plastic - parametros do plastico                                   *
c *           plastic(1) - dilatacao volumetrica plastica              *  
c *                        do passo anterior                           *
c *           plastic(2) - dilatacao volumetrica plastica              *
c *           plastic(3) - paramentro de endurecimento                 *  
c * npi     - numero de pontos de integracao                           *
c * numel   - numero de elementos                                      *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c  * plastic(1) - atualizado com a dilatacao volumetrica plastica      *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
      include 'termprop.fi'
      integer i,j,npi,numel,ix(nen+1,*),nen,ma
      real*8 plastic(3,npi,numel),alfa,e(prop,*),de
      real*8 k_plastic,e0,lambda_plastic
      logical flag_def,flag_pc
c ......................................................................
c
c ... atualizacao dos incrementos dos deslocamentos
      do i = 1, numel
        ma             = ix(nen+1,i)
        e0             = e(8,ma) 
        lambda_plastic = e(9,ma)
        k_plastic      = e(10,ma)
        alfa           = (1+e0)/(lambda_plastic-k_plastic)
        do j = 1, npi                  
          de             = plastic(2,j,i) - plastic(1,j,i) 
          plastic(1,j,i) = plastic(2,j,i)
          plastic(3,j,i) = plastic(3,j,i)*dexp(-alfa*de)  
        enddo
      enddo
c ......................................................................
c
c ...
      return
      end
c *********************************************************************
