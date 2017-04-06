c **********************************************************************
c * Data de criacao    : 28/03/2017                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * KONZEY_CARMAN : calcula a permeabilidade em funcao da porosidade   *
c * usando a relação Konzey-Carman                                     * 
c * -------------------------------------------------------------------* 
c * perm   - permebilidade inicial                                     *
c * poro   - porosidade atual                                          *
c * poro0  - porosidade inicial                                        *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * permeabilidade para porosidade atual                               *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c **********************************************************************
      real*8 function konzey_carman(perm0,poro,poro0)
      implicit none
      real*8 perm0,poro,poro0,tmp1,tmp2
      tmp1 = (poro**3)/((1.d0-poro)**2)   
      tmp2 = ((1.d0-poro0)**2)/(poro0**3) 
      konzey_carman = perm0*tmp1*tmp2
      return
      end
c **********************************************************************
c 
c **********************************************************************
c * Data de criacao    : 03/03/2017                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * KS_AND_MUS :  : calcula a permeabilidade em funcao da porosidade   *
c * usando a relação Konzey-Carman                                     * 
c * -------------------------------------------------------------------* 
c * perm   - permebilidade inicial                                     *
c * poro   - porosidade atual                                          *
c * poro0  - porosidade inicial                                        *
c * k      - modulo volumetrico homogenizado                           *
c * mu     - modulo de cisalhamento homogenizado                       *
c * ------------------------------------------------------------------ * 
c * k      - modulo volumetrico da parte solida                        *
c * mu     - modulo de cisalhamento da parte solida                    *
c * ------------------------------------------------------------------ * 
c * permeabilidade para porosidade atual                               *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c **********************************************************************
      subroutine ks_and_mus(ym,ps,poro,k,mu)
      implicit none
      real(kind=8) ym,ps,poro,k,mu,ha,hb,hc,a,b,c,mus,ks,delta
c ... calculo mus
      ha = 1.d0 - poro
      hb = 9.d0 +  6.d0*poro
      hc = 8.d0 + 12.d0*poro
c ... forma quadratica
      a = 32.d0*ha*ha
      b = 4.d0*ha*( 9.d0*k - mu*hc - 6.d0*poro*k)
      c = 3.d0*mu*hc*poro*k - 4.d0*k*mu*hb
c ... 
      delta = b*b - 4.0*a*c
      if (delta .gt. 0) then
        delta = dsqrt(delta)
        mus = 0.5d0*(delta - b)/a
        if( mus .lt. 0 ) then
          print*,'mus < 0 !!'
          call stop_mef()
        endif
      else
        print*,'Raizes imaginarias'
        call stop_mef()
      endif   
c .....................................................................
c
c ... calculo kus
      ks = 4.d0*k*mus/(4.d0*mus*ha - 3.d0*poro*k)
      if( ks .lt. 0 ) then
        call stop_mef()
      endif
c .....................................................................
c
c ...
      k  = ks
      mu = mus
c .....................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 05/04/2017                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * hashin_shtrikman : calcula das propriedades mecanicas usando as    *
c * relação de hashin-shtrikman                                        * 
c * -------------------------------------------------------------------* 
c * poro   - porosidade atual                                          *
c * ks     - modulo volumetrico da parte solida                        *
c * mus    - modulo de cisalhamento da parte solida                    *
c * k      - nao definido                                              *
c * mu     - nao definido                                              *
c * im_biot- nao definido                                              *
c * c_biot - nao definido                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * k      - modulo volumetrico homogenizado                           *
c * mu     - modulo de cisalhamento homogenizado                       *
c * im_biot- inverso do modulo de biot                                 *
c * c_biot - coeficiente de biot                                       *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c **********************************************************************
      subroutine hashin_shtrikman(poro  ,ks    ,mus
     1                           ,k     ,mu
     2                           ,im_biot,c_biot)     
      implicit none
      real(8) poro,ks,mus,k,mu,im_biot,c_biot,ha,hb,hc
c ...  
      ha = 1.d0 - poro
      hb = 9.d0 +  6.d0*poro
      hc = 8.d0 + 12.d0*poro
c ......................................................................
c
c ...
      k = (4.d0*ks*mus*ha)/(3.d0*ks*poro + 4.d0*mus)
c 
      mu = ( mus*ha*(9.0*ks + 8*mus) ) / (ks*hb + mus*hc)
c 
      c_biot = 1.d0 - k/ks
c
      im_biot = (c_biot - poro)/ks
c
      return
      end 
c **********************************************************************

