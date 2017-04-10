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
      real(kind=8) function konzey_carman(perm0,poro,poro0)
      implicit none
      real(kind=8) perm0,poro,poro0,tmp1,tmp2
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
      real(kind=8) poro,ks,mus,k,mu,im_biot,c_biot,ha,hb,hc
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
c
c **********************************************************************
c * Data de criacao    : 09/04/2017                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *
c * prop_porosity : calculo das propriedades mecanicas usando as       *
c * porosidade                                                         *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * b      - nao definido                                              *
c * c      - nao definido                                              *
c * ibiot  - nao definido                                              *
c * cbiot  - nao definido                                              *
c * vprop(*) -                                                         *
c *           1 - prop variavel                 (true|false)           *
c *           2 - konzey-Caraman                (true|false)           *
c *           3 - massa especifica homogenizada (true|false)           *
c *           4 - mecanico                      (true|false)           *
c * fluid_sw - peso especifico da agua                                 *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * a, b, c- matriz constitutiva                                       *
c * perm   - permeabilidade                                            *
c * im_biot- inverso do modulo de biot                                 *
c * c_biot - coeficiente de biot                                       *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c **********************************************************************
      subroutine prop_porosity(a,b,c,ibiot,cbiot,perm,vprop,fluid_sw)
      implicit none
      real(kind=8) a,b,c,ibiot,cbiot,perm,modK,mu,fluid_sw
      real(kind=8) vprop(*)
c ...
      real(kind=8) div23,div43
      parameter (div43 = 0.133333333333333d+01)
      parameter (div23 = 0.666666666666667d0)
c ...
      modK = vprop(4)
      mu   = vprop(5)
      ibiot = vprop(6)
      cbiot = vprop(7)
c ....................................................................
c
c ... matriz contitutiva E(modK,mu)
      a = modK  + div43*mu
      b = (modK - div23*mu)/a 
      c = mu/a 
c ....................................................................
c
c ... k = k(porosity)           
      perm = vprop(2)/fluid_sw
c ....................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 09/04/2017                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *
c * prop_porosity_plastic : calculo das propriedades plastica com      *
c * porosidade                                                         *
c * -------------------------------------------------------------------*
c * alpha  - nao definido                                              *
c * c14    - nao definido                                              *
c * g11    - nao definido                                              *
c * ps     - nao definido                                              *
c * vprop(*) -                                                         *
c *           1 - prop variavel                 (true|false)           *
c *           2 - konzey-Caraman                (true|false)           *
c *           3 - massa especifica homogenizada (true|false)           *
c *           4 - mecanico                      (true|false)           *
c * lam    - lamdba plastico                                           *
c * k      - k plastico                                                *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * alpha  -                                                           *
c * c14    -                                                           *
c * g11    -                                                           *
c * ps     - coeficiente de poisson                                    *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c **********************************************************************      
      subroutine prop_porosity_plastic(alpha,c14,g11,ps,vprop,lam,k)
      implicit none
      real(kind=8) alpha_exp,c14,g11,ps,lam,k,poro,modK,mu,alpha
      real(kind=8) vprop(*)
c ...
      real(kind=8) div23,div43
      parameter (div43 = 0.133333333333333d+01)
      parameter (div23 = 0.666666666666667d0)
c ...
      poro = vprop(1)
      modK = vprop(4)
      mu   = vprop(5)
c ....................................................................
c
c ...
      ps  = (3.d0*modK - 2.d0*mu) / (6.d0*modK + 2.d0*mu) 
      c14 = ps / ( 1.d0 - 2.d0*ps )
      g11 = mu
c
      alpha = 1.d0/( (1.d0 - poro)*(lam - k) )
c ......................................................................
c
c ...
      return
      end
c **********************************************************************
       