      subroutine plasticity3d_pm(eps1 ,eps2 ,tx
     1                          ,mcs  ,alfa ,pr 
     2                          ,c14  ,g    ,sup       
     3                          ,nel)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 17/02/2017                                    *
c * ------------------------------------------------------------------ *
c * PLASTICST3D_PM: calcula tensoes 3D, considerando plasticidade      *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * eps1         - armazena deformacao plastica                        *
c *                volumetrica do passo anterior                       *
c * eps2         - armazena pc do passo anterior                       *
c * tx(6)        - delta tensao efetiva                                *
c * alfa         - camclay                                             *
c * c14,g,pr     -  coeficientes da matriz constitutiva                *
c * sup          - antingui a superficie true| false                   *
c * nel          - numero do elemento                                  *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * tx(1:6)      - tensao efetiva atualizada (Total)                   *
c * eps1         - atualizado                                          *
c * sup          - true|false                                          *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
      integer i,iyied,j,maxit,nel
      real*8 tx(*),eps1,eps2,epe1,epe2,g
      real*8 snphi,pr,se1,steff,varj2,theta,sint3,c14,devs(6)
      real*8 stot(6),abeta,dlamb,hl1,escur1,tol,a1(6),d1(6)
      real*8 lamb,mcs,pc,sm,alfa
      real*8 ds
      logical sup
c ......................................................................
c   
c ...    
      sup   = .true.
      maxit = 500000
      tol   = 1.0e-07
c ... t0 = t:
      call aequalb(stot,tx,6)
c ...  eps1 = armazena deformacao plastica volumetrica do passo anterior      
c ...  eps2 = armazena pc do passo anterior
      epe1 = eps1
      epe2 = eps2
      pc   = epe2
c ... evaluate effective stress
      call effst3d(se1,iyied,stot,devs,steff,sm,varj2,mcs,escur1,pc) 
c ... check for yielding durint this iterations
      if(escur1 .le. 1.0e-10) then
        sup = .false.        
        return
      endif
c ... loop Newton-Raphson
      do 8 i = 1, maxit
c ... calculate vectors A and D
         call effst3d(se1,iyied,stot,devs,steff,sm,varj2,mcs,escur1,pc) 
         print*,i,dabs(escur1)
         if( dabs(escur1) .lt. tol)  goto 11
         call flow3d(a1,d1,devs,abeta,steff,theta,varj2,sint3,hl1,
     .          snphi,c14,g,pr,mcs,pc,sm,epe2,alfa)
c ... compute plastic multiplier
         dlamb = escur1 *abeta
c ... update plastic strain trace
         do j = 1, 3            
           epe1 = epe1 + a1(j)*dlamb
         enddo         
c ... update hardening parameter
         ds = epe1-eps1
         pc = epe2 - alfa*epe2*ds
c ...compute elastoplastic stresses
         do 10 j = 1,6
           stot(j) = stot(j)-dlamb*d1(j)
   10    continue
    8  continue
       write(*,*) 'numero maximo de iteracoes excedido - plasticidade'
       write(*,*) 'Elemento: ',nel
       write(*,*) 'it:',maxit,'res:',dabs(escur1),'tol:',tol
       call stop_mef()
   11 continue
c ... calculate equivalent plastic strain
      tx(1:6) = stot(1:6)
c .....................................................................
      eps1 = epe1
      eps2 = epe2   
      return
      end
      subroutine effst3d(se,iyied,tx,devs,steff,sm,varj2,mcs,escur,pc)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   Effst: calculates the deviatoric stresses, effective stress and  *
c *          invariants                                                *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer iyied
      real*8 root3,sm,tx(*),devs(*),varj2,steff,se
      real*8 mcs,escur,pc
      parameter ( root3   = 1.73205080757d0 )      
c ... sm = s11+s22+s33/3 = tr(sigma)/3
      sm      = (tx(1)+tx(2)+tx(3))/3.d0
c ... sij = sigmaij - sm*deltaij
      devs(1) = tx(1)-sm
      devs(2) = tx(2)-sm
      devs(3) = tx(3)-sm
      devs(4) = tx(4)
      devs(5) = tx(5)
      devs(6) = tx(6)
c ... J2 = 0.5*sij*sji 
      varj2 = devs(4)*devs(4)+devs(5)*devs(5)+devs(6)*devs(6)+
     . 0.5d0*(devs(1)*devs(1)+devs(2)*devs(2)+devs(3)*devs(3))
c ... J3 = det(sigma)
c     varj3 = devs(1)*devs(2)*devs(3)+2.d0*devs(4)*devs(5)*devs(6)
c    .        -devs(1)*devs(5)*devs(5)-devs(2)*devs(6)*devs(6)
c    .        -devs(3)*devs(4)*devs(4)
c ... meam  effective stress = steff = sqrt(s:s/2) = sqrt(sij*sji/2)
      steff = dsqrt(varj2)
c ... Von Mises effective stress sqrt(3(sij*sji)/2)
      se = root3*steff
c ... Modified Cam-Clays
      escur = se*se + mcs*mcs*sm*(sm+pc)
      return
      end
      subroutine flow3d(a,d,devs,abeta,steff,theta,varj2,sint3,hl,snphi,
     .                  c14,g,pr,mcs,pc,sm,epe2,alfa)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   Effst: calculates vectors A and D                                *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer j
      real*8 a1(6),a2(6),a(6),d(6),devs(6),root3,steff,varj2
      real*8 cons1,cons2,theta,snphi
      real*8 c14,xml,g,abeta,sint3,hl,denom,pr,mcs
      real*8 pc,sm,epe2,alfa
      root3 = 1.73205080757d0
c ... vector A1:
      a1(1) = 1.d0
      a1(2) = 1.d0
      a1(3) = 1.d0
      a1(4) = 0.d0
      a1(5) = 0.d0
      a1(6) = 0.d0
c ... vector A2:
      a2(1) = devs(1)!/(2.d0*steff)
      a2(2) = devs(2)!/(2.d0*steff)
      a2(3) = devs(3)!/(2.d0*steff)
      a2(4) = devs(4)!/steff
      a2(5) = devs(5)!/steff
      a2(6) = devs(6)!/steff
c ... vector A3
c     a3(1) =(2.d0/3.d0)*(devs(2)*devs(3)-devs(5)*devs(5))+(1.d0/3.d0)*
c    .       (devs(6)*devs(6)+devs(4)*devs(4)-devs(1)*(devs(3)+devs(2)))
c     a3(2) =(2.d0/3.d0)*(devs(1)*devs(3)-devs(6)*devs(6))+(1.d0/3.d0)*
c    .       (devs(5)*devs(5)+devs(4)*devs(4)-devs(2)*(devs(1)+devs(3)))
c     a3(3) =(2.d0/3.d0)*(devs(1)*devs(2)-devs(4)*devs(4))+(1.d0/3.d0)*
c    .       (devs(5)*devs(5)+devs(6)*devs(6)-devs(3)*(devs(1)+devs(2)))
c     a3(4) = 2.d0*(devs(5)*devs(6)-devs(3)*devs(4))
c     a3(5) = 2.d0*(devs(4)*devs(6)-devs(1)*devs(5))
c     a3(6) = 2.d0*(devs(4)*devs(5)-devs(2)*devs(6))
c ... Modified Cam-clay
      cons1 = (2.d0/3.d0)*mcs*mcs*sm+(1.d0/3.d0)*mcs*mcs*pc
      cons2 = 3.d0
c     cons3 = 0.d0
c ... vector A
      do 50 j = 1,6
c       a(j) = cons1*a1(j)+cons2*a2(j)+cons3*a3(j)
        a(j) = cons1*a1(j)+cons2*a2(j)
   50 continue
c ... vector D
      xml  = c14*(a(1)+a(2)+a(3))
      d(1) = 2.d0*g*(a(1)+xml)
      d(2) = 2.d0*g*(a(2)+xml)  
      d(3) = 2.d0*g*(a(3)+xml)
      d(4) = g*a(4)
      d(5) = g*a(5)
      d(6) = g*a(6)
c ... compute hardening parameter
      hl = mcs*mcs*sm*alfa*epe2
c ... compute part of plastic multiplier
      denom = hl*(a(1)+a(2)+a(3))
      do 80  j = 1,6
        denom = denom + a(j)*d(j)
   80 continue
      abeta = 1.d0/denom
      return
      end     
      subroutine principal_stress(t,ntn)
c **********************************************************************
c *                                                                    *
c *   PSTRESS3D: tensoes principais                                    *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    t(ntn,nnode)  - tensoes nodais referenciadas a xyz              *
c *    nnode - numero de nos                                           *
c *    ntn   - numero de tensoes por no                                *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *    t(ntn,nnode)  - tensoes principais                              *
c *    v1(3,nnode) - direcao principal 1                               *
c *    v2(3,nnode) - direcao principal 2                               *
c *    v3(3,nnode) - direcao principal 3                               *
c *                                                                    *
c **********************************************************************
      implicit none
      integer ntn,ip(3)
      real*8  t(ntn)
      real*8  a(3,3),x(3,3),w(3)
c ......................................................................
      a(1,1) = t(1)
      a(1,2) = t(4)
      a(1,3) = 0.d0
      a(2,1) = a(1,2)
      a(2,2) = t(2)
      a(2,3) = 0.d0
      a(3,1) = a(1,3)
      a(3,2) = a(2,3)
      a(3,3) = t(3)
      if(ntn .eq. 6) then
        a(1,3) = t(6)
        a(2,3) = t(5)
        a(3,1) = a(1,3)
        a(3,2) = a(2,3)
      endif
      call jacobi(a,x,w,ip,3,1.d-12,100)
      t(1) = w(1)
      t(2) = w(2)
      t(3) = w(3)
   20 continue
      return
      end
      subroutine radian(phi)
c **********************************************************************
c *                                                                    *
c *   RADIAN : converte angulo para radianos                           *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    phi - angulo em graus                                           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *    phi - angulo em radianos                                        *
c *                                                                    *
c **********************************************************************      
      implicit none
      real*8 phi, pi, temp
      pi = 4*datan(1.d0)
      temp = phi*pi/180.d0
      phi = temp
      return
      end
c **********************************************************************
      real*8 function cam_clay_pc(tx,mcs)
c **********************************************************************
c * Data de criacao    : 24/02/2017                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *
c * CAM_CLAY_PC: calcula pc para campo tensoes para que F(tx,pc) = 0   *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * tx           - tensao efetiva                                      *
c * mcs          - paramentro do camclay                               *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************      
      implicit none
      real*8 tx(6),devs(6),mcs,sm,q2,zero
      parameter (zero = 1.d-14)

c ... sm = s11+s22+s33/3 = tr(sigma)/3
      sm      = (tx(1)+tx(2)+tx(3))/3.d0
c ... sij = sigmaij - sm*deltaij
      devs(1) = tx(1)-sm
      devs(2) = tx(2)-sm
      devs(3) = tx(3)-sm
      devs(4) = tx(4)
      devs(5) = tx(5)
      devs(6) = tx(6)
c .....................................................................
c 
c ... q = 3/2 * s:s
      q2 = 1.5d0*(devs(1)*devs(1) + devs(2)*devs(2) + devs(3)*devs(3))
     .          + devs(4)*devs(4) + devs(5)*devs(5) + devs(6)*devs(6)
c .....................................................................
c 
c ...
      cam_clay_pc = 0.d0
      if( dabs(sm) .lt. zero ) return
      cam_clay_pc = - (q2/(mcs*mcs*sm) + sm)
c .....................................................................
      return
      end

