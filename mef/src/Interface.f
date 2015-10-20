       subroutine istress2d(d11,d22,l,u,du,tx,sy,flagt)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   ISTRESS2D: calcula tensoes 2D no elemento de interface           *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     du(8)  - deslocamento                                          *
c *     d11,d22   - coeficientes da matriz constitutiva                *
c *     sy      - tensao de escoamento                                 *
c *     dtx(2) - delta tensao linear                                   *
c *     tx(2)   - tensoes anteriores                                   *
c *     flagt  - flag que indica tracao ou compressao                  *
c *     l     - comprimento do elemento                                *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     tx(2)   - tensoes atualizadas: tx = tx + dtx                   *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical flagt
      integer i,j
      real*8  u(*),du(*),tx(*),tx0(2),dtx(2),se,d11,d22,l,epe,e2
      real*8  preys,sy,hl,dste(2),se0,escur,rfac,stot(2),reduc
c ......................................................................
c    
      flagt = .false.
      se0 = tx(3)
c ... calculate elastic stress increment
      dtx(1) = d22*((du(8)-du(2)+du(6)-du(4)))*0.5d0*l !dsigma normal
      dtx(2) = d11*((du(7)-du(1)+du(5)-du(3)))*0.5d0*l !dsigma tangencial
      preys=sy+epe*hl
      dste(1) = dtx(1)
      dste(2) = dtx(2)
c ...    t0 = t:
      call aequalb(tx0,tx,2)
c ...    Tensoes lineares, t = t0 + dt
      call vsum(tx0,dste,2,tx)
c ... evaluate effective stress
      se = tx(1)
      if ( (se0-preys) .ge. 0.) goto 3
c ... check for yielding durint this iterations
      escur = se - preys
      if(escur .le. 0.) goto 4
c ... compute factor R
      rfac = escur/(se-se0)
      goto 5
    3 escur = se - se0
      if(escur .le. 0.d0)   then
         tx(1) = d22*((u(8)-u(2)+u(6)-u(4)))*0.5d0*l !tensao normal
         tx(2) = d11*((u(7)-u(1)+u(5)-u(3)))*0.5d0*l !tensao tangencial
         if (tx(1) .gt. 0.d0)   then
            flagt = .true.
            tx(1) = 0.d0
            tx(2) = 0.d0
            tx(3) = 0.d0
         endif
         goto 12
      endif
      rfac = 1.d0
    5 continue
      flagt = .true. 
      do 6 j = 1,2
    6    stot(j) = 0.d0
      do 11 j = 1,2
   11    tx(j) = stot(j)
      tx(3) = 0.d0
      goto 12
c ... compute stresses for elastic nodes or points
    4 continue
      tx(3) = se
   12 return
      end
      subroutine istress2d_grampo(d11,d22,l,u,du,tx,E2,pr,sy,thic,
     .flagt,nlit)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   ISTRESS2D: calcula tensoes 2D no elemento de interface           *
c *   -----------     com resistencia a tracao                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     du(8)  - deslocamento                                          *
c *     d11,d22   - coeficientes da matriz constitutiva                *
c *     sy      - tensao de escoamento                                 *
c *     dtx(2) - delta tensao linear                                   *
c *     tx(2)   - tensoes anteriores                                   *
c *     flagt  - flag que indica tracao ou compressao                  *
c *     l     - comprimento do elemento                                *
c *     e2    - modulo de elasticidade em tracao                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     tx(2)   - tensoes atualizadas: tx = tx + dtx                   *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical flagt
      integer i,j,nlit
      real*8  u(*),du(*),tx(*),tx0(2),dtx(2),se,d11,d22,l,epe,e2,thic
      real*8  preys,sy,hl,dste(2),se0,escur,rfac,stot(2),pr,delta
c ......................................................................
c    
      flagt = .false.
      se0 = tx(3)
      d22 = e2*(1.d0-pr)/((1.d0+pr)*(1.d0-2.d0*pr))
      d11 = e2/(2.d0*(1.d0+pr))
c ... calculate elastic stress increment
      dtx(1) = d22*((du(8)-du(2)+du(6)-du(4)))*0.5d0 !dsigma normal
      dtx(2) = d11*((du(7)-du(1)+du(5)-du(3)))*0.5d0 !dsigma tangencial
c      preys=sy+epe*hl
      preys = 0.d0
      dste(1) = dtx(1)
      dste(2) = dtx(2)
c ...    t0 = t:
      call aequalb(tx0,tx,2)
c ...    Tensoes lineares, t = t0 + dt
      call vsum(tx0,dste,2,tx)
c ... evaluate effective stress
      se = tx(1)
      if ( (se0-preys) .ge. 0.d0) goto 3
c ... check for yielding durint this iterations
      escur = se - preys
      if(escur .le. 0.d0) goto 4
c ... compute factor R
      rfac = escur/(se-se0)
      goto 5
    3 escur = se - se0
      if(escur .le. 0.d0)   then
         tx(1) = d22*((u(8)-u(2)+u(6)-u(4)))*0.5d0 !tensao normal
         tx(2) = d11*((u(7)-u(1)+u(5)-u(3)))*0.5d0 !tensao tangencial
         if (tx(1) .gt. 0.d0)   then
            flagt = .true.
c            d22 = e2*(1.d0-pr)/((1.d0+pr)*(1.d0-2.d0*pr))
c            d11 = e2/(2.d0*(1.d0+pr))
c 100        delta = ((u(8)-u(2)+u(6)-u(4)))*0.5d0*l
c            dtx(1) = d22*((du(8)-du(2)+du(6)-du(4)))*0.5d0*l !tensao normal
c            dtx(2) = d11*((du(7)-du(1)+du(5)-du(3)))*0.5d0*l !tensao tangencial
c            tx(1) = tx0(1)+dtx(1)
c            tx(2) = tx0(2)+dtx(2)
            tx(3) = tx(1)
c            if (tx(1) .gt. sy)  then
c               tx(1) = sy
c               tx(2) = 0.d0
c               tx(3) = tx(1)
c            endif
         endif
         goto 12
      endif
      rfac = 1.d0
    5 continue
      flagt = .true. 
c      d22 = e2*(1.d0-pr)/((1.d0+pr)*(1.d0-2.d0*pr))
c      d11 = e2/(2.d0*(1.d0+pr))
      tx(1) = d22*((u(8)-u(2)+u(6)-u(4)))*0.5d0 !tensao normal
      tx(2) = d11*((u(7)-u(1)+u(5)-u(3)))*0.5d0 !tensao tangencial
      tx(3) = tx(1)
      goto 12
c ... compute stresses for elastic nodes or points
    4 continue
      tx(3) = se
   12 return
      end
       subroutine istress3d(d11,d22,d33,l,u,du,tx,sy,flagt)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   ISTRESS2D: calcula tensoes 2D no elemento de interface           *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     du(8)  - deslocamento                                          *
c *     d11,d22   - coeficientes da matriz constitutiva                *
c *     sy      - tensao de escoamento                                 *
c *     dtx(2) - delta tensao linear                                   *
c *     tx(2)   - tensoes anteriores                                   *
c *     flagt  - flag que indica tracao ou compressao                  *
c *     l     - comprimento do elemento                                *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     tx(2)   - tensoes atualizadas: tx = tx + dtx                   *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical flagt
      integer i,j,nst,m
      real*8  u(*),du(*),tx(*),tx0(2),dtx(3),se,d11,d22,l,epe,e2
      real*8  preys,sy,hl,dste(3),se0,escur,rfac,stot(3),reduc,d33
      real*8 deps(3),eps(3)
c ......................................................................
c    
      flagt = .false.
      se0 = tx(4)
      deps(1) = (du(12)+du(15)+du(18)-du(3)-du(6)-du(9))*l/3.d0
      deps(2) = (du(10)+du(13)+du(16)-du(1)-du(4)-du(7))*l/3.d0
      deps(3) = (du(11)+du(14)+du(17)-du(2)-du(5)-du(8))*l/3.d0
      eps(1) = (u(12)+u(15)+u(18)-u(3)-u(6)-u(9))*l/3.d0
      eps(2) = (u(10)+u(13)+u(16)-u(1)-u(4)-u(7))*l/3.d0
      eps(3) = (u(11)+u(14)+u(17)-u(2)-u(5)-u(8))*l/3.d0
c ... calculate elastic stress increment
      dtx(1) = d11*deps(1) ! sigmaz
      dtx(2) = d22*deps(2) ! talzx
      dtx(3) = d33*deps(3) ! talzy
      preys=sy+epe*hl
      dste(1) = dtx(1)
      dste(2) = dtx(2)
      dste(3) = dtx(3)
c ...    t0 = t:
      call aequalb(tx0,tx,3)
c ...    Tensoes lineares, t = t0 + dt
      call vsum(tx0,dste,3,tx)
c ... evaluate effective stress
      se = tx(1)
      if ((se0-preys) .ge. 0.d0) goto 3
c ... check for yielding durint this iterations
      escur = se - preys
      if(escur .le. 0.) goto 4
c ... compute factor R
      rfac = escur/(se-se0)
      goto 5
    3 escur = se - se0
      if(escur .le. 0.d0)   then
         tx(1) = d11*eps(1) ! sigmaz
         tx(2) = d22*eps(2) ! talzx
         tx(3) = d33*eps(3) ! talzy
         if (tx(1) .gt. 0.d0)   then
            flagt = .true.
            tx(1) = 0.d0
            tx(2) = 0.d0
            tx(3) = 0.d0
            tx(4) = 0.d0
         endif
         goto 12
      endif
      rfac = 1.d0
    5 continue
      flagt = .true. 
      do 6 j = 1,3
    6    stot(j) = 0.d0
      do 11 j = 1,3
   11    tx(j) = stot(j)
      tx(4) = 0.d0
      goto 12
c ... compute stresses for elastic nodes or points
    4 continue
      tx(4) = se
   12 return
      end
       subroutine istress3d_hexa(d33,d22,d11,rn,rt,l,u,du,tx,sy,
     .                                                     flagt)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   ISTRESS2D: calcula tensoes 2D no elemento de interface           *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     du(8)  - deslocamento                                          *
c *     d11,d22   - coeficientes da matriz constitutiva                *
c *     sy      - tensao de escoamento                                 *
c *     dtx(2) - delta tensao linear                                   *
c *     tx(2)   - tensoes anteriores                                   *
c *     flagt  - flag que indica tracao ou compressao                  *
c *     l     - comprimento do elemento                                *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     tx(2)   - tensoes atualizadas: tx = tx + dtx                   *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical flagt
      integer i,j,nst,m
      real*8  u(*),du(*),tx(*),tx0(2),dtx(3),se,d11,d22,l,epe,e2
      real*8  preys,sy,hl,dste(3),se0,escur,rfac,stot(3),reduc,d33
      real*8 deps(3),eps(3),rn,rt
c ......................................................................
c    
      flagt = .false.
      se0 = tx(4)
      deps(1) = -(du(15)+du(18)+du(21)+du(24)
     .            -du(3)-du(6)-du(9)-du(12))*l/4.d0
      deps(2) = -(du(13)+du(16)+du(19)+du(22)
     .          -du(1)-du(4)-du(7)-du(10))*l/4.d0
      deps(3) = -(du(14)+du(17)+du(20)+du(23)
     .-du(2)-du(5)-du(8)-du(11))*l/4.d0
      eps(1) = -(u(15)+u(18)+u(21)+u(24)-u(3)-u(6)-u(9)-u(12))*l/4.d0
      eps(2) = -(u(13)+u(16)+u(19)+u(22)-u(1)-u(4)-u(7)-u(10))*l/4.d0
      eps(3) = -(u(14)+u(17)+u(20)+u(23)-u(2)-u(5)-u(8)-u(11))*l/4.d0
c ... calculate elastic stress increment
      dtx(1) = d11*deps(1) ! sigmaz
      dtx(2) = d22*deps(2) ! talzx
      dtx(3) = d33*deps(3) ! talzy
      preys=sy+epe*hl
      dste(1) = dtx(1)
      dste(2) = dtx(2)
      dste(3) = dtx(3)
c ...    t0 = t:
      call aequalb(tx0,tx,3)
c ...    Tensoes lineares, t = t0 + dt
      call vsum(tx0,dste,3,tx)
c ... evaluate effective stress
      tx(1) = d11*eps(1) ! sigmaz
      tx(2) = d22*eps(2) ! talzx
      tx(3) = d33*eps(3) ! talzy
      se = tx(1)
      if ( (se0-preys) .ge. 0.d0) goto 3
c ... check for yielding durint this iterations
      escur = se - preys
      tx(1) = d11*eps(1) ! sigmaz
      tx(2) = d22*eps(2) ! talzx
      tx(3) = d33*eps(3) ! 
      if(escur .le. 0.d0) goto 4
c ... compute factor R
      rfac = escur/(se-se0)
      goto 5
    3 escur = se - se0
      if(escur .le. 0.d0)   then
         tx(1) = d11*eps(1) ! sigmaz
         tx(2) = d22*eps(2) ! talzx
         tx(3) = d33*eps(3) ! talzy
         if (tx(1) .gt. 0.d0)   then
            flagt = .true.
            tx(1) = rn*eps(1)/l ! sigmaz
            tx(2) = rt*eps(2)/l ! talzx
            tx(3) = rt*eps(3)/l ! talzy
            tx(4) = tx(1)
         endif
         goto 12
      endif
      rfac = 1.d0
    5 continue
      flagt = .true. 
      tx(1) = rn*eps(1)/l ! sigmaz
      tx(2) = rt*eps(2)/l ! talzx
      tx(3) = rt*eps(3)/l ! 
      tx(4) = tx(1)
      goto 12
c ... compute stresses for elastic nodes or points
    4 continue
      tx(4) = se
   12 return
      end