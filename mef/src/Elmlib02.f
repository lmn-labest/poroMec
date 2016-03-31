      subroutine elmt02_m(e,x,u,v,a,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,
     .                                       ma,nlit,isw,flaghidr)
c **********************************************************************
c *                                                                    *
c *   ELMT02: Elemento triangular linear (Estado plano de deformacao)  *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     txe(nst)    - tensoes por elemento                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'                  
      integer ndm,nst,nel,isw,ma,iyied
      integer i,j,k,l,nlit,ipl
      real*8 e(*),x(ndm,*),u(nst),v(nst),a(nst),p(*),s(nst,nst),tl(*)
      real*8 hx(3),hy(3),det,wt,xj11,xj12,xj21,xj22,ut(*)
      real*8 xji11,xji12,xji21,xji22,tm,h(3)
      real*8 ym,pr,a1,b,c,c1,c2,d11,d12,d22,d33,sm1,sm2,alpha,hi(3)
c      real*8 tx,ty,txy
      real*8 tx(*),eps(*),deps(4),yield,epsi(4),tx1,tx2,tx4,txl(5)
      real*8 hl,phi,itr,c14,g,coes,txp(*)
      logical flaghidr
c ......................................................................
      goto (100,200,200,400) isw
c ======================================================================
c
c.... Input material properties:
c
c ......................................................................
  100 continue
c     e(1) = modulo de elasticidade
c     e(2) = coeficiente de Poisson
c     e(3) = massa específica
c     e(4) = coeficiente de dilatacao termica
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
c
c     Matriz Jacobiana:
c
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) go to 1000
c
c Inversa da matriz Jacobiana:
c
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao:
c
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
c
c     Matriz constitutiva:
c
      ym    = e(1)
      pr    = e(2)
      alpha = e(4)
      iyied =int(e(5))
      yield = e(11)
      phi   = e(12)
      coes  = e(13)
      hl    = e(14)
      call radian(phi)
      if( iyied .eq. 3)   yield = coes*dcos(phi)
      if( iyied .eq. 4)   yield = 6.d0*coes*dcos(phi)/(dsqrt(3.d0)*
     .    (3.d0-dsin(phi)))
      h(1) = 1.d0/3.d0
      h(2) = h(1)
      h(3) = h(1)
      call check_var_hidr(a1,b,c,d11,d12,d22,d33,ym,pr,alpha,hi(2),
     .                                ma,ndm,flaghidr,.true.)
      call check_var_term(a1,b,c,d11,d12,d22,d33,h,ut,ym,pr,alpha,
     .                        3,ndm,ma,flaghidr,.true.)
c      tm =  (ut(1)+ut(2)+ut(3))/3.d0
c      if (flaghidr .eqv. .true.) then
c         call updhidr(ym,hi(2),ma,1)
c         call updhidr(pr,hi(2),ma,2)
c         call updhidr(alpha,hi(2),ma,4)
c      else
c         call updprop(ym,tm,ma,1)
c         call updprop(pr,tm,ma,2)
c         call updprop(alpha,tm,ma,4)
c      endif
c .........................................
c      c1 = 1.d0-pr
c      c2 = 1.d0-2.d0*pr
c      a1 = ym*c1/((1.+pr)*c2)
c      b = pr/c1
c      c = c2/(2.d0*c1)
c      d11 = a1
c      d12 = a1*b
c      d22 = a1
c      d33 = a1*c
      if (isw .eq. 3) goto 300
c........................................
      wt = 0.5d0*det
c
c     Matriz de rigidez:                     
c
      do 220 j = 1, 3
         k = (j-1)*2+1
         do 210 i = 1, 3
           l = (i-1)*2+1
c ------------------------------------------------------------------
           s(l,k)     = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
c ------------------------------------------------------------------
           s(l,k+1)   = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
c ------------------------------------------------------------------
           s(l+1,k)   = ( hy(i)*d12*hx(j) + hx(i)*d33*hy(j) ) * wt
c ------------------------------------------------------------------
           s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt
c ------------------------------------------------------------------
  210    continue
  220 continue
c ......................................................................  
c
c ... Matriz de massa:
c
c      sm1 = (1.d0/12.d0)*det*e(3)
c      sm2 = sm1*0.5d0 
c
c ... Matriz de rigidez efetiva:                     
c
c      c = beta*dt*dt
c..................
c      c = 1.d0
c      sm1 = 0.d0
c      sm2 = 0.d0
c..................      
c      s(1,1) = s(1,1)*c + sm1
c      s(1,2) = s(1,2)*c
c      s(1,3) = s(1,3)*c + sm2
c      s(1,4) = s(1,4)*c
c      s(1,5) = s(1,5)*c + sm2
c      s(1,6) = s(1,6)*c
c      s(2,2) = s(2,2)*c + sm1
c      s(2,3) = s(2,3)*c
c      s(2,4) = s(2,4)*c + sm2
c      s(2,5) = s(2,5)*c
c      s(2,6) = s(2,6)*c + sm2    
c      s(3,3) = s(3,3)*c + sm1
c      s(3,4) = s(3,4)*c
c      s(3,5) = s(3,5)*c + sm2
c      s(3,6) = s(3,6)*c
c      s(4,4) = s(4,4)*c + sm1
c      s(4,5) = s(4,5)*c
c      s(4,6) = s(4,6)*c + sm2      
c      s(5,5) = s(5,5)*c + sm1
c      s(5,6) = s(5,6)*c
c      s(6,6) = s(6,6)*c + sm1
c      s(2,1) = s(1,2)   
c      s(3,1) = s(1,3)
c      s(3,2) = s(2,3)
c      s(4,1) = s(1,4)    
c      s(4,2) = s(2,4)
c      s(4,3) = s(3,4)
c      s(5,1) = s(1,5)
c      s(5,2) = s(2,5) 
c      s(5,3) = s(3,5)
c      s(5,4) = s(4,5)
c      s(6,1) = s(1,6)
c      s(6,2) = s(2,6)
c      s(6,3) = s(3,6)
c      s(6,4) = s(4,6)
c      s(6,5) = s(5,6)            
c ... deformacoes:      
      call deform2d0(hx,hy,u,0.d0,epsi,3)
c      call deform2d(hx,hy,u,eps,3)
c ... deformacoes termicas:
      tm = (tl(1)+tl(2)+tl(3))/3.d0
      epsi(1) = epsi(1) - alpha*tm
      epsi(2) = epsi(2) - alpha*tm
      epsi(3) = epsi(3) - alpha*tm
         tm = (tl(1)+tl(2)+tl(3))/3.d0
      if ( nlit .eq. 1)  then
         tx1  = d11*epsi(1) + d12*epsi(2) - d12*alpha*tm -txp(1)
         tx2  = d12*epsi(1) + d22*epsi(2) - d12*alpha*tm -txp(3)
         tx4  = d33*epsi(4) - txp(2)
c         tx1  = d11*epsi(1) + d12*epsi(2) - d12*alpha*tm 
c         tx2  = d12*epsi(1) + d22*epsi(2) - d12*alpha*tm 
c         tx4  = d33*epsi(4)
         p(1) = (hx(1)*tx1+hy(1)*tx4) * wt
         p(2) = (hy(1)*tx2+hx(1)*tx4) * wt
         p(3) = (hx(2)*tx1+hy(2)*tx4) * wt
         p(4) = (hy(2)*tx2+hx(2)*tx4) * wt
         p(5) = (hx(3)*tx1+hy(3)*tx4) * wt
         p(6) = (hy(3)*tx2+hx(3)*tx4) * wt
         return
      endif
      deps(1) = epsi(1) - eps(1)-alpha*tm
      deps(2) = epsi(2) - eps(2)-alpha*tm
      deps(3) = epsi(3) - eps(3)-alpha*tm
      deps(4) = epsi(4) - eps(4)
c ... deformacoes plasticas:
      eps(1) = eps(1) + deps(1)
      eps(2) = eps(2) + deps(2)
      eps(3) = eps(3) + deps(3)
      eps(4) = eps(4) + deps(4)
      txl(1) = tx(1)
      txl(2) = tx(4)
      txl(3) = tx(2)
      txl(4) = tx(3)
      txl(5) = tx(5)
      c14 = pr/(1.d0-2.d0*pr)
      g = ym/(2.d0*(1.d0+pr))
      call plasticst2d(a1,b,c,deps,eps,txl,txp,yield,hl,phi,pr,
     .    iyied,c14,g,.true.)
      tx(1) = txl(1)
      tx(2) = txl(3)
      tx(3) = txl(4)
      tx(4) = txl(2)
      tx(5) = txl(5)
c      write (*,*) tx(1),tx(2),tx(3),tx(4)
c      call deform2d0(hx,hy,v,0.d0,deps,3)
c ... tensoes:
c      tx(1)  = d11*eps(1) + d12*eps(2) - d12*alpha*tm
c      tx(2)  = d12*eps(1) + d22*eps(2) - d12*alpha*tm
c      tx(4)  = d33*eps(3)
c ......................................................................
c
c ... Forcas internas:
c
c      p(1) = (hx(1)*tx + hy(1)*txy)*wt + sm1*a(1)+sm2*a(3)+sm2*a(5)
c      p(2) = (hy(1)*ty + hx(1)*txy)*wt + sm1*a(2)+sm2*a(4)+sm2*a(6)
c      p(3) = (hx(2)*tx + hy(2)*txy)*wt + sm2*a(1)+sm1*a(3)+sm2*a(5)
c      p(4) = (hy(2)*ty + hx(2)*txy)*wt + sm2*a(2)+sm1*a(4)+sm2*a(6)
c      p(5) = (hx(3)*tx + hy(3)*txy)*wt + sm2*a(1)+sm2*a(3)+sm1*a(5)
c      p(6) = (hy(3)*ty + hx(3)*txy)*wt + sm2*a(2)+sm2*a(4)+sm1*a(6)
c      p(1) = (hx(1)*tx + hy(1)*txy)*wt
c      p(2) = (hy(1)*ty + hx(1)*txy)*wt
c      p(3) = (hx(2)*tx + hy(2)*txy)*wt
c      p(4) = (hy(2)*ty + hx(2)*txy)*wt
c      p(5) = (hx(3)*tx + hy(3)*txy)*wt
c      p(6) = (hy(3)*ty + hx(3)*txy)*wt
      p(1) = (hx(1)*tx(1)+hy(1)*tx(4)) * wt
      p(2) = (hy(1)*tx(2)+hx(1)*tx(4)) * wt
      p(3) = (hx(2)*tx(1)+hy(2)*tx(4)) * wt
      p(4) = (hy(2)*tx(2)+hx(2)*tx(4)) * wt
      p(5) = (hx(3)*tx(1)+hy(3)*tx(4)) * wt
      p(6) = (hy(3)*tx(2)+hx(3)*tx(4)) * wt
c      call lku(s,u,p,nst)
      return
c ======================================================================
c
c.... Tensoes nodais:
c
c ......................................................................
  300 continue
        do 410 i = 1, 3
         k = (i-1)*5+1
	   p(k)   = tx(1)
	   p(k+1) = tx(2)
	   p(k+2) = tx(3)
	   p(k+3) = tx(4)
	   p(k+4) = tx(5)
  410 continue
c      call deform2d(hx,hy,u,eps,3)
c ... deformacoes termicas:      
c      tm = (tl(1)+tl(2)+tl(3))/3.d0
c      eps(1) = eps(1) - alpha*tm
c      eps(2) = eps(2) - alpha*tm
c      p(1) = d11*eps(1) + d12*eps(2) - d12*alpha*tm
c      p(2) = d12*eps(1) + d22*eps(2) - d12*alpha*tm
c      p(3) = d33*eps(3)
c      druckerprager(p,phi,coes,mu,lamb)
c      call stress2d(d11,d12,d22,d33,eps,p)
c      p(4) = p(1)
c      p(5) = p(2)
c      p(6) = p(3)
c      p(7) = p(1)
c      p(8) = p(2)
c      p(9) = p(3)
      return
c ......................................................................            
  400 continue
      return
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT02: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
      subroutine elmt03_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,
     .                                         ma,nlit,isw,flaghidr)
c **********************************************************************
c *                                                                    *
c *   ELMT03: Elemento triangular linear (Estado plano de tensao)      *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     v(nst)     - solucao anterior incremental                      *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     tx(ntn)    - tensoes no elemento                               *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      integer ndm,nst,nel,isw,i,j,k,l,ma,nlit,iyied
      real*8 e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*),ut(*)
      real*8 hx(3),hy(3),hz(3),det,wt,xj11,xj12,xj21,xj22
      real*8 xji11,xji12,xji21,xji22,tm,h(3)
      real*8 ym,pr,a,thic,d11,d12,d22,d33,alpha,hi(3)
      real*8 ty,txy,b,c,txl(5),hl,phi,g,c14,coes
      real*8 tx(*),eps(*),deps(4),yield,epsi(4),tx1,tx2,tx4,txp(*)
      logical flaghidr
c ......................................................................
      go to (100,200,300,400) isw
c ======================================================================
c
c.... Input material properties:
c
c ......................................................................
  100 continue
c     e(1) = modulo de elasticidade
c     e(2) = coeficiente de Poisson
c     e(3) = espessura
c     e(4) = coeficiente de dilatacao termica
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
c
c     Matriz Jacobiana:
c
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) go to 1000
c
c Inversa da matriz Jacobiana:
c
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao:
c
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
      hz(1) =  0.d0
      hz(2) =  0.d0
      hz(3) =  0.d0
c
c     Matriz constitutiva:
c
      ym    = e(1)
      pr    = e(2)
      thic  = e(3)
      alpha = e(4)
      iyied =int(e(5))
      yield = e(11)
      phi   = e(12)
      coes  = e(13)
      hl    = e(14)
      call radian(phi)
      if( iyied .eq. 3)   yield = coes*dcos(phi)
      if( iyied .eq. 4)   yield = 6.d0*coes*dcos(phi)/(dsqrt(3.d0)*
     .    (3.d0-dsin(phi)))
c
c      a = ym/(1.d0-pr*pr)
c      d11 = a
c      d12 = a*pr
c      d22 = a
c      d33 = ym/(2.d0*(1.d0+pr))
      h(1) = 1.d0/3.d0
      h(2) = h(1)
      h(3) = h(1)
      call check_var_hidr(a,b,c,d11,d12,d22,d33,ym,pr,alpha,hi(2),
     .                                ma,ndm,flaghidr,.false.)
      call check_var_term(a,b,c,d11,d12,d22,d33,h,ut,ym,pr,alpha,
     .                        3,ndm,ma,flaghidr,.false.)
c      tm =  (ut(1)+ut(2)+ut(3))/3.d0
c      if (flaghidr .eqv. .true.) then
c         call updhidr(ym,hi(2),ma,1)
c         call updhidr(pr,hi(2),ma,2)
c         call updhidr(alpha,hi(2),ma,4)
c      else
c         call updprop(ym,tm,ma,1)
c         call updprop(pr,tm,ma,2)
c         call updprop(alpha,tm,ma,4)
c      endif
c
c      a = ym/(1.d0-pr*pr)
c      d11 = a
c      d12 = a*pr
c      d22 = a
c      d33 = ym/(2.d0*(1.d0+pr))
      if (isw .eq. 3) goto 300
c ......................................................................
      wt = 0.5d0*det*thic
c
c     Matriz de rigidez:                     
c
      do 220 j = 1, 3
         k = (j-1)*2+1
         do 210 i = 1, 3
           l = (i-1)*2+1
c ------------------------------------------------------------------
           s(l,k)     = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
c ------------------------------------------------------------------
           s(l,k+1)   = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
c ------------------------------------------------------------------
           s(l+1,k)   = ( hy(i)*d12*hx(j) + hx(i)*d33*hy(j) ) * wt
c ------------------------------------------------------------------
           s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt
c ------------------------------------------------------------------
  210    continue
  220 continue
c ======================================================================
c .... Forcas Internas
c ....
c ... deformacoes:      
      call deform2d0(hx,hy,u,b,epsi,3)
      tm = (tl(1)+tl(2)+tl(3))/3.d0
      if ( nlit .eq. 1)  then
c         tm = (tl(1)+tl(2)+tl(3))/3.d0
         epsi(1) = epsi(1) - alpha*tm
         epsi(2) = epsi(2) - alpha*tm
         tx1  = d11*epsi(1) + d12*epsi(2) - txp(1)
         tx2  = d12*epsi(1) + d22*epsi(2) -txp(3)
         tx4  = d33*epsi(4) - txp(2)
         p(1) = (hx(1)*tx1+hy(1)*tx4) * wt
         p(2) = (hy(1)*tx2+hx(1)*tx4) * wt
         p(3) = (hx(2)*tx1+hy(2)*tx4) * wt
         p(4) = (hy(2)*tx2+hx(2)*tx4) * wt
         p(5) = (hx(3)*tx1+hy(3)*tx4) * wt
         p(6) = (hy(3)*tx2+hx(3)*tx4) * wt
         return
      endif
c      call deform2d(hx,hy,u,eps,3)
c ... deformacoes termicas:
c      tm = (tl(1)+tl(2)+tl(3))/3.d0
c      epsi(1) = epsi(1) - alpha*tm
c      epsi(2) = epsi(2) - alpha*tm
c      epsi(3) = epsi(3) - alpha*tm
      a = ym/(1.d0-pr*pr)
      b = pr
      c = (1.d0-pr*pr)/(2.d0+2.d0*pr)
      deps(1) = epsi(1) - eps(1)- alpha*tm
      deps(2) = epsi(2) - eps(2)- alpha*tm
      deps(3) = epsi(3) - eps(3)
      deps(4) = epsi(4) - eps(4)
c ... deformacoes plasticas:
      eps(1) = eps(1) + deps(1)
      eps(2) = eps(2) + deps(2)
      eps(3) = eps(3) + deps(3)
      eps(4) = eps(4) + deps(4)
      txl(1) = tx(1)
      txl(2) = tx(4)
      txl(3) = tx(2)
      txl(4) = tx(3)
      txl(5) = tx(5)
      c14 = pr/(1.d0-pr)
      g = ym/(2.d0*(1.d0+pr))
c      a = ym*(1.d0-pr)/((1.+pr)*(1.d0-2.d0*pr))
c      b = pr/(1.d0-pr)
c      c = (1.d0-2.d0*pr)/(2.d0*(1.d0-pr))
      call plasticst2d(a,b,c,deps,eps,txl,txp,yield,hl,phi,0.d0,iyied,
     .    c14,g,.false.)
      tx(1) = txl(1)
      tx(2) = txl(3)
      tx(3) = txl(4)
      tx(4) = txl(2)
      tx(5) = txl(5)
      p(1) = (hx(1)*tx(1)+hy(1)*tx(4)) * wt
      p(2) = (hy(1)*tx(2)+hx(1)*tx(4)) * wt
      p(3) = (hx(2)*tx(1)+hy(2)*tx(4)) * wt
      p(4) = (hy(2)*tx(2)+hx(2)*tx(4)) * wt
      p(5) = (hx(3)*tx(1)+hy(3)*tx(4)) * wt
      p(6) = (hy(3)*tx(2)+hx(3)*tx(4)) * wt
c      call lku(s,u,p,nst)
c ... deformacoes:      
c      call deform2d(hx,hy,u,eps,3)
c ... deformacoes termicas:
c      tm = (tl(1)+tl(2)+tl(3))/3.d0
c      eps(1) = eps(1) - alpha*tm
c      eps(2) = eps(2) - alpha*tm
c ... tensoes:
c      tx  = d11*eps(1) + d12*eps(2)
c      ty  = d12*eps(1) + d22*eps(2)
c      txy = d33*eps(3)
c ... Forcas internas:
c
c      p(1) = (hx(1)*tx + hy(1)*txy)*wt
c      p(2) = (hy(1)*ty + hx(1)*txy)*wt
c      p(3) = (hx(2)*tx + hy(2)*txy)*wt
c      p(4) = (hy(2)*ty + hx(2)*txy)*wt
c      p(5) = (hx(3)*tx + hy(3)*txy)*wt
c      p(6) = (hy(3)*ty + hx(3)*txy)*wt
c ......................................................................
c      call lku(s,u,p,nst)
      return
c ======================================================================
c
c.... Tensoes nodais:       
c
c ......................................................................
  300 continue
      do 310 i = 1, 3
         k = (i-1)*5+1
         p(k)   = tx(1)
         p(k+1) = tx(2)
         p(k+2) = tx(3)
         p(k+3) = tx(4)
         p(k+4) = tx(5)
  310 continue
c      call deform2d(hx,hy,u,eps,3)
cc ... deformacoes termicas:      
c      tm = (tl(1)+tl(2)+tl(3))/3.d0
c      eps(1) = eps(1) - alpha*tm
c      eps(2) = eps(2) - alpha*tm
c      call stress2d(d11,d12,d22,d33,eps,p)
c      p(4) = p(1)
c      p(5) = p(2)
c      p(6) = p(3)
c      p(7) = p(1)
c      p(8) = p(2)
c      p(9) = p(3)
      return
  400 continue
  
      return
c ======================================================================
 1000 continue
      print*, '*** Subrotina ELMT03: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
      subroutine elmt04_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,
     .                                           ma,nlit,isw,flaghidr)
c **********************************************************************
c *                                                                    *
c *   ELMT04: Elemento quadrilatero de 4 nos                           *
c *   ------  estado plano de deformacao                               *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     v(nst)     - solucao anterior incremental                      *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     tx(ntn)    - tensoes no elemento                               *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      common /gauss/ pg,wg
      integer ndm,nst,nel,isw,nint,nlit,pint(4),nen
      integer i,j,k,l,lx,ly,ma,iyied
      real*8  e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*),ut(*)
      real*8  h(4),hx(4),hy(4),xj(2,2),xji(2,2)
      real*8  ym,pr,a,b,c,b1,b2,yield,d11,d12,d22,d33,ri,si,wt,det,tm
      real*8  pg(10,10),wg(10,10)
      real*8  eps(5,*),alpha,hi(3),tx(5,*),tx1,tx2,tx4,epsi(4),deps(4)
      real*8 invar2,rn(4),sn(4),c1,c2,c3,c4,txl(5)
      real*8  phi, coes, hl,c14,g,txp(5,*)
c     real*8      tx,ty,txy
      logical flaghidr
      data    rn /1.d0,-1.d0,-1.d0,1.d0/,sn /1.d0,1.d0,-1.d0,-1.d0/
      data    nint/2/,pint/4,2,1,3/
c ......................................................................
      goto (100,200,300,400) isw
c ======================================================================
c
c.... Imput material properties:
c
c ......................................................................
  100 continue
c     e(1) = modulo de elasticidade
c     e(2) = coeficiente de Poisson
      return
c ======================================================================
c
c.... Matrix K:
c ......................................................................
  200 continue
c
c ... Matriz constitutiva:
c ......................................................................
      ym    = e(1)
      pr    = e(2)
      alpha = e(4)
      iyied =int(e(5))
      yield = e(11)
      phi   = e(12)
      coes  = e(13)
      hl    = e(14)
      tm = 0.d0
      call radian(phi)
      if( iyied .eq. 3)   yield = coes*dcos(phi)
      if( iyied .eq. 4)   yield = 6.d0*coes*dcos(phi)/(dsqrt(3.d0)*
     .    (3.d0-dsin(phi)))
      call check_var_hidr(a,b,c,d11,d12,d22,d33,ym,pr,alpha,hi(2),
     .                           ma,ndm,flaghidr,.true.)
c
c       ym = ym*t
c
c
c ......................................................................
c      b1 = 1.d0-pr
c      b2 = 1.d0-2.d0*pr
c      a = ym*b1/((1.+pr)*b2)
c      b = pr/b1
c      c = b2/(2.d0*b1)
c ......................................................................
c      d11 = a
c      d12 = a*b
c      d22 = a
c      d33 = a*c
c      if (isw .eq. 3) goto 300
c ......................................................................
c
c ... Matriz de rigidez:
c ......................................................................
      do i = 1, nst
      do j = 1, nst
         s(i,j) = 0.d0
      enddo
      enddo
      do 240 lx = 1, nint
         ri = pg(lx,nint)
         do 230 ly = 1, nint
            si = pg(ly,nint)
            call sfquad4(h,hx,hy,ri,si,.true.,.true.)
c ...........
            call check_var_term(a,b,c,d11,d12,d22,d33,h,ut,ym,pr,alpha,
     .                        4,ndm,ma,flaghidr,.true.)
c            if (flaghidr .neqv. .true.) then
c               tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c               call updprop(ym,tm,ma,1)
c               call updprop(pr,tm,ma,2)
c               b1 = 1.d0-pr
c               b2 = 1.d0-2.d0*pr
c               a = ym*b1/((1.+pr)*b2)
c               b = pr/b1
c               c = b2/(2.d0*b1)
c               d11 = a
c               d12 = a*b
c               d22 = a
c               d33 = a*c
c            endif
c ...........
            call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)
            wt = wg(lx,nint)*wg(ly,nint)*det
        do 220 j = 1, 4
               k = (j-1)*2+1
               do 210 i = 1, 4
                  l = (i-1)*2+1
c --------------------------------------------------------------------
                  s(l,k)     = s(l,k) +
     .                        (hx(i)*d11*hx(j)+hy(i)*d33*hy(j))*wt
c --------------------------------------------------------------------
                  s(l,k+1)   = s(l,k+1) +
     .                        (hx(i)*d12*hy(j)+hy(i)*d33*hx(j))*wt
c --------------------------------------------------------------------
                  s(l+1,k)   = s(l+1,k) +
     .                        (hy(i)*d12*hx(j)+hx(i)*d33*hy(j))*wt
c --------------------------------------------------------------------
                  s(l+1,k+1) = s(l+1,k+1) +
     .                        (hy(i)*d22*hy(j)+hx(i)*d33*hx(j))*wt
c --------------------------------------------------------------------
  210          continue
  220       continue
  230    continue
  240 continue
c ======================================================================
c
c ... Forcas internas:
c
	p(1) = 0.d0
	p(2) = 0.d0 
      p(3) = 0.d0 
	p(4) = 0.d0 
      p(5) = 0.d0 
	p(6) = 0.d0 
      p(7) = 0.d0 
	p(8) = 0.d0
	k = 0
c ......................................................................
      do 260 lx = 1, nint
         ri = pg(lx,nint)
         do 250 ly = 1, nint
            si = pg(ly,nint)
            k = k + 1
            call sfquad4(h,hx,hy,ri,si,.true.,.true.)
c ...........
            call check_var_term(a,b,c,d11,d12,d22,d33,h,ut,ym,pr,alpha,
     .                        4,ndm,ma,flaghidr,.true.)
c            if (flaghidr .neqv. .true.) then
c               tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c               call updprop(ym,tm,ma,1)
c               call updprop(pr,tm,ma,2)
c               call updprop(alpha,tm,ma,4)
c               b1 = 1.d0-pr
c               b2 = 1.d0-2.d0*pr
c               a = ym*b1/((1.+pr)*b2)
c               b = pr/b1
c               c = b2/(2.d0*b1)
c               d11 = a
c               d12 = a*b
c               d22 = a
c               d33 = a*c
c            endif
c ...........
            call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)
            wt = wg(lx,nint)*wg(ly,nint)*det
            if (nlit .eq. 1) then
c .....     Deformacoes
               call deform2d(hx,hy,u,epsi,4)
c .....     Deformacoes Termicas:
               tm = h(1)*tl(1)+h(2)*tl(2)+h(3)*tl(3)+h(4)*tl(4)
               epsi(1) = epsi(1) - alpha*tm
               epsi(2) = epsi(2) - alpha*tm
c .....     Tensoes
              tx1  = d11*epsi(1) + d12*epsi(2) - d12*alpha*tm - txp(1,k)
              tx2  = d12*epsi(1) + d22*epsi(2) - d12*alpha*tm - txp(3,k)
              tx4  = d33*epsi(3) - txp(2,k)
c .....     Forcas Internas:
               p(1) = p(1) + (hx(1)*tx1 + hy(1)*tx4)*wt
               p(2) = p(2) + (hy(1)*tx2 + hx(1)*tx4)*wt
               p(3) = p(3) + (hx(2)*tx1 + hy(2)*tx4)*wt
               p(4) = p(4) + (hy(2)*tx2 + hx(2)*tx4)*wt
               p(5) = p(5) + (hx(3)*tx1 + hy(3)*tx4)*wt
               p(6) = p(6) + (hy(3)*tx2 + hx(3)*tx4)*wt
               p(7) = p(7) + (hx(4)*tx1 + hy(4)*tx4)*wt
               p(8) = p(8) + (hy(4)*tx2 + hx(4)*tx4)*wt
            else
               call deform2d0(hx,hy,u,0.d0,epsi,4)
               deps(1) = epsi(1) - eps(1,k) - alpha*tm
               deps(2) = epsi(2) - eps(2,k) - alpha*tm
               deps(3) = epsi(3) - eps(3,k) - alpha*tm
               deps(4) = epsi(4) - eps(4,k)
               eps(1,k) = eps(1,k)+deps(1)
               eps(2,k) = eps(2,k)+deps(2)
               eps(3,k) = eps(3,k)+deps(3)
               eps(4,k) = eps(4,k)+deps(4)
c              
               txl(1) = tx(1,k)
               txl(2) = tx(4,k)
               txl(3) = tx(2,k)
               txl(4) = tx(3,k)
               txl(5) = tx(5,k)
               c14 = pr/(1.d0-2.d0*pr)
               g = ym/(2.d0*(1.d0+pr))
              call plasticst2d(a,b,c,deps,eps(1,k),txl,txp(1,k),yield,
     .             hl,phi,pr,iyied,c14,g,.true.)
               tx(1,k) = txl(1)
               tx(2,k) = txl(3)
               tx(3,k) = txl(4)
               tx(4,k) = txl(2)
               tx(5,k) = txl(5)
               p(1) = p(1) + (hx(1)*tx(1,k)+hy(1)*tx(4,k)) * wt
               p(2) = p(2) + (hy(1)*tx(2,k)+hx(1)*tx(4,k)) * wt
               p(3) = p(3) + (hx(2)*tx(1,k)+hy(2)*tx(4,k)) * wt
               p(4) = p(4) + (hy(2)*tx(2,k)+hx(2)*tx(4,k)) * wt
               p(5) = p(5) + (hx(3)*tx(1,k)+hy(3)*tx(4,k)) * wt
               p(6) = p(6) + (hy(3)*tx(2,k)+hx(3)*tx(4,k)) * wt
               p(7) = p(7) + (hx(4)*tx(1,k)+hy(4)*tx(4,k)) * wt
               p(8) = p(8) + (hy(4)*tx(2,k)+hx(4)*tx(4,k)) * wt
            endif
 250     continue
 260  continue
c      call lku(s,u,p,nst)  
      return
c ======================================================================
c
c.... Tensoes nodais:       
c ......................................................................
  300 continue
c      yield = e(11)
c      b = yield*yield
      do 320 j = 1, 4
c ......................................................................
c
c...  Transforma as tensoes nos pontos de integracao em tensoes nodais:
c ......................................................................
	   k = (j-1)*5+1
         do 310 i = 1, 4
            c1 =.250d0*(tx(i,1)+tx(i,2)+tx(i,3)+tx(i,4))
            c2 =.433d0*(tx(i,3)-tx(i,1)-tx(i,2)+tx(i,4))
            c3 =.433d0*(tx(i,2)-tx(i,1)-tx(i,3)+tx(i,4))
            c4 =.750d0*(tx(i,1)-tx(i,2)-tx(i,3)+tx(i,4))
            p(k) = c1 + c2*rn(j) + c3*sn(j) + c4*rn(j)*sn(j)
	      k = k + 1
  310 continue
c ......................................................................
c
c...  Verifica se as tensoes nos pontos de integracao em tensoes nodais:
c ......................................................................
c	   k = (j-1)*5+1
c	   a = 3.d0*invar2(p(k),4)
c	   if(a .gt. b) then
c            l = pint(j)
c	      c = 3.d0*invar2(tx(1,l),4)
c            if(c .ge. b) then
c               p(k  ) = tx(1,l)
c	         p(k+1) = tx(2,l)
c	         p(k+2) = tx(3,l)
c	         p(k+3) = tx(4,l)
c	      else
c               c = dsqrt(b/a)
c               p(k  ) = p(k  )*c
c	         p(k+1) = p(k+1)*c
c	         p(k+2) = p(k+2)*c
c	         p(k+3) = p(k+3)*c
c            endif
c	    endif
  320 continue
c      call sfquad4(h,hx,hy,1.d0,1.d0,.true.,.true.)
c ...........
c            if (flaghidr .neqv. .true.) then
c               tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c               call updprop(ym,tm,ma,1)
c               call updprop(pr,tm,ma,2)
c               call updprop(alpha,tm,ma,4)
c               b1 = 1.d0-pr
c               b2 = 1.d0-2.d0*pr
c               a = ym*b1/((1.+pr)*b2)
c               b = pr/b1
c               c = b2/(2.d0*b1)
c               d11 = a
c               d12 = a*b
c               d22 = a
c               d33 = a*c
c            endif
cc ...........
c      call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)  
c      call deform2d(hx,hy,u,eps,4)
c      eps(1) = eps(1) - alpha*tl(1)
c      eps(2) = eps(2) - alpha*tl(1)
c      p(1) = d11*eps(1)+d12*eps(2)-d12*alpha*tl(1)
c      p(2) = d12*eps(1)+d22*eps(2)-d12*alpha*tl(1)
c      p(3) = d33*eps(3)
cc      call stress2d(d11,d12,d22,d33,eps,p(1))
c      call sfquad4(h,hx,hy,-1.d0,1.d0,.true.,.true.)
c ...........
c            if (flaghidr .neqv. .true.) then
c               tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c               call updprop(ym,tm,ma,1)
c               call updprop(pr,tm,ma,2)
c               call updprop(alpha,tm,ma,4)
c               b1 = 1.d0-pr
c               b2 = 1.d0-2.d0*pr
c               a = ym*b1/((1.+pr)*b2)
c               b = pr/b1
c               c = b2/(2.d0*b1)
c               d11 = a
c               d12 = a*b
c               d22 = a
c               d33 = a*c
c            endif
cc ...........
c      call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)   
c      call deform2d(hx,hy,u,eps,4)
c      eps(1) = eps(1) - alpha*tl(2)
c      eps(2) = eps(2) - alpha*tl(2)
c      p(4) = d11*eps(1)+d12*eps(2)-d12*alpha*tl(2)
c      p(5) = d12*eps(1)+d22*eps(2)-d12*alpha*tl(2)
c      p(6) = d33*eps(3)
cc      call stress2d(d11,d12,d22,d33,eps,p(4))
c      call sfquad4(h,hx,hy,-1.d0,-1.d0,.true.,.true.) 
cc ...........
c            if (flaghidr .neqv. .true.) then
c               tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c               call updprop(ym,tm,ma,1)
c               call updprop(pr,tm,ma,2)
c               call updprop(alpha,tm,ma,4)
c               b1 = 1.d0-pr
c               b2 = 1.d0-2.d0*pr
c               a = ym*b1/((1.+pr)*b2)
c               b = pr/b1
c               c = b2/(2.d0*b1)
c               d11 = a
c               d12 = a*b
c               d22 = a
c               d33 = a*c
c            endif
c ...........
c      call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)   
c      call deform2d(hx,hy,u,eps,4)
c      eps(1) = eps(1) - alpha*tl(3)
c      eps(2) = eps(2) - alpha*tl(3)
c      p(7) = d11*eps(1)+d12*eps(2)-d12*alpha*tl(3)
c      p(8) = d12*eps(1)+d22*eps(2)-d12*alpha*tl(3)
c      p(9) = d33*eps(3)
cc      call stress2d(d11,d12,d22,d33,eps,p(7)) 
c      call sfquad4(h,hx,hy,1.d0,-1.d0,.true.,.true.) 
cc ...........
c            if (flaghidr .neqv. .true.) then
c               tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c               call updprop(ym,tm,ma,1)
c               call updprop(pr,tm,ma,2)
c               call updprop(alpha,tm,ma,4)
c               b1 = 1.d0-pr
c               b2 = 1.d0-2.d0*pr
c               a = ym*b1/((1.+pr)*b2)
c               b = pr/b1
c               c = b2/(2.d0*b1)
c               d11 = a
c              d12 = a*b
c               d22 = a
c               d33 = a*c
c            endif
cc ...........
c      call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)   
c      call deform2d(hx,hy,u,eps,4)
c      eps(1) = eps(1) - alpha*tl(4)
c      eps(2) = eps(2) - alpha*tl(4)
c      p(10) = d11*eps(1)+d12*eps(2)-d12*alpha*tl(4)
c      p(11) = d12*eps(1)+d22*eps(2)-d12*alpha*tl(4)
c      p(12) = d33*eps(3)
cc      call stress2d(d11,d12,d22,d33,eps,p(10))        
      return
  400 continue 
      return
c ....................................................................
      end
      subroutine elmt05_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,
     .                                           ma,nlit,isw,flaghidr)
c **********************************************************************
c *                                                                    *
c *   ELMT05: Elemento quadrilatero de 4 nos                           *
c *   ------  estado plano de tensao                                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      common /gauss/ pg,wg
      integer ndm,nst,nel,isw,iyied
      integer nint,i,j,k,l,lx,ly,ma,nlit,pint(4)
      real*8  e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*),ut(*)
      real*8  h(4),hx(4),hy(4),xj(2,2),xji(2,2),det,wt
      real*8  ym,pr,thic,a,d11,d12,d22,d33,tm,alpha,hi(3),txl(5)
      real*8  pg(10,10),wg(10,10),ri,si,b,c,yield,c1,c2,c3,c4,invar2
      real*8  eps(5,*),tx(5,*),epsi(4),deps(4),tx1,tx2,tx4,rn(4),sn(4)
      real*8  phi,coes,hl,c14,g,txp(5,*)
      real*8  steff,theta,varj2,sint3,snphi,devs(4)
      logical flaghidr
      data    rn /1.d0,-1.d0,-1.d0,1.d0/,sn /1.d0,1.d0,-1.d0,-1.d0/
      data    nint/2/,pint/4,2,1,3/
c ......................................................................
      goto (100,200,300,400) isw
c ======================================================================
c
c.... Imput material properties:
c
  100 continue
c     e(1) = modulo de elasticidade
c     e(2) = coeficiente de Poisson
c     e(3) = espessura
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
c
c ... Matriz constitutiva:
c ......................................................................
      ym    = e(1)
      pr    = e(2)
      thic  = e(3)
      alpha = e(4)
      iyied =int(e(5))
      yield = e(11)
      phi   = e(12)
      coes  = e(13)
      hl    = e(14)
      call radian(phi)
c      snphi=dsin(phi)
      if( iyied .eq. 3)   yield = coes*dcos(phi)
      if( iyied .eq. 4)   yield = 6.d0*coes*dcos(phi)/(dsqrt(3.d0)*
     .    (3.d0-dsin(phi)))
      call check_var_hidr(a,b,c,d11,d12,d22,d33,ym,pr,alpha,hi(2),
     .                           ma,ndm,flaghidr,.false.)
c      if (flaghidr .eqv. .true.) then
c         call updhidr(ym,hi(2),ma,1)
c         call updhidr(pr,hi(2),ma,2)
c         call updhidr(alpha,hi(2),ma,4)
c      endif
c ......................................................................
c      a = ym/(1.d0-pr*pr)
c      d11 = a
c      d12 = a*pr
c      d22 = a
c      d33 = ym/(2.d0*(1.d0+pr))
c      if (isw .eq. 3) goto 300
c ......................................................................
c
c ... Matriz de rigidez:
c ......................................................................
      do i = 1, nst
      do j = 1, nst
         s(i,j) = 0.d0
      enddo
      enddo
      do 240 lx = 1, nint
         ri = pg(lx,nint)
         do 230 ly = 1, nint
            si = pg(ly,nint)
            call sfquad4(h,hx,hy,ri,si,.true.,.true.)
c ...........
            call check_var_term(a,b,c,d11,d12,d22,d33,h,ut,ym,pr,alpha,
     .                        4,ndm,ma,flaghidr,.false.)
c            if (flaghidr .neqv. .true.) then
c               tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c               call updprop(ym,tm,ma,1)
c               call updprop(pr,tm,ma,2)
c               a = ym/(1.d0-pr*pr)
c               d11 = a
c               d12 = a*pr
c               d22 = a
c               d33 = ym/(2.d0*(1.d0+pr))
c            endif
c ...........
            call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)
            wt = wg(lx,nint)*wg(ly,nint)*det*thic
        do 220 j = 1, 4
               k = (j-1)*2+1
               do 210 i = 1, 4
                  l = (i-1)*2+1
c --------------------------------------------------------------------
                  s(l,k)     = s(l,k) +
     .                        (hx(i)*d11*hx(j)+hy(i)*d33*hy(j))*wt
c --------------------------------------------------------------------
                  s(l,k+1)   = s(l,k+1) +
     .                        (hx(i)*d12*hy(j)+hy(i)*d33*hx(j))*wt
c --------------------------------------------------------------------
                  s(l+1,k)   = s(l+1,k) +
     .                        (hy(i)*d12*hx(j)+hx(i)*d33*hy(j))*wt
c --------------------------------------------------------------------
                  s(l+1,k+1) = s(l+1,k+1) +
     .                        (hy(i)*d22*hy(j)+hx(i)*d33*hx(j))*wt
c --------------------------------------------------------------------
  210          continue
  220       continue
  230    continue
  240 continue
c
c ... Forcas internas:
c
      p(1) = 0.d0
      p(2) = 0.d0 
      p(3) = 0.d0 
      p(4) = 0.d0 
      p(5) = 0.d0 
      p(6) = 0.d0 
      p(7) = 0.d0 
      p(8) = 0.d0
      k = 0
c .....................................................................
      do 260 lx = 1, nint
         ri = pg(lx,nint)
         do 250 ly = 1, nint
            si = pg(ly,nint)
            call sfquad4(h,hx,hy,ri,si,.true.,.true.)
c ...........
            call check_var_term(a,b,c,d11,d12,d22,d33,h,ut,ym,pr,alpha,
     .                        4,ndm,ma,flaghidr,.false.)
c            if (flaghidr .neqv. .true.) then
c               tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c               call updprop(ym,tm,ma,1)
c               call updprop(pr,tm,ma,2)
c               call updprop(alpha,tm,ma,4)
c               a = ym/(1.d0-pr*pr)
c               d11 = a
c               d12 = a*pr
c               d22 = a
c               d33 = ym/(2.d0*(1.d0+pr))
c            endif
c ...........
            call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)
            wt = wg(lx,nint)*wg(ly,nint)*det*thic
            tm = h(1)*tl(1)+h(2)*tl(2)+h(3)*tl(3)+h(4)*tl(4)
            if ( nlit .eq. 1) then
               k = k + 1
c .....     Deformacoes
               call deform2d(hx,hy,u,epsi,4)
c .....     Deformacoes Termicas:
               epsi(1) = epsi(1) - alpha*tm
               epsi(2) = epsi(2) - alpha*tm
c .....     Tensoes
               tx1 = d11*epsi(1) + d12*epsi(2) - txp(1,k)
               tx2 = d12*epsi(1) + d22*epsi(2) - txp(3,k)
               tx4 = d33*epsi(3) - txp(2,k)
c .....     Forcas Internas:
               p(1) = p(1) + (hx(1)*tx1 + hy(1)*tx4)*wt
               p(2) = p(2) + (hy(1)*tx2 + hx(1)*tx4)*wt
               p(3) = p(3) + (hx(2)*tx1 + hy(2)*tx4)*wt
               p(4) = p(4) + (hy(2)*tx2 + hx(2)*tx4)*wt
               p(5) = p(5) + (hx(3)*tx1 + hy(3)*tx4)*wt
               p(6) = p(6) + (hy(3)*tx2 + hx(3)*tx4)*wt
               p(7) = p(7) + (hx(4)*tx1 + hy(4)*tx4)*wt
               p(8) = p(8) + (hy(4)*tx2 + hx(4)*tx4)*wt
            else 
               k = k + 1
c               a = ym*(1.d0-pr)/((1.d0+pr)*(1.d0-2.d0*pr))
c               b = pr/(1.d0-pr)
c               c = (1.d0-2.d0*pr)/(2.d0*(1.d0-pr))
                a = ym/(1.d0-pr*pr)
                b = pr
                c = (1.d0-pr*pr)/(2.d0+2.d0*pr)
               call deform2d0(hx,hy,u,0.d0,epsi,4)
               deps(1) = epsi(1) - eps(1,k) - alpha*tm
               deps(2) = epsi(2) - eps(2,k) - alpha*tm
c              deps(3) = epsi(3) - eps(3,k) - alpha*tm
               deps(4) = epsi(4) - eps(4,k)
               eps(1,k)=eps(1,k)+deps(1)
               eps(2,k)=eps(2,k)+deps(2)
               eps(3,k)=eps(3,k)+deps(3)
               eps(4,k)=eps(4,k)+deps(4)
c               
               txl(1) = tx(1,k)
               txl(2) = tx(4,k)
               txl(3) = tx(2,k)
               txl(4) = tx(3,k)
               txl(5) = tx(5,k)
               c14 = pr/(1.d0-pr)
               g = ym/(2.d0*(1.d0+pr))
              call plasticst2d(a,b,c,deps,eps(1,k),txl,txp(1,k),yield,
     .          hl,phi,0.d0,iyied,c14,g,.false.)
               tx(1,k) = txl(1)
               tx(2,k) = txl(3)
               tx(3,k) = txl(4)
               tx(4,k) = txl(2)
               tx(5,k) = txl(5)
               p(1) = p(1) + (hx(1)*tx(1,k)+hy(1)*tx(4,k)) * wt
               p(2) = p(2) + (hy(1)*tx(2,k)+hx(1)*tx(4,k)) * wt
               p(3) = p(3) + (hx(2)*tx(1,k)+hy(2)*tx(4,k)) * wt
               p(4) = p(4) + (hy(2)*tx(2,k)+hx(2)*tx(4,k)) * wt
               p(5) = p(5) + (hx(3)*tx(1,k)+hy(3)*tx(4,k)) * wt
               p(6) = p(6) + (hy(3)*tx(2,k)+hx(3)*tx(4,k)) * wt
               p(7) = p(7) + (hx(4)*tx(1,k)+hy(4)*tx(4,k)) * wt
               p(8) = p(8) + (hy(4)*tx(2,k)+hx(4)*tx(4,k)) * wt
            endif
 250     continue
 260  continue
c      call lku(s,u,p,nst)
      return  
c ......................................................................
  300 continue
c ======================================================================
c.... Tensoes nodais:  
c      iyied =int(e(5))     
c  	yield = e(11)
c      phi   = e(12)
c      call radian(phi)
c      snphi = dsin(phi)
c      b = yield*yield
      do 320 j = 1, 4
c ......................................................................
c
c...  Transforma as tensoes nos pontos de integracao em tensoes nodais:
c ......................................................................
	   k = (j-1)*5+1
         do 310 i = 1, 4
            c1 =.250d0*(tx(i,1)+tx(i,2)+tx(i,3)+tx(i,4))
            c2 =.433d0*(tx(i,3)-tx(i,1)-tx(i,2)+tx(i,4))
            c3 =.433d0*(tx(i,2)-tx(i,1)-tx(i,3)+tx(i,4))
            c4 =.750d0*(tx(i,1)-tx(i,2)-tx(i,3)+tx(i,4))
            p(k) = c1 + c2*rn(j) + c3*sn(j) + c4*rn(j)*sn(j)
	      k = k + 1
  310    continue
         if (iyied .eq. 0)  goto 320
c ......................................................................
c
c...  Verifica se as tensoes nos pontos de integracao em tensoes nodais:
c ......................................................................
c	   k = (j-1)*5+1
c         call effst(a,iyied,p(k),devs,steff,theta,varj2,sint3,snphi)
c	   a = 3.d0*invar2(p(k),4)
c	   if(a .gt. b) then
c            l = pint(j)
c          call effst(c,iyied,tx(1,l),devs,steff,theta,varj2,sint3,snphi)
cc	      c = 3.d0*invar2(tx(1,l),4)
c            if(c .ge. b) then
c               p(k  ) = tx(1,l)
c	         p(k+1) = tx(2,l)
c	         p(k+2) = tx(3,l)
c	         p(k+3) = tx(4,l)
c	      else
c               c = dsqrt(b/a)
c               p(k  ) = p(k  )*c
c	         p(k+1) = p(k+1)*c
c	         p(k+2) = p(k+2)*c
c	         p(k+3) = p(k+3)*c
c            endif
c	    endif
  320 continue
c ......................................................................
c      call sfquad4(h,hx,hy,1.d0,1.d0,.true.,.true.)
c ...........
c         if (flaghidr .neqv. .true.) then
c            tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c            call updprop(ym,tm,ma,1)
c            call updprop(pr,tm,ma,2)
c            call updprop(alpha,tm,ma,4)
c            a = ym/(1.d0-pr*pr)
c            d11 = a
c            d12 = a*pr
c            d22 = a
c            d33 = ym/(2.d0*(1.d0+pr))
c         endif
c ...........
c      call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)  
c      call deform2d(hx,hy,u,eps,4)
c      eps(1) = eps(1) -alpha*tl(1)
c      eps(2) = eps(2) -alpha*tl(1)
c      call stress2d(d11,d12,d22,d33,eps,p(1))
c      call sfquad4(h,hx,hy,-1.d0,1.d0,.true.,.true.)
cc ...........
c         if (flaghidr .neqv. .true.) then
c            tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c            call updprop(ym,tm,ma,1)
c            call updprop(pr,tm,ma,2)
c            call updprop(alpha,tm,ma,4)
c            a = ym/(1.d0-pr*pr)
c            d11 = a
c            d12 = a*pr
c            d22 = a
c            d33 = ym/(2.d0*(1.d0+pr))
c         endif
cc ...........      
c      call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)   
c      call deform2d(hx,hy,u,eps,4)
c      eps(1) = eps(1) - alpha*tl(2)
c      eps(2) = eps(2) - alpha*tl(2)
c      call stress2d(d11,d12,d22,d33,eps,p(4))
c      call sfquad4(h,hx,hy,-1.d0,-1.d0,.true.,.true.) 
cc ...........
c         if (flaghidr .neqv. .true.) then
c            tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c            call updprop(ym,tm,ma,1)
c            call updprop(pr,tm,ma,2)
c            call updprop(alpha,tm,ma,4)
c            a = ym/(1.d0-pr*pr)
c            d11 = a
c            d12 = a*pr
c            d22 = a
c            d33 = ym/(2.d0*(1.d0+pr))
c         endif
cc ...........         
c      call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)   
c      call deform2d(hx,hy,u,eps,4)
c      eps(1) = eps(1) - alpha*tl(3)
c      eps(2) = eps(2) - alpha*tl(3)
c      call stress2d(d11,d12,d22,d33,eps,p(7)) 
c      call sfquad4(h,hx,hy,1.d0,-1.d0,.true.,.true.) 
cc ...........
c         if (flaghidr .neqv. .true.) then
c            tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)
c            call updprop(ym,tm,ma,1)
c            call updprop(pr,tm,ma,2)
c            call updprop(alpha,tm,ma,4)
c            a = ym/(1.d0-pr*pr)
c            d11 = a
c            d12 = a*pr
c            d22 = a
c            d33 = ym/(2.d0*(1.d0+pr))
c         endif
c ...........          
c      call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)   
c      call deform2d(hx,hy,u,eps,4)
c      eps(1) = eps(1) - alpha*tl(4)
c      eps(2) = eps(2) - alpha*tl(4)
c      call stress2d(d11,d12,d22,d33,eps,p(10))    
      return
c ......................................................................
  400 continue
      return
c ......................................................................
      end
      subroutine elmt06_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,
     .    ma,nlit,isw,flaghidr,plastic)
c **********************************************************************
c *                                                                    *
c *   ELMT06: Elemento tetraedrico                                     *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - solucao anterior                                  *
c *     v(nst)     - solucao anterior incremental                      *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     tx(ntn)    - tensoes no elemento                               *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'            
      integer ndm,nst,nel,isw,nlit,flag,i1,i2,i3,iyied
      integer i,j,k,l,ma,druckerprager3d,rankine
      real*8 e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*),ut(*)
      real*8 hx(4),hy(4),hz(4),det,wt,alpha,hi(3)
      real*8 xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33,tm,dum
      real*8 xji11,xji12,xji13,xji21,xji22,xji23,xji31,xji32,xji33
      real*8 ym,nu,a,b,c,a1,a2,a3,sigma,d(6,6),h(4),epsi(6)
      real*8 phi,fc,ft,plastic(2),deps(6),eps(*),tx(*),txl(7),txp(*)
      real*8 yield, coes,c14,g,hl
      logical flaghidr
c ......................................................................
      go to (100,200,300,200,200) isw
c ======================================================================
c
c.... Input material properties:
c
c ......................................................................
  100 continue
c     e(1) = modulo de elasticidade
c     e(2) = coeficiente de Poisson
c     e(3) = sem uso
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
c
c ... Matriz Jacobiana:
c
      xj11 = x(1,1)-x(1,4)
      xj12 = x(2,1)-x(2,4)
      xj13 = x(3,1)-x(3,4)
      xj21 = x(1,2)-x(1,4)
      xj22 = x(2,2)-x(2,4)
      xj23 = x(3,2)-x(3,4)
      xj31 = x(1,3)-x(1,4)
      xj32 = x(2,3)-x(2,4)
      xj33 = x(3,3)-x(3,4)
      det  = xj11*xj22*xj33 + xj12*xj23*xj31 + xj13*xj21*xj32
     .      -xj31*xj22*xj13 - xj12*xj21*xj33 - xj11*xj32*xj23
c ......................................................................      
      if (det .le. 0.d0) go to 1000
c 
c ..  Inversa da matriz Jacobiana:
c
      xji11 = (xj22*xj33-xj23*xj32)/det
      xji12 = (xj13*xj32-xj12*xj33)/det
      xji13 = (xj12*xj23-xj13*xj22)/det
      xji21 = (xj23*xj31-xj21*xj33)/det
      xji22 = (xj11*xj33-xj13*xj31)/det
      xji23 = (xj13*xj21-xj11*xj23)/det
      xji31 = (xj21*xj32-xj22*xj31)/det
      xji32 = (xj12*xj31-xj11*xj32)/det
      xji33 = (xj11*xj22-xj12*xj21)/det
c
c ... Derivadas das funcoes de interpolacao:
c
      hx(1) =  xji11
      hy(1) =  xji21
      hz(1) =  xji31
      hx(2) =  xji12
      hy(2) =  xji22
      hz(2) =  xji32
      hx(3) =  xji13
      hy(3) =  xji23
      hz(3) =  xji33
      hx(4) = -xji11-xji12-xji13
      hy(4) = -xji21-xji22-xji23
      hz(4) = -xji31-xji32-xji33
      if (isw .eq. 4) goto 400
      if (isw .eq. 5) goto 500
c ......................................................................
c ... Propriedades do material:
c ......................................................................
      ym      =  e(1)
      nu      =  e(2)
      alpha   =  e(4)
      iyied =int(e(5))
      yield = e(11)
      phi   = e(12)
      coes  = e(13)
      hl    = e(14)
      call radian(phi)
      if( iyied .eq. 3)   yield = coes*dcos(phi)
      if( iyied .eq. 4)   yield = 6.d0*coes*dcos(phi)/(dsqrt(3.d0)*
     .    (3.d0-dsin(phi)))
      h(1) = 1.d0/4.d0
      h(2) = h(1)
      h(3) = h(1)
      h(4) = h(1)
      call check_var_hidr(a,b,c,dum,dum,dum,dum,ym,nu,alpha,hi(2),
     .                           ma,ndm,flaghidr,.false.)
      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,nu,alpha,
     .                        4,ndm,ma,flaghidr,.false.)
c      if (flaghidr .eqv. .true.) then
c         call updhidr(ym,hi(2),ma,1)
c         call updhidr(nu,hi(2),ma,2)
c         call updhidr(alpha,hi(2),ma,4)
c         hi(3) = ym
c      else if (flaghidr .neqv. .true.) then
c         call updprop(ym,tm,ma,1)
c         call updprop(nu,tm,ma,2)
c         call updprop(alpha,tm,ma,4)
c      endif
c
c ... Matriz constitutiva:
c
c      a1 = 1.d0+nu
c      a2 = 1.d0-nu
c      a3 = 1.d0-2.d0*nu
c
c      a = ym*a2/(a1*a3)
c      b = nu/a2
c      c = a3/(2.d0*a2)
c      if (isw .eq. 3) goto 300  
c ......................................................................
c
c ... Matriz de rigidez:
c
c      wt = a*det/6.d0
c      do 220 i = 1, 4
c         l = (i-1)*3+1
c         do 210 j = 1, 4
c         k = (j-1)*3+1
c         s(l  ,  k) = (hx(i)*hx(j) + c*hy(i)*hy(j) + c*hz(i)*hz(j))*wt
c         s(l  ,k+1) = (b*hx(i)*hy(j) + c*hy(i)*hx(j))*wt
c         s(l  ,k+2) = (b*hx(i)*hz(j) + c*hz(i)*hx(j))*wt
c         s(l+1,  k) = (b*hy(i)*hx(j) + c*hx(i)*hy(j))*wt
c         s(l+1,k+1) = (hy(i)*hy(j) + c*hx(i)*hx(j) + c*hz(i)*hz(j))*wt
c         s(l+1,k+2) = (b*hy(i)*hz(j) + c*hz(i)*hy(j))*wt
c         s(l+2,  k) = (b*hz(i)*hx(j) + c*hx(i)*hz(j))*wt
c         s(l+2,k+1) = (b*hz(i)*hy(j) + c*hy(i)*hz(j))*wt
c         s(l+2,k+2) = (hz(i)*hz(j) + c*hy(i)*hy(j) + c*hx(i)*hx(j))*wt
c  210    continue
c  220 continue
c ... deformacoes:      
      call deform3d(hx,hy,hz,u,epsi,4)
      wt = a*det/6.d0
      do 215 i = 1, 4
      l = (i-1)*3+1
         do 205 j = 1, 4
           k = (j-1)*3+1
           s(l  ,  k) = (hx(i)*hx(j) + c*hy(i)*hy(j) + c*hz(i)*hz(j))*wt
           s(l  ,k+1) = (b*hx(i)*hy(j) + c*hy(i)*hx(j))*wt
           s(l  ,k+2) = (b*hx(i)*hz(j) + c*hz(i)*hx(j))*wt
           s(l+1,  k) = (b*hy(i)*hx(j) + c*hx(i)*hy(j))*wt
           s(l+1,k+1) = (hy(i)*hy(j) + c*hx(i)*hx(j) + c*hz(i)*hz(j))*wt
           s(l+1,k+2) = (b*hy(i)*hz(j) + c*hz(i)*hy(j))*wt
           s(l+2,  k) = (b*hz(i)*hx(j) + c*hx(i)*hz(j))*wt
           s(l+2,k+1) = (b*hz(i)*hy(j) + c*hy(i)*hz(j))*wt
           s(l+2,k+2) = (hz(i)*hz(j) + c*hy(i)*hy(j) + c*hx(i)*hx(j))*wt
 205    continue
 215  continue
         tm = (tl(1)+tl(2)+tl(3)+tl(4))/4.d0
      if (nlit .eq. 1) then
c ... deformacoes termicas:
         epsi(1) = epsi(1) - alpha*tm
         epsi(2) = epsi(2) - alpha*tm
         epsi(3) = epsi(3) - alpha*tm
c ... tensoes:
         txl(1)  =  (  epsi(1) + b*epsi(2) + b*epsi(3))*a !- txp(1)
         txl(2)  =  (b*epsi(1) +   epsi(2) + b*epsi(3))*a !- txp(2)
         txl(3)  =  (b*epsi(1) + b*epsi(2) +   epsi(3))*a !- txp(3)
         txl(4) =   (c * epsi(4))*a !- txp(4)
         txl(5) =   (c * epsi(5))*a !- txp(5)
         txl(6) =   (c * epsi(6))*a !- txp(6)
c         txl(1)  = a*epsi(1) + b*epsi(2) + b*epsi(3)
c         txl(2) = b*epsi(1) + a*epsi(2) + b*epsi(3)
c         txl(3)  = b*epsi(1) + b*epsi(2) + a*epsi(3)
c         txl(4) = (1.d0-2.d0*nu)  / (2.d0 * (1.d0-nu)) * epsi(4)
c         txl(5) = (1.d0-2.d0*nu)  / (2.d0 * (1.d0-nu)) * epsi(5)
c         txl(6) = (1.d0-2.d0*nu)  / (2.d0 * (1.d0-nu)) * epsi(6)
c         do i = 1, 4
c         i1 = (i-1)*3+1
c         i2 = i1 + 1
c         i3 = i2 + 1
c         p(i1) = (hx(i)*txl(1)+hy(i)*txl(4)+hz(i)*txl(6)) * wt
c         p(i2) = (hy(i)*txl(2)+hx(i)*txl(4)+hz(i)*txl(5)) * wt
c         p(i3) = (hz(i)*txl(3)+hx(i)*txl(6)+hy(i)*txl(5)) * wt
c         enddo
c ..................................................................
      else
         deps(1) = epsi(1)-eps(1)-alpha*tm
         deps(2) = epsi(2)-eps(2)-alpha*tm
         deps(3) = epsi(3)-eps(3)-alpha*tm
         deps(4) = epsi(4)-eps(4)
         deps(5) = epsi(5)-eps(5)
         deps(6) = epsi(6)-eps(6)
         do 235 i = 1, 6
           eps(i) = eps(i)+deps(i)
  235    continue 
         do 240 i = 1, 7
            txl(i) = tx(i)
  240    continue 
         c14 = nu/(1.d0-2.d0*nu)
         g = ym/(2.d0*(1.d0+nu))      
         call plasticst3d(a,b,c,deps,eps,txl,txp,yield,hl,phi,nu,
     .          iyied,c14,g)
         do 250 i = 1, 7
            tx(i) = txl(i)
  250    continue 
      endif
c         call plasticity(deps,a,b,c,mu,lamb,yield,phi,coes,10,4,txl,
c     .                dsig,flag)
c         effect = effectivet(txl,phi,coes)
c         call matrizd(a,b,c,txl,phi,d,flag)              
c        wt = det/6.d0
c        do 220 i = 1, 4
c            l = (i-1)*3+1
c            do 210 j = 1, 4
c            k = (j-1)*3+1
c   	   s(l  ,  k) = (hx(i)*hx(j) + c*hy(i)*hy(j) + c*hz(i)*hz(j))*wt
c	   s(l  ,k+1) = (b*hx(i)*hy(j) + c*hy(i)*hx(j))*wt
c	   s(l  ,k+2) = (b*hx(i)*hz(j) + c*hz(i)*hx(j))*wt
c	   s(l+1,  k) = (b*hy(i)*hx(j) + c*hx(i)*hy(j))*wt
c	   s(l+1,k+1) = (hy(i)*hy(j) + c*hx(i)*hx(j) + c*hz(i)*hz(j))*wt
c	   s(l+1,k+2) = (b*hy(i)*hz(j) + c*hz(i)*hy(j))*wt
c	   s(l+2,  k) = (b*hz(i)*hx(j) + c*hx(i)*hz(j))*wt
c	   s(l+2,k+1) = (b*hz(i)*hy(j) + c*hy(i)*hz(j))*wt
c	   s(l+2,k+2) = (hz(i)*hz(j) + c*hy(i)*hy(j) + c*hx(i)*hx(j))*wt
c	   s(l  ,  k) =((hx(i)*d(1,1)+hy(i)*d(4,1)+hz(i)*d(6,1))*hx(j) +
c     .                (hx(i)*d(1,4)+hy(i)*d(4,4)+hz(i)*d(6,4))*hy(j) +
c     .                (hx(i)*d(1,6)+hy(i)*d(4,6)+hz(i)*d(6,6))*hz(j))*wt
c         s(l  ,k+1) =((hx(i)*d(1,2)+hy(i)*d(4,2)+hz(i)*d(6,2))*hy(j) +
c     .                (hx(i)*d(1,4)+hy(i)*d(4,4)+hz(i)*d(6,4))*hx(j) +
c     .                (hx(i)*d(1,5)+hy(i)*d(4,5)+hz(i)*d(6,5))*hz(j))*wt
c         s(l  ,k+2) =((hx(i)*d(1,3)+hy(i)*d(4,3)+hz(i)*d(6,3))*hz(j) +
c     .                (hx(i)*d(1,5)+hy(i)*d(4,5)+hz(i)*d(6,5))*hy(j) +
c     .                (hx(i)*d(1,6)+hy(i)*d(4,6)+hz(i)*d(6,6))*hx(j))*wt
c         s(l+1,  k) =((hy(i)*d(2,1)+hx(i)*d(4,1)+hz(i)*d(5,1))*hx(j) +
c     .                (hy(i)*d(2,4)+hx(i)*d(4,4)+hz(i)*d(5,4))*hy(j) +
c     .                (hy(i)*d(2,6)+hx(i)*d(4,6)+hz(i)*d(5,6))*hz(j))*wt
c	   s(l+1,k+1) =((hy(i)*d(2,2)+hx(i)*d(4,2)+hz(i)*d(5,2))*hy(j) +
c     .                (hy(i)*d(2,4)+hx(i)*d(4,4)+hz(i)*d(5,4))*hx(j) +
c     .                (hy(i)*d(2,5)+hx(i)*d(4,5)+hz(i)*d(5,5))*hz(j))*wt
c	   s(l+1,k+2) =((hy(i)*d(2,3)+hx(i)*d(4,3)+hz(i)*d(5,3))*hz(j) +
c     .                (hy(i)*d(2,5)+hx(i)*d(4,5)+hz(i)*d(5,5))*hy(j) +
c     .                (hy(i)*d(2,6)+hx(i)*d(4,6)+hz(i)*d(5,6))*hx(j))*wt
c	   s(l+2,  k) =((hz(i)*d(3,1)+hy(i)*d(5,1)+hx(i)*d(6,1))*hx(j) +
c     .                (hz(i)*d(3,4)+hy(i)*d(5,4)+hx(i)*d(6,4))*hy(j) +
c     .                (hz(i)*d(3,6)+hy(i)*d(5,6)+hx(i)*d(6,6))*hz(j))*wt
c	   s(l+2,k+1) =((hz(i)*d(3,2)+hy(i)*d(5,2)+hx(i)*d(6,2))*hy(j) +
c     .                (hz(i)*d(3,4)+hy(i)*d(5,4)+hx(i)*d(6,4))*hx(j) +
c     .                (hz(i)*d(3,5)+hy(i)*d(5,5)+hx(i)*d(6,5))*hz(j))*wt
c	   s(l+2,k+2) =((hz(i)*d(3,3)+hy(i)*d(5,3)+hx(i)*d(6,3))*hz(j) +
c     .                (hz(i)*d(3,5)+hy(i)*d(5,5)+hx(i)*d(6,5))*hy(j) +
c     .                (hz(i)*d(3,6)+hy(i)*d(5,6)+hx(i)*d(6,6))*hx(j))*wt
c  210       continue
c  220    continue
c         tx(1) = tx(1) + dsig(1)
c         tx(2) = tx(2) + dsig(2)
c         tx(3) = tx(3) + dsig(3)
c         tx(4) = tx(4) + dsig(4)
c         tx(5) = tx(5) + dsig(5)
c         tx(6) = tx(6) + dsig(6)  
c      endif
c ... deformacoes termicas:
c      tm = (tl(1)+tl(2)+tl(3)+tl(4))/4.d0
c      eps(1) = eps(1) - alpha*tm
c      eps(2) = eps(2) - alpha*tm
c      eps(3) = eps(3) - alpha*tm
c ... tensoes:
c      tx  = a*eps(1) + b*eps(2) + b*eps(3)
c      ty  = b*eps(1) + a*eps(2) + b*eps(3)
c      tz  = b*eps(1) + b*eps(2) + a*eps(3)
c      txy = a3  / (2 * a2) * eps(4)
c      tyz = a3  / (2 * a2) * eps(5)
c      txz = a3  / (2 * a2) * eps(6)
c ...
c ... Forcas internas:
c
c ======================================================================
c
c ...Forcas internas:
c
c ......................................................................
c 222   continue
       wt = det/6.0d0
      do 230 i = 1, 4
         i1 = (i-1)*3+1
         i2 = i1 + 1
         i3 = i2 + 1
         p(i1) = (hx(i)*txl(1)+hy(i)*txl(4)+hz(i)*txl(6)) * wt
         p(i2) = (hy(i)*txl(2)+hx(i)*txl(4)+hz(i)*txl(5)) * wt
         p(i3) = (hz(i)*txl(3)+hx(i)*txl(6)+hy(i)*txl(5)) * wt
c ......................................................................
  230 continue
c      p(1)  = (hx(1)*tx + hy(1)*txy + hz(1)*txz)*wt
c      p(2)  = (hy(1)*ty + hx(1)*txy + hz(1)*tyz)*wt
c      p(3)  = (hz(1)*tz + hx(1)*txz + hy(1)*tyz)*wt
c      p(4)  = (hx(2)*tx + hy(2)*txy + hz(2)*txz)*wt
c      p(5)  = (hy(2)*ty + hx(2)*txy + hz(2)*tyz)*wt
c      p(6)  = (hz(2)*tz + hx(2)*txz + hy(2)*tyz)*wt
c      p(7)  = (hx(3)*tx + hy(3)*txy + hz(3)*txz)*wt
c      p(8)  = (hy(3)*ty + hx(3)*txy + hz(3)*tyz)*wt
c      p(9)  = (hz(3)*tz + hx(3)*txz + hy(3)*tyz)*wt
c      p(10) = (hx(4)*tx + hy(4)*txy + hz(4)*txz)*wt
c      p(11) = (hy(4)*ty + hx(4)*txy + hz(4)*tyz)*wt
c      p(12) = (hz(4)*tz + hx(4)*txz + hy(4)*tyz)*wt
c      call lku(s,u,p,nst)      
      return
c ======================================================================
c
c ... Tensoes nodais:        
c
c ......................................................................
  300 continue
      do 310 i = 1, 4
         k = (i-1)*7+1
	   p(k)   = tx(1)
	   p(k+1) = tx(2)
	   p(k+2) = tx(3)
	   p(k+3) = tx(4)
	   p(k+4) = tx(5)
	   p(k+5) = tx(6)
	   p(k+6) = tx(7)
  310 continue
c     call deform3d(hx,hy,hz,u,eps,4)
c ... deformacoes termicas: 
c      tm = (tl(1)+tl(2)+tl(3)+tl(4))/4.d0
c      eps(1) = eps(1) - alpha*tm
c      eps(2) = eps(2) - alpha*tm
c      eps(3) = eps(3) - alpha*tm
c      call stress3d(a,b,c,eps,p)
c      p(7)  = p(1)
c      p(8)  = p(2)
c      p(9)  = p(3)      
c      p(10) = p(4)      
c      p(11) = p(5)      
c      p(12) = p(6)
c      p(13) = p(1)
c      p(14) = p(2)
c      p(15) = p(3)      
c      p(16) = p(4)      
c      p(17) = p(5)      
c      p(18) = p(6)                  
c      p(19) = p(1)
c      p(20) = p(2)
c      p(21) = p(3)      
c      p(22) = p(4)      
c      p(23) = p(5)      
c      p(24) = p(6)                                          
      return      
c ======================================================================
c ......................................................................  
c
c ... Cargas distribuidas no volume e no contorno:
c            
c ......................................................................
  400 continue       
c      if(iq(1) .eq. 0) return
c ....carga aplicada na face do tetraedro- definicoes de face:
c
c     iq(1) = 1 | nós 1 - 2 - 3  |  
c             2 | nós 1 - 2 - 4  |  
c             3 | nós 2 - 3 - 4  |  
c             4 | nós 3 - 1 - 4  |  
c            
c ......................................................................
      return     
c ======================================================================
c
c ... 
c
c ......................................................................
  500 continue
      phi   = e(11)
      fc  = e(12)
      ft = e(13)
      call deform3d(hx,hy,hz,u,eps,4)
c ... deformacoes termicas: 
      tm = (tl(1)+tl(2)+tl(3)+tl(4))/4.d0
      eps(1) = eps(1) - alpha*tm
      eps(2) = eps(2) - alpha*tm
      eps(3) = eps(3) - alpha*tm
      call stress3d(a,b,c,eps,p)
      call principal_stress(p,6)
      if ( plastic(1) .eq. 0.d0)   plastic(1) = real(druckerprager3d(p,
     .                        phi,fc))
      if ( plastic(2) .eq. 0.d0)   plastic(2) = real(rankine(p,ft))
      return
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT06: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
      subroutine elmt07_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,
     .                                ma,nlit,isw,flaghidr,plastic)
c **********************************************************************
c *                                                                    *
c *   ELMT07: Elemento hexaedrico                                      *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 6)           *
c *                                                                    *
c **********************************************************************
      implicit none
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw,ma,nlit,iyied,ninc,pint(8),k,n
      integer i,j,i1,i2,i3,j1,j2,j3,druckerprager3d,rankine
      integer nen,nint,lx,ly,lz
      real*8  e(*),x(ndm,8),u(*),v(*),p(*),s(nst,nst),tl(*),ut(*)
      real*8  ym,cp,a,b,c,rn(8),sn(8),tn(8),invar2
      real*8  plastic(2)
      real*8  h(8),hx(8),hy(8),hz(8),xj(3,3),xji(3,3),dum,txl(7)
      real*8  ri,si,ti,wt,det,tm,alpha,hi(3),epsi(6),deps(6),txp(7,*)
      real*8  pg(10,10),wg(10,10),c1,c2,c3,c4,c5,c6,c7,c8
      real*8 eps(7,*),tx(7,*),devs(6),steff,theta,varj2,sint3,snphi
c      real*8  ty,tz,txy,tyz,txz
      real*8 phi,fc,ft,coes,yield,hl,c14,g,fbfc,yield2
c      real*8 plastic(2)
      data rn / 1.d0,-1.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0, 1.d0/,
     .     sn / 1.d0, 1.d0,-1.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0/,
     .     tn / 1.d0, 1.d0, 1.d0, 1.d0,-1.d0,-1.d0,-1.d0,-1.d0/
      logical flaghidr
      data nint/2/,nen/8/,pint/8,7,5,6,4,3,1,2/
c ......................................................................
      goto (100,200,300,400,200) isw
c ======================================================================
c
c.... Leitura das propriedades do material:
c
c ......................................................................
  100 continue
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c
c ... Matriz constitutiva:
c ........................................................
      ym    = e(1)
      cp    = e(2)
      alpha = e(4)
      iyied =int(e(5))
      yield = e(11)
      phi   = e(12)
      coes  = e(13)
      hl    = e(14)
      call radian(phi)
      call check_var_hidr(a,b,c,yield,phi,coes,dum,ym,cp,alpha,hi(2),
     .                           ma,ndm,flaghidr,.false.)
c      if ( ma .eq. 2)  ym = ym*hi(2)
      if( iyied .eq. 3)   yield = coes*dcos(phi)
      if( iyied .eq. 4)   yield = 6.d0*coes*dcos(phi)/(dsqrt(3.d0)*
     .    (3.d0-dsin(phi)))
      if( iyied .eq. 6 .or. iyied .eq. 7)   then
         fbfc = e(12)/e(13)
         yield2 = e(13)
      endif
c      if (flaghidr .eqv. .true.) then
c         call updhidr(ym,hi(2),ma,1)
c         call updhidr(cp,hi(2),ma,2)
c         call updhidr(alpha,hi(2),ma,4)
c      endif
c      a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c      b = cp/(1.d0-cp)
c      c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c      if (isw .eq. 3) goto 300
      if(isw .eq. 5) goto 500
c ........................................................
c ... Matriz de rigidez:
c ........................................................
      do 210 i = 1, nst
      do 210 j = 1, nst
         s(i,j) = 0.d0
  210 continue
c ........................................................
      do 260 lz = 1, nint
      ti = pg(lz,nint)
      do 250 ly = 1, nint
      si = pg(ly,nint)
      do 240 lx = 1, nint
      ri = pg(lx,nint)
      call sfhexa8 (h,hx,hy,hz,ri,si,ti,.true.,.true.)
c .........
      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,alpha,
     .                        8,ndm,ma,flaghidr,.false.)
c      if (flaghidr .neqv. .true.) then
c         tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
c     .     h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
c         call updprop(ym,tm,ma,1)
c         call updprop(cp,tm,ma,2)
c         a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c         b = cp/(1.d0-cp)
c         c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c       endif
c .........      
      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
      wt = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det*a
c ..............................................................
      do 230 i = 1, 8
      i1 = (i-1)*3+1
      i2 = i1 + 1
      i3 = i2 + 1
      do 220 j = 1, 8
         j1 = (j-1)*3+1
         j2 = j1 + 1
         j3 = j2 + 1
         s(i1,j1) = s(i1,j1) + (  hx(i)*hx(j) + c*hy(i)*hy(j) +
     .                                          c*hz(i)*hz(j))*wt
         s(i1,j2) = s(i1,j2) + (b*hx(i)*hy(j) + c*hy(i)*hx(j))*wt
         s(i1,j3) = s(i1,j3) + (b*hx(i)*hz(j) + c*hz(i)*hx(j))*wt
         s(i2,j1) = s(i2,j1) + (b*hy(i)*hx(j) + c*hx(i)*hy(j))*wt
         s(i2,j2) = s(i2,j2) + (  hy(i)*hy(j) + c*hx(i)*hx(j) +
     .                                          c*hz(i)*hz(j))*wt
         s(i2,j3) = s(i2,j3) + (b*hy(i)*hz(j) + c*hz(i)*hy(j))*wt
         s(i3,j1) = s(i3,j1) + (b*hz(i)*hx(j) + c*hx(i)*hz(j))*wt
         s(i3,j2) = s(i3,j2) + (b*hz(i)*hy(j) + c*hy(i)*hz(j))*wt
         s(i3,j3) = s(i3,j3) + (  hz(i)*hz(j) + c*hy(i)*hy(j) +
     .                                          c*hx(i)*hx(j))*wt
c ................................................................
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
      k = 0
c ....................................................................
c ... Forcas internas:
      do 290 lz = 1, nint
         ti = pg(lz,nint)
         do 280 ly = 1, nint
            si = pg(ly,nint)
            do 270 lx = 1, nint
               ri = pg(lx,nint)
               call sfhexa8 (h,hx,hy,hz,ri,si,ti,.true.,.true.)
c .........
               call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,
     .                     alpha,8,ndm,ma,flaghidr,.false.)
c               if (flaghidr .neqv. .true.) then
c                tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
c     .                h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
c                  call updprop(ym,tm,ma,1)
c                  call updprop(cp,tm,ma,2)
c                  call updprop(alpha,tm,ma,4)
c                  a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c                  b = cp/(1.d0-cp)
c                  c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c               endif
c .........   
               call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
               wt = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....     Deformacoes
               call deform3d(hx,hy,hz,u,epsi,8)
               tm = h(1)*tl(1)+h(2)*tl(2)+h(3)*tl(3)+h(4)*tl(4)+
     .         h(5)*tl(5)+h(6)*tl(6)+h(7)* tl(7)+h(8)*tl(8)
            if (nlit .eq. 1)  then
               k = k + 1
c .....     Deformacoes Termicas:
               epsi(1) = epsi(1) - alpha*tm
               epsi(2) = epsi(2) - alpha*tm
               epsi(3) = epsi(3) - alpha*tm
c               do i = 1, 8
c                 do j = 1, 6
c                 eps(j,i) = 0.d0
c                 tx(j,i) = 0.d0
c                 enddo
c                 tx(7,i) = 0.d0
c               enddo
c .....     Tensoes
               txl(1) = (epsi(1) + b*epsi(2)+b*epsi(3))*a - txp(1,k)
               txl(2) = (epsi(2) + b*epsi(1)+b*epsi(3))*a - txp(2,k)
               txl(3) = (epsi(3) + b*epsi(1)+b*epsi(2))*a - txp(3,k)
               txl(4) = a*c*epsi(4) - txp(4,k)
               txl(5) = a*c*epsi(5) - txp(5,k)
               txl(6) = a*c*epsi(6) - txp(6,k)
            else
c               if ( ma .eq. 2) then
c                  continue
c               endif
               k = k + 1
               deps(1) = epsi(1) - eps(1,k)- alpha*tm
               deps(2) = epsi(2) - eps(2,k)- alpha*tm
               deps(3) = epsi(3) - eps(3,k)- alpha*tm
               deps(4) = epsi(4) - eps(4,k)
               deps(5) = epsi(5) - eps(5,k)
               deps(6) = epsi(6) - eps(6,k)
               do i = 1, 6
                  eps(i,k) = eps(i,k) + deps(i)
               enddo
               do i = 1, 7
                  txl(i) = tx(i,k)
               enddo
c               
               c14 = cp/(1.d0-2.d0*cp)
               g = ym/(2.d0*(1.d0+cp))
c               call plasticst3d(a,b,c,deps,eps(1,k),txl,txp(1,k),yield,
c     .         hl,phi,cp,iyied,c14,g)
              call plastic3d_new(a,b,c,deps,eps(1,k),txl,txp(1,k),
     .         yield,yield2,fbfc,hl,0.d0,phi,cp,iyied,c14,g)
               do i = 1, 7
                  tx(i,k) = txl(i)
               enddo
            endif
c .....     Forcas Internas:
            do 265 i = 1, 8
                  i1 = (i-1)*3+1
                  i2 = i1 + 1
                  i3 = i2 + 1
                  p(i1) = p(i1)+(hx(i)*txl(1)+hy(i)*txl(4)+hz(i)*txl(6))
     .                    *wt
                  p(i2) = p(i2)+(hy(i)*txl(2)+hx(i)*txl(4)+hz(i)*txl(5))
     .                    *wt
                  p(i3) = p(i3)+(hz(i)*txl(3)+hx(i)*txl(6)+hy(i)*txl(5))
     .                    *wt
  265       continue
c               p(1)  = p(1) + (hx(1)*tx + hy(1)*txy + hz(1)*txz)*wt
c               p(2)  = p(2) + (hy(1)*ty + hx(1)*txy + hz(1)*tyz)*wt
c               p(3)  = p(3) + (hz(1)*tz + hx(1)*txz + hy(1)*tyz)*wt
c               p(4)  = p(4) + (hx(2)*tx + hy(2)*txy + hz(2)*txz)*wt
c               p(5)  = p(5) + (hy(2)*ty + hx(2)*txy + hz(2)*tyz)*wt
c               p(6)  = p(6) + (hz(2)*tz + hx(2)*txz + hy(2)*tyz)*wt
c               p(7)  = p(7) + (hx(3)*tx + hy(3)*txy + hz(3)*txz)*wt
c               p(8)  = p(8) + (hy(3)*ty + hx(3)*txy + hz(3)*tyz)*wt
c               p(9)  = p(9) + (hz(3)*tz + hx(3)*txz + hy(3)*tyz)*wt
c               p(10) = p(10) + (hx(4)*tx + hy(4)*txy + hz(4)*txz)*wt
c               p(11) = p(11) + (hy(4)*ty + hx(4)*txy + hz(4)*tyz)*wt
c               p(12) = p(12) + (hz(4)*tz + hx(4)*txz + hy(4)*tyz)*wt
c               p(13) = p(13) + (hx(5)*tx + hy(5)*txy + hz(5)*txz)*wt
c               p(14) = p(14) + (hy(5)*ty + hx(5)*txy + hz(5)*tyz)*wt
c               p(15) = p(15) + (hz(5)*tz + hx(5)*txz + hy(5)*tyz)*wt
c               p(16) = p(16) + (hx(6)*tx + hy(6)*txy + hz(6)*txz)*wt
c               p(17) = p(17) + (hy(6)*ty + hx(6)*txy + hz(6)*tyz)*wt
c               p(18) = p(18) + (hz(6)*tz + hx(6)*txz + hy(6)*tyz)*wt
c               p(19) = p(19) + (hx(7)*tx + hy(7)*txy + hz(7)*txz)*wt
c               p(20) = p(20) + (hy(7)*ty + hx(7)*txy + hz(7)*tyz)*wt
c               p(21) = p(21) + (hz(7)*tz + hx(7)*txz + hy(7)*tyz)*wt
c               p(22) = p(22) + (hx(8)*tx + hy(8)*txy + hz(8)*txz)*wt
c               p(23) = p(23) + (hy(8)*ty + hx(8)*txy + hz(8)*tyz)*wt
c               p(24) = p(24) + (hz(8)*tz + hx(8)*txz + hy(8)*tyz)*wt
 270        continue
 280     continue
 290  continue
c          call lku(s,u,p,24)
      return  
c ======================================================================
c
c ... Tensoes e forcas internas:
c
c ......................................................................
 300  continue
c  400 continue
c      k = 1
c      do 401 i = 1, 8
c	  ri = rn(i)
c	  si = sn(i)
c	  ti = tn(i)
c        call sfhexa8(h,hx,hy,hz,ri,si,ti,.false.,.true.)
c        call jacob3d(hx,hy,hz,xj,xji,x,det,nen,nel)
c        call deform3d(hx,hy,hz,u,eps,nen)
c        call stress3d(a,b,c,eps,p(k))
c        k = k + 6
c  401 continue
c      return
c      iyied =int(e(5))
c      yield = e(11)
c      phi   = e(12)
c      call radian(phi)
c      snphi = dsin(phi)
c      b = yield*yield
      do 320 j = 1, 8
c ......................................................................
c
c...  Transforma as tensoes nos pontos de integracao em tensoes nodais:
c ......................................................................
	   k = (j-1)*7+1
         do 310 i = 1, 6
            c1 = 0.12500d0*( tx(i,1)+tx(i,2)+tx(i,3)+tx(i,4)+tx(i,5)+
     .                       tx(i,6)+tx(i,7)+tx(i,8))
            c2 = 0.21651d0*(-tx(i,1)+tx(i,2)-tx(i,3)+tx(i,4)-tx(i,5)+
     .                       tx(i,6)-tx(i,7)+tx(i,8))
            c3 = 0.21651d0*(-tx(i,1)-tx(i,2)+tx(i,3)+tx(i,4)-tx(i,5)-
     .                       tx(i,6)+tx(i,7)+tx(i,8))
            c4 = 0.21651d0*(-tx(i,1)-tx(i,2)-tx(i,3)-tx(i,4)+tx(i,5)+
     .                       tx(i,6)+tx(i,7)+tx(i,8))
            c5 = 0.37500d0*( tx(i,1)-tx(i,2)-tx(i,3)+tx(i,4)+tx(i,5)-
     .                       tx(i,6)-tx(i,7)+tx(i,8))
            c6 = 0.37500d0*( tx(i,1)+tx(i,2)-tx(i,3)-tx(i,4)-tx(i,5)-
     .                       tx(i,6)+tx(i,7)+tx(i,8))
            c7 = 0.37500d0*( tx(i,1)-tx(i,2)+tx(i,3)-tx(i,4)-tx(i,5)+
     .                       tx(i,6)-tx(i,7)+tx(i,8))
            c8 = 0.64952d0*(-tx(i,1)+tx(i,2)+tx(i,3)-tx(i,4)+tx(i,5)-
     .                       tx(i,6)-tx(i,7)+tx(i,8))
            p(k) = c1 + c2*rn(j) + c3*sn(j) + c4*tn(j) + c5*rn(j)*sn(j)
     .          + c6*sn(j)*tn(j) + c7*rn(j)*tn(j)+ c8*rn(j)*sn(j)*tn(j)
	      k = k + 1
  310    continue
         if ( iyied .eq. 0) goto 320
c ......................................................................
c
c...  Verifica se as ensoes nos pontos de integracao em tensoes nodais:
c ......................................................................
c	   k = (j-1)*7+1
c	   a = 3.d0*invar2(p(k),6)
c         call effst3d(a,iyied,p(k),devs,steff,theta,varj2,sint3,snphi) 
c	   if(a .gt. b) then
c            j1 = pint(j)
c       call effst3d(c,iyied,tx(1,j1),devs,steff,theta,varj2,sint3,snphi)
cc	      c = 3.d0*invar2(tx(1,j1),6)
c            if(c .ge. b) then
c               p(k  ) = tx(1,j1)
c	         p(k+1) = tx(2,j1)
c	         p(k+2) = tx(3,j1)
c	         p(k+3) = tx(4,j1)
c	         p(k+4) = tx(5,j1)
c	         p(k+5) = tx(6,j1)
c	      else
c               c = dsqrt(b/a)
c               p(k  ) = p(k  )*c
c	         p(k+1) = p(k+1)*c
c	         p(k+2) = p(k+2)*c
c	         p(k+3) = p(k+3)*c
c	         p(k+4) = p(k+4)*c
c	         p(k+5) = p(k+5)*c
c            endif
c	    endif
  320 continue
      return
c      call sfhexa8 (h,hx,hy,hz,1.d0,1.d0,1.d0,.true.,.true.)
c .........
c      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,alpha,
c     .                        8,ndm,ma,flaghidr,.false.)
c               if (flaghidr .neqv. .true.) then
c                tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
c     .                h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
c                  call updprop(ym,tm,ma,1)
c                  call updprop(cp,tm,ma,2)
c                  call updprop(alpha,tm,ma,4)
c                  a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c                  b = cp/(1.d0-cp)
c                  c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c               endif
c .........  
c      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
c      call deform3d(hx,hy,hz,u,eps,8)
c      eps(1) = eps(1) - alpha*tl(1)
c      eps(2) = eps(2) - alpha*tl(1)
c      eps(3) = eps(3) - alpha*tl(1)
c      call stress3d(a,b,c,eps,p(1))
c      call sfhexa8 (h,hx,hy,hz,-1.d0,1.d0,1.d0,.true.,.true.)
cc .........
c      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,alpha,
c     .                        8,ndm,ma,flaghidr,.false.)
cc               if (flaghidr .neqv. .true.) then
cc                tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
cc     .                h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
cc                  call updprop(ym,tm,ma,1)
cc                  call updprop(cp,tm,ma,2)
c                  call updprop(alpha,tm,ma,4)
c                  a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c                  b = cp/(1.d0-cp)
c                  c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c               endif
c .........  
c      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
c      call deform3d(hx,hy,hz,u,eps,8)
c      eps(1) = eps(1) - alpha*tl(2)
c      eps(2) = eps(2) - alpha*tl(2)
c      eps(3) = eps(3) - alpha*tl(2)
c      call stress3d(a,b,c,eps,p(7))
c      call sfhexa8 (h,hx,hy,hz,-1.d0,-1.d0,1.d0,.true.,.true.)
cc .........
c      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,alpha,
c     .                        8,ndm,ma,flaghidr,.false.)
cc               if (flaghidr .neqv. .true.) then
cc                tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
cc     .                h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
cc                  call updprop(ym,tm,ma,1)
cc                  call updprop(cp,tm,ma,2)
c                  call updprop(alpha,tm,ma,4)
c                  a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c                  b = cp/(1.d0-cp)
c                  c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c               endif
c .........  
c      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
c      call deform3d(hx,hy,hz,u,eps,8)
c      eps(1) = eps(1) - alpha*tl(3)
c      eps(2) = eps(2) - alpha*tl(3)
c      eps(3) = eps(3) - alpha*tl(3)
c      call stress3d(a,b,c,eps,p(13))
c      call sfhexa8 (h,hx,hy,hz,1.d0,-1.d0,1.d0,.true.,.true.)
cc .........
c      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,alpha,
c     .                        8,ndm,ma,flaghidr,.false.)
cc               if (flaghidr .neqv. .true.) then
cc                tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
cc     .                h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
cc                  call updprop(ym,tm,ma,1)
cc                  call updprop(cp,tm,ma,2)
c                  call updprop(alpha,tm,ma,4)
c                  a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c                  b = cp/(1.d0-cp)
c                  c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c               endif
c .........  
c      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
c      call deform3d(hx,hy,hz,u,eps,8)
c      eps(1) = eps(1) - alpha*tl(4)
c      eps(2) = eps(2) - alpha*tl(4)
c      eps(3) = eps(3) - alpha*tl(4)
c      call stress3d(a,b,c,eps,p(19))
c      call sfhexa8 (h,hx,hy,hz,1.d0,1.d0,-1.d0,.true.,.true.)
cc .........
c      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,alpha,
c     .                        8,ndm,ma,flaghidr,.false.)
cc               if (flaghidr .neqv. .true.) then
cc                tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
cc     .                h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
cc                  call updprop(ym,tm,ma,1)
cc                  call updprop(cp,tm,ma,2)
c                  call updprop(alpha,tm,ma,4)
c                  a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c                  b = cp/(1.d0-cp)
c                  c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c               endif
c .........  
c      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
c      call deform3d(hx,hy,hz,u,eps,8)
c      eps(1) = eps(1) - alpha*tl(5)
c      eps(2) = eps(2) - alpha*tl(5)
c      eps(3) = eps(3) - alpha*tl(5)
c      call stress3d(a,b,c,eps,p(25))
c      call sfhexa8 (h,hx,hy,hz,-1.d0,1.d0,-1.d0,.true.,.true.)
cc .........
c      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,alpha,
c     .                        8,ndm,ma,flaghidr,.false.)
cc               if (flaghidr .neqv. .true.) then
cc                tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
cc     .                h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
cc                  call updprop(ym,tm,ma,1)
cc                  call updprop(cp,tm,ma,2)
c                  call updprop(alpha,tm,ma,4)
c                  a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c                  b = cp/(1.d0-cp)
c                  c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c               endif
c .........  
c      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
c      call deform3d(hx,hy,hz,u,eps,8)
c      eps(1) = eps(1) - alpha*tl(6)
c      eps(2) = eps(2) - alpha*tl(6)
c      eps(3) = eps(3) - alpha*tl(6)
c      call stress3d(a,b,c,eps,p(31))
c      call sfhexa8 (h,hx,hy,hz,-1.d0,-1.d0,-1.d0,.true.,.true.)
cc .........
c      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,alpha,
c     .                        8,ndm,ma,flaghidr,.false.)
cc               if (flaghidr .neqv. .true.) then
cc                tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
cc     .                h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
cc                  call updprop(ym,tm,ma,1)
cc                  call updprop(cp,tm,ma,2)
c                  call updprop(alpha,tm,ma,4)
c                  a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
c                  b = cp/(1.d0-cp)
c                  c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
c               endif
c .........  
c      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
c      call deform3d(hx,hy,hz,u,eps,8)
c      eps(1) = eps(1) - alpha*tl(7)
c      eps(2) = eps(2) - alpha*tl(7)
c      eps(3) = eps(3) - alpha*tl(7)
c      call stress3d(a,b,c,eps,p(37))
c      call sfhexa8 (h,hx,hy,hz,1.d0,-1.d0,-1.d0,.true.,.true.)
cc .........
c      call check_var_term(a,b,c,dum,dum,dum,dum,h,ut,ym,cp,alpha,
c     .                        8,ndm,ma,flaghidr,.false.)
cc               if (flaghidr .neqv. .true.) then
cc                tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
cc     .                h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
cc                  call updprop(ym,tm,ma,1)
cc                  call updprop(cp,tm,ma,2)
c                 call updprop(alpha,tm,ma,4)
cc                  a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
cc                  b = cp/(1.d0-cp)
cc                  c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
cc               endif
cc .........  
c      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)    
c      call deform3d(hx,hy,hz,u,eps,8)
c      eps(1) = eps(1) - alpha*tl(8)
c      eps(2) = eps(2) - alpha*tl(8)
c      eps(3) = eps(3) - alpha*tl(8)
c      call stress3d(a,b,c,eps,p(43)) 
c      return
cc ======================================================================
cc
cc ... Forcas internas:
c
c ......................................................................
  400 continue
c call lku(s,u,p,nst)
c     Peso proprio:
c     1 ponto de integracao:
      if(e(15) .ne. 0)   then
         call sfhexa8(h,hx,hy,hz,0.d0,0.d0,0.d0,.true.,.true.)
         call jacob3d(hx,hy,hz,xj,xji,x,det,8,nel)
         wt  = 8.d0*det*e(15)/8.d0  
         p(3) = wt
         p(6) = p(3)
         p(9) = p(3)
         p(12)= p(3)
         p(15)= p(3)
         p(18)= p(3)
         p(21)= p(3)
         p(24)= p(3)
      endif
      return
c ======================================================================
c
c
c ......................................................................
  500 continue
c ... Criterios de ruptura:
c
c ......................................................................
c ......................................................................
      i = 7
      plastic(1) =  0.12500d0*( tx(i,1)+tx(i,2)+tx(i,3)+tx(i,4)+tx(i,5)+
     .                       tx(i,6)+tx(i,7)+tx(i,8))
      if(plastic(2) .ne. 1.d0)   plastic(2) = 0.d0
      if( plastic(1) .ge. yield)   plastic(2) = 1.d0
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      print*,'Matriz geometrica nao disponivel p/ o elemento tipo 7 !'
      stop
c ....................................................................
      end
      subroutine elmt07quad_m(e,x,u,v,p,s,tl,ut,ndm,nst,nel,ma,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT07: Elemento hexaedrico de 20 nos                            *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 6)           *
c *                                                                    *
c **********************************************************************
      implicit none
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw,ma
      integer i,j,i1,i2,i3,j1,j2,j3
      integer nen,nint,lx,ly,lz
      real*8  e(*),x(ndm,8),u(*),v(*),p(*),s(nst,nst),tl(*),ut(*)
      real*8  ym,cp,a,b,c
      real*8  h(20),hx(20),hy(20),hz(20),xj(3,3),xji(3,3)
      real*8  ri,si,ti,wt,det,tm,alpha
      real*8  pg(10,10),wg(10,10)
      real*8 eps(6),tx,ty,tz,txy,tyz,txz
      data nint/4/,nen/20/
c      data nint/2/,nen/8/
c ......................................................................
      goto (100,200,200,400) isw
c ======================================================================
c
c.... Leitura das propriedades do material:
c
c ......................................................................
  100 continue
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c
c ... Matriz constitutiva:
c ........................................................
      ym    = e(1)
      cp    = e(2)
      alpha = e(4)
      call sfhexa8 (h,hx,hy,hz,ri,si,ti,.true.,.false.)
      tm = h(1)*ut(1)+h(2)*ut(2)+h(3)*ut(3)+h(4)*ut(4)+
     .    h(5)*ut(5)+h(6)*ut(6)+h(7)*ut(7)+h(8)*ut(8)
c      call varmec1(tm,ym,ma)
c      call varmec2(tm,alpha,ma)
      a = ym*(1.d0-cp)/((1.d0+cp)*(1.d0-2.d0*cp))
      b = cp/(1.d0-cp)
      c = (1.d0-2.d0*cp)/(2.d0*(1.d0-cp))
      if (isw .eq. 3) goto 300
c ........................................................
c ... Matriz de rigidez:
c ........................................................
      do 210 i = 1, nst
         do 210 j = 1, nst
            s(i,j) = 0.d0
  210 continue
c ........................................................
      do 260 lz = 1, nint
      ti = pg(lz,nint)
      do 250 ly = 1, nint
      si = pg(ly,nint)
      do 240 lx = 1, nint
      ri = pg(lx,nint)
      call sfhexa8 (h,hx,hy,hz,ri,si,ti,.false.,.true.)
c      call sfhexa20 (h,hx,hy,hz,ri,si,ti,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,8,nel)
      wt = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det*a
c ..............................................................
      call sfhexa20 (h,hx,hy,hz,ri,si,ti,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,20,nel)      
      do 230 i = 1, 20
      i1 = (i-1)*3+1
      i2 = i1 + 1
      i3 = i2 + 1
      do 220 j = 1, 20
         j1 = (j-1)*3+1
         j2 = j1 + 1
         j3 = j2 + 1
         s(i1,j1) = s(i1,j1) + (  hx(i)*hx(j) + c*hy(i)*hy(j) +
     .                                          c*hz(i)*hz(j))*wt
c
         s(i1,j2) = s(i1,j2) + (b*hx(i)*hy(j) + c*hy(i)*hx(j))*wt
c
         s(i1,j3) = s(i1,j3) + (b*hx(i)*hz(j) + c*hz(i)*hx(j))*wt
c
         s(i2,j1) = s(i2,j1) + (b*hy(i)*hx(j) + c*hx(i)*hy(j))*wt
c
         s(i2,j2) = s(i2,j2) + (  hy(i)*hy(j) + c*hx(i)*hx(j) +
     .                                          c*hz(i)*hz(j))*wt
c
         s(i2,j3) = s(i2,j3) + (b*hy(i)*hz(j) + c*hz(i)*hy(j))*wt
c
         s(i3,j1) = s(i3,j1) + (b*hz(i)*hx(j) + c*hx(i)*hz(j))*wt
c
         s(i3,j2) = s(i3,j2) + (b*hz(i)*hy(j) + c*hy(i)*hz(j))*wt
c
         s(i3,j3) = s(i3,j3) + (  hz(i)*hz(j) + c*hy(i)*hy(j) +
     .                                          c*hx(i)*hx(j))*wt
c ................................................................
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
c ....................................................................
c ... Forcas internas:
      do 290 lz = 1, nint
         ti = pg(lz,nint)
         do 280 ly = 1, nint
            si = pg(ly,nint)
            do 270 lx = 1, nint
               ri = pg(lx,nint)
               call sfhexa8 (h,hx,hy,hz,ri,si,ti,.false.,.true.)
c               call sfhexa20 (h,hx,hy,hz,ri,si,ti,.false.,.true.)
               call jacob3d (hx,hy,hz,xj,xji,x,det,8,nel)
               wt = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....     Deformacoes
               call sfhexa20 (h,hx,hy,hz,ri,si,ti,.false.,.true.)
                call jacob3d (hx,hy,hz,xj,xji,x,det,20,nel)              
               call deform3d(hx,hy,hz,u,eps,20)
c .....     Deformacoes Termicas:
               call sfhexa8 (h,hx,hy,hz,ri,si,ti,.true.,.false.)
               tm = h(1)*tl(1)+h(2)*tl(2)+h(3)*tl(3)+h(4)*tl(4)+
     .         h(5)*tl(5)+h(6)*tl(6)+h(7)* tl(7)+h(8)*tl(8)
               eps(1) = eps(1) - alpha*tm
               eps(2) = eps(2) - alpha*tm
               eps(3) = eps(3) - alpha*tm
c .....     Tensoes
               tx = (eps(1) + b*eps(2)+b*eps(3))*a
               ty = (eps(2) + b*eps(1)+b*eps(3))*a
               tz = (eps(3) + b*eps(1)+b*eps(2))*a 
               txy = a*c*eps(4)
               tyz = a*c*eps(5)
               txz = a*c*eps(6)
c .....     Forcas Internas:
c               call sfhexa20 (h,hx,hy,hz,ri,si,ti,.false.,.true.)
               p(1)  = p(1) + (hx(1)*tx + hy(1)*txy + hz(1)*txz)*wt
               p(2)  = p(2) + (hy(1)*ty + hx(1)*txy + hz(1)*tyz)*wt
               p(3)  = p(3) + (hz(1)*tz + hx(1)*txz + hy(1)*tyz)*wt
               p(4)  = p(4) + (hx(2)*tx + hy(2)*txy + hz(2)*txz)*wt
               p(5)  = p(5) + (hy(2)*ty + hx(2)*txy + hz(2)*tyz)*wt
               p(6)  = p(6) + (hz(2)*tz + hx(2)*txz + hy(2)*tyz)*wt
               p(7)  = p(7) + (hx(3)*tx + hy(3)*txy + hz(3)*txz)*wt
               p(8)  = p(8) + (hy(3)*ty + hx(3)*txy + hz(3)*tyz)*wt
               p(9)  = p(9) + (hz(3)*tz + hx(3)*txz + hy(3)*tyz)*wt
               p(10) = p(10) + (hx(4)*tx + hy(4)*txy + hz(4)*txz)*wt
               p(11) = p(11) + (hy(4)*ty + hx(4)*txy + hz(4)*tyz)*wt
               p(12) = p(12) + (hz(4)*tz + hx(4)*txz + hy(4)*tyz)*wt
               p(13) = p(13) + (hx(5)*tx + hy(5)*txy + hz(5)*txz)*wt
               p(14) = p(14) + (hy(5)*ty + hx(5)*txy + hz(5)*tyz)*wt
               p(15) = p(15) + (hz(5)*tz + hx(5)*txz + hy(5)*tyz)*wt
               p(16) = p(16) + (hx(6)*tx + hy(6)*txy + hz(6)*txz)*wt
               p(17) = p(17) + (hy(6)*ty + hx(6)*txy + hz(6)*tyz)*wt
               p(18) = p(18) + (hz(6)*tz + hx(6)*txz + hy(6)*tyz)*wt
               p(19) = p(19) + (hx(7)*tx + hy(7)*txy + hz(7)*txz)*wt
               p(20) = p(20) + (hy(7)*ty + hx(7)*txy + hz(7)*tyz)*wt
               p(21) = p(21) + (hz(7)*tz + hx(7)*txz + hy(7)*tyz)*wt
               p(22) = p(22) + (hx(8)*tx + hy(8)*txy + hz(8)*txz)*wt
               p(23) = p(23) + (hy(8)*ty + hx(8)*txy + hz(8)*tyz)*wt
               p(24) = p(24) + (hz(8)*tz + hx(8)*txz + hy(8)*tyz)*wt
               p(25)  = p(25) + (hx(9)*tx + hy(9)*txy + hz(9)*txz)*wt
               p(26)  = p(26) + (hy(9)*ty + hx(9)*txy + hz(9)*tyz)*wt
               p(27)  = p(27) + (hz(9)*tz + hx(9)*txz + hy(9)*tyz)*wt
               p(28)  = p(28) + (hx(10)*tx + hy(10)*txy + hz(10)*txz)*wt
               p(29)  = p(29) + (hy(10)*ty + hx(10)*txy + hz(10)*tyz)*wt
               p(30)  = p(30) + (hz(10)*tz + hx(10)*txz + hy(10)*tyz)*wt
               p(31)  = p(31) + (hx(11)*tx + hy(11)*txy + hz(11)*txz)*wt
               p(32)  = p(32) + (hy(11)*ty + hx(11)*txy + hz(11)*tyz)*wt
               p(33)  = p(33) + (hz(11)*tz + hx(11)*txz + hy(11)*tyz)*wt
               p(34)  = p(34) + (hx(12)*tx + hy(12)*txy + hz(12)*txz)*wt
               p(35)  = p(35) + (hy(12)*ty + hx(12)*txy + hz(12)*tyz)*wt
               p(36)  = p(36) + (hz(12)*tz + hx(12)*txz + hy(12)*tyz)*wt
               p(37)  = p(37) + (hx(13)*tx + hy(13)*txy + hz(13)*txz)*wt
               p(38)  = p(38) + (hy(13)*ty + hx(13)*txy + hz(13)*tyz)*wt
               p(39)  = p(39) + (hz(13)*tz + hx(13)*txz + hy(13)*tyz)*wt
               p(40)  = p(40) + (hx(14)*tx + hy(14)*txy + hz(14)*txz)*wt
               p(41)  = p(41) + (hy(14)*ty + hx(14)*txy + hz(14)*tyz)*wt
               p(42)  = p(42) + (hz(14)*tz + hx(14)*txz + hy(14)*tyz)*wt
               p(43)  = p(43) + (hx(15)*tx + hy(15)*txy + hz(15)*txz)*wt
               p(44)  = p(44) + (hy(15)*ty + hx(15)*txy + hz(15)*tyz)*wt
               p(45)  = p(45) + (hz(15)*tz + hx(15)*txz + hy(15)*tyz)*wt
               p(46)  = p(46) + (hx(16)*tx + hy(16)*txy + hz(16)*txz)*wt
               p(47)  = p(47) + (hy(16)*ty + hx(16)*txy + hz(16)*tyz)*wt
               p(48)  = p(48) + (hz(16)*tz + hx(16)*txz + hy(16)*tyz)*wt
               p(49)  = p(49) + (hx(17)*tx + hy(17)*txy + hz(17)*txz)*wt
               p(50)  = p(50) + (hy(17)*ty + hx(17)*txy + hz(17)*tyz)*wt
               p(51)  = p(51) + (hz(17)*tz + hx(17)*txz + hy(17)*tyz)*wt
               p(52)  = p(52) + (hx(18)*tx + hy(18)*txy + hz(18)*txz)*wt
               p(53)  = p(53) + (hy(18)*ty + hx(18)*txy + hz(18)*tyz)*wt
               p(54)  = p(54) + (hz(18)*tz + hx(18)*txz + hy(18)*tyz)*wt
               p(55)  = p(55) + (hx(19)*tx + hy(19)*txy + hz(19)*txz)*wt
               p(56)  = p(56) + (hy(19)*ty + hx(19)*txy + hz(19)*tyz)*wt
               p(57)  = p(57) + (hz(19)*tz + hx(19)*txz + hy(19)*tyz)*wt
               p(58)  = p(58) + (hx(20)*tx + hy(20)*txy + hz(20)*txz)*wt
               p(59)  = p(59) + (hy(20)*ty + hx(20)*txy + hz(20)*tyz)*wt
               p(60)  = p(60) + (hz(20)*tz + hx(20)*txz + hy(20)*tyz)*wt
 270        continue
 280     continue
 290  continue
c          call lku(s,u,p,nst)
      return  
c ======================================================================
c
c ... Tensoes e forcas internas:
c
c ......................................................................
 300  continue
      call sfhexa8 (h,hx,hy,hz,1.d0,1.d0,1.d0,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
      call deform3d(hx,hy,hz,u,eps,8)
      eps(1) = eps(1) - e(4)*tl(1)
      eps(2) = eps(2) - e(4)*tl(1)
      eps(3) = eps(3) - e(4)*tl(1)
      call stress3d(a,b,c,eps,p(1))
      call sfhexa8 (h,hx,hy,hz,-1.d0,1.d0,1.d0,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
      call deform3d(hx,hy,hz,u,eps,8)
      eps(1) = eps(1) - e(4)*tl(2)
      eps(2) = eps(2) - e(4)*tl(2)
      eps(3) = eps(3) - e(4)*tl(2)
      call stress3d(a,b,c,eps,p(7))
      call sfhexa8 (h,hx,hy,hz,-1.d0,-1.d0,1.d0,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
      call deform3d(hx,hy,hz,u,eps,8)
      eps(1) = eps(1) - e(4)*tl(3)
      eps(2) = eps(2) - e(4)*tl(3)
      eps(3) = eps(3) - e(4)*tl(3)
      call stress3d(a,b,c,eps,p(13))
      call sfhexa8 (h,hx,hy,hz,1.d0,-1.d0,1.d0,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
      call deform3d(hx,hy,hz,u,eps,8)
      eps(1) = eps(1) - e(4)*tl(4)
      eps(2) = eps(2) - e(4)*tl(4)
      eps(3) = eps(3) - e(4)*tl(4)
      call stress3d(a,b,c,eps,p(19))
      call sfhexa8 (h,hx,hy,hz,1.d0,1.d0,-1.d0,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
      call deform3d(hx,hy,hz,u,eps,8)
      eps(1) = eps(1) - e(4)*tl(5)
      eps(2) = eps(2) - e(4)*tl(5)
      eps(3) = eps(3) - e(4)*tl(5)
      call stress3d(a,b,c,eps,p(25))
      call sfhexa8 (h,hx,hy,hz,-1.d0,1.d0,-1.d0,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
      call deform3d(hx,hy,hz,u,eps,8)
      eps(1) = eps(1) - e(4)*tl(6)
      eps(2) = eps(2) - e(4)*tl(6)
      eps(3) = eps(3) - e(4)*tl(6)
      call stress3d(a,b,c,eps,p(31))
      call sfhexa8 (h,hx,hy,hz,-1.d0,-1.d0,-1.d0,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)
      call deform3d(hx,hy,hz,u,eps,8)
      eps(1) = eps(1) - e(4)*tl(7)
      eps(2) = eps(2) - e(4)*tl(7)
      eps(3) = eps(3) - e(4)*tl(7)
      call stress3d(a,b,c,eps,p(37))
      call sfhexa8 (h,hx,hy,hz,1.d0,-1.d0,-1.d0,.false.,.true.)
      call jacob3d (hx,hy,hz,xj,xji,x,det,nen,nel)    
      call deform3d(hx,hy,hz,u,eps,8)
      eps(1) = eps(1) - e(4)*tl(8)
      eps(2) = eps(2) - e(4)*tl(8)
      eps(3) = eps(3) - e(4)*tl(8)
      call stress3d(a,b,c,eps,p(43)) 
      return
c ======================================================================
c
c ... Forcas internas:
c
c ......................................................................
  400 continue
c call lku(s,u,p,nst)
      return
c ======================================================================
c
c ... Matriz de massa:
c
c ......................................................................
  500 continue
      print*,'Matriz de massa nao disponivel para o elemento tipo 6 !'
      stop
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      print*,'Matriz geometrica nao disponivel p/ o elemento tipo 6 !'
      stop
c ....................................................................
      end
      subroutine elmt08_m(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,ma,nlit,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT01_I: Elemento de interface de 4 nos                         *
c *   ------  Para o caso 2D                                           * 
c *                                                                    *
c *    4----------3     |                                              *
c *    |          |     h -> espessura                                 *
c *    1----------2    _|_                                             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
c      common /gauss/ pg,wg
      integer ndm,nst,nel,isw,ma,nlit
      integer nint,i,j,k,l,lx,ly,pi,NPI
      real*8  e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*)
      real*8  h(4),hx(4),hy(4),xj(2,2),xji(2,2),det,wt,w(3)
      real*8  ym,pr,thic,a,d11,d12,d22,d33,tm,ponto,peso
      real*8  N(2,8),Lc1,Lc2,auxt1,auxt2
      real*8  tx(*),aux,ZERO,d11Aux,d22Aux,r(2,2),du(8),u0(8)
      real*8 N1,N2,N3,N4,tx1,tx2,esp
      logical flagt
      parameter (ZERO = 1.d-15)
      data    nint/2/
c ......................................................................
      goto (100,200,300,400) isw
c ======================================================================
c
c.... Imput material properties:
c
  100 continue
c     e(1) = parametro de penalidade para direcao normal
c     e(2) = parametro de penalidade para direcao tangencial
c     e(3) = espessura do elemento de interface
      return
c ======================================================================
c.... Matrix K:
c
c ......................................................................
  200 continue
      flagt = .false.
c
c ... Matriz constitutiva:
c ......................................................................
      d11Aux    = e(1) !d11 parametro de penalidade para direcao tangencial
      d22Aux    = e(2) !d22 parametro de penalidade para direcao normal
      esp   = e(3)    !h
c ... Matriz de rigidez do Elemento de interface:
c ......................................................................
      do i = 1, nst
         do j = 1, nst
            s(i,j) = 0.d0
         enddo
      enddo     
c
c Pontos e pesos para esquema de integracao de Gauss - 3 pontos
c
      hx(1) = -SQRT(0.6) !GAUSS3(1) = -SQRT(0.6)
      hx(2) = 0.0
      hx(3) = SQRT(0.6)
      w(1) = 0.555555555555555556 !PESO3(1) = 0.55555555555555
      w(2) = 0.888888888888888889
      w(3) = 0.555555555555555556
      NPI=3
c
c comprimentos das linhas inferior e superior como se segue
c	  
      Lc1 = DSQRT( (x(1,2)-x(1,1))*(x(1,2)-x(1,1)) +
     .	    (x(2,2)-x(2,1))*(x(2,2)-x(2,1)) )
      det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
c ...........................................................
      Lc2 = DSQRT( (x(1,4)-x(1,3))*(x(1,4)-x(1,3)) +
     .       (x(2,4)-x(2,3))*(x(2,4)-x(2,3)) )
c       1000000000.0
      thic = Lc1/1.e+09        
      d11 = ((d11Aux)/(2.0*(1.0+d22Aux)))/thic
      d22 = (d11Aux)/thic
c
c       d11 = d11Aux
c       d22 = d22Aux
      if (isw .eq. 3) goto 300
      do 240 i=1,8
         do 230 j=1,8      
         aux=0.0
         
c     Loop sobre os pontos de integracao   
            do 220 pi=1,NPI
              ponto = hx(pi)
              peso = w(pi)
c     Funcao de interpolacao 
              h(1) = 0.5*(1.0 - ponto)
              h(2) = 0.5*(1.0 + ponto)
c     Matriz das funcoes de interpolacao
          N(1,1) = -h(1)*Lc1
          N(2,1) =  0.0
          N(1,2) =  0.0
          N(2,2) = -h(1)*Lc1
          N(1,3) = -h(2)*Lc1
          N(2,3) =  0.0
          N(1,4) =  0.0
          N(2,4) = -h(2)*Lc1     
          N(1,5) =  h(2)*Lc2
          N(2,5) =  0.0
          N(1,6) =  0.0
          N(2,6) =  h(2)*Lc2
          N(1,7) =  h(1)*Lc2
          N(2,7) =  0.0
          N(1,8) =  0.0
          N(2,8) =  h(1)*Lc2
c
          wt = peso*0.5d0
          if (esp .ne. 0.d0)   wt = wt*esp          
c
c         Matriz de rigidez:
c
         aux = aux + (N(1,i)*d11*N(1,j) + N(2,i)*d22*N(2,j))*wt        
220      continue 
c
c Matriz de rigidez:               
c
         s(i,j) = aux
c       
230      continue
240   continue
      call irotmatrix(x,r)
      call irotatek(s,r)
c      call irotmatrix(x,r)
c
c ..... Forcas Internas
c
      call irotate(u,r,u,.false.) 
      call irotate(u(3),r,u(3),.false.)
      call irotate(u(5),r,u(5),.false.)
      call irotate(u(7),r,u(7),.false.)
      if (nlit .eq. 1)  then
          do i = 1, 8
c             u0(i) = 0.d0
             u(i) = 0.d0
          enddo
c          tx(1) = 0.d0
c          tx(2) = 0
c          tx1 = d22*((u(8)-u(2)+u(6)-u(4)))*0.5d0*Lc1 !y
c          tx2 = d11*((u(7)-u(1)+u(5)-u(3)))*0.5d0*Lc1 !xy
      else 
         do 241 i = 1, 8
            du(i) = u(i) - u0(i)
            u0(i) = u(i)
  241    continue
         call istress2d(d11,d22,Lc1,u,du,tx,0.d0,flagt)
         tx1 = tx(1)
         tx2 = tx(2)
cc       tx = 0.d0  !x
cc          txl(1) = d22*((u(8)-u(2)+u(6)-u(4)))*0.5d0*Lc1 !y
cc          txl(2) = d11*((u(7)-u(1)+u(5)-u(3)))*0.5d0*Lc1 !xy
      endif
      call irotate(u,r,u,.true.)
      call irotate(u(3),r,u(3),.true.)
      call irotate(u(5),r,u(5),.true.)
      call irotate(u(7),r,u(7),.true.)
      do 245 i = 1, 8
         p(i) = 0.d0
  245 continue
c ... Modifica a matriz de rigidez no caso de tracao
      if (flagt .eqv. .true.)   then
         do i = 1, nst
            do j = 1, nst
               s(i,j) = 0.d0
            enddo
         enddo 
c         d11 = ((d11Aux)/(2.0*(1.0+d22Aux)))*thic
c         d22 = (d11Aux)*thic
c         do 270 i=1,8
c            do 260 j=1,8      
c               aux=0.0
c         
c     Loop sobre os pontos de integracao   
c               do 250 pi=1,NPI
c                  ponto = hx(pi)
c                  peso = w(pi)
c     Funcao de interpolacao 
c                  h(1) = 0.5*(1.0 - ponto)
c                  h(2) = 0.5*(1.0 + ponto)
c     Matriz das funcoes de interpolacao
c                  N(1,1) = -h(1)*Lc1
c                  N(2,1) =  0.d0
c                  N(1,2) =  0.0d0
c                  N(2,2) = -h(1)*Lc1
c                  N(1,3) = -h(2)*Lc1
c                  N(2,3) =  0.0d0
c                  N(1,4) =  0.0d0
c                  N(2,4) = -h(2)*Lc1     
c                  N(1,5) =  h(2)*Lc2
c                  N(2,5) =  0.0d0
c                  N(1,6) =  0.0d0
c                  N(2,6) =  h(2)*Lc2
c                  N(1,7) =  h(1)*Lc2
c                  N(2,7) =  0.0d0
c                  N(1,8) =  0.0d0
c                  N(2,8) =  h(1)*Lc2
c
c                  wt = peso*0.5
c                  if (esp .ne. 0.d0)   wt = wt*esp 
c
c         Matriz de rigidez:
c
c                  aux = aux + (N(1,i)*d11*N(1,j) + N(2,i)*d22*N(2,j))*wt
c250            continue 
c
c Matriz de rigidez:               
c
c               s(i,j) = aux
c       
c260         continue
c270      continue
c         call irotmatrix(x,r)
c         call irotatek(s,r)
         if ( nlit .eq. 2) goto 280
      endif
       p(1) = s(1,1)*u(1)+s(1,3)*u(3)+s(1,5)*u(5)+s(1,7)*u(7)
       p(2) = s(2,2)*u(2)+s(2,4)*u(4)+s(2,6)*u(6)+s(2,8)*u(8)
       p(3) = s(3,1)*u(1)+s(3,3)*u(3)+s(3,5)*u(5)+s(3,7)*u(7)
       p(4) = s(4,2)*u(2)+s(4,4)*u(4)+s(4,6)*u(6)+s(4,8)*u(8)
       p(5) = s(5,1)*u(1)+s(5,3)*u(3)+s(5,5)*u(5)+s(5,7)*u(7)
       p(6) = s(6,2)*u(2)+s(6,4)*u(4)+s(6,6)*u(6)+s(6,8)*u(8)
       p(7) = s(7,1)*u(1)+s(7,3)*u(3)+s(7,5)*u(5)+s(7,7)*u(7)
       p(8) = s(8,2)*u(2)+s(8,4)*u(4)+s(8,6)*u(6)+s(8,8)*u(8)
       return
c
c     Loop sobre os pontos de integracao
c      do 250 pi=1,NPI
c         ponto = hx(pi)
c         peso = w(pi)
c     Funcao de interpolacao 
c         h(1) = 0.5*(1.0 - ponto)
c         h(2) = 0.5*(1.0 + ponto)
c     Matriz das funcoes de interpolacao
c        N1 = -h(1)*Lc1
c         N2 = -h(2)*Lc1
c         N3 =  h(2)*Lc2
c         N4 =  h(1)*Lc2
c         wt =  peso*0.5d0
c          p(1) = p(1) + N1*tx2*wt
c          p(2) = p(2) + N1*tx1*wt
c          p(3) = p(3) + N2*tx2*wt
c          p(4) = p(4) + N2*tx1*wt
c          p(5) = p(5) + N3*tx2*wt
c          p(6) = p(6) + N3*tx1*wt
c          p(7) = p(7) + N4*tx2*wt
c          p(8) = p(8) + N4*tx1*wt
c250   continue 
c      do  255 i = 1, 4
c          j = 2*i-1
c          call irotate(p(j),r,p(j),.true.)
c         p(1) = 0.d0
c         p(2) = -0.50d0
c         p(3) = 0.d00
c         p(4) = -0.50d0
c        p(5) = 0.d0
c         p(6) = 0.50d0
c         p(7) = 0.d0
c         p(8) = 0.50d0
280   continue
      return 
c ======================================================================
c.... Tensoes nodais e Deslocamentos Relativos:       
c ......................................................................
  300 continue
c
c .... Sistema local de coordenadas:
c
      call irotmatrix(x,r)
c      call irotate(u,r,u,.false.) 
c      call irotate(u(3),r,u(3),.false.)
c      call irotate(u(5),r,u(5),.false.)
c      call irotate(u(7),r,u(7),.false.)      
      p(1) = 0.d0
      p(2) = tx(1)
      p(3) = tx(2)
      call irotate2(p(1),r,p(1),.true.) 
      p(4) = p(3)
      p(3) = 0.d0
c       p(1) = 0.0  !x
c       p(3) = d11*((u(7)-u(1)+u(5)-u(3)))*0.5*Lc1 !xy      
c       p(2) = d22*((u(8)-u(2)+u(6)-u(4)))*0.5*Lc1 !y
c      call irotate2(p(1),r,p(1),.true.)
c      p(4) = p(3)
c      p(3) = 0.d0
      do 310 i = 2, 4
         k = (i-1)*5+1
         p(k)   = p(1)
         p(k+1) = p(2)
         p(k+2) = p(3)
         p(k+3) = p(4)
  310 continue
      return
  400 continue
      return
c ......................................................................
      end
      subroutine elmt08_grampo(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,
     .                                              ma,nlit,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT01_I: Elemento de interface de 4 nos                         *
c *   ------  Para o caso 2D                                           * 
c *                                                                    *
c *    4----------3     |                                              *
c *    |          |     h -> espessura                                 *
c *    1----------2    _|_                                             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
c      common /gauss/ pg,wg
      integer ndm,nst,nel,isw,ma,nlit
      integer nint,i,j,k,l,lx,ly,pi,NPI
      real*8  e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*)
      real*8  h(4),hx(4),hy(4),xj(2,2),xji(2,2),det,wt,w(3)
      real*8  ym,pr,thic,a,d11,d12,d22,d33,tm,ponto,peso
      real*8  N(2,8),Lc1,Lc2,auxt1,auxt2
      real*8  tx(*),aux,ZERO,d11Aux,d22Aux,r(2,2),du(8),u0(8)
      real*8 N1,N2,N3,N4,tx1,tx2,e2,sy,esp
      logical flagt
      parameter (ZERO = 1.d-15)
      data    nint/2/
c ......................................................................
      goto (100,200,300,400) isw
c ======================================================================
c
c.... Imput material properties:
c
  100 continue
c     e(1) = parametro de penalidade para direcao normal
c     e(2) = parametro de penalidade para direcao tangencial
c     e(3) = espessura do elemento de interface
      return
c ======================================================================
c.... Matrix K:
c
c ......................................................................
  200 continue
      flagt = .false.
c
c ... Matriz constitutiva:
c ......................................................................
      ym    = e(1) !d11 parametro de penalidade para direcao tangencial
      pr    = e(2) !d22 parametro de penalidade para direcao normal
      esp   = e(3)    !h
      e2    = e(4)
      sy    = e(5)
c ... Matriz de rigidez do Elemento de interface:
c ......................................................................
      do i = 1, nst
         do j = 1, nst
            s(i,j) = 0.d0
         enddo
      enddo
c
c Pontos e pesos para esquema de integracao de Gauss - 3 pontos
c
      hx(1) = -SQRT(0.6) !GAUSS3(1) = -SQRT(0.6)
      hx(2) = 0.0
      hx(3) = SQRT(0.6)
      w(1) = 0.555555555555555556 !PESO3(1) = 0.55555555555555
      w(2) = 0.888888888888888889
      w(3) = 0.555555555555555556
      NPI=3
c
c comprimentos das linhas inferior e superior como se segue
c	  
      Lc1 = DSQRT( (x(1,2)-x(1,1))*(x(1,2)-x(1,1)) +
     .           (x(2,2)-x(2,1))*(x(2,2)-x(2,1)) )
      det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
c ...........................................................
      Lc2 = DSQRT( (x(1,4)-x(1,3))*(x(1,4)-x(1,3)) +
     .       (x(2,4)-x(2,3))*(x(2,4)-x(2,3)) )
c       1000000000.0
      thic = Lc1/1.e+09        
      d11 = ((ym)/(2.0d0*(1.0d0+pr)))/thic
      d22 = (ym)/thic
c
c       d11 = d11Aux
c       d22 = d22Aux
      if (isw .eq. 3) goto 300
      do 240 i=1,8
         do 230 j=1,8      
            aux=0.0
         
c     Loop sobre os pontos de integracao   
            do 220 pi=1,NPI
              ponto = hx(pi)
              peso = w(pi)
c     Funcao de interpolacao 
              h(1) = 0.5*(1.0 - ponto)
              h(2) = 0.5*(1.0 + ponto)
c     Matriz das funcoes de interpolacao
          N(1,1) = -h(1)*Lc1
          N(2,1) =  0.0
          N(1,2) =  0.0
          N(2,2) = -h(1)*Lc1
          N(1,3) = -h(2)*Lc1
          N(2,3) =  0.0
          N(1,4) =  0.0
          N(2,4) = -h(2)*Lc1     
          N(1,5) =  h(2)*Lc2
          N(2,5) =  0.0
          N(1,6) =  0.0
          N(2,6) =  h(2)*Lc2
          N(1,7) =  h(1)*Lc2
          N(2,7) =  0.0
          N(1,8) =  0.0
          N(2,8) =  h(1)*Lc2
c
          wt = peso*0.5d0
          if (esp .ne. 0.d0)   wt = wt*esp          
c
c         Matriz de rigidez:
c
         aux = aux + (N(1,i)*d11*N(1,j) + N(2,i)*d22*N(2,j))*wt        
220      continue 
c
c Matriz de rigidez:               
c
         s(i,j) = aux
c       
230      continue
240   continue
      call irotmatrix(x,r)
      call irotatek(s,r)
c      call irotmatrix(x,r)
c
c ..... Forcas Internas
c
      call irotate(u,r,u,.false.) 
      call irotate(u(3),r,u(3),.false.)
      call irotate(u(5),r,u(5),.false.)
      call irotate(u(7),r,u(7),.false.)
      if (nlit .eq. 1)  then
          do i = 1, 8
             u0(i) = 0.d0
             u(i) = 0.d0
          enddo
c          tx(1) = 0.d0
c          tx(2) = 0.d0
c          tx(3) = 0.d0
          tx1 = d22*((u(8)-u(2)+u(6)-u(4)))*0.5d0*Lc1 !y
          tx2 = d11*((u(7)-u(1)+u(5)-u(3)))*0.5d0*Lc1 !xy
      else
         do 241 i = 1, 8
            du(i) = u(i) - u0(i)
            u0(i) = u(i)
  241    continue
c         if (nlit .eq. 2)   then
c            call istress2d(d11,d22,Lc1,u,du,tx,0.d0,flagt)
c         else 
            call istress2d_grampo(d11,d22,Lc1,u,du,tx,e2,pr,sy,thic,
     .                                                       flagt,nlit)
c            tx1 = (tx(1)+tx(4))/2.d0
c            tx2 = (tx(2)+tx(5))/2.d0
            tx1=tx(1)
            tx2=tx(2)
c         endif
      endif
      call irotate(u,r,u,.true.)
      call irotate(u(3),r,u(3),.true.)
      call irotate(u(5),r,u(5),.true.)
      call irotate(u(7),r,u(7),.true.)
      do 245 i = 1, 8
         p(i) = 0.d0
  245 continue
c ... Modifica a matriz de rigidez no caso de tracao
      if (flagt)   then
         do i = 1, nst
            do j = 1, nst
               s(i,j) = 0.d0
            enddo
         enddo 
         d11 = e2*(1.d0-pr)/((1.d0+pr)*(1.d0-2.d0*pr))
         d22 = e2/(2.d0*(1.d0+pr))
c         d11 = ((ym)/(2.0*(1.0+pr)))*thic
c         d22 = (ym)*thic
         do 270 i=1,8
            do 260 j=1,8      
               aux=0.0
         
c     Loop sobre os pontos de integracao   
               do 250 pi=1,NPI
                  ponto = hx(pi)
                  peso = w(pi)
c     Funcao de interpolacao 
                  h(1) = 0.5*(1.0 - ponto)
                  h(2) = 0.5*(1.0 + ponto)
c     Matriz das funcoes de interpolacao
                  N(1,1) = -h(1)*Lc1
                  N(2,1) =  0.0
                  N(1,2) =  0.0
                  N(2,2) = -h(1)*Lc1
                  N(1,3) = -h(2)*Lc1
                  N(2,3) =  0.0
                  N(1,4) =  0.0
                  N(2,4) = -h(2)*Lc1     
                  N(1,5) =  h(2)*Lc2
                  N(2,5) =  0.0
                  N(1,6) =  0.0
                  N(2,6) =  h(2)*Lc2
                  N(1,7) =  h(1)*Lc2
                  N(2,7) =  0.0
                  N(1,8) =  0.0
                  N(2,8) =  h(1)*Lc2

                  wt = peso*0.5
                  if (esp .ne. 0.d0)   wt = wt*esp 
c
c         Matriz de rigidez:
c
                  aux = aux + (N(1,i)*d11*N(1,j) + N(2,i)*d22*N(2,j))*wt
250            continue 
c
c Matriz de rigidez:               
c
               s(i,j) = aux
c       
260         continue
270      continue
         call irotmatrix(x,r)
         call irotatek(s,r)
         if ( nlit .eq. 2) goto 280
c         p(1) = -tx2*Lc1*0.5d0
c         p(2) = -tx1*Lc1*0.5d0
c         p(3) = -tx2*Lc1*0.5d0
c         p(4) = -tx1*Lc1*0.5d0
c         p(5) =  tx2*Lc1*0.5d0
c         p(6) =  tx1*Lc1*0.5d0
c         p(7) =  tx2*Lc1*0.5d0
c         p(8) =  tx1*Lc1*0.5d0
c         do  255 i = 1, 4
c            j = 2*i-1
c            call irotate(p(j),r,p(j),.true.)
c255      continue
      else
         if (nlit .ge. 2) call istress2d(d11,d22,Lc1,u,du,tx,0.d0,flagt)
      endif
c         p(1) = -tx2*Lc1*0.5d0
c         p(2) = -tx1*Lc1*0.5d0
c         p(3) = -tx2*Lc1*0.5d0
c         p(4) = -tx1*Lc1*0.5d0
c         p(5) =  tx2*Lc1*0.5d0
c         p(6) =  tx1*Lc1*0.5d0
c         p(7) =  tx2*Lc1*0.5d0
c         p(8) =  tx1*Lc1*0.5d0
       p(1) = p(1)+s(1,1)*u(1)+s(1,3)*u(3)+s(1,5)*u(5)+s(1,7)*u(7)
       p(2) = p(2)+s(2,2)*u(2)+s(2,4)*u(4)+s(2,6)*u(6)+s(2,8)*u(8)
       p(3) = p(3)+s(3,1)*u(1)+s(3,3)*u(3)+s(3,5)*u(5)+s(3,7)*u(7)
       p(4) = p(4)+s(4,2)*u(2)+s(4,4)*u(4)+s(4,6)*u(6)+s(4,8)*u(8)
       p(5) = p(5)+s(5,1)*u(1)+s(5,3)*u(3)+s(5,5)*u(5)+s(5,7)*u(7)
       p(6) = p(6)+s(6,2)*u(2)+s(6,4)*u(4)+s(6,6)*u(6)+s(6,8)*u(8)
       p(7) = p(7)+s(7,1)*u(1)+s(7,3)*u(3)+s(7,5)*u(5)+s(7,7)*u(7)
       p(8) = p(8)+s(8,2)*u(2)+s(8,4)*u(4)+s(8,6)*u(6)+s(8,8)*u(8)
c     Loop sobre os pontos de integracao
c      do 250 pi=1,NPI
c         ponto = hx(pi)
c         peso = w(pi)
c     Funcao de interpolacao 
c         h(1) = 0.5*(1.0 - ponto)
c         h(2) = 0.5*(1.0 + ponto)
c     Matriz das funcoes de interpolacao
c        N1 = -h(1)*Lc1
c         N2 = -h(2)*Lc1
c         N3 =  h(2)*Lc2
c         N4 =  h(1)*Lc2
c         wt =  peso*0.5d0
c          p(1) = p(1) + N1*tx2*wt
c          p(2) = p(2) + N1*tx1*wt
c          p(3) = p(3) + N2*tx2*wt
c          p(4) = p(4) + N2*tx1*wt
c          p(5) = p(5) + N3*tx2*wt
c          p(6) = p(6) + N3*tx1*wt
c          p(7) = p(7) + N4*tx2*wt
c          p(8) = p(8) + N4*tx1*wt
c250   continue 
c      do  255 i = 1, 4
c          j = 2*i-1
c          call irotate(p(j),r,p(j),.true.)
c       if (nlit .eq. 1) then
c         p(1) = 0.d0
c         p(2) = -50d0
c         p(3) = 0.d00
c         p(4) = -50d0
c         p(5) = 0.d0
c         p(6) = 50d0
c         p(7) = 0.d0
c         p(8) = 50d0
c       endif
280   continue
      return 
c ======================================================================
c.... Tensoes nodais e Deslocamentos Relativos:       
c ......................................................................
  300 continue
c
c .... Sistema local de coordenadas:
c
      call irotmatrix(x,r)
c      call irotate(u,r,u,.false.) 
c      call irotate(u(3),r,u(3),.false.)
c      call irotate(u(5),r,u(5),.false.)
c      call irotate(u(7),r,u(7),.false.)   
c          e2 = 0.5d0
cc          pr = 0.3d0
c         d22 = e2*(1.d0-pr)/((1.d0+pr)*(1.d0-2.d0*pr))
c         ym = ((u(8)-u(2)+u(6)-u(4)))*0.5d0
      p(1) = 0.d0
      p(2) = tx(1)
      p(3) = tx(2)
      call irotate2(p(1),r,p(1),.true.) 
      p(4) = p(3)
      p(3) = 0.d0
c       p(1) = 0.0  !x
c       p(3) = d11*((u(7)-u(1)+u(5)-u(3)))*0.5*Lc1 !xy      
c       p(2) = d22*((u(8)-u(2)+u(6)-u(4)))*0.5*Lc1 !y
c      call irotate2(p(1),r,p(1),.true.)
c      p(4) = p(3)
c      p(3) = 0.d0
      do 310 i = 2, 4
         k = (i-1)*5+1
         p(k)   = p(1)
         p(k+1) = p(2)
         p(k+2) = p(3)
         p(k+3) = p(4)
  310 continue
      return
  400 continue
      return
c ......................................................................
      end
      subroutine elmt08_mLRH(e,x,u,v,p,s,tl,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT01_I: Elemento de interface (mola) de 4 nos                  *
c *   ------  Para o caso 2D                                           * 
c *                                                                    *
c *    4----------3     |                                              *
c *    %          %     h -> espessura                                 *
c *    1----------2    _|_                                             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      common /gauss/ pg,wg
      integer ndm,nst,nel,isw
      integer nint,i,j,k,l,lx,ly,pi,NPI
      real*8  e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*)
      real*8  h(4),hx(4),hy(4),xj(2,2),xji(2,2),det,wt
      real*8  ym,pr,thic,a,d11,d12,d22,d33,tm,ponto,peso
      real*8  pg(10,10),wg(10,10),N(2,8),ri,si,Lc1,Lc2
      real*8  eps(3),tx,ty,txy,aux,ZERO,aux1,aux2,d11Aux,d22Aux
      parameter (ZERO = 1.d-14)
      data    nint/2/
c ......................................................................
      goto (100,200,200,400) isw
c ======================================================================
c
c.... Imput material properties:
c
  100 continue
c     e(1) = parametro de penalidade para direcao normal
c     e(2) = parametro de penalidade para direcao tangencial
c     e(3) = espessura do elemento de interface
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
c
c ... Matriz constitutiva:
c ......................................................................
      d11Aux    = e(1) !d11
      d22Aux    = e(2) !d22
      thic   = e(3) !h
c
c ... Matriz de rigidez do Elemento de interface:
c ......................................................................
      do i = 1, nst
         do j = 1, nst
           s(i,j) = 0.d0
         enddo
      enddo
c
c comprimentos das linhas inferior e superior como se segue
c	  
            Lc1 = DSQRT( (x(1,2)-x(1,1))*(x(1,2)-x(1,1)) +
     #	 (x(2,2)-x(2,1))*(x(2,2)-x(2,1)) )*0.5

            Lc2 = DSQRT( (x(1,4)-x(1,3))*(x(1,4)-x(1,3)) +
     #	 (x(2,4)-x(2,3))*(x(2,4)-x(2,3)) )*0.5
        
        thic = Lc1/10000000000.0         
        d11 = ((d11Aux)/(2.0*(1.0+d22Aux)))/thic
        d22 = (d11Aux)/thic        
           
        if (isw .eq. 3) goto 300 !para calculo das tensões
        print*,'d11 = ',d11
        print*,'d22 = ',d22
c
c	Matriz Kij:
c
c
c	 Matriz de rigidez:               
c		
        s(1,1)=d11*Lc1*0.5
 	s(2,2)=d22*Lc1*0.5
 	s(3,3)=d11*Lc1*0.5
 	s(4,4)=d22*Lc1*0.5
	s(5,5)=d11*Lc2*0.5
 	s(6,6)=d22*Lc2*0.5
 	s(7,7)=d11*Lc2*0.5
 	s(8,8)=d22*Lc2*0.5
 	s(1,7)=-d11*Lc2*0.5
 	s(2,8)=-d22*Lc2*0.5
 	s(3,5)=-d11*Lc2*0.5
	s(4,6)=-d22*Lc2*0.5
 	s(5,3)=-d11*Lc1*0.5
 	s(6,4)=-d22*Lc1*0.5
        s(7,1)=-d11*Lc1*0.5
 	s(8,2)=-d22*Lc1*0.5

c        s(3,1)=-d11*Lc1*0.5
c 	s(1,3)=-d11*Lc1*0.5

c        s(7,5)=d11*Lc1*0.5
c 	s(5,7)=d11*Lc1*0.5


c ======================================================================
c
c.... Tensoes nodais:       
c ......................................................................
  300 continue
      aux1 = ABS(u(8)-u(2))
      aux2 = ABS(u(6)-u(4))
      if ((aux1.LE.1.0e-10).OR.(aux2.LE.1.0e-10))then
      else 
        print*, 'Penetração no Elemento ',nel  
      endif
      aux1 = ABS(u(7)-u(1))
      aux2 = ABS(u(5)-u(3))
      if ((aux1.LE.1.0e-10).OR.(aux2.LE.1.0e-10))then
      else
        print*, 'Escorregamento no Elemento ',nel 
      endif
      p(1) = 0.0
      p(3) = d11*((u(7)-u(1)+u(5)-u(3)))*Lc1 !xy      
      p(2) = d22*((u(8)-u(2)+u(6)-u(4)))*Lc1 !y
      print*,'p2 y = ',p(2)
      print*,'p3 xy = ',p(3)
      print*,' ---'
      p(4) = p(1)
      p(5) = p(2) 
      p(6) = p(3)

      p(7) = p(1)  
      p(8) = p(2)
      p(9) = p(3)

      p(10) = p(1)
      p(11) = p(2)
      p(12) = p(3)

c      p(4) = d11*(u(3)-u(5))/Lc1
c      p(5) = d22*(u(4)-u(6))/Lc1 
c      p(6) = 0.0

c      p(7) = d11*(u(3)-u(5))/Lc1 
c      p(8) = d22*(u(4)-u(6))/Lc1
c      p(9) = 0.0

c      p(10) = d11*(u(1)-u(7))/Lc1 
c      p(11) = d22*(u(2)-u(8))/Lc1
c      p(12) = 0.0

      return
  400 continue
      return
c ......................................................................
      end
c IE5 - ALGA
      subroutine elmt08_miE_5(e,x,u,v,p,s,tl,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT01_I: Elemento de interface (mola) de 4 nos                  *
c *   ------  Para o caso 2D                                           * 
c *                                                                    *
c *    4----------3     |                                              *
c *    %     %    %     h -> espessura                                 *
c *    1----------2    _|_                                             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      common /gauss/ pg,wg
      integer ndm,nst,nel,isw
      integer nint,i,j,k,l,lx,ly,pi,NPI
      real*8  e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*)
      real*8  h(4),hx(4),hy(4),xj(2,2),xji(2,2),det,wt
      real*8  ym,pr,thic,a,d11,d12,d22,d33,tm,ponto,peso
      real*8  pg(10,10),wg(10,10),N(2,8),ri,si,Lc1,Lc2,c1,c2,c4
      real*8  eps(3),tx,ty,txy,aux,ZERO,aux1,aux2,d11Aux,d22Aux
      parameter (ZERO = 1.d-14)
      data    nint/2/
c ......................................................................
      goto (100,200,200,400) isw
c ======================================================================
c
c.... Imput material properties:
c
  100 continue
c     e(1) = parametro de penalidade para direcao normal
c     e(2) = parametro de penalidade para direcao tangencial
c     e(3) = espessura do elemento de interface
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
c
c ... Matriz constitutiva:
c ......................................................................
      d11Aux    = e(1) !d11
      d22Aux    = e(2) !d22
      thic   = e(3) !h
c
c ... Matriz de rigidez do Elemento de interface:
c ......................................................................
      do i = 1, nst
         do j = 1, nst
           s(i,j) = 0.d0
         enddo
      enddo
c
c comprimentos das linhas inferior e superior como se segue
c	  
            Lc1 = DSQRT( (x(1,2)-x(1,1))*(x(1,2)-x(1,1)) +
     #	 (x(2,2)-x(2,1))*(x(2,2)-x(2,1)) )*0.5

            Lc2 = DSQRT( (x(1,4)-x(1,3))*(x(1,4)-x(1,3)) +
     #	 (x(2,4)-x(2,3))*(x(2,4)-x(2,3)) )*0.5
        
        thic = Lc1/1000000000.0         
        d11 = ((d11Aux)/(2.0*(1.0+d22Aux)))/thic
        d22 = (d11Aux)/thic        
        c1 = d11*Lc1*0.5
        c2 = d22*Lc1/6.0
        c4 = 4.0*c2	  

        if (isw .eq. 3) goto 300 !para calculo das tensões
        print*,'d11 = ',d11
        print*,'d22 = ',d22
c
c	Matriz Kij:
c
c
c	 Matriz de rigidez:               
c	
        s(1,1)=c1
 	s(2,2)=c2+c4*0.25
 	s(3,3)=c1
 	s(4,4)=c2+c4*0.25
	s(5,5)=c1
 	s(6,6)=c2+c4*0.25
 	s(7,7)=c1
 	s(8,8)=c2+c4*0.25
 	s(1,7)=-c1
 	s(3,5)=-c1
 	s(5,3)=-c1
 	s(6,4)=-(c2+c4*0.25)
        s(7,1)=-c1
 	s(8,2)=-(c2+c4*0.25)
        s(2,4)= c4*0.25
        s(4,2)= c4*0.25
        s(6,2)= -c4*0.25
        s(8,4)= -c4*0.25
        s(2,6)= -c4*0.25
 	s(2,8)=-(c2+c4*0.25)
 	s(4,6)=-(c2+c4*0.25)
 	s(4,8)=-c4*0.25
 	s(6,8)= c4*0.25
 	s(8,6)= c4*0.25


        do i=1,8 
           write(*,'(a,7e15.5e3,7e15.5e3,7e15.5e3,7e15.5e3,
     #7e15.5e3,7e15.5e3,7e15.5e3,7e15.5e3)')
     #' {', s(i,1), s(i,2), s(i,3), s(i,4), s(i,5), s(i,6),
     #s(i,7), s(i,8)
           print*,' '
        enddo
c ======================================================================
c
c.... Tensoes nodais:       
c ......................................................................
  300 continue
      aux1 = ABS(u(8)-u(2))
      aux2 = ABS(u(6)-u(4))
      if ((aux1.LE.1.0e-10).OR.(aux2.LE.1.0e-10))then
      else 
        print*, 'Penetração no Elemento ',nel  
      endif
      aux1 = ABS(u(7)-u(1))
      aux2 = ABS(u(5)-u(3))
      if ((aux1.LE.1.0e-10).OR.(aux2.LE.1.0e-10))then
      else
        print*, 'Escorregamento no Elemento ',nel 
      endif
      p(1) = 0.0
      p(3) = c1*((u(7)-u(1)+u(5)-u(3))*0.5)/Lc1 !xy      
      p(2) = (c2+c4*0.5)*((u(8)-u(2)+u(6)-u(4)))/Lc1 !y
      p(4) = p(1)
      p(5) = p(2) 
      p(6) = p(3)
      p(7) = p(1)  
      p(8) = p(2)
      p(9) = p(3)
      p(10) = p(1)
      p(11) = p(2)
      p(12) = p(3)

      return
  400 continue
      return
c ......................................................................
      end
      subroutine elmt08_3d(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,ma,
     .                          nlit,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT02_I: Elemento de Interface para o Elemento tetraedrico      *
c *   ------                                                           *
c *  1--3                                                              *
c *  |\/ \  h -> espessura = 0.0                                       *
c *  |/4--6                                                            *
c *  2 | /                                                             *
c *   \|/                                                              *
c *    5                                                               *
c *	  !nst = 18                                                        *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - solucao anterior                                  *
c *     v(nst)     - solucao anterior incremental                      *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     tx(ntn)    - tensoes no elemento                               *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'            
      integer ndm,nst,nel,isw,ma,nlit
      integer i,j,k,l,q,i1,i2,i3,j1,j2,j3
      real*8 e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*),r(3,3)
      real*8 hx(4),hy(4),hz(4),eps(6),det,wt,espes(3),N(3,18)
      real*8 xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33,tm
      real*8 xji11,xji12,xji13,xji21,xji22,xji23,xji31,xji32,xji33
      real*8 ym,nu,a,b,c,a1,a2,a3,d11,d22,d33
      real*8 lDxDksi,lDyDksi,lDxDeta,lDyDeta,lDzDksi,lDzDeta
      real*8 uDxDksi,uDyDksi,uDxDeta,uDyDeta,uDzDksi,uDzDeta
      real*8 udet,ldet,uJ1,uJ2,uJ3,lJ1,lJ2,lJ3,peso,auxdet,paux(9)
      real*8 lxj11,lxj12,lxj21,lxj22,uxj11,uxj12,uxj21,uxj22,thic
      real*8 u0(*),du(18),tx1,tx2,tx3,tx(*)
      logical flagt    
c ......................................................................
      go to (100,200,300,200,200) isw
c ======================================================================
c
c.... Input material properties:
c
c ......................................................................
  100 continue
c     e(1) = parametro de penalidade para direçao normal
c     e(2) = parametro de penalidade para direçao xz
c     e(3) = parametro de penalidade para direçao yz
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
      flagt = .false.
c ......................................................................
c ... Propriedades do material:
c ......................................................................
c      d11 = e(1) !d11
c     d22 = e(2) !d22
c     d33 = e(3) !d33
c
c ... Calculo das espessuras entre nós homólogos
c 
      espes(1) = SQRT( (x(1,1)-x(1,4))*(x(1,1)-x(1,4)) +
     #	              (x(2,1)-x(2,4))*(x(2,1)-x(2,4)) +
     #                  (x(3,1)-x(3,4))*(x(3,1)-x(3,4)) )
      espes(2) = SQRT( (x(1,2)-x(1,5))*(x(1,2)-x(1,5)) +
     #	              (x(2,2)-x(2,5))*(x(2,2)-x(2,5)) +
     #                  (x(3,2)-x(3,5))*(x(3,2)-x(3,5)) )
      espes(3) = SQRT( (x(1,3)-x(1,6))*(x(1,3)-x(1,6)) +
     #	              (x(2,3)-x(2,6))*(x(2,3)-x(2,6)) +
     #                  (x(3,3)-x(3,6))*(x(3,3)-x(3,6)) )
      peso = 1.d0/3.d0
c
c ... Matriz Jacobiana:
c
c     upper
c
      uDxDksi = x(2,2)-x(2,1) !xj11
      uDyDksi = x(3,3)-x(3,1) !xj12
      uDzDksi = x(3,2)-x(3,1) !xj13

      uDxDeta = x(2,3)-x(2,1) !xj21
      uDyDeta = x(1,3)-x(1,1) !xj22
      uDzDeta = x(1,2)-x(1,1) !xj23
     	
      uJ1 = uDxDksi*uDyDksi-uDzDksi*uDxDeta
      uJ2 = uDzDksi*uDyDeta-uDzDeta*uDyDksi
      uJ3 = uDzDeta*uDxDeta-uDxDksi*uDyDeta
      udet = (uJ1 + uJ2 + uJ3)
      if (udet .le. 0.d0) goto 1000
c
c
c     lower
c
      lDxDksi = x(2,5)-x(2,4) !xj11
      lDyDksi = x(3,6)-x(3,4) !xj12
      lDzDksi = x(3,5)-x(3,4) !xj13

      lDxDeta = x(2,6)-x(2,4) !xj21
      lDyDeta = x(1,6)-x(1,4) !xj22
      lDzDeta = x(1,5)-x(1,4) !xj23
     	
      lJ1 = lDxDksi*lDyDksi-lDzDksi*lDxDeta
      lJ2 = lDzDksi*lDyDeta-lDzDeta*lDyDksi
      lJ3 = lDzDeta*lDxDeta-lDxDksi*lDyDeta
      ldet = (lJ1 + lJ2 + lJ3)
	if (ldet .le. 0.d0)  goto 1000
c ...................................................................... 
c ... Propriedades do material: 
c ...................................................................... 
      thic = e(4)                     
      if (thic .eq. 0.d0) thic = 0.5d0*udet/1000000000.d0
      d11 = (e(1)/(2.0*(1.0+e(2))))/thic   !d11   
      d22 = (e(1)/(2.0*(1.0+e(3))))/thic   !d22 
      d33 = e(1)/thic                      !d33

c      if (isw .eq. 3) goto 300  
c
c     Matriz de rigidez:                     
c
	wt = udet/24.d0
      do i = 1, 18
         do j = 1, 18
            s(i,j) = 0.d0   
         enddo
      enddo
      s(1,1) = 2.d0*d22*wt
      s(2,2) = 2.d0*d33*wt
      s(3,3) = 2.d0*d11*wt
      do i = 1, 6
         do j = 1, 6
            i1=3*(i-1)+1
            i2=3*(i-1)+2
            i3=3*(i-1)+3
            j1=3*(j-1)+1
            j2=3*(j-1)+2
            j3=3*(j-1)+3 
            if ((i .gt. 3 .and. j .le. 3) .or. 
     .      (j .gt. 3 .and. i .le. 3))  then
               s(i1,j1) = -s(1,1)/2.d0
               s(i2,j2) = -s(2,2)/2.d0
               s(i3,j3) = -s(3,3)/2.d0  
              goto 205  
            endif
            if (i .eq. j)   then
               s(i1,j1) = s(1,1)
               s(i2,j2) = s(2,2)
               s(i3,j3) = s(3,3)
               goto 205          
            endif                 
            s(i1,j1) = s(1,1)/2.d0
            s(i2,j2) = s(2,2)/2.d0
            s(i3,j3) = s(3,3)/2.d0
  205       continue
         enddo
      enddo
      j = 0
      do i = 10, 18
          j = j+1
          s(i,j) = -s(i,i)
          s(j,i) = -s(i,i)
      enddo
      call rotmatrix(x,r)
      call irotatek3d(s,r,nst)
c
c ..... Forcas Internas
c
      call rotate(u,r,u,.false.) 
      call rotate(u(4),r,u(4),.false.)
      call rotate(u(7),r,u(7),.false.)
      call rotate(u(10),r,u(10),.false.)
      call rotate(u(13),r,u(13),.false.)
      call rotate(u(16),r,u(16),.false.)
      if (nlit .eq. 1)  then
          do i = 1, 18
             u(i) = 0.d0
          enddo
      else 
         do 241 i = 1, 18
            du(i) = u(i) - u0(i)
            u0(i) = u(i)
  241    continue
         call istress3d(d11,d22,d33,udet,u,du,tx,0.d0,flagt)
         tx1 = tx(1)
         tx2 = tx(2)
         tx3 = tx(3)
      endif
      call rotate(u,r,u,.true.)
      call rotate(u(4),r,u(4),.true.)
      call rotate(u(7),r,u(7),.true.)
      call rotate(u(10),r,u(10),.true.)
      call rotate(u(13),r,u(13),.true.)
      call rotate(u(16),r,u(16),.true.)
      do 245 i = 1, 18
         p(i) = 0.d0
  245 continue
c ... Modifica a matriz de rigidez no caso de tracao
      if (flagt)   then
         do i = 1, nst
            do j = 1, nst
               s(i,j) = 0.d0
            enddo
         enddo
c         d11 = (e(1)/(2.0*(1.0+e(2))))*thic   !d11   
c         d22 = (e(1)/(2.0*(1.0+e(3))))*thic   !d22 
c         d33 = e(1)*thic                      !d33
c
c     Matriz de rigidez:                     
c
c	   wt = udet/24.d0
c         do i = 1, 18
c            do j = 1, 18
c               s(i,j) = 0.d0   
c            enddo
c         enddo
c         s(1,1) = 2*d22*wt
c         s(2,2) = 2*d33*wt
c         s(3,3) = 2*d11*wt
c         do i = 1, 6
c            do j = 1, 6
c               i1=3*(i-1)+1
c               i2=3*(i-1)+2
c               i3=3*(i-1)+3
c               j1=3*(j-1)+1
c               j2=3*(j-1)+2
c               j3=3*(j-1)+3 
c               if ((i .gt. 3 .and. j .le. 3) .or. 
c     .        (j .gt. 3 .and. i .le. 3))  then
c                  s(i1,j1) = -s(1,1)/2.d0
c                  s(i2,j2) = -s(2,2)/2.d0
c                  s(i3,j3) = -s(3,3)/2.d0  
c                  goto 270  
c               endif
c               if (i .eq. j)   then
c                  s(i1,j1) = s(1,1)
c                  s(i2,j2) = s(2,2)
c                  s(i3,j3) = s(3,3)
c                  goto 270          
c               endif                 
c               s(i1,j1) = s(1,1)/2.d0
c               s(i2,j2) = s(2,2)/2.d0
c               s(i3,j3) = s(3,3)/2.d0
c  270          continue
c            enddo
c         enddo
c         j = 0
c         do i = 10, 18
c             j = j+1
c             s(i,j) = -s(i,i)
c             s(j,i) = -s(i,i)
c         enddo
c         call rotmatrix(x,r)
c         call irotatek3d(s,r,nst)
         if ( nlit .eq. 2) goto 280
      endif
c      do i = 1, 18
c         p(i) = 0.d0
c      enddo
      do i = 1, 6
         do j = 1, 6
         p(3*i-2) = p(3*i-2)  + s(3*i-2,3*j-2)*u(3*j-2)
         p(3*i-1) = p(3*i-1)  + s(3*i-1,3*j-1)*u(3*j-1)
         p(3*i)   = p(3*i)    + s(3*i ,3*j)*u(3*j)
         enddo
      enddo
280   continue
      return
c ======================================================================
c
c ... Tensoes nodais:        
c
c ......................................................................
  300 continue
      do i = 1, 36
         p(i) = 0.d0
      enddo
c      do i = 1, 6
c         eps(i) = 0.d0
c      enddo
      call rotmatrix(x,r)
c      call rotate(u,r,u,.false.) 
c      call rotate(u(4),r,u(4),.false.)
c      call rotate(u(7),r,u(7),.false.)
c      call rotate(u(10),r,u(10),.false.)  
c      call rotate(u(13),r,u(13),.false.) 
c      call rotate(u(16),r,u(16),.false.)  

c      eps(3) = (u(12)+u(15)+u(18)-u(3)-u(6)-u(9))/3.d0
c      eps(6) = (u(10)+u(13)+u(16)-u(1)-u(4)-u(7))/3.d0
c      eps(5) = (u(11)+u(14)+u(17)-u(2)-u(5)-u(8))/3.d0
c      p(3) = d11*eps(3) ! sigmaz
c      p(6) = d22*eps(6) ! talzx
c      p(5) = d33*eps(5) ! talzy
       p(3) = tx(1)
       p(6) = tx(2)
       p(5) = tx(3)
c      write(*,*)'*********'
c      write(*,*)d11*eps(3)
c      write(*,*)'*********'

      paux(1) = p(1)
      paux(2) = p(4)
      paux(3) = p(6)
      paux(4) = p(4)
      paux(5) = p(2)
      paux(6) = p(5)
      paux(7) = p(6)
      paux(8) = p(5)
      paux(9) = p(3)
      call irotatek3d(paux,r,3)  
      p(1) = paux(1)
      p(2) = paux(5) 
      p(3) = paux(9)
      p(4) = paux(2)
      p(5) = paux(6) 
      p(6) = paux(3)
      do i = 1, 5
         p(6*i+1) = p(1)
         p(6*i+2) = p(2)
         p(6*i+3) = p(3)
         p(6*i+4) = p(4)
         p(6*i+5) = p(5)
         p(6*i+6) = p(6)
      enddo
      return      
c ======================================================================
c
c ... Forcas nodais equivalentes:        
c
c ......................................................................
  400 continue
c
c ... Matriz Jacobiana:
c 
c      p(1)  =  0.d0
c      p(2)  =  0.d0
c      p(3)  =  e(3)*det/24.d0
c      p(4)  =  0.d0
c      p(5)  =  0.d0
c      p(6)  =  e(3)*det/24.d0
c      p(7)  =  0.d0
c      p(8)  =  0.d0
c      p(9)  =  e(3)*det/24.d0
c      p(10) =  0.d0
c      p(11) =  0.d0
c      p(12) =  e(3)*det/24.d0           
      return
c ======================================================================
c
c ... 
c
c ......................................................................
  500 continue
      return
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT09: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
      subroutine elmt08_3dhexa(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,ma,
     .                          nlit,isw,plastic)
c **********************************************************************
c *                                                                    *
c *   ELMT02_I: Elemento de Interface hexaedro mecanico                *
c *   ------                                                           *
c *  1--4                                                              *
c *  |\ |\  h -> espessura = 0.0                                       *
c *  | 2--3                                                            *
c *  5-|-8|                                                            *
c *   \|  \                                                            *
c *    6--7                                                            *
c *	  !nst = 24                                                      *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - solucao anterior                                  *
c *     v(nst)     - solucao anterior incremental                      *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     tx(ntn)    - tensoes no elemento                               *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'   
      common /gauss/ pg,wg       
      integer ndm,nst,nel,isw,ma,nlit,nint
      integer i,j,k,l,q,i1,i2,i3,j1,j2,j3,lx,ly
      real*8 e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*),r(3,3)
      real*8 hx(4),hy(4),hz(4),eps(6),det,wt,h(8)
      real*8 xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33,tm
      real*8 xji11,xji12,xji13,xji21,xji22,xji23,xji31,xji32,xji33
      real*8 ym,nu,a,b,c,a1,a2,a3,d11,d22,d33
      real*8 lDxDksi,lDyDksi,lDxDeta,lDyDeta,lDzDksi,lDzDeta
      real*8 uDxDksi,uDyDksi,uDxDeta,uDyDeta,uDzDksi,uDzDeta
      real*8 udet,ldet,uJ1,uJ2,uJ3,lJ1,lJ2,lJ3,peso,auxdet,paux(9)
      real*8 lxj11,lxj12,lxj21,lxj22,uxj11,uxj12,uxj21,uxj22,thic
      real*8 u0(*),du(nst),tx1,tx2,tx3,tx(*),xj(2,2),xji(2,2)
      real*8 pg(10,10),wg(10,10),ri,si,plastic(2)
      real*8 x1,x2,x3,area
      logical flagt    
      data    nint/2/
c ......................................................................
      go to (100,200,300,200,500) isw
c ======================================================================
c
c.... Input material properties:
c
c ......................................................................
  100 continue
c     e(1) = parametro de penalidade para direçao normal
c     e(2) = parametro de penalidade para direçao xz
c     e(3) = parametro de penalidade para direçao yz
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
      flagt = .false.
c ......................................................................
c ... Propriedades do material:
c ......................................................................
c      d11 = e(1) !d11
c     d22 = e(2) !d22
c     d33 = e(3) !d33
c
c ...................................................................... 
c      if (isw .eq. 3) goto 300  
c
c     Matriz de rigidez:                     
c
      call rotmatrix(x,r)
      call rotate(x(1,1),r,x(1,1),.false.) 
      call rotate(x(1,2),r,x(1,2),.false.)
      call rotate(x(1,3),r,x(1,3),.false.)
      call rotate(x(1,4),r,x(1,4),.false.)

c ... Propriedades do material: 
c ...................................................................... 
c      thic = e(4)                     
c      if (thic .eq. 0.d0) thic = 1.e-08
c      thic = 1.e-08
c     Determinacao da area do quadrilatero da interface
c      v1x = x(1,1) - x(1,2)
c      v1y = x(2,1) - x(2,2)
c      v1z = x(3,1) - x(3,2)
c      v2x = x(1,3) - x(1,2)
c      v2y = x(2,3) - x(2,2)
c      v2z = x(3,3) - x(3,2)
      x1 = (x(3,3) - x(3,2))*(x(2,1) - x(2,2))-(x(2,3) - 
     .x(2,2))*( x(3,1) - x(3,2))
      x2 = (x(3,1) - x(3,2))*(x(1,3) - x(1,2)) -
     .(x(1,1) - x(1,2))*(x(3,3) - x(3,2))
      x3 = (x(1,1) - x(1,2))*(x(2,3) - x(2,2)) -
     .(x(1,3) - x(1,2))*(x(2,1) - x(2,2))
      area = sqrt(x1*x1+x2*x2+x3*x3)
      thic = area*1.0e-04
      d11 = (e(1)/(2.0*(1.0+e(2))))/thic   !d11   
      d22 = (e(1)/(2.0*(1.0+e(2))))/thic   !d22 
      d33 = e(1)/thic                      !d33

      do i = 1, nst
      do j = 1, nst
         s(i,j) = 0.d0
      enddo
      enddo
      do 260 lx = 1, nint
         ri = pg(lx,nint)
         do 250 ly = 1, nint
            si = pg(ly,nint)
            call sfquad4(h,hx,hy,ri,si,.true.,.true.)
            call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)
            wt = wg(lx,nint)*wg(ly,nint)*det
            h(5) = -h(1)
            h(6) = -h(2)
            h(7) = -h(3)
            h(8) = -h(4)
           do 230 i = 1, 8
           i1 = (i-1)*3+1
           i2 = i1 + 1
           i3 = i2 + 1
           do 220 j = 1, 8
              j1 = (j-1)*3+1
              j2 = j1 + 1
              j3 = j2 + 1
              s(i1,j1) = s(i1,j1) + h(i)*d22*h(j)*wt
              s(i2,j2) = s(i2,j2) + h(i)*d33*h(j)*wt
              s(i3,j3) = s(i3,j3) + h(i)*d11*h(j)*wt
c ................................................................
  220      continue
  230      continue
  250    continue
  260 continue
c      call rotmatrix(x,r)
      call irotatek3d(s,r,nst)
c
c ..... Forcas Internas
c
      call rotate(u(1),r,u(1),.false.) 
      call rotate(u(4),r,u(4),.false.)
      call rotate(u(7),r,u(7),.false.)
      call rotate(u(10),r,u(10),.false.)
      call rotate(u(13),r,u(13),.false.)
      call rotate(u(16),r,u(16),.false.)
      call rotate(u(19),r,u(19),.false.)
      call rotate(u(22),r,u(22),.false.)
      if (nlit .eq. 1)  then
c      if (nlit .eq. 1 .or. nlit .eq. 2)  then
          do i = 1, nst
             u(i) = 0.d0
          enddo
      else 
         do 241 i = 1, nst
            du(i) = u(i) - u0(i)
            u0(i) = u(i)
  241    continue
         call istress3d_hexa(d11,d22,d33,0.d0,0.d0,area,u,du,tx,
     .                                                  0.d0,flagt)
         tx1 = tx(1)
         tx2 = tx(2)
         tx3 = tx(3)
      endif
      call rotate(u,r,u,.true.)
      call rotate(u(4),r,u(4),.true.)
      call rotate(u(7),r,u(7),.true.)
      call rotate(u(10),r,u(10),.true.)
      call rotate(u(13),r,u(13),.true.)
      call rotate(u(16),r,u(16),.true.)
      call rotate(u(19),r,u(19),.true.)
      call rotate(u(22),r,u(22),.true.)
      do 245 i = 1, nst
         p(i) = 0.d0
  245 continue
c ... Modifica a matriz de rigidez no caso de tracao
      if (flagt)   then
         do i = 1, nst
            do j = 1, nst
               s(i,j) = 0.d0
            enddo
         enddo
         return
      endif
c ......................................................................      
      do i = 1, 8
         do j = 1, 8
         p(3*i-2) = p(3*i-2)  + s(3*i-2,3*j-2)*u(3*j-2)
         p(3*i-1) = p(3*i-1)  + s(3*i-1,3*j-1)*u(3*j-1)
         p(3*i)   = p(3*i)    + s(3*i ,3*j)*u(3*j)
         enddo
      enddo
      return
c ======================================================================
c
c ... Tensoes nodais:        
c
c ......................................................................
  300 continue
      do i = 1, 36
         p(i) = 0.d0
      enddo
c      do i = 1, 6
c         eps(i) = 0.d0
c      enddo
      call rotmatrix(x,r)
c      call rotate(u,r,u,.false.) 
c      call rotate(u(4),r,u(4),.false.)
c      call rotate(u(7),r,u(7),.false.)
c      call rotate(u(10),r,u(10),.false.)  
c      call rotate(u(13),r,u(13),.false.) 
c      call rotate(u(16),r,u(16),.false.)  

c      eps(3) = (u(12)+u(15)+u(18)-u(3)-u(6)-u(9))/3.d0
c      eps(6) = (u(10)+u(13)+u(16)-u(1)-u(4)-u(7))/3.d0
c      eps(5) = (u(11)+u(14)+u(17)-u(2)-u(5)-u(8))/3.d0
c      p(3) = d11*eps(3) ! sigmaz
c      p(6) = d22*eps(6) ! talzx
c      p(5) = d33*eps(5) ! talzy
       p(3) = tx(1)
       p(6) = tx(2)
       p(5) = tx(3)
c      write(*,*)'*********'
c      write(*,*)d11*eps(3)
c      write(*,*)'*********'

      paux(1) = p(1)
      paux(2) = p(4)
      paux(3) = p(6)
      paux(4) = p(4)
      paux(5) = p(2)
      paux(6) = p(5)
      paux(7) = p(6)
      paux(8) = p(5)
      paux(9) = p(3)
      call irotatek3d(paux,r,3)  
      p(1) = paux(1)
      p(2) = paux(5) 
      p(3) = paux(9)
      p(4) = paux(2)
      p(5) = paux(6) 
      p(6) = paux(3)
      do i = 1, 7
         p(6*i+1) = p(1)
         p(6*i+2) = p(2)
         p(6*i+3) = p(3)
         p(6*i+4) = p(4)
         p(6*i+5) = p(5)
         p(6*i+6) = p(6)
      enddo
      return      
c ======================================================================
c
c ... Forcas nodais equivalentes:        
c
c ......................................................................
  400 continue
c
c ... Matriz Jacobiana:
c 
c      p(1)  =  0.d0
c      p(2)  =  0.d0
c      p(3)  =  e(3)*det/24.d0
c      p(4)  =  0.d0
c      p(5)  =  0.d0
c      p(6)  =  e(3)*det/24.d0
c      p(7)  =  0.d0
c      p(8)  =  0.d0
c      p(9)  =  e(3)*det/24.d0
c      p(10) =  0.d0
c      p(11) =  0.d0
c      p(12) =  e(3)*det/24.d0           
      return
c ======================================================================
c
c ... 
c
c ......................................................................
  500 continue
c .... Define a tensao local normal no elemento de interface e se houve
c      rompimento
      plastic(1) = tx(1)
      plastic(2) = 0.d0
      if (tx(1) .ge. 0.d0)  plastic(2) = 1.d0
      return
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT11: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
      subroutine elmt08_3dhexa_grampo(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,
     .                          ma,nlit,isw,plastic)
c **********************************************************************
c *                                                                    *
c *   ELMT08_I: Elemento de Interface hexaedro mecanico -grampo        *
c *   ------                                                           *
c *  1--4                                                              *
c *  |\ |\  h -> espessura = 0.0                                       *
c *  | 2--3                                                            *
c *  5-|-8|                                                            *
c *   \|  \                                                            *
c *    6--7                                                            *
c *	  !nst = 24                                                      *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - solucao anterior                                  *
c *     v(nst)     - solucao anterior incremental                      *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     tx(ntn)    - tensoes no elemento                               *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'   
      common /gauss/ pg,wg       
      integer ndm,nst,nel,isw,ma,nlit,nint
      integer i,j,k,l,q,i1,i2,i3,j1,j2,j3,lx,ly
      real*8 e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*),r(3,3)
      real*8 hx(4),hy(4),hz(4),eps(6),det,wt,h(8)
      real*8 xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33,tm
      real*8 xji11,xji12,xji13,xji21,xji22,xji23,xji31,xji32,xji33
      real*8 ym,nu,a,b,c,a1,a2,a3,d11,d22,d33
      real*8 lDxDksi,lDyDksi,lDxDeta,lDyDeta,lDzDksi,lDzDeta
      real*8 uDxDksi,uDyDksi,uDxDeta,uDyDeta,uDzDksi,uDzDeta
      real*8 udet,ldet,uJ1,uJ2,uJ3,lJ1,lJ2,lJ3,peso,auxdet,paux(9)
      real*8 lxj11,lxj12,lxj21,lxj22,uxj11,uxj12,uxj21,uxj22,thic
      real*8 u0(*),du(nst),tx1,tx2,tx3,tx(*),xj(2,2),xji(2,2)
      real*8 pg(10,10),wg(10,10),ri,si,plastic(2)
      real*8 x1,x2,x3,area
      logical flagt    
      data    nint/2/
c ......................................................................
      go to (100,200,300,200,500) isw
c ======================================================================
c
c.... Input material properties:
c
c ......................................................................
  100 continue
c     e(1) = parametro de penalidade para direçao normal
c     e(2) = parametro de penalidade para direçao xz
c     e(3) = parametro de penalidade para direçao yz
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
      flagt = .false.
c ......................................................................
c ... Propriedades do material:
c ......................................................................
c      d11 = e(1) !d11
c     d22 = e(2) !d22
c     d33 = e(3) !d33
c

c     Matriz de rigidez:                     
c
      call rotmatrix(x,r)
      call rotate(x(1,1),r,x(1,1),.false.) 
      call rotate(x(1,2),r,x(1,2),.false.)
      call rotate(x(1,3),r,x(1,3),.false.)
      call rotate(x(1,4),r,x(1,4),.false.)
c     Determinacao da area do quadrilatero da interface
c      v1x = x(1,1) - x(1,2)
c      v1y = x(2,1) - x(2,2)
c      v1z = x(3,1) - x(3,2)
c      v2x = x(1,3) - x(1,2)
c      v2y = x(2,3) - x(2,2)
c      v2z = x(3,3) - x(3,2)
      x1 = (x(3,3) - x(3,2))*(x(2,1) - x(2,2))-(x(2,3) - 
     .x(2,2))*( x(3,1) - x(3,2))
      x2 = (x(3,1) - x(3,2))*(x(1,3) - x(1,2)) -
     .(x(1,1) - x(1,2))*(x(3,3) - x(3,2))
      x3 = (x(1,1) - x(1,2))*(x(2,3) - x(2,2)) -
     .(x(1,3) - x(1,2))*(x(2,1) - x(2,2))
      area = sqrt(x1*x1+x2*x2+x3*x3)
      thic = area/1.0e+4
c ...................................................................... 
c ... Propriedades do material: 
c ...................................................................... 
c      thic = e(4)                     
c      if (thic .eq. 0.d0) thic = 1.e-08
c      thic = 1.e-08
      d11 = (e(1)/(2.0*(1.0+e(2))))/thic   !d11   
      d22 = (e(1)/(2.0*(1.0+e(2))))/thic   !d22 
      d33 = e(1)/thic                      !d33

c      if (isw .eq. 3) goto 300  
c
      do i = 1, nst
      do j = 1, nst
         s(i,j) = 0.d0
      enddo
      enddo
      do 260 lx = 1, nint
         ri = pg(lx,nint)
         do 250 ly = 1, nint
            si = pg(ly,nint)
            call sfquad4(h,hx,hy,ri,si,.true.,.true.)
            call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)
            wt = wg(lx,nint)*wg(ly,nint)*det
            h(5) = -h(1)
            h(6) = -h(2)
            h(7) = -h(3)
            h(8) = -h(4)
           do 230 i = 1, 8
           i1 = (i-1)*3+1
           i2 = i1 + 1
           i3 = i2 + 1
           do 220 j = 1, 8
              j1 = (j-1)*3+1
              j2 = j1 + 1
              j3 = j2 + 1
              s(i1,j1) = s(i1,j1) + h(i)*d22*h(j)*wt
              s(i2,j2) = s(i2,j2) + h(i)*d33*h(j)*wt
              s(i3,j3) = s(i3,j3) + h(i)*d11*h(j)*wt
c ................................................................
  220      continue
  230      continue
  250    continue
  260 continue
c      call rotmatrix(x,r)
      call irotatek3d(s,r,nst)
c
c ..... Forcas Internas
c
      call rotate(u(1),r,u(1),.false.) 
      call rotate(u(4),r,u(4),.false.)
      call rotate(u(7),r,u(7),.false.)
      call rotate(u(10),r,u(10),.false.)
      call rotate(u(13),r,u(13),.false.)
      call rotate(u(16),r,u(16),.false.)
      call rotate(u(19),r,u(19),.false.)
      call rotate(u(22),r,u(22),.false.)
      if (nlit .eq. 1)  then
c      if (nlit .eq. 1 .or. nlit .eq. 2)  then
          do i = 1, nst
             u(i) = 0.d0
          enddo
      else 
         do 241 i = 1, nst
            du(i) = u(i) - u0(i)
            u0(i) = u(i)
  241    continue
         call istress3d_hexa(d11,d22,d33,e(3),e(4),area,u,
     .                                         du,tx,0.d0,flagt)
         tx1 = tx(1)
         tx2 = tx(2)
         tx3 = tx(3)
      endif
      call rotate(u,r,u,.true.)
      call rotate(u(4),r,u(4),.true.)
      call rotate(u(7),r,u(7),.true.)
      call rotate(u(10),r,u(10),.true.)
      call rotate(u(13),r,u(13),.true.)
      call rotate(u(16),r,u(16),.true.)
      call rotate(u(19),r,u(19),.true.)
      call rotate(u(22),r,u(22),.true.)
      do 245 i = 1, nst
         p(i) = 0.d0
  245 continue
c ... Modifica a matriz de rigidez no caso de tracao
      if (flagt)   then
         d11 = e(3)
         d22 = e(4)
         d33 = e(4)
c
c     Matriz de rigidez:                     
c
         do i = 1, nst
         do j = 1, nst
            s(i,j) = 0.d0
         enddo
         enddo
         do 295 lx = 1, nint
            ri = pg(lx,nint)
            do 290 ly = 1, nint
               si = pg(ly,nint)
               call sfquad4(h,hx,hy,ri,si,.true.,.true.)
               call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)
               wt = wg(lx,nint)*wg(ly,nint)*det
               h(5) = -h(1)
               h(6) = -h(2)
               h(7) = -h(3)
               h(8) = -h(4)
               do 280 i = 1, 8
                  i1 = (i-1)*3+1
                  i2 = i1 + 1
                  i3 = i2 + 1
                  do 270 j = 1, 8
                     j1 = (j-1)*3+1
                     j2 = j1 + 1
                     j3 = j2 + 1
                     s(i1,j1) = s(i1,j1) + h(i)*d22*h(j)*wt
                     s(i2,j2) = s(i2,j2) + h(i)*d33*h(j)*wt
                     s(i3,j3) = s(i3,j3) + h(i)*d11*h(j)*wt
c ................................................................
  270             continue
  280          continue
  290       continue
  295    continue
         call irotatek3d(s,r,nst)
      endif
c ......................................................................      
      do i = 1, 8
         do j = 1, 8
         p(3*i-2) = p(3*i-2)  + s(3*i-2,3*j-2)*u(3*j-2)
         p(3*i-1) = p(3*i-1)  + s(3*i-1,3*j-1)*u(3*j-1)
         p(3*i)   = p(3*i)    + s(3*i ,3*j)*u(3*j)
         enddo
      enddo
      return
c ======================================================================
c
c ... Tensoes nodais:        
c
c ......................................................................
  300 continue
      do i = 1, 36
         p(i) = 0.d0
      enddo
c      do i = 1, 6
c         eps(i) = 0.d0
c      enddo
      call rotmatrix(x,r)
       p(3) = tx(1)
       p(6) = tx(2)
       p(5) = tx(3)
      paux(1) = p(1)
      paux(2) = p(4)
      paux(3) = p(6)
      paux(4) = p(4)
      paux(5) = p(2)
      paux(6) = p(5)
      paux(7) = p(6)
      paux(8) = p(5)
      paux(9) = p(3)
      call irotatek3d(paux,r,3)  
      p(1) = paux(1)
      p(2) = paux(5) 
      p(3) = paux(9)
      p(4) = paux(2)
      p(5) = paux(6) 
      p(6) = paux(3)
      do i = 1, 7
         p(6*i+1) = p(1)
         p(6*i+2) = p(2)
         p(6*i+3) = p(3)
         p(6*i+4) = p(4)
         p(6*i+5) = p(5)
         p(6*i+6) = p(6)
      enddo
      return      
c ======================================================================
c
c ... Forcas nodais equivalentes:        
c
c ......................................................................
  400 continue
c
c ... Matriz Jacobiana:
c 
c      p(1)  =  0.d0
c      p(2)  =  0.d0
c      p(3)  =  e(3)*det/24.d0
c      p(4)  =  0.d0
c      p(5)  =  0.d0
c      p(6)  =  e(3)*det/24.d0
c      p(7)  =  0.d0
c      p(8)  =  0.d0
c      p(9)  =  e(3)*det/24.d0
c      p(10) =  0.d0
c      p(11) =  0.d0
c      p(12) =  e(3)*det/24.d0           
      return
c ======================================================================
c
c ... 
c
c ......................................................................
  500 continue
c .... Define a tensao local normal no elemento de interface e se houve
c      rompimento
      plastic(1) = tx(1)
      plastic(2) = 0.d0
      if (tx(1) .ge. e(5))  plastic(2) = 1.d0
      return
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT11: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
      subroutine elmt09_mItetra(e,x,u,v,p,s,tl,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT02_I: Elemento de Interface para o Elemento tetraedrico      *
c *   ------                                                           *
c *  1--3                                                              *
c *  |\/ \  h -> espessura = 0.0                                       *
c *  |/4--6                                                            *
c *  2 | /                                                             *
c *   \|/                                                              *
c *    5                                                               *
c *	  !nst = 18                                                      *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - solucao anterior                                  *
c *     v(nst)     - solucao anterior incremental                      *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     tx(ntn)    - tensoes no elemento                               *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'            
      integer ndm,nst,nel,isw
      integer i,j,k,l,q,i1,i2,i3,j1,j2,j3
      real*8 e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),tl(*),r(3,3)
      real*8 hx(4),hy(4),hz(4),eps(6),det,wt,espes(3),N(3,18)
      real*8 xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33,tm
      real*8 xji11,xji12,xji13,xji21,xji22,xji23,xji31,xji32,xji33
      real*8 ym,nu,a,b,c,a1,a2,a3,tx,ty,tz,txy,txz,tyz,d11,d22,d33
      real*8 lDxDksi,lDyDksi,lDxDeta,lDyDeta,lDzDksi,lDzDeta
      real*8 uDxDksi,uDyDksi,uDxDeta,uDyDeta,uDzDksi,uDzDeta
      real*8 udet,ldet,uJ1,uJ2,uJ3,lJ1,lJ2,lJ3,peso,auxdet,paux(9)
      real*8 lxj11,lxj12,lxj21,lxj22,uxj11,uxj12,uxj21,uxj22,thic      
c ......................................................................
      go to (100,200,200,200,200) isw
c ======================================================================
c
c.... Input material properties:
c
c ......................................................................
  100 continue
c     e(1) = parametro de penalidade para direçao normal
c     e(2) = parametro de penalidade para direçao xz
c     e(3) = parametro de penalidade para direçao yz
      return
c ======================================================================
c
c.... Matrix K:
c
c ......................................................................
  200 continue
c ......................................................................
c ... Propriedades do material:
c ......................................................................
c      d11 = e(1) !d11
c     d22 = e(2) !d22
c     d33 = e(3) !d33
c
c ... Calculo das espessuras entre nós homólogos
c 
      espes(1) = SQRT( (x(1,1)-x(1,4))*(x(1,1)-x(1,4)) +
     #	              (x(2,1)-x(2,4))*(x(2,1)-x(2,4)) +
     #                  (x(3,1)-x(3,4))*(x(3,1)-x(3,4)) )
      espes(2) = SQRT( (x(1,2)-x(1,5))*(x(1,2)-x(1,5)) +
     #	              (x(2,2)-x(2,5))*(x(2,2)-x(2,5)) +
     #                  (x(3,2)-x(3,5))*(x(3,2)-x(3,5)) )
      espes(3) = SQRT( (x(1,3)-x(1,6))*(x(1,3)-x(1,6)) +
     #	              (x(2,3)-x(2,6))*(x(2,3)-x(2,6)) +
     #                  (x(3,3)-x(3,6))*(x(3,3)-x(3,6)) )
      peso = 1.d0/3.d0
c
c ... Matriz Jacobiana:
c
c     upper
c
      uDxDksi = x(2,2)-x(2,1) !xj11
      uDyDksi = x(3,3)-x(3,1) !xj12
      uDzDksi = x(3,2)-x(3,1) !xj13

      uDxDeta = x(2,3)-x(2,1) !xj21
      uDyDeta = x(1,3)-x(1,1) !xj22
      uDzDeta = x(1,2)-x(1,1) !xj23
     	
      uJ1 = uDxDksi*uDyDksi-uDzDksi*uDxDeta
      uJ2 = uDzDksi*uDyDeta-uDzDeta*uDyDksi
      uJ3 = uDzDeta*uDxDeta-uDxDksi*uDyDeta
      udet = (uJ1 + uJ2 + uJ3)
      if (udet .le. 0.d0) goto 1000
c
c
c     lower
c
      lDxDksi = x(2,5)-x(2,4) !xj11
      lDyDksi = x(3,6)-x(3,4) !xj12
      lDzDksi = x(3,5)-x(3,4) !xj13

      lDxDeta = x(2,6)-x(2,4) !xj21
      lDyDeta = x(1,6)-x(1,4) !xj22
      lDzDeta = x(1,5)-x(1,4) !xj23
     	
      lJ1 = lDxDksi*lDyDksi-lDzDksi*lDxDeta
      lJ2 = lDzDksi*lDyDeta-lDzDeta*lDyDksi
      lJ3 = lDzDeta*lDxDeta-lDxDksi*lDyDeta
      ldet = (lJ1 + lJ2 + lJ3)
	if (ldet .le. 0.d0)  goto 1000
c ...................................................................... 
c ... Propriedades do material: 
c ...................................................................... 
      thic = e(4)                     
      if (thic.eq.0.0) thic = 0.5*udet/1000000000.0
      d11 = (e(1)/(2.0*(1.0+e(2))))/thic   !d11   
      d22 = (e(1)/(2.0*(1.0+e(3))))/thic   !d22 
      d33 = e(1)/thic                      !d33

      if (isw .eq. 3) goto 300  
c
c     Matriz de rigidez:                     
c
	wt = udet/24.d0
      do i = 1, 18
         do j = 1, 18
            s(i,j) = 0.d0   
         enddo
      enddo
      s(1,1) = 2*d22*wt
      s(2,2) = 2*d33*wt
      s(3,3) = 2*d11*wt
      do i = 1, 6
         do j = 1, 6
            i1=3*(i-1)+1
            i2=3*(i-1)+2
            i3=3*(i-1)+3
            j1=3*(j-1)+1
            j2=3*(j-1)+2
            j3=3*(j-1)+3 
            if ((i .gt. 3 .and. j .le. 3) .or. 
     .      (j .gt. 3 .and. i .le. 3))  then
               s(i1,j1) = -s(1,1)/2.d0
               s(i2,j2) = -s(2,2)/2.d0
               s(i3,j3) = -s(3,3)/2.d0  
              goto 205  
            endif
            if (i .eq. j)   then
               s(i1,j1) = s(1,1)
               s(i2,j2) = s(2,2)
               s(i3,j3) = s(3,3)
               goto 205          
            endif                 
            s(i1,j1) = s(1,1)/2.d0
            s(i2,j2) = s(2,2)/2.d0
            s(i3,j3) = s(3,3)/2.d0
  205       continue
         enddo
      enddo
      j = 0
      do i = 10, 18
          j = j+1
          s(i,j) = -s(i,i)
          s(j,i) = -s(i,i)
      enddo
      call rotmatrix(x,r)
      call irotatek3d(s,r,nst)
      return
c ======================================================================
c
c ... Tensoes nodais:        
c
c ......................................................................
  300 continue
      do i = 1, 36
         p(i) = 0.d0
      enddo
      do i = 1, 6
         eps(i) = 0.d0
      enddo
      call rotmatrix(x,r)
      call rotate(u,r,u,.false.) 
      call rotate(u(4),r,u(4),.false.)
      call rotate(u(7),r,u(7),.false.)
      call rotate(u(10),r,u(10),.false.)  
      call rotate(u(13),r,u(13),.false.) 
      call rotate(u(16),r,u(16),.false.)  

      eps(3) = (u(12)+u(15)+u(18)-u(3)-u(6)-u(9))*udet/3.d0
      eps(6) = (u(10)+u(13)+u(16)-u(1)-u(4)-u(7))*udet/3.d0
      eps(5) = (u(11)+u(14)+u(17)-u(2)-u(5)-u(8))*udet/3.d0
      p(3) = d11*eps(3) ! sigmaz
      p(6) = d22*eps(6) ! talzx
      p(5) = d33*eps(5) ! talzy
c     
c      write(*,*)'*********'
c      write(*,*)d11*eps(3)
c      write(*,*)'*********'

      paux(1) = p(1)
      paux(2) = p(4)
      paux(3) = p(6)
      paux(4) = p(4)
      paux(5) = p(2)
      paux(6) = p(5)
      paux(7) = p(6)
      paux(8) = p(5)
      paux(9) = p(3)
      call irotatek3d(paux,r,3)  
      p(1) = paux(1)
      p(2) = paux(5) 
      p(3) = paux(9)
      p(4) = paux(2)
      p(5) = paux(6) 
      p(6) = paux(3)
      do i = 1, 5
         p(6*i+1) = p(1)
         p(6*i+2) = p(2)
         p(6*i+3) = p(3)
         p(6*i+4) = p(4)
         p(6*i+5) = p(5)
         p(6*i+6) = p(6)
      enddo
      return      
c ======================================================================
c
c ... Forcas nodais equivalentes:        
c
c ......................................................................
  400 continue
c
c ... Matriz Jacobiana:
c 
c      p(1)  =  0.d0
c      p(2)  =  0.d0
c      p(3)  =  e(3)*det/24.d0
c      p(4)  =  0.d0
c      p(5)  =  0.d0
c      p(6)  =  e(3)*det/24.d0
c      p(7)  =  0.d0
c      p(8)  =  0.d0
c      p(9)  =  e(3)*det/24.d0
c      p(10) =  0.d0
c      p(11) =  0.d0
c      p(12) =  e(3)*det/24.d0           
      return
c ======================================================================
c
c ... 
c
c ......................................................................
  500 continue
      return
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT09: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
c      subroutine elmt08_m(e,x,u,v,p,s,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT08: Elemento cunha                                           *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - deslocamento nodais do elemento                   *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 = leitura das constantes fisicas                       *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 = forcas internas ku                                   *
c *           4 = tensoes nodais                                       *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2,3) ou tensoes (isw = 4)           *
c *                                                                    *
c **********************************************************************
c      implicit none
c      common /gauss/ pg, wg
c      integer ndm,nst,nel,isw,nen
c      integer nint,lz,i,j,i1,i2,i3,j1,j2,j3
c      real*8 e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst)
c      real*8 h(6),hx(6),hy(6),hz(6),xj(3,3),xji(3,3)
c      real*8 ri,si,ti,wt,det
c      real*8 ym,cp,a,b,c
c      real*8 pg(10,10),wg(10,10)
c      real*8 eps(6)
c      data   nint/2/,nen/6/
cc ......................................................................
c      go to (100,200,200,400) isw
cc ======================================================================
cc
cc.... Input material properties:
cc
cc ......................................................................
c  100 continue
cc     e(1) = modulo de elasticidade
cc     e(2) = coeficiente de Poisson
c      return
cc ======================================================================
cc
cc.... Matrix K:
cc ......................................................................
c  200 continue
c
c ... Propriedades do material:
c ......................................................................
c      ym    = e(1)
c      cp    = e(2)
cc ......................................................................
cc ... Matriz constitutiva:
cc ......................................................................
c      a = ym*(1.-cp)/((1.+cp)*(1.-2.d0*cp))
c      b = cp/(1.-cp)
c      c = (1.-2.d0*cp)/(2.d0*(1.-cp))
c      if (isw .eq. 3) goto 300
c ......................................................................
c ... Matriz de rigidez:
c ......................................................................
c      do 210 i = 1, nst
c      do 210 j = 1, nst
c         s(i,j) = 0.d0
c  210 continue
cc ......................................................................
c      do 240 lz = 1, nint
c      ri = 1.d0/3.d0
c      si = 1.d0/3.d0
c      ti = pg(lz,nint)
c      call sfwedge(h,hx,hy,hz,ri,si,ti,.false.,.true.)
c      call jacob3d(hx,hy,hz,xj,xji,x,det,nen,nel)
c      wt = wg(lz,nint)*0.5d0*det*a
c      do 230 i = 1, 6
c      i1 = (i-1)*3+1
c      i2 = i1 + 1
c      i3 = i2 + 1
c      do 220 j = 1, 6
c        j1 = (j-1)*3+1
c        j2 = j1 + 1
c        j3 = j2 + 1
c        s(i1,j1) = s(i1,j1) + (  hx(i)*hx(j) + c*hy(i)*hy(j) +
c     .                                             c*hz(i)*hz(j))*wt
c        s(i1,j2) = s(i1,j2) + (b*hx(i)*hy(j) + c*hy(i)*hx(j))*wt
c        s(i1,j3) = s(i1,j3) + (b*hx(i)*hz(j) + c*hz(i)*hx(j))*wt
c        s(i2,j1) = s(i2,j1) + (b*hy(i)*hx(j) + c*hx(i)*hy(j))*wt
c        s(i2,j2) = s(i2,j2) + (  hy(i)*hy(j) + c*hx(i)*hx(j) +
c     .                                             c*hz(i)*hz(j))*wt
c        s(i2,j3) = s(i2,j3) + (b*hy(i)*hz(j) + c*hz(i)*hy(j))*wt
c        s(i3,j1) = s(i3,j1) + (b*hz(i)*hx(j) + c*hx(i)*hz(j))*wt
c        s(i3,j2) = s(i3,j2) + (b*hz(i)*hy(j) + c*hy(i)*hz(j))*wt
c        s(i3,j3) = s(i3,j3) + (  hz(i)*hz(j) + c*hy(i)*hy(j) +
c     .                                             c*hx(i)*hx(j))*wt
c  220    continue
c  230 continue
c  240 continue
cc ======================================================================
c
c ... Forcas internas:
c ......................................................................
c      call lku(s,u,p,nst)
c      return
cc ======================================================================
cc
cc ... Tensoes nodais:        
cc
cc ......................................................................
c  300 continue
c      call sfwedge(h,hx,hy,hz,1.d0/3.d0,1.d0/3.d0,1.d0,.false.,.true.)
c      call jacob3d(hx,hy,hz,xj,xji,x,det,nen,nel)
c      call deform3d(hx,hy,hz,u,eps,nen)
c      call stress3d(a,b,c,eps,p(1))
c      p(7)  = p(1)
c      p(8)  = p(2)
c      p(9)  = p(3)
c      p(10) = p(4)
c      p(11) = p(5)
c      p(12) = p(6)
c      p(13) = p(1)
c      p(14) = p(2)
c      p(15) = p(3)
c      p(16) = p(4)
c      p(17) = p(5)
c      p(18) = p(6)
c      call sfwedge(h,hx,hy,hz,1.d0/3.d0,1.d0/3.d0,-1.d0,.false.,.true.)
c      call jacob3d(hx,hy,hz,xj,xji,x,det,nen,nel)
c      call deform3d(hx,hy,hz,u,eps,nen)
c      call stress3d(a,b,c,eps,p(19))   
c      p(25) = p(19)
c      p(26) = p(20)
c      p(27) = p(21)
c      p(28) = p(22)
c      p(29) = p(23)
c      p(30) = p(24)
c      p(31) = p(19)
c      p(32) = p(20)
c      p(33) = p(21)
c      p(34) = p(22)
c      p(35) = p(23)
c      p(36) = p(24)       
c      return
cc ======================================================================
cc
cc ... Matriz de massa:
cc
c ......................................................................
c  400 continue
cc      print*,'WARNING : Cargas de Volume nao disponiveis !'
c      return
cc ======================================================================
cc
cc ... Matriz geometrica:
cc
cc ......................................................................
c  600 continue
c      print*,'Matriz geometrica nao disponivel p/ o elemento tipo 5 !'
c      stop
cc ......................................................................
c 1000 continue
c      print*, '*** Subrotina ELMT08: determinante nulo ou negativo do el
c     .emento ',nel
c      stop
c      end
