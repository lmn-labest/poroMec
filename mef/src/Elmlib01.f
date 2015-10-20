      subroutine elmt01_t(e,iq,xl,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT01: Elemento unidimensional (difusao)                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nem)- coordenadas nodais locais                         *
c *     u(nst)     - valores prescritos                                *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     p   = nodal forces                                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'            
      integer iq(*),ndm,nst,nel,isw,i,j,k,ma,nlit
      real*8 e(*),xl(ndm,*),u(*),v(*),p(*),s(nst,*),difus,capac,tm
      real*8 hx(2),det,wt,m11,m12,a(3),ea,clat,atil,hi(3),dhi,ksi
c ......................................................................
      goto (100,200,200,500) isw
c ......................................................................
c
c ... Propriedades do material:
c
c ......................................................................
  100 continue
c     e(1) = difusao
c     e(2) = coeficiente da derivada no tempo
      return
c ......................................................................
c
c.... Matrix K:
c
c ......................................................................
  200 continue
      difus = e(1)
      capac = e(2)
      tm = (u(1)+u(2))/2
      call updprop(difus,tm,ma,6)
      call updprop(capac,tm,ma,7)
c     Propriedades Termo - quimico : 
      ea    = e(3)
      clat  = e(4)
      atil  = e(5)
      if (nlit .eq. 1)    hi(1) = hi(2)
      call calc_hidr(hi(1),hi(2),dhi,hi(3),tm,atil,ea,ma,nlit)
c
      if (isw .eq. 3) goto 300
c
      det = 0.d0
      do 210 i = 1, ndm
         det = det + (xl(i,2)-xl(i,1))*(xl(i,2)-xl(i,1))
  210 continue
       det = dsqrt(det)
       if (det .le. 0.d0) goto 1000
c
c ... Derivadas das funcoes de interpolacao:
c
      hx(1) = -1.d0/det
      hx(2) =  1.d0/det
c
c ... Matriz de rigidez:                     
c
      wt = det*difus
      s(1,1) = (hx(1)*hx(1))*wt
      s(1,2) = (hx(1)*hx(2))*wt
      s(2,2) = (hx(2)*hx(2))*wt
      s(2,1) =  s(1,2)
c
c ... Matriz de massa consistente:
c
      m11 = (1.d0/3.d0)*det*capac
      m12 = m11*0.5d0
c
c ... Forcas internas:
c      
      p(1) = s(1,1)*u(1)+s(1,2)*u(2) + m11*v(1)+m12*v(2) 
     . - dhi*clat*det*0.5d0
      p(2) = s(1,2)*u(1)+s(2,2)*u(2) + m12*v(1)+m11*v(2) 
     . - dhi*clat*det*0.5d0
c
c ... Matriz de rigidez efetiva:                     
c
      s(1,1) = m11 + s(1,1)*alfa*dt
      s(1,2) = m12 + s(1,2)*alfa*dt
      s(2,2) = m11 + s(2,2)*alfa*dt
      s(2,1) =  s(1,2)            
      return
  300 continue
c ......................................................................
c
c.... Fluxos nos vertices:
c
c ......................................................................  
      det = 0.d0
      do 310 i = 1, ndm
         a(i) = xl(i,2)-xl(i,1)
         det = det + a(i)*a(i)
  310 continue  
       det = dsqrt(det)
       if (det .le. 0.d0) goto 1000
c
c ... Derivadas das funcoes de interpolacao:
c
      hx(1) = -1.d0/det
      hx(2) =  1.d0/det
c
c ... Fluxo:
c
      m11 = -(hx(1)*u(1)+hx(2)*u(2))*difus
      do 330 i = 1, ndm
         a(i) = a(i)/det
         do 320 j = 1, 2
            k = (j-1)*ndm+i
            p(k) = m11*a(i)
  320    continue
  330 continue
      return
c ......................................................................
c
c.... Cargas de volume:
c
c ......................................................................
  500 continue
      if(iq(1) .eq. 0) return
      call tload(iq(1),t,u,m11)
      det = 0.d0
      do 510 i = 1, ndm
         a(i) = xl(i,2)-xl(i,1)
         det = det + a(i)*a(i)
  510 continue  
      det = dsqrt(det)
      if (det .le. 0.d0) goto 1000
      p(1) = -det*.5d0*m11
      p(2) = -det*.5d0*m11
      return      
c ......................................................................      
 1000 continue
      print*, '*** Subrotina ELMT01: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end            
      subroutine elmt01c_t(e,iq,xl,u,v,p,s,a,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT01: Elemento unidimensional (conveccao-difusao, SUPG)        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nem)- coordenadas nodais locais                         *
c *     u(nst)     - valores prescritos                                *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     p   = nodal forces                                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'            
      integer iq(*),ndm,nst,nel,isw,i,j,k
      real*8 e(*),xl(ndm,*),u(*),v(*),p(*),s(nst,*),a(ndm,*),q
      real*8 det,wt,m11,m12,m21,m22,vx,vx1,vx2,vmod,peclet,alpha,xk,kr
c ......................................................................
      goto (100,200,300,500) isw
c ......................................................................
c
c ... Propriedades do material:
c
c ......................................................................
  100 continue
c     e(1) = difusao
c     e(2) = coeficiente da derivada no tempo
      return
c ......................................................................
c
c.... Matrix K:
c
c ......................................................................
  200 continue
      det = 0.d0
      do 210 i = 1, ndm
         det = det + (xl(i,2)-xl(i,1))*(xl(i,2)-xl(i,1))
  210 continue
       det = dsqrt(det)
       if (det .le. 0.d0) goto 1000
c
c     Velocidade media:
c
      vx1 = a(1,1)
      vx2 = a(1,2)
      vx = (vx1+vx2)/2.d0
c
c ... difusao artificial (xk):
c
      vmod = dabs(vx)
      xk    = 0.d0
      alpha = 1.d0
      if(e(1) .gt. 0.d0) then
         peclet= (det*vmod)/(2.d0*e(1))
         alpha = min((1.d0/3.d0)*peclet,1.d0)
      endif
      if(vmod .gt. 0.d0) xk = alpha*det*0.5d0/vmod
c
c ... Matriz de rigidez:                     
c
      vx = xk * (vx1*vx1 + vx1*vx2 + vx2*vx2)/(3.d0*det)
      wt = e(1)/det
      s(1,1) =  wt - (2.d0*vx1+vx2)/6.d0 + vx
      s(1,2) = -wt + (2.d0*vx1+vx2)/6.d0 - vx
      s(2,1) = -wt - (2.d0*vx2+vx1)/6.d0 - vx
      s(2,2) =  wt + (2.d0*vx2+vx1)/6.d0 + vx
c
c ... Matriz de massa consistente:
c
      m11 = (1.d0/3.d0)*det - xk*(2.d0*vx1+vx2)/6.d0
      m12 = (1.d0/6.d0)*det - xk*(2.d0*vx2+vx1)/6.d0
      m21 = (1.d0/6.d0)*det + xk*(2.d0*vx1+vx2)/6.d0            
      m22 = (1.d0/3.d0)*det + xk*(2.d0*vx2+vx1)/6.d0
c
c    Adiciona termos de decaimento:
c
      kr = 0.d0
      s(1,1) = s(1,1) + kr*m11
      s(1,2) = s(1,2) + kr*m12
      s(2,1) = s(2,1) + kr*m21    
      s(2,2) = s(2,2) + kr*m22
c
c ... Matriz de massa consistente x e(2):
c
      m11 = m11*e(2)
      m12 = m12*e(2)
      m21 = m21*e(2)          
      m22 = m22*e(2)
c
c ... Forcas internas:
c      
      p(1) = s(1,1)*u(1)+s(1,2)*u(2) + m11*v(1)+m12*v(2)
      p(2) = s(2,1)*u(1)+s(2,2)*u(2) + m21*v(1)+m22*v(2)
c
c ... Matriz de rigidez efetiva:                     
c
      s(1,1) = m11 + s(1,1)*alfa*dt
      s(1,2) = m12 + s(1,2)*alfa*dt
      s(2,1) = m21 + s(2,1)*alfa*dt            
      s(2,2) = m22 + s(2,2)*alfa*dt
      return
c ......................................................................
c
c.... Fluxos nos vertices:
c
c ......................................................................  
  300 continue
      det = 0.d0
      do 310 i = 1, ndm
         det = det + (xl(i,2)-xl(i,1))*(xl(i,2)-xl(i,1))
  310 continue  
      det = dsqrt(det)
      if (det .le. 0.d0) goto 1000
c ... Fluxo:
      m11 = (u(2)-u(1))*e(1)/det
      do 330 i = 1, ndm
      do 320 j = 1, 2
         k = (j-1)*ndm+i
         p(k) = -m11*(xl(i,2)-xl(i,1))/det
  320 continue
  330 continue
      return
c ......................................................................
c
c.... Cargas de volume:
c
c ......................................................................
  500 continue
      if(iq(1) .eq. 0) return
      call tload(iq(1),t,u,q)
      det = 0.d0
      do 510 i = 1, ndm
         det = det + (xl(i,2)-xl(i,1))*(xl(i,2)-xl(i,1))
  510 continue  
      det = dsqrt(det)
      if (det .le. 0.d0) goto 1000
c ......................................................................      
c     Velocidade media:
      vx1 = a(1,1)
      vx2 = a(1,2)
      vx = (vx1+vx2)/2.d0
c ... difusao artificial:
      vmod = dabs(vx)
      xk    = 0.d0
      alpha = 1.d0
      if(e(1) .gt. 0.d0) then
         peclet = (det*vmod)/(2.d0*e(1))
         alpha  =  min((1.d0/3.d0)*peclet,1.d0)
      end if
      if(vmod .gt. 0.d0) xk = alpha*det*0.5d0/vmod 
c ......................................................................                 
      p(1) = -q*(det*.5d0 - xk*vx)
      p(2) = -q*(det*.5d0 + xk*vx)
      return      
c ......................................................................      
 1000 continue
      print*, '*** Subrotina ELMT01: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
      subroutine elmt03_t(e,iq,xl,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c **********************************************************************
c *                                                                    *
c *   elmt03: elemento triangular linear (difusao)                     *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nem)- coordenadas nodais locais                         *
c *     u(nst)     - valores prescritos                                *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     p   = nodal forces                                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'                  
      integer iq(*),ndm,nst,nel,isw,i,k1,k2,k3,ma,nlit
      real*8 e(*),xl(ndm,*),u(*),v(*),p(*),s(nst,*),tm,difus,capac
      real*8 hx(2),hy(2),det,wt,xj11,xj12,xj21,xj22,sm1,sm2,r(3,3)
      real*8 hi(3),ea,clat,atil,dhi,a1,a2,b1,c1
c ......................................................................
      goto (100,200,200,500) isw
c ......................................................................
c
c ... Propriedades do material:
c
c ......................................................................
  100 continue
c     e(1) = difusao
c     e(2) = coeficiente da derivada no tempo
      return
c ......................................................................
c
c.... Matrix K:
c
c ......................................................................
  200 continue
      difus = e(1)
      capac = e(2)
      tm = (u(1)+u(2)+u(3))/3.d0
      call updprop(difus,tm,ma,6)
      call updprop(capac,tm,ma,7)
c 
c     Propriedades Termo - quimico : 
      ea    = e(3)
      atil  = e(4)
      clat  = e(5)
      if (nlit .eq. 1)    then
         hi(1) = hi(2)
         hi(3) = tm
      endif
      call calc_hidr(hi(1),hi(2),dhi,hi(3),tm,atil,ea,ma,nlit)
c     
      if (isw .eq. 3) goto 300
c    
      if(ndm .eq. 3) then
         call rotmatrix(xl,r)
         call rotate(xl(1,1),r,xl(1,1),.false.)
         call rotate(xl(1,2),r,xl(1,2),.false.)
         call rotate(xl(1,3),r,xl(1,3),.false.)
      endif  
c
c ... Matriz Jacobiana:
c
      xj11 = xl(1,1)-xl(1,3)
      xj12 = xl(2,1)-xl(2,3)
      xj21 = xl(1,2)-xl(1,3)
      xj22 = xl(2,2)-xl(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) goto 1000
c
c ... Derivadas das funcoes de interpolacao:
c
      hx(1) =  xj22/det
      hx(2) = -xj12/det
      hy(1) = -xj21/det
      hy(2) =  xj11/det
c
c ... Matriz de rigidez:                     
c
      wt = 0.5d0*det*difus
      s(1,1) = (hx(1)*hx(1)+hy(1)*hy(1))*wt
      s(1,2) = (hx(1)*hx(2)+hy(1)*hy(2))*wt
      s(2,2) = (hx(2)*hx(2)+hy(2)*hy(2))*wt
      s(1,3) = -s(1,1)-s(1,2)
      s(2,3) = -s(1,2)-s(2,2)
      s(3,3) = -s(1,3)-s(2,3)
      s(2,1) =  s(1,2)
      s(3,1) =  s(1,3)
      s(3,2) =  s(2,3)
c
c ... Matriz de massa:
c
      sm1 = (1.d0/12.d0)*det*capac
      sm2 = sm1*0.5d0
c
c ... Forcas internas:
c      
      p(1)= s(1,1)*u(1)+s(1,2)*u(2)+s(1,3)*u(3)+sm1*v(1)+sm2*(v(2)+v(3))
     .       - dhi*clat*capac*det/6.d0
      p(2)= s(1,2)*u(1)+s(2,2)*u(2)+s(2,3)*u(3)+sm1*v(2)+sm2*(v(1)+v(3))
     .       - dhi*clat*capac*det/6.d0
      p(3)= s(1,3)*u(1)+s(2,3)*u(2)+s(3,3)*u(3)+sm1*v(3)+sm2*(v(1)+v(2))
     .       - dhi*clat*capac*det/6.d0
c
c ... Matriz de rigidez efetiva:                     
c
      s(1,1) = sm1 + s(1,1)*alfa*dt
      s(1,2) = sm2 + s(1,2)*alfa*dt
      s(1,3) = sm2 + s(1,3)*alfa*dt
      s(2,2) = sm1 + s(2,2)*alfa*dt
      s(2,3) = sm2 + s(2,3)*alfa*dt
      s(3,3) = sm1 + s(3,3)*alfa*dt
      s(2,1) = s(1,2)
      s(3,1) = s(1,3)
      s(3,2) = s(2,3)    
      return
  300 continue
c ......................................................................
c
c.... Fluxos nos vertices:
c
c ......................................................................
c
c ... sistema local de coordenadas, se ndm = 3:
c
      if(ndm .eq. 3) then
         call rotmatrix(xl,r)
         call rotate(xl(1,1),r,xl(1,1),.false.)
         call rotate(xl(1,2),r,xl(1,2),.false.)
         call rotate(xl(1,3),r,xl(1,3),.false.)
      endif
c ......................................................................        
c
c ... Derivadas da solucao:
c
c ... Matriz Jacobiana:
c
      xj11 = xl(1,1)-xl(1,3)
      xj12 = xl(2,1)-xl(2,3)
      xj21 = xl(1,2)-xl(1,3)
      xj22 = xl(2,2)-xl(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) goto 1000
c
c ... Derivadas das funcoes de interpolacao:
c
      hx(1) =  xj22/det
      hx(2) = -xj12/det
      hy(1) = -xj21/det
      hy(2) =  xj11/det
c ......................................................................              
c
c ... Fluxo:
c
      do 310 i = 1, 3
         k1 = (i-1)*ndm+1
         k2 = k1 + 1
         p(k1) = (-hx(1)*u(1)-hx(2)*u(2)+(hx(1)+hx(2))*u(3))*difus
         p(k2) = (-hy(1)*u(1)-hy(2)*u(2)+(hy(1)+hy(2))*u(3))*difus
         if(ndm .eq. 3) then
            k3 = k2 + 1
            p(k3) = 0.d0
            call rotate(p(k1),r,p(k1),.true.)
         endif 
  310 continue        
      return
c ......................................................................  
c
c ... Cargas distribuidas no volume e no contorno:
c            
c ......................................................................
  500 continue
      p(1) = 0.d0
      p(2) = 0.d0
      p(3) = 0.d0
      if(iq(1) .eq. 0)   return
      call tload(iq(1),t,u,sm1)
      xj11 = xl(1,1)-xl(1,3)
      xj12 = xl(2,1)-xl(2,3)
      xj21 = xl(1,2)-xl(1,3)
      xj22 = xl(2,2)-xl(2,3)
      det  = xj11*xj22-xj12*xj21
      sm1 = sm1*det/6.d0      
      p(1) = -sm1
      p(2) = -sm1
      p(3) = -sm1
      return                  
c ......................................................................      
 1000 continue
      print*, '*** Subrotina ELMT03: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
      subroutine elmt03c_t(e,iq,xl,u,v,p,s,vn,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT03: Elemento triangular linear, conveccao-difusao (SUPG)     *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nen)- coordenadas nodais locais                         *
c *     u(nst)     - solucao nodal                                     *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     v   - derivadas temporais de u                                 *
c *     vn  - velocidades locais                                       *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      integer iq(*),ndm,nst,nel,isw
      real*8 e(10),xl(ndm,*),u(*),v(*),p(*),s(nst,*),vn(ndm,*)
      real*8 hx(3),hy(3),det,area,wt,xj11,xj12,xj21,xj22,area_3,area_6
      real*8 hpg1,hpg2,hpg3,xk,s11,s12,s13,s21,s22,s23,s31,s32,s33
      real*8 sm11,sm12,sm13,sm21,sm22,sm23,sm31,sm32,sm33,sm1,sm2
      real*8 hc,peclet,alpha,a1,kd,ks,kr,div
      real*8 vx,vy,vmod,vx1,vy1,vx2,vy2,vx3,vy3,conv11,conv22,conv33
      real*8 conv1,conv2,conv3,conv12,conv13,conv21,conv23,conv31,conv32
c ......................................................................
      goto (100,200,300,500) isw
c ......................................................................
c
c.... Input material properties:
c
c ......................................................................
  100 continue
c     e(1) = coeficiente de difusao
c     e(2) =                         
c     e(3) =                          
c     e(4) =                                  
      return
c ......................................................................
c
c....Matriz de Condutividade   
c
c ......................................................................
  200 continue
c
c     Matriz Jacobiana:
c
      xj11 = xl(1,1)-xl(1,3)
      xj12 = xl(2,1)-xl(2,3)
      xj21 = xl(1,2)-xl(1,3)
      xj22 = xl(2,2)-xl(2,3)
      det  = xj11*xj22-xj12*xj21
      area = det*0.5d0
      if (det .le. 0.d0) goto 1000
      area_3 = area/3.d0
      area_6 = area/6.d0 
c
c     Derivadas das funcoes de interpolacao:
c
      hx(1) =  xj22/det
      hx(2) = -xj12/det
      hx(3) = (xj12-xj22)/det
      hy(1) = -xj21/det
      hy(2) =  xj11/det
      hy(3) = (xj21-xj11)/det
c
c     Velocidade media:
c
      vx = (vn(1,1)+vn(1,2)+vn(1,3))/3.d0
      vy = (vn(2,1)+vn(2,2)+vn(2,3))/3.d0
c
c     Velocidades nos lados:
c
      vx1 = (vn(1,1)+vn(1,2))*0.5d0
      vy1 = (vn(2,1)+vn(2,2))*0.5d0
      vx2 = (vn(1,2)+vn(1,3))*0.5d0
      vy2 = (vn(2,2)+vn(2,3))*0.5d0
      vx3 = (vn(1,3)+vn(1,1))*0.5d0
      vy3 = (vn(2,3)+vn(2,1))*0.5d0
c
c ... Divergente do campo de velocidades:      
c
      div = vn(1,1)*hx(1)+vn(1,2)*hx(2)+vn(1,3)*hx(3)+
     .      vn(2,1)*hy(1)+vn(2,2)*hy(2)+vn(2,3)*hy(3)
c
c ... Termos de conveccao (Galerkin):
c
      conv11 = (vx1*hx(1)+vy1*hy(1)+vx3*hx(1)+vy3*hy(1))*area_6
      conv22 = (vx1*hx(2)+vy1*hy(2)+vx2*hx(2)+vy2*hy(2))*area_6
      conv33 = (vx2*hx(3)+vy2*hy(3)+vx3*hx(3)+vy3*hy(3))*area_6
      conv12 = (vx1*hx(2)+vy1*hy(2)+vx3*hx(2)+vy3*hy(2))*area_6
      conv13 = (vx1*hx(3)+vy1*hy(3)+vx3*hx(3)+vy3*hy(3))*area_6
      conv21 = (vx1*hx(1)+vy1*hy(1)+vx2*hx(1)+vy2*hy(1))*area_6
      conv23 = (vx1*hx(3)+vy1*hy(3)+vx2*hx(3)+vy2*hy(3))*area_6
      conv31 = (vx2*hx(1)+vy2*hy(1)+vx3*hx(1)+vy3*hy(1))*area_6
      conv32 = (vx2*hx(2)+vy2*hy(2)+vx3*hx(2)+vy3*hy(2))*area_6
c
      vmod=dsqrt(vx*vx+vy*vy)
c
c ... tamanho caracteristico do elememto (raiz da area):
c
      hc = dsqrt(area)
c
      if(e(1) .gt. 0.d0) then
           peclet= (hc*vmod)/(2.d0*e(1))
           alpha = min((1.d0/3.d0)*peclet,1.0d0)
      else
           alpha=1.0d0
      end if
c
c ... difusao artificial:
c
      xk = 0.d0
      if(vmod .gt. 0.d0) xk = alpha*hc*0.5d0/vmod
c
c ... funcoes peso de Petrov-Galerkin:
c
      conv1 = vx*hx(1)+vy*hy(1)
      conv2 = vx*hx(2)+vy*hy(2)
      conv3 = vx*hx(3)+vy*hy(3)
      hpg1  = xk*conv1
      hpg2  = xk*conv2
      hpg3  = xk*conv3
c
c ... Concentracao e profundidade medias:
c
c      um = (u(1)+u(2)+u(3))/3.d0
c      hm = (xl(3,1)+xl(3,2)+xl(3,3))/3.d0
c
c ... coeficientes de decaimento:
c
c      kd = e(2)
c      ks = e(3)/hm
c      ks = e(3)
c      kr = kd+ks+div
      kr = 0.d0
c     
c     Matriz de rigidez:                     
c
c ... Termos de difusao:
c
      wt  = area*e(1)
      s11 =(hx(1)*hx(1)+hy(1)*hy(1))*wt
      s12 =(hx(1)*hx(2)+hy(1)*hy(2))*wt
      s22 =(hx(2)*hx(2)+hy(2)*hy(2))*wt
      s13 = -s11-s12
      s23 = -s12-s22
      s33 = -s13-s23
      s21 = s12
      s31 = s13
      s32 = s23 
c
c     Matriz de massa:                     
c
      sm1  = area_6        
      sm2  = sm1*0.5d0        
      sm11 = (sm1+hpg1*area_3)
      sm22 = (sm1+hpg2*area_3)
      sm33 = (sm1+hpg3*area_3)
      sm12 = (sm2+hpg1*area_3)
      sm21 = (sm2+hpg2*area_3)
      sm31 = (sm2+hpg3*area_3)
      sm13 = sm12               
      sm23 = sm21
      sm32 = sm31
c
c ... termos difusao e conveccao (Galerkin) mais 
c     termos de conveccao (Petrov/Galerkin):
c
      s11 = s11 + hpg1*conv1*area + kr*sm11 + conv11
      s22 = s22 + hpg2*conv2*area + kr*sm22 + conv22       
      s33 = s33 + hpg3*conv3*area + kr*sm33 + conv33
      s21 = s12 + hpg2*conv1*area + kr*sm21 + conv21
      s31 = s13 + hpg3*conv1*area + kr*sm31 + conv31
      s32 = s23 + hpg3*conv2*area + kr*sm32 + conv32
      s12 = s12 + hpg1*conv2*area + kr*sm12 + conv12
      s13 = s13 + hpg1*conv3*area + kr*sm13 + conv13
      s23 = s23 + hpg2*conv3*area + kr*sm23 + conv23
c
c     Matriz de massa * porosidade:
c      
      sm11 = sm11*e(2)
      sm22 = sm22*e(2)
      sm33 = sm33*e(2)
      sm12 = sm12*e(2)
      sm21 = sm21*e(2)
      sm31 = sm31*e(2)
      sm13 = sm13*e(2)      
      sm32 = sm32*e(2)
      sm23 = sm23*e(2) 
c           
c     Residuo: R = ku + Mv
c
      p(1) = s11*u(1)+s12*u(2)+s13*u(3)+sm11*v(1)+sm12*(v(2)+v(3))
      p(2) = s21*u(1)+s22*u(2)+s23*u(3)+sm22*v(2)+sm21*(v(1)+v(3))
      p(3) = s31*u(1)+s32*u(2)+s33*u(3)+sm33*v(3)+sm31*(v(1)+v(2))
c
c     Matriz efetiva:  M + alfa.dt.K
c
      a1 = alfa*dt
      s(1,1) = s11*a1 + sm11
      s(2,2) = s22*a1 + sm22
      s(3,3) = s33*a1 + sm33
      s(2,1) = s21*a1 + sm21
      s(3,1) = s31*a1 + sm31
      s(3,2) = s32*a1 + sm32
      s(1,2) = s12*a1 + sm12
      s(1,3) = s13*a1 + sm13
      s(2,3) = s23*a1 + sm23
      return
c.......................................................................
  300 continue
      return
c ......................................................................  
c
c ... Cargas distribuidas no volume e no contorno:
c            
c ......................................................................
  500 continue
      if(iq(1) .eq. 0) return
      call tload(iq(1),t,u,a1)
      xj11 = xl(1,1)-xl(1,3)
      xj12 = xl(2,1)-xl(2,3)
      xj21 = xl(1,2)-xl(1,3)
      xj22 = xl(2,2)-xl(2,3)
      det  = xj11*xj22-xj12*xj21
      vx = (vn(1,1)+vn(1,2)+vn(1,3))/3.d0
      vy = (vn(2,1)+vn(2,2)+vn(2,3))/3.d0            
      vmod=dsqrt(vx*vx+vy*vy)
      area = det*.5d0
      hc = dsqrt(area)
      if(e(1) .gt. 0.d0) then
           peclet= (hc*vmod)/(2.d0*e(1))
           alpha = min((1.d0/3.d0)*peclet,1.0d0)
      else
           alpha=1.0d0
      end if
      xk = 0.d0
      if(vmod .gt. 0.d0) xk = alpha*hc*0.5d0/vmod 
      hx(1) =  xj22/det
      hx(2) = -xj12/det
      hx(3) = (xj12-xj22)/det
      hy(1) = -xj21/det
      hy(2) =  xj11/det
      hy(3) = (xj21-xj11)/det
      hpg1 = (vx*hx(1)+vy*hy(1))*xk
      hpg2 = (vx*hx(2)+vy*hy(2))*xk
      hpg3 = (vx*hx(3)+vy*hy(3))*xk            
      p(1) = -a1*(1.d0/3.d0 + hpg1)*area
      p(2) = -a1*(1.d0/3.d0 + hpg2)*area
      p(3) = -a1*(1.d0/3.d0 + hpg3)*area           
      return            
c ......................................................................           
 1000 continue
      print*, '*** Subrotina ELMT03: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end      
      subroutine elmt03b_t(e,iq,xl,u,v,p,s,vn,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT03: Elemento triangular linear, conveccao-difusao (CG)       *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nen)- coordenadas nodais locais                         *
c *     u(nst)     - solucao nodal                                     *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     v   - derivadas temporais de u                                 *
c *     vn  - velocidades locais                                       *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *     nin - arquivo de entrada                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      integer iq(*),ndm,nst,nel,isw
      real*8 e(10),xl(ndm,*),u(*),v(*),p(*),s(nst,*),vn(ndm,*)
      real*8 hx(3),hy(3),det,area,wt,xj11,xj12,xj21,xj22,area_3,area_6
      real*8 hpg1,hpg2,hpg3,xk,s11,s12,s13,s21,s22,s23,s31,s32,s33
      real*8 a1,kd,ks,kr,div,sm1,sm2
      real*8 vx,vy,vx1,vy1,vx2,vy2,vx3,vy3,conv11,conv22,conv33
      real*8 conv12,conv13,conv21,conv23,conv31,conv32
c ......................................................................
      goto (100,200,300,500) isw
c ......................................................................
c
c.... Input material properties:
c
c ......................................................................
  100 continue
c     e(1) = coeficiente de difusao
c     e(2) = porosidade              
c     e(3) =                          
c     e(4) =                                  
      return
c ......................................................................
c
c....Matriz de Condutividade   
c
c ......................................................................
  200 continue
c
c     Matriz Jacobiana:
c
      xj11 = xl(1,1)-xl(1,3)
      xj12 = xl(2,1)-xl(2,3)
      xj21 = xl(1,2)-xl(1,3)
      xj22 = xl(2,2)-xl(2,3)
      det  = xj11*xj22-xj12*xj21
      area = det*0.5d0
      if (det .le. 0.d0) goto 1000
      area_3 = area/3.d0
      area_6 = area/6.d0 
c
c     Derivadas das funcoes de interpolacao:
c
      hx(1) =  xj22/det
      hx(2) = -xj12/det
      hx(3) = (xj12-xj22)/det
      hy(1) = -xj21/det
      hy(2) =  xj11/det
      hy(3) = (xj21-xj11)/det
c
c     Velocidade media:
c
      vx = (vn(1,1)+vn(1,2)+vn(1,3))/3.d0
      vy = (vn(2,1)+vn(2,2)+vn(2,3))/3.d0
c
c     Velocidades nos lados:
c
      vx1 = (vn(1,1)+vn(1,2))*0.5d0
      vy1 = (vn(2,1)+vn(2,2))*0.5d0
      vx2 = (vn(1,2)+vn(1,3))*0.5d0
      vy2 = (vn(2,2)+vn(2,3))*0.5d0
      vx3 = (vn(1,3)+vn(1,1))*0.5d0
      vy3 = (vn(2,3)+vn(2,1))*0.5d0
c
c ... Divergente do campo de velocidades:      
c
      div = vn(1,1)*hx(1)+vn(1,2)*hx(2)+vn(1,3)*hx(3)+
     .      vn(2,1)*hy(1)+vn(2,2)*hy(2)+vn(2,3)*hy(3)
c
c ... Termos de conveccao (Galerkin):
c
      conv11 = (vx1*hx(1)+vy1*hy(1)+vx3*hx(1)+vy3*hy(1))*area_6
      conv22 = (vx1*hx(2)+vy1*hy(2)+vx2*hx(2)+vy2*hy(2))*area_6
      conv33 = (vx2*hx(3)+vy2*hy(3)+vx3*hx(3)+vy3*hy(3))*area_6
      conv12 = (vx1*hx(2)+vy1*hy(2)+vx3*hx(2)+vy3*hy(2))*area_6
      conv13 = (vx1*hx(3)+vy1*hy(3)+vx3*hx(3)+vy3*hy(3))*area_6
      conv21 = (vx1*hx(1)+vy1*hy(1)+vx2*hx(1)+vy2*hy(1))*area_6
      conv23 = (vx1*hx(3)+vy1*hy(3)+vx2*hx(3)+vy2*hy(3))*area_6
      conv31 = (vx2*hx(1)+vy2*hy(1)+vx3*hx(1)+vy3*hy(1))*area_6
      conv32 = (vx2*hx(2)+vy2*hy(2)+vx3*hx(2)+vy3*hy(2))*area_6
c
c ... difusao artificial:
c
      xk = dt*.5d0
c
c ... termos de estabilizacao:
c
      hpg1 = vx*hx(1)+vy*hy(1)
      hpg2 = vx*hx(2)+vy*hy(2)
      hpg3 = vx*hx(3)+vy*hy(3)
c
c ... Concentracao e profundidade medias:
c
c      um = (u(1)+u(2)+u(3))/3.d0
c      hm = (xl(3,1)+xl(3,2)+xl(3,3))/3.d0
c
c ... coeficientes de decaimento:
c
c      kd = e(2)
c      ks = e(3)/hm
c      ks = e(3)
c      kr = kd+ks+div
      kr = 0.d0
c     
c     Matriz de rigidez:                     
c
c ... Termos de difusao:
c
      wt  = area*e(1)
      s11 =(hx(1)*hx(1)+hy(1)*hy(1))*wt
      s12 =(hx(1)*hx(2)+hy(1)*hy(2))*wt
      s22 =(hx(2)*hx(2)+hy(2)*hy(2))*wt
      s13 = -s11-s12
      s23 = -s12-s22
      s33 = -s13-s23
      s21 = s12
      s31 = s13
      s32 = s23 
c
c     Matriz de massa:                     
c
      sm1  = area_6        
      sm2  = sm1*0.5d0
c
c ... termos difusao e conveccao (Galerkin) mais 
c     termos de estabilizacao:
c
      s11 = s11 + (conv11 + hpg1*hpg1*area*xk)*(1.d0+kr) + kr*sm1 
      s22 = s22 + (conv22 + hpg2*hpg2*area*xk)*(1.d0+kr) + kr*sm1
      s33 = s33 + (conv33 + hpg3*hpg3*area*xk)*(1.d0+kr) + kr*sm1 
      s21 = s12 + (conv21 + hpg2*hpg1*area*xk)*(1.d0+kr) + kr*sm2 
      s31 = s13 + (conv31 + hpg3*hpg1*area*xk)*(1.d0+kr) + kr*sm2 
      s32 = s23 + (conv32 + hpg3*hpg2*area*xk)*(1.d0+kr) + kr*sm2 
      s12 = s12 + (conv12 + hpg1*hpg2*area*xk)*(1.d0+kr) + kr*sm2 
      s13 = s13 + (conv13 + hpg1*hpg3*area*xk)*(1.d0+kr) + kr*sm2 
      s23 = s23 + (conv23 + hpg2*hpg3*area*xk)*(1.d0+kr) + kr*sm2 
c
c     Matriz de massa * porosidade:
c      
      sm1 = sm1*e(2)
      sm2 = sm2*e(2)
c           
c     Residuo: R = ku + Mv
c
      p(1) = s11*u(1)+s12*u(2)+s13*u(3)+sm1*v(1)+sm2*(v(2)+v(3))
      p(2) = s21*u(1)+s22*u(2)+s23*u(3)+sm1*v(2)+sm2*(v(1)+v(3))
      p(3) = s31*u(1)+s32*u(2)+s33*u(3)+sm1*v(3)+sm2*(v(1)+v(2))
c
c     Matriz efetiva:  M + alfa.dt.K
c
      a1 = alfa*dt
      s(1,1) = s11*a1 + sm1
      s(2,2) = s22*a1 + sm1
      s(3,3) = s33*a1 + sm1
      s(2,1) = s21*a1 + sm2
      s(3,1) = s31*a1 + sm2
      s(3,2) = s32*a1 + sm2
      s(1,2) = s12*a1 + sm2
      s(1,3) = s13*a1 + sm2
      s(2,3) = s23*a1 + sm2
      return
c.......................................................................
  300 continue
c
c....
c      xj11 = xl(1,1)-xl(1,3)
c      xj12 = xl(2,1)-xl(2,3)
c      xj21 = xl(1,2)-xl(1,3)
c      xj22 = xl(2,2)-xl(2,3)
c      det  = xj11*xj22-xj12*xj21
c      area = det*0.5d0 
c      soma = soma + ((u(1)+u(2)+u(3))/3.d0)*area
      return
c ......................................................................  
c
c ... Cargas distribuidas no volume e no contorno:
c            
c ......................................................................
  500 continue
      if(iq(1) .eq. 0) return
      call tload(iq(1),t,u,a1)
      vx = (vn(1,1)+vn(1,2)+vn(1,3))/3.d0
      vy = (vn(2,1)+vn(2,2)+vn(2,3))/3.d0      
      xj11 = xl(1,1)-xl(1,3)
      xj12 = xl(2,1)-xl(2,3)
      xj21 = xl(1,2)-xl(1,3)
      xj22 = xl(2,2)-xl(2,3)
      det  = xj11*xj22-xj12*xj21
      hx(1) =  xj22/det
      hx(2) = -xj12/det
      hx(3) = (xj12-xj22)/det
      hy(1) = -xj21/det
      hy(2) =  xj11/det
      hy(3) = (xj21-xj11)/det
      hpg1 = (vx*hx(1)+vy*hy(1))*dt
      hpg2 = (vx*hx(2)+vy*hy(2))*dt
      hpg3 = (vx*hx(3)+vy*hy(3))*dt            
      area = det*0.5
      p(1) = -a1*(1.d0/3.d0 + hpg1)*area
      p(2) = -a1*(1.d0/3.d0 + hpg2)*area
      p(3) = -a1*(1.d0/3.d0 + hpg3)*area
      return            
c ......................................................................           
 1000 continue
      print*, '*** Subrotina elmt03: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end
      subroutine elmt04_t(e,iq,xl,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c **********************************************************************
c *                                                                    *
c *   elmt04: Elemento quadrilatero bilinear (difusao)                 *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nem)- coordenadas nodais locais                         *
c *     u(nst)     - valores prescritos                                *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     p   = nodal forces                                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi' 
      include 'gauss.fi'                 
      integer iq(*),ndm,nst,nel,isw,nint,lx,ly,i,j,i1,i2,i3,ma,nlit
      real*8 e(*),xl(ndm,*),u(*),v(*),p(*),s(nst,*),m(4,4),r(3,3)
      real*8 h(4),hx(4),hy(4),xj(2,2),xji(2,2),det,wt,w,hi(3)
      real*8 rn(4),sn(4),ri,si,q,tm,difus,capac,ea,atil,clat,dhi
      data   rn / 1.d0,-1.d0,-1.d0, 1.d0/,sn / 1.d0, 1.d0,-1.d0,-1.d0/
      data   nint / 2 /
c ......................................................................
      go to (100,200,200,500) isw
c ......................................................................
c
c ... Propriedades do material:
c
c ......................................................................
  100 continue
c ... e(1) = difusao
c     e(2) = coeficiente da derivada no tempo
      return
c ......................................................................
c
c.... Matrix K:
c
c ......................................................................
  200 continue
      difus = e(1)
      capac = e(2)
c  
c     Propriedades Termo - quimico : 
      ea    = e(3)
      atil  = e(4)
      clat  = e(5)
      call sfquad4(h,hx,hy,0.0d0,0.0d0,.true.,.false.)
      tm = h(1)*u(1)+h(2)*u(2)+h(3)*u(3)+h(4)*u(4)
      if (nlit .eq. 1)    then
         hi(1) = hi(2)
         hi(3) = tm
      endif
      call calc_hidr(hi(1),hi(2),dhi,hi(3),tm,atil,ea,ma,nlit)
      if (isw .eq. 3) goto 300
c
      if(ndm .eq. 3) then
         call rotmatrix(xl,r)
         call rotate(xl(1,1),r,xl(1,1),.false.)
         call rotate(xl(1,2),r,xl(1,2),.false.)
         call rotate(xl(1,3),r,xl(1,3),.false.)
         call rotate(xl(1,4),r,xl(1,4),.false.)
      endif
      do 206 i = 1, nst
         do 205 j = 1, nst
            s(i,j) = 0.d0
            m(i,j) = 0.d0
  205    continue
         p(i) = 0.d0
  206 continue
      do 240 ly = 1, nint
         si = pg(ly,nint)
         do 230 lx = 1, nint
            ri = pg(lx,nint)
            call sfquad4(h,hx,hy,ri,si,.true.,.true.)
c ......
            tm = h(1)*u(1)+h(2)*u(2)+h(3)*u(3)+h(4)*u(4)
            call updprop(difus,tm,ma,6)
            call updprop(capac,tm,ma,7)
c ......
            call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
            w  = wg(lx,nint)*wg(ly,nint)*det            
            wt = w*difus
            do 220 j = 1, 4
            do 210 i = 1, 4
c ------------------------------------------------------------------
               s(i,j) = s(i,j) + (hx(i)*hx(j)+ hy(i)*hy(j)) * wt
c ------------------------------------------------------------------
               m(i,j) = m(i,j) + h(i)*h(j) * w * capac
c ------------------------------------------------------------------
  210      continue
           p(j) = p(j) - capac*dhi*clat*w*h(j)
  220      continue
  230 continue
  240 continue
c ......................................................................  
c
c ... Forcas internas:
c  
c      call lku(s,u,p,nst)
      do 245 i = 1,nst
      do 245 j = 1, nst
         p(i) = p(i) + s(i,j)*u(j)
  245 continue
      call addlku(m,v,p,nst)
c ......................................................................        
c
c ... Matriz de rigidez efetiva:                     
c      
      do 260 i = 1, nst
      do 250 j = 1, nst
         s(i,j) = m(i,j)+alfa*dt*s(i,j)
  250 continue
  260 continue
      return
c ......................................................................
c
c.... Fluxos nos vertices:
c
c ... sistema local de coordenadas, se ndm = 3:
c
c ......................................................................
  300 continue
      if(ndm .eq. 3) then
         call rotmatrix(xl,r)
         call rotate(xl(1,1),r,xl(1,1),.false.)
         call rotate(xl(1,2),r,xl(1,2),.false.)
         call rotate(xl(1,3),r,xl(1,3),.false.)
         call rotate(xl(1,4),r,xl(1,4),.false.)
      endif
      do 310 i = 1, 4
         i1 = (i-1)*ndm+1
         i2 =  i1+1
         call sfquad4(h,hx,hy,rn(i),sn(i),.true.,.true.)
c ......
         tm = h(1)*u(1)+h(2)*u(2)+h(3)*u(3)+h(4)*u(4)
         call updprop(difus,tm,ma,6)
c ......
         call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
         p(i1) = -(hx(1)*u(1)+hx(2)*u(2)+hx(3)*u(3)+hx(4)*u(4))*difus
         p(i2) = -(hy(1)*u(1)+hy(2)*u(2)+hy(3)*u(3)+hy(4)*u(4))*difus
         if(ndm .eq. 3) then
            i3 = i2 + 1
            p(i3) = 0.d0
            call rotate(p(i1),r,p(i1),.true.)
         endif
  310 continue
      return
c ......................................................................  
c
c ... Cargas distribuidas no volume e no contorno:
c            
c ......................................................................
  500 continue
      if(iq(1) .eq. 0) return
      call tload(iq(1),t,u,q)
c     1 ponto de integracao:                                   
      call sfquad4(h,hx,hy,0.d0,0.d0,.true.,.true.)
      call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
      p(1) = -q*det
      p(2) = -q*det
      p(3) = -q*det
      p(4) = -q*det             
      return 
c ......................................................................                  
      stop
      end                  
      subroutine elmt04c_t(e,iq,xl,u,v,p,s,a,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT04: Elemento quadrilatero bilinear (conveccao-difusao-SUPG)  *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nem)- coordenadas nodais locais                         *
c *     u(nst)     - valores prescritos                                *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     p   = nodal forces                                             *
c *     a(ndm,*) - velocidade nos vertices                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'      
      include 'gauss.fi'                 
      integer iq(*),ndm,nst,nel,isw,nint,lx,ly,i,j,k,i1,i2
      real*8 e(*),xl(ndm,*),u(*),v(*),p(*),s(nst,*),a(ndm,*),m(4,4)
      real*8 h(4),hx(4),hy(4),hpg(4),xj(2,2),xji(2,2),det,wt,w
      real*8 rn(4),sn(4),ri,si,ax,ay,hc,q
      data   rn / 1.d0,-1.d0,-1.d0, 1.d0/,sn / 1.d0, 1.d0,-1.d0,-1.d0/
      data   nint / 2 /
c ......................................................................
      go to (100,200,300,500) isw
c ......................................................................
c
c ... Propriedades do material:
c
c ......................................................................
  100 continue
c ... e(1) = difusao
c     e(2) = coeficiente da derivada no tempo
      return
c ......................................................................
c
c.... Matrix K:
c
c ......................................................................
  200 continue
      do 205 i = 1, nst
      do 205 j = 1, nst
      s(i,j) = 0.d0
      m(i,j) = 0.d0
  205 continue
c
c ... tamanho caracteristico do elemento:
c
      call sfquad4(h,hx,hy,0.d0,0.d0,.true.,.true.)
      call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
      hc = dsqrt(det*4.d0)
c
c ... Velocidade no centro do elemento:
c
      ax = (a(1,1)+a(1,2)+a(1,3)+a(1,4))*0.25d0
      ay = (a(2,1)+a(2,2)+a(2,3)+a(2,4))*0.25d0
c ..................................................
      k = 0
      do 240 ly = 1, nint
         si = pg(ly,nint)
         do 230 lx = 1, nint
            ri = pg(lx,nint)
            k = k + 1
            i1 = (k-1)*ndm+1
            i2 = i1 + 1
            call sfquad4(h,hx,hy,ri,si,.true.,.true.)
            call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
            call hpgquad4(h,hx,hy,ax,ay,e(1),hc,hpg)
            w  = wg(lx,nint)*wg(ly,nint)*det            
            wt = w*e(1)
            do 220 j = 1, 4
            do 210 i = 1, 4
c ------------------------------------------------------------------
               s(i,j) = s(i,j) + (hx(i)*hx(j)+ hy(i)*hy(j))*wt + 
     .                            hpg(i)*(ax*hx(j)+ay*hy(j))*w
c ------------------------------------------------------------------
                m(i,j) = m(i,j) + hpg(i)*h(j) * w * e(2)
c ------------------------------------------------------------------
  210      continue
  220      continue
  230 continue
  240 continue
c ......................................................................  
c
c ... Forcas internas:
c
      call lku(s,u,p,nst)      
      call addlku(m,v,p,nst)
c ......................................................................        
c
c ... Matriz de rigidez efetiva:                     
c      
      do 260 i = 1, nst
      do 250 j = 1, nst
         s(i,j) = m(i,j)+alfa*dt*s(i,j)
  250 continue
  260 continue
      return
  300 continue
c ......................................................................
c
c.... Fluxos nos vertices:
c
c ......................................................................
c
c ... sistema local de coordenadas, se ndm = 3:
c
      do 310 i = 1, 4
         i1 = (i-1)*ndm+1
         i2 =  i1+1
         call sfquad4(h,hx,hy,rn(i),sn(i),.false.,.true.)
         call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
         p(i1) = -(hx(1)*u(1)+hx(2)*u(2)+hx(3)*u(3)+hx(4)*u(4))*e(1)
         p(i2) = -(hy(1)*u(1)+hy(2)*u(2)+hy(3)*u(3)+hy(4)*u(4))*e(1)
  310 continue
c      soma = soma + (u(1)+u(2)+u(3)+u(4))*25.d0  
      return
  400 continue
c ......................................................................
c
c.... Fluxos nos pontos de integracao:
c
c ......................................................................
c
c ... sistema local de coordenadas, se ndm = 3:
c
      i = 0
      do 420 ly = 1, nint      
         si = pg(ly,nint)
         do 410 lx = 1, nint
            i = i + 1
            ri = pg(lx,nint)
            i1 = (i-1)*ndm+1
            i2 =  i1+1
            call sfquad4(h,hx,hy,ri,si,.false.,.true.)
            call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
            p(i1) = -(hx(1)*u(1)+hx(2)*u(2)+hx(3)*u(3)+hx(4)*u(4))*e(1)
            p(i2) = -(hy(1)*u(1)+hy(2)*u(2)+hy(3)*u(3)+hy(4)*u(4))*e(1)
  410    continue
  420 continue
      return
c ......................................................................  
c
c ... Cargas distribuidas no volume e no contorno:
c            
c ......................................................................
  500 continue       
      if(iq(1) .eq. 0) return
      call tload(iq(1),t,u,q)
      ax = (a(1,1)+a(1,2)+a(1,3)+a(1,4))*0.25d0
      ay = (a(2,1)+a(2,2)+a(2,3)+a(2,4))*0.25d0 
c     1 ponto de integracao:                                   
      call sfquad4(h,hx,hy,0.d0,0.d0,.true.,.true.)
      call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
      w  = 4.d0*det                        
      hc = dsqrt(w)                  
      call hpgquad4(h,hx,hy,ax,ay,e(1),hc,hpg)            
      p(1) = -q*hpg(1)*w
      p(2) = -q*hpg(2)*w
      p(3) = -q*hpg(3)*w
      p(4) = -q*hpg(4)*w            
      return            
c ......................................................................           
      end
      subroutine elmt04b_t(e,iq,xl,u,v,p,s,a,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   elmt04: Elemento quadrilatero bilinear (conveccao-difusao,       *
c *           formulacao Characteristic-Galerkin)                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nem)- coordenadas nodais locais                         *
c *     u(nst)     - valores prescritos                                *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     p   = nodal forces                                             *
c *     a(ndm,*) - velocidade nos vertices                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'      
      include 'gauss.fi'                 
      integer iq(*),ndm,nst,nel,isw,nint,lx,ly,i,j,k,i1,i2
      real*8 e(*),xl(ndm,*),u(*),v(*),p(*),s(nst,*),a(ndm,*),m(4,4)
      real*8 h(4),hx(4),hy(4),hpg(4),xj(2,2),xji(2,2),det,wt,w,q
      real*8 rn(4),sn(4),ri,si,ax,ay
      data   rn / 1.d0,-1.d0,-1.d0, 1.d0/,sn / 1.d0, 1.d0,-1.d0,-1.d0/
      data   nint / 2 /
c ......................................................................
      go to (100,200,300,500) isw
c ......................................................................
c
c ... Propriedades do material:
c
c ......................................................................
  100 continue
c ... e(1) = difusao
c     e(2) = coeficiente da derivada no tempo
      return
c ......................................................................
c
c.... Matrix K:
c
c ......................................................................
  200 continue
      do 206 i = 1, nst
         p(i) = 0.d0
         do 205 j = 1, nst
            s(i,j)  = 0.d0
            m(i,j)  = 0.d0
  205    continue
  206 continue
c
c ... Velocidade no centro do elemento:
c
      ax = (a(1,1)+a(1,2)+a(1,3)+a(1,4))*0.25d0
      ay = (a(2,1)+a(2,2)+a(2,3)+a(2,4))*0.25d0
c ..................................................
      k = 0
      do 240 ly = 1, nint
         si = pg(ly,nint)
         do 230 lx = 1, nint
            ri = pg(lx,nint)
            k = k + 1
            i1 = (k-1)*ndm+1
            i2 = i1 + 1
            call sfquad4(h,hx,hy,ri,si,.true.,.true.)
            call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
            call uhx(hx,hy,ax,ay,hpg,4)            
            w  = wg(lx,nint)*wg(ly,nint)*det            
            wt = w*e(1)
            do 220 j = 1, 4
               do 210 i = 1, 4
c ----------------------------------------------------------------------
                  s(i,j) = s(i,j) + (hx(i)*hx(j)+ hy(i)*hy(j))*wt +
     .                    (h(i)*hpg(j) + 0.5d0*dt*hpg(i)*hpg(j))*w     
c ----------------------------------------------------------------------
                  m(i,j) = m(i,j) + h(i)*h(j) * w * e(2)
c ----------------------------------------------------------------------
  210         continue
  220      continue
  230 continue
  240 continue
c ......................................................................  
c
c ... Forcas internas:
c 
      call lku(s,u,p,nst)
      call addlku(m,v,p,nst)
c ......................................................................        
c
c ... Matriz de rigidez efetiva:                     
c
      do 260 i = 1, nst
      do 250 j = 1, nst
         s(i,j) = m(i,j)+alfa*dt*s(i,j)         
  250 continue
  260 continue
      return
c ......................................................................
c
c.... Fluxos nos vertices:
c
c ......................................................................
  300 continue
c
c ... sistema local de coordenadas, se ndm = 3:
c
      do 310 i = 1, 4
         i1 = (i-1)*ndm+1
         i2 =  i1+1
         call sfquad4(h,hx,hy,rn(i),sn(i),.false.,.true.)
         call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
         p(i1) = -(hx(1)*u(1)+hx(2)*u(2)+hx(3)*u(3)+hx(4)*u(4))*e(1)
         p(i2) = -(hy(1)*u(1)+hy(2)*u(2)+hy(3)*u(3)+hy(4)*u(4))*e(1)
  310 continue
c      soma = soma + (u(1)+u(2)+u(3)+u(4))*25.d0
      return
  400 continue
c ......................................................................
c
c.... Fluxos nos pontos de integracao:
c
c ......................................................................
c
c ... sistema local de coordenadas, se ndm = 3:
c
      i = 0
      do 420 ly = 1, nint      
         si = pg(ly,nint)
         do 410 lx = 1, nint
            i = i + 1
            ri = pg(lx,nint)
            i1 = (i-1)*ndm+1
            i2 =  i1+1
            call sfquad4(h,hx,hy,ri,si,.false.,.true.)
            call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
            p(i1) = -(hx(1)*u(1)+hx(2)*u(2)+hx(3)*u(3)+hx(4)*u(4))*e(1)
            p(i2) = -(hy(1)*u(1)+hy(2)*u(2)+hy(3)*u(3)+hy(4)*u(4))*e(1)
  410    continue
  420 continue
      return 
c ......................................................................  
c
c ... Cargas distribuidas no volume e no contorno:
c            
c ......................................................................
  500 continue       
      if(iq(1) .eq. 0) return
      call tload(iq(1),t,u,q)
      ax = (a(1,1)+a(1,2)+a(1,3)+a(1,4))*0.25d0
      ay = (a(2,1)+a(2,2)+a(2,3)+a(2,4))*0.25d0
c     1 ponto de integracao:            
      call sfquad4(h,hx,hy,0.d0,0.d0,.true.,.true.)
      call jacob2d(hx,hy,xj,xji,xl,det,4,ndm,nel)
      call uhx(hx,hy,ax,ay,hpg,4)                        
      w  = 4.d0*det                        
      p(1) = -q*(h(1)+hpg(1))*w
      p(2) = -q*(h(2)+hpg(2))*w
      p(3) = -q*(h(3)+hpg(3))*w
      p(4) = -q*(h(4)+hpg(4))*w            
      return
c ......................................................................  
      end             
      subroutine elmt06_t(e,iq,xl,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c **********************************************************************
c *                                                                    *
c *   elmt06: Elemento tetraedrico (difusao)                           *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nem)- coordenadas nodais locais                         *
c *     u(nst)     - valores prescritos                                *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     p   = nodal forces                                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'                        
      integer iq(*),ndm,nst,nel,isw,i,j,ma,nlit
      real*8 e(*),xl(ndm,*),u(nst),v(*),p(*),s(nst,nst)
      real*8 hx(4),hy(4),hz(4),det,wt,tm,difus,capac,cp,ro
      real*8 xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33,sm1,sm2
      real*8 r(3,3),x1(3),x2(3),x3(3)
      real*8 ea,clat,hi(3),atil,dhi
c ......................................................................
      go to (100,200,200,400) isw
c ......................................................................
c
c ... Propriedades do material:
c
c ......................................................................
  100 continue
c     e(1) = coeficiente de difusao
c     e(2) = calor especifico
c     e(3) = massa especifica
c
c     e(2) = coeficiente da derivada no tempo
      return
c ......................................................................
c
c.... Matrix K:
c
c ......................................................................
  200 continue
c
      difus = e(1)
      capac = e(2)
c      cp = e(2)
c      ro = e(3)
      tm = (u(1)+u(2)+u(3)+u(4))/4.d0
      call updprop(difus,tm,ma,6)
      call updprop(capac,tm,ma,7)
c      call updprop(ro,tm,ma,8)
c      capac = cp*ro
c     Propriedades Termo - quimico : 
      ea    = e(3)
      atil  = e(4)
      clat  = e(5)
      if (ea .eq. 0.d0) goto 205
      if (nlit .eq. 1)    then
         hi(1) = hi(2)
         hi(3) = tm
      endif
      call calc_hidr(hi(1),hi(2),dhi,hi(3),tm,atil,ea,ma,nlit)
c 
  205 continue
c      
      if (isw .eq. 3) goto 300
c
c ... Matriz Jacobiana:
c
      xj11 = xl(1,1)-xl(1,4)
      xj12 = xl(2,1)-xl(2,4)
      xj13 = xl(3,1)-xl(3,4)
      xj21 = xl(1,2)-xl(1,4)
      xj22 = xl(2,2)-xl(2,4)
      xj23 = xl(3,2)-xl(3,4)
      xj31 = xl(1,3)-xl(1,4)
      xj32 = xl(2,3)-xl(2,4)
      xj33 = xl(3,3)-xl(3,4)
      det  = xj11*xj22*xj33 + xj12*xj23*xj31 + xj13*xj21*xj32
     .      -xj31*xj22*xj13 - xj12*xj21*xj33 - xj11*xj32*xj23 
      if (det .le. 0.d0) go to 1000
c
c ... Derivadas das funcoes de interpolacao:
c
      hx(1) =  (xj22*xj33-xj23*xj32)/det
      hy(1) =  (xj23*xj31-xj21*xj33)/det
      hz(1) =  (xj21*xj32-xj22*xj31)/det
      hx(2) =  (xj13*xj32-xj12*xj33)/det
      hy(2) =  (xj11*xj33-xj13*xj31)/det
      hz(2) =  (xj12*xj31-xj11*xj32)/det
      hx(3) =  (xj12*xj23-xj13*xj22)/det
      hy(3) =  (xj13*xj21-xj11*xj23)/det
      hz(3) =  (xj11*xj22-xj12*xj21)/det
      hx(4) = -hx(1)-hx(2)-hx(3)
      hy(4) = -hy(1)-hy(2)-hy(3)
      hz(4) = -hz(1)-hz(2)-hz(3)
c ......................................................................        
c
c ... Matriz K:                     
c
      wt = det*difus/6.d0
      do 220 i = 1, 4
      do 210 j = 1, 4
         s(i,j) = (hx(i)*hx(j) + hy(i)*hy(j) + hz(i)*hz(j))*wt
  210 continue
  220 continue
c
c ... Matriz de massa:
c
      sm1 = (1.d0/60.d0)*det*capac
      sm2 = sm1*0.5d0  
c ......................................................................  
c
c ... Forcas internas:
c      
      p(1)= s(1,1)*u(1) + s(1,2)*u(2) + s(1,3)*u(3) + s(1,4)*u(4) +
     .         sm1*v(1) + sm2*(v(2)+v(3)+v(4)) 
     .      - capac*dhi*clat*det/24.d0
      p(2)= s(1,2)*u(1) + s(2,2)*u(2) + s(2,3)*u(3) + s(2,4)*u(4) +
     .         sm1*v(2) + sm2*(v(1)+v(3)+v(4))
     .      - capac*dhi*clat*det/24.d0
      p(3)= s(1,3)*u(1) + s(2,3)*u(2) + s(3,3)*u(3) + s(3,4)*u(4) +
     .         sm1*v(3) + sm2*(v(1)+v(2)+v(4))
     .      - capac*dhi*clat*det/24.d0
      p(4)= s(1,4)*u(1) + s(2,4)*u(2) + s(3,4)*u(3) + s(4,4)*u(4) +
     .         sm1*v(4) + sm2*(v(1)+v(2)+v(3))
     .      - capac*dhi*clat*det/24.d0
c
c ... Matriz de rigidez efetiva:                     
c
      s(1,1) = sm1 + s(1,1)*alfa*dt
      s(1,2) = sm2 + s(1,2)*alfa*dt
      s(1,3) = sm2 + s(1,3)*alfa*dt
      s(1,4) = sm2 + s(1,4)*alfa*dt      
      s(2,2) = sm1 + s(2,2)*alfa*dt
      s(2,3) = sm2 + s(2,3)*alfa*dt
      s(2,4) = sm2 + s(2,4)*alfa*dt      
      s(3,3) = sm1 + s(3,3)*alfa*dt
      s(3,4) = sm2 + s(3,4)*alfa*dt
      s(4,4) = sm1 + s(4,4)*alfa*dt            
      s(2,1) = s(1,2)
      s(3,1) = s(1,3)
      s(3,2) = s(2,3)
      s(4,1) = s(1,4)
      s(4,2) = s(2,4)
      s(4,3) = s(3,4)   
      return
c ......................................................................
c
c.... Fluxos nos vertices:
c
c ......................................................................
  300 continue
c
c ... Matriz Jacobiana:
c
      xj11 = xl(1,1)-xl(1,4)
      xj12 = xl(2,1)-xl(2,4)
      xj13 = xl(3,1)-xl(3,4)
      xj21 = xl(1,2)-xl(1,4)
      xj22 = xl(2,2)-xl(2,4)
      xj23 = xl(3,2)-xl(3,4)
      xj31 = xl(1,3)-xl(1,4)
      xj32 = xl(2,3)-xl(2,4)
      xj33 = xl(3,3)-xl(3,4)
      det  = xj11*xj22*xj33 + xj12*xj23*xj31 + xj13*xj21*xj32
     .      -xj31*xj22*xj13 - xj12*xj21*xj33 - xj11*xj32*xj23 
      if (det .le. 0.d0) go to 1000
c
c ... Derivadas das funcoes de interpolacao:
c
      hx(1) =  (xj22*xj33-xj23*xj32)/det
      hy(1) =  (xj23*xj31-xj21*xj33)/det
      hz(1) =  (xj21*xj32-xj22*xj31)/det
      hx(2) =  (xj13*xj32-xj12*xj33)/det
      hy(2) =  (xj11*xj33-xj13*xj31)/det
      hz(2) =  (xj12*xj31-xj11*xj32)/det
      hx(3) =  (xj12*xj23-xj13*xj22)/det
      hy(3) =  (xj13*xj21-xj11*xj23)/det
      hz(3) =  (xj11*xj22-xj12*xj21)/det
      hx(4) = -hx(1)-hx(2)-hx(3)
      hy(4) = -hy(1)-hy(2)-hy(3)
      hz(4) = -hz(1)-hz(2)-hz(3)
c
c ... Fluxo:
c
      p(1)  = -(hx(1)*u(1)+hx(2)*u(2)+hx(3)*u(3)+hx(4)*u(4))*difus
      p(2)  = -(hy(1)*u(1)+hy(2)*u(2)+hy(3)*u(3)+hy(4)*u(4))*difus
      p(3)  = -(hz(1)*u(1)+hz(2)*u(2)+hz(3)*u(3)+hz(4)*u(4))*difus
      p(4)  = p(1)
      p(5)  = p(2)
      p(6)  = p(3)
      p(7)  = p(1)
      p(8)  = p(2)
      p(9)  = p(3)
      p(10) = p(1)
      p(11) = p(2)
      p(12) = p(3)            
      return
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
      if (iq(1) .ne. 0) then
         call rotmatrix2(xl(1,1),xl(1,3),xl(1,2),r)
         call rotate(xl(1,1),r,x1(1),.false.)
         call rotate(xl(1,3),r,x2(1),.false.)
         call rotate(xl(1,2),r,x3(1),.false.) 
         xj11 = x1(1)-x3(1)
         xj12 = x1(2)-x3(2)
         xj21 = x2(1)-x3(1)
         xj22 = x2(2)-x3(2)
         det  = xj11*xj22-xj12*xj21
         if (det .le. 0.d0) go to 1000
         tm = (u(1)+u(2)+u(3))/3.d0
         call tload(iq(1),t,tm,wt)
         p(1)=-wt*det/6.d0
         p(2)=-wt*det/6.d0
         p(3)=-wt*det/6.d0
      endif
      if (iq(2) .ne. 0) then
         call rotmatrix2(xl(1,1),xl(1,2),xl(1,4),r)
         call rotate(xl(1,1),r,x1(1),.false.)
         call rotate(xl(1,2),r,x2(1),.false.)
         call rotate(xl(1,4),r,x3(1),.false.) 
         xj11 = x1(1)-x3(1)
         xj12 = x1(2)-x3(2)
         xj21 = x2(1)-x3(1)
         xj22 = x2(2)-x3(2)
         det  = xj11*xj22-xj12*xj21
         if (det .le. 0.d0) go to 1000
         tm = (u(1)+u(2)+u(4))/3.d0
         call tload(iq(2),t,tm,wt)
         p(1)=-wt*det/6.d0
         p(2)=-wt*det/6.d0
         p(4)=-wt*det/6.d0
      endif
      if (iq(3) .ne. 0) then
         call rotmatrix2(xl(1,2),xl(1,3),xl(1,4),r)
         call rotate(xl(1,2),r,x1(1),.false.)
         call rotate(xl(1,3),r,x2(1),.false.)
         call rotate(xl(1,4),r,x3(1),.false.) 
         xj11 = x1(1)-x3(1)
         xj12 = x1(2)-x3(2)
         xj21 = x2(1)-x3(1)
         xj22 = x2(2)-x3(2)
         det  = xj11*xj22-xj12*xj21
         if (det .le. 0.d0) go to 1000
         tm = (u(2)+u(3)+u(4))/3.d0
         call tload(iq(3),t,tm,wt)
         p(2)=-wt*det/6.d0
         p(3)=-wt*det/6.d0
         p(4)=-wt*det/6.d0
      endif
      if (iq(4) .ne. 0) then     
         call rotmatrix2(xl(1,3),xl(1,1),xl(1,4),r)
         call rotate(xl(1,3),r,x1(1),.false.)
         call rotate(xl(1,1),r,x2(1),.false.)
         call rotate(xl(1,4),r,x3(1),.false.) 
         xj11 = x1(1)-x3(1)
         xj12 = x1(2)-x3(2)
         xj21 = x2(1)-x3(1)
         xj22 = x2(2)-x3(2)
         det  = xj11*xj22-xj12*xj21 
         if (det .le. 0.d0) go to 1000
         tm = (u(1)+u(3)+u(4))/3.d0
         call tload(iq(4),t,tm,wt)
         p(3)=-wt*det/6.d0
         p(1)=-wt*det/6.d0
         p(4)=-wt*det/6.d0
      endif
      return       
c ...................................................................... 
     
 1000 continue
      print*, '*** Subrotina ELMT06: determinante nulo ou negativo do el
     .emento ',nel
      stop
      end                                         
      subroutine elmt07_t(e,iq,xl,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c **********************************************************************
c *                                                                    *
c *   elmt07: Elemento hexaedrico (difusao)                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nem)- coordenadas nodais locais                         *
c *     u(nst)     - valores prescritos                                *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     p   = nodal forces                                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'                  
      include 'gauss.fi'                 
      integer iq(*),ndm,nst,nel,isw,nint,lx,ly,lz,i,j,i1,i2,i3,ma
      integer nlit
      real*8 e(*),xl(ndm,*),u(*),v(*),p(*),s(nst,*),m(8,8),q
      real*8 h(8),hx(8),hy(8),hz(8),xj(3,3),xji(3,3),det,w,wt
      real*8 rn(8),sn(8),tn(8),ri,si,ti,ea,atil,clat,hi(3),dhi
      real*8 tm,difus,capac
      data   rn / 1.d0,-1.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0, 1.d0 /,
     .       sn / 1.d0, 1.d0,-1.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0 /,
     .       tn / 1.d0, 1.d0, 1.d0, 1.d0,-1.d0,-1.d0,-1.d0,-1.d0 /
      data   nint / 2 /
c ......................................................................
      go to (100,200,200,400) isw
c ......................................................................
c
c ... Propriedades do material:
c
c ......................................................................
  100 continue
c ... e(1) = coeficiente de difusao
c     e(2) = coeficiente da derivada no tempo
      return
c ......................................................................
  200 continue
c
c ... definicao das propriedades
c
      difus = e(1)
      capac = e(2)
      call sfhexa8(h,hx,hy,hz,0.0,0.0,0.0,.true.,.false.)
      tm = h(1)*u(1)+h(2)*u(2)+h(3)*u(3)+h(4)*u(4)+
     .     h(5)*u(5)+h(6)*u(6)+h(7)*u(7)+h(8)*u(8)
      call updprop(difus,tm,ma,6)
      call updprop(capac,tm,ma,7)
c     Propriedades Termo - quimico : 
      ea    = e(3)
      atil  = e(4)
      clat  = e(5)
      if (nlit .eq. 1)    then
         hi(1) = hi(2)
         hi(3) = tm
      endif
      call calc_hidr(hi(1),hi(2),dhi,hi(3),tm,atil,ea,ma,nlit)
      if (isw .eq. 3) goto 300
c
c ... Matriz K:                     
c
      do 205 i = 1, nst
      do 206 j = 1, nst
      s(i,j) = 0.d0
      m(i,j) = 0.d0
  206 continue
      p(i) = 0.d0
  205 continue
      do 250 lz = 1, nint
         ti = pg(lz,nint)
         do 240 ly = 1, nint
            si = pg(ly,nint)
            do 230 lx = 1, nint
               ri = pg(lx,nint)
               call sfhexa8(h,hx,hy,hz,ri,si,ti,.true.,.true.)
               call jacob3d(hx,hy,hz,xj,xji,xl,det,8,nel)
               tm = h(1)*u(1)+h(2)*u(2)+h(3)*u(3)+h(4)*u(4)+
     .              h(5)*u(5)+h(6)*u(6)+h(7)*u(7)+h(8)*u(8)
               call updprop(difus,tm,ma,6)
               call updprop(capac,tm,ma,7)
               w  = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
               wt = w*difus
               do 220 j = 1, 8
               do 210 i = 1, 8
c ---------------------------------------------------------------------
               s(i,j) = s(i,j)+(hx(i)*hx(j)+hy(i)*hy(j)+hz(i)*hz(j))*wt
c ---------------------------------------------------------------------
               m(i,j) = m(i,j) + h(i)*h(j)*w*capac
c ---------------------------------------------------------------------               
  210         continue
              p(j) = p(j) - capac*dhi*clat*w*h(j)
  220         continue
  230       continue
  240    continue
  250 continue
c ......................................................................
c
c ... Forcas internas:
c  
      call addlku(s,u,p,nst)
      call addlku(m,v,p,nst)
c ......................................................................      
c
c ... Matriz de rigidez efetiva:                     
c      
      do 270 i = 1, nst
      do 260 j = 1, nst
         s(i,j) = m(i,j)+alfa*dt*s(i,j)
  260 continue
  270 continue
c ......................................................................  
      return
c ......................................................................
c
c.... Fluxos nos vertices:
c
c ......................................................................
  300 continue
      do 310 i = 1, 8
         i1 = (i-1)*3+1
         i2 = i1+1
         i3 = i2+1
         call sfhexa8(h,hx,hy,hz,rn(i),sn(i),tn(i),.true.,.true.)
         call jacob3d(hx,hy,hz,xj,xji,xl,det,8,nel)    
         tm = h(1)*u(1)+h(2)*u(2)+h(3)*u(3)+h(4)*u(4)+
     .        h(5)*u(5)+h(6)*u(6)+h(7)*u(7)+h(8)*u(8)
         call updprop(difus,tm,ma,6)   
         p(i1) = -(hx(1)*u(1)+hx(2)*u(2)+hx(3)*u(3)+hx(4)*u(4)+
     .             hx(5)*u(5)+hx(6)*u(6)+hx(7)*u(7)+hx(8)*u(8))*difus
         p(i2) = -(hy(1)*u(1)+hy(2)*u(2)+hy(3)*u(3)+hy(4)*u(4)+
     .             hy(5)*u(5)+hy(6)*u(6)+hy(7)*u(7)+hy(8)*u(8))*difus
         p(i3) = -(hz(1)*u(1)+hz(2)*u(2)+hz(3)*u(3)+hz(4)*u(4)+
     .             hz(5)*u(5)+hz(6)*u(6)+hz(7)*u(7)+hz(8)*u(8))*difus
  310 continue
      return
c ......................................................................
  400 continue
      if(iq(1) .eq. 0) return
c      call tload(iq(1),t,u,q)
c     1 ponto de integracao:
c      call sfhexa8(h,hx,hy,hz,0.d0,0.d0,0.d0,.true.,.true.)
c      call jacob3d(hx,hy,hz,xj,xji,xl,det,8,nel)
c      w  = 8.d0*det                        
c      p(1) = -q*w
c      p(2) = -q*w
c      p(3) = -q*w
c      p(4) = -q*w
c      p(5) = -q*w
c      p(6) = -q*w
c      p(7) = -q*w
c      p(8) = -q*w
c
c ....carga aplicada na face do hexaedro- definicoes de face:
c
c     iq(1) = 1 | nós 1 - 2 - 3 - 4 |  
c             2 | nós 2 - 6 - 7 - 3 |  
c             3 | nós 5 - 8 - 7 - 6 |  
c             4 | nós 8 - 5 - 1 - 4 |  
c             5 | nós 5 - 6 - 2 - 1 |  
c             6 | nós 7 - 8 - 4 - 3 |  
c            
c ......................................................................
c      if (iq(1) .eq. 1) then
c         xj11 = xl(1,1)-xl(1,3)
c         xj12 = xl(2,1)-xl(2,3)
c         xj21 = xl(1,2)-xl(1,3)
c         xj22 = xl(2,2)-xl(2,3)
c         det  = xj11*xj22-xj12*xj21
c         tm = (u(1)+u(2)+u(3))/3.d0
c         call tload(iq(2),t,tm,wt)
c         p(1)=-wt*2/det
c         p(2)=-wt*2/det
c         p(3)=-wt*2/det
c      elseif (iq(1) .eq. 2) then
c         xj11 = xl(1,1)-xl(1,4)
c         xj12 = xl(2,1)-xl(2,4)
c         xj21 = xl(1,2)-xl(1,4)
c         xj22 = xl(2,2)-xl(2,4)
c         det  = xj11*xj22-xj12*xj21
c         tm = (u(1)+u(2)+u(4))/3.d0
c         call tload(iq(2),t,tm,wt)
c         p(1)=-wt*2/det
c         p(2)=-wt*2/det
c         p(4)=-wt*2/det
c      elseif (iq(1) .eq. 3) then
c         xj11 = xl(1,2)-xl(1,4)
c         xj12 = xl(2,2)-xl(2,4)
c         xj21 = xl(1,3)-xl(1,4)
c         xj22 = xl(2,3)-xl(2,4)
c         det  = xj11*xj22-xj12*xj21
c         tm = (u(2)+u(3)+u(4))/3.d0
c         call tload(iq(2),t,tm,wt)
c         p(2)=-wt*2/det
c         p(3)=-wt*2/det
c        p(4)=-wt*2/det
c      elseif (iq(1) .eq. 4) then     
c         xj11 = xl(1,3)-xl(1,4)
c         xj12 = xl(2,3)-xl(2,4)
c         xj21 = xl(1,1)-xl(1,4)
c        xj22 = xl(2,1)-xl(2,4)
c        det  = xj11*xj22-xj12*xj21 
c        tm = (u(1)+u(3)+u(4))/3.d0
c         call tload(iq(2),t,tm,wt)
c         p(3)=-wt*2/det
c        p(1)=-wt*2/det
c         p(4)=-wt*2/det
c      endif
      return  
c ......................................................................                 
      stop
      end
      subroutine elmt07c_t(e,iq,xl,u,v,p,s,a,ndm,nst,nel,isw)
c **********************************************************************
c *                                                                    *
c *   elmt07: Elemento hexaedrico (conveccao-difusao, SUPG)            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     xl(ndm,nem)- coordenadas nodais locais                         *
c *     u(nst)     - valores prescritos                                *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     p   = nodal forces                                             *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - correcao do vetor de forcas                                *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'                  
      include 'gauss.fi'                 
      integer iq(*),ndm,nst,nel,isw,nint,lx,ly,lz,i,j,i1,i2,i3
      real*8 e(*),xl(ndm,*),u(*),v(*),p(*),s(nst,*),a(ndm,*),m(8,8)
      real*8 h(8),hx(8),hy(8),hz(8),hpg(8),xj(3,3),xji(3,3),det,w,wt
      real*8 rn(8),sn(8),tn(8),ri,si,ti,ax,ay,az,hc,q
      data   rn / 1.d0,-1.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0, 1.d0 /,
     .       sn / 1.d0, 1.d0,-1.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0 /,
     .       tn / 1.d0, 1.d0, 1.d0, 1.d0,-1.d0,-1.d0,-1.d0,-1.d0 /
      data   nint / 2 /
c ......................................................................
      go to (100,200,300,500) isw
c ......................................................................
c
c ... Propriedades do material:
c
c ......................................................................
  100 continue
c ... e(1) = coeficiente de difusao
c     e(2) = coeficiente da derivada no tempo
      return
c ......................................................................
  200 continue
c
c ... Matriz K:                     
c
      do 205 i = 1, nst
      do 205 j = 1, nst
      s(i,j) = 0.d0
      m(i,j) = 0.d0
  205 continue
c
c ... tamanho caracteristico do elemento:
c
      call sfhexa8(h,hx,hy,hz,0.d0,0.d0,0.d0,.true.,.true.)
      call jacob3d(hx,hy,hz,xj,xji,xl,det,8,nel)
      hc = (det*8.d0)**(1.d0/3.d0)
c
c ... Velocidade no centro do elemento:
c
      ax =(a(1,1)+a(1,2)+a(1,3)+a(1,4)+a(1,5)+a(1,6)+a(1,7)+a(1,8))/8.d0
      ay =(a(2,1)+a(2,2)+a(2,3)+a(2,4)+a(2,5)+a(2,6)+a(2,7)+a(2,8))/8.d0
      az =(a(3,1)+a(3,2)+a(3,3)+a(3,4)+a(3,5)+a(3,6)+a(3,7)+a(3,8))/8.d0
c ......................................................................  
      do 250 lz = 1, nint
         ti = pg(lz,nint)
         do 240 ly = 1, nint
            si = pg(ly,nint)
            do 230 lx = 1, nint
               ri = pg(lx,nint)
               call sfhexa8(h,hx,hy,hz,ri,si,ti,.true.,.true.)
               call jacob3d(hx,hy,hz,xj,xji,xl,det,8,nel)
               call hpghexa8(h,hx,hy,hz,ax,ay,az,e(1),hc,hpg)
               w  = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
               wt = w*e(1)
               do 220 j = 1, 8
               do 210 i = 1, 8
c ---------------------------------------------------------------------
               s(i,j) = s(i,j)+(hx(i)*hx(j)+hy(i)*hy(j)+hz(i)*hz(j))*wt
     .                + hpg(i)*(ax*hx(j)+ay*hy(j)+az*hz(j))*w
c ---------------------------------------------------------------------
               m(i,j) = m(i,j) + hpg(i)*h(j)*w*e(2)
c ---------------------------------------------------------------------               
  210         continue
  220         continue
  230       continue
  240    continue
  250 continue
c ......................................................................
c
c ... Forcas internas:
c  
      call lku(s,u,p,nst)
      call addlku(m,v,p,nst)
c ......................................................................      
c
c ... Matriz de rigidez efetiva:                     
c      
      do 270 i = 1, nst
      do 260 j = 1, nst
         s(i,j) = m(i,j)+alfa*dt*s(i,j)
  260 continue
  270 continue
c ......................................................................  
      return
c ......................................................................
c
c.... Fluxos nos vertices:
c
c ......................................................................
  300 continue
      do 310 i = 1, 8
         i1 = (i-1)*3+1
         i2 = i1+1
         i3 = i2+1
         call sfhexa8(h,hx,hy,hz,rn(i),sn(i),tn(i),.false.,.true.)
         call jacob3d(hx,hy,hz,xj,xji,xl,det,8,nel)         
         p(i1) = -(hx(1)*u(1)+hx(2)*u(2)+hx(3)*u(3)+hx(4)*u(4)+
     .             hx(5)*u(5)+hx(6)*u(6)+hx(7)*u(7)+hx(8)*u(8))*e(1)
         p(i2) = -(hy(1)*u(1)+hy(2)*u(2)+hy(3)*u(3)+hy(4)*u(4)+
     .             hy(5)*u(5)+hy(6)*u(6)+hy(7)*u(7)+hy(8)*u(8))*e(1)
         p(i3) = -(hz(1)*u(1)+hz(2)*u(2)+hz(3)*u(3)+hz(4)*u(4)+
     .             hz(5)*u(5)+hz(6)*u(6)+hz(7)*u(7)+hz(8)*u(8))*e(1)
  310 continue
      return
c ......................................................................
  500 continue
      if(iq(1) .eq. 0) return
      call tload(iq(1),t,u,q)
      ax =(a(1,1)+a(1,2)+a(1,3)+a(1,4)+a(1,5)+a(1,6)+a(1,7)+a(1,8))/8.d0
      ay =(a(2,1)+a(2,2)+a(2,3)+a(2,4)+a(2,5)+a(2,6)+a(2,7)+a(2,8))/8.d0
      az =(a(3,1)+a(3,2)+a(3,3)+a(3,4)+a(3,5)+a(3,6)+a(3,7)+a(3,8))/8.d0
c ......................................................................        
c     1 ponto de integracao:
      call sfhexa8(h,hx,hy,hz,0.d0,0.d0,0.d0,.true.,.true.)
      call jacob3d(hx,hy,hz,xj,xji,xl,det,8,nel)
      w  = 8.d0*det                        
      hc = w**(1.d0/3.d0)
      call hpghexa8(h,hx,hy,hz,ax,ay,az,e(1),hc,hpg)
      p(1) = -q*hpg(1)*w
      p(2) = -q*hpg(2)*w
      p(3) = -q*hpg(3)*w
      p(4) = -q*hpg(4)*w
      p(5) = -q*hpg(5)*w
      p(6) = -q*hpg(6)*w
      p(7) = -q*hpg(7)*w
      p(8) = -q*hpg(8)*w            
      return
c ......................................................................                 
      stop
      end    
      subroutine elmt08_t3d(e,iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT02_I: Elemento de Interface prisma termico                   *
c *   ------                                                           *
c *  1--3                                                              *
c *  |\/ \  h -> espessura = 0.0                                       *
c *  |/4--6                                                            *
c *  2 | /                                                             *
c *   \|/                                                              *
c *    5                                                               *
c *	  !nst = 6                                                         *
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
c *     p - forcas internas (isw = 2,3) ou fluxo (isw = 4)             *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'            
      integer ndm,nst,nel,isw,ma,nlit
      integer i,j,k,l,q,i1,i2,i3,j1,j2,j3
      real*8 e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),r(3,3)
      real*8 det,wt,espes(3)
      real*8 xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33,tm
      real*8 xji11,xji12,xji13,xji21,xji22,xji23,xji31,xji32,xji33
      real*8 ym,nu,a,b,c,a1,a2,a3,d11,d22,d33
      real*8 lDxDksi,lDyDksi,lDxDeta,lDyDeta,lDzDksi,lDzDeta
      real*8 uDxDksi,uDyDksi,uDxDeta,uDyDeta,uDzDksi,uDzDeta
      real*8 udet,ldet,uJ1,uJ2,uJ3,lJ1,lJ2,lJ3,peso,auxdet,paux(9)
      real*8 lxj11,lxj12,lxj21,lxj22,uxj11,uxj12,uxj21,uxj22,thic
      real*8 difus,capac,sm1,sm2,iq,hi(*)
      logical flagt    
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
c     d11 = e(1) !d11
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
      if (thic .eq. 0.d0) thic = 0.5d0*udet/1e+09
      difus = e(6)/thic
      capac = e(7)

      if (isw .eq. 3) goto 300  
c
c     Matriz de rigidez:                     
c
      wt = udet*difus/12.d0
      do i = 1, 6
         do j = 1, 6
            if (i .eq. j)   then
               s(i,j) = wt
            elseif(i .gt. 3 .and. j .le. 3 .or. 
     .         i .le. 3 .and. j .gt. 3 ) then
               s(i,j) = - wt*0.5d0
            else
               s(i,j) = wt*0.5d0
            endif
         enddo
      enddo
      s(4,1) = -wt
      s(5,2) = -wt
      s(6,3) = -wt
      s(1,4) = -wt
      s(2,5) = -wt
      s(3,6) = -wt
c
c ..... Matriz de Massa:
c
c      sm1 = (1.d0/12.d0)*udet*capac
c      sm2 = sm1*0.5d0
c
c      Forcas Internas:
c
c      p(1)=s(1,1)*u(1)+s(1,2)*u(2)+s(1,3)*u(3)+sm1*v(1)+sm2*(v(2)+v(3))
c      p(2)=s(2,1)*u(1)+s(2,2)*u(2)+s(2,3)*u(3)+sm1*v(2)+sm2*(v(1)+v(3))
c      p(3)=s(3,1)*u(1)+s(3,2)*u(2)+s(3,3)*u(3)+sm1*v(3)+sm2*(v(1)+v(2))
c      p(4)=s(4,1)*u(1)+s(4,2)*u(2)+s(4,3)*u(3)+sm1*v(4)+sm2*(v(5)+v(6))
c      p(5)=s(5,1)*u(1)+s(5,2)*u(2)+s(5,3)*u(3)+sm1*v(5)+sm2*(v(4)+v(6))
c      p(6)=s(6,1)*u(1)+s(6,2)*u(2)+s(6,3)*u(3)+sm1*v(6)+sm2*(v(4)+v(5))
      p(1)=s(1,1)*u(1)+s(1,2)*u(2)+s(1,3)*u(3)+s(1,4)*u(4)+s(1,5)*u(5)+
     .s(1,6)*u(6)
      p(2)=s(2,1)*u(1)+s(2,2)*u(2)+s(2,3)*u(3)+s(2,4)*u(4)+s(2,5)*u(5)+
     .s(2,6)*u(6)
      p(3)=s(3,1)*u(1)+s(3,2)*u(2)+s(3,3)*u(3)+s(3,4)*u(4)+s(3,5)*u(5)+
     .s(3,6)*u(6)
      p(4)=s(4,1)*u(1)+s(4,2)*u(2)+s(4,3)*u(3)+s(4,4)*u(4)+s(4,5)*u(5)+
     .s(4,6)*u(6)
      p(5)=s(5,1)*u(1)+s(5,2)*u(2)+s(5,3)*u(3)+s(5,4)*u(4)+s(5,5)*u(5)+
     .s(5,6)*u(6)
      p(6)=s(6,1)*u(1)+s(6,2)*u(2)+s(6,3)*u(3)+s(6,4)*u(4)+s(6,5)*u(5)+
     .s(6,6)*u(6)
c
c ..... Matriz de rigidez efetiva:
c
      do i = 1,6
         do j = 1,6
            s(i,j) = alfa*dt*s(i,j)
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
      subroutine elmt08_t3dhexa(e,iq,x,u,v,p,s,hi,ndm,nst,nel,ma,
     .                                                  nlit,isw)
c **********************************************************************
c *                                                                    *
c *   ELMT02_I: Elemento de Interface hexaedro termico                 *
c *   ------                                                           *
c *  1--4                                                              *
c *  |\ |\  h -> espessura = 0.0                                       *
c *  | 2--3                                                            *
c *  5-|-8|                                                            *
c *   \|  \                                                            *
c *    6--7                                                            *
c *	  !nst = 8                                                       *
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
c *     p - forcas internas (isw = 2,3) ou fluxo (isw = 4)             *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'     
      common /gauss/ pg,wg       
      integer ndm,nst,nel,isw,ma,nlit,nint
      integer i,j,k,l,q,i1,i2,i3,j1,j2,j3,lx,ly
      real*8 e(*),x(ndm,*),u(*),v(*),p(*),s(nst,nst),r(3,3)
      real*8 det,wt,espes(3),sl(8,8)
      real*8 xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33,tm
      real*8 xji11,xji12,xji13,xji21,xji22,xji23,xji31,xji32,xji33
      real*8 ym,nu,a,b,c,a1,a2,a3,d11,d22,d33
      real*8 h(8),hx(4),hy(4),xj(2,2),xji(2,2)
      real*8  pg(10,10),wg(10,10),ri,si
      real*8 udet,ldet,uJ1,uJ2,uJ3,lJ1,lJ2,lJ3,peso,auxdet,paux(9)
      real*8 lxj11,lxj12,lxj21,lxj22,uxj11,uxj12,uxj21,uxj22,thic
      real*8 difus,capac,sm1,sm2,iq,hi(*)
      logical flagt    
      data    nint/2/
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
c     d11 = e(1) !d11
c     d22 = e(2) !d22
c     d33 = e(3) !d33
c
c ...................................................................... 
c ... Propriedades do material: 
c ...................................................................... 
c      thic = e(4)                     
c      if (thic .eq. 0.d0) thic = 0.5d0*udet/1e+09
      difus = e(6)
      capac = e(7)
      thic = 0.5d0/1e+09
      difus = difus/thic
c ... Matriz de rigidez
c 
      call rotmatrix(x,r)
      call rotate(x,r,x,.false.) 
      call rotate(x(1,2),r,x(1,2),.false.)
      call rotate(x(1,3),r,x(1,3),.false.)
      call rotate(x(1,4),r,x(1,4),.false.)
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
            call jacob2d(hx,hy,xj,xji,x,det,4,ndm,nel)
            wt = wg(lx,nint)*wg(ly,nint)*det*difus
            h(5) = -h(1)
            h(6) = -h(2)
            h(7) = -h(3)
            h(8) = -h(4)
            do 220 j = 1, 8
               do 210 i = 1, 8
                  s(i,j)     = s(i,j) + h(i)*h(j)*wt
  210          continue
  220       continue
  230    continue
  240 continue
      if (isw .eq. 3) goto 300 
c      call rotmatrix(x,r)
c      call irotatek3d(s,r,nst)
c
c ..... Forcas Internas
c
c
c      Forcas Internas:
c
c      p(1)=s(1,1)*u(1)+s(1,2)*u(2)+s(1,3)*u(3)+sm1*v(1)+sm2*(v(2)+v(3))
c      p(2)=s(2,1)*u(1)+s(2,2)*u(2)+s(2,3)*u(3)+sm1*v(2)+sm2*(v(1)+v(3))
c      p(3)=s(3,1)*u(1)+s(3,2)*u(2)+s(3,3)*u(3)+sm1*v(3)+sm2*(v(1)+v(2))
c      p(4)=s(4,1)*u(1)+s(4,2)*u(2)+s(4,3)*u(3)+sm1*v(4)+sm2*(v(5)+v(6))
c      p(5)=s(5,1)*u(1)+s(5,2)*u(2)+s(5,3)*u(3)+sm1*v(5)+sm2*(v(4)+v(6))
c      p(6)=s(6,1)*u(1)+s(6,2)*u(2)+s(6,3)*u(3)+sm1*v(6)+sm2*(v(4)+v(5))
      p(1)=s(1,1)*u(1)+s(1,2)*u(2)+s(1,3)*u(3)+s(1,4)*u(4)+s(1,5)*u(5)+
     .s(1,6)*u(6)+s(1,7)*u(7)+s(1,8)*u(8)
      p(2)=s(2,1)*u(1)+s(2,2)*u(2)+s(2,3)*u(3)+s(2,4)*u(4)+s(2,5)*u(5)+
     .s(2,6)*u(6)+s(2,7)*u(7)+s(2,8)*u(8)
      p(3)=s(3,1)*u(1)+s(3,2)*u(2)+s(3,3)*u(3)+s(3,4)*u(4)+s(3,5)*u(5)+
     .s(3,6)*u(6)+s(3,7)*u(7)+s(3,8)*u(8)
      p(4)=s(4,1)*u(1)+s(4,2)*u(2)+s(4,3)*u(3)+s(4,4)*u(4)+s(4,5)*u(5)+
     .s(4,6)*u(6)+s(4,7)*u(7)+s(4,8)*u(8)
      p(5)=s(5,1)*u(1)+s(5,2)*u(2)+s(5,3)*u(3)+s(5,4)*u(4)+s(5,5)*u(5)+
     .s(5,6)*u(6)+s(5,7)*u(7)+s(5,8)*u(8)
      p(6)=s(6,1)*u(1)+s(6,2)*u(2)+s(6,3)*u(3)+s(6,4)*u(4)+s(6,5)*u(5)+
     .s(6,6)*u(6)+s(6,7)*u(7)+s(6,8)*u(8)
      p(7)=s(7,1)*u(1)+s(7,2)*u(2)+s(7,3)*u(3)+s(7,4)*u(4)+s(7,5)*u(5)+
     .s(7,6)*u(6)+s(7,7)*u(7)+s(7,8)*u(8)
      p(8)=s(8,1)*u(1)+s(8,2)*u(2)+s(8,3)*u(3)+s(8,4)*u(4)+s(8,5)*u(5)+
     .s(8,6)*u(6)+s(8,7)*u(7)+s(8,8)*u(8)
c      do i = 1, 8
c         p(i) = 0.d0
c      enddo
c ... Matriz de rigidez efetiva
      do i = 1, nst
      do j = 1, nst
         s(i,j) = alfa*dt*s(i,j)
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

