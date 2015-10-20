c*****************************Svn***************************************      
c*$Date: 2013-04-19 11:09:37 -0300 (Fri, 19 Apr 2013) $                 
c*$Rev: 967 $                                                           
c*$Author: ana $                                                   
c***********************************************************************
      subroutine elmlibpmec(e,iq,x,u,p,s,dt,ndm,nst,nel,iel,isw,ma,nlit
     .                     ,ilib)
c **********************************************************************
c *                                                                    *
c *   ELMLIB: biblioteca de elementos.                                 *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     iq(7) - cargas nos elementos                                   *
c *     x(ndm,nen)- coordenadas nodais locais                          *
c *     u(nst)     - solucao anterior                                  *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     iel - tipo   do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *     ma   -  numero de material do elemento                         *
c *     nlit -  numero da iteracao nao linear                          *
c *     ilib -  codigo da biblioteca                                   *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - residuo                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer iq(*),iel,nel,ndm,nst,isw,ilib,ma,nlit
      real*8 e(*),x(*),u(*),p(*),s(nst,*)
      real*8 dt
c ......................................................................
      goto (100,200,300,400,500,600,700,800,900,1000,1100,1200) iel
   10 write(*,2000) iel,nel
      stop
c ......................................................................
  100 continue
c     if (ilib .eq. 1) then
c     endif 
      return
c ......................................................................            
  200 continue
c     if (ilib .eq. 1) then
c     endif 
      return
c ......................................................................
  300 continue
c     if (ilib .eq. 1) then
c     endif  
      return
c ......................................................................
  400 continue
c     if (ilib .eq. 1) then
c     endif 
      return
c ......................................................................
  500 continue
c     if (ilib .eq. 1) then
c     endif  
      return
c ......................................................................
  600 continue
c      if (ilib .eq. 1) then
c      endif    
      return
c ......................................................................
  700 continue
      if (ilib .eq. 1) then  
c     Elemento hexaedrico de 20 nos (poromec)
        call elmt10(e,x,u,p,s,dt,ndm,nst,nel,isw)
      endif
      return
c ......................................................................
  800 continue
c      if (ilib .eq. 1) then
c      endif    
      return
c ......................................................................            
  900 continue
c     if (ilib .eq. 1) then
c     endif 
      return 
c ......................................................................
 1000 continue
c     if (ilib .eq. 1) then
c     endif 
      return 
c ......................................................................
 1100 continue
c     if (ilib .eq. 1) then
c     endif 
      return       
c ......................................................................
 1200 continue
c     if (ilib .eq. 1) then
c     endif 
      return       
c ......................................................................        
 2000 format(1x,'SUBROUTINE ELMLIBPMEC:'
     .,/,5x,'tipo de elemento ',i2,' nao existente, elemento ',i6,' !')
       end
c **********************************************************************
c
c **********************************************************************       
      subroutine elmlib(e,iq,x,u,v,a,w,p,s,tl,ut,hi,u0,tx,txp,eps,ndm,
     .                nst,nel,iel,isw,ma,nlit,ilib,flaghidr,plastic)
c **********************************************************************
c *                                                                    *
c *   ELMLIB: biblioteca de elementos.                                 *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     iq(7) - cargas nos elementos                                   *
c *     x(ndm,nem)- coordenadas nodais locais                          *
c *     u(nst)     - solucao anterior                                  *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     iel - tipo   do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *     ilib -  codigo da biblioteca                                   *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - residuo                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer iq(*),iel,nel,ndm,nst,isw,ilib,ma,nlit
      real*8 plastic(*)
      real*8 e(*),x(*),u(*),v(*),a(*),w(*),p(*),s(nst,*),tl(*),ut(*)
      real*8 hi(*), tx(*),eps(*),u0(*),txp(*)
      logical flaghidr
c ......................................................................
      goto (100,200,300,400,500,600,700,800,900,1000,1100,1200) iel
   10 write(*,2000) iel,nel
      stop
c ......................................................................
  100 continue
      if (ilib .eq. 1) then
c     Elemento unidimensional (difusao):
      call elmt01_t(e(6),iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)     
c     Elemento unidimensional (conveccao-difusao, SUPG):  
c     call elmt01c_t(e(6),iq,x,u,v,p,s,w,ndm,nst,nel,isw)
      elseif (ilib .eq. 2) then
      goto 10
      endif
      return
c ......................................................................            
  200 continue
      if (ilib .eq. 1) then
c     Elemento triangular linear (difusao):
      call elmt03_t(e(6),iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c     Elemento triangular linear (conveccao-difusao,SUPG):
c     call elmt03c_t(e(6),iq,x,u,v,p,s,w,ndm,nst,nel,isw)
c     Elemento triangular linear (conveccao-difusao,CG):
c     call elmt03b_t(e(6),iq,x,u,v,p,s,w,ndm,nst,nel,isw)        
      elseif (ilib .eq. 2) then
c ... Estado plano de deformacao - triangulo linear:
      call elmt02_m(e,x,u,v,a,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,ma,
     .nlit,isw,flaghidr)
      endif
      return
c ......................................................................
  300 continue
      if (ilib .eq. 1) then
c     Elemento triangular linear (difusao):
      call elmt03_t(e(6),iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
      elseif (ilib .eq. 2) then
c ... Estado plano de tensao - triangulo linear:
      call elmt03_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,ma,nlit,
     . isw,flaghidr)
      endif
      return
c ......................................................................
  400 continue
      if (ilib .eq. 1) then  
c     Elemento quadrilatero bilinear (difusao)
      call elmt04_t(e(6),iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c     Elemento quadrialtero linear (conveccao-difusao,SUPG)
c     call elmt04c_t(e(6),iq,x,u,v,p,s,w,ndm,nst,nel,isw)
c     Elemento quadrialtero linear (conveccao-difusao,CG)          
c     call elmt04b_t(e(6),iq,x,u,v,p,s,w,ndm,nst,nel,isw)  
      elseif (ilib .eq. 2) then
c ... Estado plano de deformaçao - quadrilatero linear:
      call elmt04_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,ma,nlit,
     .                                                 isw,flaghidr)
      endif
      return
c ......................................................................
  500 continue
      if (ilib .eq. 1) then  
c     Elemento quadrilatero bilinear (difusao)
      call elmt04_t(e(6),iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
      elseif (ilib .eq. 2) then
c ... Estado plano de tensao - quadrilatero linear:
      call elmt05_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,ma,nlit,
     .                                                isw,flaghidr)
      endif
      return
c ......................................................................
  600 continue
      if (ilib .eq. 1) then  
c     Elemento tetraedrico (difusao)
      call elmt06_t(e(6),iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
      elseif (ilib .eq. 2) then
c ... Elasticidade 3D - tetraedro de 4 nos:
      call elmt06_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,ma,nlit,
     .      isw,flaghidr,plastic)
      endif      
      return
c ......................................................................
  700 continue
      if (ilib .eq. 1) then  
c     Elemento hexaedrico de 8 nos (difusao)
      call elmt07_t(e(6),iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c     Elemento hexaedrico de 8 nos (conveccao-difusao, SUPG)
c     call elmt07c_t(e(6),iq,x,u,v,p,s,w,ndm,nst,nel,isw)
      elseif (ilib .eq. 2) then
c ... Elasticidade 3D - hexaedro de 8 nos:
      call elmt07_m(e,x,u,v,p,s,hi,tx,txp,eps,tl,ut,ndm,nst,nel,ma,nlit,
     . isw,flaghidr,plastic)
c      call elmt07quad_m(e,x,u,v,p,s,tl,ut,ndm,nst,nel,ma,isw)
      endif      
      return
c ......................................................................
  800 continue
      if (ilib .eq. 1) then  
c ... nao implementado:  
c ... Elasticidade 3D - elemento cunha:
c      call elmt08_m(e,x,u,v,p,s,ndm,nst,nel,isw)
      elseif (ilib .eq. 2) then
c ... Elasticidade 3D - elemento interface:
      call elmt08_m(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,ma,nlit,isw)
      endif        
      return
c ......................................................................            
  900 continue
      if (ilib .eq. 1) then  
c ... nao implementado:  
      write (*,*) 'elemento nao implementado'
      stop
c ... Elasticidade 3D - elemento cunha:
c      call elmt08_m(e,x,u,v,p,s,ndm,nst,nel,isw)
      elseif (ilib .eq. 2) then
c ... Elasticidade 2D - elemento interface - grampo:
      call elmt08_grampo(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,ma,nlit,isw)
      endif 
      return 
c ......................................................................
 1000 continue
       if (ilib .eq. 1) then  
c ... Termico 3D - elemento interface Prisma:
      call elmt08_t3d(e,iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
c      stop
      elseif (ilib .eq. 2) then
c ... Elasticidade 3D - elemento interface:
      call elmt08_3d(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,ma,nlit,isw)
c      call elmt09_mItetra(e,x,u,v,p,s,tl,ndm,nst,nel,isw)
      endif 
      return 
c ......................................................................
 1100 continue
       if (ilib .eq. 1) then  
c ... Termico 3D - elemento interface Hexa:
      call elmt08_t3dhexa(e,iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
      elseif (ilib .eq. 2) then
c ... Elasticidade 3D - elemento interface Hexa:
      call elmt08_3dhexa(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,ma,nlit,isw,
     .   plastic)
      endif 
      return       
c ......................................................................
 1200 continue
       if (ilib .eq. 1) then  
c ... Termico 3D - elemento interface Hexa:
      call elmt08_t3dhexa(e,iq,x,u,v,p,s,hi,ndm,nst,nel,ma,nlit,isw)
      elseif (ilib .eq. 2) then
c ... Elasticidade 3D - elemento interface Hexa - grampo:
      call elmt08_3dhexa_grampo(e,x,u,v,p,s,tx,u0,tl,ndm,nst,nel,ma,
     .     nlit,isw,plastic)
      endif 
      return       
c ......................................................................        
 2000 format(1x,'SUBROUTINE ELMLIB:',/,5x,'tipo de elemento ',i2,' nao e
     .xistente, elemento ',i6,' !')
      end
c **********************************************************************
c
c **********************************************************************
      subroutine uhx(hx,hy,vx,vy,hxu,nen)
c **********************************************************************
c *                                                                    *
c *   Termos de estabilzacao CG (characteristic-Galerkin)              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   hx(nen) - derivadas de h                                         *
c *   hy(nen) - derivadas de h                                         *
c *   vx,vy - componentes da velocidade                                *
c *   kd    - coeficiente de difusao                                   *
c *   hc    - tamanho caracteristico do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   hpg(nen) - termos de estabilizacao: hpg = v.grad(h)              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,i
      real*8  hx(*),hy(*),vx,vy,hxu(*)
c ......................................................................
      do 100 i = 1, nen
         hxu(i) = vx*hx(i) + vy*hy(i)
  100 continue
c ......................................................................      
      return
      end           
      subroutine hpgquad4(h,hx,hy,vx,vy,kd,hc,hpg)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   HPGQUAD4: funcoes de interpolacao de Petrov-Galerkin             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   h(4)  - funcoes de interpolacao                                  *
c *   hx(4) - derivadas de h                                           *
c *   hy(4) - derivadas de h                                           *
c *   vx,vy - componeentes da velocodade                               *
c *   kd    - coeficiente de difusao                                   *
c *   hc    - tamanho caracteristico do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   hpg(4) - funcoes de Petrov-Galerkin                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      real*8  h(*),hx(*),hy(*),vx,vy,kd,hc,hpg(*)
      real*8 zero,peclet,alpha,xk,vm
      parameter (zero = 1.d-14)
c ......................................................................
c
c ... Modulo da velocidade:
c
      vm = dsqrt(vx*vx + vy*vy)
c
c ... Parametro de upwind:
c
      if(kd .gt. zero) then
           peclet = (hc*vm)/(2.d0*kd)
           alpha  = min((1.d0/3.d0)*peclet,1.d0)
      else
           alpha = 1.d0
      end if
c
c ... difusao artificial:
c
      xk = 0.d0
      if(vm .gt. zero) xk = alpha*hc*0.5d0/vm
c
c ... funcoes de Petrov-Galerkin:
c
      do 100 i = 1, 4
         hpg(i) = h(i) + xk*(vx*hx(i) + vy*hy(i))
  100 continue
c ......................................................................      
      return
      end
      subroutine hpghexa8(h,hx,hy,hz,vx,vy,vz,kd,hc,hpg)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   HPGHEXA8: funcoes de interpolacao de Petrov-Galerkin             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   h(8)  - funcoes de interpolacao                                  *
c *   hx(8) - derivadas de h                                           *
c *   hy(8) - derivadas de h                                           *
c *   hz(8) - derivadas de h                                           *
c *   vx,vy,vz - componeentes da velocodade                            *
c *   kd    - coeficiente de difusao                                   *
c *   hc    - tamanho caracteristico do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   hpg(4) - funcoes de Petrov-Galerkin                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      real*8  h(*),hx(*),hy(*),hz(*),vx,vy,vz,kd,hc,hpg(*)
      real*8 zero,peclet,alpha,xk,vm
      parameter (zero = 1.d-14)
c ......................................................................
c
c ... Modulo da velocidade:
c
      vm = dsqrt(vx*vx + vy*vy + vz*vz)
c
c ... Parametro de upwind:
c
      if(kd .gt. zero) then
           peclet = (hc*vm)/(2.d0*kd)
           alpha  = min((1.d0/3.d0)*peclet,1.d0)
      else
           alpha = 1.d0
      end if
c
c ... difusao artificial:
c
      xk = 0.d0
      if(vm .gt. zero) xk = alpha*hc*0.5d0/vm
c
c ... funcoes de Petrov-Galerkin:
c
      do 100 i = 1, 8
         hpg(i) = h(i) + xk*(vx*hx(i) + vy*hy(i) + vz*hz(i))
  100 continue
c ......................................................................      
      return
      end      
      subroutine lku(s,u,p,nst)
c **********************************************************************
c *                                                                    *
c *   LKU: forcas internas K.U                                         *
c *   ---                                                              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     s    - matriz de elemento                                      *
c *     u    - deslocamentos                                           *
c *     nst  - numero de graus de liberdade por elemento               *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     p = s.u - forcas internas                                      *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nst,i,j
      real*8  s(nst,*),u(*),p(*)
c ......................................................................
      do 200 i = 1, nst
         p(i) = 0.d0
         do 100 j = 1, nst
            p(i) = p(i) + s(i,j)*u(j)
  100    continue
  200 continue
c ......................................................................  
      return
      end
      subroutine addlku(s,u,p,nst)
c **********************************************************************
c *                                                                    *
c *   ADDLKU: forcas internas K.U                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     s    - matriz de elemento                                      *
c *     u    - deslocamentos                                           *
c *     nst  - numero de graus de liberdade por elemento               *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     p = s.u - forcas internas                                      *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nst,i,j
      real*8  s(nst,*),u(*),p(*)
c ......................................................................
      do 200 i = 1, nst
         do 100 j = 1, nst
            p(i) = p(i) + s(i,j)*u(j)
  100    continue
  200 continue
c ......................................................................  
      return
      end      
      subroutine sftria3(h,hr,hs,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    22/12/01        *
c *   SFTRIA3:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do triangulo de 3 nos.                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s - coordenadas naturais                                     *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(6)    - funcoes de interolocao no ponto (r,s)                *
c *     hr(6)   - derivadas de h em relacao a r                        *
c *     hs(6)   - derivadas de h em relacao a s                        *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  h(*),hr(*),hs(*),r,s,t
      logical afl,bfl
c ......................................................................
      t = 1.d0 - r - s
      if (afl) then
c
c ...... Funcoes de interpolacao quadraticas standard:
c
         h(1) = r
         h(2) = s
         h(3) = t
      endif
      if (bfl) then
c
c ...... Derivadas em relacao a r :
c
         hr(1) =  1.d0
         hr(2) =  0.d0
         hr(3) = -1.d0
c
c ...... Derivadas em relacao a s :
c
         hs(1) =  0.
         hs(2) =  1.d0
         hs(3) = -1.d0
      endif
c ......................................................................                    
      return
      end
      subroutine sftria6(h,hr,hs,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    22/12/01        *
c *   SFTRIA6:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do triangulo de 6 nos.                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s - coordenadas naturais                                     *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(6)    - funcoes de interolocao no ponto (r,s)                *
c *     hr(6)   - derivadas de h em relacao a r                        *
c *     hs(6)   - derivadas de h em relacao a s                        *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  h(*),hr(*),hs(*),r,s,t
      logical afl,bfl
c ......................................................................
      t = 1.d0 - r - s
      if (afl) then
c
c ...... Funcoes de interpolacao quadraticas standard:
c
         h(1) = r * (2.d0*r-1)
         h(2) = s * (2.d0*s-1)
         h(3) = t * (2.d0*t-1)
         h(4) = 4.d0 * r * s
         h(5) = 4.d0 * s * t
         h(6) = 4.d0 * t * r
      endif
      if (bfl) then
c
c ...... Derivadas em relacao a r :
c
         hr(1) =  4.d0* r - 1.d0
         hr(2) =  0.d0
         hr(3) =  1.d0- 4.d0* t
         hr(4) =  4.d0* s
         hr(5) = -4.d0* s
         hr(6) =  4.d0*(t-r)
c
c ...... Derivadas em relacao a s :
c
         hs(1) =  0.d0
         hs(2) =  4.d0* s - 1.d0
         hs(3) =  1.d0 - 4.d0* t
         hs(4) =  4.d0* r
         hs(5) =  4.d0*(t-s)
         hs(6) = -4.d0* r
      endif
c ......................................................................      
      return
      end
      subroutine sftetra10(h,hr,hs,ht,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    22/12/01        *
c *   SFTETRA10:                                                       *
c *   ---------                                                        *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t) do tetraedro de 10 nos.                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   r,s,t - coordenadas naturais                                     *
c *   afl = true ( calcula funcoes )                                   *
c *   bfl = true ( calcula as derivadas )                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   h(10)    - funcoes de interolocao no ponto (r,s)                 *
c *   hr(10)   - derivadas de h em relacao a r                         *
c *   hs(10)   - derivadas de h em relacao a s                         *
c *   ht(10)   - derivadas de h em relacao a t                         *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  h(*),hr(*),hs(*),ht(*),r,s,t,u
      logical afl,bfl
c ......................................................................
      u = 1.d0 - r - s - t
      if (afl) then
c
c ...... Funcoes de interpolacao quadraticas standard:
c
         h(1) = r * (2.d0*r-1)
         h(2) = s * (2.d0*s-1)
         h(3) = t * (2.d0*t-1)
         h(4) = u * (2.d0*u-1)         
         h(5) = 4.d0 * r * s
         h(6) = 4.d0 * r * t
         h(7) = 4.d0 * r * u
         h(8) = 4.d0 * s * t
         h(9) = 4.d0 * t * u
         h(10)= 4.d0 * s * u         
      endif
      if (bfl) then
c
c ...... Derivadas em relacao a r :
c
         hr(1) =  4.d0*r - 1.d0
         hr(2) =  0.d0
         hr(3) =  0.d0
         hr(4) =  1.d0 - 4.d0*u
         hr(5) =  4.d0*s
         hr(6) =  4.d0*t
         hr(7) =  4.d0*(u-r)
         hr(8) =  0.d0
         hr(9) = -4.d0*t 
         hr(10)= -4.d0*s         
c
c ...... Derivadas em relacao a s :
c
         hs(1) =  0.d0
         hs(2) =  4.d0*s - 1.d0
         hs(3) =  0.d0
         hs(4) =  1.d0 - 4.d0*u
         hs(5) =  4.d0*r
         hs(6) =  0.d0
         hs(7) = -4.d0*r
         hs(8) =  4.d0*t
         hs(9) = -4.d0*t 
         hs(10)=  4.d0*(u-s)
c
c ...... Derivadas em relacao a t :
c
         ht(1) =  0.d0
         ht(2) =  0.d0
         ht(3) =  4.d0*t - 1.d0
         ht(4) =  1.d0 - 4.d0*u
         ht(5) =  0.d0
         ht(6) =  4.d0*r
         ht(7) = -4.d0*r
         ht(8) =  4.d0*s
         ht(9) =  4.d0*(u-t)
         ht(10)= -4.d0*s
      endif
c ......................................................................      
      return
      end      
      subroutine sfquad4(h,hx,hy,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    27/01/89        *
c *   SFQUAD4:                                                         *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do quadrilatero de 4 nos.                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   r,s - coordenadas naturais                                       *
c *   afl = true ( calcula funcoes )                                   *
c *   bfl = true ( calcula as derivadas )                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   h(4)    - funcoes de interolocao no ponto (r,s)                  *
c *   hx(4)   - derivadas de h em relacao a r                          *
c *   hy(4)   - derivadas de h em relacao a s                          *
c *                                                                    *
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8  h(*),hx(*),hy(*),r,s
c ......................................................................
      if ( afl ) then
c
c     funcoes de interpolacao
c
      h(1)  = ( 1.d0+r ) * ( 1.d0+s ) / 4.d0
      h(2)  = ( 1.d0-r ) * ( 1.d0+s ) / 4.d0
      h(3)  = ( 1.d0-r ) * ( 1.d0-s ) / 4.d0
      h(4)  = ( 1.d0+r ) * ( 1.d0-s ) / 4.d0
c
      endif
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
      hx(1)  =   (1.d0+s) / 4.d0
      hx(2)  = - (1.d0+s) / 4.d0
      hx(3)  = - (1.d0-s) / 4.d0
      hx(4)  =   (1.d0-s) / 4.d0
c
c     derivadas em relacao a s :
c
      hy(1)  =  (1.d0+r) / 4.d0
      hy(2)  =  (1.d0-r) / 4.d0
      hy(3)  = -(1.d0-r) / 4.d0
      hy(4)  = -(1.d0+r) / 4.d0
c
      endif
c ......................................................................      
      return
      end
      subroutine sfquad8(h,hx,hy,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    27/01/89        *
c *   SFQUAD8:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do quadrilatero quadratico de 8 nos.              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s - coordenadas naturais                                     *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(8)    - funcoes de interolocao no ponto (r,s)                *
c *     hx(8)   - derivadas de h em relacao a r                        *
c *     hy(8)   - derivadas de h em relacao a s                        *
c *                                                                    *
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8  h(*),hx(*),hy(*),r,s
c ......................................................................
      if ( afl ) then
c
c     funcoes de interpolacao
c
      h(5)  = ( 1.d0-r*r ) * ( 1.d0+s ) / 2.
      h(6)  = ( 1.d0-s*s ) * ( 1.d0-r ) / 2.
      h(7)  = ( 1.d0-r*r ) * ( 1.d0-s ) / 2.
      h(8)  = ( 1.d0-s*s ) * ( 1.d0+r ) / 2.
      h(1)  = ( 1.d0+r ) * ( 1.d0+s ) / 4.d0 - h(5) / 2.d0 - h(8) / 2.d0
      h(2)  = ( 1.d0-r ) * ( 1.d0+s ) / 4.d0 - h(5) / 2.d0 - h(6) / 2.d0
      h(3)  = ( 1.d0-r ) * ( 1.d0-s ) / 4.d0 - h(6) / 2.d0 - h(7) / 2.d0
      h(4)  = ( 1.d0+r ) * ( 1.d0-s ) / 4.d0 - h(7) / 2.d0 - h(8) / 2.d0
c
      endif
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
      hx(5)  = - r * (1.d0+s)
      hx(6)  = -.5d0 * (1.d0-s*s)
      hx(7)  = - r * (1.d0-s)
      hx(8)  =  .5d0 * (1.d0-s*s)
      hx(1)  =   (1.d0+s) / 4.d0 - hx(5) / 2.d0 - hx(8) / 2.d0
      hx(2)  = - (1.d0+s) / 4.d0 - hx(5) / 2.d0 - hx(6) / 2.d0
      hx(3)  = - (1.d0-s) / 4.d0 - hx(6) / 2.d0 - hx(7) / 2.d0
      hx(4)  =   (1.d0-s) / 4.d0 - hx(7) / 2.d0 - hx(8) / 2.d0
c
c     derivadas em relacao a s :
c
      hy(5)  =  (1.d0-r*r) / 2.d0
      hy(6)  = - s * (1.d0-r)
      hy(7)  = -(1.d0-r*r) / 2.d0
      hy(8)  = - s * (1.+r)
      hy(1)  =  (1.d0+r) / 4.d0 - hy(5) / 2.d0 - hy(8) / 2.d0
      hy(2)  =  (1.d0-r) / 4.d0 - hy(5) / 2.d0 - hy(6) / 2.d0
      hy(3)  = -(1.d0-r) / 4.d0 - hy(6) / 2.d0 - hy(7) / 2.d0
      hy(4)  = -(1.d0+r) / 4.d0 - hy(7) / 2.d0 - hy(8) / 2.d0
c
      endif
c ......................................................................      
      return
      end
      subroutine sfhexa8(h,hx,hy,hz,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    27/01/89        *
c *   SFHEXA8:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t), do hexaedro linear de 8 nos.                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s,t - coordenadas naturais                                   *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(8)    - funcoes de interolocao no ponto (r,s,t)              *
c *     hx(8)   - derivadas de h em relacao a r                        *
c *     hy(8)   - derivadas de h em relacao a s                        *
c *     hz(8)   - derivadas de h em relacao a t                        *
c *                                                                    *
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8  h(*),hx(*),hy(*),hz(*),r,s,t
c ......................................................................
      if ( afl ) then
c
c     funcoes de interpolacao
c
      h(1)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
      h(2)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
      h(3)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
      h(4)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
      h(5)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
      h(6)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
      h(7)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
      h(8)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
c
      endif
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
      hx(1)  =   (1.d0+s) * (1.d0+t) / 8.d0
      hx(2)  = - (1.d0+s) * (1.d0+t) / 8.d0
      hx(3)  = - (1.d0-s) * (1.d0+t) / 8.d0
      hx(4)  =   (1.d0-s) * (1.d0+t) / 8.d0
      hx(5)  =   (1.d0+s) * (1.d0-t) / 8.d0
      hx(6)  = - (1.d0+s) * (1.d0-t) / 8.d0
      hx(7)  = - (1.d0-s) * (1.d0-t) / 8.d0
      hx(8)  =   (1.d0-s) * (1.d0-t) / 8.d0
c
c     derivadas em relacao a s :
c
      hy(1)  =  (1.d0+r) * (1.d0+t) / 8.d0
      hy(2)  =  (1.d0-r) * (1.d0+t) / 8.d0
      hy(3)  = -(1.d0-r) * (1.d0+t) / 8.d0
      hy(4)  = -(1.d0+r) * (1.d0+t) / 8.d0
      hy(5)  =  (1.d0+r) * (1.d0-t) / 8.d0
      hy(6)  =  (1.d0-r) * (1.d0-t) / 8.d0
      hy(7)  = -(1.d0-r) * (1.d0-t) / 8.d0
      hy(8)  = -(1.d0+r) * (1.d0-t) / 8.d0
c
c     derivadas em relacao a t :
c
      hz(1)  =  ( 1.d0+r ) * ( 1.d0+s ) / 8.d0
      hz(2)  =  ( 1.d0-r ) * ( 1.d0+s ) / 8.d0
      hz(3)  =  ( 1.d0-r ) * ( 1.d0-s ) / 8.d0
      hz(4)  =  ( 1.d0+r ) * ( 1.d0-s ) / 8.d0
      hz(5)  = -( 1.d0+r ) * ( 1.d0+s ) / 8.d0
      hz(6)  = -( 1.d0-r ) * ( 1.d0+s ) / 8.d0
      hz(7)  = -( 1.d0-r ) * ( 1.d0-s ) / 8.d0
      hz(8)  = -( 1.d0+r ) * ( 1.d0-s ) / 8.d0
      endif
c ......................................................................      
      return
      end
      subroutine sfhexa20(h,hx,hy,hz,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    30/09/05        *
c *   SFHEXA20:                                                        *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t), do hexaedro quadratico de 20 nos.              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s,t - coordenadas naturais                                   *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(2)    - funcoes de interolocao no ponto (r,s,t)              *
c *     hx(20)  - derivadas de h em relacao a r                        *
c *     hy(20)  - derivadas de h em relacao a s                        *
c *     hz(20)  - derivadas de h em relacao a t                        *
c *                                                                    *
c * OBS:                                                               *
c * no  ( r, s, t)                                                     *
c * no1 ( 1, 1, 1)                                                     * 
c * no2 (-1, 1, 1)                                                     * 
c * no3 (-1,-1, 1)                                                     * 
c * no4 ( 1,-1, 1)                                                     * 
c * no5 ( 1, 1,-1)                                                     * 
c * no6 (-1, 1,-1)                                                     * 
c * no7 (-1,-1,-1)                                                     * 
c * no8 ( 1,-1,-1)                                                     * 
c * no9 ( 0, 1, 1)                                                     * 
c * no10(-1, 0, 1)                                                     * 
c * no11( 0,-1, 1)                                                     * 
c * no12( 1, 0, 1)                                                     * 
c * no13( 0, 1,-1)                                                     * 
c * no14(-1, 0,-1)                                                     * 
c * no15( 0,-1,-1)                                                     * 
c * no16( 1, 0,-1)                                                     * 
c * no17( 1, 1, 0)                                                     * 
c * no18(-1, 1, 0)                                                     * 
c * no19(-1,-1, 0)                                                     * 
c * no20( 1,-1, 0)                                                     * 
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8 onePlusR,onePlusS,onePlusT
      real*8 oneMinusR,oneMinusS,oneMinusT
      real*8 oneMinusRr,oneMinusSs,oneMinusTt
      real*8  h(*),hx(*),hy(*),hz(*),r,s,t
c ......................................................................
c
c ...
      onePlusR  = 1.d0 + r
      onePlusS  = 1.d0 + s
      onePlusT  = 1.d0 + t
      oneMinusR = 1.d0 - r
      oneMinusS = 1.d0 - s
      oneMinusT = 1.d0 - t
      oneMinusRr= 1.d0-r*r
      oneMinusSs= 1.d0-s*s
      oneMinusTt= 1.d0-t*t
c ......................................................................
c
c ... funcoes de interpolacao
c
      if ( afl ) then
c     Nos dos Vertices
        h(1)  = onePlusR *onePlusS *onePlusT *( r+s+t-2.d0)*0.125d0
        h(2)  = oneMinusR*onePlusS *onePlusT *(-r+s+t-2.d0)*0.125d0
        h(3)  = oneMinusR*oneMinusS*onePlusT *(-r-s+t-2.d0)*0.125d0
        h(4)  = onePlusR *oneMinusS*onePlusT *( r-s+t-2.d0)*0.125d0
        h(5)  = onePlusR *onePlusS *oneMinusT*( r+s-t-2.d0)*0.125d0
        h(6)  = oneMinusR*onePlusS *oneMinusT*(-r+s-t-2.d0)*0.125d0
        h(7)  = oneMinusR*oneMinusS*oneMinusT*(-r-s-t-2.d0)*0.125d0
        h(8)  = onePlusR *oneMinusS*oneMinusT*( r-s-t-2.d0)*0.125d0
c     Nos pontos médios das arestas      
        h(9)  = oneMinusRr*onePlusS  *onePlusT  *0.25d0
        h(10) = oneMinusR *oneMinusSs*onePlusT  *0.25d0
        h(11) = oneMinusRr*oneMinusS *onePlusT  *0.25d0
        h(12) =  onePlusR *oneMinusSs*onePlusT  *0.25d0
        h(13) = oneMinusRr*onePlusS  *oneMinusT *0.25d0
        h(14) = oneMinusR *oneMinusSs*oneMinusT *0.25d0
        h(15) = oneMinusRr*oneMinusS *oneMinusT *0.25d0
        h(16) =  onePlusR *oneMinusSs*oneMinusT *0.25d0    
        h(17) =  onePlusR *onePlusS  *oneMinusTt*0.25d0
        h(18) = oneMinusR *onePlusS  *oneMinusTt*0.25d0
        h(19) = oneMinusR *oneMinusS *oneMinusTt*0.25d0
        h(20) =  onePlusR *oneMinusS *oneMinusTt*0.25d0
c
      endif
c .....................................................................
c
c ... derivadas
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
c     Nos dos Vertices
        hx(1)  = onePlusS *onePlusT *(2.d0*r+s+t-1.0d0)*0.125d0     
        hx(2)  = onePlusS *onePlusT *(2.d0*r-s-t+1.0d0)*0.125d0
        hx(3)  = oneMinusS*onePlusT *(2.d0*r+s-t+1.0d0)*0.125d0
        hx(4)  = oneMinusS*onePlusT *(2.d0*r-s+t-1.0d0)*0.125d0
        hx(5)  = onePlusS *oneMinusT*(2.d0*r+s-t-1.0d0)*0.125d0
        hx(6)  = onePlusS *oneMinusT*(2.d0*r-s+t+1.0d0)*0.125d0
        hx(7)  = oneMinusS*oneMinusT*(2.d0*r+s+t+1.0d0)*0.125d0
        hx(8)  = oneMinusS*oneMinusT*(2.d0*r-s-t-1.0d0)*0.125d0
c     Nos pontos médios das arestas      
        hx(9)  = -r*onePlusS  *onePlusT  *0.50d0
        hx(10) =   -oneMinusSs*onePlusT  *0.25d0
        hx(11) = -r*oneMinusS *onePlusT  *0.50d0
        hx(12) =    oneMinusSs*onePlusT  *0.25d0
        hx(13) = -r*onePlusS  *oneMinusT *0.50d0
        hx(14) =   -oneMinusSs*oneMinusT *0.25d0
        hx(15) = -r*oneMinusS *oneMinusT *0.50d0
        hx(16) =    oneMinusSs*oneMinusT *0.25d0                     
        hx(17) =    onePlusS  *oneMinusTt*0.25d0                     
        hx(18) =   -onePlusS  *oneMinusTt*0.25d0                     
        hx(19) =   -oneMinusS *oneMinusTt*0.25d0                     
        hx(20) =    oneMinusS *oneMinusTt*0.25d0                     

c
c     derivadas em relacao a s :
        hy(1)  = onePlusR  *onePlusT *( r+2.d0*s+t-1.0d0)*0.125d0     
        hy(2)  = oneMinusR *onePlusT *(-r+2.d0*s+t-1.0d0)*0.125d0     
        hy(3)  = oneMinusR *onePlusT *( r+2.d0*s-t+1.0d0)*0.125d0     
        hy(4)  = oneMinusR *onePlusT *(-r+2.d0*s-t+1.0d0)*0.125d0     
        hy(5)  = onePlusR  *oneMinusT*( r+2.d0*s-t-1.0d0)*0.125d0     
        hy(6)  = oneMinusR *oneMinusT*(-r+2.d0*s-t-1.0d0)*0.125d0     
        hy(7)  = oneMinusR *oneMinusT*( r+2.d0*s+t+1.0d0)*0.125d0     
        hy(8)  = onePlusR  *oneMinusT*(-r+2.d0*s+t+1.0d0)*0.125d0     
c     Nos pontos médios das arestas      
        hy(9)  =   oneMinusRr*onePlusT  *0.25d0
        hy(10) =-s*oneMinusR *onePlusT  *0.50d0
        hy(11) =  -oneMinusRr*onePlusT  *0.25d0
        hy(12) =-s*onePlusR  *onePlusT  *0.50d0
        hy(13) =   oneMinusRr*oneMinusT *0.25d0
        hy(14) =-s*oneMinusR *oneMinusT *0.50d0
        hy(15) =  -oneMinusRr*oneMinusT *0.25d0
        hy(16) =-s*onePlusR  *oneMinusT *0.50d0
        hy(17) =   onePlusR  *oneMinusTt*0.25d0
        hy(18) =   oneMinusR *oneMinusTt*0.25d0
        hy(19) =  -oneMinusR *oneMinusTt*0.25d0
        hy(20) =  -onePlusR  *oneMinusTt*0.25d0
c
c
c     derivadas em relacao a t :
c
        hz(1)  = onePlusR *onePlusS  *( r+s+2.d0*t-1.0d0)*0.125d0     
        hz(2)  = oneMinusR*onePlusS  *(-r+s+2.d0*t-1.0d0)*0.125d0     
        hz(3)  = oneMinusR*oneMinusS *(-r+s+2.d0*t-1.0d0)*0.125d0     
        hz(4)  = onePlusR *oneMinusS *(-r-s+2.d0*t-1.0d0)*0.125d0     
        hz(5)  = onePlusR *onePlusS  *(-r-s+2.d0*t+1.0d0)*0.125d0     
        hz(6)  = oneMinusR*onePlusS  *( r-s+2.d0*t+1.0d0)*0.125d0     
        hz(7)  = oneMinusR*oneMinusS *( r+s+2.d0*t+1.0d0)*0.125d0     
        hz(8)  = oneMinusR*oneMinusS *(-r+s+2.d0*t+1.0d0)*0.125d0     
c     Nos pontos médios das arestas      
        hz(9)  = oneMinusRr *onePlusS  *0.25d0
        hz(10) = oneMinusR  *oneMinusSs*0.25d0
        hz(11) = oneMinusRr *oneMinusS *0.25d0
        hz(12) =  onePlusR  *oneMinusSs*0.25d0
        hz(13) = -oneMinusRr*onePlusS  *0.25d0
        hz(14) = -oneMinusR *oneMinusSs*0.25d0
        hz(15) = -oneMinusRr*oneMinusS *0.25d0
        hz(16) =  -onePlusR *oneMinusSs*0.25d0    
        hz(17) =-t*onePlusR *onePlusS  *0.50d0
        hz(18) =-t*oneMinusR*onePlusS  *0.50d0
        hz(19) =-t*oneMinusR*oneMinusS *0.50d0
        hz(20) =-t*onePlusR *oneMinusS *0.50d0
      endif
c ......................................................................      
      return
      end
c **********************************************************************
c
c     subroutine sfhexa20(h,hx,hy,hz,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    27/01/11        *
c *   SFHEXA20:                                                        *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t), do hexaedro quadratico de 20 nos.              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s,t - coordenadas naturais                                   *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(8)    - funcoes de interolocao no ponto (r,s,t)              *
c *     hx(8)   - derivadas de h em relacao a r                        *
c *     hy(8)   - derivadas de h em relacao a s                        *
c *     hz(8)   - derivadas de h em relacao a t                        *
c *                                                                    *
c **********************************************************************
c     implicit none
c     logical afl,bfl
c     real*8  h(*),hx(*),hy(*),hz(*),r,s,t
c ......................................................................
c     if ( afl ) then
c
c     funcoes de interpolacao
c
c     Nos dos Vertices
c     h(1)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) *(r+s+t-2)/ 8.d0
c     h(2)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t) *(-r+s+t-2)/ 8.d0
c     h(3)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) *(-r-s+t-2)/ 8.d0
c     h(4)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t) *(r-s+t-2)/ 8.d0
c     h(5)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t) *(r+s-t-2)/ 8.d0
c     h(6)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t) *(-r+s-t-2)/ 8.d0
c     h(7)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) *(-r-s-t-2)/ 8.d0
c     h(8)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t) *(r-s-t-2)/ 8.d0
c     Nos pontos médios das arestas      
c     h(9)  =  ( 1.d0-r**2 ) * ( 1.d0+s ) * (1.d0+t) / 4.d0
c     h(10)  = ( 1.d0-r ) * ( 1.d0-s**2 ) * (1.d0+t) / 4.d0
c     h(11)  = ( 1.d0-r**2 ) * ( 1.d0-s ) * (1.d0+t) / 4.d0
c     h(12)  = ( 1.d0+r ) * ( 1.d0-s**2 ) * (1.d0+t) / 4.d0
c     h(13)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t**2) / 4.d0
c     h(14)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t**2) / 4.d0
c     h(15)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t**2) / 4.d0
c     h(16)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t**2) / 4.d0
c     h(17)  =  ( 1.d0-r**2 ) * ( 1.d0+s ) * (1.d0-t) / 4.d0
c     h(18)  =  ( 1.d0-r ) * ( 1.d0-s**2 ) * (1.d0-t) / 4.d0
c     h(19)  =  ( 1.d0-r**2 ) * ( 1.d0-s ) * (1.d0-t) / 4.d0
c     h(20)  =  ( 1.d0+r ) * ( 1.d0-s**2 ) * (1.d0-t) / 4.d0    
c
c     endif
c     if ( bfl ) then
c
c     derivadas em relacao a r :
c
c     hx(1) = ( 1.d0+s ) * (1.d0+t) *(r+s+t-2)/ 8.d0+
c    . ( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0     
c     hx(2) = -( 1.d0+s ) * (1.d0+t) *(-r+s+t-2)/ 8.d0-
c    .( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
c     hx(3) = - ( 1.d0-s ) * (1.d0+t) *(-r-s+t-2)/ 8.d0-
c    .( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
c     hx(4) =  ( 1.d0-s ) * (1.d0+t) *(r-s+t-2)/ 8.d0+
c    .( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
c     hx(5) = ( 1.d0+s ) * (1.d0-t) *(r+s-t-2)/ 8.d0+
c    . ( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
c     hx(6) = -( 1.d0+s ) * (1.d0-t) *(-r+s-t-2)/ 8.d0-
c    .( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t)/ 8.d0
c     hx(7) = -( 1.d0-s ) * (1.d0-t) *(-r-s-t-2)/ 8.d0-
c    .( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
c     hx(8) =  ( 1.d0-s ) * (1.d0-t) *(r-s-t-2)/ 8.d0+
c    .( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t)/ 8.d0
c     hx(9)  =  -2*r * ( 1.d0+s ) * (1.d0+t) / 4.d0
c     hx(10)  = -( 1.d0-s**2 ) * (1.d0+t) / 4.d0
c     hx(11)  = -2*r * ( 1.d0-s ) * (1.d0+t) / 4.d0
c     hx(12)  =  ( 1.d0-s**2 ) * (1.d0+t) / 4.d0
c     hx(13)  =  ( 1.d0+s ) * (1.d0-t**2) / 4.d0
c     hx(14)  = -( 1.d0+s ) * (1.d0-t**2) / 4.d0
c     hx(15)  = -( 1.d0-s ) * (1.d0-t**2) / 4.d0
c     hx(16)  =  ( 1.d0-s ) * (1.d0-t**2) / 4.d0
c     hx(17)  = -2*r * ( 1.d0+s ) * (1.d0-t) / 4.d0
c     hx(18)  = -( 1.d0-s**2 ) * (1.d0-t) / 4.d0
c     hx(19)  = -2*r * ( 1.d0-s ) * (1.d0-t) / 4.d0
c     hx(20)  =  ( 1.d0-s**2 ) * (1.d0-t) / 4.d0  
c
c     derivadas em relacao a s :
c
c     hy(1)  = ( 1.d0+r ) *  (1.d0+t) *(r+s+t-2)/ 8.d0+
c    . ( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
c     hy(2)  = ( 1.d0-r ) * (1.d0+t) *(-r+s+t-2)/ 8.d0+
c    .( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
c     hy(3)  = -( 1.d0-r ) *  (1.d0+t) *(-r-s+t-2)/ 8.d0-
c    .( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
c     hy(4)  = -( 1.d0+r ) *  (1.d0+t) *(r-s+t-2)/ 8.d0-
c    .( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t)/ 8.d0
c     hy(5)  = ( 1.d0+r ) * (1.d0-t) *(r+s-t-2)/ 8.d0+
c    .( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
c     hy(6)  = ( 1.d0-r ) * (1.d0-t) *(-r+s-t-2)/ 8.d0+
c    .( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
c     hy(7)  = -( 1.d0-r ) * (1.d0-t) *(-r-s-t-2)/ 8.d0-
c    .( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
c     hy(8)  = -( 1.d0+r ) *  (1.d0-t) *(r-s-t-2)/ 8.d0-
c    .( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
c     hy(9)  =  ( 1.d0-r**2 ) * (1.d0+t) / 4.d0
c     hy(10)  = ( 1.d0-r ) * ( -2*s ) * (1.d0+t) / 4.d0
c     hy(11)  = -( 1.d0-r**2 )  * (1.d0+t) / 4.d0
c     hy(12)  = ( 1.d0+r ) * (-2*s ) * (1.d0+t) / 4.d0
c     hy(13)  = ( 1.d0+r )  * (1.d0-t**2) / 4.d0
c     hy(14)  = ( 1.d0-r )  * (1.d0-t**2) / 4.d0
c     hy(15)  = -( 1.d0-r )  * (1.d0-t**2) / 4.d0
c     hy(16)  = -( 1.d0+r )  * (1.d0-t**2) / 4.d0
c     hy(17)  =  ( 1.d0-r**2 ) * (1.d0-t) / 4.d0
c     hy(18)  =  ( 1.d0-r ) * ( -2*s ) * (1.d0-t) / 4.d0
c     hy(19)  = -( 1.d0-r**2 )  * (1.d0-t) / 4.d0
c     hy(20)  =  ( 1.d0+r ) * ( -2*s ) * (1.d0-t) / 4.d0  
c
c     derivadas em relacao a t :
c
c     hz(1)  = ( 1.d0+r ) * ( 1.d0+s ) * (r+s+t-2)/ 8.d0+
c    .( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
c     hz(2)  = ( 1.d0-r ) * ( 1.d0+s ) *(-r+s+t-2)/ 8.d0+
c    .( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t)/ 8.d0
c     hz(3)  = ( 1.d0-r ) * ( 1.d0-s ) *(-r-s+t-2)/ 8.d0+
c    . ( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
c     hz(4)  = ( 1.d0+r ) * ( 1.d0-s ) * (r-s+t-2)/ 8.d0+
c    .( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
c     hz(5)  = -( 1.d0+r ) * ( 1.d0+s ) * (r+s-t-2)/ 8.d0-
c    .( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t)/ 8.d0
c     hz(6)  = -( 1.d0-r ) * ( 1.d0+s ) * (-r+s-t-2)/ 8.d0-
c    .( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t)/ 8.d0
c     hz(7)  = -( 1.d0-r ) * ( 1.d0-s ) * (-r-s-t-2)/ 8.d0-
c    .( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
c     hz(8)  = -( 1.d0+r ) * ( 1.d0-s )  *(r-s-t-2)/ 8.d0-
c    .( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
c     hz(9)  =  ( 1.d0-r**2 ) * ( 1.d0+s )  / 4.d0
c     hz(10)  = ( 1.d0-r ) * ( 1.d0-s**2 )  / 4.d0
c     hz(11)  = ( 1.d0-r**2 ) * ( 1.d0-s )  / 4.d0
c     hz(12)  = ( 1.d0+r ) * ( 1.d0-s**2 )  / 4.d0
c     hz(13)  = ( 1.d0+r ) * ( 1.d0+s ) * (-2*t) / 4.d0
c     hz(14)  = ( 1.d0-r ) * ( 1.d0+s ) * (-2*t) / 4.d0
c     hz(15)  = ( 1.d0-r ) * ( 1.d0-s ) * (-2*t) / 4.d0
c     hz(16)  = ( 1.d0+r ) * ( 1.d0-s ) * (-2*t) / 4.d0
c     hz(17)  =  -( 1.d0-r**2 ) * ( 1.d0+s )  / 4.d0
c     hz(18)  =  -( 1.d0-r ) * ( 1.d0-s**2 )  / 4.d0
c     hz(19)  =  -( 1.d0-r**2 ) * ( 1.d0-s )  / 4.d0
c     hz(20)  =  -( 1.d0+r ) * ( 1.d0-s**2 )  / 4.d0 
c     endif
c ......................................................................      
c     return
c     end
c **********************************************************************
      subroutine sfwedge(h,hx,hy,hz,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    26/03/03        *
c *   SFWEDGE:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t) para o elemento cunha                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s,t - coordenadas naturais                                   *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(6)    - funcoes de interolocao no ponto (r,s,t)              *
c *     hr(6)   - derivadas de h em relacao a r                        *
c *     hs(6)   - derivadas de h em relacao a s                        *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  h(*),hx(*),hy(*),hz(*),r,s,t
      logical afl,bfl
c ......................................................................
      if (afl) then
c
c ...... Funcoes de interpolacao:
c
         h(1) = r * (1.d0+t) * 0.5d0
         h(2) = s * (1.d0+t) * 0.5d0
         h(3) =(1.d0-r-s) * (1.d0+t) * 0.5d0
         h(4) = r * (1.d0-t) * 0.5d0
         h(5) = s * (1.d0-t) * 0.5d0
         h(6) =(1.d0-r-s) * (1.d0-t) * 0.5d0
      endif
      if (bfl) then
c
c ...... Derivadas em relacao a r:
c
         hx(1) = (1.d0+t) * 0.5d0
         hx(2) =  0.d0
         hx(3) =-(1.d0+t) * 0.5d0
         hx(4) = (1.d0-t) * 0.5d0
         hx(5) =  0.d0
         hx(6) = (t-1.d0) * 0.5d0
c
c ...... Derivadas em relacao a s:
c
         hy(1) =  0.d0
         hy(2) = (1.d0+t) * 0.5d0
         hy(3) =-(1.d0+t) * 0.5d0
         hy(4) =  0.d0
         hy(5) = (1.d0-t) * 0.5d0
         hy(6) = (t-1.d0) * 0.5d0
c
c ...... Derivadas em relacao a t:
c
         hz(1) =  r * 0.5d0
         hz(2) =  s * 0.5d0
         hz(3) = (1.d0-r-s) * 0.5d0
         hz(4) = -r * 0.5d0
         hz(5) = -s * 0.5d0
         hz(6) = (r+s-1.d0) * 0.5d0
      endif
c ......................................................................      
      return
      end
      subroutine jacob2d(hx,hy,xj,xji,x,det,nen,ndm,nel)
c **********************************************************************
c *                                                                    *
c *                                                    20/03/03        *
c *                                                                    *
c *   JACOB2D: determinante do Jacobiano no ponto (r,s)                *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     hx(nen)   - derivadas de h em relacao a r                      *
c *     hy(nen)   - derivadas de h em relacao a s                      *
c *     x(ndm,nen)- coordenadas nodais do elemento                     *
c *     nen       - numero de nos do elemento                          *
c *     nel       - numero do elemento                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     hx(nen)  - derivadas de h em relacao a x                       *
c *     hy(nen)  - derivadas de h em relacao a y                       *
c *     xj(2,2)  - matriz jacobiana no ponto (r,s,t)                   *
c *     xji(2,2) - inversa de xj                                       *
c *     det      - determinante da matriz jacobiana no ponto (r,s)     *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,ndm,nel,j,k
      real*8  ZERO
      real*8  hx(*),hy(*),xj(2,*),xji(2,*),x(ndm,*)
      real*8  det,hxk,hyk
      parameter (ZERO = 1.d-14)
c ......................................................................
c
c ... Matriz Jacobiana:
c
      do 200 j = 1 , 2
         xj(1,j) = 0.d0
         xj(2,j) = 0.d0
         do 100 k = 1 , nen
            xj(1,j) = xj(1,j) + hx(k) * x(j,k)
            xj(2,j) = xj(2,j) + hy(k) * x(j,k)
  100    continue
  200 continue
c
c ... Determinante da matriz Jacobiana:  
c
c ......................................................................              
      det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
      if (det .le. ZERO) then        
         print*,'*** Subrotina ELMT: determinante <= 0 ',nel
         stop
      endif
c ......................................................................                    
c
c ... Inversa da matriz Jacobiana:  
c
      xji(1,1) =  xj(2,2) / det
      xji(1,2) = -xj(1,2) / det
      xji(2,1) = -xj(2,1) / det
      xji(2,2) =  xj(1,1) / det
c
c ... Derivadas das funcoes de interpolacao:
c
      do 300 k = 1, nen
         hxk = hx(k)
         hyk = hy(k)
         hx(k) = xji(1,1)*hxk + xji(1,2)*hyk
         hy(k) = xji(2,1)*hxk + xji(2,2)*hyk
  300 continue
c ......................................................................                
      return
      end
      subroutine jacob3d(hx,hy,hz,xj,xji,x,det,nen,nel,afl)
c **********************************************************************
c *                                                                    *
c *                                                    10/10/01        *
c *                                                                    *
c *   JACOB3D: determinante do Jacobiano no ponto (r,s,t)              *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     hx(nen)   - derivadas de h em relacao a r                      *
c *     hy(nen)   - derivadas de h em relacao a s                      *
c *     hz(nen)   - derivadas de h em relacao a s                      *
c *     x(ndm,nen)- coordenadas nodais do elemento                     *
c *     nen       - numero de nos do elemento                          *
c *     nel       - numero do elemento                                 *
c *     afl       - true - calcula a inversa da matriz Jacobiana       *
c *                 false- calcula somente as derivadas das funcoes de *
c *                 interpolacao                                       *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     hx(nen)  - derivadas de h em relacao a x                       *
c *     hy(nen)  - derivadas de h em relacao a y                       *
c *     hz(nen)  - derivadas de h em relacao a z                       *
c *     xj(3,3)  - matriz jacobiana no ponto (r,s,t)                   *
c *     xji(3,3) - inversa de xj                                       *
c *     det      - determinante da matriz jacobiana no ponto (r,s,t)   *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,nel,j,k
      real*8  hx(*),hy(*),hz(*),xj(3,*),xji(3,*),x(3,*)
      real*8  det,hxk,hyk,hzk
      logical afl
      real*8  ZERO
      parameter (ZERO = 1.d-14)
c ......................................................................
c
c
      if (afl) then 
c ... Matria Jacobiana:
c
        do 200 j = 1 , 3
          xj(1,j) = 0.
          xj(2,j) = 0.
          xj(3,j) = 0.
          do 100 k = 1 , nen
            xj(1,j) = xj(1,j) + hx(k) * x(j,k)
            xj(2,j) = xj(2,j) + hy(k) * x(j,k)
            xj(3,j) = xj(3,j) + hz(k) * x(j,k)
  100     continue
  200   continue
c
c ... Determinante da matriz Jacobiana:  
c
        det = xj(1,1)*xj(2,2)*xj(3,3) + xj(1,2)*xj(2,3)*xj(3,1)+xj(1,3)*
     .        xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2)*xj(1,3)-xj(2,1)*xj(1,2)*
     .        xj(3,3) - xj(1,1)*xj(3,2)*xj(2,3)
c ......................................................................                   
        if (det .le. ZERO) then
          print*,'*** Subrotina ELMT__: determinante <= 0 ',nel
          stop
        endif
c ......................................................................                    
c     
c ... Inversa da matriz Jacobiana:  
c
        xji(1,1) =  ( xj(2,2) * xj(3,3) - xj(2,3) * xj(3,2) ) / det
        xji(2,1) = -( xj(2,1) * xj(3,3) - xj(2,3) * xj(3,1) ) / det
        xji(3,1) =  ( xj(2,1) * xj(3,2) - xj(2,2) * xj(3,1) ) / det
        xji(1,2) = -( xj(1,2) * xj(3,3) - xj(1,3) * xj(3,2) ) / det
        xji(2,2) =  ( xj(1,1) * xj(3,3) - xj(1,3) * xj(3,1) ) / det
        xji(3,2) = -( xj(1,1) * xj(3,2) - xj(1,2) * xj(3,1) ) / det
        xji(1,3) =  ( xj(1,2) * xj(2,3) - xj(1,3) * xj(2,2) ) / det
        xji(2,3) = -( xj(1,1) * xj(2,3) - xj(1,3) * xj(2,1) ) / det
        xji(3,3) =  ( xj(1,1) * xj(2,2) - xj(1,2) * xj(2,1) ) / det
c ......................................................................          
      endif
c ......................................................................          
c
c ... Derivadas das funcoes de interpolacao:
c
      do 300 k = 1, nen
         hxk = hx(k)
         hyk = hy(k)
         hzk = hz(k)
         hx(k) = xji(1,1)*hxk + xji(1,2)*hyk + xji(1,3)*hzk
         hy(k) = xji(2,1)*hxk + xji(2,2)*hyk + xji(2,3)*hzk
         hz(k) = xji(3,1)*hxk + xji(3,2)*hyk + xji(3,3)*hzk
  300 continue
c ......................................................................              
      return
      end
      subroutine rotmatrix(x,r)
c **********************************************************************
c *                                                                    *
c *                                                    26/03/03        *
c *   ROTACAO:                                                         *
c *   -------                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(3,*),r(3,3)
      real*8 v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,v1,v2,v3
c ......................................................................
      v1x = x(1,2)-x(1,1)
      v1y = x(2,2)-x(2,1)
      v1z = x(3,2)-x(3,1)
      v2x = x(1,3)-x(1,1)
      v2y = x(2,3)-x(2,1)
      v2z = x(3,3)-x(3,1)
      v3x = v1y*v2z - v1z*v2y
      v3y = v1z*v2x - v1x*v2z
      v3z = v1x*v2y - v1y*v2x
      v2x = v1z*v3y - v1y*v3z
      v2y = v1x*v3z - v1z*v3x
      v2z = v1y*v3x - v1x*v3y
      v1 = dsqrt(v1x*v1x+v1y*v1y+v1z*v1z)
      v2 = dsqrt(v2x*v2x+v2y*v2y+v2z*v2z)
      v3 = dsqrt(v3x*v3x+v3y*v3y+v3z*v3z)
      r(1,1) = v1x/v1
      r(1,2) = v1y/v1
      r(1,3) = v1z/v1
      r(2,1) = v2x/v2
      r(2,2) = v2y/v2
      r(2,3) = v2z/v2
      r(3,1) = v3x/v3
      r(3,2) = v3y/v3
      r(3,3) = v3z/v3
c ......................................................................      
      return
      end
      subroutine rotmatrix2(x1,x2,x3,r)
c **********************************************************************
c *                                                                    *
c *                                                    26/03/03        *
c *   ROTACAO:                                                         *
c *   -------                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x1(3),x2(3),x3(3),r(3,3)
      real*8 v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,v1,v2,v3
c ......................................................................
      
      v1x = x2(1)-x1(1)
      v1y = x2(2)-x1(2)
      v1z = x2(3)-x1(3)
      v2x =x3(1)-x1(1)
      v2y = x3(2)-x1(2)
      v2z = x3(3)-x1(3)
      v3x = v1y*v2z - v1z*v2y
      v3y = v1z*v2x - v1x*v2z
      v3z = v1x*v2y - v1y*v2x
      v2x = v1z*v3y - v1y*v3z
      v2y = v1x*v3z - v1z*v3x
      v2z = v1y*v3x - v1x*v3y
      v1 = dsqrt(v1x*v1x+v1y*v1y+v1z*v1z)
      v2 = dsqrt(v2x*v2x+v2y*v2y+v2z*v2z)
      v3 = dsqrt(v3x*v3x+v3y*v3y+v3z*v3z)
      r(1,1) = v1x/v1
      r(1,2) = v1y/v1
      r(1,3) = v1z/v1
      r(2,1) = v2x/v2
      r(2,2) = v2y/v2
      r(2,3) = v2z/v2
      r(3,1) = v3x/v3
      r(3,2) = v3y/v3
      r(3,3) = v3z/v3
c ......................................................................      
      return
      end
      subroutine irotmatrix(x,r)
c **********************************************************************
c *                                                                    *
c *                                                    21/11/2010      *
c *   ROTACAO:                                                         *
c *   calcula matriz de rotacao sistema local para 2D                  *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(2,*),r(2,2)
      real*8 v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,v1,v2,v3
c ......................................................................
      v1x = x(1,2)-x(1,1)
      v1y = x(2,2)-x(2,1)
      v2x = -v1y
      v2y = v1x
      v1 = dsqrt(v1x*v1x+v1y*v1y)
      v2 = dsqrt(v2x*v2x+v2y*v2y)
      r(1,1) = v1x/v1
      r(1,2) = v1y/v1
      r(2,1) = v2x/v2
      r(2,2) = v2y/v2
c ......................................................................      
      return
      end
      subroutine rotate(x,r,xl,transp)
c **********************************************************************
c *                                                                    *
c *                                                    26/03/03        *
c *   ROTACAO:                                                         *
c *   -------                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(3),r(3,3),xl(3),x1,x2,x3
      logical transp
c ......................................................................
      if(transp) then
         x1 = r(1,1)*x(1)+r(2,1)*x(2)+r(3,1)*x(3)
         x2 = r(1,2)*x(1)+r(2,2)*x(2)+r(3,2)*x(3)
         x3 = r(1,3)*x(1)+r(2,3)*x(2)+r(3,3)*x(3)      
      else
         x1 = r(1,1)*x(1)+r(1,2)*x(2)+r(1,3)*x(3)
         x2 = r(2,1)*x(1)+r(2,2)*x(2)+r(2,3)*x(3)
         x3 = r(3,1)*x(1)+r(3,2)*x(2)+r(3,3)*x(3)
      endif
      xl(1) = x1
      xl(2) = x2
      xl(3) = x3
c ......................................................................      
      return
      end
      subroutine irotate(x,r,xl,transp)
c **********************************************************************
c *                                                                    *
c *                                                   21/11/2010       *
c *   ROTACAO:                                                         *
c *   Rotaciona as coordenadas                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(2),r(2,2),xl(2),x1,x2
      logical transp
c ......................................................................
      if(transp) then
         x1 = r(1,1)*x(1)+r(2,1)*x(2)
         x2 = r(1,2)*x(1)+r(2,2)*x(2)
      else
         x1 = r(1,1)*x(1)+r(1,2)*x(2)
         x2 = r(2,1)*x(1)+r(2,2)*x(2)
      endif
      xl(1) = x1
      xl(2) = x2
c ......................................................................      
      return
      end
      subroutine irotate2(x,r,xl,inv)
c **********************************************************************
c *                                                                    *
c *                                                   21/11/2010       *
c *   ROTACAO:                                                         *
c *   Rotaciona tensor de segunda ordem                                *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(3),r(2,2),xl(3),rt11,rt12,rt21,rt22,x1,x2,x3
      logical inv
c ......................................................................
      if(inv) then
         rt11 = r(1,1)*x(1)+r(2,1)*x(3)
         rt12 = r(1,1)*x(3)+r(2,1)*x(2)
         rt21 = r(1,2)*x(1)+r(2,2)*x(3)
         rt22 = r(1,2)*x(3)+r(2,2)*x(2)
         x1 = rt11*r(1,1)+rt12*r(2,1)
         x2 = rt21*r(1,2)+rt22*r(2,2)
         x3 = rt11*r(1,2)+rt12*r(2,2)
      else
         rt11 = r(1,1)*x(1)+r(1,2)*x(3)
         rt12 = r(1,1)*x(3)+r(1,2)*x(2)
         rt21 = r(2,1)*x(1)+r(2,2)*x(3)
         rt22 = r(2,1)*x(3)+r(2,2)*x(2)
         x1 = rt11*r(1,1)+rt12*r(1,2)
         x2 = rt21*r(2,1)+rt22*r(2,2)
         x3 = rt11*r(2,1)+rt12*r(2,2)
      endif
      xl(1) = x1
      xl(2) = x2
      xl(3) = x3
c ......................................................................      
      return
      end
      subroutine irotatek(s,r) 
c **********************************************************************
c *                                                                    *
c *                                                   21/11/2010       *
c *   ROTACAO:                                                         *
c *   Rotaciona matriz de rigidez                                      *
c *                                                                    *
c **********************************************************************
      implicit none 
      real*8 s(8,8),r(2,2),sl(8,8) 
      integer i,j 
c ......................................................................
      do j = 1,4 
         do i =1,4 
            sl(2*j-1,2*i-1)=r(1,1)*s(2*j-1,2*i-1)*r(1,1)+ 
     .                  r(2,1)*s(2*j,2*i)*r(2,1) 
            sl(2*j-1,2*i)=r(1,1)*s(2*j-1,2*i-1)*r(1,2)+ 
     .                  r(2,1)*s(2*j,2*i)*r(2,2) 
            sl(2*j,2*i-1)=r(1,2)*s(2*j-1,2*i-1)*r(1,1)+ 
     .                  r(2,2)*s(2*j,2*i)*r(2,1) 
            sl(2*j,2*i)=r(1,2)*s(2*j-1,2*i-1)*r(1,2)+ 
     .                  r(2,2)*s(2*j,2*i)*r(2,2) 
      enddo 
      enddo 
      do i = 1,8 
         do j = 1,8 
            s(i,j) = sl(i,j) 
         enddo 
      enddo 
c ...................................................................... 
      return 
      end
      subroutine irotatek3d(s,r,nst) 
c **********************************************************************
c *                                                                    *
c *                                                   21/11/2010       *
c *   ROTACAO:                                                         *
c *   Rotaciona matriz de rigidez 3d                                   *
c *                                                                    *
c **********************************************************************
      implicit none 
      integer i,j,k,l,i1,i2,i3,j1,j2,j3,nst
      real*8 s(nst,nst),r(3,3),sl(nst,nst),aux(nst,nst),raux(nst,nst)
c ......................................................................
      l = nst/3
      do i = 1, nst
         do j = 1, nst
            sl(i,j) = 0.d0
            aux(i,j) = 0.d0
            raux(i,j) = 0.d0
         enddo
      enddo
c
c      do i = 1, l 
c      i1 = 3*(i-1)+1
c      i2 = 3*(i-1)+2
c      i3 = 3*(i-1)+3
c         do j =1, l 
c            j1 = 3*(j-1)+1
c            j2 = 3*(j-1)+2
c            j3 = 3*(j-1)+3
c            do k = 1, 18
c               
c               aux(i1,j1) = aux(i1,j1) + r(1,k)*s(k,j1)
c               aux(i2,j2) = aux(i2,j2) + r(2,k)*s(k,j2)
c               aux(i3,j3) = aux(i3,j3) + r(3,k)*s(k,j3) 
c            enddo  
c         enddo 
c      enddo 
c      do i = 1, l 
c      i1 = 3*(i-1)+1
c      i2 = 3*(i-1)+2
c      i3 = 3*(i-1)+3
c         do j =1, l 
c            j1 = 3*(j-1)+1
c            j2 = 3*(j-1)+2
c            j3 = 3*(j-1)+3
c            do k = 1, 3
c               sl(i1,j1) = sl(i1,j1) + aux(i1,k)*r(1,k) 
c               sl(i2,j2) = sl(i2,j2) + aux(i2,k)*r(2,k) 
c               sl(i3,j3) = sl(i3,j3) + aux(i3,k)*r(3,k) 
c            enddo  
c         enddo 
c      enddo 
      do i = 1, l
         i1 = 3*(i-1)+1
         i2 = 3*(i-1)+2
         i3 = 3*(i-1)+3
            raux(i1,i1) = r(1,1)
            raux(i1,i2) = r(1,2)
            raux(i1,i3) = r(1,3)
            raux(i2,i1) = r(2,1)
            raux(i2,i2) = r(2,2)
            raux(i2,i3) = r(2,3)
            raux(i3,i1) = r(3,1)
            raux(i3,i2) = r(3,2)
            raux(i3,i3) = r(3,3)            
      enddo
      do i = 1, nst
         do j = 1, nst
            do k  = 1, nst
               aux(i,j) = aux(i,j) + raux(k,i)*s(k,j)
            enddo
         enddo
      enddo
      do i = 1, nst
         do j = 1, nst
            do k  = 1, nst
               sl(i,j) = sl(i,j) + aux(i,k)*raux(k,j)
            enddo
         enddo
      enddo      
      do i = 1,nst 
         do j = 1,nst 
            s(i,j) = sl(i,j) 
         enddo 
      enddo 
c ...................................................................... 
      return 
      end
c **********************************************************************
c
c **********************************************************************
      block data
c **********************************************************************
c *                                                                    *
c *                                                    25/01/89        *
c *   - Descricao:                                                     *
c *                                                                    *
c *     Inicializa os pontos de integracao.                            *
c *                                                                    *
c *     Gauss-Legendre:                                                *
c *     --------------                                                 *
c *     do lz = 1, nint                                                *
c *       ti = pg(lz,nint)                                             *
c *       do ly = 1, nint                                              *
c *          si = pg(ly,nint)                                          *
c *          do lx = 1, nint                                           *
c *             ri = pg(lx,nint)                                       *
c *             wt = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det           *
c *          enddo                                                     *
c *       enddo                                                        *
c *     enddo                                                          *
c *                                                                    *
c *     Integracao em triangulos:                                      *
c *     -------------------------                                      *
c *     nint = npint(igrau)                                            *
c *     do lx = 1, nint                                                *
c *        ri = pri(lx,igrau)                                          *
c *        si = psi(lx,igrau)                                          *
c *        ti = 1-ri-si                                                *
c *        call dgh(node,xl,ri,si,ti,det,noel,hfl)                     *
c *        wt = wf(lx,igrau) * 0.5d0 * det * thic                      *
c *     end do                                                         *
c *                                                                    *
c *     Integracao em tetraedros:                                      *
c *     -------------------------                                      *
c *     nint = npint4(grau)                                            *
c *     do lx = 1, nint                                                *
c *        ri = pri4(lx,igrau)                                         *
c *        si = psi4(lx,igrau)                                         *
c *        ti = pti4(lx,igrau)                                         *
c *        ui = 1 - ri - si - ti                                       *
c *        call dgh(node,xl,ri,si,ti,ui,det,noel,hfl)                  *
c *        wt = wf4(lx,igrau)                                          *
c *     end do                                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      common /gauss/ pg, wg
      common /pint / pri,psi,wf,npint
      common /pint4/ pri4,psi4,pti4,wf4,npint4      
      real*8  pg(10,10), wg(10,10)
      real*8  pri(12,5),psi(12,5),wf(12,5)
      real*8  pri4(4,2),psi4(4,2),pti4(4,2),wf4(4,2)
      integer npint(5),npint4(2)
c ======================================================================
c
c ... Pontos de integracao de Gauss-Legendre:
c
      data pg /      0.d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     .               0.d0,    0.d0, 0.d0,
     .-.577350269189626d0,   .577350269189626d0,    0.d0, 0.d0,  0.d0,
     .  0.d0, 0.d0,  0.d0,    0.d0, 0.d0,
     .-.774596669241483d0,    0.d0,                .774596669241483d0,
     .  0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,    0.d0,
     .-.861136311594053d0,  -.339981043584856d0,   .339981043584856d0,
     . .861136311594053d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     .-.906179845938664d0,  -.538469310105683d0,    0.d0,
     . .538469310105683d0,   .906179845938664d0,    0.d0, 0.d0,  0.d0,
     .  0.d0, 0.d0,
     .-.932469514203152d0,  -.661209386466265d0,  -.238619186083197d0,
     . .238619186083197d0,   .661209386466265d0,   .932469514203152d0,
     .  0.d0, 0.d0,  0.d0,    0.d0,
     .-.949107912342759d0,  -.741531185599394d0,  -.405845151377397d0,
     .  0.d0,                .405845151377397d0,   .741531185599394d0,
     . .949107912342759d0,    0.d0, 0.d0,  0.d0,
     .-.960289856497536d0,  -.796666477413627d0,  -.525532409916329d0,
     .-.183434642495650d0,   .183434642495650d0,   .525532409916329d0,
     . .796666477413627d0,   .960289856497536d0,    0.d0, 0.d0,
     .-.968160239507626d0,  -.836031107326636d0,  -.613371432700590d0,
     .-.324253423403809d0,    0.d0,                .324253423403809d0,
     . .613371432700590d0,   .836031107326636d0,   .968160239507626d0,
     .  0.d0,
     .-.973906528517172d0,  -.865063366688985d0,  -.679409568299024d0,
     .-.433395394129247d0,  -.148874338981631d0,   .148874338981631d0,
     . .433395394129247d0,   .679409568299024d0,   .865063366688985d0,
     . .973906528517172d0 /
c
c ... Fatores de ponderacao de Gauss-Legendre:
c
      data wg /      2.d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     .               0.d0,    0.d0, 0.d0,
     .  1.d0, 1.d0,  0.d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     .  0.d0,
     . .555555555555556d0,   .888888888888889d0,   .555555555555556d0,
     .  0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,    0.d0,
     . .347854845137454d0,   .652145154862546d0,   .652145154862546d0,
     . .347854845137454d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     . .236926885056189d0,   .478628670499366d0,   .568888888888889d0,
     . .478628670499366d0,   .236926885056189d0,    0.d0, 0.d0,  0.d0,
     .  0.d0, 0.d0,
     . .171324492379170d0,   .360761573048139d0,   .467913934572691d0,
     . .467913934572691d0,   .360761573048139d0,   .171324492379170d0,
     .  0.d0, 0.d0,  0.d0,    0.d0,
     . .129484966168870d0,   .279705391489277d0,   .381830050505119d0,
     . .417959183673469d0,   .381830050505119d0,   .279705391489277d0,
     . .129484966168870d0,    0.d0, 0.d0,  0.d0,
     . .101228532690376d0,   .222381034453374d0,   .313706645877887d0,
     . .362683783378362d0,   .362683783378362d0,   .313706645877887d0,
     . .222381034453374d0,   .101228536290376d0,    0.d0, 0.d0,
     . .081274388361574d0,   .180648160694857d0,   .260610696402935d0,
     . .312347077040003d0,   .330239355001260d0,   .312347077040003d0,
     . .260610696402935d0,   .180648160694857d0,   .081274388361574d0,
     .  0.d0,
     . .066671344308688d0,   .149451349150581d0,   .219086362515982d0,
     . .269266719309996d0,   .295524224714753d0,   .295524224714753d0,
     . .269266719309996d0,   .219086362515982d0,   .149451349150581d0,
     . .066671344308688d0 /
c ======================================================================     
c
c ... Pontos de integracao em triangulos:
c ... Numero de pontos de integracao (grau = 1,..., 5):
      data npint /1,3,4,7,12/      
c ... Coordenadas r:
      data pri /.333333333333333d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     .                      .5d0,.0,.5d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     .    .333333333333333d0,.6d0,.2d0,.2d0,.0,.0,.0,.0,.0,.0,.0,.0,
     ..333333333333333d0,.059715871789770d0,.470142064105115d0,
     ..470142064105115d0,.797426985353087d0,.101286507323456d0,
     ..101286507323456d0,.0,.0,.0,.0,.0,
     ..873821971016996d0,.063089014491502d0,.063089014491502d0,
     ..501426509658179d0,.249286745170910d0,.249286745170910d0,
     ..636502499121399d0,.636502499121399d0,.310352451033785d0,
     ..310352451033785d0,.053145049844816d0,.053145049844816d0/
c ... Coordenadas s:
      data psi /.333333333333333d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     .                      .5d0,.5d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     .    .333333333333333d0,.2d0,.6d0,.2d0,.0,.0,.0,.0,.0,.0,.0,.0,
     ..333333333333333d0,.470142064105115d0,.059715871789770d0,
     ..470142064105115d0,.101286507323456d0,.797426985353087d0,
     ..101286507323456d0,.0,.0,.0,.0,.0,
     ..063089014491502d0,.873821971016996d0,.063089014491502d0,
     ..249286745170910d0,.501426509658179d0,.249286745170910d0,
     ..310352451033785d0,.053145049844816d0,.636502499121399d0,
     ..053145049844816d0,.636502499121399d0,.310352451033785d0/
c ... Pesos para triangulos:
      data wf /                 1.d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     . .333333333333333d0,.333333333333333d0,.333333333333333d0,.0,.0,
     ..0,.0,.0,.0,.0,.0,.0,
     .-.562500000000000d0,.520833333333333d0,.520833333333333d0,
     . .520833333333333d0,.0,.0,.0,.0,.0,.0,.0,.0,
     . .225030000300000d0,.132394152788506d0,.132394152788506d0,
     . .132394152788506d0,.125939180544827d0,.125939180544827d0,
     . .125939180544827d0,.0,.0,.0,.0,.0,
     . .050844906370207d0,.050844906370207d0,.050844906370207d0,
     . .116786275726379d0,.116786275726379d0,.116786275726379d0,
     . .082851075618374d0,.082851075618374d0,.082851075618374d0,
     . .082851075618374d0,.082851075618374d0,.082851075618374d0/
c=======================================================================     
c
c ... Pontos de integracao em tetraedros:
c
c ... Numero de pontos de integracao (grau = 1, 2):
      data npint4 /1,4/    
c ... Coordenadas r:
      data pri4 / 0.25d0,      0.d0,        0.d0,        0.d0,
     .            0.58541020d0,0.13819660d0,0.13819660d0,0.13819660d0/
c ... Coordenadas s:
      data psi4 / 0.25d0,      0.d0,        0.d0,        0.d0,
     .            0.13819660d0,0.58541020d0,0.13819660d0,0.13819660d0/
c ... Coordenadas t:
      data pti4 / 0.25d0,      0.d0,        0.d0,        0.d0,
     .            0.13819660d0,0.13819660d0,0.58541020d0,0.13819660d0/
c ... Pesos para tetraedros:
      data wf4 /1.d0, 1.d0, 1.d0, 1.d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0/
c ......................................................................      
      end
