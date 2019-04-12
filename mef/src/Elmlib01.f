c *********************************************************************
c * Biblioteca de elementos termicos                                  * 
c *********************************************************************
c * ------------------ Elementos lineares --------------------------- *
c *                                                                   *
c * ELMT07_PM - hexaedros de  8 nos para o problema termico           *
c *                                                                   *
c * ------------------ Elementos quadraticos ------------------------ *
c *                                                                   *
c * ----------------------------------------------------------------- *
c *                                                                   *
c *********************************************************************       
       subroutine elmt07_t(e,iq,xl,u,v,p,s,ndm,nst,nel,ma,isw)
c **********************************************************************
c * Data de criacao    : 11/04/2019                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *       
c * elmt07: Elemento hexaedrico (difusao)                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = coeficiente de difusao                            *
c *           e(2) = capacidade termica (rocp)                         *
c * iq (7) - cargas por elemento                                       *
c * xl(ndm,nem)- coordenadas nodais locais                             *
c * u(nst)     - valores prescritos                                    *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * p   - nodal forces                                                 *
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento                    *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = fluxos                                                        *
c *  4 = forcas de volume e superficies                                *
c *  5 =                                                               *
c *  6 =                                                               *
c *  7 =                                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  fluxo                                                 *
c *     isw = 4  cargas de superfice e volume                          * 
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
       
c **********************************************************************
      implicit none
      include 'transiente.fi'                  
      include 'gauss.fi'                 
      integer iq(*),ndm,nst,nel,isw,nint,lx,ly,lz,i,j,i1,i2,i3,ma
      real*8 e(*),xl(ndm,*),u(*),v(*),p(*),s(nst,*),m(8,8),q
      real*8 h(8),hx(8),hy(8),hz(8),xj(3,3),xji(3,3),det,w,wt
      real*8 rn(8),sn(8),tn(8),ri,si,ti
      real*8 tm,difus,capac
      data rn / 1.d0,-1.d0,-1.d0, 1.d0     ! r1, r2, r3, r4
     .        , 1.d0,-1.d0,-1.d0, 1.d0 /   ! r5, r6, r7, r8
      data sn / 1.d0, 1.d0,-1.d0,-1.d0     ! s1, s2, s3, s4
     .        , 1.d0, 1.d0,-1.d0,-1.d0 /   ! s5, s6, s7, s8
      data tn / 1.d0, 1.d0, 1.d0, 1.d0     ! t1, t2, t3, t4
     .        ,-1.d0,-1.d0,-1.d0,-1.d0 /   ! t5, t6, r7, t8
      data nint / 2 /
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
      if (isw .eq. 3) goto 300
c
c ... Matriz K:                     
c
      do i = 1, nst
        do j = 1, nst
          s(i,j) = 0.d0
          m(i,j) = 0.d0
  205   enddo
      enddo
c .....................................................................
c
c ...
      do 250 lz = 1, nint
         ti = pg(lz,nint)
         do 240 ly = 1, nint
            si = pg(ly,nint)
            do 230 lx = 1, nint
               ri = pg(lx,nint)
c ...
               call sfhexa8_m(h,hx,hy,hz,ri,si,ti,.true.,.true.)
               call jacob3d_m(hx,hy,hz,xj,xji,xl,det,8,nel,.true.)
c ......................................................................
c
c ...
               w  = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c ......................................................................
c
c ...            
               wt = w*difus
               do 220 j = 1, 8
                 do 210 i = 1, 8
c ......................................................................
                      s(i,j) = s(i,j)
     .                       +(hx(i)*hx(j)+hy(i)*hy(j)+hz(i)*hz(j))*wt
c ......................................................................
                       m(i,j) = m(i,j) + h(i)*h(j)*w*capac
c ......................................................................            
  210            continue
c ......................................................................
  220          continue
c ......................................................................   
  230       continue
c ......................................................................
  240    continue
c ......................................................................
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
      do j = 1, nst
        do i = 1, nst
          s(i,j) = m(i,j)+alfa*dt*s(i,j)
        enddo
      enddo
c ......................................................................  
      return
c ======================================================================
c
c ... Tensoes nodais e fluxo nodais:
c
c ......................................................................
  300 continue
      do 310 i = 1, 8
         i1 = (i-1)*3+1
         i2 = i1+1
         i3 = i2+1
         call sfhexa8(h,hx,hy,hz,rn(i),sn(i),tn(i),.true.,.true.)
         call jacob3d(hx,hy,hz,xj,xji,xl,det,8,nel)    
         p(i1) = -(hx(1)*u(1)+hx(2)*u(2)+hx(3)*u(3)+hx(4)*u(4)+
     .             hx(5)*u(5)+hx(6)*u(6)+hx(7)*u(7)+hx(8)*u(8))*difus
         p(i2) = -(hy(1)*u(1)+hy(2)*u(2)+hy(3)*u(3)+hy(4)*u(4)+
     .             hy(5)*u(5)+hy(6)*u(6)+hy(7)*u(7)+hy(8)*u(8))*difus
         p(i3) = -(hz(1)*u(1)+hz(2)*u(2)+hz(3)*u(3)+hz(4)*u(4)+
     .             hz(5)*u(5)+hz(6)*u(6)+hz(7)*u(7)+hz(8)*u(8))*difus
  310 continue
      return
c ======================================================================
c
c ... Cargas distribuidas no volume, no contorno e variaveis no passo
c     anterior:
c
c ......................................................................
  400 continue
      if(iq(1) .eq. 0) return
      return  
c ......................................................................                 
      stop
      end
c *********************************************************************