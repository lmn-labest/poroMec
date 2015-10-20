      subroutine elmt10(e,x,u,p,s,dt,ndm,nst,nel,isw)
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
c *           e(3) = coeficiente de Darcy                              *
c *           e(4) = coeficiente de M                                  *
c *           e(5) = coeficiente de Biot                               *
c *           e(6) = massa especifica do esquelto solido               *
c *           e(7) = massa especifica do fluido                        *
c *           e(8) = acelaracao da gravidade gx                        *
c *           e(8) = acelaracao da gravidade gy                        *
c *           e(10)= acelaracao da gravidade gz                        *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)       graus de liberade por elemento (u + p)            *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento ( 3*60 + 1*8)  *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 =                                                      *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 =                                                      *
c *           4 =                                                      *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - forcas internas (isw = 2) ou tensoes (isw = 6)             *
c *                                                                    *
c *   OBS:                                                             *
c *   ------------------                                               *
c *   u(1:60) - deslocamentos                                          *
c *   u(61:8) - diferencao de presso                                   *
c **********************************************************************
      implicit none
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k
      integer nen,nint,lx,ly,lz
      real*8 u(*)
      real*8 hu(20),hux(20),huy(20),huz(20)
      real*8 hp(8),hpx(8),hpy(8),hpz(8)
      real*8 xj(3,3),xji(3,3),g(3)
      real*8 ri,si,ti,w,wt1,wt2,det
      real*8 rn(20),sn(20),tn(20)
      real*8 pi,epsi(6),txi(6)
      real*8 dt
      real*8 pg(10,10),wg(10,10)
      real*8 e(*),x(ndm,*),p(*),s(nst,nst)
      real*8 darcy,a,b,c,m,ym,ps
      real*8 a1,a2
      real*8 dFluid,dtDarcy,dSolid
      real*8 dtFluid,dtFluidDarcy,lambda,mi,biot
c 
      data rn / 1.d0,-1.d0,-1.d0,1.d0, ! r1, r2, r3, r4
     .          1.d0,-1.d0,-1.d0,1.d0, ! r5, r6, r7, r8              
     .          0.d0,-1.d0, 0.d0,1.d0, ! r9,r10,r11,r12              
     .          0.d0,-1.d0, 0.d0,1.d0, !r13,r14,r15,r16               
     .          1.d0,-1.d0,-1.d0,1.d0/ !r17,r18,r19,r20  
c
      data sn / 1.d0, 1.d0,-1.d0,-1.d0, ! s1, s2, s3, s4           
     .          1.d0, 1.d0,-1.d0,-1.d0, ! s5, s6, s7, s8             
     .          1.d0, 0.d0,-1.d0, 0.d0, ! s9,s10,s11,s12             
     .          1.d0, 0.d0,-1.d0, 0.d0, !s13,s14,s15,s16             
     .          1.d0, 1.d0,-1.d0,-1.d0/ !s17,s18,s19,s20  
c
      data tn / 1.d0, 1.d0,-1.d0,-1.d0, ! t1, t2, t3, t4
     .          1.d0, 1.d0,-1.d0,-1.d0, ! t5, t6, t7, t8
     .          1.d0, 0.d0,-1.d0, 0.d0, ! t9,t10,t11,t12
     .          1.d0, 0.d0,-1.d0, 0.d0, !t13,t14,t15,t16      
     .          1.d0, 1.d0,-1.d0,-1.d0/ !t17,t18,t19,t20       
c
      data nen/20/
c ......................................................................
      goto (100,200,300,400) isw
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
      ym      = e(1)
      ps      = e(2)
      darcy   = e(3)
      m       = 1.d0/e(4)
      biot    = e(5)
c ... 
      dtDarcy = darcy*dt
c ...
      a1      = 1.d0 - ps
      a2      = 1.d0 - 2.d0*ps
c ...
      a       = ym*a1/((1.d0+ps)*a2)
      b       = ps/a1
      c       = 0.5d0*a2/a1 
c .....................................................................
c
c ... Matriz de rigidez:
      do 205 i = 1, nst
         do 205 j = 1, nst
            s(i,j) = 0.d0
  205 continue
c .....................................................................
c
c ... 
      nint = 3
      do 211 lz = 1, nint
        ti = pg(lz,nint)
        do 212 ly = 1, nint
          si = pg(ly,nint)
          do 213 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8 (hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
            call jacob3d (hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
            wt1 = w*a
c .....................................................................
c
c ...
            call sfhexa20 (hu,hux,huy,huz,ri,si,ti,.false.,.true.)
            call jacob3d (hux,huy,huz,xj,xji,x,det,20,nel,.false.)
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
            do 214 i = 1, 20
c             l1 = (i-1)*3+1
              l1 = 3*i-2
              l2 = l1 + 1
              l3 = l2 + 1
              do 215 j = 1, 20
c               k1 = (j-1)*3+1
                k1 = 3*j-2
                k2 = k1 + 1
                k3 = k2 + 1
                s(l1,k1) = s(l1,k1) +
     .          (hux(i)*hux(j) + c*huy(i)*huy(j) + c*huz(i)*huz(j))*wt1
                s(l1,k2) = s(l1,k2) +
     .          (b*hux(i)*huy(j) + c*huy(i)*hux(j))*wt1
                s(l1,k3) = s(l1,k3) + 
     .          (b*hux(i)*huz(j) + c*huz(i)*hux(j))*wt1
                s(l2,k1) = s(l2,k1) + 
     .          (b*huy(i)*hux(j) + c*hux(i)*huy(j))*wt1
                s(l2,k2) = s(l2,k2) + 
     .          (huy(i)*huy(j) + c*hux(i)*hux(j) + c*huz(i)*huz(j))*wt1
                s(l2,k3) = s(l2,k3) +
     .          (b*huy(i)*huz(j) + c*huz(i)*huy(j))*wt1
                s(l3,k1) = s(l3,k1) +
     .          (b*huz(i)*hux(j) + c*hux(i)*huz(j))*wt1
                s(l3,k2) = s(l3,k2) + 
     .          (b*huz(i)*huy(j) + c*huy(i)*huz(j))*wt1
                s(l3,k3) = s(l3,k3) +
     .          (huz(i)*huz(j) + c*huy(i)*huy(j) + c*hux(i)*hux(j))*wt1
  215         continue
  214       continue
c .....................................................................
c
c ... Kpu ( Int((b*(N~t)*(B)*dV) )  
            do 217 i = 1, 8
              l   = i + 60
              wt1 = b*hp(i)*w
              do 218 j = 1, 20
c               k1 = (j-1)*3+1
                k1 = 3*j-2
                k2 = k1 + 1
                k3 = k2 + 1
                s(l,k1) = s(l,k1) + hux(j)*wt1
                s(l,k2) = s(l,k2) + huy(j)*wt1
                s(l,k3) = s(l,k3) + huz(j)*wt1
  218         continue
  217       continue
c ................................................................ 
  213     continue
  212   continue
  211 continue
c .....................................................................
c
c ... Kpp ( Int((1/M)*(N~t)*(N~)*dV) )    
      nint = 2
      do 219 lz = 1, nint
        ti = pg(lz,nint)
        do 220 ly = 1, nint
          si = pg(ly,nint)
          do 221 lx = 1, nint
            ri = pg(lx,nint)
c ...
            call sfhexa8 (hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
            call jacob3d (hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c ..............................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
            wt1 = w*m
            wt2 = w*dtDarcy
c ..............................................................
c
c ...    
            do 222 i = 1, 8
              l = i + 60
              do 223 j = 1, 8
                k = j + 60
                s(l,k) = s(l,k) + wt1*(hp(i)*hp(j))
     .          + wt2*(hpx(i)*hpx(j) + hpy(i)*hpy(j) + hpz(i)*hpz(j))
  223         continue
  222       continue
c .....................................................................           
  221     continue
  220   continue
  219 continue 
c .....................................................................
c
c ... Forcas internas:
      call lku(s,u,p,nst)
c .....................................................................
c
c ... 
      call printmatrixs(s,nst)
c .....................................................................      
      return  
c ======================================================================
c
c ... Forcas internas:
c
c ......................................................................
 300  continue
c
c ... Matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
      darcy     = e(3)
      dSolid    = e(6)
      dFluid    = e(7)
      g(1)      = e(8)
      g(2)      = e(9)
      g(3)      = e(10)
c .....................................................................
c
c ...
      dtDarcy      = darcy*dt
      dtFluid      = dFluid*dt
      dtFluidDarcy = darcy*dFluid*dt
c .....................................................................
c
c ...
      a1      = 1.d0 - ps
      a2      = 1.d0 - 2.d0*ps
c ...
      a       = ym*a1/((1.d0+ps)*a2)
      b       = ps/a1
      c       = 0.5d0*a2/a1 
c .....................................................................
c
c ...      
      p(1:nst) = 0.d0
c .....................................................................
c
c ... Fu = int(BeT*sigma*dv)
      nint = 2
      do 301 lz = 1, nint
        ti = pg(lz,nint)
        do 302 ly = 1, nint
          si = pg(ly,nint)
          do 303 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8 (hp,hpx,hpy,hpz,ri,si,ti,.false.,.true.)
            call jacob3d (hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
            wt1 = w  
c .....................................................................
c
c ...
            call sfhexa20 (hu,hux,huy,huz,ri,si,ti,.false.,.true.)
            call jacob3d (hux,huy,huz,xj,xji,x,det,20,nel,.false.)
c .....................................................................
c
c ... deformacoes nos pontos de integracao
            call deform3d(hux,huy,huz,u,epsi,20)
            call stress3d(a,b,c,epsi,txi)
c .....................................................................            
c ...
            do 391 i = 1, 20
c             l1  = (i-1)*3+1
              l1  = 3*i-2
              l2  = l1 + 1
              l3  = l2 + 1
              p(l1) = p(l1) 
     .        - wt1*(hux(i)*txi(1) + huy(i)*txi(4) + huz(i)*txi(6)) 
              p(l2) = p(l2) 
     .        - wt1*(hux(i)*txi(4) + huy(i)*txi(2) + huz(i)*txi(5)) 
              p(l1) = p(l1) 
     .        - wt1*(hux(i)*txi(6) + huy(i)*txi(5) + huz(i)*txi(3)) 
  391       continue
c .....................................................................
  303     continue
  302   continue
  301 continue
c .....................................................................
c
c ... Fu = int(dSold*Ne*ge*dv)                                       
      nint = 2
      do 304 lz = 1, nint
        ti = pg(lz,nint)
        do 305 ly = 1, nint
          si = pg(ly,nint)
          do 306 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8 (hp,hpx,hpy,hpz,ri,si,ti,.false.,.true.)
            call jacob3d (hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
            wt1 = w*dSolid 
c .....................................................................
c
c ...
            call sfhexa20 (hu,hux,huy,huz,ri,si,ti,.true.,.false.)
c .....................................................................
c
c ...
            do 307 i = 1, 20
              wt2 = hu(i)*wt1  
c             l1  = (i-1)*3+1
              l1  = 3*i-2
              l2  = l1 + 1
              l3  = l2 + 1
              p(l1) = p(l1) + wt2*g(1)
              p(l2) = p(l2) + wt2*g(2)
              p(l3) = p(l3) + wt2*g(3)
  307       continue
c .....................................................................
  306     continue
  305   continue
  304 continue
c ....................................................................
c
c ... Fp = int(dt*darcy*dFluid*(B~T)*ge*dv)  
      nint = 2
      do 308 lz = 1, nint
        ti = pg(lz,nint)
        do 309 ly = 1, nint
          si = pg(ly,nint)
          do 310 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8 (hp,hpx,hpy,hpz,ri,si,ti,.false.,.true.)
            call jacob3d (hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
            wt1 = w*dtFluidDarcy 
c .....................................................................
c
c ...
            print*,ri,si,ti,hpz(5),det
            do 311 i = 1, 8 
              l    = i + 60
              p(l) = p(l) 
     .             + wt1*(hpx(i)*g(1) + hpy(i)*g(2) + hpz(i)*g(3))
 
  311       continue
c .....................................................................
  310     continue
  309   continue
  308 continue
c ....................................................................
c
c ... Fp = int(dt*k*(B~T)*(B~)*pe*dv)                                       
      nint = 2
      do 312 lz = 1, nint
        ti = pg(lz,nint)
        do 313 ly = 1, nint
          si = pg(ly,nint)
          do 314 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8 (hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
            call jacob3d (hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
            wt1 = w*dtDarcy 
c .....................................................................
c
c ... interpola p no ponto de integracao
            pi = 0.0d0
            do i = 1, 8
              l  = i + 60  
              pi = pi + hp(i)*u(l)    
            enddo      
c .....................................................................
c
c .....................................................................            
            do 315 i = 1, 8 
              l    = i + 60
              p(l) = p(l) 
     .        - wt1*pi*(hpx(i)*hpx(i) + hpy(i)*hpy(i) + hpz(i)*hpz(i))
  315       continue
c .....................................................................
  314     continue
  313   continue
  312 continue
c ....................................................................
      return
c ======================================================================
c
c ... Forcas internas:
c
c ......................................................................
  400 continue
      return
c ======================================================================
c
c ... Matriz de massa:
c
c ......................................................................
  500 continue
      print*,'Matriz de massa nao disponivel para o elemento tipo 10 !'
      stop
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      print*,'Matriz geometrica nao disponivel p/ o elemento tipo 10 !'
      stop
c ....................................................................
      end
c *********************************************************************
      subroutine printmatrixs(s,nst)
      implicit none
      integer i,j,nst
      real*8 s(nst,*)
c ...  
      open(14,file='sEl.dat',action='write')
       do i = 1, nst
        write(14,*),i,i,s(i,i)          
      enddo
      do i = 1, nst
        do j = 1, i-1
          write(14,*),i,j,s(i,j)          
        enddo    
      enddo
      close(14)
c .....................................................................      
      return
      end
