      subroutine elmt7_pm(e,iq,x,u,dp,p,s,dt,ndm,nst,nel,isw,block_pu)
c **********************************************************************
c *                                                                    *
c *   ELMT07_pm: Elemento hexaedrico de 20 nos                         *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = coeficiente de Darcy                              *
c *           e(4) = modulo de Biot                                    *
c *           e(5) = coeficiente de Biot                               *
c *           e(6) = massa especifica homogenizada do meio poroso      *
c *           e(7) = massa especifica do fluido                        *
c *     x(ndm,nem) - coordenadas nodais locais                         *
c *     u(nst)     - graus de liberade por elemento (u + p)            *
c *     dp(*)      - delta p ( p(n  ,0  ) - p(0) )                     *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst) - nao definido                                      *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento (3*60 + 1*8)   *
c *     nel - numero do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *           1 =                                                      *
c *           2 = matriz K e forcas internas Ku                        *
c *           3 =                                                      *
c *           4 = forcas de volume, superficies e integrais do passo   *
c *             de tempo anterior                                      *
c *           5 =                                                      *
c *           6 =                                                      *
c *           7 =                                                      *
c *    block_pu    - true - armazenamento em blocos Kuu,Kpp e kpu      *
c *                  false- aramzenamento em unico bloco               *      
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
      include 'poro_mec_prop.fi'
      include 'load.fi'
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k,tp
      integer nen,nint,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_u(8),face_f(3),n_flux
      integer carg
c ...
      real*8 u(*),dp(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 hu(20),hux(20),huy(20),huz(20)
      real*8 hp(8),hpx(8),hpy(8),hpz(8)
c ...
      real*8 xj(3,3),xji(3,3),xj2D(2,2),xji2D(2,2),r(3,3)
      real*8 ri,si,ti,w,wt1,wt2,wt3,det
      real*8 rn(20),sn(20),tn(20)
      real*8 pi,epsi(6),txi(6),dpm
      real*8 pg(10,10),wg(10,10)
c ...
      real*8 dt,dt_c
      real*8 e(*),x(ndm,*),xf(3,4)
c ...
      real*8 l_c
      real*8 perm,a,b,c,ym,ps
      real*8 fluid_d,dt_perm,pm_d
      real*8 dt_fluid,dt_fluid_perm,lambda,mi,gl(3)
      real*8 imod_biot,coef_biot
      real*8 fluid_sw
      real*8 a1,a2,a3,tmp
c ... 
      integer hexa_face_node20(8,6),no
c ...
      logical block_pu
      data hexa_face_node20 / 1, 2, 3, 4, 9,10,11,12
     .                      , 5, 6, 7, 8,13,14,15,16
     .                      , 5, 6, 2, 1,13,18, 9,17
     .                      , 8, 4, 3, 7,20,11,12,15
     .                      , 1, 4, 8, 5,12,20,16,17
     .                      , 2, 6, 7, 3,18,17,12,10/
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
      goto (100,200,300,400,500) isw
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
c ...
      gl(1)     =  gravity(1)
      gl(2)     =  gravity(2)
      gl(3)     =  gravity(3)
c ... 
      fluid_d   = e(7)

c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod

c ... Matriz constitutiva:
      ym       = e(1)
      ps       = e(2)
c ... 
      perm     = e(3)/fluid_sw
c ...
      imod_biot= 1.d0/e(4)
      coef_biot= e(5)
c ... 
      dt_perm  = perm*dt
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - 2.d0*ps
      a3       = 1.d0 + ps
c ...
      a        = (ym*a1)/(a3*a2)
      b        = ps/a1
      c        = 0.5d0*(a2/a1) 
c .....................................................................
c
c ...
c     l_c        = 0.1d0*0.1d0
c     dt_c       = ((l_c)/perm) 
c    .           * ( imod_biot+( coef_biot*coef_biot*a3*a2 )/( ym*a1 ) )
c     print*,dt_c
c     l_c        = 1.d0
c     dt_c       = ((l_c)/perm) 
c    .           * ( imod_biot+( coef_biot*coef_biot*a3*a2 )/( ym*a1 ) )
c     print*,dt_c
c     if (dt .gt. dt_c) then
c       print*,'numero de elemento',nel,'delta critico',dt_c,dt
c     endif
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... 
      nint = 3
      do lz = 1, nint
        ti = pg(lz,nint)
        do ly = 1, nint
          si = pg(ly,nint)
          do lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
            call jacob3d (hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....................................................................
c
c ...
            call sfhexa20_m(hu,hux,huy,huz,ri,si,ti,.false.,.true.)
            call jacob3d (hux,huy,huz,xj,xji,x,det,20,nel,.false.)
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
            wt1 = w*a
            do i = 1, 20
c             l1 = (i-1)*3+1
              l1 = 3*i-2
              l2 = l1 + 1
              l3 = l2 + 1
              do j = 1, 20
c               k1 = (j-1)*3+1
                k1 = 3*j-2
                k2 = k1 + 1
                k3 = k2 + 1
                s(l1,k1) = s(l1,k1) +
     .          (hux(i)*hux(j) + c*(huy(i)*huy(j) + huz(i)*huz(j)))*wt1
c
                s(l1,k2) = s(l1,k2) +
     .          (b*hux(i)*huy(j) + c*huy(i)*hux(j))*wt1
c
                s(l1,k3) = s(l1,k3) + 
     .          (b*hux(i)*huz(j) + c*huz(i)*hux(j))*wt1
c
                s(l2,k1) = s(l2,k1) + 
     .          (b*huy(i)*hux(j) + c*hux(i)*huy(j))*wt1
c
                s(l2,k2) = s(l2,k2) + 
     .          (huy(i)*huy(j) + c*(hux(i)*hux(j) + huz(i)*huz(j)))*wt1
c
                s(l2,k3) = s(l2,k3) +
     .          (b*huy(i)*huz(j) + c*huz(i)*huy(j))*wt1
c
                s(l3,k1) = s(l3,k1) +
     .          (b*huz(i)*hux(j) + c*hux(i)*huz(j))*wt1
c
                s(l3,k2) = s(l3,k2) + 
     .          (b*huz(i)*huy(j) + c*huy(i)*huz(j))*wt1
c
                s(l3,k3) = s(l3,k3) +
     .          (huz(i)*huz(j) + c*(huy(i)*huy(j) + hux(i)*hux(j)))*wt1
              enddo
            enddo   
c .....................................................................
c
c ...-Kpu ( Int((b*(N~t)*(1t)(B)*dV) )  
            wt2 = w*coef_biot
            do i = 1, 8
              l   = i + 60
              wt1 = hp(i)*wt2
              do j = 1, 20
c               k1      = (j-1)*3+1
                k1      = 3*j-2
                k2      = k1 + 1
                k3      = k2 + 1
                s(l,k1) = s(l,k1) - hux(j)*wt1
                s(l,k2) = s(l,k2) - huy(j)*wt1
                s(l,k3) = s(l,k3) - huz(j)*wt1
              enddo
            enddo
c ................................................................ 
c
c ...    
            wt1 = w*imod_biot
            wt2 = w*dt_perm
c ...................................................................
c
c ... Kpp ( Int((1/M)*(N~t)*(N~)*dV) )    
            do i = 1, 8
              l = i + 60
              do j = 1, 8
                k = j + 60
                s(l,k) = s(l,k) - wt1*hp(i)*hp(j)
     .          - wt2*(hpx(i)*hpx(j) + hpy(i)*hpy(j) + hpz(i)*hpz(j))
              enddo   
            enddo   
c .....................................................................
          enddo
        enddo
      enddo
c .....................................................................
c
c ...Kup = kpu  
      do i = 1, 8 
        l   = i + 60 
        do j = 1, 20 
c         k1 = (j-1)*3+1 
          k1      = 3*j-2 
          k2      = k1 + 1  
          k3      = k2 + 1
          s(k1,l) = s(l,k1)
          s(k2,l) = s(l,k2)
          s(k3,l) = s(l,k3)
        enddo
      enddo
c ................................................................ 
c
c ...
      if(block_pu) then
c ... -Kpp
        do i = 1, 8
          l = i + 60
          do j = 1, 8
            k = j + 60
            s(l,k) = -s(l,k)
          enddo  
        enddo   
c ....................................................................
c
c ... Kpu = -Kpu
        do i = 1, 8
          l   = i + 60
          do j = 1, 20
c           k1 = (j-1)*3+1
            k1      = 3*j-2
            k2      = k1 + 1
            k3      = k2 + 1
            s(l,k1) = -s(l,k1)
            s(l,k2) = -s(l,k2)
            s(l,k3) = -s(l,k3)
          enddo
        enddo
c .....................................................................
      endif
c .....................................................................
c
c ... Forcas internas:
      call lku(s,u,p,nst)
c .....................................................................
c
c .....................................................................      
      return  
c ======================================================================
c
c ... Tensoes nodais e fluxo:
c
c ......................................................................
  300 continue
      return
c ======================================================================
c
c ... Cargas distribuidas no volume, no contorno e variaveis no passo
c     anterior:
c
c ......................................................................
 400  continue
c
c ... Matriz constitutiva:
      gl(1)     =  gravity(1)
      gl(2)     =  gravity(2)
      gl(3)     =  gravity(3)
c ...
      pm_d      = e(6)
      fluid_d   = e(7)

c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod
c ...      
      ym        = e(1)
      ps        = e(2)
c ...
      coef_biot = e(5)
      perm      = e(3)/fluid_sw
c .....................................................................
c
c ...
      dt_perm       = perm*dt
      dt_fluid      = fluid_d*dt
      dt_fluid_perm = fluid_d*dt_perm
c .....................................................................
c
c ...
      a1      = 1.d0 - ps
      a2      = 1.d0 - 2.d0*ps
      a3      = 1.d0 + ps
c ...
      a       = (ym*a1)/(a3*a2)
      b       = ps/a1
      c       = 0.5d0*(a2/a1) 
c .....................................................................
c
c ...                            
      nint = 3
      do lz = 1, nint
        ti = pg(lz,nint)
        do ly = 1, nint
          si = pg(ly,nint)
          do lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
            call jacob3d (hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....................................................................
c
c ...
            call sfhexa20_m(hu,hux,huy,huz,ri,si,ti,.true.,.true.)
            call jacob3d (hux,huy,huz,xj,xji,x,det,20,nel,.false.)
c .....................................................................
c
c ... deformacoes nos pontos de integracao
            call deform3d(hux,huy,huz,u,epsi,20)
c ... tensao nos pontos de integracao
            call stress3d(a,b,c,epsi,txi)
c .....................................................................            
c
c ... delta p nos pontos de integracao
            dpm = 0.0d0
            do i = 1, 8 
              dpm = dpm + hp(i)*dp(i)
            enddo    
c .....................................................................
c
c ...
            dpm    = coef_biot*dpm
            txi(1) = txi(1) - dpm
            txi(2) = txi(2) - dpm
            txi(3) = txi(3) - dpm
c .....................................................................
c
c ...
            wt1 = w  
            wt2 = w*pm_d 
c .....................................................................
c
c ...
            do i = 1, 20
c             l1  = (i-1)*3+1
              l1  = 3*i-2
              l2  = l1 + 1
              l3  = l2 + 1
c ... Fu = int(BeT*sigma*dv)
              p(l1) = p(l1) 
     .        + wt1*(hux(i)*txi(1) + huy(i)*txi(4) + huz(i)*txi(6)) 
              p(l2) = p(l2) 
     .        + wt1*(hux(i)*txi(4) + huy(i)*txi(2) + huz(i)*txi(5)) 
              p(l3) = p(l3) 
     .        + wt1*(hux(i)*txi(6) + huy(i)*txi(5) + huz(i)*txi(3))
c .....................................................................
c
c ... Fp = int(dt*perm*fuild_d*(B~T)*ge*dv)  
              wt3 = hu(i)*wt2  
              p(l1) = p(l1) - wt3*gl(1)
              p(l2) = p(l2) - wt3*gl(2)
              p(l3) = p(l3) - wt3*gl(3)
            enddo   
c .....................................................................
c
c ...
            wt1 = w*dt_fluid_perm 
            wt2 = w*dt_perm 
c .....................................................................
c
c ...
            do i = 1, 8 
              l    = i + 60
c ... Fp = int(dt*perm*fuild_d*(B~T)*ge*dv)  
              p(l) = p(l) 
     .             + wt1*(hpx(i)*gl(1) + hpy(i)*gl(2) + hpz(i)*gl(3))
c .....................................................................
c
c ... Fp = int(dt*k*(B~T)*(B~)*pe*dv)
              tmp = 0.0d0
              do j = 1, 8
                k    = j + 60
                tmp  = tmp       
     .          + wt2*u(k)*(hpx(i)*hpx(j)+hpy(i)*hpy(j)+hpz(i)*hpz(j))
              enddo 
              p(l) = p(l) - tmp
c .....................................................................
            enddo 
c .....................................................................
          enddo   
        enddo   
      enddo   
c .....................................................................
c
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 1 2 3 4  9 10 11 12 |
c             2 | no 5 6 7 8 13 14 15 16 |
c             3 | no 5 6 2 1 13 18  9 17 |     
c             4 | no 8 4 3 7 20 11 12 15 |
c             5 | no 1 4 8 5 12 20 16 17 |
c             6 | no 2 6 7 3 18 17 12 10 |
      tp = 0
      do i = 1, 6
        tp    = tp + iq(i)
      enddo
      if( tp .gt. 0 ) then
c ... face 1
        if(iq(1) .gt. 0 ) then
c ...
          do i = 1, 4
            no       = hexa_face_node20(i,1)
            xf(1,no) = x(1,no) 
            xf(2,no) = x(2,no) 
            xf(3,no) = x(3,no)
          enddo
c ...     
          call rotmatrix(xf,r)
          call rotate(xf(1,1),r,xf(1,1),.false.)
          call rotate(xf(1,2),r,xf(1,2),.false.)
          call rotate(xf(1,3),r,xf(1,3),.false.)
          call rotate(xf(1,4),r,xf(1,4),.false.)
c ...
          carg = load(1,iq(1))
c ... forcas
          if( carg .eq. 40 ) then
            call tload(iq(1),0.d0,face_u,n_flux,face_f)
            nint = 2
            do ly = 1, nint
              si = pg(ly,nint)
              do lx = 1, nint
                ri = pg(lx,nint)
c ...                                               
                call sfquad4_m(hp,hpx,hpy,ri,si,.false.,.true.)
                call jacob2d(hpx,hpy,xj2D,xji2D,xf,det,4,ndm
     .                      ,nel,.true.)
c .....................................................................
c
c ...                                               
                w   = wg(lx,nint)*wg(ly,nint)*det
c .....................................................................
c
c ...                                               
                call sfquad8_m(hu,hux,huy,ri,si,.true.,.false.)
c .....................................................................
c
c ...
                do i = 1, 8
                  no  = hexa_face_node20(i,1)
c                 l1  = (no-1)*3+1
                  wt1 = hu(i)*w
                  l1  = 3*no-2
                  l2  = l1 + 1
                  l3  = l2 + 1
                  p(l1) = p(l1) - face_f(1)*wt1
                  p(l2) = p(l2) - face_f(2)*wt1
                  p(l3) = p(l3) - face_f(3)*wt1  
                enddo
c .....................................................................
              enddo
            enddo
c .....................................................................
c
c .....................................................................
c
c ... fluxo 
          elseif( carg .eq. 41) then
          endif
c .....................................................................
c
c ... face 2
        else if(iq(2) .gt. 0 ) then
          carg = load(1,iq(2))
c .....................................................................
        endif
      endif
c .....................................................................
c
c ... bloco Fp = -Fp
      if(block_pu) then
        p(61:nst) =-p(61:nst)
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Forcas internas:
c
c ......................................................................
  500 continue
      stop
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      stop
c ....................................................................
      end
c *********************************************************************
c
c *********************************************************************
c     subroutine printmatrixs(s,nst)
c     implicit none
c     integer i,j,nst
c     real*8 s(nst,*)
c ...  
c     open(14,file='sEl.dat',action='write')
c     do i = 1, nst
c       write(14,*),i,i,s(i,i)          
c     enddo
c     do i = 1, nst
c       do j = 1, i-1
c         write(14,*),i,j,s(i,j)          
c       enddo    
c     enddo
c     do i = 1, 61
c       do j = 61, 68
c         write(14,*),i,j,s(i,j)          
c       enddo    
c     enddo
c     do i = 61, 68
c       do j = 1, 60
c         write(14,*),i,j,s(i,j)          
c       enddo    
c     enddo
c     do i = 1, nst
c       write(14,'(100e10.3)')(s(i,j),j=1,nst)          
c     enddo
c     do i = 61, nst
c       do j = 61, nst
c         write(14,*),i-60,j-60,s(i,j)          
c       enddo    
c     enddo
c     close(14)
c .....................................................................      
c     return
c     end
c *********************************************************************
