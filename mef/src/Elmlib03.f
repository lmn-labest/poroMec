      subroutine elmt6_pm(e,iq,x,u,dp,p,s,txn,dt,ndm,nst,nel,isw
     .                   ,block_pu)
c **********************************************************************
c * Data de criacao    : 27/03/2016                                    *
c * Data de modificaco : 28/03/2016                                    * 
c * ------------------------------------------------------------------ *      
c * ELMT06_PM: Elemento tetraedrico de 10 nos para problemas           *  
c * poromecanico elasticos                                             *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = coeficiente de Darcy                              *
c *           e(4) = modulo de Biot                                    *
c *           e(5) = coeficiente de Biot                               *
c *           e(6) = massa especifica homogenizada do meio poroso      *
c *           e(7) = massa especifica do fluido                        *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u + p)                *
c * dp(*)      - delta p ( p(n  ,0  ) - p(0) )                         *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(6,nen) - tensoes nodais                                        *
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (3*10 + 1*4)       *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 = delta t critico                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tesnsoes e fluxos                                             *
c *  4 = forcas de volume, superficies e integrais do passo            *
c *    de tempo anterior                                               *
c *  5 = Tensoes iniciais                                              *
c *  6 =                                                               *
c *  7 =                                                               *
c * block_pu - true - armazenamento em blocos Kuu,Kpp e kpu            *
c *            false- aramzenamento em unico bloco                     *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao e fluxo                                        *
c *     isw = 4  cargas de superfice, volume e integras do passo       *
c *         de tempo anterior                                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * u(1:30) - isw = 2 diferenca de deslocamentos                       *
c * u(1:30) - isw = 3 e isw = 4 deslocamentos                          *
c * u(31:34)- isw = 2 diferenca de pressao                             *
c * u(31:34)- isw = 3 e isw = 4 pressao                                *
c **********************************************************************
      implicit none
      include 'poro_mec_prop.fi'
      include 'load.fi'
      common /pint4/ pri4,psi4,pti4,wf4,npint4
      common /pint / pri,psi,wf,npint
      real*8  div6         
      parameter ( div6 = 0.166666666666667d0)
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k,tp,tp1
      integer nen,nint,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_u(3),face_f(3),n_flux
      integer carg
c ...
      real*8 u(*),dp(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 hu(10),hux(10),huy(10),huz(10)
      real*8 hp(4),hpx(4),hpy(4),hpz(4)
c ...
      real*8 xj(3,3),xji(3,3),xj2D(2,2),xji2D(2,2),r(3,3)
      real*8 ri,si,ti,w,wt1,wt2,wt3,det
      real*8 rn(10),sn(10),tn(10)
      real*8 pi,epsi(6),txi(6),txn(6,*),dpm,pm
c ... integracao numerica de tetraedros     
      real*8 pri4(5,3),psi4(5,3),pti4(5,3),wf4(5,3)
      real*8  pri(12,5),psi(12,5),wf(12,5)
      integer igrau,igrau_vol,igrau_face,npint(5),npint4(3)
c ...
      real*8 dt,dt_c
      real*8 e(*),x(ndm,*),xf(3,3)
c ...
      real*8 l_c
      real*8 perm,a,b,c,ym,ps
      real*8 fluid_d,dt_perm,pm_d
      real*8 dt_fluid,dt_fluid_perm,lambda,mi,gl(3)
      real*8 imod_biot,coef_biot
      real*8 fluid_sw
      real*8 a1,a2,a3,tmp
c ... 
      integer tetra_face_node10(6,4),no
      real*8 tetra_vol,volum
c ...
      logical block_pu
      data tetra_face_node10 / 2, 3, 4, 8, 9,10
     .                       , 1, 4, 3, 7, 9, 6
     .                       , 1, 2, 4, 5,10, 7
     .                       , 1, 3, 2, 6, 8, 5/
c
      data rn / 1.0d0, 0.0d0, 0.0d0 ! r1, r2, r3
     .       ,  0.0d0, 0.5d0, 0.5d0 ! r4, r5, r6              
     .       ,  0.5d0, 0.0d0, 0.0d0 ! r7, r8, r9              
     .       ,  0.0d0/              ! r10              
c
      data sn / 0.0d0, 1.0d0, 0.0d0 ! s1, s2, s3
     .        , 0.0d0, 0.5d0, 0.0d0 ! s4, s5, s6              
     .        , 0.0d0, 0.5d0, 0.0d0 ! s7, s8, s9              
     .        , 0.5d0/              ! s10              
c
      data tn / 0.0d0, 0.0d0, 0.0d0 ! t1, t2, t3
     .        , 1.0d0, 0.0d0, 0.0d0 ! t4, t5, t6              
     .        , 0.5d0, 0.0d0, 0.0d0 ! t7, t8, t9              
     .        , 0.5d0/              ! t10              
c
      data nen/10/,igrau_vol/2/,igrau_face/2/    
c ......................................................................
      goto (100,200,300,400,500) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
c ...
      fluid_d   = e(7)*1.0d-06

c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod

c ... matriz constitutiva:
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
c .....................................................................
c
c ...
      volum = tetra_vol(x)
      l_c   = volum**(1.0d0/3.d0)
c ...
      dt_c  = ((l_c)/perm) 
     .      * ( imod_biot+( coef_biot*coef_biot*a3*a2 )/( ym*a1 ) )
      p(1)  = dt_c
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ...
      gl(1)     =  gravity(1)
      gl(2)     =  gravity(2)
      gl(3)     =  gravity(3)
c ... 
      fluid_d   = e(7)*1.0d-06
c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod

c ... matriz constitutiva:
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
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... 
      igrau = igrau_vol 
      nint = npint4(igrau) 
      do 205 lx = 1, nint
c ...
        ti = pti4(lx,igrau)
        si = psi4(lx,igrau)
        ri = pri4(lx,igrau)
c ...                                               
        call sftetra4(hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,4,nel,.true.)
c .....................................................................
c
c ...
        w   = wf4(lx,igrau)*det
        w   = w*div6
c .....................................................................
c
c ...
        call sftetra10(hu,hux,huy,huz,ri,si,ti,.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,10,nel,.false.)
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
        wt1 = w*a
        do 210 i = 1, 10
c         l1 = (i-1)*3+1
          l1 = 3*i-2
          l2 = l1 + 1
          l3 = l2 + 1
          do 215 j = 1, 10
c           k1 = (j-1)*3+1
            k1 = 3*j-2
            k2 = k1 + 1
            k3 = k2 + 1
            s(l1,k1) = s(l1,k1) +
     .      (hux(i)*hux(j) + c*(huy(i)*huy(j) + huz(i)*huz(j)))*wt1
c
            s(l1,k2) = s(l1,k2) +
     .      (b*hux(i)*huy(j) + c*huy(i)*hux(j))*wt1
c
            s(l1,k3) = s(l1,k3) + 
     .      (b*hux(i)*huz(j) + c*huz(i)*hux(j))*wt1
c
            s(l2,k1) = s(l2,k1) + 
     .      (b*huy(i)*hux(j) + c*hux(i)*huy(j))*wt1
c
            s(l2,k2) = s(l2,k2) + 
     .      (huy(i)*huy(j) + c*(hux(i)*hux(j) + huz(i)*huz(j)))*wt1
c
            s(l2,k3) = s(l2,k3) +
     .      (b*huy(i)*huz(j) + c*huz(i)*huy(j))*wt1
c
            s(l3,k1) = s(l3,k1) +
     .      (b*huz(i)*hux(j) + c*hux(i)*huz(j))*wt1
c
            s(l3,k2) = s(l3,k2) + 
     .      (b*huz(i)*huy(j) + c*huy(i)*huz(j))*wt1
c
            s(l3,k3) = s(l3,k3) +
     .      (huz(i)*huz(j) + c*(huy(i)*huy(j) + hux(i)*hux(j)))*wt1
c .....................................................................            
  215     continue
c .....................................................................     
  210   continue 
c .....................................................................
c
c ...-Kpu ( Int((b*(N~t)*(1t)(B)*dV) )  
        wt2 = w*coef_biot
        do 220 i = 1, 4
          l   = i + 30
          wt1 = hp(i)*wt2
          do 225 j = 1, 10
c           k1      = (j-1)*3+1
            k1      = 3*j-2
            k2      = k1 + 1
            k3      = k2 + 1
            s(l,k1) = s(l,k1) - hux(j)*wt1
            s(l,k2) = s(l,k2) - huy(j)*wt1
            s(l,k3) = s(l,k3) - huz(j)*wt1
c .....................................................................
  225     continue
c .....................................................................
  220   continue    
c .....................................................................
c
c ...    
        wt1 = w*imod_biot
        wt2 = w*dt_perm
c ...................................................................
c
c ... Kpp ( Int((1/M)*(N~t)*(N~)*dV) )    
        do 230 i = 1, 4
          l = i + 30
          do 235 j = 1, 4
            k = j + 30
            s(l,k) = s(l,k) - wt1*hp(i)*hp(j)
     .      - wt2*(hpx(i)*hpx(j) + hpy(i)*hpy(j) + hpz(i)*hpz(j))
c .....................................................................
 235      continue
c .....................................................................
 230    continue    
c ....................................................................
 205  continue 
c .....................................................................
c
c ...Kup = kpu  
      do 240 i = 1, 4 
        l   = i + 30 
        do 245 j = 1, 10 
c         k1 = (j-1)*3+1 
          k1      = 3*j-2 
          k2      = k1 + 1  
          k3      = k2 + 1
          s(k1,l) = s(l,k1)
          s(k2,l) = s(l,k2)
          s(k3,l) = s(l,k3)
c ................................................................          
  245   continue 
c ................................................................
  240 continue
c ................................................................ 
c
c ...
      if(block_pu) then
c ... -Kpp
        do 250 i = 1, 4
          l = i + 30
          do 255 j = 1, 4
            k = j + 30
            s(l,k) = -s(l,k)
c ................................................................            
  255     continue 
c ....................................................................
  250   continue
c ....................................................................
c
c ... Kpu = -Kpu
        do 260 i = 1, 4
          l   = i + 30
          do 265 j = 1, 10
c           k1 = (j-1)*3+1
            k1      = 3*j-2
            k2      = k1 + 1
            k3      = k2 + 1
            s(l,k1) = -s(l,k1)
            s(l,k2) = -s(l,k2)
            s(l,k3) = -s(l,k3)
c .....................................................................       
  265     continue 
c .....................................................................
  260   continue
c .....................................................................
      endif
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c     if(block_pu) then
c       call lku_m(s,u,p,nst)
c     else
c       call lku_sym(s,u,p,nst)    
c     endif   
c ......................................................................
c
c ......................................................................   
      return  
c ======================================================================
c
c ... Tensoes nodais e fluxo nodais:
c
c ......................................................................
  300 continue

c ...                        
      gl(1)     =  gravity(1)
      gl(2)     =  gravity(2)
      gl(3)     =  gravity(3)
c ... 
      fluid_d   = e(7)*1.0d-06

c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod

c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ... 
      coef_biot = e(5)
      perm      = e(3)/fluid_sw
c ...
      a1        = 1.d0 - ps
      a2        = 1.d0 - 2.d0*ps
      a3        = 1.d0 + ps
c ...
      a         = (ym*a1)/(a3*a2)
      b         = ps/a1
      c         = 0.5d0*(a2/a1)
c .....................................................................
c
c ... tensao nodal total
      do 310 i = 1, 4
c       tp = (i-1)*6 + 1
        tp  = 6*i - 5
c ... tensao efetiva de biot (4x6+1=25)     
c       tp1 = (i-1)*6 + 25 
        tp1 = 6*i +19
c ...   p(1...4) = u(31...34)
        l   = i  + 30 
c ... calculo do terminante
        call sftetra4(hp,hpx,hpy,hpz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,4,nel ,.true.)
c ... calculo da derivadas das funcoes de interpolacao
        call sftetra10(hu,hux,huy,huz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,10,nel ,.false.)
c .....................................................................
        call deform3d(hux,huy,huz,u,epsi,10)
        call stress3d(a,b,c,epsi,p(tp))
c ... Tensao = D*deformacao elastica - biot*dp
        dpm     = coef_biot*dp(i)
         pm     = coef_biot*u(l)
c ... tensao total
        p(tp)   = p(tp)   - dpm
        p(tp+1) = p(tp+1) - dpm
        p(tp+2) = p(tp+2) - dpm
c ... tensao efetiva de biot
        p(tp1)  = p(tp)   + pm
        p(tp1+1)= p(tp+1) + pm
        p(tp1+2)= p(tp+2) + pm
        p(tp1+3)= p(tp+3) 
        p(tp1+4)= p(tp+4) 
        p(tp1+5)= p(tp+5) 
c .....................................................................
  310 continue
c .....................................................................
c
c ... fluxo nodal 
      do 320 i = 1, 4
c ... fuxo (4x6 + 4x6=49)            
c       tp = (i-1)*3 + 49
        tp = 3*i + 46
c ...
        call sftetra4(hp,hpx,hpy,hpz,rn(i),sn(i),tn(i),.true.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,4,nel,.true.)
c ... p nos pontos de integracao
        call darcy_flux(perm,gl,fluid_d,hpx,hpy,hpz,u(31),4,p(tp))  
c .....................................................................
  320 enddo
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume, no contorno e variaveis no passo
c     anterior:
c
c ......................................................................
 400  continue
c ...                        
      gl(1)     =  gravity(1)
      gl(2)     =  gravity(2)
      gl(3)     =  gravity(3)
c ...
      pm_d      = e(6)*1.0d-06           
      fluid_d   = e(7)*1.0d-06

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
      igrau = igrau_vol 
      nint  = npint4(igrau) 
      do 410 lx = 1, nint
c ...
        ti = pti4(lx,igrau)
        si = psi4(lx,igrau)
        ri = pri4(lx,igrau)
c ...                                               
        call sftetra4(hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,4,nel,.true.)
c .....................................................................
c
c ...
        w   = wf4(lx,igrau)*det
        w   = w*div6
c .....................................................................
c
c ...
       call sftetra10(hu,hux,huy,huz,ri,si,ti,.true.,.true.)
       call jacob3d_m(hux,huy,huz,xj,xji,x,det,10,nel,.false.)
c .....................................................................
c
c ... deformacoes nos pontos de integracao
       call deform3d(hux,huy,huz,u,epsi,10)
c ... tensao nos pontos de integracao
       call stress3d(a,b,c,epsi,txi)
c .....................................................................            
c
c ... delta p nos pontos de integracao
       dpm = hp(1)*dp(1) + hp(2)*dp(2) + hp(3)*dp(3) + hp(4)*dp(4) 
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
        do 420 i = 1, 10
c         l1  = (i-1)*3+1
          l1  = 3*i-2
          l2  = l1 + 1
          l3  = l2 + 1
c ... Fu = int(BeT*sigma*dv)
          p(l1) = p(l1) 
     .    + wt1*(hux(i)*txi(1) + huy(i)*txi(4) + huz(i)*txi(6))
c          
          p(l2) = p(l2) 
     .    + wt1*(hux(i)*txi(4) + huy(i)*txi(2) + huz(i)*txi(5)) 
c          
          p(l3) = p(l3) 
     .    + wt1*(hux(i)*txi(6) + huy(i)*txi(5) + huz(i)*txi(3))
c .....................................................................
c
c ... Fu = int(pm_d*(N~T)*ge*dv)  
          wt3   = hu(i)*wt2  
          p(l1) = p(l1) - wt3*gl(1)
          p(l2) = p(l2) - wt3*gl(2)
          p(l3) = p(l3) - wt3*gl(3)
  420   continue    
c .....................................................................
c
c ...
        wt1 = w*dt_fluid_perm 
        wt2 = w*dt_perm 
c .....................................................................
c
c ...
        do 425 i = 1, 4 
          l    = i + 30
c ... Fp = int(dt*perm*fuild_d*(B~T)*ge*dv)  
          p(l) = p(l) 
     .    + wt1*(hpx(i)*gl(1) + hpy(i)*gl(2) + hpz(i)*gl(3))
c .....................................................................
c
c ... Fp = int(dt*k*(B~T)*(B~)*pe*dv)
          tmp = 0.0d0
          do 430 j = 1, 4
            k    = j + 30
            tmp  = tmp       
     .      + wt2*u(k)*(hpx(i)*hpx(j)+hpy(i)*hpy(j)+hpz(i)*hpz(j))
 430      continue 
          p(l) = p(l) - tmp
c .....................................................................
 425    continue  
c .....................................................................
 410  continue   
c .....................................................................
c
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 2 3 4  8  9 10 |
c             2 | no 1 4 3  7  9  6 |
c             3 | no 1 2 4  5 10  7 |     
c             4 | no 1 3 2  6  8  5 |
c ... verifica se ha alguma face com carga      
      tp = 0
      do 435 i = 1, 4
        tp  = tp + iq(i)
  435 continue
c .....................................................................
c
c ...      
      if( tp .gt. 0 ) then
        do 440 j = 1, 4
c ... face  
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 3
              no       = tetra_face_node10(i,j)
              xf(1,i) = x(1,no) 
              xf(2,i) = x(2,no) 
              xf(3,i) = x(3,no)
  445      continue
c ...    
            call rotmatrix(xf,r)
            call rotate(xf(1,1),r,xf(1,1),.false.)
            call rotate(xf(1,2),r,xf(1,2),.false.)
            call rotate(xf(1,3),r,xf(1,3),.false.)
c ...
            carg = load(1,iq(j))
c ... forcas
            if( carg .eq. 40 ) then
              call tload(iq(j),0.d0,face_u,n_flux,face_f)
              igrau = igrau_face 
              nint = npint(igrau)
              do 450 lx = 1, nint
                ri = pri(lx,igrau)
                si = psi(lx,igrau)
c ...                                 
                call sftria3(hp,hpx,hpy,ri,si,.false.,.true.)
                call jacob2d_m(hpx,hpy,xj2D,xji2D,xf,det,3,ndm
     .                        ,nel,.true.)
c .....................................................................
c
c ...                                               
                w   = wf(lx,igrau)*det
                w   = 0.5d0*w
c .....................................................................
c
c ...                   
                call sftria6(hu,hux,huy,ri,si,.true.,.false.)
c .....................................................................
c
c ...
                do 455 i = 1, 6 
                  no  = tetra_face_node10(i,j)
c                 l1  = (no-1)*3+1
                  wt1   = hu(i)*w
                  l1    = 3*no-2
                  l2    = l1 + 1
                  l3    = l2 + 1
                  p(l1) = p(l1) - face_f(1)*wt1
                  p(l2) = p(l2) - face_f(2)*wt1
                  p(l3) = p(l3) - face_f(3)*wt1
  455           continue
c .....................................................................
  450         continue
c .....................................................................
c
c ... fluxo 
            elseif( carg .eq. 41) then
            endif
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
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
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
c ...     
      igrau = igrau_vol 
      nint  = npint4(igrau) 
      do 510 lx = 1, nint
c ...
        ti = pti4(lx,igrau)
        si = psi4(lx,igrau)
        ri = pri4(lx,igrau)
c ...                                               
        call sftetra4(hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,4,nel,.true.)
c .....................................................................
c
c ...
        w   = wf4(lx,igrau)*det
        wt1 = w*div6
c .....................................................................
c
c ...
        call sftetra10(hu,hux,huy,huz,ri,si,ti,.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,10,nel,.false.)
c .....................................................................
c  
c ...
        txi(1:6)= 0.d0
        do 515 j = 1,   4
          txi(1) = txi(1) + hp(j)*txn(1,j)
          txi(2) = txi(2) + hp(j)*txn(2,j) 
          txi(3) = txi(3) + hp(j)*txn(3,j) 
          txi(4) = txi(4) + hp(j)*txn(4,j) 
          txi(5) = txi(5) + hp(j)*txn(5,j)
          txi(6) = txi(6) + hp(j)*txn(6,j) 
  515   continue
c .....................................................................
c        
c ...
        do 520 i = 1, 10
c         l1  = (i-1)*3+1
          l1  = 3*i-2
          l2  = l1 + 1
          l3  = l2 + 1
c ... Fu = int(BeT*sigma*dv)
          p(l1) = p(l1) 
     .    + wt1*(hux(i)*txi(1) + huy(i)*txi(4) + huz(i)*txi(6))
c         
          p(l2) = p(l2) 
     .    + wt1*(hux(i)*txi(4) + huy(i)*txi(2) + huz(i)*txi(5))
c          
          p(l3) = p(l3) 
     .    + wt1*(hux(i)*txi(6) + huy(i)*txi(5) + huz(i)*txi(3))
c .....................................................................
  520   continue  
c .....................................................................
  510 continue 
c .....................................................................
      return
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
      subroutine elmt7_pm(e,iq,x,u,dp,p,s,txn,dt,ndm,nst,nel,isw
     .                   ,block_pu)
c **********************************************************************
c * Data de criacao    : 10/12/2015                                    *
c * Data de modificaco : 02/04/2016                                    * 
c * ------------------------------------------------------------------ *       
c * ELMT07_PM: Elemento hexaedricos de 20 nos para problemas           *  
c * poromecanico elasticos                                             *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = coeficiente de Darcy                              *
c *           e(4) = modulo de Biot                                    *
c *           e(5) = coeficiente de Biot                               *
c *           e(6) = massa especifica homogenizada do meio poroso      *
c *           e(7) = massa especifica do fluido                        *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u + p)                *
c * dp(*)      - delta p ( p(n  ,0  ) - p(0) )                         *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(6,nen) - tensoes nodais                                        *      
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (3*20 + 1*8)       *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 = delta t critico                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tesnsoes e fluxos                                             *
c *  4 = forcas de volume, superficies e integrais do passo            *
c *    de tempo anterior                                               *
c *  5 = Tensoes iniciais                                              *
c *  6 =                                                               *
c *  7 =                                                               *
c * block_pu - true - armazenamento em blocos Kuu,Kpp e kpu            *
c *            false- aramzenamento em unico bloco                     *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao e fluxo                                        *
c *     isw = 4  cargas de superfice, volume e integras do passo       *
c *     de tempo anterior                                              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * u(1:60) - isw = 2 diferenca de deslocamentos                       *
c * u(1:60) - isw = 3 e isw = 4 deslocamentos                          *
c * u(61:68)- isw = 2 diferenca de pressao                             *
c * u(61:68)- isw = 3 e isw = 4 pressao                                *
c **********************************************************************
      implicit none
      include 'poro_mec_prop.fi'
      include 'load.fi'
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k,tp,tp1
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
      real*8 pi,epsi(6),txi(6),txn(6,*),dpm,pm
c ... integracao numerica de tetraedros          
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
      real*8 hexa_vol,volum
c ...
      logical block_pu
      data hexa_face_node20 / 1, 2, 3, 4, 9,10,11,12
     .                      , 5, 6, 7, 8,13,14,15,16
     .                      , 1, 5, 6, 2,17,13,18, 9
     .                      , 4, 3, 7, 8,11,19,15,20
     .                      , 1, 4, 8, 5,12,20,16,17
     .                      , 2, 6, 7, 3,18,14,19,10/
c
      data rn / 1.d0,-1.d0,-1.d0, 1.d0, ! r1, r2, r3, r4
     .          1.d0,-1.d0,-1.d0, 1.d0, ! r5, r6, r7, r8              
     .          0.d0,-1.d0, 0.d0, 1.d0, ! r9,r10,r11,r12              
     .          0.d0,-1.d0, 0.d0, 1.d0, !r13,r14,r15,r16               
     .          1.d0,-1.d0,-1.d0, 1.d0/ !r17,r18,r19,r20  
c
      data sn / 1.d0, 1.d0,-1.d0,-1.d0, ! s1, s2, s3, s4           
     .          1.d0, 1.d0,-1.d0,-1.d0, ! s5, s6, s7, s8             
     .          1.d0, 0.d0,-1.d0, 0.d0, ! s9,s10,s11,s12             
     .          1.d0, 0.d0,-1.d0, 0.d0, !s13,s14,s15,s16             
     .          1.d0, 1.d0,-1.d0,-1.d0/ !s17,s18,s19,s20  
c
      data tn / 1.d0, 1.d0, 1.d0, 1.d0, ! t1, t2, t3, t4
     .         -1.d0,-1.d0,-1.d0,-1.d0, ! t5, t6, t7, t8
     .          1.d0, 1.d0, 1.d0, 1.d0, ! t9,t10,t11,t12
     .         -1.d0,-1.d0,-1.d0,-1.d0, !t13,t14,t15,t16      
     .          0.d0, 0.d0, 0.d0, 0.d0/ !t17,t18,t19,t20       
c
      data nen/20/
c ......................................................................
      goto (100,200,300,400,500) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
c ...
      fluid_d   = e(7)*1.0d-06

c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod

c ... matriz constitutiva:
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
c .....................................................................
c
c ...
      volum = hexa_vol(x)
      l_c   = volum**(1.0d0/3.d0)
c ...
      dt_c  = ((l_c)/perm) 
     .      * ( imod_biot+( coef_biot*coef_biot*a3*a2 )/( ym*a1 ) )
      p(1)  = dt_c
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ...
      gl(1)     =  gravity(1)
      gl(2)     =  gravity(2)
      gl(3)     =  gravity(3)
c ... 
      fluid_d   = e(7)*1.0d-06
c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod

c ... matriz constitutiva:
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
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... 
      nint = 4 
      do 205 lz = 1, nint
        ti = pg(lz,nint)
        do 210 ly = 1, nint
          si = pg(ly,nint)
          do 215 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
            call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....................................................................
c
c ...
            call sfhexa20_m(hu,hux,huy,huz,ri,si,ti,.false.,.true.)
            call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel,.false.)
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
            wt1 = w*a
            do 220 i = 1, 20
c             l1 = (i-1)*3+1
              l1 = 3*i-2
              l2 = l1 + 1
              l3 = l2 + 1
              do 225 j = 1, 20
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
c .....................................................................                   
  225         continue
c .....................................................................              
  220       continue    
c .....................................................................
c
c ...-Kpu ( Int((b*(N~t)*(1t)(B)*dV) )  
            wt2 = w*coef_biot
            do 230 i = 1, 8
              l   = i + 60
              wt1 = hp(i)*wt2
              do 235 j = 1, 20
c               k1      = (j-1)*3+1
                k1      = 3*j-2
                k2      = k1 + 1
                k3      = k2 + 1
                s(l,k1) = s(l,k1) - hux(j)*wt1
                s(l,k2) = s(l,k2) - huy(j)*wt1
                s(l,k3) = s(l,k3) - huz(j)*wt1
c .....................................................................
  235         continue
c .....................................................................
  230       continue    
c .....................................................................
c
c ...    
            wt1 = w*imod_biot
            wt2 = w*dt_perm
c ...................................................................
c
c ... Kpp ( Int((1/M)*(N~t)*(N~)*dV) )    
            do 240 i = 1, 8
              l = i + 60
              do 245 j = 1, 8
                k = j + 60
                s(l,k) = s(l,k) - wt1*hp(i)*hp(j)
     .          - wt2*(hpx(i)*hpx(j) + hpy(i)*hpy(j) + hpz(i)*hpz(j))
c .....................................................................
  245         continue
c .....................................................................
  240       continue    
c .....................................................................
  215     continue
c .....................................................................                
  210   continue 
c .....................................................................              
  205 continue
c .....................................................................
c
c ...Kup = kpu  
      do 250 i = 1, 8 
        l   = i + 60 
        do 255 j = 1, 20 
c         k1 = (j-1)*3+1 
          k1      = 3*j-2 
          k2      = k1 + 1  
          k3      = k2 + 1
          s(k1,l) = s(l,k1)
          s(k2,l) = s(l,k2)
          s(k3,l) = s(l,k3)
c ................................................................          
  255   continue 
c ................................................................
  250 continue
c ................................................................ 
c
c ...
      if(block_pu) then
c ... -Kpp
        do 260 i = 1, 8
          l = i + 60
          do 265 j = 1, 8
            k = j + 60
            s(l,k) = -s(l,k)
c ................................................................            
  265     continue 
c ....................................................................
  260   continue
c ....................................................................
c
c ... Kpu = -Kpu
        do 270 i = 1, 8
          l   = i + 60
          do 275 j = 1, 20
c           k1 = (j-1)*3+1
            k1      = 3*j-2
            k2      = k1 + 1
            k3      = k2 + 1
            s(l,k1) = -s(l,k1)
            s(l,k2) = -s(l,k2)
            s(l,k3) = -s(l,k3)
c .....................................................................       
  275     continue 
c .....................................................................
  270   continue
c .....................................................................
      endif
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c     if(block_pu) then
c       call lku_m(s,u,p,nst)
c     else
c       call lku_sym(s,u,p,nst)    
c     endif    
c .....................................................................
c
c .....................................................................      
      return  
c ======================================================================
c
c ... Tensoes nodais e fluxo nodais:
c
c ......................................................................
  300 continue
c ...                        
      gl(1)     =  gravity(1)
      gl(2)     =  gravity(2)
      gl(3)     =  gravity(3)
c ... 
      fluid_d   = e(7)*1.0d-06

c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod

c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ... 
      coef_biot = e(5)
      perm      = e(3)/fluid_sw
c ...
      a1        = 1.d0 - ps
      a2        = 1.d0 - 2.d0*ps
      a3        = 1.d0 + ps
c ...
      a         = (ym*a1)/(a3*a2)
      b         = ps/a1
      c         = 0.5d0*(a2/a1)
c .....................................................................
c
c ... tensao nodal total
      do 310 i = 1, 8
c       tp = (i-1)*6 + 1
        tp  = 6*i - 5
c ... tensao efetiva de biot (8x6+1=49)         
c       tp1 = (i-1)*6 + 49  
        tp1 = 6*i +43
c ...   p(1...8) = u(61...68) 
        l   = i  + 60 
c ... calculo do terminante
        call sfhexa8_m(hp,hpx,hpy,hpz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,8,nel ,.true.)
c ... calculo da derivadas das funcoes de interpolacao
        call sfhexa20_m(hu,hux,huy,huz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel ,.false.)
c .....................................................................
        call deform3d(hux,huy,huz,u,epsi,20)
        call stress3d(a,b,c,epsi,p(tp))
c ... Tensao = D*deformacao elastica - biot*dp
        dpm     = coef_biot*dp(i)
         pm     = coef_biot*u(l)
c ... tensao total
        p(tp)   = p(tp)   - dpm
        p(tp+1) = p(tp+1) - dpm
        p(tp+2) = p(tp+2) - dpm
c ... tensao efetiva de biot
        p(tp1)  = p(tp)   + pm
        p(tp1+1)= p(tp+1) + pm
        p(tp1+2)= p(tp+2) + pm
        p(tp1+3)= p(tp+3) 
        p(tp1+4)= p(tp+4) 
        p(tp1+5)= p(tp+5) 
c .....................................................................
  310 continue
c .....................................................................
c
c ... fluxo nodal 
      do 320 i = 1, 8
c ... fuxo (8x6 + 8x6=96)                 
c       tp = (i-1)*3 + 97
        tp = 3*i + 94
c ...
        call sfhexa8_m(hp,hpx,hpy,hpz,rn(i),sn(i),tn(i),.true.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c ... p nos pontos de integracao
        call darcy_flux(perm,gl,fluid_d,hpx,hpy,hpz,u(61),8,p(tp))  
c .....................................................................
  320 enddo
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume, no contorno e variaveis no passo
c     anterior:
c
c ......................................................................
 400  continue
c ...                        
      gl(1)     =  gravity(1)
      gl(2)     =  gravity(2)
      gl(3)     =  gravity(3)
c ...
      pm_d      = e(6)*1.0d-06             
      fluid_d   = e(7)*1.0d-06

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
      nint = 4 
      do 405 lz = 1, nint
        ti = pg(lz,nint)
        do 410 ly = 1, nint
          si = pg(ly,nint)
          do 415 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
            call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....................................................................
c
c ...
            call sfhexa20_m(hu,hux,huy,huz,ri,si,ti,.true.,.true.)
            call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel,.false.)
c .....................................................................
c
c ... deformacoes nos pontos de integracao
            call deform3d(hux,huy,huz,u,epsi,20)
c ... tensao nos pontos de integracao
            call stress3d(a,b,c,epsi,txi)
c .....................................................................            
c
c ... delta p nos pontos de integracao
            dpm = hp(1)*dp(1) + hp(2)*dp(2) + hp(3)*dp(3) + hp(4)*dp(4) 
     .          + hp(5)*dp(5) + hp(6)*dp(6) + hp(7)*dp(7) + hp(8)*dp(8) 
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
            do 420 i = 1, 20
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
              wt3   = hu(i)*wt2  
              p(l1) = p(l1) - wt3*gl(1)
              p(l2) = p(l2) - wt3*gl(2)
              p(l3) = p(l3) - wt3*gl(3)
  420       continue    
c .....................................................................
c
c ...
            wt1 = w*dt_fluid_perm 
            wt2 = w*dt_perm 
c .....................................................................
c
c ...
            do 425 i = 1, 8 
              l    = i + 60
c ... Fp = int(dt*perm*fuild_d*(B~T)*ge*dv)  
              p(l) = p(l) 
     .             + wt1*(hpx(i)*gl(1) + hpy(i)*gl(2) + hpz(i)*gl(3))
c .....................................................................
c
c ... Fp = int(dt*k*(B~T)*(B~)*pe*dv)
              tmp = 0.0d0
              do 430 j = 1, 8
                k    = j + 60
                tmp  = tmp       
     .          + wt2*u(k)*(hpx(i)*hpx(j)+hpy(i)*hpy(j)+hpz(i)*hpz(j))
  430         continue
              p(l) = p(l) - tmp
c .....................................................................
  425       continue
c .....................................................................
  415     continue
c .....................................................................
  410   continue
c .....................................................................
  405 continue
c .....................................................................
c
c .....................................................................
c
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 1 2 3 4  9 10 11 12 |
c             2 | no 5 6 7 8 13 14 15 16 |
c             3 | no 1 5 6 2 17 13 18  9 |     
c             4 | no 4 3 7 8 11 19 15 20 |
c             5 | no 1 4 8 5 12 20 16 17 |
c             6 | no 2 6 7 3 18 14 19 10 |
c
c ... verifica se ha alguma face com carga
      tp = 0
      do 435 i = 1, 6
        tp    = tp + iq(i)
  435 continue
c .....................................................................
c
c ...
      if( tp .gt. 0 ) then
        do 440 j = 1, 6
c ... face 
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 4
              no       = hexa_face_node20(i,j)
              xf(1,i) = x(1,no) 
              xf(2,i) = x(2,no) 
              xf(3,i) = x(3,no)
  445       continue
c ...    
            call rotmatrix(xf,r)
            call rotate(xf(1,1),r,xf(1,1),.false.)
            call rotate(xf(1,2),r,xf(1,2),.false.)
            call rotate(xf(1,3),r,xf(1,3),.false.)
            call rotate(xf(1,4),r,xf(1,4),.false.)
c ...
            carg = load(1,iq(j))
c ... forcas
            if( carg .eq. 40 ) then
              call tload(iq(j),0.d0,face_u,n_flux,face_f)
              nint = 3
              do 450 ly = 1, nint
                si = pg(ly,nint)
                do 455 lx = 1, nint
                  ri = pg(lx,nint)
c ...                                               
                  call sfquad4_m(hp,hpx,hpy,ri,si,.false.,.true.)
                  call jacob2d_m(hpx,hpy,xj2D,xji2D,xf,det,4,ndm
     .                        ,nel,.true.)
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
                  do 460 i = 1, 8
                    no  = hexa_face_node20(i,j)
c                   l1  = (no-1)*3+1
                    wt1 = hu(i)*w
                    l1    = 3*no-2
                    l2    = l1 + 1
                    l3    = l2 + 1
                    p(l1) = p(l1) - face_f(1)*wt1
                    p(l2) = p(l2) - face_f(2)*wt1
                    p(l3) = p(l3) - face_f(3)*wt1  
  460             continue
c .....................................................................
  455           continue
c .....................................................................             
  450         continue
c .....................................................................
c
c ... fluxo 
            elseif( carg .eq. 41) then
            endif
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
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
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
      nint = 4 
      do 510 lz = 1, nint
        ti = pg(lz,nint)
        do 520 ly = 1, nint
          si = pg(ly,nint)
          do 530 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(hp,hpx,hpy,hpz,ri,si,ti,.true.,.true.)
            call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
            wt1 = w  
c .....................................................................
c
c ...
            call sfhexa20_m(hu,hux,huy,huz,ri,si,ti,.false.,.true.)
            call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel,.false.)
c .....................................................................
c
c ...
            txi(1:6)= 0.d0
            do 540 j = 1,   8
              txi(1) = txi(1) + hp(j)*txn(1,j)
              txi(2) = txi(2) + hp(j)*txn(2,j) 
              txi(3) = txi(3) + hp(j)*txn(3,j) 
              txi(4) = txi(4) + hp(j)*txn(4,j) 
              txi(5) = txi(5) + hp(j)*txn(5,j)
              txi(6) = txi(6) + hp(j)*txn(6,j) 
  540       continue 
c ..................................................................... 
c
c ...
            do 560 i = 1, 20
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
c ...................................................................... 
  560       continue  
c .....................................................................
  530     continue
c .....................................................................
  520   continue
c .....................................................................
  510 continue
c .....................................................................
c
c ....
      return
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
