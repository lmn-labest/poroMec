c *********************************************************************
c * Biblioteca de elementos poro-mecanicos                            *
c * ----------------------------------------------------------------- *
c * Elastico:                                                         *
c * ----------------------------------------------------------------- *
c * ------------------ Elementos lineares --------------------------- *
c *                                                                   *
c * ELMT02_PM- triangulo de 3 nos para o problema poromecanico        *
c * elastico estatico (Estado plano de deformacao) (NAO IMPLEMENTADO) *
c *                                                                   *
c * ELMT03_PM - triangulo de 3 nos para o problema poromecanico       *
c * elastico estatico (Estado plano de deformacao) (NAO IMPLEMENTADO) *
c *                                                                   *
c * ELMT04_PM - quadrilateros de 4 nos para o problema poromecanico   * 
c * elastico estatico (Estado plano de deformacao) (NAO IMPLEMENTADO) *
c *                                                                   *
c * ELMT05_PM - quadrilateros de 4 nos para o problema poromecanico   *
c * elastico estatico (Estado plano de tensao) (NAO IMPLEMENTADO)     *
c *                                                                   *
c * ELMT06_PM - tetraedros de 4 nos para o problema poromecanico      * 
c * elastico estatico(NAO IMPLEMENTADO)                               *
c *                                                                   *
c * ELMT07_PM - hexaedros de  8 nos para o problema poromecanico      *
c * elastico estatico (NAO IMPLEMENTADO)                              *
c *                                                                   *
c * ------------------ Elementos quadraticos ------------------------ *
c *                                                                   *
c * ELMT16_PM - tetraedros de 10 nos para o problema poro mecanico    *
c * elastico                                                          *
c *                                                                   *
c * ELMT17_PM - hexaedros de 20 nos para o problema poro mecanico     *
c * elastico                                                          *
c *                                                                   *
c * ----------------------------------------------------------------- *
c * Plastico-estatico:                                                *
c * ----------------------------------------------------------------- *
c * ------------------ Elementos lineares --------------------------- *
c *                                                                   *
c * ------------------ Elementos quadraticos ------------------------ *
c *                                                                   *
c * ELMT36_PM - tetraedros de 10 nos para o problema poro mecanico    *
c * plastico (Nao implementado)                                       *
c *                                                                   *
c * ELMT37_PM - hexaedros de 20 nos para o problema poro mecanico     *
c * plastico (Nao implementado)                                       *
c *                                                                   *
c *********************************************************************  
      subroutine elmt16_pm(e,iq,x,u,dp,p,s,txn,ndm,nst,nel,isw
     .                    ,block_pu)
c **********************************************************************
c * Data de criacao    : 27/03/2016                                    *
c * Data de modificaco : 10/03/2017                                    * 
c * ------------------------------------------------------------------ *      
c * ELMT16_PM: Elemento tetraedrico de 10 nos para problemas           *  
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
      include 'transiente.fi'
      include 'gravity.fi'
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
      real*8 face_u(3),face_f(3),n_carg,ddum
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
      real*8 dt_c
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
      real*8 scale
      parameter (scale = 1.d-06)
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
      goto (100,200,300,400,500,600,700) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
c ...
      fluid_d   = e(7)*scale
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
      dt_c  = ((l_c*l_c)/perm) 
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
      fluid_d   = e(7)*scale
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
        p(i) = 0.d0
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
c ... -Kpp ( Int((1/M)*(N~t)*(N~)*dV) )    
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
c ... kpp=-Kpp
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
      fluid_d   = e(7)*scale
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
      pm_d      = e(6)*scale          
      fluid_d   = e(7)*scale
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
c .....................................................................
c
c ...
            carg = load(1,iq(j))
c ... forcas qualquer direcao ou normal
            if( carg .eq. 40 .or. carg .eq. 41) then
c ... calculo do verto normal externo a face o elemento
              if( carg .eq. 41) call face_normal_vector(xf,face_f,ndm)
c .....................................................................
c
c ... rotacionando os eixos    
              call rotmatrix(xf,r)
              call rotate(xf(1,1),r,xf(1,1),.false.)
              call rotate(xf(1,2),r,xf(1,2),.false.)
              call rotate(xf(1,3),r,xf(1,3),.false.)
c .....................................................................
c
c ... carga normal ao elemento
              if( carg .eq. 41) then
                call tload(iq(j),t,face_u,n_carg,ddum) 
                face_f(1:ndm) = n_carg*face_f(1:ndm)
c ... carga distribuida
              elseif ( carg .eq. 40) then
                call tload(iq(j),t,face_u,ddum,face_f)
              endif
c ......................................................................
c
c ...
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
            elseif( carg .eq. 42) then
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
        p(31:nst) =-p(31:nst)
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
c
c ... Variacao da porosidade             
c
c ......................................................................
c 
c ...  
  700 continue                      
c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ... 
      coef_biot = e(5)
      imod_biot = 1.d0/e(4)
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
c ... variacao de porosidade
      do 710 i = 1, 4
c ... calculo do terminante
        call sftetra4(hp,hpx,hpy,hpz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,4,nel ,.true.)
c ... calculo da derivadas das funcoes de interpolacao
        call sftetra10(hu,hux,huy,huz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,10,nel ,.false.)
c .....................................................................
        call deform3d(hux,huy,huz,u,epsi,10)
c .....................................................................
c
c ... variacao da porosidade 
c       p(1...4) = u(31...34) 
        l    = i  + 30 
        pi   = imod_biot*u(l)
c .....................................................................
c
c ... variacao da porosidade = biot*tr(deps)total + dp/M 
        p(i) = coef_biot*(epsi(1) + epsi(2) + epsi(3) ) +  pi
c ..................................................................... 
  710 continue
c .....................................................................
      return
      end
c *********************************************************************
      subroutine elmt17_pm(e,iq,x,u,dp,p,s,txn,ndm,nst,nel,isw
     .                    ,block_pu)
c **********************************************************************
c * Data de criacao    : 10/12/2015                                    *
c * Data de modificaco : 04/03/2017                                    * 
c * ------------------------------------------------------------------ *       
c * ELMT17_PM: Elemento hexaedricos de 20 nos para problemas           *  
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
c *  7 = variacao da porosidade                                        *
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
c *     isw = 7  porosidade                                            *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * u(1:60) - isw = 2 diferenca de deslocamentos                       *
c * u(1:60) - isw = 3 e isw = 4 deslocamentos                          *
c * u(61:68)- isw = 2 diferenca de pressao                             *
c * u(61:68)- isw = 3 e isw = 4 pressao                                *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      include 'gravity.fi'
      include 'load.fi'
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k,tp,tp1
      integer nen,nint,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_u(8),face_f(3),n_carg,ddum
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
      real*8 dt_c
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
      real*8 scale
      parameter (scale = 1.d-06)
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
      goto (100,200,300,400,500,600,700) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
c ...
      fluid_d   = e(7)*scale
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
      dt_c  = ((l_c*l_c)/perm) 
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
      fluid_d   = e(7)*scale
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
        p(i) = 0.d0
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
      fluid_d   = e(7)*scale

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
c ... fuxo (8x6 + 8x6 + 1 = 97)                 
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
      pm_d      = e(6)*scale            
      fluid_d   = e(7)*scale

c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod
c ...      
      ym        = e(1)
      ps        = e(2)
c ...
      coef_biot = e(5)
      perm      =  e(3)/fluid_sw
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
c ... Fu = int(pm_d*(N~T)*ge*dv)
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
c     iq(1) = 1 | no 1 2 3 4  9 10 11 12 | normal externa
c             2 | no 5 6 7 8 13 14 15 16 | normal interna
c             3 | no 1 5 6 2 17 13 18  9 | normal externa    
c             4 | no 4 3 7 8 11 19 15 20 | normal externa
c             5 | no 1 4 8 5 12 20 16 17 | normal externa
c             6 | no 2 6 7 3 18 14 19 10 | normal externa
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
c .....................................................................
c
c ...
            carg = load(1,iq(j))
c ... forcas
            if( carg .eq. 40 .or. carg .eq. 41) then
c ... calculo do verto normal externo a face o elemento
              if( carg .eq. 41) then
                call face_normal_vector(xf,face_f,ndm)
                if (j .ne. 2) face_f(1:ndm) = -1.0d0*face_f(1:ndm)
              endif
c .....................................................................
c          
c ... rotacionando os eixos 
              call rotmatrix(xf,r)
              call rotate(xf(1,1),r,xf(1,1),.false.)
              call rotate(xf(1,2),r,xf(1,2),.false.)
              call rotate(xf(1,3),r,xf(1,3),.false.)
              call rotate(xf(1,4),r,xf(1,4),.false.)
c ...................................................................
c
c ... carga normal ao elemento
              if( carg .eq. 41) then
                call tload(iq(j),t,face_u,n_carg,ddum) 
                face_f(1:ndm) = n_carg*face_f(1:ndm)
c ... carga distribuida
              elseif ( carg .eq. 40) then
                call tload(iq(j),t,face_u,ddum,face_f)
              endif
c ......................................................................
c
c ...
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
c ======================================================================
c
c ... Variacao da porosidade             
c
c ......................................................................
c 
c ...  
  700 continue                      
c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ... 
      coef_biot = e(5)
      imod_biot = 1.d0/e(4)
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
c ... variacao de porosidade
      do 710 i = 1, 8
c ... calculo do terminante
        call sfhexa8_m(hp,hpx,hpy,hpz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,8,nel ,.true.)
c .....................................................................
c
c ... calculo da derivadas das funcoes de interpolacao
        call sfhexa20_m(hu,hux,huy,huz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel ,.false.)
c .....................................................................
c
c ... variacao da deformacao
        call deform3d(hux,huy,huz,u,epsi,20)
c ......................................................................
c
c ... variacao da porosidade = biot (eps11 + eps22 + eps33) + dp/M
c ...   p(1...8) = u(61...68) 
        l    = i  + 60 
        pi   = imod_biot*u(l)
c .....................................................................
c
c ... variacao da porosidade = biot*tr(deps)total + dp/M 
        p(i) = coef_biot*(epsi(1) + epsi(2) + epsi(3) ) + pi
c ..................................................................... 
  710 continue
c .....................................................................
      return
c ======================================================================
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine elmt36_pm(e       ,iq       ,x       ,u    ,p0      
     1                    ,p       ,s        ,tx0     ,tx   ,depsi
     2                    ,vplastic,elplastic,ndm     ,nst  ,nel
     3                    ,isw     ,block_pu ,nlit)
c **********************************************************************
c * Data de criacao    : 22/12/2016                                    *
c * Data de modificaco : 10/03/2017                                    * 
c * ------------------------------------------------------------------ *      
c * ELMT36_PM: Elemento tetraedrico de 10 nos para problemas           *  
c * poromecanico plastico                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(12) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = coeficiente de Darcy                              *
c *           e(4) = modulo de Biot                                    *
c *           e(5) = coeficiente de Biot                               *
c *           e(6) = massa especifica homogenizada do meio poroso      *
c *           e(7) = massa especifica do fluido                        *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u + p)                *
c * p0(*)      - poropressao do passo de tempo anterior                *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * tx0(6,*)   - tensao nos pontos de integracao no passo de tempo     *
c *              anterior                                              *
c * tx(6,*)    - tensao nos pontos de integracao                       *
c * depsi(7,*) - incremento de deformascoe e poropressoes              *
c * vplastic(3,*)- deformacao volumetricas plasticas no passo de tempo *
c *                anterior                                            *
c *                deformacao volumetricas plasticas                   *
c *                paramentro de endurecimento nos pontos de integracao*
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
c *  7 = variacao da porosidades                                       *
c *  8 = extrapolacao das pressoes de consolidacao nos vertices        *
c *  9 = incializacao da pressao de consolidacao nos pontos de         *
c *  integracao                                                        *
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
c * vplastic(1,*) - dilatacao volumetrica plastica do passo de tempo   *
c *                  anterior                                          *
c * vplastic(2,*) - dilatacao volumetrica plastica                     *
c * vplastic(3,*)  - paramentro de endurecimento                       *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      include 'gravity.fi'
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
      real*8 face_u(3),face_f(3),n_carg
      real*8 ddum
      integer carg
c ...
      real*8 u(*),p0(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 hu(10),hux(10),huy(10),huz(10)
      real*8 hp(4),hpx(4),hpy(4),hpz(4)
c ...
      real*8 xj(3,3),xji(3,3),xj2D(2,2),xji2D(2,2),r(3,3)
      real*8 ri,si,ti,w,wt1,wt2,wt3,det
      real*8 rn(10),sn(10),tn(10)
      real*8 ddepsi(6),epsi(6),txi(6)
      real*8 vplastic(3,*),tx0(6,*),tx(6,*),depsi(7,*)
      real*8 pi,dpi,ddpi
c ...
      real*8 etx(6,4),pce(4),pci(4)
c ... integracao numerica de tetraedros     
      real*8 pri4(5,3),psi4(5,3),pti4(5,3),wf4(5,3)
      real*8  pri(12,5),psi(12,5),wf(12,5)
      integer igrau,igrau_vol,igrau_face,npint(5),npint4(3)
c ...
      real*8 dt_c
      real*8 e(*),x(ndm,*),xf(3,3)
c ...
      real*8 l_c
      real*8 perm,a,b,c,ym,ps,e0
      real*8 fluid_d,dt_perm,pm_d
      real*8 dt_fluid,dt_fluid_perm,lambda,mi,gl(3)
      real*8 imod_biot,coef_biot
      real*8 fluid_sw
      real*8 a1,a2,a3,tmp
c ... plasticidade
      real*8 alpha_exp,pc0,mcs,c14,g11,def_vol_plast,lambda_plastic
      real*8 k_plastic,cam_clay_pc
      integer elplastic
      logical sup 
c ... 
      integer tetra_face_node10(6,4),no
      real*8 tetra_vol,volum
c ...
      real*8 scale
      parameter (scale = 1.d-06)
c ... 
      integer nlit
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
      data tn / 0.0d0, 0.0d0, 1.0d0 ! t1, t2, t3
     .        , 0.0d0, 0.5d0, 0.0d0 ! t4, t5, t6              
     .        , 0.5d0, 0.5d0, 0.0d0 ! t7, t8, t9              
     .        , 0.0d0/              ! t10              
c
      data nen/10/,igrau_vol/2/,igrau_face/2/    
c ......................................................................
      goto (100,200,300,400,500,600,700,800,900) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
c ...
      fluid_d   = e(7)*scale
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
      dt_c  = ((l_c*l_c)/perm) 
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
      fluid_d   = e(7)*scale
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
c ... plasticidade
      e0             = e(8) 
      lambda_plastic = e(9)
      k_plastic      = e(10)
      mcs            = e(11)
      alpha_exp = (1+e0)/(lambda_plastic-k_plastic)
c ...
      c14 = ps/a2
      g11 = 0.5d0*(ym/a3)
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
        p(i) = 0.d0
      enddo
c .....................................................................
c
c ...
      elPlastic = 0
      sup       =.false.
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
c ... -Kpp ( Int((1/M)*(N~t)*(N~)*dV) )    
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
c ... Kpp = -kpp
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
      if(nlit .eq. 1 ) then
        call lku_m(s,u,p,nst)
c
        igrau = igrau_vol 
        nint  = npint4(igrau) 
        do lx = 1, nint
c ... zera o delta de deformacoes e o delta de pressoes 
c     nos pontos de integracoes
          depsi(1:7,lx) = 0.d0
        enddo
        return  
      endif
c ......................................................................
c
c ...    
      igrau = igrau_vol 
      nint  = npint4(igrau) 
      do 270 lx = 1, nint
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
c ... incremetno de deformacoes nos pontos de integracao
        call deform3d(hux,huy,huz,u,epsi,10)
c .....................................................................           
c
c ... incremetno de p nos pontos de integracao
        dpi = hp(1)*u(31) + hp(2)*u(32) + hp(3)*u(33) + hp(4)*u(34) 
c .....................................................................
c
c ... p do passo anterior nos pontos de integracao
        pi = hp(1)*p0(1) + hp(2)*p0(2) + hp(3)*p0(3) + hp(4)*p0(4) 
c .....................................................................
c
c ... ddeps(i,n+1) = deps(i,n+1) -  deps(i-1,n+1) 
        ddepsi(1:6) = epsi(1:6) - depsi(1:6,lx)
c ... ddpi(i,n+1) = dp(i,n+1) -  dp(i-1,n+1) 
        ddpi       = dpi       - depsi(7,lx)   
c .....................................................................
c
c ... deps(i,n+1)
        depsi(1:6,lx) = epsi(1:6)
c ... dp(i,n+1)
        depsi(7,lx)   = dpi       
c .....................................................................
c
c ... tx = D(ddeps)
        call stress3d(a,b,c,ddepsi,txi)
c .....................................................................
c
c ... incremento do incremento de tensao efetiva nos pontos de integracao
        tmp      = (1.0d0-coef_biot)*ddpi
        txi(1:3) = txi(1:3) + tmp        
c .....................................................................
c
c ... tensao efetiva:
        if(nlit .eq. 2) then
          txi(1:3) = tx(1:3,lx) + pi + txi(1:3) 
        else
          txi(1:3) = tx(1:3,lx) + pi + dpi + txi(1:3) 
        endif
        txi(4:6) = tx(4:6,lx) + txi(4:6)
c .....................................................................
c
c ... 
        call plasticity3d_pm(vplastic(2,lx) ,vplastic(3,lx),txi
     1                          ,mcs              ,alpha_exp       ,ps 
     2                          ,c14              ,g11             ,sup 
     3                          ,nel)  
c .....................................................................
c
c ...         
        if(sup) elPlastic = 1
c .....................................................................
c
c ... tensao total 
        if(nlit .eq. 2) then
          txi(1:3) = txi(1:3) - (pi+ddpi)
        else 
          txi(1:3) = txi(1:3) - (pi+dpi+ddpi)
        endif  
c ......................................................................
c
c ... tx(i,n+1)  
        tx(1:6,lx) = txi(1:6)
c .....................................................................
c
c ... incremento da tensao total ( sigma(n+1,i) - sigma(n) )
        txi(1:6) = tx(1:6,lx) - tx0(1:6,lx)
c .....................................................................
c
c ...
        wt1 = w  
c .....................................................................
c
c ...
        do 275 i = 1, 10
c         l1  = (i-1)*3+1
          l1  = 3*i-2
          l2  = l1 + 1
          l3  = l2 + 1
c ... Fu = Fu - p = Fu - int(BeT*sigma*dv) 
          p(l1) = p(l1) 
     .    + wt1*(hux(i)*txi(1) + huy(i)*txi(4) + huz(i)*txi(6))
c          
          p(l2) = p(l2) 
     .    + wt1*(hux(i)*txi(4) + huy(i)*txi(2) + huz(i)*txi(5)) 
c          
          p(l3) = p(l3) 
     .    + wt1*(hux(i)*txi(6) + huy(i)*txi(5) + huz(i)*txi(3))
c .....................................................................
  275   continue   
c ..................................................................... 
c
c ...
        wt1 = w*(1.0d0-coef_biot)*(vplastic(2,lx)-vplastic(1,lx))       
c .....................................................................
c
c ...
        do i = 1, 4 
          l    = i + 30
c ... Fp =  Fp - int((1-coef_biot)*trace(delta_def_plast)*(N~)dv)  
          p(l) = p(l) - wt1*hp(i)
c .....................................................................
        enddo
c .....................................................................
 270  continue   
c .....................................................................
c
c ... p = [kpu kpp ] [u p] 
      do j = 1, nst
        p(31) = p(31) + s(31,j)*u(j)
        p(32) = p(32) + s(32,j)*u(j)
        p(33) = p(33) + s(33,j)*u(j)
        p(34) = p(34) + s(34,j)*u(j)
      enddo   
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
      fluid_d   = e(7)*scale
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
c ... extrapolacao da tensao do ponto de integracao 
c     para o nos dos vertice
      call extrapol_tetra10(tx,etx,rn,sn,tn,6)
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
c .....................................................................
c
c ... tensao total
       p(tp)   = etx(1,i)           
       p(tp+1) = etx(2,i) 
       p(tp+2) = etx(3,i) 
       p(tp+3) = etx(4,i) 
       p(tp+4) = etx(5,i) 
       p(tp+5) = etx(6,i)  
c .....................................................................
c
c ...
        pi      = coef_biot*u(l)
c ... tensao efetiva de biot
        p(tp1)  = p(tp)   + pi 
        p(tp1+1)= p(tp+1) + pi
        p(tp1+2)= p(tp+2) + pi
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
      pm_d      = e(6)*scale         
      fluid_d   = e(7)*scale
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
c ... tensao nos pontos de integracao
        txi(1:6) = tx0(1:6,lx)
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
c .....................................................................
c
c ...
            carg = load(1,iq(j))
c ... forcas qualquer direcao ou normal
            if( carg .eq. 40 .or. carg .eq. 41) then
c ... calculo do verto normal externo a face o elemento
              if( carg .eq. 41) call face_normal_vector(xf,face_f,ndm)
c .....................................................................
c
c ... rotacionando o eixos    
              call rotmatrix(xf,r)
              call rotate(xf(1,1),r,xf(1,1),.false.)
              call rotate(xf(1,2),r,xf(1,2),.false.)
              call rotate(xf(1,3),r,xf(1,3),.false.)
c .....................................................................
c
c ... carga normal ao elemento
              if( carg .eq. 41) then
                call tload(iq(j),t,face_u,n_carg,ddum) 
                face_f(1:ndm) = n_carg*face_f(1:ndm)
c ... carga distribuida
              elseif ( carg .eq. 40) then
                call tload(iq(j),t,face_u,ddum,face_f)
              endif
c ......................................................................
c
c ...
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
c ... fluxo de massa
            elseif( carg .eq. 42) then             
c .....................................................................
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
        p(31:nst) =-p(31:nst)
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
        call sftetra4(hp,hpx,hpy,hpz,ri,si,ti,.true.,.false.)
c .....................................................................
c  
c ...
        txi(1:6)= 0.d0
        do 515 j = 1,   4
          txi(1) = txi(1) + hp(j)*tx0(1,j)
          txi(2) = txi(2) + hp(j)*tx0(2,j) 
          txi(3) = txi(3) + hp(j)*tx0(3,j) 
          txi(4) = txi(4) + hp(j)*tx0(4,j) 
          txi(5) = txi(5) + hp(j)*tx0(5,j)
          txi(6) = txi(6) + hp(j)*tx0(6,j) 
  515   continue
c .....................................................................
c 
c ...
        tx(1:6,lx) = txi(1:6)       
c .....................................................................
  510 continue 
c ....................................................................
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      stop
c ======================================================================
c
c ... Variacao da porosidade             
c
c ......................................................................
c 
c ...  
  700 continue                      
c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ... 
      coef_biot = e(5)
      imod_biot = 1.d0/e(4)
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
c ... delta da variacao volumetrica plastica nos nos
      igrau       = igrau_vol 
      nint        = npint4(igrau) 
      pci(1:nint) = vplastic(2,1:nint) - vplastic(1,1:nint)
c ... pressao de consolidacao nodal total
      call extrapol_tetra10(pci,pce,rn,sn,tn,1)
c .....................................................................
c
c ... variacao de porosidade
      do 710 i = 1, 4
c ... calculo do terminante
        call sftetra4(hp,hpx,hpy,hpz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,4,nel ,.true.)
c ... calculo da derivadas das funcoes de interpolacao
        call sftetra10(hu,hux,huy,huz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,10,nel ,.false.)
c .....................................................................
        call deform3d(hux,huy,huz,u,epsi,10)
c ....................................................................
c
c ... variacao da porosidade 
c       p(1...4) = u(31...34) 
        l    = i  + 30 
        pi   = imod_biot*u(l)
c .....................................................................
c
c ... variacao da porosidade = biot*tr(deps)total +
c     (1-b)*tr(deps)pastic + dp/M 
        p(i) = coef_biot*(epsi(1) + epsi(2) + epsi(3) ) 
     .       + (1.d0-coef_biot)*pce(i) +  pi
c ..................................................................... 
  710 continue
c .....................................................................
      return
c ======================================================================
c
c ... extrapola das pressoes de consolidacao para os nos
c
c ......................................................................
  800 continue
      igrau       = igrau_vol 
      nint        = npint4(igrau) 
      pci(1:nint) = vplastic(3,1:nint)
c ... pressao de consolidacao nodal total
      call extrapol_tetra10(pci,pce,rn,sn,tn,1)
      p(1:4)      = pce(1:4)    
c .....................................................................
      return
c ======================================================================
c
c ======================================================================
c
c ... calculo do parametro de encruamento                
c
c ......................................................................
  900 continue
c ...     
      pc0  = e(12)
      mcs  = e(11)
      igrau = igrau_vol 
      nint  = npint4(igrau) 
      do 910 lx = 1, nint
c ...
        ti = pti4(lx,igrau)
        si = psi4(lx,igrau)
        ri = pri4(lx,igrau)
c ...                                               
        call sftetra4(hp,hpx,hpy,hpz,ri,si,ti,.true.,.false.)
c .....................................................................
c  
c ...
        pi       = 0.d0
        txi(1:6) = 0.d0
        do 915 j = 1,   4
          l      = j  + 30 
          txi(1) = txi(1) + hp(j)*tx0(1,j)
          txi(2) = txi(2) + hp(j)*tx0(2,j) 
          txi(3) = txi(3) + hp(j)*tx0(3,j) 
          txi(4) = txi(4) + hp(j)*tx0(4,j) 
          txi(5) = txi(5) + hp(j)*tx0(5,j)
          txi(6) = txi(6) + hp(j)*tx0(6,j) 
          pi     = pi     + hp(j)*u(l)
  915   continue
c .....................................................................
c
c ... tensao efetiva
        txi(1:3) =  txi(1:3) + pi
c ..................................................................... 
c
c ...
        vplastic(3,lx) = max(pc0,cam_clay_pc(txi,mcs))
c ..................................................................... 
  910 continue 
c .....................................................................
      return
c =====================================================================
      end
c *********************************************************************
c
c *********************************************************************
      subroutine elmt37_pm(e    ,iq      ,x       ,u    ,p0      
     1                    ,p    ,s       ,tx0     ,tx   ,depsi
     2                    ,vplastic,elplastic,ndm ,nst  ,nel
     3                    ,isw     ,block_pu ,nlit)
c **********************************************************************
c * Data de criacao    : 10/12/2015                                    *
c * Data de modificaco : 20/03/2017                                    * 
c * ------------------------------------------------------------------ *       
c * ELMT37_PM: Elemento hexaedricos de 20 nos para problemas           *  
c * poromecanico plastico                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(12) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = coeficiente de Darcy                              *
c *           e(4) = modulo de Biot                                    *
c *           e(5) = coeficiente de Biot                               *
c *           e(6) = massa especifica homogenizada do meio poroso      *
c *           e(7) = massa especifica do fluido                        *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u + p)                *
c * p0(*)      - poropressao do passo de tempo anterior                *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * tx0(6,*)   - tensao nos pontos de integracao no passo de tempo     *
c *              anterior                                              *
c * tx(6,*)    - tensao nos pontos de integracao                       *
c * depsi(7,*) - incremento de deformascoe e poropressoes              *
c * vplastic(3,*)- deformacao volumetricas plasticas no passo de tempo *
c *                anterior                                            *
c *                deformacao volumetricas plasticas                   *
c *                paramentro de endurecimento nos pontos de integracao*
c * eplastic - identificao se o elemento plastificou ou nao (0 ou 1)   *
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
c * eplastic - identificao se o elemento plastificou ou nao (0 ou 1)   *
c * nlit     - iteracao nao linear                                     *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao e fluxo                                        *
c *     isw = 4  cargas de superfice, volume e integras do passo       *
c *     de tempo anterior                                              *
c *     isw = 7  porosidade                                            *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * u(1:60) - isw = 2 diferenca de deslocamentos                       *
c * u(1:60) - isw = 3 e isw = 4 deslocamentos                          *
c * u(61:68)- isw = 2 diferenca de pressao                             *
c * u(61:68)- isw = 3 e isw = 4 pressao                                *
c * vplastic(1,*) - dilatacao volumetrica plastica do passo de tempo   *
c *                  anterior                                          *
c * vplastic(2,*) - dilatacao volumetrica plastica                     *
c * vplastic(3,*)  - paramentro de endurecimento                       *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      include 'gravity.fi'
      include 'load.fi'
      common /gauss/ pg, wg
      integer nint
      parameter (nint = 4)
c ...
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k,tp,tp1,inpi
      integer nen,nints,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_u(8),face_f(3),n_carg,ddum
      integer carg
c ...
      real*8 u(*),p0(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 hu(20),hux(20),huy(20),huz(20)
      real*8 hp(8),hpx(8),hpy(8),hpz(8)
c ...
      real*8 xj(3,3),xji(3,3),xj2D(2,2),xji2D(2,2),r(3,3)
      real*8 ri,si,ti,w,wt1,wt2,wt3,det
      real*8 rn(20),sn(20),tn(20)
      real*8 ddepsi(6),epsi(6),txi(6)
      real*8 vplastic(3,*),tx0(6,*),tx(6,*),depsi(7,*)
      real*8 pi,dpi,ddpi
c ...
      real*8 etx(6,8),pce(8),pci(64)
c ... integracao numerica de tetraedros          
      real*8 pg(10,10),wg(10,10)
c ...
      real*8 dt_c
      real*8 e(*),x(ndm,*),xf(3,4)
c ...
      real*8 l_c
      real*8 perm,a,b,c,ym,ps,e0
      real*8 fluid_d,dt_perm,pm_d
      real*8 dt_fluid,dt_fluid_perm,lambda,mi,gl(3)
      real*8 imod_biot,coef_biot
      real*8 fluid_sw
      real*8 a1,a2,a3,tmp
c ... plasticidade
      real*8 alpha_exp,pc0,mcs,c14,g11,def_vol_plast,lambda_plastic
      real*8 k_plastic,cam_clay_pc
      integer elplastic
      logical sup 
c ... 
      integer hexa_face_node20(8,6),no
      real*8 hexa_vol,volum
c ...
      real*8 scale
      parameter (scale = 1.d-06)
c ...
      integer nlit
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
c 
c ......................................................................
      goto (100,200,300,400,500,600,700,800,900) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
c ...
      fluid_d   = e(7)*scale

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
      dt_c  = ((l_c*l_c)/perm) 
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
      fluid_d   = e(7)*scale  
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
c ... plasticidade
      e0             = e(8) 
      lambda_plastic = e(9)
      k_plastic      = e(10)
      mcs            = e(11)
      alpha_exp      = (1+e0)/(lambda_plastic-k_plastic)
c ...
      c14 = ps/a2
      g11 = 0.5d0*(ym/a3)
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
        p(i) = 0.d0
      enddo
c .....................................................................
c
c ...
      elPlastic = 0
      sup       =.false.
c .....................................................................
c
c ... 
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
      if(nlit .eq. 1 ) then
        call lku_m(s,u,p,nst)
c        
        do i = 1, nint*nint*nint
c ... zera o delta de deformacoes e o delta de pressoes 
c     nos pontos de integracoes
          depsi(1:7,i) = 0.d0          
        enddo
        return
      endif
c .....................................................................
c
c ... 
      inpi = 0 
      do 276 lz = 1, nint
        ti = pg(lz,nint)
        do 280 ly = 1, nint
          si = pg(ly,nint)
          do 285 lx = 1, nint
            inpi = inpi + 1
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
c ... incremetno de deformacoes nos pontos de integracao
            call deform3d(hux,huy,huz,u,epsi,20)
c .....................................................................           
c
c ... incremetno de p nos pontos de integracao
            dpi = hp(1)*u(61) + hp(2)*u(62) 
     .          + hp(3)*u(63) + hp(4)*u(64) 
     .          + hp(5)*u(65) + hp(6)*u(66) 
     .          + hp(7)*u(67) + hp(8)*u(68) 
c .....................................................................
c
c ... p do passo anterior nos pontos de integracao
            pi = hp(1)*p0(1) + hp(2)*p0(2) 
     .         + hp(3)*p0(3) + hp(4)*p0(4) 
     .         + hp(5)*p0(5) + hp(6)*p0(6) 
     .         + hp(7)*p0(7) + hp(8)*p0(8) 
c .....................................................................
c
c ... ddeps(i,n+1) = deps(i,n+1) -  deps(i-1,n+1) 
            ddepsi(1:6) = epsi(1:6) - depsi(1:6,inpi)
c ... ddpi(i,n+1) = dp(i,n+1) -  dp(i-1,n+1) 
            ddpi        = dpi       - depsi(7,inpi)   
c .....................................................................
c
c ... deps(i,n+1)
            depsi(1:6,inpi) = epsi(1:6)
c ... dp(i,n+1)
            depsi(7,inpi)   = dpi       
c .....................................................................
c
c ... tx = D(ddeps)
            call stress3d(a,b,c,ddepsi,txi)
c .....................................................................
c
c ... incremento do incremento de tensao efetiva nos pontos de integracao
            tmp      = (1.0d0-coef_biot)*ddpi
            txi(1:3) = txi(1:3) + tmp        
c .....................................................................
c
c ... tensao efetiva:
            if(nlit .eq. 2) then
              txi(1:3) = tx(1:3,inpi) + pi + txi(1:3) 
            else
              txi(1:3) = tx(1:3,inpi) + pi + dpi + txi(1:3) 
            endif
            txi(4:6) = tx(4:6,inpi) + txi(4:6)
c .....................................................................
c
c ... 
            call plasticity3d_pm(vplastic(2,inpi) ,vplastic(3,inpi),txi
     1                          ,mcs              ,alpha_exp       ,ps 
     2                          ,c14              ,g11             ,sup 
     3                          ,nel)  
c .....................................................................
c
c ...         
            if(sup) elPlastic = 1
c .....................................................................
c
c ... tensao total 
            if(nlit .eq. 2) then
              txi(1:3) = txi(1:3) - (pi+ddpi)
            else 
              txi(1:3) = txi(1:3) - (pi+dpi+ddpi)
            endif  
c ......................................................................
c
c ... tx(i,n+1)  
            tx(1:6,inpi) = txi(1:6)
c .....................................................................
c
c ... incremento da tensao total ( sigma(n+1,i) - sigma(n) )
            txi(1:6) = tx(1:6,inpi) - tx0(1:6,inpi)
c .....................................................................
c
c ...
            wt1 = w
c .....................................................................
c
c ...
            do 290 i = 1, 20
c             l1  = (i-1)*3+1
              l1  = 3*i-2
              l2  = l1 + 1
              l3  = l2 + 1
c ... Fu = int(BeT*sigma*dv)
              p(l1) = p(l1) 
     .        + wt1*(hux(i)*txi(1) + huy(i)*txi(4) + huz(i)*txi(6)) 
c
              p(l2) = p(l2) 
     .        + wt1*(hux(i)*txi(4) + huy(i)*txi(2) + huz(i)*txi(5)) 
c
              p(l3) = p(l3) 
     .        + wt1*(hux(i)*txi(6) + huy(i)*txi(5) + huz(i)*txi(3))
c .....................................................................
  290       continue    
c .....................................................................
c
c ...
            wt1 =w*(1.0d0-coef_biot)*(vplastic(2,inpi)-vplastic(1,inpi))
c .....................................................................
c
c ...
            do i = 1, 8 
              l    = i + 60
c ... Fp =  Fp - int((1-coef_biot)*trace(delta_def_plast)*(N~)dv)  
              p(l) = p(l) - wt1*hp(i)
c .....................................................................
            enddo
c .....................................................................            
  285     continue 
c .....................................................................            
  280   continue 
c .....................................................................              
  276 continue
c .....................................................................     
c
c ... p = [kpu kpp ] [u p] 
      do j = 1, nst
        p(61) = p(61) + s(61,j)*u(j)
        p(62) = p(62) + s(62,j)*u(j)
        p(63) = p(63) + s(63,j)*u(j)
        p(64) = p(64) + s(64,j)*u(j)
        p(65) = p(65) + s(65,j)*u(j)
        p(66) = p(66) + s(66,j)*u(j)
        p(67) = p(67) + s(67,j)*u(j)
        p(68) = p(68) + s(68,j)*u(j)
      enddo   
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
      fluid_d   = e(7)*scale

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
c ...    
      call extrapol_hexa20_v2(tx,etx,rn,sn,tn,6)
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
c ... 
        pi     = coef_biot*u(l)
c ... tensao total no nos
        p(tp)   = etx(1,i)           
        p(tp+1) = etx(2,i) 
        p(tp+2) = etx(3,i) 
        p(tp+3) = etx(4,i) 
        p(tp+4) = etx(5,i) 
        p(tp+5) = etx(6,i)  
c ... tensao efetiva de biot
        p(tp1)  = p(tp)   + pi
        p(tp1+1)= p(tp+1) + pi
        p(tp1+2)= p(tp+2) + pi
        p(tp1+3)= p(tp+3) 
        p(tp1+4)= p(tp+4) 
        p(tp1+5)= p(tp+5) 
c .....................................................................
  310 continue
c .....................................................................
c
c ... fluxo nodal 
      do 320 i = 1, 8
c ... fuxo (8x6 + 8x6 + 1 = 97)                 
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
      pm_d      = e(6)*scale            
      fluid_d   = e(7)*scale

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
      inpi = 0                          
      do 405 lz = 1, nint
        ti = pg(lz,nint)
        do 410 ly = 1, nint
          si = pg(ly,nint)
          do 415 lx = 1, nint
            ri = pg(lx,nint)
c ...
            inpi = inpi + 1
c .....................................................................
c
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
c ... tensao nos pontos de integracao
            txi(1:6) = tx0(1:6,inpi)
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
c
              p(l2) = p(l2) 
     .        + wt1*(hux(i)*txi(4) + huy(i)*txi(2) + huz(i)*txi(5)) 
c
              p(l3) = p(l3) 
     .        + wt1*(hux(i)*txi(6) + huy(i)*txi(5) + huz(i)*txi(3))
c .....................................................................
c
c ... Fu = int(pm_d*(N~T)*ge*dv)
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
c .....................................................................
c
c ...
            carg = load(1,iq(j))
c ... forcas
            if( carg .eq. 40 .or. carg .eq. 41) then
c ... calculo do verto normal externo a face o elemento
              if( carg .eq. 41) then
                call face_normal_vector(xf,face_f,ndm)
                if (j .ne. 2) face_f(1:ndm) = -1.0d0*face_f(1:ndm)
              endif
c .....................................................................
c          
c ... rotacionando os eixos 
              call rotmatrix(xf,r)
              call rotate(xf(1,1),r,xf(1,1),.false.)
              call rotate(xf(1,2),r,xf(1,2),.false.)
              call rotate(xf(1,3),r,xf(1,3),.false.)
              call rotate(xf(1,4),r,xf(1,4),.false.)
c ...................................................................
c
c ... carga normal ao elemento
              if( carg .eq. 41) then
                call tload(iq(j),t,face_u,n_carg,ddum) 
                face_f(1:ndm) = n_carg*face_f(1:ndm)
c ... carga distribuida
              elseif ( carg .eq. 40) then
                call tload(iq(j),t,face_u,ddum,face_f)
              endif
c ......................................................................
c
c ...
              nints = 3
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
                    no    = hexa_face_node20(i,j)
c                   l1  = (no-1)*3+1
                    wt1   = hu(i)*w
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
            elseif( carg .eq. 42) then
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
      inpi  = 0
      do 510 lz = 1, nint
        ti = pg(lz,nint)
        do 520 ly = 1, nint
          si = pg(ly,nint)
          do 530 lx = 1, nint
            inpi = inpi + 1
            ri   = pg(lx,nint)
c ...                                               
            call sfhexa8_m(hp,hpx,hpy,hpz,ri,si,ti,.true.,.false.)
c .....................................................................
c
c ...
            txi(1:6)= 0.d0
            do 540 j = 1,   8
              txi(1) = txi(1) + hp(j)*tx0(1,j)
              txi(2) = txi(2) + hp(j)*tx0(2,j) 
              txi(3) = txi(3) + hp(j)*tx0(3,j) 
              txi(4) = txi(4) + hp(j)*tx0(4,j) 
              txi(5) = txi(5) + hp(j)*tx0(5,j)
              txi(6) = txi(6) + hp(j)*tx0(6,j) 
  540       continue 
c ..................................................................... 
c
c ...
            tx(1:6,inpi) = txi(1:6) 
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
c ======================================================================
c
c ... Variacao da porosidade             
c
c ......................................................................
c 
c ...  
  700 continue                      
c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ... 
      coef_biot = e(5)
      imod_biot = 1.d0/e(4)
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
c ... delta da variacao volumetrica plastica nos nos
      inpi        = nint*nint*nint
      pci(1:inpi) = vplastic(2,1:inpi) - vplastic(1,1:inpi)
      call extrapol_hexa20_v2(pci,pce,rn,sn,tn,1)
c .....................................................................
c
c ... variacao de porosidade
      do 710 i = 1, 8
c ... calculo do terminante
        call sfhexa8_m(hp,hpx,hpy,hpz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,8,nel ,.true.)
c .....................................................................
c
c ... calculo da derivadas das funcoes de interpolacao
        call sfhexa20_m(hu,hux,huy,huz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel ,.false.)
c .....................................................................
c
c ... variacao da deformacao
        call deform3d(hux,huy,huz,u,epsi,20)
c ......................................................................
c
c ... variacao da porosidade
c       p(1...8) = u(61...68) 
        l    = i  + 60 
        pi   = imod_biot*u(l)
c .....................................................................
c
c ... variacao da porosidade = biot*tr(deps)total +
c     (1-b)*tr(deps)pastic + dp/M 
        p(i) = coef_biot*(epsi(1) + epsi(2) + epsi(3) ) 
     .       + (1.d0-coef_biot)*pce(i) +  pi
c ..................................................................... 
  710 continue
c .....................................................................
      return
c ======================================================================
c
c ======================================================================
c
c ... extrapola das pressoes de consolidacao para os nos
c
c ......................................................................
  800 continue
      inpi        = nint*nint*nint
      pci(1:inpi) = vplastic(3,1:inpi)
c ... pressao de consolidacao nodal total
      call extrapol_hexa20_v2(pci,pce,rn,sn,tn,1)
c .....................................................................
c
c ...
      p(1:8) = pce(1:8) 
      return
c ======================================================================
c
c ======================================================================
c
c ... calculo do parametro de encruamento                
c
c ......................................................................
  900 continue
      pc0  = e(12)
      mcs  = e(11)
      inpi = 0  
      do 910 lz = 1, nint
        ti = pg(lz,nint)
        do 920 ly = 1, nint
          si = pg(ly,nint)
          do 930 lx = 1, nint
            ri   = pg(lx,nint)
c ...
            inpi = inpi + 1
c .....................................................................
c
c ...
            call sfhexa8_m(hp,hpx,hpy,hpz,ri,si,ti,.true.,.false.)
c .....................................................................
c
c ... interpolacao da tensoes e pressoes nos pontos de integracoes
            txi(1:6)= 0.d0
            pi      = 0.d0
            do 940 j = 1,   8
              l      = j  + 60 
              txi(1) = txi(1) + hp(j)*tx0(1,j)
              txi(2) = txi(2) + hp(j)*tx0(2,j) 
              txi(3) = txi(3) + hp(j)*tx0(3,j) 
              txi(4) = txi(4) + hp(j)*tx0(4,j) 
              txi(5) = txi(5) + hp(j)*tx0(5,j)
              txi(6) = txi(6) + hp(j)*tx0(6,j) 
              pi     = pi     + hp(j)*u(l)
  940       continue 
c .....................................................................
c
c ... tensao efetiva
             txi(1:3) =  txi(1:3) + pi
c ..................................................................... 
c
c ...
             vplastic(3,inpi) = max(pc0,cam_clay_pc(txi,mcs))
c ..................................................................... 
c
  930     continue
  920   continue
  910 continue
c ......................................................................
      return
c ======================================================================
      end
c *********************************************************************
