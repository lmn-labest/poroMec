      subroutine calc_hidr(hi,hi1,dhi,tm,tm1,atil,ea,ma,nlit)
c **********************************************************************
c *                                                                    *
c *   calc_hidr:                                                       *
c *   -----                                                            *
c *                                                                    *
c *   Calculo valores de ksi(grau de hidratacao e passa para vetor b   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    atil - curva de afinidade                                       *
c *    tm   - temperatura no centro do elemento (em Kelvin)            *
c *    hi - grau de hidratacao                                         *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    hi1- grau de hidratacao                                         *
c *                                                                    *
c *                                                                    *
c **********************************************************************  
      implicit none
      include 'termprop.fi'
      include 'transiente.fi'
      integer n,i,maxit,ma,nlit,j,naf
      real*8 hi,hi1,dhi,dhi1,tm,atil,a,b,c1,c2,fa,fb,ea,r,beta1
      real*8 x,fx,eps,temp,temp1,tm1,r1,delta
      logical ok
c.......................................................................
      if (hi .gt. 0.98) then
        hi1 = 1.d0
        dhi = 0.d0
        return
      endif
c.......................................................................
c     parametros necessarios:
      naf = 1000
      maxit = 1000
      eps = 1.0e-3
      beta1 = 0.95
      temp = tm + 273.d0
      temp1 = tm1 + 273.d0
      n = eprop(9,ma)
      if (n .ne. 0 ) then
         a = min(hi+0.01,1.d0)
         b = min(a+0.01,1.d0)
         r = -ea/temp
         r1 = -ea/temp1
         c1 = exp(r1)*beta1*dt
         j =(int(hi*dfloat(naf-1))+1)
         atil = nprop(j,9,ma)
         c2 = atil*(1-beta1)*exp(r)*dt+hi
         j =(int(a*dfloat(naf-1))+1)
         atil = nprop(j,9,ma)
         fa = c2+c1*atil-a
         j =(int(b*dfloat(naf-1))+1)
         atil = nprop(j,9,ma)
         fb = c2 + c1*atil-b
c
c     busca tipo regula falsi:
c
         i = 0
         ok   = .false.
         do 100 while (.not. ok)
            i = i+1
            hi1 = (a*fb - b*fa)/(fb-fa)
            a    = b
            b    = hi1
            fa   = fb
            j =(int(b*dfloat(naf-1))+1)
            if((j.gt.0).and.(j.le.naf))then
               atil = nprop(j,9,ma)
               fb   = c2 + c1*atil - b
               if (dabs(fb).lt. eps) ok = .true.
               if (i    .eq.maxit) ok = .true.
            else
               if (hi .gt. 0.98) then
                  hi1 = 1.d0
                  dhi = 0.d0
                  return
               endif
               print*, '** erro em grauhid'
               print*, '  hyd1=',hi1
               print*, '  nit =',i
               stop
            endif
  100    continue 
         if (i .ge. maxit) then
            write (*,*) 'o método da falsa posicao nao atingiu 
     .      convergencia ate  a iteracao ', maxit, 'hi = ', hi1
            stop
         endif
         goto 1000
      endif
      return
 1000 continue
      hi1 = b
      j =(int(hi1*dfloat(naf-1))+1)
      atil = nprop(j,9,ma)
      dhi = atil*exp(r1)
      return
      end 
      subroutine check_var_hidr(a,b,c,d11,d12,d22,d33,ym,pr,alpha,hi,
     .                           ma,ndm,flaghidr,epd)
c **********************************************************************
c *                                                                    *
c *   CHECK_VAR_HIDR:                                                  *
c *   -----                                                            *
c *                                                                    *
c *   Verifica se as propriedades variam com o grau                    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    ma - material                                                   *
c *    hi - grau de hidratacao                                         *
c *    e - numero da propriedade                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    c - valor da propriedade                                        *
c *                                                                    *
c *                                                                    *
c **********************************************************************    
      implicit none
      integer ma,ndm
      real*8 ym,pr,alpha,hi,a,b,c,d11,d12,d22,d33,b1,b2,a1,a2,a3
      logical flaghidr,epd
      if (flaghidr)  then
         call updhidr(ym,hi,ma,1)
         call updhidr(pr,hi,ma,2)
         call updhidr(alpha,hi,ma,4)
         if( ndm .eq. 2) then
            if(epd)  then
               b1 = 1.d0-pr
               b2 = 1.d0-2.d0*pr
               a = ym*b1/((1.d0+pr)*b2)
               b = pr/b1
               c = b2/(2.d0*b1)
               d11 = a
               d12 = a*b
               d22 = a
               d33 = a*c  
            else
               a = ym/(1.d0-pr*pr)
               d11 = a
               d12 = a*pr
               d22 = a
               d33 = ym/(2.d0*(1.d0+pr))
            endif         
         elseif(ndm .eq. 3) then
            call updhidr(d11,hi,ma,11)
            call updhidr(d12,hi,ma,12)
            call updhidr(d22,hi,ma,13)
            a1 = 1.d0+pr
            a2 = 1.d0-pr
            a3 = 1.d0-2.d0*pr
            a = ym*a2/(a1*a3)
            b = pr/a2
            c = a3/(2.d0*a2)
         endif
      endif
      return
      end      
      subroutine updhidr(c,hi,ma,e)
c **********************************************************************
c *                                                                    *
c *   UPDHIDR:                                                         *
c *   -----                                                            *
c *                                                                    *
c *   Interpola valores com base em uma poligonal                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    ma - material                                                   *
c *    hi - grau de hidratacao                                         *
c *    e - numero da propriedade                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    c - valor da propriedade                                        *
c *                                                                    *
c *                                                                    *
c **********************************************************************    
      implicit none
      include 'termprop.fi'
      integer ma,e,n
      real*8 c,hi,a
      n = eprop(e,ma)
      if (n .ne. 0) then
          call interpol(nprop(1,e,ma),nprop(n+1,e,ma),hi,n,a)
          c = a * c
      endif
      return
      end
      subroutine teste_ksi0(hi,ix,nen,numel,nconc,flagm)
c **********************************************************************
c *                                                                    *
c *   UPDHIDR:                                                         *
c *   -----                                                            *
c *                                                                    *
c *   Interpola valores com base em uma poligonal                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    ma - material                                                   *
c *    hi - grau de hidratacao                                         *
c *    e - numero da propriedade                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    c - valor da propriedade                                        *
c *                                                                    *
c *                                                                    *
c **********************************************************************    
      implicit none
      integer i, ix(nen+1,*),nconc,nen,numel,ma,flagm,k
      real*8 hi(3,*),ksi0
      flagm = 1
      ksi0 = 0.1d0
      k = 0
      do 100 i = 1, numel
         ma = ix(nen+1,i)
         if (ma .eq. nconc)   then
            k = 1
            if(hi(2,i) .lt. ksi0)  flagm = 0
         endif
 100  continue
      if (k .eq. 0)  flagm = 1
      return
      end