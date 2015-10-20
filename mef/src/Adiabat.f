c **********************************************************************
c *                                                                    *
c *                          AFINIDADE                                 *
c *                          ---------                                 *
c *                                                                    *
c *   Determina curva de afinidade normalizada a partir de curva       *
c *   de temperatura adiabatica                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   dt          -   amplitude dos intervalos de tempo                *
c *   afin(npt)   -   lista de temperaturas adiabaticas (C)            *
c *   npt         -   numero de pontos                                 *
c *   Ea          -   energia de ativacao                              *
c *                                                                    *
c *   Parametro de saida:                                              *
c *   ------------------                                               *
c *                                                                    *
c *   hyd (npt)   -   grau de hidratacao                               *
c *   afin(npt)   -   afinidade normalizada                            *
c *                                                                    *
c **********************************************************************
      subroutine afinidade(t,Tad,ea,ma,npt)
      implicit none
      include 'termprop.fi'
      integer i,j,k,numat,npt,ma
      real*8  t(*),Tad(*)
      real*8  Ea,DTad,Tad0
      real*8  x1,x2,x3,y1,y2,y3,a,b
c.......................................................................
c
c     constantes
c
      Tad0 = Tad(1)
      DTad = Tad(npt)-Tad0
c.......................................................................
c
c     primeiro ponto
c
      x1 = t(1)
      x2 = t(2)
      x3 = t(3)
      y1 = Tad(1)
      y2 = Tad(2)
      y3 = Tad(3)
c
c     parábola (Tad x t)
c
      a = ((y1-y3)*(x2-x3)-(y2-y3)*(x1-x3))/(-(x2-x1)*(x1-x3)*(x2-x3))
      b = (y2-y3)/(x2-x3)-a*(x2+x3)
c      c = y3-x3*(a*x3+b)
c
c     derivada da parábola d(Tad)/dt em t=x1
c
      y3 = 2*a*x1 + b
c
c     afinidade
c
c      nprop(npt+1,9,ma) = exp( Ea/y1 ) * y3 / DTad
      nprop(npt+1,9,ma) = 0.d0
c
c     grau de hidratacao
c
      nprop(1,9,ma) = 0.d0
c.......................................................................
c
c     loop nos pontos intermediários
c
      do i=2,npt-1
c
c       pontos auxiliares
c
        x1 = t(i-1)
        x2 = t(i)
        x3 = t(i+1)
        y1 = Tad(i-1)
        y2 = Tad(i)
        y3 = Tad(i+1)
c
c       parábola (Tad x t)
c
        a = ((y1-y3)*(x2-x3)-(y2-y3)*(x1-x3))/((x2-x3)*(x1-x3)*(x2-x3))
        b = (y2-y3)/(x2-x3) - a*(x2+x3)
c        c = y3-x3*(a*x3+b)
c
c       derivada da parábola d(Tad)/dt em t=x2
c
        y3 = 2*a*x2 + b
c
c       afinidade
c
        nprop(npt+i,9,ma) = exp( Ea/y2 ) * y3 / DTad
c
c       grau de hidratacao
c
        nprop(i,9,ma) = (y2-Tad0)/DTad
        
c        write(125,*) hyd(i) , afin(i)
      enddo
c.......................................................................
c
c     ultimo ponto
c
      nprop(npt,9,ma) = 1.d0
      nprop(npt+npt,9,ma) = 0.d0
c.......................................................................
      Ea = DTad
c.......................................................................
      return
      end
c **********************************************************************
c *                                                                    *
c *   UNIFORMIZA:                                                      *
c *   ----------                                                       *
c *                                                                    *
c *   uniformiza uma curva dada por pontos usando interpolacao         *
c *   quadratica e intervalos uniformes na abcissa                     *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   x0(n1)       -  lista de valores na abcissa  [0,n1]              *
c *   y0(n1)       -  lista de valores na ordenada [0,n1]              *
c *   n1           -  numero de pontos original                        *
c *   n2           -  numero de pontos final                           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   yu(n2)       -  lista de valores na ordenada [0,n2]              *
c *                                                                    *
c **********************************************************************
      subroutine uniformiza(ma,n1,n2)
      implicit none
      include 'termprop.fi'
      integer n1,n2,i,j,k,m,nint,ma
      real yu(1000)
      real*8  x,ampl,x1,x2,x3,y1,y2,y3,a,b,c,dx,x23,x13,y23
c.......................................................................
c
c     k = numeracao da sequencia uniforme
      k = 1
c
c     m = numeracao da sequencia original
      m = 1
c
c     x = abscissa uniformizada atual
      x = nprop(1,9,ma)
c
c     yu(*) = lista de ordenadas uniformizadas
      yu(1) = nprop(n1+1,9,ma)
c
c     ampl = amplitude uniforme de intervalo
      ampl = (nprop(n1,9,ma)-x)/(n2-1)
c
c     loop nos intervalos originais
      do while (m .le. n1-2)
         m  = m + 1
         x1 = x
         y1 = yu(k)
         x2 = nprop(m,9,ma)
         y2 = nprop(n1+m,9,ma)
         dx = x2-x1
         if (dx .gt. ampl) then
            x3 = nprop(m+1,9,ma)
            y3 = nprop(n1+m+1,9,ma)
            x23 = x2-x3
            x13 = x1-x3
            y23 = y2-y3
            a = ((y1-y3)*x23-y23*x13)/(-dx*x13*x23)
            b = y23/x23-a*(x2+x3)
            c = y3-x3*(a*x3+b)
            nint = dx/ampl
            do j = 1,nint
               x = x + ampl
               k = k + 1
               yu(k) = x*(a*x + b) + c
            enddo
         else
            do while (x2 .lt. x+ampl)
               m = m + 1
               x2 = nprop(m,9,ma)
            enddo
            y2 = nprop(n1+m,9,ma)
            a = (y2-y1)/(x2-x1)
            b = y1 - a*x1
            x = x + ampl
            k = k + 1
            yu(k) = a*x + b
            if (x2-x .gt. ampl) m = m - 1
         endif
      enddo
c.......................................................................
c
c     Ultimo intervalo
c
      if (k .lt. n2) then
         x1 = nprop(n1-2,9,ma)
         y1 = nprop(n1+n1-2,9,ma)
         x2 = x
         y2 = yu(k)
         x3 = nprop(n1,9,ma)
         y3 = nprop(n1+n1,9,ma)
         dx = x3-x2
         if (dx .gt. ampl) then
            x23 = x2-x3
            x13 = x1-x3
            y23 = y2-y3
            a = ((y1-y3)*x23-y23*x13)/(-dx*x13*x23)
            b = y23/x23-a*(x2+x3)
            c = y3-x3*(a*x3+b)
            nint = dx/ampl
            do j = 1,nint
               x = x + ampl
               k = k + 1
               yu(k) = x*(a*x + b) + c
            enddo
         else
            k = k + 1
            yu(k) = y3
         endif
         if (k.lt.n2) then
            yu(n2) = y3
         endif
      endif
      do i = 1, n2
         nprop(i,9,ma) = yu(i)
      enddo
c.......................................................................
      return
      end