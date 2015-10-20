c*****************************Svn***************************************      
c*$Date: 2013-06-21 12:34:21 -0300 (Fri, 21 Jun 2013) $                 
c*$Rev: 969 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************          
      subroutine pload(x,id,f,u,v,b,nload,nnode,ndf)
c **********************************************************************
c *                                                                    *
c *   PLOAD: Monta o vetor de forcas                                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    id(ndf,nnode)   - numeracao das equacoes                        *
c *    f(ndf,nnode)    - cargas concentradas e incognitas prescritas   *
c *    b(neq) - nao definido                                           *
c *    nload(ndf,nnode) - numero da carga nodal                        *
c *    nnode  - numero de nos acessados na particao                    *
c *    ndf    - numero de graus de liberdade por no                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    b - vetor de forcas                                             *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      integer id(ndf,*),nload(ndf,*),nnode,ndf,i,j,k,nc
      real*8  x(3,*),f(ndf,*),u(ndf,*),v(ndf,*),b(*),c,a
c.......................................................................
c
c ... Cargas nodais e desloc. prescritos variaveis no tempo:
c
      do 110 i = 1, nnode
	   do 100 j = 1, ndf
	      nc = nload(j,i)
	      if(nc .gt. 0) then
	         call tload(nc,t,u(j,i),c)
c ===========================================
c               if (j .eq. 3) then
c                  gH=(x(3,i)+c)*9.81
c                  f(ndf,i) = 2.d0*dsqrt(gH)
c               endif
c ===========================================
               f(j,i) = c
	      endif
  100    continue
  110 continue
c.......................................................................  
c
c ... Cargas nodais e deslocamentos prescritos no tempo t:
c
      a = alfa*dt
      do 200 i = 1, nnode
	do 200 j = 1, ndf
	   k = id(j,i)
         if (k .gt. 0) then
            b(k) = f(j,i)
         else
c            v(j,i) = (f(j,i)-u(j,i))/a - v(j,i)*((1.d0-alfa)/alfa)
c            v(j,i) = 0.d0            
            v(j,i) = (f(j,i)-u(j,i))/dt
            u(j,i) =  f(j,i)
         endif
  200 continue
c.......................................................................  
      return
      end
      subroutine tload(nc,t,u,c)
c **********************************************************************
c *                                                                    *
c *   TLOAD:                                                           *
c *   -----                                                            *
c *                                                                    *
c *   Funcoes de cargas nodais                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    load(2,numload)  - load(1,j) = tipo da carga j                  *
c *                     - load(2,j) = numero de parcelas da carga j    *
c *    fload(numpar,numsum,numload) - definicao das cargas             *
c *                 fload(i,j,k) = parametro i da parcela j da carga j *
c *    nterms  - numero maximo de parcelas                             *
c *    npar    - numero maximo de parametros                           *
c *    nc      - numero da carga                                       *
c *    t       - tempo                                                 *
c *    c       - nao definido                                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    c - valor da carga nodal no tempo t                             *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'load.fi'
      integer nc,itype,nparc,i
      real*8  t,u,c,a,w,t0,ti,DOISPI,tm,uext,b,h
      parameter (DOISPI = 6.283185307d0)
c.......................................................................
c      
c ... Cargas nodais variaveis no tempo:
c
	c  = 0.d0
	if(nc .le. 0) return
      itype = load(1,nc)
	nparc = load(2,nc)
      if(itype .eq. 1) then
c
c ...    Carga constante no tempo:
c
         c = fload(1,1,nc)
c         if(t .gt. 50.d0) then
c            c = 0.d0
c         endif
      elseif(itype .eq. 2) then
c
c ...    Variacao linear:            
c
         c = fload(1,1,nc) + fload(2,1,nc)*t
c ....   "gancho":
c         if (t .gt. 28800.0) then
c            if (t .le. 43200.0) then
c               c = -0.044444d0*t+1950.d0
c            else
c               c = 30.d0
c            endif
c        endif
c        c = 40
      elseif(itype .eq. 3) then
c
c ...    troca de calor de Newton:            
c
         c = (fload(1,1,nc)-u)*fload(2,1,nc)
	elseif(itype .eq. 4) then
c
c ...    Variacao senoidal:            
c
	   do 200 i = 1, nparc
c ...       Amplitude:
	      a = fload(1,i,nc)
c ...       Frequencia:
		   w = 1.d0/fload(2,i,nc)
c ...       Fase:
		   t0 = fload(3,i,nc)
	      ti = max(t-t0,0.d0)
	      c  = c + a*dsin(DOISPI*w*ti) 
  200    continue
	elseif(itype .eq. 5) then
c ...    Variacao cosssenoidal:            
	   do 300 i = 1, nparc
c ...       Amplitude:
	      a = fload(1,i,nc)
c ...       Frequencia:
		   w = 1.d0/fload(2,i,nc)
c ...       Fase:
		   t0 = fload(3,i,nc)
	      ti = max(t-t0,0.d0)
	      c  = c + a*dcos(DOISPI*w*ti)
  300    continue
	elseif(itype .eq. 6) then
	   do 400 i = 1, nparc	  
	  	   a  = fload(1,i,nc)
c           periodo:
		   w  = 1.d0/fload(2,i,nc)
c           frequencia:
c		    w  = fload(2,i,nc)
		   t0 = fload(3,i,nc)
c	      ti = max(t-t0,0.d0)
c           ti = t + t0
	      c  = c + a*dcos(DOISPI*w*t+t0)
  400    continue
	elseif(itype .eq. 7) then
c ...       kdu/dx = emiss * const(Stef-Boltz) *(uext4-u4)+H(uext-u)
c         call interpol(fload(1,2,nc),fload(1,3,nc),t,nparc,uext)
           call interpol(fload(1,2,nc),fload(1,3,nc),t,nparc,c)
c         a = fload(1,1,nc)
c         b = fload(2,1,nc)
c         h = fload(3,1,nc)
c         c = a*b*(uext**4-u**4)+h*(uext-u)
 500	elseif(itype .eq. 8) then
c ...       kdu/dx = emiss * const(Stef-Boltz) *(uext4-u4)+H(uext-u)
c         call interpol(fload(1,2,nc),fload(1,3,nc),t,nparc,uext)
c           if( t .ge. 172800.d0)   then
              call interpol(fload(1,1,nc),fload(1,2,nc),t,nparc,c)
c           else  
c              c = 0.d0
c           endif
c         a = fload(1,1,nc)
c         b = fload(2,1,nc)
c         h = fload(3,1,nc)
c         c = a*b*(uext**4-u**4)+h*(uext-u)
	endif
      return
      end
      subroutine global(id,f,x,u,nnode,ndf)
c **********************************************************************
c *                                                                    *
c *   Subroutine GLOBAL                                                *
c *   -----------------                                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8  f(ndf,*),x(*),u(ndf,*)
c ......................................................................
      do 100 i = 1, nnode
      do 100 j = 1, ndf
         k = id(j,i)
         if (k .gt. 0) then
            u(j,i) = x(k)
         else
            u(j,i) = f(j,i)
         endif
  100 continue
      return
      end
      subroutine updgeo(id,f,v,u,x,nnode,ndm,ndf)
c **********************************************************************
c *                                                                    *
c *   Subroutine UPDGEO                                                *
c *   -----------------                                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndm,ndf,id(ndf,*),i,j,k
      real*8  x(ndm,*),u(ndf,*),v(*),f(ndf,*)
c ......................................................................
      do 100 i = 1, nnode
      do 100 j = 1, ndm
         k = id(j,i)
         if (k .gt. 0) then
            x(j,i) = x(j,i) + v(k)
         else
            x(j,i) = x(j,i) + f(j,i) - u(j,i)
         endif
  100 continue
      return
      end
