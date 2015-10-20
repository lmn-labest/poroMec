c*****************************Svn***************************************      
c*$Date: 2010-09-06 18:52:24 -0300 (Mon, 06 Sep 2010) $                 
c*$Rev: 847 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************    
c ======================================================================
c      Preditor-multicorretor para derivada de ordem 1.                |
c ====================================================================== 
      subroutine predict(nnode,ndf,id,u,v)
c **********************************************************************
c *                                                                    *
c *   PREDICT                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 u(ndf,*),v(ndf,*),b
c ......................................................................      
      b = dt*(1.d0-alfa)
      do 110 i = 1, nnode
      do 100 j = 1, ndf
         k = id(j,i)
         if(k .gt. 0) then         
            u(j,i) = u(j,i) + v(j,i)*b
            v(j,i) = 0.d0
         endif
  100 continue
  110 continue
c ......................................................................        
      return
      end
      subroutine update(nnode,ndf,id,u,v,x)
c **********************************************************************
c *                                                                    *
c *   UPDATE                                                           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'      
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 u(ndf,*),v(ndf,*),x(*),b
c ......................................................................            
      b = alfa*dt
      do 110 i = 1, nnode
         do 100 j = 1, ndf
            k = id(j,i)
            if(k .gt. 0) then
                u(j,i) =  u(j,i) + x(k)*b
                v(j,i) =  v(j,i) + x(k)
            endif
  100    continue
  110 continue
c ......................................................................  
      return
      end
c ======================================================================
c      Preditor-multicorretor para derivada de ordem 2.                |
c ====================================================================== 
	subroutine predict2(nnode,ndf,id,u,v,a,dt,alfa,beta)
c **********************************************************************
c *                                                                    *
c *   PREDICT2                                                         *
c *                                                                    *
c **********************************************************************
	integer nnode,ndf,k
	integer id(ndf,*)
	real*8  u(ndf,*),v(ndf,*),a(ndf,*)
      real*8  dt,alfa,beta,c,d
c ......................................................................
      c = dt*dt*.5d0*(1.d0-2.d0*beta)
	d = dt*(1.d0-alfa)
	do 110 i = 1, nnode
	   do 100 j = 1, ndf
	      k = id(j,i)
	      if (k .gt. 0) then
	         u(j,i) = u(j,i) + dt*v(j,i) + c*a(j,i)
	         v(j,i) = v(j,i) + d*a(j,i)
               a(j,i) = 0.d0
            endif
  100    continue
  110 continue
      return
	end
	subroutine update2(nnode,ndf,id,u,v,a,x,dt,alfa,beta)
c **********************************************************************
c *                                                                    *
c *   UPDATE2                                                          *
c *                                                                    *
c **********************************************************************
	integer nnode,ndf,k
	integer id(ndf,*)
	real*8  u(ndf,*),v(ndf,*),a(ndf,*),x(*)
      real*8  dt,alfa,beta,c,d
c ......................................................................
      c = beta*dt*dt
	d = alfa*dt
	do 110 i = 1, nnode
	   do 100 j = 1, ndf
            k = id(j,i)
            if (k .gt. 0) then		
	         a(j,i) = a(j,i) + x(k)
	         u(j,i) = u(j,i) + c*x(k)
	         v(j,i) = v(j,i) + d*x(k)
	      endif
  100    continue
  110 continue
      return
	end
