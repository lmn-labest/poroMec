c*****************************Svn***************************************      
c*$Date: 2011-08-08 15:20:15 -0300 (Mon, 08 Aug 2011) $                 
c*$Rev: 939 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************       
      subroutine profil(ix,id,jd,nnode,numel,nen,ndf,neq,nad,ptr)
c **********************************************************************
c *                                                                    *
c *   Subroutine PROFIL                                                *
c *                                                                    *
c *   Computes equation profile                                        *
c *                                                                    *
c *   INPUT :                                                          *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    id(ndf,nnode)   - nodal equation numbering                      *
c *    jd(neq)- not defined                                            *
c *    nnode - number of nodes                                         *
c *    numel - number of elements                                      *
c *    ndf   - number of degrees of freedom per node                   *
c *    neq   - number of equations                                     *
c *    nad   - not defined                                             *
c *                                                                    *
c *   OUTPUT :                                                         *
c *                                                                    *
c *    jd    - profile pointer                                         *
c *    nad   - number of terms in profile                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer ix(nen+1,*),id(ndf,*),jd(*),nnode,numel,nen,ndf,neq,nad
      integer i,n,ii,jj,mm,k,no
      logical ptr
c ......................................................................
      do 100 i = 1, neq
        jd(i) = 0
  100 continue
c ......................................................................  
c
c.... compute column heights
c
      do 220 n = 1, numel
         mm = 0
         do 200 i = 1, nen
            no = ix(i,n)
            if(no .eq. 0) goto 200
            do k = 1, ndf
               jj = id(k,no)
               if (jj .gt. 0) then
                  if(mm .eq. 0) mm = jj
                  mm = min(mm,jj)
               endif 
            enddo
  200    continue
         do 210 i = 1, nen
            no = ix(i,n)
            if(no .eq. 0) goto 210
            do k = 1, ndf
               ii = id(k,no)
                   if(ii .gt. 0) then
                      jj = jd(ii)
                      jd(ii) = max(jj,ii-mm)
                   endif 
            enddo
  210    continue
  220 continue
c ......................................................................  
c
c.... compute diagonal pointers for profile
c
      nad   = 0
      jd(1) = 0
      if (neq .gt. 1) then
        do 300 n = 2, neq
          jd(n) = jd(n) + jd(n-1)
  300   continue
        nad = jd(neq)
      endif
c
c.... equation summary
c
      mm = 0
      if (neq .gt. 0) mm = (nad+neq)/neq
      if (ptr) write(*,2001) neq,nnode,mm,numel,nad,(nad+neq)*8+neq*4,
     .                      (nad+nad+neq)*8+neq*4
c ......................................................................     
      return
 2001 format(/' (PROFIL) Equation Problem Summary:'//
     1 5x,'Number of equations  =',i9,5x,'Number nodes      =',i5/
     2 5x,'Average col. height  =',i9,5x,'Number elements   =',i5/
     3 5x,'No. terms in profile =',i9/
     4 5x,'Memory required  =  ',i11,' bytes (symmetric)'/26x,i10,
     5 ' bytes (unsymmetric)'/)
      end
      subroutine dtri(ad,au,al,jp,neq,flg)
c **********************************************************************
c *                                                                    *
c *   Subroutine DTRI                                                  *
c *   ---------------                                                  *
c *                                                                    *
c *   Decomposition of coefficient matrix into its triangular          *
c *   factors (LDU).                                                   *
c *                                                                    *
c *                                                                    *
c *   INPUT :                                                          *
c *   -------                                                          *
c *                                                                    *
c *    jp   - pointer to last element in each column/row of au/al      *
c *           respectively                                             *
c *           dimension: jp(neq)                                       *
c *    ad   - diagonal coefficients                                    *
c *           dimension: ad(neq)                                       *
c *    au   - coefficients above diagonal                              *
c *           dimension: au(nad), nad = jp(neq)                        *
c *    al   - coefficients below diagonal                              *
c *           dimension: al(nad), nad = jp(neq)                        *
c *           if flg = .true., al = au                                 *
c *    neq  - number of equations                                      *
c *    flg  - flag                                                     *
c *           if true equations are treated as unsymmetric and         *
c *           separate storage must be provided for au and al          *
c *                                                                    *
c *   OUTPUT :                                                         *
c *   -------                                                          *
c *                                                                    *
c *    jp   - unchanged                                                *
c *    ad   - factor D                                                 *
c *    au   - factor U                                                 *
c *    al   - factor L                                                 *
c *    neq  - unchanged                                                *
c *    flg  - unchanged                                                *
c *                                                                    *
c **********************************************************************
      implicit none
      integer jp(*),neq,i,id,ie,ih,is,idh,ifig,j,jd,jh,jr,jrh
      real*8 al(*),au(*),ad(*),dot,zero,one,tol,dimx,dimn,dfig,daval,dd
      logical flg
c ......................................................................
c.... n.b.  tol should be set to approximate half-word precision.
      data zero,one/0.0d0,1.0d0/, tol/0.5d-07/
c
c.... set initial values for conditioning check
c
      dimx = zero
      dimn = zero
      do 50 j = 1,neq
      dimn = max(dimn,abs(ad(j)))
50    continue
      dfig = zero
      dd   = zero
c
c.... loop through the columns to perform the triangular decomposition
c
      jd = 1
      do 200 j = 1,neq
        jr = jd + 1
        jd = jp(j)
        jh = jd - jr
        if(jh.gt.0) then
          is = j - jh
          ie = j - 1
c
c.... if diagonal is zero compute a norm for singularity test
c
          if(ad(j).eq.zero) call dtest(au(jr),jh,daval)
          do 100 i = is,ie
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              jrh = jr - ih
              idh = id - ih + 1
              au(jr) = au(jr) - dot(au(jrh),al(idh),ih)
              if(flg) al(jr) = al(jr) - dot(al(jrh),au(idh),ih)
            endif
100       continue
        endif
c
c.... reduce the diagonal
c
        if(jh.ge.0) then
          dd = ad(j)
          jr = jd - jh
          jrh = j - jh - 1
          call dredu(al(jr),au(jr),ad(jrh),jh+1,flg,ad(j))
c
c.... check for possible errors and print warnings
c
          if(abs(ad(j)).lt.tol*abs(dd))    write(*,2000) j
          if(dd.lt.zero.and.ad(j).gt.zero) write(*,2001) j
          if(dd.gt.zero.and.ad(j).lt.zero) write(*,2001) j
          if(ad(j) .eq.  zero)             write(*,2002) j
          if(dd.eq.zero.and.jh.gt.0) then
            if(abs(ad(j)) .lt. tol*daval)  write(*,2003) j
          endif
        endif
c
c.... store reciprocal of diagonal, compute condition checks
c
        if(ad(j).ne.zero) then
          dimx  = max(dimx,abs(ad(j)))
          dimn  = min(dimn,abs(ad(j)))
          dfig  = max(dfig,abs(dd/ad(j)))
          ad(j) = one/ad(j)
        endif
200   continue
c
c.... print conditioning information
c
      dd = zero
      if(dimn.ne.zero) dd = dimx/dimn
      ifig = dlog10(dfig) + 0.6
      write(*,2004) dimx,dimn,dd,ifig
      return
 2000 format(' ***DTRI WARNING 1*** Loss of at least 7 digits in',
     1 ' reducing diagonal of equation',i5)
 2001 format(' ***DTRI WARNING 2*** Sign of diagonal changed when',
     1 ' reducing equation',i5)
 2002 format(' ***DTRI WARNING 3*** Reduced diagonal is zero for',
     1 ' equation',i5)
 2003 format(' ***DTRI WARNING 4*** Rank failure for zero unreduced',
     1 ' diagonal in equation',i5)
 2004 format(/,' (DTRI) Condition check:',/,
     .  ' D-max',e11.4,'; D-min',e11.4,
     1 '; Ratio',e11.4/' Maximum no. diagonal digits lost:',i3,/)
 2005 format('Cond ck: Dmax',1p1e9.2,'; Dmin',1p1e9.2,'; Ratio',1p1e9.2)
      end
      subroutine dredu(al,au,ad,jh,flg,dj)
c **********************************************************************
c *                                                                    *
c *   Subroutine DREDU                                                 *
c *   ----------------                                                 *
c *                                                                    *
c *   Reduce diagonal element in triangular decomposition              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer jh,j
      real*8 al(jh),au(jh),ad(jh),ud,dj
      logical flg
c ......................................................................
      do 100 j = 1,jh
        ud = au(j)*ad(j)
        dj = dj - al(j)*ud
        au(j) = ud
100   continue
c
c.... finish computation of column of al for unsymmetric matrices
c
      if(flg) then
        do 200 j = 1,jh
          al(j) = al(j)*ad(j)
200     continue
      endif
      return
      end
      subroutine dtest(au,jh,daval)
c **********************************************************************
c *                                                                    *
c *   Subroutine DTEST                                                 *
c *   ----------------                                                 *
c *                                                                    *
c *   Test for rank                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer jh,j
      real*8 au(jh),daval
c ......................................................................
      daval = 0.0d0
      do 100 j = 1,jh
         daval = daval + abs(au(j))
100   continue
      return
      end
      subroutine dsolv(ad,au,al,b,jp,neq,energy,ptr)
c **********************************************************************
c *                                                                    *
c *   Subroutine DSOLV                                                 *
c *   ----------------                                                 *
c *                                                                    *
c *   Direct solver for symmetric or unsymmetric systems of            *
c *   equations (LDU decomposition).                                   *
c *   Coefficient matrix must be decomposed into its triangular        *
c *   factors using subroutine DTRI.                                   *
c *                                                                    *
c *                                                                    *
c *   INPUT :                                                          *
c *   -------                                                          *
c *                                                                    *
c *    jp   - pointer to last element in each column/row of au/al      *
c *           respectively                                             *
c *           dimension: jp(neq)                                       *
c *    ad   - diagonal elements (factor D)                             *
c *           dimension: ad(neq)                                       *
c *    au   - Factor matrix U                                          *
c *           dimension: au(nad), nad = jp(neq)                        *
c *    al   - Factor matrix L                                          *
c *           dimension: al(nad), nad = jp(neq)                        *
c *    b    - right side of equations                                  *
c *           dimension: b(neq)                                        *
c *    neq  - number of equations                                      *
c *    energy - not defined                                            *
c *    ptr    - if .true. -> print energy                              *
c *                                                                    *
c *                                                                    *
c *   OUTPUT :                                                         *
c *   -------                                                          *
c *                                                                    *
c *    jp   - unchanged                                                *
c *    ad   - unchanged                                                *
c *    au   - unchanged                                                *
c *    al   - unchanged                                                *
c *    b    - solution                                                 *
c *    neq  - unchanged                                                *
c *    energy - energy norm of solution                                *
c *                                                                    *
c **********************************************************************
      implicit none
      logical ptr
      integer jp(*),neq,is,j,jh,jr
      real*8  al(*),au(*),ad(*),b(*),dot,energy,zero,bd
      data zero/1.d-15/
c ......................................................................
c
c.... find the first nonzero entry in the right hand side
c
      do 100 is = 1, neq
        if (dabs(b(is)) .gt. zero) go to 200
100   continue
      write(*,2000)
      return
200   if (is .lt. neq) then
c
c.... reduce the right hand side
c
        do 300 j = is+1, neq
          jr = jp(j-1)
          jh = jp(j) - jr
          if (jh .gt. 0) then
            b(j) = b(j) - dot(al(jr+1),b(j-jh),jh)
          endif
300     continue
      endif
c
c.... multiply by inverse of diagonal elements
c
      energy = zero
      do 400 j = is, neq
        bd = b(j)
        b(j) = b(j)*ad(j)
        energy = energy + bd*b(j)
400   continue
c
c.... backsubstitution
c
      if (neq .gt. 1) then
        do 500 j = neq, 2, -1
          jr = jp(j-1)
          jh = jp(j) - jr
          if (jh .gt. 0) then
            call saxpb(au(jr+1),b(j-jh),-b(j),jh, b(j-jh))
          endif
500     continue
      endif
      if (ptr) write(*,2100) energy
      return
 2000 format(' ***DSOLV WARNING 1*** Zero right-hand-side vector')
 2100 format(1x,'(DSOLV) Energy norm = ',d13.6,/)
      end
