
      subroutine gauleg(x,w)
      integer j
      REAL*8 x(16),w(16)

      x(1)= -0.9894009349916499d0
      x(2)= -0.9445750230732326d0
      x(3)= -0.8656312023878318d0
      x(4)= -0.7554044083550030d0
      x(5)= -0.6178762444026438d0
      x(6)= -0.4580167776572274d0
      x(7)= -0.2816035507792589d0
      x(8)= -9.501250983763744d-02

      w(1)= 2.715245941175185d-02
      w(2)= 6.225352393864778d-02
      w(3)= 9.515851168249290d-02
      w(4)= 0.1246289712555339d0
      w(5)= 0.1495959888165733d0
      w(6)= 0.1691565193950024d0
      w(7)= 0.1826034150449236d0
      w(8)= 0.1894506104550685d0

      do j=1,8
      x(17-j)=-x(j)
      w(17-j)=w(j)
      enddo

      end subroutine gauleg

!     c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE gaulegmaker(x,w,n)
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      x1=-1.d0
      x2=1.d0
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
         z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
 1       continue
      p1=1.d0
      p2=0.d0
      do 11 j=1,n
         p3=p2
         p2=p1
         p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
 11      continue
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
         
         z=z1-p1/pp
         if(abs(z-z1).gt.EPS)goto 1
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
 12   continue
!      write(6,*)x

      return
      END

      
!+++++++++++++++++++++++++++++++++++++++++++++++++++
      
      SUBROUTINE orderint(n,arr)
      INTEGER n,M,NSTACK
      integer arr(n), a,temp
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M)then
         
         do j=l+1,ir

            a=arr(j)
            do i=j-1,l,-1

               if(arr(i).le.a)goto 2
               arr(i+1)=arr(i)
            enddo
            i=l-1
 2          arr(i+1)=a

         enddo

         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      Else
         k=(l+ir)/2
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
         if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l).gt.arr(l+1))then
            temp=arr(l)
            arr(l)=arr(l+1)
            arr(l+1)=temp
         endif
         i=l+1
         j=ir
         a=arr(l+1)

 3       continue
         i=i+1
         if(arr(i).lt.a)goto 3
 4       continue
         j=j-1
         if(arr(j).gt.a)goto 4
         if(j.lt.i)goto 5
         temp=arr(i)
         arr(i)=arr(j)
         arr(j)=temp
         goto 3
 5       arr(l+1)=arr(j)
         arr(j)=a
         jstack=jstack+2

         if(jstack.gt.NSTACK)write(6,*)'NSTACK small in sort, stop'
         if(ir-i+1.ge.j-l)then

            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else

            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1

      END SUBROUTINE orderint
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=500)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL *8 a

      do j=1,n
         indx(j)=j
      enddo
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M)then
         do j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do i=j-1,l,-1
               if(arr(indx(i)).le.a)goto 2
               indx(i+1)=indx(i)
            enddo
            i=l-1
 2          indx(i+1)=indxt
         enddo

         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l)).gt.arr(indx(ir)))then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif

         if(arr(indx(l+1)).gt.arr(indx(ir)))then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif

         if(arr(indx(l)).gt.arr(indx(l+1)))then
            itemp=indx(l)
            indx(l)=indx(l+1)
            indx(l+1)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l+1)
         a=arr(indxt)
 3       continue
         i=i+1
         if(arr(indx(i)).lt.a)goto 3
 4       continue
         j=j-1
         if(arr(indx(j)).gt.a)goto 4
         if(j.lt.i)goto 5
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         goto 3
 5       indx(l+1)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
         if(jstack.gt.NSTACK)write(6,*)'NSTACK too small in indexx'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      END SUBROUTINE indexx

!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ---------------------------------------------------------------
      SUBROUTINE GCLOCK(XTIME)
      DOUBLE PRECISION XTIME
      REAL*4 TIME,TTIME
      DIMENSION TTIME(2)
!
!C     THIS ROUTINE IS MACHINE-DEPENDENT.
!C     IT SHOULD RETURN THE ELAPSED CPU TIME IN UNITS OF SECONDS.
!C     ONLY DIFFERENCES ARE USED, SO IT NEED NOT BE AN ABSOLUTE VALUE.
!C
!C     DUMMY RESULT FOR VANILLA DISTRIBUTION
      XTIME=0.D0
!C
!C!     CODE BELOW CALLS THE BSD UNIX TIMING ROUTINE.
!C     A C VERSION OF etime FOR MOST OTHER UNIX SYSTEMS IS AVAILABLE FROMJMH.
      TIME=etime(TTIME)
      XTIME=DBLE(TIME)
!C
!C     CODE BELOW IS THE GISS ROUTINE
!C     CALL CLOCKS(ITIME)
!C     XTIME=-ITIME
!C     XTIME=XTIME*1.D-2
      RETURN
      END SUBROUTINE GCLOCK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
      subroutine inv(a,o,n,np)
      implicit real*8 (a-h,o-z)
!     indx(nmax)
      dimension a(np,np),indx(n),o(np,np)
      call ludcmp(a,n,np,indx,d)
      do 12 i = 1,n
         do 11 j = 1,n
            o(i,j) = 0.0d0
 11      continue
         o(i,i) = 1.0d0
 12   continue
      do 13 j = 1,n
         call lubksb(a,n,np,indx,o(1,j))
 13   continue
      return
      end
!     c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (TINY=1.0D-20)
      DIMENSION A(NP,NP),INDX(N),VV(N)
      D=1.
      DO 12 I=1,N
         AAMAX=0.
         DO 11 J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
 11      CONTINUE
         IF (AAMAX.EQ.0.)then
            write(6,*)'Singular matrix.'
            return
         endif
         VV(I)=1./AAMAX
 12   CONTINUE
      DO 19 J=1,N
         DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
               DO 13 K=1,I-1
                  SUM=SUM-A(I,K)*A(K,J)
 13            CONTINUE
               A(I,J)=SUM
            ENDIF
 14      CONTINUE
         AAMAX=0.
         DO 16 I=J,N
            SUM=A(I,J)
            DO 15 K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
 15         CONTINUE
            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            ENDIF
 16      CONTINUE
         IF (J.NE.IMAX)THEN
            DO 17 K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
 17         CONTINUE
            D=-D
            VV(IMAX)=VV(J)
         ENDIF
         INDX(J)=IMAX
         IF(A(J,J).EQ.0.)A(J,J)=TINY
         IF(J.NE.N)THEN
            DUM=1./A(J,J)
            DO 18 I=J+1,N
               A(I,J)=A(I,J)*DUM
 18         CONTINUE
         ENDIF
 19   CONTINUE
      RETURN
      END
!     c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0)THEN
            DO 11 J=II,I-1
               SUM=SUM-A(I,J)*B(J)
 11         CONTINUE
         ELSE IF (SUM.NE.0.) THEN
            II=I
         ENDIF
         B(I)=SUM
 12   CONTINUE
      DO 14 I=N,1,-1
         SUM=B(I)
         IF(I.LT.N)THEN
            DO 13 J=I+1,N
               SUM=SUM-A(I,J)*B(J)
 13         CONTINUE
         ENDIF
         B(I)=SUM/A(I,I)
 14   CONTINUE
      RETURN
      END

      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK) write(6,*) 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
!     (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
      
      SUBROUTINE zbrak(fx,x1,x2,n,xb1,xb2,nb)
!     Given a function fx defined on the interval from x1-x2 subdivide the interval into n equally spaced segments,
!and search for zero crossings of the function. nb is input as the maxi- mum number of roots sought, and is reset to the number of bracketing pairs xb1(1:nb), xb2(1:nb) that are found.
      INTEGER n,nb
      REAL x1,x2,xb1(nb),xb2(nb),fx
      EXTERNAL fx
      INTEGER i,nbb
      REAL dx,fc,fp,x
      nbb=0
      x=x1
      dx=(x2-x1)/n
      fp=fx(x)
      do 11 i=1,n
        x=x+dx
        fc=fx(x)
        if(fc*fp.lt.0.) then
          nbb=nbb+1
          xb1(nbb)=x-dx
          xb2(nbb)=x
          if(nbb.eq.nb)goto 1
        endif
        fp=fc
11    continue
1     continue
      nb=nbb
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.

      SUBROUTINE zbrac(func,x1,x2,succes)
      INTEGER NTRY
      REAL x1,x2,func,FACTOR
      EXTERNAL func
      PARAMETER (FACTOR=1.4d0,NTRY=50)
      INTEGER j
      REAL f1,f2
      LOGICAL succes
      if(x1.eq.x2) then
         write(6,*) 'you have to guess an initial range in zbrac'
      endif
      f1=func(x1)
      f2=func(x2)
      succes=.true.
      do 11 j=1,NTRY
        if(f1*f2.lt.0.)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+FACTOR*(x1-x2)
          f1=func(x1)
        else
          x2=x2+FACTOR*(x2-x1)
          f2=func(x2)
        endif
11    continue
      succes=.false.
      return
      END

      SUBROUTINE MyZbrac(func,x1,x2,succes)
      INTEGER NTRY
      REAL x1,x2,func,FACTOR
      EXTERNAL func
      PARAMETER (FACTOR=1.5d0,NTRY=50)
      INTEGER j
      REAL f1,f2
      LOGICAL succes
      if(x1.eq.x2) then
         write(6,*) 'you have to guess an initial range in zbrac'
      endif
      f1=func(x1)
      f2=func(x2)
      succes=.true.
      do 11 j=1,NTRY
        if(f1*f2.lt.0.)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+FACTOR*(x1-x2)
          f1=func(x1)
        else
          x2=x2+FACTOR*(x2-x1)
          f2=func(x2)
        endif
11    continue
      succes=.false.
      return
      END
C     (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
      
      FUNCTION rtbis(func,x1,x2,xacc)
      INTEGER JMAX
      REAL*8 rtbis,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      REAL*8 dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.) write(6,*) 'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      write(6,*) 'too many bisections in rtbis'
      END
C     (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
      
      FUNCTION zriddr(func,x1,x2,xacc)
      INTEGER MAXIT
      REAL zriddr,x1,x2,xacc,func,UNUSED
      PARAMETER (MAXIT=60,UNUSED=-1.11E30)
      EXTERNAL func
CU    USES func
      INTEGER j
      REAL fh,fl,fm,fnew,s,xh,xl,xm,xnew
      fl=func(x1)
      fh=func(x2)
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
        xl=x1
        xh=x2
        zriddr=UNUSED
        do 11 j=1,MAXIT
          xm=0.5*(xl+xh)
          fm=func(xm)
          s=sqrt(fm**2-fl*fh)
          if(s.eq.0.)return
          xnew=xm+(xm-xl)*(sign(1.,fl-fh)*fm/s)
          if (abs(xnew-zriddr).le.xacc) return
          zriddr=xnew
          fnew=func(zriddr)
          if (fnew.eq.0.) return
          if(sign(fm,fnew).ne.fm) then
            xl=xm
            fl=fm
            xh=zriddr
            fh=fnew
          else if(sign(fl,fnew).ne.fl) then
            xh=zriddr
            fh=fnew
          else if(sign(fh,fnew).ne.fh) then
            xl=zriddr
            fl=fnew
          else
            write(6,*) 'never get here in zriddr'
          endif
          if(abs(xh-xl).le.xacc) return
11      continue
        write(6,*) 'zriddr exceed maximum iterations'
      else if (fl.eq.0.) then
        zriddr=x1
      else if (fh.eq.0.) then
        zriddr=x2
      else
        write(6,*) 'root must be bracketed in zriddr'
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
