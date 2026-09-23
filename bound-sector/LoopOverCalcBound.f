c THIS PROGRAM CALCULATES THE BOUND STATES SUPPORTED BY A MULTICHANNEL SET OF 
C POTENTIALS READ IN FROM FILES FORT.100.  THE (ADIABATIC) COUPLINGS ARE READ IN FROM FORT.101, FORT.102, FORT.103
c NB: Because this is a bound-state code, it does NOT do R-matrix poropagation and therefore
c     we must alway shave iSector=1.
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      program LoopOverCalcBound
      implicit none
      call CalcBound
      end

      subroutine CalcBound
      use bspline
      implicit none
      integer LegPoints,RNumPoints,NumChan,PotInterpPoints,NumTot,PotInterpOrder ! 
      integer Order, RDimMin,Num2BodyChan,NumSectors,n
      integer Left,Right,RDim
      double precision, allocatable :: R(:), SectorR(:)
      integer i,j,k,iR,Left0,iSector,iEnergy
      double precision RMin, RMax,alpha,delR,mass,EnergyStart,EnergyEnd,deltaEnergy ! 
      double precision mu, mu12, Pi,energy,x,RMaxRun
      double precision, allocatable :: RMatNew(:,:), RMatOld(:,:)
      double precision, allocatable :: xLeg(:),wLeg(:)
      double precision, allocatable :: KK(:,:),kchan(:),BigPot(:,:,:,:,:),BigAPot(:,:,:,:,:) ! 
      double complex, allocatable ::  SS(:,:)
      double complex, allocatable :: tmp1(:,:), tmp2(:,:)
      double precision, allocatable :: Ucurves(:,:), Potknots(:), Pbcoef(:,:,:),Qbcoef(:,:,:),Ebcoef(:,:) ! 
      double precision, allocatable :: P(:,:,:),QTemp(:,:,:),Q(:,:,:),dP(:,:,:) ! 
      character*64 LegendreFile
      double precision hbarc
      double complex II
      parameter(II=(0.0d0,1.0d0))
      common hbarc

      PotInterpOrder=3
c      hbarc=197.32858d0
      hbarc=1.0d0

c     read in number of channels and order of splines
      print*, '#NumTot, NumChan, Order'
      read(5,*)
      read(5,*) NumTot, NumChan, order
      write(6,*) NumTot, NumChan, order
      print*, '#Legendre File'
      read(5,*)
      read(5,*)
      read(5,1002) LegendreFile
      write(6,1002) LegendreFile

      print*, '#alpha,  LegPoints'
      read(5,*)
      read(5,*)
      read(5,*) alpha, LegPoints
      print*, alpha, LegPoints
      Pi=dacos(-1.d0)
      write(6,*) '#Pi=',Pi

c     read in domain information
      read(5,*)
      read(5,*)
      read(5,*) PotInterpPoints,RMin,RMax, mass
      write(6,*) PotInterpPoints,RMin,RMax, mass
      read(5,*)
      read(5,*)
      read(5,*) RNumPoints, NumSectors ! 
      write(6,*) 'RnumPoints = ',RNumPoints,'NumSectors=',NumSectors ! 

      mu=mass/dsqrt(3.0d0)

      print*, '# Reading in Gauss-Legendre points and weights'
      allocate(xLeg(LegPoints),wLeg(LegPoints))
      call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)
      
      allocate(R(PotInterpPoints))
      
      allocate(Ucurves(NumChan,PotInterpPoints))
      allocate(P(NumChan,NumChan,PotInterpPoints),
     >     Q(NumChan,NumChan,PotInterpPoints), ! 
     >     QTemp(NumChan,NumChan,PotInterpPoints),
     >     dP(NumChan,NumChan,PotInterpPoints)) ! 

      write(6,*), '#Reading in the adiabatic potentials and couplings'

      call ReadCurves(PotInterpPoints,R,Ucurves,P,QTemp,Q,NumChan,NumTot) ! 

      allocate(Potknots(PotInterpPoints+PotInterpOrder))
      allocate(Pbcoef(NumChan,NumChan,PotInterpPoints),Qbcoef(NumChan,NumChan,PotInterpPoints)) ! 
      allocate(Ebcoef(NumChan,PotInterpPoints))
      Num2BodyChan=0
      
      do n=1,NumChan
         if(Ebcoef(n,PotInterpPoints).lt.0.0d0) Num2BodyChan=Num2BodyChan+1 ! 
      enddo

      call setup_interp(R,PotInterpOrder,Potknots,PotInterpPoints,P,Q,Ucurves,Ebcoef,Pbcoef,Qbcoef,NumChan) ! 
c      do k=1,NumChan
c      x=R(1)
c         do while(x.lt.RMax)
         do i=1,PotInterpPoints
c            write(20,*) R(i), Ucurves(1,i) !dbsval(R(i),PotInterpOrder,Potknots,PotInterpPoints,Ebcoef(1,:)) ! 
c            write(20,*) x, dbsval(x,PotInterpOrder,Potknots,PotInterpPoints,Ebcoef(k,:))! 
c            write(20,*) x, dbsval(x,PotInterpOrder,Potknots,PotInterpPoints,Qbcoef(1,1,:)), ! 
c     >           dbsval(x,PotInterpOrder,Potknots,PotInterpPoints,Qbcoef(1,2,:)) !  ! 
c            write(20,*) R(i), dbsval(R(i),PotInterpOrder,Potknots,PotInterpPoints,Qbcoef(1,1,:)), ! 
c     >           dbsval(R(i),PotInterpOrder,Potknots,PotInterpPoints,Qbcoef(1,2,:)) !  ! 
            write(20,*) R(i), dbsval(R(i),PotInterpOrder,Potknots,PotInterpPoints,Pbcoef(1,1,:)), ! 
     >           dbsval(R(i),PotInterpOrder,Potknots,PotInterpPoints,Pbcoef(1,2,:)), 
     >           dbsval(R(i),PotInterpOrder,Potknots,PotInterpPoints,Pbcoef(2,1,:)) 
c            x=x+0.05d0
         enddo
c         write(20,*)
c      enddo
      deallocate(P,Q,QTemp,dP)
      deallocate(Ucurves)

      allocate(SectorR(NumSectors+1))



      SectorR(1)=RMin
      SectorR(2)=RMax
      print*,'SectorR(',1,') = ',SectorR(1)
      print*,'SectorR(',2,') = ',SectorR(2)


      allocate(BigPot(NumChan,NumChan,LegPoints,RNumPoints,NumSectors))
      allocate(BigAPot(NumChan,NumChan,LegPoints,RNumPoints,NumSectors))
      iEnergy=1
      

      Left=0
      Right=0
      RDimMin=RNumPoints+order-3
      RDim=RDimMin
      if (Left .eq. 2) RDim = RDim + 1
      if (Right .eq. 2) RDim = RDim + 1
      iSector=1
      print*, 'iSector = ',iSector,'  #Left BC  Right BC = ', Left, Right, 'RNumPoints = ', RNumPoints ! 
      print*, '# RDimMin = ',RDimMin,' RDim = ',RDim
      
      
      call CalcBoundStates(NumSectors,alpha,mu,iSector,SectorR(iSector),SectorR(iSector+1),! 
     >     Left,Right,RDim,NumChan,RNumPoints,xLeg,wLeg,LegPoints,Order,  ! 
     >     Ebcoef,Pbcoef,Qbcoef,Potknots,PotInterpPoints,PotInterpOrder,BigPot,BigAPot) ! 
      
      deallocate(BigPot,BigAPot)
c 10   format(1P,100e25.15)
 20   format(1P,100e16.8)
      
 1002 format(a64)

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CalcBoundStates(NumSectors,alpha,mu,iSector,RMin,RMax,Left,Right,
     >     RDim,NumChan,RNumPoints,xLeg,wLeg,LegPoints,Order,
     >     Ebcoef,Pbcoef,Qbcoef,Potknots,PotInterpPoints,PotInterpOrder,BigPot,BigAPot) ! 
      implicit none
      double precision, external :: phirecon, kdelta
      integer Left,Right,beta,betaMax,NumChan,MatrixDim,RDim,RNumPoints,Order, NumSectors,! 
     >     LegPoints,PotInterpPoints,Num2BodyChan,i,j,k,nch,mch,iSector,PotInterpOrder,iEnergy ! 
      double precision RMin,RMax,RMatOld(NumChan,NumChan),RMatNew(NumChan,NumChan),! 
     >     xLeg(LegPoints),wLeg(LegPoints),alpha,energy,mu,x
      double precision Ebcoef(NumChan,PotInterpPoints),Pbcoef(NumChan,NumChan,PotInterpPoints), ! 
     >     Qbcoef(NumChan,NumChan,PotInterpPoints),Potknots(PotInterpPoints+PotInterpOrder) ! 
      double precision, allocatable :: u(:,:,:),ux(:,:,:),uxx(:,:,:),u0(:,:,:),u0x(:,:,:),RPoints(:),asymptoticE(:) ! 
      double precision, allocatable :: H(:,:),S(:,:),evec(:,:),eval(:),RMatLocal(:,:),Kmat(:,:),Lambda(:,:) ! 
      double precision hbarc,BigPot(NumChan,NumChan,LegPoints,RNumPoints,NumSectors),
     >     BigAPot(NumChan,NumChan,LegPoints,RNumPoints,NumSectors) ! 
      common hbarc

      
      integer, allocatable :: RBounds(:),ikeep(:)


      MatrixDim = RDim*NumChan
      allocate(u(LegPoints,RNumPoints,RDim),ux(LegPoints,RNumPoints,RDim),uxx(LegPoints,RNumPoints,RDim)) ! 
      allocate(RPoints(RNumPoints))
      allocate(RBounds(RNumPoints+2*Order))
      print*, '#Calling GridMaker'
      
      call GridMaker(RNumPoints,RMin,RMax,RPoints)
      print*, '#Calculating basis functions'

      call setup_potential_matrix(BigPot,BigAPot,mu,NumChan,NumSectors,LegPoints,RNumPoints,PotInterpOrder,Potknots, ! 
     >     PotInterpPoints,Ebcoef,Pbcoef,Qbcoef,RPoints,iSector,xLeg)

      call CalcBasisFuncs(Left,Right,Order,RPoints,LegPoints,xLeg,
     >     RDim,RBounds,RNumPoints,0,u)
      call CalcBasisFuncs(Left,Right,Order,RPoints,LegPoints,xLeg,
     >     RDim,RBounds,RNumPoints,1,ux)
      call CalcBasisFuncs(Left,Right,Order,RPoints,LegPoints,xLeg,
     >     RDim,RBounds,RNumPoints,2,uxx)

c      call CheckBasis(u,RDim,RNumPoints,LegPoints,xLeg,RPoints,999)

      allocate(S(MatrixDim,MatrixDim),H(MatrixDim,MatrixDim)) ! 
      Num2BodyChan=1
      print*, 'RDim=',RDim
      print*, '#Calculating Hamiltonian matrix'


      call CalcHamiltonian(RMin,RMax,mu,NumSectors,alpha,energy,Left,Right,Order,Num2BodyChan,RPoints,LegPoints,xLeg,wLeg,RDim, ! 
     >     RNumPoints,u,ux,uxx,RBounds,MatrixDim,NumChan,H,S,Ebcoef,Pbcoef,Qbcoef,Potknots,PotInterpPoints,PotInterpOrder, ! 
     >     BigPot,BigAPot,iSector)                    ! 

      deallocate(u,ux,uxx)
      deallocate(RBounds)
      
      print*, '#Diagonalizing System'
      
      allocate(evec(MatrixDim,MatrixDim),eval(MatrixDim))
      call Mydggev(MatrixDim,H,MatrixDim,S,MatrixDim,eval,evec)

      print*, '# done diagonalizing... deallocating memory' ! 
      
c      allocate(ikeep(2*NumChan))
c      call CheckBasisPhi(RMin,RMax,Left,Right,RDim,RNumPoints,RPoints,0,Order,666)

      j=1
      do i = 1, 10
         if(eval(i).le.0.0d0) then
            print*,'# ****************  eval(',i,')=',eval(i)
         endif
      enddo

      x=RMin
      do while (x.lt.RMax)
         write(777,*), x, phirecon(x,1,1,evec,Left,Right,RDim,MatrixDim,RNumPoints,RPoints,order) ! 
         x=x+0.05
      enddo

      deallocate(H,S)


      deallocate(evec,eval)      
      deallocate(RPoints)

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine setup_interp(R,PotInterpOrder,Potknots,PotInterpPoints,P,Q,Ucurves,Ebcoef,Pbcoef,Qbcoef, Numchan)
      use bspline
      implicit none
      integer PotInterpPoints, NumChan,PotInterpOrder
      double precision Ucurves(NumChan,PotInterpPoints),P(NumChan,NumChan,PotInterpPoints) ! 
      double precision Q(NumChan,NumChan,PotInterpPoints),R(PotInterpPoints) ! 
      double precision Potknots(PotInterpPoints+PotInterpOrder), Ebcoef(NumChan,PotInterpPoints) ! 
      double precision Pbcoef(NumChan,NumChan,PotInterpPoints),Qbcoef(NumChan,NumChan,PotInterpPoints) ! 
      integer nch, mch


      call dbsnak(PotInterpPoints, R, PotInterpOrder, Potknots)      

      do nch=1, NumChan
         call dbsint(PotInterpPoints,R,Ucurves(nch,:),PotInterpOrder,Potknots,Ebcoef(nch,:)) ! 
         do mch=1, NumChan
            call dbsint(PotInterpPoints,R,P(nch,mch,:),PotInterpOrder,Potknots,Pbcoef(nch,mch,:)) ! 
            call dbsint(PotInterpPoints,R,Q(nch,mch,:),PotInterpOrder,Potknots,Qbcoef(nch,mch,:)) ! 
         enddo
      enddo
            
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine ReadCurves(PotInterpPoints,R,Ucurves,P,QTemp,Q,NumChan,NumTot) ! 
      implicit none
      integer PotInterpPoints,NumChan,NumTot,i,j,iR
      double precision R(PotInterpPoints),Ucurves(NumChan,PotInterpPoints) ! 
      double precision P(NumChan,NumChan,PotInterpPoints),Q(NumChan,NumChan,PotInterpPoints) ! 
      double precision QTemp(NumChan,NumChan,PotInterpPoints)
      double precision, allocatable :: dP(:,:,:)

      allocate(dP(NumChan,NumChan,PotInterpPoints))

      do iR = 1, PotInterpPoints
         read(101,*) R(iR)
         read(102,*) R(iR)
         read(103,*) R(iR)
         read(200,*) R(iR), (Ucurves(i,iR),i=1,NumChan)
         do i = 1, NumChan
            read(101,*) (P(i,j,iR), j = 1,NumChan)
c            read(101,*)
            read(102,*) (QTemp(i,j,iR), j = 1,NumChan)
c            read(102,*)
            read(103,*) (dP(i,j,iR), j = 1,NumChan)
         enddo
         do i=NumChan+1,NumTot
            read(101,*)
            read(102,*)
            read(103,*)
         enddo
      enddo
      print*, '# ...Done'
c      print*, '# symeterizing the Q matrix '
      do iR = 1, PotInterpPoints
         do i = 1, NumChan
            do j = 1, NumChan
               
c     Q(i,j,iR) = 0.5d0*(QTemp(i,j,iR)+QTemp(j,i,iR))
c     dP(i,j,iR)=0.0d0
c     Q(i,j,iR)=0.0d0
c     THIS IS THE ACTUAL 2ND DERIV COUPLING IN TERMS OF THE SYM AND ASYM PIECES <phi|d2(phi)/dR2>
               Q(i,j,iR) = (QTemp(i,j,iR)-dP(i,j,iR))
c     Q(i,j,iR) = -(dP(i,j,iR))
c     Q(i,j,iR) = QTemp(i,j,iR)
c     Q(i,j,iR) = 0.0d0
c     P(i,j,iR) = 0.0d0
c     Ucurves(i,iR) = 0.0d0
               
               
            enddo
         enddo
c     USE THIS BIT OF CODE TO FORCE THE P-MATRIX TO BE ANTISYMMETRIC:
         do i = 1, NumChan
            do j = i, NumChan
               P(i,j,iR) = 0.5d0*(P(i,j,iR) - P(j,i,iR))
               P(j,i,iR) = -P(i,j,iR)
            enddo
         enddo

      enddo


      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine GridMaker(xNumPoints,xMin,xMax,xPoints)
      implicit none
      integer xNumPoints
      double precision xMin,xMax,xPoints(xNumPoints)

      integer i,j,k,OPGRID
      double precision Pi
      double precision r0New
      double precision xRswitch
      double precision xDelt,x0,x1,x2


      Pi = 3.1415926535897932385d0


      x0 = xMin
      x1 = xMax
      k = 1
      xDelt = (x1-x0)/dble(xNumPoints-1)
      do i = 1,xNumPoints
c         xPoints(k) = x1*((i-1)*xDelt/(x1-x0))**2 + x0
         xPoints(k) = (i-1)*xDelt + x0
c         print*, k, xPoints(k)
         k = k + 1

      enddo
      OPGRID=1
      if(OPGRID.eq.1) then

         r0New=100.0d0
         x0 = xMin
         x1 = r0New
         x2 = xMax
         print*, x0,x1,x2
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints/2)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x0
c            print*, k, xPoints(k), xDelt
            k = k + 1
         enddo
         xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x1
c            print*, k, xPoints(k), xDelt
            k = k + 1
         enddo
      endif


 15   format(6(1x,1pd12.5))
      

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      subroutine CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,
     >     xNumPoints,u,xBounds,MatrixDim,NumChan,S)
      implicit none
      integer Order,LegPoints,xDim,xNumPoints,xBounds(xNumPoints+2*Order),MatrixDim,NumChan
      double precision xPoints(*),xLeg(*),wLeg(*)
      double precision S(MatrixDim,MatrixDim)
      double precision u(LegPoints,xNumPoints,xDim)

      integer ix,ixp,kx,lx, nch, mch
      integer i1,i1p
      integer Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:)
      double precision a,b,m
      double precision xTempS
      double precision ax,bx
      double precision, allocatable :: xIntScale(:)
   
      allocate(xIntScale(xNumPoints))
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))

      S = 0.0d0

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xIntScale(kx) = 0.5d0*(bx-ax)
      enddo

      do ix = 1,xDim
         do ixp = 1,xDim
            kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
            kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
         enddo
      enddo

      do nch = 1,NumChan
         do mch = 1,NumChan
            do ix = 1,xDim
               do ixp = max(1,ix-Order),min(xDim,ix+Order)
c                  xS(ixp,ix) = 0.0d0
                  S((NumChan-1)*xDim+ix,(NumChan-1)*xDim+ixp) = 0.0d0
                  do kx = kxMin(ixp,ix),kxMax(ixp,ix)
                     xTempS = 0.0d0

                     do lx = 1,LegPoints
                        a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
                        b = a*u(lx,kx,ixp)
                        xTempS = xTempS + b
                     enddo
c                     xS(ixp,ix) = xS(ixp,ix) +   xTempS
                     S((NumChan-1)*xDim+ix,(NumChan-1)*xDim+ixp) = S((NumChan-1)*xDim+ix,(NumChan-1)*xDim+ixp)
     >                    + xTempS
                  enddo
               enddo
            enddo
         enddo
      enddo


      deallocate(xIntScale)
      deallocate(kxMin,kxMax)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CalcHamiltonian(RMin,RMax,mu,NumSectors,alpha,energy,Left,Right, ! 
     >     Order,Num2BodyChan,RPoints,LegPoints,xLeg,wLeg,RDim, ! 
     >     RNumPoints,u,ux,uxx,RBounds,MatrixDim,NumChan,H,S,Ebcoef,Pbcoef,Qbcoef,Potknots,! 
     >     PotInterpPoints,PotInterpOrder, ! 
     >     BigPot,BigAPot,iSector)                    ! 
      use bspline
      implicit none
      double precision, external :: kdelta
      integer PotInterpPoints, MatrixDim, RNumPoints,kx,lx,RDim,Order,ix,ixp,Left,Right,LegPoints,NumChan,mch,nch ! 
      integer PotInterpOrder,iSector,NumSectors
      double precision RPoints(RNumPoints),xLeg(LegPoints),wLeg(LegPoints),alpha,mu,energy,a,prefact ! 
      double precision Ebcoef(NumChan,PotInterpPoints),Pbcoef(NumChan,NumChan,PotInterpPoints), ! 
     >     Qbcoef(NumChan,NumChan,PotInterpPoints),Potknots(PotInterpPoints+PotInterpOrder) ! 
      double precision u(LegPoints,RNumPoints,RDim), uxx(LegPoints,RNumPoints,RDim), ! 
     >     ux(LegPoints,RNumPoints,RDim),RMax,RMin ! 
      double precision x, xScaledZero,Znet,hbarc
      double precision, allocatable :: Pot(:,:,:,:),APot(:,:,:,:)
      double precision H(MatrixDim,MatrixDim), S(MatrixDim,MatrixDim),ax, bx,potvalue ! 
      double precision, allocatable :: xIntScale(:),xIntPoints(:,:)
      double precision, allocatable ::  threshold(:) 
      integer, allocatable :: kxMin(:,:),kxMax(:,:)
      common hbarc
      double precision TempH, TempS, test,BigPot(NumChan,NumChan,LegPoints,RNumPoints,NumSectors), ! 
     >     BigAPot(NumChan,NumChan,LegPoints,RNumPoints,NumSectors)                     ! 
      integer RBounds(RNumPoints+2*Order),Num2BodyChan,i
c      common PotInterpOrder
      print*, 'NumChan=',NumChan

      allocate(xIntScale(RNumPoints))
      allocate(kxMin(RDim,RDim),kxMax(RDim,RDim)) ! 
      allocate(xIntPoints(LegPoints,RNumPoints))

      print*,'# NumChan = ', NumChan
      print*,'# MatrixDim = ',MatrixDim

      prefact=0.5d0/mu      

c      print*,'#seting up kxMin and kxMax'
      do ix = 1,RDim
         do ixp = 1,RDim
            kxMin(ixp,ix) = max(RBounds(ix),RBounds(ixp))
            kxMax(ixp,ix) = min(RBounds(ix+Order+1),RBounds(ixp+Order+1))-1 ! 
c            print*,ixp,ix,kxMin(ixp,ix),kxMax(ixp,ix)
         enddo
      enddo
      
      do mch = 1, NumChan
         do nch = 1, NumChan
            do kx = 1,RNumPoints-1
               ax = RPoints(kx)
               bx = RPoints(kx+1)
               xIntScale(kx) = 0.5d0*(bx-ax)
               xScaledZero = 0.5d0*(bx+ax)
               do lx = 1,LegPoints
                  xIntPoints(lx,kx) = xIntScale(kx)*xLeg(lx)+xScaledZero
               enddo
            enddo

            do ix = 1,RDim
               do ixp = 1,RDim
                  S((mch-1)*RDim+ix,(nch-1)*RDim+ixp)=0.0d0
                  H((mch-1)*RDim+ix,(nch-1)*RDim+ixp)=0.0d0
               enddo
            enddo
         enddo
      enddo
      
c      call CheckPot(BigPot,RNumPoints,NumChan,LegPoints,xLeg,RPoints,iSector) ! 
      
      print*, '#Calculating Hamiltonian Matrix elements'
      
      do ix = 1,RDim
         do ixp = max(1,ix-Order),min(RDim,ix+Order)
            
            do mch = 1,NumChan
               do nch = 1, NumChan
c                  print*, 'mch, nch = ',mch, nch
                  do kx = kxMin(ixp,ix),kxMax(ixp,ix)
                     TempH = 0.0d0
                     TempS = 0.0d0
                     do lx = 1,LegPoints
                        a = xIntPoints(lx,kx)*wLeg(lx)*xIntScale(kx)
                        TempS = TempS + a*kdelta(mch,nch)*u(lx,kx,ix)*u(lx,kx,ixp) ! 
                        TempH = TempH - a*prefact*kdelta(mch,nch)*u(lx,kx,ix)*(uxx(lx,kx,ixp) + ! 
     >                       ux(lx,kx,ixp)/xIntPoints(lx,kx)) ! 
                        TempH = TempH + a*prefact*u(lx,kx,ix)*(BigPot(mch,nch,lx,kx,iSector))*u(lx,kx,ixp) ! sym piece
                        TempH = TempH - a*prefact*BigAPot(mch,nch,lx,kx,iSector)* ! antisym P matrix 
     >                       (2.0*u(lx,kx,ix)*ux(lx,kx,ixp)+u(lx,kx,ix)*u(lx,kx,ixp)/xIntPoints(lx,kx)) ! cont.
c                    print*, lx
                     enddo
                     H((mch-1)*RDim+ix,(nch-1)*RDim+ixp) = H((mch-1)*RDim+ix,(nch-1)*RDim+ixp) + TempH ! 
                     S((mch-1)*RDim+ix,(nch-1)*RDim+ixp) = S((mch-1)*RDim+ix,(nch-1)*RDim+ixp) + TempS ! 
                  enddo
               enddo
            enddo
            
         enddo
      enddo
      

      print*, '#...Done'
      deallocate(kxMax,kxMin)
      deallocate(xIntScale,xIntPoints)

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function kdelta(mch,nch)
      
      integer mch,nch
      if (mch.eq.nch) then
         kdelta = 1.0d0
      else 
         kdelta = 0.0d0
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Mydggev(N,H,LDH,L,LDL,eval,evec)
      implicit none
      integer LDH,N,LDL,info
      double precision H(LDH,N),L(LDL,N),eval(N),evec(N,N)
      double precision, allocatable :: alphar(:),alphai(:),beta(:),work(:),VL(:,:)
      integer lwork,i,im,in,j

      allocate(alphar(N),alphai(N),beta(N))


      info = 0
      lwork = -1
      allocate(work(1))
      call dggev('N','V',N,H,LDH,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info) ! 
      do im = 1,N
         alphar(im)=0.0d0
         alphai(im)=0.0d0
         beta(im)=0.0d0
         eval(im)=0.0d0
         do in = 1,N
            evec(im,in)=0.0d0
         enddo
      enddo

      lwork=work(1)
      deallocate(work)
      allocate(work(lwork))
      call dggev('N','V',N,H,LDH,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info) ! 
      
      do i = 1, N
         if (abs(alphai(i)).ge.1e-12) then
            print*, '#eigenvalue may be complex! alphai(',i,')=',alphai(i) ! 
         endif
         if(abs(beta(i)).ge.1e-12) then
            eval(i) = -alphar(i)/beta(i)
         endif
      enddo

      call deigsrt(eval,evec,N,N)

      do i = 1,N
         eval(i)=-eval(i)
      enddo

      deallocate(alphar,alphai,beta)
      deallocate(work)

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CompSqrMatInv(A, N)
      implicit none
      integer N,info,lwk
      integer, allocatable :: ipiv(:)
      double precision, allocatable :: work(:)
      double complex A(N,N)
      allocate(ipiv(N))
      call zgetrf(N, N, A, N, ipiv, info)
      allocate(work(1))
      lwk = -1
      call zgetri(N, A, N, ipiv, work, lwk, info)
      lwk = work(1)
      deallocate(work)
      allocate(work(lwk))
      call zgetri(N, A, N, ipiv, work, lwk, info)
      deallocate(ipiv,work)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CheckPot(Pot,RNumPoints,NumChan,LegPoints,xLeg,RPoints,iSector) ! 
      implicit none
      integer RNumPoints,NumChan,LegPoints,mch,nch,kx,lx,iSector
      double precision Pot(NumChan,NumChan,LegPoints,RNumPoints,1),ax,bx,xLeg(LegPoints)
      double precision RPoints(RNumPoints),xScaledZero,x
      double precision, allocatable :: xIntScale(:)
      allocate(xIntScale(RNumPoints))
      do mch = 1,NumChan
         do nch = 1,NumChan
            do kx = 1,RNumPoints-1
               ax = RPoints(kx)
               bx = RPoints(kx+1)
               xIntScale(kx) = 0.5d0*(bx-ax)
               xScaledZero = 0.5d0*(bx+ax)
               do lx = 1,LegPoints
                  x = xIntScale(kx)*xLeg(lx)+xScaledZero
                  write(888,*), x, Pot(mch,nch,lx,kx,iSector)
               enddo
            enddo
            write(888,*), ' '
         enddo
      enddo
      
      deallocate(xIntScale)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CheckBasisPhi(RMin,RMax,Left,Right,RDim,RNumPoints,RPoints,Deriv,order,file)
      double precision, external :: BasisPhi
      integer MatrixDim,RDim,nch,beta,i,RNumPoints,Left,Right,Deriv,order,file
      double precision R,RMax,RPoints(RNumPoints)

      do ix=1,RDim
         R=RMin
         do while (R.le.RMax)
            write(file,*) R, BasisPhi(R,Left,Right,order,RDim,RPoints,RNumPoints,Deriv,ix)      
            R = R+0.0001d0
         enddo
         write(file,*) ' '
      enddo

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CheckBasis(u,RDim,RNumPoints,LegPoints,xLeg,RPoints,file)
      implicit none
      integer RNumPoints,NumChan,LegPoints,kx,lx,ix,RDim,file
      double precision u(LegPoints,RNumPoints,RDim),ax,bx,xLeg(LegPoints)
      double precision RPoints(RNumPoints),xScaledZero,x
      double precision, allocatable :: xIntScale(:)
      allocate(xIntScale(RNumPoints))

      do ix=1,RDim
         do kx = 1,RNumPoints-1
            ax = RPoints(kx)
            bx = RPoints(kx+1)
            xIntScale(kx) = 0.5d0*(bx-ax)
            xScaledZero = 0.5d0*(bx+ax)
            do lx = 1,LegPoints
               x = xIntScale(kx)*xLeg(lx)+xScaledZero
               write(file,*), x, u(lx,kx,ix)
            enddo
         enddo
         write(file,*), ' '
      enddo      
      
      deallocate(xIntScale)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function phirecon(R,beta,nch,evec,left,right,RDim,MatrixDim,RNumPoints,RPoints,order)
      implicit none
      double precision, external :: BasisPhi
      integer MatrixDim,RDim,nch,beta,i,RNumPoints,left,right,order
      double precision R,evec(MatrixDim,MatrixDim),RPoints(RNumPoints)
      phirecon = 0.0d0
      do i = 1,RDim
         phirecon = phirecon + evec((nch-1)*RDim+i,beta)*BasisPhi(R,left,right,order,RDim,RPoints,
     >        RNumPoints,0,i)
      enddo
c      if(R.ne.0.0d0) then
c         phirecon = phirecon/dsqrt(R) 
c      endif
c      print*, R, phirecon
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine SetMultipoleCoup(MultipoleCoupMat,threshold,NumChan,NumLam)
      implicit none
      integer NumChan,NumLam
      double precision MultipoleCoupMat(NumChan,NumChan,NumLam), threshold(NumChan)

      threshold(1)=0.0d0
c      threshold(2)=0.0d0
      threshold(2)=-0.2d0
c      threshold(3)=0.5d0*2.169d0
      threshold(3)=-0.1d0
      

      MultipoleCoupMat(1,1,1)=0.0d0
      MultipoleCoupMat(1,2,1)=0.0d0
      MultipoleCoupMat(1,3,1)=-0.5682977d0
      MultipoleCoupMat(2,1,1)=0.0d0
      MultipoleCoupMat(2,2,1)=0.0d0
      MultipoleCoupMat(2,3,1)=-0.8036944d0
      MultipoleCoupMat(3,1,1)=-0.5682977d0
      MultipoleCoupMat(3,2,1)=-0.8036944d0
      MultipoleCoupMat(3,3,1)=0.0d0

      MultipoleCoupMat(1,1,2)=0.0d0
      MultipoleCoupMat(1,2,2)=-0.5554260d0
      MultipoleCoupMat(1,3,2)=0.0d0
      MultipoleCoupMat(2,1,2)=-0.5554260d0
      MultipoleCoupMat(2,2,2)=-0.3927455d0
      MultipoleCoupMat(2,3,2)=0.0d0
      MultipoleCoupMat(3,1,2)=0.0d0
      MultipoleCoupMat(3,2,2)=0.0d0
      MultipoleCoupMat(3,3,2)=0.0d0


      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine printmatrix(M,nr,nc,file)
      implicit none
      integer nr,nc,file,j,k
      double precision M(nr,nc)
      
      do j = 1,nr
         write(file,*) (M(j,k), k = 1,nc)
      enddo
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copymatrix(A,B,nr,nc)
      implicit none
      integer nr,nc,i,j
      double precision A(nr,nc),B(nr,nc)

      do i=1,nr
         do j=1,nc
            B(i,j)=A(i,j)
         enddo
      enddo
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ModifyBasis(u,ux,RDim,RNumPoints,LegPoints,xLeg,RPoints,Order)
      implicit none
      double precision, external  :: BasisPhi
      integer RNumPoints,NumChan,LegPoints,kx,lx,ix,RDim,Order
      double precision u(LegPoints,RNumPoints,RDim),ux(LegPoints,RNumPoints,RDim),
     >     u0(LegPoints,RNumPoints,RDim),u0x(LegPoints,RNumPoints,RDim),ax,bx,xLeg(LegPoints)
      double precision RPoints(RNumPoints),xScaledZero,x
      double precision, allocatable :: xIntScale(:)
      allocate(xIntScale(RNumPoints))

      do ix=1,RDim!Order-1
         do kx = 1,RNumPoints-1
            ax = RPoints(kx)
            bx = RPoints(kx+1)
            xIntScale(kx) = 0.5d0*(bx-ax)
            xScaledZero = 0.5d0*(bx+ax)
            do lx = 1,LegPoints
               x = xIntScale(kx)*xLeg(lx)+xScaledZero
               ux(lx,kx,ix)=u(lx,kx,ix)/(2.0d0*dsqrt(x))+ux(lx,kx,ix)*dsqrt(x) ! 
               u(lx,kx,ix)=dsqrt(x)*u(lx,kx,ix)
            enddo
         enddo
      enddo      
      deallocate(xIntScale)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setup_potential_matrix(Pot,APot,mu,NumChan,NumSectors,LegPoints,RNumPoints,PotInterpOrder,Potknots,
     >     PotInterpPoints,Ebcoef,Pbcoef,Qbcoef,RPoints,iSector,xLeg)
      use bspline
      implicit none
      integer NumChan,NumSectors,LegPoints,RNumPoints,PotInterpOrder,lx,kx,mch,nch,PotInterpPoints, ! 
     >     iSector                     ! 
      double precision, external :: kdelta
      double precision Ebcoef(NumChan,PotInterpPoints),Qbcoef(NumChan,NumChan,PotInterpPoints), mu,! 
     >     Pbcoef(NumChan,NumChan,PotInterpPoints)
c      double precision xIntPoints(LegPoints,RNumPoints)
      double precision, allocatable :: xIntScale(:),xIntPoints(:,:)
      double precision Pot(NumChan,NumChan,LegPoints,RNumPoints,NumSectors),ax,bx,xScaledZero,xLeg(LegPoints) ! 
      double precision APot(NumChan,NumChan,LegPoints,RNumPoints,NumSectors),RPoints(RNumPoints) ! 
      double precision Potknots(PotInterpPoints+PotInterpOrder)
      allocate(xIntPoints(LegPoints,RNumPoints))
      allocate(xIntScale(RNumPoints))
      print*,'#setting up potential matrix,iSector = ',iSector

      do mch = 1, NumChan
         do nch = 1, NumChan
c            print*,'kx,  lx, ax,              bx,            xIntScale,     xScaledZero,     x' ! 
            do kx = 1,RNumPoints-1
               ax = RPoints(kx)
               bx = RPoints(kx+1)
c               print*, 'RPoints(1)=',RPoints(kx), 'RPoints(2)=',RPoints(kx+1)
               xIntScale(kx) = 0.5d0*(bx-ax)
               xScaledZero = 0.5d0*(bx+ax)
               do lx = 1,LegPoints

                  Pot(mch,nch,lx,kx,iSector)=0.0d0
                  APot(mch,nch,lx,kx,iSector)=0.0d0

                  xIntPoints(lx,kx) = xIntScale(kx)*xLeg(lx)+xScaledZero
                  
c     print*, kx, lx, ax, bx, xIntScale, xScaledZero, xIntPoints(lx,kx) ! 
                  
                  Pot(mch,nch,lx,kx,iSector) = kdelta(mch,nch)*2.0d0*mu* ! 
     >                 dbsval(xIntPoints(lx,kx),PotInterpOrder,Potknots,PotInterpPoints,Ebcoef(mch,:)) ! 
                  
                  
                  Pot(mch,nch,lx,kx,iSector) = Pot(mch,nch,lx,kx,iSector) - ! 
     >                 dbsval(xIntPoints(lx,kx),PotInterpOrder,Potknots,PotInterpPoints,Qbcoef(mch,nch,:)) ! 
                  
                  APot(mch,nch,lx,kx,iSector) = dbsval(xIntPoints(lx,kx),
     >                 PotInterpOrder,Potknots,PotInterpPoints,Pbcoef(mch,nch,:)) ! set the P matrix

                  
               enddo
            enddo
            
         enddo
      enddo
      deallocate(xIntScale,xIntPoints)
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE deigsrt(d,v,n,np)
      INTEGER n,np
      double precision d(np),v(np,np)
      INTEGER i,j,k
      double precision p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
