!c THIS PROGRAM CALCULATES THE BOUND STATES SUPPORTED BY A MULTICHANNEL SET OF 
!C POTENTIALS THAT ARE SPECIFIED BY A SUBROUTINE.
!c NB: Because this is a bound-state code, it does NOT do R-matrix poropagation and therefore
!C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Module modb
  implicit none
  real*8 pi,t0,xm0,b0,ee0,d0,f0,w0,sig0,rtol,bwidth
  integer nleg
  Parameter(pi=3.1415926535897932d0,nleg=16,ee0 = 6.579685d6)
  Parameter(t0=3.15934d5,xm0=5.48579903d-4,b0 = 4.70110d9/2.d0)
  PARAMETER(d0=2.541d0,f0=5.1456d9,w0=2.195d5)
  PARAMETER(sig0=2.8002852e-17)
  PARAMETER(rtol=1d-5, bwidth=1d-3)
end Module modb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program MCBound
  implicit none
  call Setup
end program MCBound

subroutine Setup
  implicit none
  integer LegPoints,RNumPoints,NumChan
  integer Order, RDimMin,n
  integer Left,Right,RDim
  integer i,j,k,iR,Left0
  double precision RMin, RMax,alpha,delR,mass,EnergyStart,EnergyEnd,deltaEnergy ! 
  double precision mu, mu12, Pi,energy,x,RMaxRun
  double precision, allocatable :: xLeg(:),wLeg(:)
  double precision, allocatable :: kchan(:),BigPot(:,:,:,:)


  character*64 LegendreFile
  double precision hbarc
  double complex II
  parameter(II=(0.0d0,1.0d0))
  common hbarc

  !c      hbarc=197.32858d0
  hbarc=1.0d0

  !c     read in number of channels and order of splines
  print*, '#NumChan, Order'
  read(5,*)
  read(5,*)  NumChan, order
  write(6,*) NumChan, order
  print*, '#Legendre File'
  read(5,*)
  read(5,*)
  read(5,1002) LegendreFile
  write(6,1002) LegendreFile

  read(5,*)
  read(5,*)
  read(5,*) alpha, LegPoints
  write(6,*) "alpha, LegPoints"
  print*, alpha, LegPoints
!  Pi=dacos(-1.d0)
!  write(6,*) '#Pi=',Pi

  !c     read in domain information
  read(5,*)
  read(5,*)
  read(5,*)  left,right, mass
  write(6,*) "Left, Right, mass"
  write(6,*) left,right, mass
  read(5,*)
  read(5,*)
  read(5,*)  RMin,RMax, RNumPoints
  write(6,*) "RMin, RMax, RNumPoints"
  write(6,*) RMin,RMax, RNumPoints


  !mu=mass/dsqrt(3.0d0)
  mu = mass/2.d0
  print*, '# Reading in Gauss-Legendre points and weights'
  allocate(xLeg(LegPoints),wLeg(LegPoints))
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  allocate(BigPot(NumChan,NumChan,LegPoints,RNumPoints))

!  Left=0
!  Right=0
  RDimMin=RNumPoints+order-3
  RDim=RDimMin
  if (Left .eq. 2) RDim = RDim + 1
  if (Right .eq. 2) RDim = RDim + 1

  print*,  '#Left BC  Right BC = ', Left, Right, 'RNumPoints = ', RNumPoints ! 
  print*, '# RDimMin = ',RDimMin,' RDim = ',RDim

  call CalcBoundStates(alpha,mu,RMin,RMax, &! 
       Left,Right,RDim,NumChan,RNumPoints,xLeg,wLeg,LegPoints,Order,BigPot)
  
  deallocate(BigPot)
  !c     10   format(1P,100e25.15)
20 format(1P,100e16.8)

1002 format(a64)

end subroutine Setup

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcBoundStates(alpha,mu,RMin,RMax,Left,Right,&
     RDim,NumChan,RNumPoints,xLeg,wLeg,LegPoints,Order,BigPot)
  implicit none
  double precision, external :: phirecon, kdelta
  integer Left,Right,beta,betaMax,NumChan,MatrixDim,RDim,RNumPoints,Order,&
       LegPoints,Num2BodyChan,i,j,k,nch,mch
  double precision RMin,RMax,xLeg(LegPoints),wLeg(LegPoints),alpha,energy,mu,x
  double precision, allocatable :: u(:,:,:),ux(:,:,:),uxx(:,:,:),u0(:,:,:),u0x(:,:,:),RPoints(:),Ethresh(:)
  double precision, allocatable :: H(:,:),S(:,:),evec(:,:),eval(:)
  double precision hbarc,BigPot(NumChan,NumChan,LegPoints,RNumPoints)
  common hbarc
  integer, allocatable :: RBounds(:),ikeep(:)

  MatrixDim = RDim*NumChan
  allocate(u(LegPoints,RNumPoints,RDim),ux(LegPoints,RNumPoints,RDim),uxx(LegPoints,RNumPoints,RDim)) ! 
  allocate(RPoints(RNumPoints))
  allocate(RBounds(RNumPoints+2*Order))
  print*, '#Calling GridMaker'

  call GridMaker(RNumPoints,RMin,RMax,RPoints)
  print*, '#Calculating basis functions'

  call setup_potential_matrix(BigPot,mu,NumChan,LegPoints,RNumPoints,RPoints,xLeg)

  call CalcBasisFuncs(Left,Right,Order,RPoints,LegPoints,xLeg,RDim,RBounds,RNumPoints,0,u)
  call CalcBasisFuncs(Left,Right,Order,RPoints,LegPoints,xLeg,RDim,RBounds,RNumPoints,1,ux)
  call CalcBasisFuncs(Left,Right,Order,RPoints,LegPoints,xLeg,RDim,RBounds,RNumPoints,2,uxx)

  !c      call CheckBasis(u,RDim,RNumPoints,LegPoints,xLeg,RPoints,999)

  allocate(S(MatrixDim,MatrixDim),H(MatrixDim,MatrixDim)) ! 

  print*, 'RDim=',RDim
  print*, '#Calculating Hamiltonian matrix'

  call CalcHamiltonian(RMin,RMax,mu,alpha,energy,Left,Right,Order,RPoints,LegPoints,xLeg,wLeg,RDim,&
       RNumPoints,u,ux,uxx,RBounds,MatrixDim,NumChan,H,S,BigPot) ! 

  deallocate(u,ux,uxx)
  deallocate(RBounds)

  print*, '#Diagonalizing System'

  allocate(evec(MatrixDim,MatrixDim),eval(MatrixDim))
  call Mydggev(MatrixDim,H,MatrixDim,S,MatrixDim,eval,evec)

  print*, '# done diagonalizing... deallocating memory' ! 

  !c      allocate(ikeep(2*NumChan))
  !c      call CheckBasisPhi(RMin,RMax,Left,Right,RDim,RNumPoints,RPoints,0,Order,666)

  j=1
  do i = 1, MatrixDim
!     if(eval(i).le.0.0d0) then
        print*,'# ****************  eval(',i,')=',eval(i)
!     endif
  enddo

  x=RMin
  do while (x.lt.RMax)
     write(777,*) x, phirecon(x,1,1,evec,Left,Right,RDim,MatrixDim,RNumPoints,RPoints,order) ! 
     x=x+0.05
  enddo

  deallocate(H,S)

  deallocate(evec,eval)      
  deallocate(RPoints)

end subroutine CalcBoundStates
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
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
     !c         xPoints(k) = x1*((i-1)*xDelt/(x1-x0))**2 + x0
     xPoints(k) = (i-1)*xDelt + x0
     !c         print*, k, xPoints(k)
     k = k + 1

  enddo
  OPGRID=1
  if(OPGRID.eq.1) then

     r0New=1.0d0
     x0 = xMin
     x1 = r0New
     x2 = xMax
     print*, x0,x1,x2
     k = 1
     xDelt = (x1-x0)/dfloat(xNumPoints/2)
     do i = 1,xNumPoints/2
        xPoints(k) = (i-1)*xDelt + x0
        !c            print*, k, xPoints(k), xDelt
        k = k + 1
     enddo
     xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
     do i = 1,xNumPoints/2
        xPoints(k) = (i-1)*xDelt + x1
        !c            print*, k, xPoints(k), xDelt
        k = k + 1
     enddo
  endif

15 format(6(1x,1pd12.5))

  return
end subroutine GridMaker
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
subroutine CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,xNumPoints,u,xBounds,MatrixDim,NumChan,S)
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
              !c                  xS(ixp,ix) = 0.0d0
              S((NumChan-1)*xDim+ix,(NumChan-1)*xDim+ixp) = 0.0d0
              do kx = kxMin(ixp,ix),kxMax(ixp,ix)
                 xTempS = 0.0d0

                 do lx = 1,LegPoints
                    a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
                    b = a*u(lx,kx,ixp)
                    xTempS = xTempS + b
                 enddo
                 !c                     xS(ixp,ix) = xS(ixp,ix) +   xTempS
                 S((NumChan-1)*xDim+ix,(NumChan-1)*xDim+ixp) = S((NumChan-1)*xDim+ix,(NumChan-1)*xDim+ixp) &
                      + xTempS
              enddo
           enddo
        enddo
     enddo
  enddo


  deallocate(xIntScale)
  deallocate(kxMin,kxMax)

  return
end subroutine CalcOverlap
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcHamiltonian(RMin,RMax,mu,alpha,energy,Left,Right,Order,RPoints,LegPoints,xLeg,wLeg,RDim,& 
     RNumPoints,u,ux,uxx,RBounds,MatrixDim,NumChan,H,S,BigPot)

  implicit none
  double precision, external :: kdelta
  integer MatrixDim, RNumPoints,kx,lx,RDim,Order,ix,ixp,Left,Right,LegPoints,NumChan,mch,nch ! 

  double precision RPoints(RNumPoints),xLeg(LegPoints),wLeg(LegPoints),alpha,mu,energy,a,prefact ! 
  double precision u(LegPoints,RNumPoints,RDim), uxx(LegPoints,RNumPoints,RDim),&
       ux(LegPoints,RNumPoints,RDim),RMax,RMin ! 
  double precision x, xScaledZero,Znet,hbarc

  double precision H(MatrixDim,MatrixDim),S(MatrixDim,MatrixDim),ax,bx,potvalue ! 
  double precision, allocatable :: xIntScale(:),xIntPoints(:,:)
  double precision, allocatable ::  threshold(:) 
  integer, allocatable :: kxMin(:,:),kxMax(:,:)
  common hbarc
  double precision TempH, TempS, test, BigPot(NumChan,NumChan,LegPoints,RNumPoints)
  integer RBounds(RNumPoints+2*Order),i

  print*, 'NumChan=',NumChan

  allocate(xIntScale(RNumPoints))
  allocate(kxMin(RDim,RDim),kxMax(RDim,RDim)) ! 
  allocate(xIntPoints(LegPoints,RNumPoints))

  print*,'# NumChan = ', NumChan
  print*,'# MatrixDim = ',MatrixDim

  prefact=0.5d0/mu      

  !c      print*,'#seting up kxMin and kxMax'
  do ix = 1,RDim
     do ixp = 1,RDim
        kxMin(ixp,ix) = max(RBounds(ix),RBounds(ixp))
        kxMax(ixp,ix) = min(RBounds(ix+Order+1),RBounds(ixp+Order+1))-1 ! 
        !c            print*,ixp,ix,kxMin(ixp,ix),kxMax(ixp,ix)
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

!  call CheckPot(BigPot,RNumPoints,NumChan,LegPoints,xLeg,RPoints) ! 

  print*, '#Calculating Hamiltonian Matrix elements'

  do ix = 1,RDim
     do ixp = max(1,ix-Order),min(RDim,ix+Order)

        do mch = 1,NumChan
           do nch = 1, NumChan
              !c                  print*, 'mch, nch = ',mch, nch
              do kx = kxMin(ixp,ix),kxMax(ixp,ix)
                 TempH = 0.0d0
                 TempS = 0.0d0
                 do lx = 1,LegPoints
                    a = xIntPoints(lx,kx)*wLeg(lx)*xIntScale(kx)
                    TempS = TempS + a*kdelta(mch,nch)*u(lx,kx,ix)*u(lx,kx,ixp) ! 
                    TempH = TempH - a*prefact*kdelta(mch,nch)*u(lx,kx,ix)*uxx(lx,kx,ixp) ! 
                    TempH = TempH + a*u(lx,kx,ix)*BigPot(mch,nch,lx,kx)*u(lx,kx,ixp) ! sym piece
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

end subroutine CalcHamiltonian
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function kdelta(mch,nch)

  integer mch,nch
  if (mch.eq.nch) then
     kdelta = 1.0d0
  else 
     kdelta = 0.0d0
  endif
  return
end function kdelta
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

end subroutine Mydggev

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
end subroutine CompSqrMatInv
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CheckPot(Pot,RNumPoints,NumChan,LegPoints,xLeg,RPoints)
  implicit none
  integer RNumPoints,NumChan,LegPoints,mch,nch,kx,lx,iSector
  double precision Pot(NumChan,NumChan,LegPoints,RNumPoints),ax,bx,xLeg(LegPoints)
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
              write(888,*) x, Pot(mch,nch,lx,kx)
           enddo
        enddo
        write(888,*) ' '
     enddo
  enddo

  deallocate(xIntScale)
end subroutine CheckPot
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

end subroutine CheckBasisPhi
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
           write(file,*) x, u(lx,kx,ix)
        enddo
     enddo
     write(file,*) ' '
  enddo

  deallocate(xIntScale)
end subroutine CheckBasis
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function phirecon(R,beta,nch,evec,left,right,RDim,MatrixDim,RNumPoints,RPoints,order)
  implicit none
  double precision, external :: BasisPhi
  integer MatrixDim,RDim,nch,beta,i,RNumPoints,left,right,order
  double precision R,evec(MatrixDim,MatrixDim),RPoints(RNumPoints)
  phirecon = 0.0d0
  do i = 1,RDim
     phirecon = phirecon + evec((nch-1)*RDim+i,beta)*BasisPhi(R,left,right,order,RDim,RPoints,RNumPoints,0,i)
  enddo
  !c      if(R.ne.0.0d0) then
  !c         phirecon = phirecon/dsqrt(R) 
  !c      endif
  !c      print*, R, phirecon
  return
end function phirecon

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printmatrix(M,nr,nc,file)
  implicit none
  integer nr,nc,file,j,k
  double precision M(nr,nc)

  do j = 1,nr
     write(file,*) (M(j,k), k = 1,nc)
  enddo
end subroutine printmatrix
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine setup_potential_matrix(Pot,mu,NumChan,LegPoints,RNumPoints,RPoints,xLeg)

  implicit none
  integer NumChan,LegPoints,RNumPoints,lx,kx,mch,nch
  double precision Vsquare, r0, mu, D0, Cij,Delta
  double precision, external :: kdelta
  double precision, allocatable :: xIntScale(:),xIntPoints(:,:)
  double precision Pot(NumChan,NumChan,LegPoints,RNumPoints),ax,bx,xScaledZero,xLeg(LegPoints) ! 
  double precision RPoints(RNumPoints), ethresh(NumChan), D1,vmat(NumChan,NumChan)! 

  allocate(xIntPoints(LegPoints,RNumPoints))
  allocate(xIntScale(RNumPoints))

  ! Some square well parameters:
  r0 = 1.d0 !width of well
  D1 = 10d0
  Cij = 3d0
  Delta = 50d0
  call makeVTridiag(NumChan,1,D1,Delta,Cij,vmat,ethresh)
  
  do mch = 1, NumChan
     do nch = 1, NumChan
        !c     print*,'kx,  lx, ax,              bx,            xIntScale,     xScaledZero,     x' ! 
        do kx = 1,RNumPoints-1
           ax = RPoints(kx)
           bx = RPoints(kx+1)
           !c     print*, 'RPoints(1)=',RPoints(kx), 'RPoints(2)=',RPoints(kx+1)
           xIntScale(kx) = 0.5d0*(bx-ax)
           xScaledZero = 0.5d0*(bx+ax)
           do lx = 1,LegPoints
              Pot(mch,nch,lx,kx)=0.0d0
              xIntPoints(lx,kx) = xIntScale(kx)*xLeg(lx)+xScaledZero
              
              !c     print*, kx, lx, ax, bx, xIntScale, xScaledZero, xIntPoints(lx,kx) !
              
              if(xIntPoints(lx,kx).lt.r0) then
                 Pot(mch,nch,lx,kx) = vmat(mch,nch)
              else
                 Pot(mch,nch,lx,kx) = kdelta(mch,nch)*ethresh(nch)
              endif

           enddo
        enddo

     enddo
  enddo
  deallocate(xIntScale,xIntPoints)
end subroutine setup_potential_matrix
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine makeVTridiag(NumChan,NumOpen,D,Delta,Cij,Vmat,ethresh)
  use modb
  implicit none
  integer NumChan,NumOpen,i,j
  double precision Vmat(NumChan,NumChan),eoffset(NumChan - NumOpen),ethresh(NumChan)
  double precision D,Delta,Cij
  
  Vmat = 0.d0 !initialize the potential matrix to zero
  do i = 1, NumChan
     ethresh(i) = Delta
     do j = 1,NumChan
        if(i.eq.j)    Vmat(i,j) = -D
        if((i.eq.j+1).or.(i.eq.j-1)) Vmat(i,j) = Cij
     enddo
  enddo
  ethresh(NumChan) = 0d0
  
end subroutine makeVTridiag
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE deigsrt(d,v,n,np)
  INTEGER n,np
  double precision d(np),v(np,np)
  INTEGER i,j,k
  double precision p
  do i=1,n-1
     k=i
     p=d(i)
     do j=i+1,n
        if(d(j).ge.p)then
           k=j
           p=d(j)
        endif
     enddo
     if(k.ne.i)then
        d(k)=d(i)
        d(i)=p
        do j=1,n
           p=v(j,i)
           v(j,i)=v(j,k)
           v(j,k)=p
        enddo
     endif
  enddo
  return
END SUBROUTINE deigsrt
        !C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
