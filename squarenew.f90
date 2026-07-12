!  magnetic Feshbach resonances with square wells
!  compile with make -f MagneticF.mak  (see the make file for dependencies)
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
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Module potmod
  implicit none
  real*8 mass,Rmatch
  integer L
end Module potmod
!****************************************************************************************************
module DataStructures
  implicit none
  !-----------------
  type SolutionData
     real*8 ascat,alsqr,zclose,BSdet,bgphase,phase,sigma,phaseopi,timedelay
     real*8 abg,dtaudE,d2tau
  end type SolutionData
  !-----------------
  type PotentialData
     real*8, allocatable :: V(:,:),Vcc(:,:),Lambda(:),LambdaC(:),UC(:,:),U(:,:)
     real*8, allocatable :: ethresh(:),eoffset(:)
     real*8  BCE, BCO
  end type PotentialData
  !-----------------
  type ResonanceData
     integer iel,ier ! array position for bracketing energies
     real*8 el,er, eres ! the left and right energy values bracketing the resonance eres
     real*8 al,ar  !the scattering length at the left and right points
     real*8 zclose,deltabg, width, fanoq, timedelay, bgampsqr
  end type ResonanceData
  !--------------------
  type BinData
     !     real*8, allocatable :: FitBins(:), Bins(:)
     integer, allocatable :: BinC(:)
     real*8, allocatable :: pbrody(:)
     real*8 binsum, wparam(1),stdev
  end type BinData
  !--------------------
  type StatisticalData
     real*8, allocatable :: ls(:)
     real*8 avgd, avgwidth
     integer, allocatable :: BinC(:)
  end type StatisticalData
  !-----------------
  contains
    subroutine AllocatePotData(Vdata,NumChan,NumClosed)
      !allocates memory for the potential matrix and corresponding transformation matrices
      implicit none
      type(PotentialData) :: Vdata
      integer NumChan,NumClosed

      allocate(Vdata%V(NumChan,NumChan),Vdata%Vcc(NumClosed,NumClosed),Vdata%ethresh(NumChan),Vdata%eoffset(NumClosed))
      allocate(Vdata%Lambda(NumChan),Vdata%LambdaC(NumClosed),Vdata%U(NumChan,NumChan),Vdata%UC(NumClosed,NumClosed))

      Vdata%V=0d0
      Vdata%Vcc=0d0
      Vdata%ethresh=0d0
      Vdata%eoffset=0d0
      Vdata%Lambda=0d0
      Vdata%LambdaC=0d0
      Vdata%U=0d0
      Vdata%UC=0d0
      
    end subroutine AllocatePotData
  end module DataStructures
  !****************************************************************************************************
Program main
  use modb
  use potmod
  use DataStructures
  implicit none
  external SBrody
  double precision de, dB,blo,bhi,elo, ehi, xx1, time1, time2, tloop,ei,ef
  double precision BCE,BCO,cdlength,gocmin,gocmax,dgoc,goc,gccmin,gccmax,vopen,Boffset,diffE,atl,atr,abgr,abgl
  double precision DBClosedWindow,DBOpenRes,DVClosedWindow,goctemp,dBins,te1,te2,timedelayold
  double precision, allocatable :: Bins(:),FitBins(:), ensavgd(:,:),ensavgwidth(:,:)
  integer i,j,ne,nb,iB, NumChan, NumClosed, NumOpen, NumEven, ie,ifile,nens,ncoup,iens,icoup,idum
  integer numres,NumResMax,ires,NumBins,NumIter,WriteFlag
  integer, allocatable :: nres(:,:,:)
  character*11 filename
  double precision, allocatable :: energy(:), Bfield(:), gocscale(:),gccscale(:)
  type(SolutionData) :: sol
  type(PotentialData) :: V0,Vshifted
  type(PotentialData), allocatable :: VArray(:,:)
  type(SolutionData), allocatable :: solarray(:,:,:,:)
  type(ResonanceData), allocatable :: Res(:,:,:,:)
  type(StatisticalData), allocatable :: SD(:,:,:), SDtot(:,:)
  type(BinData), allocatable :: BC(:,:)
  type(ResonanceData) :: ResTemp
!  allocate(data1%ascat(10))
  
  call CPU_TIME(time1)
  
  !  open(unit=910,file="background.dat")
  !open(unit=920,file="result.dat")
  !open(unit=930,file="ampn.dat")

  ! Read in some info from the input file
  open(unit=890,file="squarenew.par")
  read(890,*)
  read(890,*) elo,ehi,ne,blo,bhi,nb,nens,ncoup  !Energy range and number of energy steps
  read(890,*)
  read(890,*) NumChan,NumEven,mass,Rmatch,L, NumBins  !number of channels, mass, matching radius, partial wave
  read(890,*)
  read(890,*) BCE, BCO, gocmin,gocmax, gccmin,gccmax, cdlength, vopen !threshold spacing, open-closed coupling, decay length of Vmat, number of gcc's to try
  read(890,*)
  read(890,*) DVClosedWindow, DBClosedWindow,DBOpenRes, Boffset
  read(890,*)
  read(890,*) WriteFlag
  write(6,*) "Rmatch = ", Rmatch
  NumOpen = 1
  NumClosed = NumChan-NumOpen ! number of closed channels
!  NumEven = NumClosed/2
  allocate(energy(ne))
  allocate(gocscale(ncoup),gccscale(ncoup))
  allocate(nres(nens,ncoup,nb))

  if(ncoup.gt.1) then
     dgoc = (gocmax - gocmin)/dble(ncoup-1)
     do icoup = 1, ncoup
        gocscale(icoup) = gocmin + dgoc*(icoup-1)
        gccscale(icoup) = gccmin + (gccmax - gccmin)/(dble(ncoup-1))*(icoup-1)
     enddo
  else
     ncoup = 1
     gocscale(1) = gocmax
     gccscale(1) = gccmax
  endif
  if(ne.gt.1) then
     ei = log(elo)
     ef = log(ehi)
     de=(ehi-elo)/(ne-1) !for linear grid
     !de=(ef-ei)/(ne-1) !for log grid
     diffE = de/100d0
    
     do ie = 1,ne
       ! energy(ie) = exp(ei+de*(ie-1)) ! log grid
        energy(ie) = elo+de*(ie-1) !linear grid
!        write(6,*) ie, energy(ie)
     end do
  else
     de=0d0
     diffE = 1d-5
     ne = 1
     energy(1) = elo
!     write(6,*) "1", energy(1)
  endif
  allocate(Bfield(nb))
  if(nb.gt.1) then
     dB = (bhi-blo)/(nb-1)
     do iB = 1, nB
        Bfield(iB) = blo + dB*(iB-1)
 !       write(6,*) iB, Bfield(iB)
     enddo
  else
     db = 0d0
     nb = 1
     bfield(1) = blo
!     write(6,*) "1", bfield(1)
  end if
  NumResMax = 200
  allocate(Res(nens,ncoup,nb,NumResMax))   ! allocate memory for the resonance data

  !----------------------------------------------------------------------------------------------------
  ! allocate the memory for the array of potentials
  allocate(VArray(nens,ncoup))
  allocate(solarray(nens,ncoup,nb,ne))
  do iens = 1, nens
     do icoup = 1, ncoup
        call AllocatePotData(VArray(iens,icoup),NumChan,NumClosed)
!        call printmatrix(VArray(iens,icoup)%V,NumChan,NumChan,6)
     enddo
  enddo
  !----------------------------------------------------------------------------------------------------
  ! set the initial open channel depth. v0 = (n-1/2)**2 * pi**2 gives the nth bound state at threshold,
  ! so we'll make the well a bit deeper than than needed for one state.
!  vopen = 100d0!((1.d0 - 0.5d0)*pi)**2 + 10.d0  
  
  call AllocatePotData(V0,NumChan,NumClosed)
  call AllocatePotData(Vshifted,NumChan,NumClosed)
  call ranset_()
  write(6,*) "DVClosedWindow = ", DVClosedWindow
  write(6,*) "initialize the positions of the uncoupled resonances for each member of the ensemble by setting eoffset..."
  do iens = 1, nens
     call rgnf_lux(V0%eoffset,NumClosed)
     call sort(NumClosed,V0%eoffset)
     !V0%eoffset = log(elo) + V0%eoffset*log(ehi) !for log scale
     !V0%eoffset = exp(V0%eoffset) ! for log scale
     V0%eoffset = V0%eoffset*DVClosedWindow
     call makeVNew(NumChan,NumOpen,NumEven,vopen,BCE,BCO,V0%eoffset,&
          V0%V,V0%ethresh)
     do icoup = 1, ncoup
        Varray(iens,icoup)=V0
        do i = 1, NumChan-1
           do j = i+1,NumChan-1
              Varray(iens,icoup)%V(i,j) = gccscale(icoup)*V0%V(i,j)
              Varray(iens,icoup)%V(j,i) = gccscale(icoup)*V0%V(j,i)
           enddo
           Varray(iens,icoup)%V(i,NumChan) = gocscale(icoup)*V0%V(i,NumChan)
           Varray(iens,icoup)%V(NumChan,i) = gocscale(icoup)*V0%V(NumChan,i)
        enddo
!        call rgnf_lux(Varray(iens,icoup)%eoffset,NumClosed)
        !do i = 1, NumClosed
           !Varray(iens,icoup)%eoffset(i) = (dble(i))/dble(NumClosed)
           !Varray(iens,icoup)%eoffset(i) = (dble(i)*3.18d0)
        !enddo
!        call sort(NumClosed,Varray(iens,icoup)%eoffset)
!        Varray(iens,icoup)%eoffset = Varray(iens,icoup)%eoffset*DVClosedWindow
!        call printmatrix(Varray(iens,icoup)%eoffset,1,NumClosed,6)
        
!        call makeVOld(NumChan,NumOpen,NumEven,vopen,BCE,BCO,gocscale(icoup),gccscale(icoup),VArray(iens,icoup)%eoffset,&
!             VArray(iens,icoup)%V,VArray(iens,icoup)%ethresh)
!        call makeVTridiag(NumChan,NumOpen,1000d0,1000d0,300d0,VArray(iens,icoup)%V,VArray(iens,icoup)%ethresh)
! The general GOE-ensemble path above (scaling V0%V by gccscale/gocscale) is what produced the
! published N=41 Brody-statistics figures (Sec. V of PRA 98, 062703).
! To instead reproduce the Sec. IV pedagogical N=3 tridiagonal example (Figs. 2-3), uncomment:
!        call makeVTridiag(NumChan,NumOpen,50d0,200d0,5d0,VArray(iens,icoup)%V,VArray(iens,icoup)%ethresh)
! and loosen the resonance-selection threshold below from 4.0d0 to 0.0d0 (the N=3 case has fewer,
! less sharp resonances than the N=41 ensemble).
     enddo
  enddo

  ! Calculate the solutions

  do iens = 1, nens
     do icoup = 1, ncoup
        call CPU_TIME(te1)
        do iB = 1,nB
           call ShiftV(Bfield(iB),NumChan,NumOpen,DVClosedWindow,DBClosedWindow,DBOpenRes,VArray(iens,icoup),Vshifted)
           do ie=1,ne     ! ENERGY LOOP
              call EnergyPoint(rmatch,L,mass,energy(ie),NumChan,NumClosed,Vshifted,solarray(iens,icoup,iB,ie))
!              write(6,'(4(A10,I8),5(A15,ES16.5))') 'iens = ',iens,'icoup = ',icoup,'iB =',iB,&
!                   'ie = ',ie,'goc = ',gocscale(icoup),'gcc = ',gccscale(icoup),'B=',Bfield(iB),'E =',energy(ie)
           enddo ! ENERGY LOOP
        enddo ! B FIELD LOOP
        call CPU_TIME(te2)
        write(6,'(A10,I10,A10,I10,A10,e14.5,A10,e14.5,A10,e14.5,A1)') "iens =",iens, "icoup = ",&
             icoup, "goc = ",gocscale(icoup), "gcc = ", gccscale(icoup), "time = ", (te2-te1)/60d0,"m"
     enddo ! COUPLING LOOP
  enddo ! ENSEMBLE LOOP
  write(6,*) "Calling fixphasejumps..."
  call fixphasejumps(nens,ncoup,nB,ne,solarray,energy,Bfield,WriteFlag) ! this fixes the phase jumps and prints to a file.
  
  
  ! go back and do a search for the peaks of the time delay

  write(800,'(4A10,100A14)') "iB","iens","icoup","numres", "eres","zclose","width","fanoq","bgampsqr"
  do iB = 1,nb
     do iens = 1, nens
        do icoup = 1,ncoup
           call ShiftV(Bfield(iB),NumChan,NumOpen,DVClosedWindow,DBClosedWindow,DBOpenRes,VArray(iens,icoup),Vshifted)
           write(6,'(A50,3I10)',advance='no') "Finding resonance positions for :",iB,iens,icoup
           !write(6,'(A50,3I10)') "Finding resonance positions for:",iB,iens,icoup
           numres = 0
           do ie = 2, ne
!!$              if( (solarray(iens,icoup,iB,ie-1)%timedelay.gt.solarray(iens,icoup,iB,ie-2)%timedelay)&
!!$                   .and.(solarray(iens,icoup,iB,ie-1)%timedelay.gt.solarray(iens,icoup,iB,ie)%timedelay).and.&
!!$                   (log(solarray(iens,icoup,iB,ie-1)%timedelay).gt.0d0) ) then
              
              if( ((solarray(iens,icoup,iB,ie-1)%dtaude*solarray(iens,icoup,iB,ie)%dtaude).lt.0d0)&
                   .and. ((solarray(iens,icoup,iB,ie)%dtaude-solarray(iens,icoup,iB,ie-1)%dtaude).lt.0d0)) then

                 ResTemp%iel = ie-1
                 ResTemp%ier = ie
                 ResTemp%el = energy(ie-1)
                 ResTemp%er = energy(ie)
                 ResTemp%eres =  energy(ie-1)
                 ResTemp%timedelay = solarray(iens,icoup,iB,ie-1)%timedelay
                 ResTemp%zclose = solarray(iens,icoup,iB,ie-1)%zclose
                 ResTemp%width = 4d0/solarray(iens,icoup,iB,ie-1)%timedelay
                 ResTemp%bgampsqr = sin(solarray(iens,icoup,iB,ie-1)%bgphase)**2
                 ResTemp%fanoq = -1d0/tan(solarray(iens,icoup,iB,ie-1)%bgphase)

                 ! go and refine the resonance positions and values
                 call FindResPos(iens,Bfield(ib),mass,rmatch,L,NumChan,NumClosed,&
                      ResTemp,Vshifted,solarray(iens,icoup,iB,ie-1))
                 ! select only real resonances (peaks in the time delay, not minima, and not weak maxima)
                 ! threshold 4.0d0 is tuned for the N=41 GOE ensemble; use 0.0d0 for the N=3
                 ! tridiagonal pedagogical example (see comment near makeVTridiag above)
                 if (log(ResTemp%timedelay).gt.4.0d0) then
                    numres = numres + 1
                    Res(iens,icoup,iB,numres)=ResTemp
                    write(800,14) iB,iens,icoup,numres, Res(iens,icoup,iB,numres)%eres, Res(iens,icoup,iB,numres)%zclose,&
                         Res(iens,icoup,iB,numres)%width,Res(iens,icoup,iB,numres)%fanoq,Res(iens,icoup,iB,numres)%bgampsqr
!                    write(6,*) "NumRes = ", numres, Res(iens,icoup,iB,numres)%eres
                 endif

                 nres(iens,icoup,iB) = numres
!                 write(6,'(A10,I10)')"numres = ", numres
              endif
           enddo
           write(6,'(A10,I10)')"numres = ", numres
        enddo
     enddo
  enddo
14 format(4I10,100e14.5)
  ! Collect statistical data for each member of the ensemble at each coupling and Bfield
  dbins = 5.d0/NumBins
  allocate(Bins(NumBins+1),Fitbins(NumBins))

!  write(6,*) "---------------"
!  write(6,*) "Bins ..."
  do i = 1,NumBins+1
     Bins(i) = 0.d0 + (i-1)*dbins
!     write(6,*) i,Bins(i)
  enddo
  write(6,*) "---------------"
  do i = 1,NumBins
     FitBins(i) = 0.5d0*(Bins(i) + Bins(i+1))
     !FitBins(i) = (Bins(i+1))
  enddo

  allocate(SD(nens,ncoup,nb))
  allocate(ensavgd(ncoup,nb),ensavgwidth(ncoup,nb))
  ensavgd = 0d0
  ensavgwidth = 0d0
  do iB = 1, nB
     do icoup = 1, ncoup
        do iens = 1, nens
           numres=nres(iens,icoup,iB)
           allocate(SD(iens,icoup,ib)%ls(numres-1))
           allocate(SD(iens,icoup,ib)%BinC(NumBins))
           do ires = 1, numres - 1
              SD(iens,icoup,iB)%ls(ires) = Res(iens,icoup,iB,ires+1)%eres - Res(iens,icoup,iB,ires)%eres
!              write(6,*) iens,icoup,iB,ires,SD(iens,icoup,iB)%ls(ires)
           enddo
           SD(iens,icoup,iB)%avgd=sum(SD(iens,icoup,iB)%ls(1:numres-1))/dble(numres-1) ! average level spacing for this member of the ensemble
           SD(iens,icoup,iB)%avgwidth = sum(Res(iens,icoup,iB,1:numres)%width)/dble(numres) !average width for this member of the ensemble
           write(6,*) iB, icoup, iens, "avgd =",SD(iens,icoup,iB)%avgd, "avgwidth = ", SD(iens,icoup,iB)%avgwidth
           if(SD(iens,icoup,iB)%avgwidth.lt.0d0) then
              write(6,*) "note width negative for icoup,iens = ",icoup,iens, (Res(iens,icoup,iB,ires)%width, ires=1,numres)
              write(6,*) "press enter to continue..."
              read(5,*)
           endif
           ensavgd(icoup,iB) = ensavgd(icoup,iB) + SD(iens,icoup,iB)%avgd
           ensavgwidth(icoup,iB) = ensavgwidth(icoup,iB) + SD(iens,icoup,iB)%avgwidth
!           SD(iens,icoup,iB)%ls(:) = SD(iens,icoup,iB)%ls(:)/SD(iens,icoup,iB)%avgd ! scale the level spacings by the average
!           call bincount(NumBins,numres-1,Bins,SD(iens,icoup,iB)%ls(:),SD(iens,icoup,iB)%BinC)
!           write(200,'(100I5)') iB, icoup, iens, SD(iens,icoup,iB)%BinC
!           write(6,'(I4,I4,I4,A1,100I5)') iB, icoup, iens, "|", SD(iens,icoup,iB)%BinC
        enddo
        ensavgd(icoup,iB) = ensavgd(icoup,iB)/dble(nens)
        ensavgwidth(icoup,iB) = ensavgwidth(icoup,iB)/dble(nens)
        write(6,*) iB, icoup,"ensavgd =",ensavgd(icoup,iB), "ensavgwidth = ", ensavgwidth(icoup,iB)
     enddo
  enddo

  allocate(BC(ncoup,nb))
  NumIter = 10
  if(NumChan.gt.2) then
     write(700,'(2A10,100A14)') "iB","icoup","<D>","goc/<D>","gcc/<D>","<Gamma>/<D>",&
          "brody-w","stdev-w"
     write(500,'(2A10,6A14,A50)') "iB","icoup","<D>","goc/<D>","gcc/<D>","<Gamma>/<D>",&
          "brody-w","stdev-w","level-spacing-bin-counts."
     write(400,'(2A10,6A14,A50)') "iB","icoup","<D>","goc/<D>","gcc/<D>","<Gamma>/<D>",&
          "brody-w","stdev-w","Brody-dist.-probabilities"
     do iB = 1,nB
        do icoup = 1,ncoup
           allocate(BC(icoup,ib)%BinC(NumBins))
           allocate(BC(icoup,ib)%pbrody(NumBins))
           BC(icoup,iB)%BinC=0d0
           do iens = 1, nens
              numres=nres(iens,icoup,iB)
              SD(iens,icoup,iB)%ls(:) = SD(iens,icoup,iB)%ls(:)/ensavgd(icoup,iB) ! scale the level spacings by the ensemble average
              call bincount(NumBins,numres-1,Bins,SD(iens,icoup,iB)%ls(:),SD(iens,icoup,iB)%BinC)
              !           write(6,*) "Bincounts = ",SD(iens,icoup,iB)%BinC
              BC(icoup,iB)%BinC(:) = BC(icoup,iB)%BinC + SD(iens,icoup,iB)%BinC(:)
           enddo
           BC(icoup,iB)%binsum = dble(sum(BC(icoup,iB)%BinC))
           BC(icoup,iB)%pbrody = dble(BC(icoup,iB)%BinC)/BC(icoup,iB)%binsum/dbins
           !        write(6,*) "Total BC = ", BC(icoup,iB)%BinC
           !        write(6,*) "Brody Dist = ",BC(icoup,iB)%pbrody
           BC(icoup,iB)%wparam(1) = 0.5d0
           call Param_LS(1,NumIter,NumBins,BC(icoup,iB)%stdev,0.001d0,0.1d0,BC(icoup,iB)%wparam,Fitbins,BC(icoup,iB)%pbrody,SBrody)
           write(400,11) iB, icoup,ensavgd(icoup,iB), gocscale(icoup)/ensavgd(icoup,iB), gccscale(icoup)/ensavgd(icoup,iB), &
                ensavgwidth(icoup,iB)/ensavgd(icoup,iB),BC(icoup,iB)%wparam(1),BC(icoup,iB)%stdev, &
                (BC(icoup,iB)%BinC(i), i = 1, NumBins)
           write(500,12) iB, icoup, ensavgd(icoup,iB), gocscale(icoup)/ensavgd(icoup,iB), gccscale(icoup)/ensavgd(icoup,iB), &
                ensavgwidth(icoup,iB)/ensavgd(icoup,iB),BC(icoup,iB)%wparam(1),BC(icoup,iB)%stdev, &
                (BC(icoup,iB)%pbrody(i), i = 1, NumBins)
           write(700,12) iB, icoup,ensavgd(icoup,iB), gocscale(icoup)/ensavgd(icoup,iB), &
                gccscale(icoup)/ensavgd(icoup,iB), &
                !gccscale(icoup), &
                ensavgwidth(icoup,iB)/ensavgd(icoup,iB), BC(icoup,iB)%wparam(1), BC(icoup,iB)%stdev
        enddo
     enddo
  endif
11 format(I10,I10, e14.5,e14.5,e14.5,e14.5,e14.5,e14.5,1000I5)
12 format(I10,I10, 100e14.5)
        
  call CPU_TIME(time2)
  xx1=(time2-time1)/60d0
  write(6,*)'Done, total time for this calculation:',xx1,"m"
  !c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !        call MPI_FINALIZE(ierr)
10 format(100e24.14)
End Program main

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fixphasejumps(nens,ncoup,nB,ne,solarray,energy,Bfield,WriteFlag)
  use DataStructures
  use modb
  implicit none
  integer nens, ncoup, nB, ne,iens,icoup,ifile,iB,ie,idum,WriteFlag
  double precision energy(ne),timedelayold,Bfield(nB),dtdeold
  type(SolutionData) solarray(nens,ncoup,nB,ne)
  
  do iens = 1, nens
     do icoup = 1, ncoup
        ifile = 1000*iens + icoup
        if(WriteFlag.eq.1) open(ifile)
        do iB = 1,nB
           write(6,'(A50,3I10)') "Fixing the phase and finding finding time delay :",iB,iens,icoup
           do ie=2,ne     ! ENERGY LOOP
              ! first go through and fix the phase jump due to the background reaching -pi/2
              if(( solarray(iens,icoup,iB,ie)%phaseopi - solarray(iens,icoup,iB,ie-1)%phaseopi).gt.0.05d0) then
!                 write(6,*) "Fixing phase at:", energy(ie)
                 do idum = ie, ne
                    solarray(iens,icoup,iB,idum)%phaseopi = solarray(iens,icoup,iB,idum)%phaseopi - 1.d0
                    solarray(iens,icoup,iB,idum)%phase = solarray(iens,icoup,iB,idum)%phase - pi 
                 enddo
              endif
              ! This searches for discontinuities in the phase...
              if ( (solarray(iens,icoup,iB,ie)%phaseopi - solarray(iens,icoup,iB,ie-1)%phaseopi ) .lt. -0.05d0) then
                 do idum = ie, ne
                    solarray(iens,icoup,iB,idum)%phaseopi = solarray(iens,icoup,iB,idum)%phaseopi + 1.d0
                    solarray(iens,icoup,iB,idum)%phase = solarray(iens,icoup,iB,idum)%phase + pi
                 enddo
              endif
              !calculate the time delay
!              timedelayold = 2d0*(solarray(iens,icoup,iB,ie)%phase - solarray(iens,icoup,iB,ie-1)%phase)&
!                   /(energy(ie)-energy(ie-1))
!              dtdeold = (solarray(iens,icoup,iB,ie)%timedelay - solarray(iens,icoup,iB,ie-1)%timedelay)&
              !                   /(energy(ie)-energy(ie-1))
              if(WriteFlag.eq.1) write(ifile,13) iens,energy(ie-1),Bfield(iB),solarray(iens,icoup,iB,ie-1)

           enddo ! ENERGY LOOP
        enddo  ! B FIELD LOOP
        if(WriteFlag.eq.1) write(ifile,*)
        if(WriteFlag.eq.1) flush(ifile)
        if(WriteFlag.eq.1) close(ifile)
     enddo ! COUPLING LOOP
  enddo ! ENSEMBLE LOOP
13 format(I8,100e14.5)

end subroutine fixphasejumps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bincount(NumBins,NumPoints,Bins, A, Counts)
  integer NumBins, NumPoints,ip,ib,Counts(NumBins)
  real*8 Bins(NumBins+1), A(NumPoints)
  Counts=0
  do ip = 1,NumPoints
     do ib = 1,NumBins
        if( (A(ip).lt.Bins(ib+1)) .and. (A(ip).gt.Bins(ib)) ) Counts(ib) = Counts(ib)+1
     enddo
  enddo
  
end subroutine bincount
!****************************************************************************************************
subroutine FixedBEnergyRun(EPoints,NumEPoints,NumChan,NumOpen,NumEven,Bfield)
  use modb
  integer NumChan,NumOpen,NumEven,NumEPoints, i,j
  double precision EPoints(NumEPoints),Bfield
  
  
end subroutine FixedBEnergyRun
!****************************************************************************************************
subroutine EnergyPoint(rmatch,L,mass,energy,n,nc,V,sol)
  !------------------------------------------ 
  use modb
  use DataStructures
  implicit none
  real*8 energy,mass,rmatch 
  integer i,ip,j,k,nc,n,ipiv(2*n+1),n2,n21,L
  real*8 Amat(2*n,n),Bmat(2*n,n+1),d1amat(2*n,n),d1bmat(2*n,n+1)
  real*8 d2amat(2*n,n),d2bmat(2*n,n+1)
  real*8 d3amat(2*n,n),d3bmat(2*n,n+1),d3abmat(2*n+1,2*n+1)
  real*8 PhiC(nc,nc),dPhiC(nc,nc)
  real*8 Gmat(nc,nc),ABmat(2*n+1,2*n+1)
  real*8 d1abmat(2*n+1,2*n+1),ABmatTEMP(2*n+1,2*n+1)
  real*8 d2abmat(2*n+1,2*n+1)
  real*8 fc(nc),dfc(nc), dfcde(nc), ddfcde(nc),d2fc(nc),d2dfc(nc),d3fc(nc),d3dfc(nc)
  real*8 bvec(2*n+1),bcoef(n),aclose(nc),amp(nc),cd,sd,BSmat(nc,nc),cvec(2*n+1),dsde,dcde
  real*8 dvec(2*n+1),d3vec(2*n+1)
  real*8 d2s, d2c,d3c,d3s
  type(PotentialData), intent(in) :: V  ! input
  type(SolutionData) :: sol ! output
  
  n2=2*n
  n21=2*n+1

  Amat=0d0
  d2amat=0d0; d3amat=0d0
  d2bmat=0d0
  d3bmat=0d0
!  write(6,*) "make Amat, Bmat..."
  call makeAmat(energy,V%Lambda,V%U,Amat,d1amat,d2amat,d3amat,n,V%ethresh) !build the interior block structure
  call makeBmat(energy,Bmat,d1bmat,d2bmat,d3bmat,Gmat,nc,n,V%ethresh) !build the exterior block structure

  ABmat=0d0
  bvec=0d0
  d1abmat=0d0
  cvec=0d0
  d2abmat=0d0
  dvec=0d0
  d3abmat=0d0
  d3vec=0d0
  !put it all together into a big AB matrix
  do i=1,n2
     do ip=1,n+1
        if(ip.le.n) then
           ABmat(i,ip)=Amat(i,ip)
           d1abmat(i,ip)=d1amat(i,ip)
           d2abmat(i,ip)=d2amat(i,ip)
           d3abmat(i,ip)=d3amat(i,ip)
        endif
        ABmat(i,ip+n)=(-1d0)*Bmat(i,ip)
        d1abmat(i,ip+n)=(-1d0)*d1bmat(i,ip)
        d2abmat(i,ip+n)=(-1d0)*d2bmat(i,ip)
        d3abmat(i,ip+n)=(-1d0)*d3bmat(i,ip)
     enddo
  enddo
  ABmat(n21,n+1)=1d0; bvec(n21)=1d0 ! pick aclose(1)=1
  d1abmat(n21,n+1)=0d0
  d2abmat(n21,n+1)=0d0
!!$  ABmat(n21,n2)=1d0; bvec(n21)=1d0 ! pick cd=1
!!$  d1abmat(n21,n2)=0d0
  ABmatTEMP = ABmat
  call DGESV(N21,1,ABmatTEMP,N21,IPIV,bvec,N21,i)  !returns the solution vector in bvec (bvec overwritten)
  cvec = -matmul(d1abmat,bvec)

  ABmatTEMP = ABmat
  call DGESV(N21,1,ABmatTEMP,N21,IPIV,cvec,N21,i)  !returns the solution vector in cvec (cvec overwritten)

  ABmatTEMP = ABmat
  dvec = - 2d0*matmul(d1abmat,cvec) - matmul(d2abmat,bvec) 
!  call printmatrix(d2abmat,n21,n21,6)
  call DGESV(N21,1,ABmatTEMP,N21,IPIV,dvec,N21,i)  !returns the solution vector in dvec (dvec overwritten)

  ABmatTEMP = ABmat
  d3vec = -3d0*(matmul(d1abmat,dvec) + matmul(d2abmat,cvec)) - matmul(d3abmat,bvec)
  call DGESV(N21,1,ABmatTEMP,N21,IPIV,d3vec,N21,i)  !returns the solution vector in dvec (dvec overwritten)
  
  do i=1,n
     bcoef(i)=bvec(i) !interior coefficients
  enddo
  do i=1,nc
     aclose(i)=bvec(n+i) !exterior coefficients
  enddo
  
  cd = bvec(n2) 
  sd = bvec(n21) 
  
  sol%phase = atan(sd/cd) !delta
  sol%phaseopi=sol%phase/pi
  sol%alsqr = sin(sol%phase)**2
  sol%sigma = 2d0*pi*sin(sol%phase)**2/(mass*energy) !partial wave cross section
  sol%ascat = -tan(sol%phase)/sqrt(2d0*mass*energy)

  dcde = cvec(n2)
  dsde = cvec(n21)

  d2c = dvec(n2) 
  d2s = dvec(n21)

  d3c = d3vec(n2)
  d3s = d3vec(n21)
  
  sol%timedelay = 2d0*cos(sol%phase)**2*( -cd**(-2d0) * dcde * sd + dsde/cd)
  sol%dtaude = 2d0*(cd**2*d2s*cos(sol%phase)**2 - 2*cd*dcde*dsde*cos(sol%phase)**2 - &
       cd*d2c*sd*cos(sol%phase)**2 + 2*dcde**2*sd*cos(sol%phase)**2 - &
       (cd**3*sol%timedelay**2*tan(sol%phase))/2.d0)/cd**3
  sol%d2tau =  ((-3*cd**4*sol%timedelay**3)/4. + &
       (cd**4*sol%timedelay**3*Cos(sol%phase)**2)/2. + cd**3*d3s*Cos(sol%phase)**4 - &
       3*cd**2*d2s*dcde*Cos(sol%phase)**4 - &
       3*cd**2*d2c*dsde*Cos(sol%phase)**4 + &
       6*cd*dcde**2*dsde*Cos(sol%phase)**4 - cd**2*d3c*sd*Cos(sol%phase)**4 + &
       6*cd*d2c*dcde*sd*Cos(sol%phase)**4 - 6*dcde**3*sd*Cos(sol%phase)**4 - &
       (3*cd**4*sol%dtaude*sol%timedelay*Cos(sol%phase)*Sin(sol%phase))/2.)/(cos(sol%phase)**2*cd**4)


  !--------------------------------------------------------
  ! comment/uncomment for zclose
    call makePSI(mass,L,energy,V%Lambda,V%U,amp,aclose,bcoef,cd,sd,n,V%ethresh)  !Comment out if you want to skip the lengthy zclose calculation
     sol%zclose=sum(amp) !sum up the closed channel amplitudes
  !--------------------------------------------------------
  call background(energy,sol%bgphase,n,V%V) !this just performs a single-channel calcualtion on the open channel (neglecting all closed channels.)
  sol%abg = -tan(sol%bgphase)/sqrt(2d0*mass*energy)  
 
  
  !--------------------------------------
  !Closed Channel Sector...
!  write(6,*) "calling makeF for closed channel sector...energy =", energy
  call makeF(nc,rmatch,0,mass,energy,V%LambdaC,fc,dfc,dfcdE,ddfcdE,d2fc,d2dfc,d3fc,d3dfc) !make the diagonal f and df matrices (L=0 inside well)
  do k = 1, nc
     PhiC(1:nc,k) = V%UC(1:nc,k)*fc(k)!matmul(V%UC,PhiC)
     dPhiC(1:nc,k) = V%UC(1:nc,k)*dfc(k) !dPhiC = matmul(V%UC,dPhiC)
  enddo
  BSMat = dPhiC - matmul(Gmat,PhiC)
!!$  write(13,*) "Energy = ", Energy, "dPhiC = "
!!$  call printmatrix(dPhiC,nc,nc,13)
!!$  
!!$  write(13,*) "Energy = ", Energy, "Gmat = "
!!$  call printmatrix(Gmat,nc,nc,13)
!!$  
!!$  write(13,*) "Energy = ", Energy, "PhiC = "
!!$  call printmatrix(PhiC,nc,nc,13)
 
  call MyCalcDet(BSMat,nc,sol%BSdet)
  !--------------------------------------
end subroutine EnergyPoint
!****************************************************************************************************
subroutine FindResPos(iens,bfield,mass,rmatch,L,NumChan,NumClosed,Res,V,sol)
  use DataStructures
  use modb
  implicit none
  double precision rmatch,eleft,eright,emid,tol,mass,de,diffE,bfield
  integer L,i,imax,NumChan,NumClosed,iefine,nefine,iens,count
  double precision, allocatable :: efine(:), dtaude(:)
  type(PotentialData) :: Vtemp
  type(PotentialData), intent(in) :: V  ! input
  type(SolutionData), intent(in) :: sol ! input
  type(SolutionData) :: solmid,solleft,solright
  type(SolutionData), allocatable :: sola(:)
  type(ResonanceData) :: Res
!  logical condition
  tol = 1d-12
  eleft = Res%el
  eright = Res%er
  emid = 0.5d0*(eleft+eright)
!  write(6,*) "(el,er) = ", Res%el, Res%er

  nefine = 10
  allocate(sola(nefine))
  allocate(efine(nefine))!,dtaude(nefine))

  do count = 1, 10
     de = (eright-eleft)/dble(nefine-1)
     do iefine = 1, nefine
        efine(iefine) = eleft + de*(iefine - 1)
     enddo
     do iefine = 1, nefine
        !write(6,*) "calling energypoint in findrespos...", efine(iefine)
        call EnergyPoint(rmatch, L, mass, efine(iefine), NumChan, NumClosed,V,solmid)
        sola(iefine) = solmid
     enddo
     
     do iefine = 1, nefine-1
        if( (sola(iefine)%dtaude*sola(iefine+1)%dtaude).lt.0d0 ) then
           eright = efine(iefine+1)
           eleft = efine(iefine)
           Res%el = efine(iefine)
           Res%er = efine(iefine+1)
           !Res%eres = (efine(iefine) + efine(iefine+1))*0.5d0
           Res%eres = efine(iefine)
           Res%timedelay = sola(iefine)%timedelay!2d0*(sola(iefine+1)%phase - sola(iefine)%phase)/(efine(iefine+1) - efine(iefine))
           Res%width = 4d0/Res%timedelay
           Res%Fanoq = -1d0/tan(sola(iefine)%bgphase)
           Res%zclose = sola(iefine)%zclose
           Res%bgampsqr = sin(sola(iefine)%bgphase)**2
!           write(6,*) "(er-el) = ",eright - eleft
           !write(6,*) "(el,er) = ",iefine, Res%el, Res%er
           exit
        endif
     enddo
     if ((eright-eleft).lt.tol) exit
  enddo
  
end subroutine FindResPos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fdfgdg(L,k,r,f,df,g,dg) ! returns the oscillatory Riccati functions and their radial derivatives.
  use modb
  implicit none
  integer L
  double precision k, r, f, df, g, dg, nu, x, factor
  double precision J, dJ, RTPIO2, Y, dY
  PARAMETER (RTPIO2=1.25331413731550d0)

  x = k*r
  if(L.lt.0.d0.or.x.le.0.d0) write(6,*) 'bad arguments in sphbesjy'

  nu=L+0.5d0
!  write(6,*) "k, r:", k, r
  call bessjy(x,nu,J,Y,dJ,dY)
  factor=RTPIO2*sqrt(x)
  f=factor*J
  g=factor*Y
  df=k*factor*(dJ+J/(2.d0*x))
  dg=k*factor*(dY+Y/(2.d0*x))
  
end subroutine fdfgdg
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine bfdfgdg(L,k,r,xscale,f,df,g,dg,ldi,ldk) ! returns the exponential-type Riccati functions and their radial derivatives.
  use modb
  implicit none
  integer L
  double precision k, r, f, df, g, dg, nu,xscale
  double precision RTPIO2, RT2OPI, factor, x
  double precision alpha, beta, alphap, betap, ldi, ldk
  PARAMETER (RTPIO2=1.25331413731550d0) 
  PARAMETER (RT2OPI=0.7978845608028654d0)
  x=k*r
  if(L.lt.0.d0.or.x.le.0.d0) write(6,*) 'bad arguments in sphbesik, (L, x) = ', L, x
  nu = L + 0.5d0
  call MyScaledBessIK(x, nu, alpha, beta, alphap, betap, ldi,ldk)
  factor=RT2OPI*sqrt(x)

  f = factor*exp(x-xscale)*alpha
  df = k*factor*exp(x-xscale)*( alpha*(1.d0 + 1.d0/(2d0*x)) + alphap )

  g = factor*exp(xscale-x)*beta
  dg = k*factor*exp(xscale-x)*( beta*(-1.d0 + 1.d0/(2.d0*x)) + betap )
  !  ldi = (alphap/alpha + 1.d0 + 1.d0/(2.d0*x))*k
  !  ldk = (betap/beta - 1.d0 + 1.d0/(2.d0*x))*k

  ldi = df/f 
  ldk = dg/g
  
end subroutine bfdfgdg
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MyCalcDet(A,NA,det)
  implicit none
  integer NA,LDA,IPIV(NA),info,i
  double precision A(NA,NA),det,sgn
  double precision ALocal(NA,NA)
  ALocal=A
  LDA=NA
  call DGETRF( NA, NA, ALocal, LDA, IPIV, INFO )
  if(info < 0) write(6,*) "DGETRF illegal argument i = ",-info
  if(info > 0) write(6,*) "DGETRF U(i,i) = 0 for i = ",info
  sgn = 1.d0
  det = 1.d0
  do i=1,NA
     det = det*ALocal(i,i)
     if(IPIV(i).ne.i) then
        sgn = -sgn
     else
        sgn = sgn
     endif
  enddo
  det = sgn*det
end subroutine MyCalcDet
!****************************************************************************************************
double precision function kdelta(i,j)
  integer i, j

  if(i.eq.j) then
     kdelta=1.d0
  else
     kdelta=0.0d0
  endif
  return
end function kdelta
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine makeF(numf,R,L,mass,energy,Lambda,f,df,dfde,ddfde,d2f,d2df,d3f,d3df)
 ! use potmod
  use modb
  implicit none
  integer i,j,k,lm,numf,L
  real*8 ek,ee,energy,x,R,xscale,ldk,ldi,mass,twomu,dek,d2ek
  real*8 Lambda(numf)
  real*8 f(numf), df(numf),dfde(numf), ddfde(numf)
  real*8 d2f(numf),d2df(numf),d3f(numf),d3df(numf)
  real*8 rj,ry,drj,dry,pre,rk,ri,drk,dri,fpp,fp
  
  f=0d0; df=0d0;
  dfde=0d0; ddfde=0d0;
  d2f=0d0; d2df=0d0;
  d3f=0d0; d3df=0d0
  
  twomu=2.d0*mass
  
  do k=1,numf
     ee=energy-Lambda(k)
     if(ee.gt.0d0)then
        ek = sqrt(twomu*ee) 
        x=ek*R
        dek = mass/ek
        d2ek = -mass**2/ek**3

!!$  The following for testing purposes only:
!!$        write(6,*) "in make F; x = ", x
!!$        write(6,*) "in make F; ek = ", ek
!!$        write(6,*) "in make F; R = ", R
!!$        write(6,*) "in make F; mass = ", mass
!!$        write(6,*) "in make F; ee = ", ee
!!$        write(6,*) "in make F; Lambda (",k,") = ", Lambda(k)
        
        call fdfgdg(L,ek,r,rj,drj,ry,dry)
        !        write(6,*) "done with call from makeF..."
        !NOTE Here, drj indicates the derivative of the Riccati function with r, while fp indicates the derivative with kr
        f(k) = rj
        df(k) = drj
        fp = drj/ek
        fpp = (L*(L+1)/x**2 - 1)*rj
        
        dfde(k) = drj*r/(2.d0*ee)
        ddfde(k) = drj/(2.d0*ee) + mass*r*fpp
        
        d2f(k) = r*(fp*d2ek + dek**2*r*fpp)
        d2df(k) = (mass**2*(-2*rj*L*(1 + L) + &
             ek*r*(ek*fpp*r + fp*(-1 + L + L**2 - ek**2*r**2))))/(ek**4*r)

        d3f(k) = -((-(ek*fp*(3 + L + L**2)*r) + ek**3*fp*r**3 + &
             ek**2*r**2*(fpp - 2*rj) + 4*L*(1 + L)*rj)*twomu**3)/(8.d0*ek**6)
        d3df(k) = ((ek*fp*(3 - 4*L*(1 + L))*r - ek**4*fpp*r**4 + 4*L*(1 + L)*rj + &
             ek**2*r**2*(fpp*(-1 + L + L**2) + 2*rj))*twomu**3)/(8.d0*ek**6*r) 
        
     else
        ek = sqrt(-2.d0*mass*ee) 
!        write(6,*) "Calling exponentially growing function in makeF"
!        read(*,*)
        x=ek*R
        call bfdfgdg(L,ek,r,ek,ri,dri,rk,drk,ldi,ldk) ! really x*sphbes
        f(k)=ri
        df(k)=dri
        fp = dri/ek
        fpp = (L*(L+1)/x**2 + 1)*ri
        
        dfde(k)=dri*r/(2.d0*ee)
        ddfde(k) = dri/(2.d0*ee) - mass*r*ri*( (L*(L + 1))/(x*x) + 1)
        d2f(k) = ((-(dri*ek*r) + (L + L**2 + ek**2*r**2)*ri)*twomu**2)/(4.d0*ek**4)
        d2df(k) = ((dri*ek*r*(-1 + L + L**2 + ek**2*r**2) - &
             (L + L**2 - ek**2*r**2)*ri)*twomu**2)/(4.d0*ek**4*r)

        d3f(k) = -((3*ek*fp*r + ek*fp*L*r + ek*fp*L**2*r - ek**2*fpp*r**2 + &
            ek**3*fp*r**3 - 4*L*ri - 4*L**2*ri - 2*ek**2*r**2*ri)*twomu**3)/ &
            (8.d0*ek**6)
        d3df(k) =  -((3*ek*fp*r - 4*ek*fp*L*r - 4*ek*fp*L**2*r - ek**2*fpp*r**2 + &
            ek**2*fpp*L*r**2 + ek**2*fpp*L**2*r**2 + ek**4*fpp*r**4 + &
            4*L*ri + 4*L**2*ri - 2*ek**2*r**2*ri)*twomu**3)/(8.d0*ek**6*r)
        
     endif
     !write(6,*)k,pre,rj,drj
  enddo
end subroutine makeF
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module teqparams
  double precision v0,r0
end module teqparams
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function teq0(energy)
  use teqparams
  double precision energy
  teq0 = sqrt(v0+energy)*cos(r0*sqrt(v0+energy)) + sqrt(-energy)*sin(r0*sqrt(v0+energy))
  return
end function teq0
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module teqparamsv0
  double precision energy,r0
end module teqparamsv0
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function teqv0(v0)
  use teqparamsv0
  double precision v0
  teqv0 = sqrt(v0+energy)*cos(r0*sqrt(v0+energy))+sqrt(-energy)*sin(r0*sqrt(v0+energy))
  return
end function teqv0
!****************************************************************************************************
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
!****************************************************************************************************
! Make the potential matrix model for Dysprosium (or Dy like atom)
! NumEven = number of even-parity (ground state) closed-channel resonances.
! DBOpenRes = separation in magnetic field (in gauss) of the broad open-channel features
! DVClosedWindow = window in energy over which the Feshbach resonances should live.
! DBClosedWindow = Window in B field over which the Feshbach resonances should live
! dmu = range of magnetic moments = DVClosedWindow/DBClosedWindow
! BC = binding of lowest-energy closed channel resonance.
! ei0 = Initial energy position of all closed channel resonances,
!       the first NumEven of which are the corresponding ground state of the closed channel
! nopen = number of bound states supported by the open channel
!****************************************************************************************************
!****************************************************************************************************
subroutine makeVOld(NumChan,NumOpen,NumEven,vopen,BCE,BCO,gco,gcc,eoffset,Vmat,ethresh)
  use modb
  implicit none
  integer NumChan,NumOpen,NumEven,nopen,i,j
  double precision Vmat(NumChan,NumChan),eoffset(NumChan - NumOpen),ethresh(NumChan),VOD(NumChan,NumChan)
  double precision BCE,BCO,ECE,ECO,gco,gcc,cdlength
  double precision vcE,vcO,tol
  double precision, intent(in) :: vopen
  !  double precision, allocatable :: eopen(:)
  integer mode
  external rgnf_lux

  cdlength = 1.d0
  Vmat = 0.d0 !initialize the potential matrix to zero
  mode = 2 !potential depth search
  tol = 1d-14
  ECE = -BCE
  ECO = -BCO
  call SquareWellEnergy(vcE,ECE,mode,1,tol)
  call SquareWellEnergy(vcO,ECO,mode,2,tol)

  nopen = floor(sqrt(vopen)/pi - 0.5d0)
!  allocate(eopen(nopen))
!!$  mode = 1 !switch to energy search
!!$  do i = 1, nopen
!!$     call SquareWellEnergy(vcO,eopen(i),mode,i,tol)
!!$  end do

!!$  write(6,*) "NumEven = ", NumEven
!!$  write(6,*) "NumOpen = ", NumOpen
!!$  write(6,*) "NumChan = ", NumChan
  
  call GOE(NumChan,VOD,rgnf_lux) ! call a gaussian matrix
  vmat = gcc*VOD  ! assign the gaussian matrix to the potential 
!    vmat = VOD  ! assign the gaussian matrix to the potential 

  ! now write over the diagonal parts so only the off-diagonals are gaussian
  do i = 1, NumEven
     Vmat(i,i) = BCE - vcE + eoffset(i)
     ethresh(i) = BCE + eoffset(i)
     write(6,*) "Vmat-diag-even",i,Vmat(i,i), ethresh(i), eoffset(i)
!!$
!!$          Vmat(i,i) = BCO - vcO + eoffset(i) 
!!$     ethresh(i) = BCO + eoffset(i)
!!$     write(6,*) "Vmat-diag-odd",i,Vmat(i,i),ethresh(i), eoffset(i)

  enddo
  do i = NumEven + 1, NumChan - NumOpen
     Vmat(i,i) = BCO - vcO + eoffset(i) 
     ethresh(i) = BCO + eoffset(i)
     write(6,*) "Vmat-diag-odd",i,Vmat(i,i),ethresh(i), eoffset(i)

!!$     Vmat(i,i) = BCE - vcE + eoffset(i)
!!$     ethresh(i) = BCE + eoffset(i)
!!$     write(6,*) "Vmat-diag-even",i,Vmat(i,i), ethresh(i), eoffset(i)

  enddo
  
  do i = NumChan - NumOpen + 1, NumChan
     
     Vmat(i,i) = -vopen
     ethresh(i) = 0.d0
     write(6,*) "Vmat-diag-open",i,Vmat(i,i),ethresh(i)!, eoffset(i)
  enddo
!  write(6,*) "press enter to continue..."
!  read(*,*)
  !Now set the open-closed coupling separately
!  do i = 1, NumChan-NumOpen
!     vmat(i,NumChan) = gco
!     vmat(NumChan,i) = gco
!  enddo

  !!$  do i = 1,NumChan - NumOpen
!!$     do j = 1,NumChan - NumOpen
!!$        if(i.ne.j) Vmat(i,j) = gcc*exp(-dble(abs(i-j)-1)/cdlength) 
!!$     enddo
!!$  enddo
!!$  do i = 1,NumChan - NumOpen
!!$     do j = 1,NumChan - NumOpen
!!$        if(i.ne.j) Vmat(i,j) = gcc*exp(-dble(abs(i-j)-1)/cdlength) 
!!$     enddo
!!$  enddo
!  deallocate(eopen)
end subroutine makeVOld
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine makeVNew(NumChan,NumOpen,NumEven,vopen,BCE,BCO,eoffset,Vmat,ethresh)
  use modb
  implicit none
  integer NumChan,NumOpen,NumEven,nopen,i,j
  double precision Vmat(NumChan,NumChan),eoffset(NumChan - NumOpen),ethresh(NumChan),VOD(NumChan,NumChan)
  double precision BCE,BCO,ECE,ECO,cdlength
  double precision vcE,vcO,tol
  double precision, intent(in) :: vopen
  !  double precision, allocatable :: eopen(:)
  integer mode
  external rgnf_lux

  cdlength = 1.d0
  Vmat = 0.d0 !initialize the potential matrix to zero
  mode = 2 !potential depth search
  tol = 1d-14
  ECE = -BCE
  ECO = -BCO
  call SquareWellEnergy(vcE,ECE,mode,1,tol)
  call SquareWellEnergy(vcO,ECO,mode,2,tol)

  nopen = floor(sqrt(vopen)/pi - 0.5d0)
!  allocate(eopen(nopen))
!!$  mode = 1 !switch to energy search
!!$  do i = 1, nopen
!!$     call SquareWellEnergy(vcO,eopen(i),mode,i,tol)
!!$  end do

!!$  write(6,*) "NumEven = ", NumEven
!!$  write(6,*) "NumOpen = ", NumOpen
!!$  write(6,*) "NumChan = ", NumChan
  
  call GOE(NumChan,VOD,rgnf_lux) ! call a gaussian matrix
!  vmat = gcc*VOD  ! assign the gaussian matrix to the potential
    vmat = VOD  ! assign the gaussian matrix to the potential 

  ! now write over the diagonal parts so only the off-diagonals are gaussian
  do i = 1, NumEven
     Vmat(i,i) = BCE - vcE + eoffset(i)
     ethresh(i) = BCE + eoffset(i)
     write(6,*) "Vmat-diag-even",i,Vmat(i,i), ethresh(i), eoffset(i)
!!$
!!$          Vmat(i,i) = BCO - vcO + eoffset(i) 
!!$     ethresh(i) = BCO + eoffset(i)
!!$     write(6,*) "Vmat-diag-odd",i,Vmat(i,i),ethresh(i), eoffset(i)

  enddo
  do i = NumEven + 1, NumChan - NumOpen
     Vmat(i,i) = BCO - vcO + eoffset(i) 
     ethresh(i) = BCO + eoffset(i)
     write(6,*) "Vmat-diag-odd",i,Vmat(i,i),ethresh(i), eoffset(i)

!!$     Vmat(i,i) = BCE - vcE + eoffset(i)
!!$     ethresh(i) = BCE + eoffset(i)
!!$     write(6,*) "Vmat-diag-even",i,Vmat(i,i), ethresh(i), eoffset(i)

  enddo
  
  do i = NumChan - NumOpen + 1, NumChan
     
     Vmat(i,i) = -vopen
     ethresh(i) = 0.d0
     write(6,*) "Vmat-diag-open",i,Vmat(i,i),ethresh(i)!, eoffset(i)
  enddo

!  write(6,*) "press enter to continue..."
!  read(*,*)
  !Now set the open-closed coupling separately
!!$  do i = 1, NumChan-NumOpen
!!$     vmat(i,NumChan) = gco
!!$     vmat(NumChan,i) = gco
!!$  enddo

  !!$  do i = 1,NumChan - NumOpen
!!$     do j = 1,NumChan - NumOpen
!!$        if(i.ne.j) Vmat(i,j) = gcc*exp(-dble(abs(i-j)-1)/cdlength) 
!!$     enddo
!!$  enddo
!!$  do i = 1,NumChan - NumOpen
!!$     do j = 1,NumChan - NumOpen
!!$        if(i.ne.j) Vmat(i,j) = gcc*exp(-dble(abs(i-j)-1)/cdlength) 
!!$     enddo
!!$  enddo
!  deallocate(eopen)
end subroutine makeVNew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! This next subroutine modifies the potential by the magnetic field
subroutine ShiftV(Bfield,NumChan,NumOpen,DVClosedWindow,DBClosedWindow,DBOpenRes,V0,Vshifted)
  use modb
  use DataStructures
  implicit none
  integer NumChan,NumOpen,i,NumClosed,j
!  double precision ethresh(NumChan),ethresh0(NumChan),Vmat0(NumChan,NumChan),Vmat(NumChan,NumChan)
  double precision dmu, DVClosedWindow,DBClosedWindow,DBOpenRes,vo1,vo2,dOpenDepth
  double precision, intent(in) :: Bfield
  type(PotentialData) :: V0, Vshifted
!!$  write(6,*) "*V0%V:"
!  call printmatrix(V0%V,NumChan,NumChan,6)
  Vshifted = V0
  NumClosed = NumChan - NumOpen
  dmu = DVClosedWindow/DBClosedWindow ! range of magnetic moments needed
!  write(6,*) "dmu = ", dmu, Bfield
  vo1 = (0.5d0*pi)**2  ! depth for which the ground state is at threshold
  vo2 = (1.5d0*pi)**2  ! depth for which the first excited state is at threshold
  dOpenDepth = (vo2 - vo1)/DBOpenRes
  do i = 1, NumChan - NumOpen
     Vshifted%V(i,i) = V0%V(i,i) + dmu*Bfield
     Vshifted%ethresh(i) = V0%ethresh(i) + dmu*Bfield
  enddo
  do i = NumChan - NumOpen + 1, NumChan
     Vshifted%V(i,i) = V0%V(i,i) -  dOpenDepth * Bfield
  end do

  do i = 1,NumClosed
     do j = 1,NumClosed
        Vshifted%Vcc(i,j) = Vshifted%V(i,j)!(1::NumClosed,1::NumClosed)
     enddo
  enddo
  call MyDSYEV(Vshifted%V,NumChan,Vshifted%Lambda,Vshifted%U) !diagonalize the interior hamiltonian
  call MyDSYEV(Vshifted%Vcc,NumClosed,Vshifted%LambdaC,Vshifted%UC) !diagonalize the closed-channel sector only
  
!!$  write(6,*) "*Vshifted%V:"
!!$  call printmatrix(Vshifted%V,NumChan,NumChan,6)
!!$
!!$  write(6,*) "Vshifted%LambdaC:"
!!$  call printmatrix(Vshifted%LambdaC,1,NumClosed,6)
!!$
!!$  write(6,*) "Vshifted%Vcc:"
!!$  call printmatrix(Vshifted%Vcc,NumClosed,NumClosed,6)
!!$
!!$  write(6,*) "V0%Vcc:"
!!$  call printmatrix(V0%Vcc,NumClosed,NumClosed,6)
end subroutine ShiftV
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine testbessik
  implicit none
  real*8 ri,rk,dri,drk, x, dx, ldi, ldk,xscale,nu0, r, k
  real*8 si, sk, sip, skp, alpha, beta, alphap, betap
  real*8 v, vm,pi,xinit,xfinal
  complex*16, allocatable :: YCI(:), YCK(:), YDI(:), YDK(:)
!  complex*16 z
  integer ixtot, L, ix
  pi = acos(-1.d0)
  !write(6,*) 'pi = ',pi

  L=0
  nu0 = 0.5d0
  xinit=0.001d0
  xfinal=100.d0
  ixtot=100
  dx=(xfinal-xinit)/dble(ixtot-1)

  xscale = 0.d0
  
!  write(6,*) dx, ixtot
!  read(*,*)
  !  x=0.0d0
  k = 1.d0

  do ix = 1, ixtot
     x = xinit + (ix-1)*dx
     r = x/k
!     xscale = x
     !     call Mysphbesik(L,x,xscale,si,sk,sip,skp,ldk,ldi) ! really x*sphbes
     call bfdfgdg(L,k,r,xscale,ri,dri,rk,drk,ldi,ldk) ! really x*sphbes
!     call sphbesik(L,x,si,sk,sip,skp) ! really x*sphbes
     write(6,10) x, ri, dri, rk, drk, ldi, ldk
     
  end do

10 format(100e24.14)
end subroutine testbessik
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine background(energy,phase,n,V)
  use potmod
  use modb
  implicit none
  integer i,j,k,nc,nr,ir,n
  real*8 ek,ee,energy,x,ko,kc,xx,V(n,n)
  real*8 y,g,f,df,dg,kmat,phase
  real*8 pvec(n),rj,ry,rjp,ryp,norm
  ko=sqrt(2d0*mass*energy)
  kc=sqrt(2d0*mass*(energy-V(n,n))) 
  xx=kc*Rmatch
  call sphbesjy(L,xx,rj,ry,rjp,ryp)
  y=kc*rjp/rj
  xx=ko*Rmatch
  call sphbesjy(L,xx,f,g,df,dg) ! psi=f-k g
  kmat=(y*f-ko*df)/(y*g-ko*dg)
  phase=atan(kmat)
end subroutine background

!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine makePSI(mass,L,energy,Lambda,U,Amp,aclose,bcoef,cd,sd,n,ethresh)
  !this makes the wavefunctions for a given energy and U, but doesn't save them.
  !Uses primitive integration to calculate closed-channel amplitude
  !should replace with gauss-legendre integration or some such.
!  use potmod
  use modb
  implicit none
  integer i,j,k,nc,nr,ir,n,L
  real*8, allocatable :: kappa(:)
  real*8 ek,ee,energy,x,ko,kc,sd,cd,rinit,rfin,dr,dl,xx,r,mass,Rmatch
  real*8 Lambda(n),amp(n-1),U(n,n),aclose(n-1),bcoef(n)
  real*8 f(n),df(n),dfde(n), ddfde(n),d2f(n),d2df(n),d3f(n),d3df(n),psi(n),ethresh(n)
  real*8 pvec(n),rj,ry,drj,dry,norm, ri,rk,dri,drk,ldi,ldk,xscale
  Rmatch = 1d0
  nc=n-1
  allocate(kappa(nc))
!  kc=sqrt(2d0*mass*(eth1-energy))
  ko=sqrt(2d0*mass*energy)
  norm=1d0/sqrt(cd*cd+sd*sd) !normalize the open channel amplitude
  amp=0d0
  rinit=0.01d0
  rfin=Rmatch*2.5d0 ! increase if there are very weakly closed channels
  !  xx=-minval(Lambda)+eth1
  xx=-minval(Lambda)+ethresh(1)
  dr=6.3d0/sqrt(2d0*mass*xx)/20d0
  nr=int((rfin-rinit)/dr)
  do ir=1,nr
     r=rinit+dr*(ir-1)
     if(r.lt.rmatch) then
!        write(6,*) "Calling makeF from makePsi..."
!        write(6,*) "in makePsi, r = ", r
!        write(6,*) "in makePsi, energy = ", energy
!        write(6,*) "in makePsi, Lambda = ", Lambda
        call makeF(n,r,0,mass,energy,Lambda,f,df,dfde,ddfde,d2f,d2df,d3f,d3df) !Set L = 0 inside the well
        pvec=f*bcoef!matmul(f,bcoef)
        psi=norm*matmul(U,pvec)
     else
        do i=1,nc
           kappa(i) = sqrt(2d0*mass*(ethresh(i)-energy))
           !psi(i)=norm*aclose(i)*exp(-kappa(i)*(r-rmatch+1d0))  !  CLOSED PSI r>Rmatch
           !psi(i)=norm*aclose(i)*exp(-kc*(r-rmatch+1d0))  !  CLOSED PSI r>Rmatch

           !replace the above with the sqrt(x)*besselK
           xscale=kappa(i)*rmatch
           x=kappa(i)*r
           if(x.le.0.d0) then
              write(6,*) 'sending bad argument to Mysphbesik, i, kappa(i), r, x = ', i,kappa(i),r,x
           endif
!           write(6,*) "in makePSI calling bound bessel..."
           call bfdfgdg(L,kappa(i),r,xscale,ri,dri,rk,drk,ldi,ldk) ! really x*sphbes
           psi(i)=norm*aclose(i)*rk
           
        enddo
        
        xx=ko*r
        !call sphbesjy(L,xx,rj,ry,rjp,ryp)
!        write(6,*) "in makePSI calling open channel fdfgdg... ko = ", ko
!        write(6,*) "in makePSI calling open channel fdfgdg... r = ", r
        call fdfgdg(L,ko,r,rj,drj,ry,dry)
!        write(6,*) "done with call..."

        psi(n)=norm*(cd*rj-sd*ry)
     endif
     !write(100,*)r,(psi(j),j=1,n)
     do i=1,nc
        amp(i)=amp(i)+dr*psi(i)**2 !closed channel amplitude in channel i
     enddo
  enddo
!  deallocate(kappa)
end subroutine makePSI
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine makeAmat(energy,Lambda,U,Amat,d1amat,d2amat,d3amat,n,eth)
  use potmod
  use modb
  implicit none
  integer i,j,k,n
  real*8 ek,ee,energy,x
  real*8 Lambda(n),Amat(2*n,n),d1amat(2*n,n),U(n,n),eth(n),d2amat(2*n,n),d3amat(2*n,n)
  real*8 Phi(n,n), dPhi(n,n),dPhidE(n,n),ddPhidE(n,n),d2phi(n,n),d2dphi(n,n),d3phi(n,n),d3dphi(n,n)
  real*8 f(n), df(n), dfde(n),ddfde(n),d2f(n),d2df(n),d3f(n),d3df(n)

  Amat=0d0
  d1amat=0d0 ! involves energy derivative
  d2amat=0d0 ! second energy derivative matrix
  d3amat=0d0 ! third energy derivative matrix
  
!  write(6,*) "in makeAmat, energy = ", energy
!  write(6,*) "Calling makeF from makeAmat..."
  call makeF(n,rmatch,0,mass,energy,Lambda,f,df,dfde,ddfde,d2f,d2df,d3f,d3df) !make the diagonal f and df matrices (L=0 inside well)
!  write(6,*) "the f matrix:"
!  call printmatrix(f,n,n,6)
!  stop
!!$  Phi=matmul(U,f) ! transform to open-closed basis; make the interior PHI_(i,alpha) matrix
!!$  dPhi=matmul(U,df)
!!$  dPhidE=matmul(U,dfde)
!!$  ddPhidE=matmul(U,ddfde)
!!$  d2phi = matmul(U,d2f)
!!$  d2dphi = matmul(U,d2df)
  do k = 1, n
     Phi(1:n,k) = U(1:n,k)*f(k)
     dPhi(1:n,k) = U(1:n,k)*df(k)
     dPhidE(1:n,k) = U(1:n,k)*dfde(k)
     ddPhidE(1:n,k) = U(1:n,k)*ddfde(k)
     d2Phi(1:n,k) = U(1:n,k)*d2f(k)
     d2dPhi(1:n,k) = U(1:n,k)*d2df(k)
     d3phi(1:n,k) = U(1:n,k)*d3f(k)
     d3dphi(1:n,k) = U(1:n,k)*d3df(k)
  enddo

!  call printmatrix(d2dphi,n,n,6)
  ! Load the interior functions into the A matrix
  do i=1,n
     do j=1,n
        Amat(2*i-1,j) = Phi(i,j)
        Amat(2*i  ,j) = dPhi(i,j)
        d1amat(2*i-1,j) = dPhidE(i,j)
        d1amat(2*i,  j) = ddPhidE(i,j)
        d2amat(2*i-1,j) = d2phi(i,j)
        d2amat(2*i,  j) = d2dphi(i,j)
        d3amat(2*i-1,j) = d3phi(i,j)
        d3amat(2*i,  j) = d3dphi(i,j) 
     enddo
  enddo
  
end subroutine makeAmat
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine makeBmat(energy,Bmat,d1bmat,d2bmat,d3bmat,Gmat,nc,n,eth)
  use potmod
  use modb
  implicit none
  integer i,j,k,nc,n
  real*8 ek,ee,energy,x,eth(n)
  real*8 Lambda(n),Bmat(2*n,n+1),d1bmat(2*n,n+1),Gmat(nc,nc),d2bmat(2*n,n+1),d3bmat(2*n,n+1)
  real*8 f,df,g,dg,fp,gp,fpp,gpp
  real*8 rj,ry,drj,dry,r
  real*8 ri,rk,dri,drk,xscale,ldi,ldk,twomu
  twomu = 2.d0*mass
  Bmat=0d0
  d1bmat=0d0
  d2bmat=0d0
  r=rmatch
  nc=n-1
  !closed channels
  Gmat=0.d0
  do i=1,nc
     ee = eth(i)-energy !NPM
     ek = sqrt(twomu*ee)
     xscale = ek
     x = ek*r
     call bfdfgdg(L,ek,r,xscale,ri,dri,rk,drk,ldi,ldk) ! really x*sphbesik (the modified Riccati-Bessel function)
     g = rk
     dg = drk
     gp = drk/ek
     gpp = g*( L*(L+1)/(x*x) + 1)

     Bmat(2*i-1,i) = g
     Bmat(2*i  ,i) = dg

     d1bmat(2*i-1,i) = -dg*r/(2*ee)
     d1bmat(2*i  ,i) = -dg/(2.d0*ee) - mass*r*gpp

     d2bmat(2*i-1,i) =-((drk*r - L*rk - L**2*rk - ek**2*r**2*rk)*twomu**2)/(4.d0*ek**4)
     d2bmat(2*i,i) = ((-(drk*r) + drk*L*r + drk*L**2*r + drk*ek**2*r**3 - L*rk - &
          L**2*rk + ek**2*r**2*rk)*twomu**2)/(4.d0*ek**4*r)

     d3bmat(2*i-1,i) = -((3*ek*gp*r + ek*gp*L*r + ek*gp*L**2*r - ek**2*gpp*r**2 + &
          ek**3*gp*r**3 - 4*L*rk - 4*L**2*rk - 2*ek**2*r**2*rk)*twomu**3)/ &
          (8.d0*ek**6)
     d3bmat(2*i,i) = -((3*ek*gp*r - 4*ek*gp*L*r - 4*ek*gp*L**2*r - ek**2*gpp*r**2 + &
            ek**2*gpp*L*r**2 + ek**2*gpp*L**2*r**2 + ek**4*gpp*r**4 + &
            4*L*rk + 4*L**2*rk - 2*ek**2*r**2*rk)*twomu**3)/(8.d0*ek**6*r)
     
     Gmat(i,i)=ldk
  enddo

  !open channel
  ee = energy-eth(n) !NPM
  ek = sqrt(2.d0*mass*ee)
  x = ek*r
  call fdfgdg(L,ek,r,f,df,g,dg)
  fp = f/ek
  fpp = f*( L*(L+1)/(x*x) - 1)
  gp = g/ek
  gpp = g*( L*(L+1)/(x*x) - 1)
  
  Bmat(2*n-1,n)=f
  Bmat(2*n  ,n)=df

  Bmat(2*n-1,n+1)=-g
  Bmat(2*n  ,n+1)=-dg

  d1bmat(2*n-1,n) = df*r/(2.d0*ee)
  d1bmat(2*n  ,n) = df/(2.d0*ee) + mass*r*fpp

  d1bmat(2*n-1,n+1) = -dg*r/(2.d0*ee)
  d1bmat(2*n  ,n+1) = -dg/(2.d0*ee) - mass*r*gpp

  d2bmat(2*n-1,n) = -((-(f*L) - f*L**2 + df*r + ek**2*f*r**2)*twomu**2)/(4.d0*ek**4)
  d2bmat(2*n, n) =  -((f*L + f*L**2 + df*r - df*L*r - df*L**2*r + ek**2*f*r**2 + &
         df*ek**2*r**3)*twomu**2)/(4.d0*ek**4*r)
  
  d2bmat(2*n-1,n+1)  = ((-(g*L) - g*L**2 + dg*r + ek**2*g*r**2)*twomu**2)/(4.d0*ek**4)
  d2bmat(2*n, n+1) =   ((g*L + g*L**2 + dg*r - dg*L*r - dg*L**2*r + ek**2*g*r**2 + &
            dg*ek**2*r**3)*twomu**2)/(4.d0*ek**4*r)

  d3bmat(2*n-1,n) =  -((-(ek*fp*(3 + L + L**2)*r) + ek**3*fp*r**3 + &
       ek**2*r**2*(fpp - 2*f) + 4*L*(1 + L)*f)*twomu**3)/(8.d0*ek**6)
  d3bmat(2*n,n) = ((ek*fp*(3 - 4*L*(1 + L))*r - ek**4*fpp*r**4 + 4*L*(1 + L)*f + &
       ek**2*r**2*(fpp*(-1 + L + L**2) + 2*f))*twomu**3)/(8.d0*ek**6*r)

  d3bmat(2*n-1,n+1) =  -((-(ek*gp*(3 + L + L**2)*r) + ek**3*gp*r**3 + &
       ek**2*r**2*(gpp - 2*g) + 4*L*(1 + L)*g)*twomu**3)/(8.d0*ek**6)
  d3bmat(2*n,n+1) = ((ek*gp*(3 - 4*L*(1 + L))*r - ek**4*gpp*r**4 + 4*L*(1 + L)*g + &
       ek**2*r**2*(gpp*(-1 + L + L**2) + 2*g))*twomu**3)/(8.d0*ek**6*r)

  
end subroutine makeBmat
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine makeABmatNew

end subroutine makeABmatNew
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine printmatrix(M,nr,nc,file)
  implicit none
  integer nr,nc,file,j,k
  double precision M(nr,nc)

  do j = 1,nr
     write(file,10) (M(j,k), k = 1,nc)
  enddo

20 format(1P,100e16.8)
10 format(100e14.5)
30 format(100F14.5)
  
end subroutine printmatrix
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine diagH(n,Lambda,U,H)
  implicit none
  integer n,k,i,ip
  real*8 Lambda(n),U(n,n),UT(n,n),work(3*n),H(n,n),xx,sgn(n),Phi(n,n)
  UT=H
  call DSYEV('V','U',N,UT,N,Lambda,work,3*N,k)
  do i=1,n
     sgn(i)=UT(i,i)/abs(UT(i,i))
  enddo
  do i=1,n
     xx=0d0
     do ip=1,n
        xx=xx+UT(ip,i)**2
     enddo
     U(:,i)=sgn(i)*UT(:,i)/sqrt(xx) ! use for psi = U.phi project dressed on undressed
  enddo

end subroutine diagH
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE qgaus(t,ss,nleg,wleg,xleg)
  REAL*8 ss,t(nLeg),wLeg(nLeg),xLeg(nLeg)
  Integer j,n
  ss=0.d0
  do j=1,nLeg/2
     ss=ss+wLeg(j)*(t(j)+t(nleg-j+1))
  enddo
  ss=ss
  return
END SUBROUTINE qgaus
!****************************************************************************************************
! subroutine SquareWellEnergy will calculate either the potential depth required to give a particular energy, or the energy for a particular potential depth.
! This requires solving a transcendental equation using Newton's method.
! Use mode = 1 for energy given v0
! Use mode = 2 for v0 given energy
!****************************************************************************************************
subroutine SquareWellEnergy(v0,energy,mode,n,tol)
  use modb
  implicit none
  integer, intent(in) :: n, mode
  integer i,j,imax
  double precision v0,energy,res,t1,t2,r0,tol,B,teq,teqprime
  double precision vmin,vmax,vmid,emin,emid,emax,teqa,teqb,teqc

  r0=1.d0
  imax = 100

  i=1;
!  write(6,*) "# mode = ", mode
  if (mode.eq.2) then
     ! make the initial guess for v0 in mode=2
     if (energy.lt.0.d0) then
        v0=0.5d0*( (n - 0.5d0)**2 + n**2)*pi**2 - energy
     else
        goto 200
     endif
     
100  teq=sqrt(v0+energy)*cos(r0*sqrt(v0+energy)) + sqrt(-energy)*sin(r0*sqrt(v0+energy))
     teqprime=((1.d0 + Sqrt(-energy))*Cos(Sqrt(energy + v0)) - Sqrt(energy + v0)*Sin(Sqrt(energy + v0)))/(2.d0*Sqrt(energy + v0))
     ! Update estimate
     v0=v0-(teq/teqprime)
     i=i+1
     if (i>=imax) goto 200
     if (dabs(teq/teqprime)>tol) goto 100
     
  else if (mode.eq.1) then
     ! make the initial guess for the energy in mode = 1
     if (v0.lt.((n - 0.5)**2*pi**2)) goto 200  ! require the potential be deep enough to support the state
     !energy = 0.5d0*( (n - 0.5d0)**2 + n**2)*pi**2 - v0
     emin = ((n - 0.5d0)*pi)**2 - v0
     emax = min((n*pi)**2 - v0,0.d0)
     teqa = sqrt(v0+emin)*cos(r0*sqrt(v0+emin)) + sqrt(-emin)*sin(r0*sqrt(v0+emin))
     teqb = sqrt(v0+emax)*cos(r0*sqrt(v0+emax)) + sqrt(-emax)*sin(r0*sqrt(v0+emax))

     if (teqa*teqb.gt.0.d0) then  ! require the root to be bracketed
        write(6,*) "root not bracketed in mode 2"
        goto 200
     end if
     !This section uses the bisection method.
     emid = 0.5d0*(emax + emin)
102  teqc = sqrt(v0+emid)*cos(r0*sqrt(v0+emid)) + sqrt(-emid)*sin(r0*sqrt(v0+emid))
!     write(6,*) emin, emax, teqa, teqb, teqa*teqb
     if (teqa*teqc.lt.0) then
        emax = emid
        emid = 0.5*(emin + emax)
     else if (teqb*teqc.lt.0) then
        emin = emid
        emid = 0.5*(emin + emax)
     end if
     teqa = sqrt(v0+emin)*cos(r0*sqrt(v0+emin)) + sqrt(-emin)*sin(r0*sqrt(v0+emin))
     teqb = sqrt(v0+emax)*cos(r0*sqrt(v0+emax)) + sqrt(-emax)*sin(r0*sqrt(v0+emax))
     energy = emid
     i = i+1
     if (i.ge.imax) goto 200
     if(dabs(teqc).gt.tol) goto 102
     
  else
     write(6,*) "potential depth not sufficient to support bound state with n = ", n
     energy = 0.d0
     goto 200
  endif
  
200 return
end subroutine SquareWellEnergy
!***********************************************
!*         Newton"s method subroutine          *
!* ------------------------------------------- *
!* This routine calculates the zeros of a      *
!* function Y(x) by Newton"s method.           *
!* The routine requires an initial guess, x0,  *
!* and a convergence factor, e. Also required  *
!* is a limit on the number of iterations, m.  *
!***********************************************
Subroutine Newton(m,n,e,x0,yy,Y)
  implicit none
  real*8 e,x0,yy,y1,Y
  external Y
  
  integer m,n
  n=0
  ! Get y and y1
100 yy=Y(x0,y1)
  ! Update estimate
  x0=x0-(yy/y1)
  n=n+1
  if (n>=m) return
  if (dabs(yy/y1)>e) goto 100
  return
end Subroutine Newton
!***************************************************************
!* Parametric least squares curve fit subroutine. This program *
!* least squares fits a function to a set of data values by    *
!* successively reducing the variance. Convergence depends on  *
!* the initial values and is not assured.                      *
!* n pairs of data values, X(i), Y(i), are given. There are l  *
!* parameters, A(j), to be optimized across.                   *
!* Required are initial values for the A(l) and e. Another     *
!* important parameter which affects stability is e1, which is *
!* initially converted to e1(l) for the first intervals.       *
!* The parameters are multiplied by (1 - e1(i)) on each pass.  *
!***************************************************************
Subroutine Param_LS(l,m,n,d,e,ee1,A,X,Y,func)  
  !Labels: 50,100
  implicit none
  external func
  integer SIZE
  parameter(SIZE=25)
  real*8 d,e,ee1,A(l),X(SIZE),Y(SIZE)
  real*8 E1(l)
  integer i,l,m,n
  real*8 a0,l1,l2,m0,m1
  do i = 1, l
     E1(i) = ee1
  end do;
  !Set up test residual
  l1 = 1.d6
  !Make sweep through all parameters
50 do i = 1, l
     a0 = A(i)
     !Get value of residual
     A(i) = a0
100  call S200(l,l2,n,d,A,X,Y,func)
     !Store result in m0
     m0 = l2
     !Repeat for m1
     A(i) = a0 * (1.d0 - E1(i))
     call S200(l,l2,n,d,A,X,Y,func)
     m1 = l2
     !Change interval size if called for
     !If variance was increased, halve E1(i) 
     if (m1 > m0)  then
        E1(i) = -E1(i) / 2.d0
     end if
     !If variance was reduced, increase step size by increasing E1(i)
     if (m1 < m0)  then
        E1(i) = 1.2d0 * E1(i)
     end if
     !If variance was increased, try to reduce it
     if (m1 > m0) then  
        A(i) = a0
     end if
     if (m1 > m0)  then
        goto 100
     end if
  end do !i loop
  !End of a complete pass
  !Test for convergence 
  m = m + 1
  if (l2.eq.0.d0) then
     return
  end if
  if (dabs((l1 - l2) / l2) < e)  then
     return
  end if
  !If this point is reached, another pass is called for
  l1 = l2
  goto 50
  
end Subroutine Param_LS ! Param_LS()
!****************************************************************************************************
Subroutine S200(l,l2,n,d,A,X,Y,func)
  implicit none
  external func
  integer SIZE
  parameter(SIZE=25)
  real*8 d,l2,A(l),X(SIZE),Y(SIZE) 
  integer l,n,j
  real*8 xx,yy
  l2 = 0.d0
  do j = 1, n
     xx = X(j)
     ! Obtain function
     call func(l,xx,yy,A)
     l2 = l2 + (Y(j) - yy) * (Y(j) - yy)
  end do
  d = dsqrt(l2 / dfloat(n - l))
  return
end Subroutine S200
!****************************************************************************************************
subroutine SBrody(l,x,y,A)
  implicit none
  double precision, external :: mygamma
  integer l
  real*8 x,y,A(l), w, alpha

  w = A(1)
  alpha = mygamma( (2.d0 + w)/(1.d0 + w) )
  !A(1) = w = Brody parameter
  y = (1 + w) * (x**w) * alpha**(1.d0 + w)  * exp(- (alpha * x)**(1.d0 + w)) 
  
  return
end subroutine SBrody
!****************************************************************************************************
subroutine MeanStDev(Data,n,Mean,Variance,StdDev)
  implicit none
  integer n, i
  double precision Data(n), Mean, Variance, StdDev
  
  Mean = 0.0                           ! compute mean
  DO i = 1, n
     Mean = Mean + Data(i)
  END DO
  Mean = Mean / dble(n)
  
  Variance = 0.0                       ! compute variance
  DO i = 1, n
     Variance = Variance + (Data(i) - Mean)**2
  END DO
  StdDev = Variance / dble(n - 1)
  StdDev   = SQRT(StdDev)            ! compute standard deviation
  
end subroutine MeanStDev
