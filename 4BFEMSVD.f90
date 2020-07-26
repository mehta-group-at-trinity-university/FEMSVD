module datastructures
  type irix
     integer ir, ix
     logical bridge
  end type irix
end module datastructures
program SVD
  use datastructures
  implicit none
  real*8 n,w,m,Rmin,Rmax,res,temp,overlap,tempA,tempB,epsout,res1,res2
  real*8 Shift,alpha,r0diatom,etaOVERpi,xMin,xMax,yMin,yMax,RDerivDelt,RFirst,RLast,V2Depth
  real*8 e,eMin,eMax,m1,m2,m3,m4,mu4,massarray(4)
  real*8, allocatable :: nodesR(:),weightsR(:),nds(:),wts(:),wmn(:,:),xmn(:,:)
  real*8, allocatable :: eigenVals(:),RSectors(:)
  real*8, allocatable :: EVals(:),EVecs(:,:),beta(:)
  real*8, allocatable :: dXrAll(:,:),dXr(:,:),fp(:,:,:),TmatTemp(:,:),TmatR(:,:),weightsDMR(:,:),S(:,:),eigenVecs(:,:)
  real*8, allocatable :: iPsi(:,:),jPsi(:,:),H(:,:),B(:,:)
  real*8, allocatable :: Uad(:,:,:),Psi(:,:,:),O(:,:), sectorT(:,:,:),wtemp(:,:,:), AllT(:,:,:,:)
  integer LP,NRSectors
  integer i,j,k,l,ivec,ifilenodes,ix,sR,sRc,row,mu,nu,inu,jmu,numN0A,numN0B,probSize,count
  integer xDim,yDim,ncv,matrixdim,HalfBandWidth,set
  integer NumStates,PsiFlag,CouplingFlag,LegPoints,order,Left,Right,Bottom,Top,xNumPoints,yNumPoints
  integer fparam(64),info,loop,Tstart,Tend,rate,iRS,jRS,in,im,jn,NRGridPoints
  integer, allocatable :: Hrow(:),Hcol(:),Brow(:),Bcol(:),indexOf(:,:), lindex(:,:)
  logical bridgebra, bridgeket
  type(irix), allocatable :: rxindex(:)
  character*64 LegendreFile
  character*64 Outputfile
  CHARACTER*64 LobattoFile
  character*3 wfnfile
  character*3 wfnsuffix
  character*64 wfnout
  integer iwfn, l1, l2
  double precision, parameter :: Pi = 3.1415926535897932384626433832795
  integer, external :: kdelta
  read(5,*)
  read(5,*) xNumPoints, yNumPoints
  read(5,*)
  read(5,*)
  read(5,*) NRSectors, LP, NumStates, Order
  read(5,*)
  read(5,*)
  read(5,*) V2Depth, Rmin, Rmax, alpha
  read(5,*)
  read(5,*)
  read(5,*) m1,m2,m3,m4
  read(5,*)
  read(5,*)
  read(5,*) xMin, xMax, yMin, yMax
  read(5,*)
  read(5,*)
  read(5,*) Left, Right, Bottom, Top

  xMin = xMin*Pi
  xMax = xMax*Pi
  yMin = yMin*Pi
  yMax = yMax*Pi
  write(6,*) "Bounds = ", xMin, xMax, yMin, yMax
  write(6,*) "V2Depth = ", V2Depth

  !mu4=(m1*m2*m3*m4/(m1+m2+m3+m4))**(1.0d0/3.0d0)
  mu4=(m1*m2*m3*m4)**(0.25d0)
  write(6,*) "mu4 = ", mu4
  massarray = (/m1,m2,m3,m4/)
  Outputfile = "Eigenvals.dat"
  open(unit=1,file=Outputfile)

  call system_clock(Tstart)

  NRGridPoints=NRSectors+1
  sR=(NRGridPoints-1)*LP-NRGridPoints+2
  sRc=sR-2
  write(6,*) "sR = ", sR
  write(6,*) "sRc = ", sRc
  allocate(rxindex(sR))
  allocate(nds(LP),wts(LP))
  allocate(nodesR(sR),weightsR(sR),dXrAll(sR,sR),weightsDMR(sR,sR))
  allocate(wmn(LP,NRSectors),xmn(LP,NRSectors),lindex(LP,NRSectors))
  allocate(RSectors(NRSectors+1))
  allocate(SectorT(LP,LP,NRSectors))
  allocate(fp(LP,LP,NRSectors))
  allocate(wtemp(LP,LP,NRSectors))

  call makeRSectors(RSectors,NRSectors, Rmin, Rmax)
  ! read in the Gauss-Lobatto nodes and weights.  Be sure LP is equal to one of
  ! the options in the input following file:
  LobattoFile = 'GaussLobatto.dat'
  call GetGaussLobattoFactors(LobattoFile,LP,nds,wts)
  
!  call rescale(sR,Rmin,Rmax,nds,wts,nodesR,weightsR)
  i = 0

  write(6,*) "setting all nodes:"
  weightsDMR=0d0
  wtemp=0d0
! count = 0
  do iRS = 1, NRSectors
     write(6,*) "iRS = ", iRS, RSectors(iRS) 
     call rescale(LP,RSectors(iRS),RSectors(iRS+1),nds,wts,xmn(:,iRS),wmn(:,iRS))
     do jn = 1, LP
        wtemp(jn,jn,iRS) = wmn(jn,iRS)
        if((iRS.lt.NRSectors).and.(jn.eq.LP)) cycle
!        count = count+1
        i = (iRS-1)*LP-iRS+1+jn
        lindex(jn,iRS) = i
        rxindex(i)%ir = iRS
        rxindex(i)%ix = jn
        rxindex(i)%bridge = .false.
        nodesR(i) = xmn(jn,iRS)
        weightsR(i) = wmn(jn,iRS)
 
        if((iRS.gt.1).and.(jn.eq.1)) then
!           count = count - 1
           weightsR(i) = wmn(LP,iRS-1) + wmn(1,iRS)
           rxindex(i)%bridge = .true.
        endif
        weightsDMR(i,i) = weightsR(i)
        !
!        write(6,*) iRS, jn, xmn(jn,iRS), wmn(jn,iRS), i,lindex(jn,iRS)
     enddo
  enddo
       write(6,*) "iRS = ", iRS, RSectors(iRS) 
  write(6,*) "all nodes and weights:"
  do i = 1, sR
     write(6,*) i, nodesR(i), weightsR(i), weightsDMR(i,i)
  enddo
  write(6,*) "do a test integral \int{x^2} from 0 to 10"
  res1=0d0
  do i = 1, sR
     res1 = res1 + weightsR(i)*nodesR(i)**2
  enddo
  write(6,*) "... result is:", res1
  write(6,*) "try it again with matrix math:"
  write(6,*) dot_product(nodesR,MATMUL(weightsDMR,nodesR))

  fp(:,:,:)=0d0

  do iRS = 1, NRSectors
     do in = 1, LP
        do jn = 1, LP
           call Calcfmnprime(xmn, wmn, iRS, iRS, in, jn, fp(jn,in,iRS), NRSectors, LP)
           
!           write(6,*) iRS, in, jn, fp(jn,in,iRS)
        enddo
!        fp(:,in,iRS) = fp(:,in,iRS)/sqrt(wmn(in,iRS))
!        write(6,*)
     enddo
     
     SectorT(:,:,iRS) = matmul(Transpose(fp(:,:,iRS)),matmul(wtemp(:,:,iRS),fp(:,:,iRS)))
!     call printmatrix(0.5d0*SectorT(:,:,iRS),LP,LP,6)
!     write(6,*)
  enddo
  
  ! now construct the basis derivatives (Suno eq. 15)
!  dXrAll(:,:)=0d0
  allocate(AllT(LP-1,LP-1,NRSectors,NRSectors))
  AllT=0d0
  allocate(TmatTemp(sr-1,sr-1))
  allocate(TmatR(src,src))
  TmatR=0d0
  TmatTemp=0d0
  do iRS = 1, NRSectors
     do jRS = 1, NRSectors
        do im = 1, LP-1
           do in = 1, LP-1
              if ( (in.gt.1) .and. (im.gt.1) ) then
                 AllT(im,in,iRS,jRS) = kdelta(iRS,jRS)*SectorT(im,in,iRS)*(wmn(im,iRS)*wmn(in,iRS))**(-0.5d0)
              else if( (im.eq.1).and.(iRS.gt.1).and.(in.gt.1) ) then
                 AllT(im,in,iRS,jRS) = (kdelta(iRS-1,jRS)*SectorT(LP,in,jRS) + kdelta(iRS,jRS)*SectorT(1,in,jRS)) &
                      *(wmn(in,jRS)*(wmn(LP,iRS-1)+wmn(1,jRS)))**(-0.5d0)
              else if( (in.eq.1) .and. (jRS.gt.1) .and. (im.gt.1) ) then
                 AllT(im,in,iRS,jRS) = (kdelta(iRS,jRS-1)*SectorT(im,LP,iRS) + kdelta(iRS,jRS)*SectorT(im,1,iRS)) &
                      *(wmn(im,iRS)*(wmn(LP,jRS-1)+wmn(1,iRS)))**(-0.5d0)
              else if ( (im.eq.1) .and. (in.eq.1) .and. (iRS.gt.1) .and. (jRS.gt.1)) then
                 AllT(im,in,iRS,jRS) = kdelta(iRS,jRS)*SectorT(LP,LP,iRS-1)
                 AllT(im,in,iRS,jRS) = AllT(im,in,iRS,jRS) + kdelta(iRS-1,jRS)*SectorT(LP,1,iRS-1)
                 AllT(im,in,iRS,jRS) = AllT(im,in,iRS,jRS) + kdelta(iRS,jRS-1)*SectorT(1,LP,iRS)
                 AllT(im,in,iRS,jRS) = AllT(im,in,iRS,jRS) + kdelta(iRS,jRS)*SectorT(1,1,iRS)
                 AllT(im,in,iRS,jRS) = AllT(im,in,iRS,jRS)/sqrt(wmn(LP,iRS-1)+wmn(1,iRS))/sqrt(wmn(LP,jRS-1)+wmn(1,jRS))
              endif
              TmatTemp(lindex(im,iRS),lindex(in,jRS)) = AllT(im,in,iRS,jRS)
              
           enddo
        enddo
     enddo
  enddo
  TmatR(:,:) = 0.5d0*TmatTemp(2:sr-1,2:sr-1)
  deallocate(TmatTemp)
!!$  allocate(dXr(sR,sRc))
!!$  dXr(:,:) = 0d0
!!$  do i = 1, sRc
!!$     dXr(:,i) = dXrAll(:,i+1)
!!$  enddo
!!$  TmatR=0d0
!!$  TmatR = (0.5d0/mu4)*MATMUL(Transpose(dXr),MATMUL(weightsDMR,dXr))

  !TmatR = (0.5d0/mu4)*MATMUL(dXr,MATMUL(weightsDMR,TRANSPOSE(dXr)))
  !  write(6,*) "Kinetic energy matrix:"
  
  !call printmatrix(TmatR,sRc,sRc,6)
  
  !  call testbasis(nodesR, weightsR, sr)
  call test1D(TmatR,sr,weightsR,nodesR)
  

  !obtain the adiabatic eigenvalues and eigenvectors 
  LegendreFile='Legendre.dat'
  LegPoints=10

  ! Set the "shift" for the DSBAND solver to -5*depth
  Shift = 0.d0!-5d0*V2Depth

  RDerivDelt=0.0001d0
  RFirst=nodesR(2)
  RLast=nodesR(sR-1)

  !ncv=2*NumStates
  xDim=xNumPoints+Order-3
  yDim=yNumPoints+Order-3+1
  MatrixDim=xDim*yDim
  HalfBandWidth=yDim*Order+Order
  print *, "MatDIM: ", MatrixDim
  probSize=sRc*numStates
  PsiFlag=0

  write(6,*) NUmStates
  write(6,*) HalfBandWidth, MatrixDim
  allocate(Uad(sRc,numStates,2),Psi(MatrixDim,numStates,sRc),S(HalfBandWidth+1,MatrixDim),&
       indexOf(sRc,NumStates),O(probSize,probSize),&
       iPsi(MatrixDim,numStates),jPsi(MatrixDim,numStates))
  call adiabaticSolver(NumStates,PsiFlag,0,LegendreFile,LegPoints,Shift,Order,Left,Right,top,Bottom,alpha,massarray,&
       xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,sRc,RDerivDelt,RFirst,RLast,V2Depth,&
       nodesR(2:sR-1),Uad,Psi,numstates,MatrixDim,S,HalfBandWidth+1)

  write(6,*)"About to calculate the overlap matrix..."
  ! create an array for the collective index
  row=1
  do i=1,sRc
     do j=1,numStates
        indexOf(i,j)=row
        row=row+1
     end do
  end do
  ! Begin calculation of overlap matrix
  O=0d0

  do i=1,sRc
     write(6,*) i
     iPsi(:,:) = Psi(:,:,i)
     do j=1,sRc
        jPsi(:,:) = Psi(:,:,j)
        call CalcOMatrix(Numstates,HalfBandWidth,MatrixDim,iPsi,jPsi,&
             S,O,sRc,i,j,indexOf,probSize,rxindex)
     end do
  end do

  allocate(H(probsize,probsize),EVals(probsize),EVecs(probsize,probsize))
  H=0d0

  write(6,*) "Constructing the Full Hamiltonian..."
  do i=1,sRc
     do nu=1,numStates
        do j=1,sRc
           do mu=1,numStates
              inu=indexOf(i,nu)
              jmu=indexOf(j,mu)

              !H(inu,jmu) = H(inu,jmu) + TmatR(i,j)*O(inu,jmu)
              H(inu,jmu) = TmatR(i,j)*O(inu,jmu)
              if(inu.eq.jmu) then
                 H(inu,jmu) = H(inu,jmu) + Uad(i,nu,1)
              endif

           end do
        end do
     end do
  end do

!!!!!!
  !Hamiltonian symmetry check
!!$     k=0
!!$     do i=1,sRc
!!$        do nu=1,numStates
!!$           do mu=1,numStates
!!$              inu=indexOf(i,nu)
!!$              jmu=indexOf(i,mu)
!!$              if(O(inu,jmu) .ne. O(inu,jmu) )  then
!!$                 k=k+1
!!$                 write(6,*) O(inu,jmu)
!!$              endif
!!$           end do
!!$        end do
!!$     end do
!!$     write(6,*) "DIFFS"
!!$     do i=1,sRc
!!$        do nu=1,numStates
!!$           !do j=1,sRc
!!$           !do mu=1,numStates
!!$           inu=indexOf(i,nu)
!!$           jmu=indexOf(i+1,nu)
!!$           write(6,*) O(inu,jmu)-O(jmu,inu)
!!$           !if (O(inu,jmu) .ne. O(jmu,inu)) k=k+1
!!$           !end do
!!$           !end do
!!$        end do
!!$     end do
!!$     write(6,*) k
!!$!!!!!

  write(6,*) "About to diagonalize the Hamiltonian. problem size = ", probsize
  call Mydsyev(H,probsize,EVals,EVecs)
  write(6,*) "done diagonalization"
  call system_clock(Tend,rate)
!  write(6,*) info
  write(6,*) dble((Tend-Tstart)/rate)/60.d0, "minutes"
  write(1,*) '#',dble((Tend-Tstart)/rate)/60.d0, numstates, Uad(sRc,1,1)

  do i=1,probsize
     write(1,*) EVals(i)
     !write(301,*) (EVecs(i,j),j=1,probsize)
  end do
  wfnfile = "wfn"
  do ivec = 1, 20
     write(wfnsuffix,'(I3)') 100+ivec
     wfnout = wfnfile // "-" // wfnsuffix // ".dat"
     write(6,*) wfnout
     !        open(9999,"wfnout")
     !       do iR = 1, nR
     !          write(file,10) nodesR(iR+1), (evecs(indexOf(iR,nu),ivec)/sqrt(weightsR(iR+1)), nu = 1, NumStates)
     !      enddo
     !        close(9999+ifile)
     call Plotwfn(ivec,indexOf,nodesR,weightsR,sr,sRc,NumStates,probsize,EVecs,wfnout)

  enddo
10 format(1P,2000e16.8)     
  close(1)
  !close(2)
  !close(3)
  !close(4)
  !close(5)
  !close(6)
  !close(8)
  close(10)
  close(20)
  !close(15)
  !close(100)
  !close(101)
  !close(102)
  !close(200)
  !close(201)
  !close(300)
  !close(301)
  !close(400)
  !close(401)


  deallocate(nodesR,weightsR,weightsDMR,TmatR,Uad,Psi,S,&
       H,EVals,EVecs,indexOf,O,iPsi,jPsi)



end program SVD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Plotwfn(ivec,indexOf,nodesR,weightsR,ngl,nR,NumStates,probsize,evecs,wfnfile)
  implicit none
  integer ivec,nu,iR,ngl,NumStates,nR,probsize,inu,indexOf(nR,NumStates)
  real*8 nodesR(ngl),weightsR(ngl),evecs(probsize,probsize)
  character*64 wfnfile
  !  write(6,*) file
  open(unit=8,file=wfnfile)
  do iR = 1, nR
     write(8,10) nodesR(iR+1), (evecs(indexOf(iR,nu),ivec)*sqrt(weightsR(iR+1)), nu = 1, NumStates)
  enddo
  close(8)
10 format(1P,2000e16.8)

end subroutine Plotwfn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rescale(n,a,b,nodes,weights,xn,wn)
  implicit none
  integer  n,i
  real*8  nodes(n), xn(n), wn(n), weights(n)
  real*8  a,b
  do i=1,n
     xn(i) = 0.5d0*(a+b) + nodes(i)*0.5d0*(b-a)
     wn(i) = 0.5d0*(b-a)*weights(i)
  end do
end subroutine rescale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find the derivative of the DVR basis function using the Eq. (24) in Suno, J. Chem. Phys 134, 064318 (2011)
! NRS = number of radial sectors
! ML = number of Lobatto points in each sector
! returns the derivative of the ith shape function at the jth point
subroutine Calcfmnprime(x, w, n, np, i, j, res, NRS, ML)
  implicit none
  integer n, np, i, j, k,NRS, ML, kd1, kdM
  double precision res, x(ML,NRS), w(ML,NRS),holder(ML)
!  call printmatrix(x,ML,NRS,6)
!  call printmatrix(w,ML,NRS,6)
  res = 0d0

  if(n.eq.np) then
     kd1 = 0
     kdM = 0
     if (i.eq.ML) kdM = 1
     if (j.eq.1) kd1 = 1
     holder = 0d0
     if(i.eq.j) then
        res = 0.5d0*(kdM - kd1)/w(i,n) 
     else
        do k = 1, ML
           if((k .ne. i).and.(k.ne.j)) then
              holder(k)=(x(j,n)-x(k,n))/(x(i,n)-x(k,n))
           else
              holder(k)=1d0
           endif
        enddo
        res=PRODUCT(holder)/((x(i,n)-x(j,n)))
     end if
  endif
  
end subroutine Calcfmnprime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function kdelta(i,j)
  integer i, j
  kdelta=0
  if(i.eq.j) kdelta=1
  
end function kdelta

!find derivitave of the DVR basis at the DVR points
! Returns the derivative of the ith basis function at the lth node
subroutine CalcDX(nodes,i,l,sz,weights,res)
  implicit none
  integer  i,sz,j,k,l
  real*8 nodes(sz),holder(sz),weights(sz)
  real*8 val,res
  
  res=0

  do k=1,sz
     holder=0d0
     if(k .ne. i) then
        do j=1,sz
           if(j.ne.i) holder(j)=(nodes(l)-nodes(j))/(nodes(i)-nodes(j))
        end do
        holder(i)=1d0
        holder(k)=1d0

        val=PRODUCT(holder)/((nodes(i)-nodes(k)))
        res = res + val
     end if
  end do
  res=res*(weights(i)**(-.5d0))
end subroutine CalcDX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine makeRSectors(RSectors,NRSectors, Rmin, Rmax)
  integer NRSectors,iR
  double precision RSectors(NRSectors+1), Rmin, Rmax

  do iR = 1, NRSectors+1
     RSectors(iR) = Rmin + (dble(iR-1)/dble(NRSectors)) * (Rmax - Rmin)
  enddo
end subroutine makeRSectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function CalcXfun(nodes,i,R,sz,weights)
  implicit none
  integer  i,sz,j,k
  real*8 nodes(sz),holder(sz),weights(sz)
  real*8 val,res,R
  res=0
  do j=1,i-1
     holder(j)=(R-nodes(j))/(nodes(i)-nodes(j))
  end do
  holder(i)=1d0
  do j=i+1,sz
     holder(j)=(R-nodes(j))/(nodes(i)-nodes(j))
  end do
  res=PRODUCT(holder)/sqrt(weights(i))
  CalcXfun=res
end function CalcXfun
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function CalcDXfun(nodes,i,R,sz,weights)
  implicit none
  integer i,sz,j,k
  real*8 nodes(sz),holder(sz),weights(sz)
  real*8 val,res,R
  res=0

  do k=1,sz
     holder=0d0
     if(k .ne. i) then
        do j=1,sz
           if(j.ne.i) holder(j)=(R-nodes(j))/(nodes(i)-nodes(j))
        end do
        holder(i)=1d0
        holder(k)=1d0

        val=PRODUCT(holder)/((nodes(i)-nodes(k)))
        res = res + val
     end if
  end do
  res=res*(weights(i)**(-.5d0))
  CalcDXfun=res
  
end function CalcDXfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testbasis(nodes, weights, sr)
  implicit none
  integer i,j,sr,sRc,np
  double precision, allocatable :: R(:)
  double precision nodes(sr), weights(sr)
  double precision, external :: CalcXfun, CalcDXfun
  np=1000
  allocate(R(np))
  do i = 1, np
     R(i) = nodes(1) + (dble(i)-1d0)/(dble(np)-1d0) * (nodes(sr) - nodes(1))
  enddo
  sRc = sr-2
  do i = 1, np
     write(2000,*) R(i), (CalcXfun(nodes,j,R(i),sr,weights), j=1,sr)
  enddo

  do i = 1, np
     write(2001,*) R(i), (CalcDXfun(nodes,j,R(i),sr,weights), j=1,sr)
  enddo

end subroutine testbasis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test1D(Tmat,sr,weights,nodes)
  implicit none
  integer sr, sRc, i, j, numvecs
  double precision Tmat(sr-2,sr-2),V(sr-2,sr-2)
  double precision weights(sr),nodes(sr)
  double precision evals(sr-2), evecs(sr-2,sr-2)
  numvecs=10
  sRc = sr-2
  V=0d0
  do i=1,sRc
     !V(i,i) = 6d0/(2d0*nodes(i+1)**2)
     V(i,i) = 0.5d0*nodes(i+1)**2
  enddo
  call Mydsyev(Tmat+V,sRc,evals,evecs)
  open(unit = 9, file = "HO-Evals.dat")
  do i = 1, sRc
     write(9,*) evals(i)
  enddo
  close(9)
  open(unit=10, file = "HO-Evecs.dat")
  write(10,*) nodes(1),(0d0, j=1,numvecs)
  do i = 1, sRc
     write(10,*) nodes(i+1), (evecs(i,j)/sqrt(weights(i+1)), j=1,numvecs)
  enddo
  write(10,*) nodes(sR),(0d0, j=1,numvecs)
  close(10)
end subroutine test1D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcOMatrix(NumStates,HalfBandWidth,MatrixDim,iPsi,jPsi,S,O,sRc,i,j,&
     indexOf,probSize,rxindex)
  use datastructures
  implicit none
  integer NumStates,HalfBandWidth,MatrixDIm,i,j,k,R,sRc,nu,mu
  integer indexOf(sRc,numStates),inu,jmu,probSize,irs,jrs
  real*8 ddot
  double precision iPsi(MatrixDim,NumStates),jPsi(MatrixDim,numStates),S(HalfBandWidth+1,MatrixDim)
  double precision O(probSize,probSize),TempPsi(MatrixDim)
  type(irix) rxindex(sRc + 2)
  
  do nu=1,numStates
     call dsbmv('U',MatrixDim,Halfbandwidth,1.0d0,S,Halfbandwidth+1,iPsi(:,nu),1,0.0d0,tempPsi,1)
     do mu=1,numStates
        irs = rxindex(i)%ir
        jrs = rxindex(j)%ir
        inu=indexOf(i,nu)
        jmu=indexOf(j,mu)
        if ( abs( irs - jrs  ) .le. 1) then
           O(inu,jmu) = ddot(MatrixDim,TempPsi,1,jPsi(1,mu),1)
           !write(6,*) irs, jrs, i, j, nu, mu, O(inu,jmu)
        else
           O(inu,jmu) = 0d0
        endif
        !               write(401,*) O(inu,jmu)
     end do
  end do
end subroutine CalcOMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE printmatrix(M,nr,nc,file)
  IMPLICIT NONE
  INTEGER nr,nc,file,j,k
  DOUBLE PRECISION M(nr,nc)
  
  DO j = 1,nr
     WRITE(file,30) (M(j,k), k = 1,nc)
  ENDDO
  
20 FORMAT(1P,100D16.8)
30 FORMAT(100D14.4)
END SUBROUTINE printmatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
