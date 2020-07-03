program SVD
  implicit none
  real*8 n,w,m,Rmax,res,temp,overlap,tempA,tempB,epsout, &
       Shift,alpha,r0diatom,etaOVERpi,xMin,xMax,yMin,yMax,RDerivDelt,RFirst,RLast,V2Depth,&
       e,eMin,eMax,m1,m2,m3,m4,mu4
  real*8, allocatable :: nodesR(:),weightsR(:),holder(:),Hsparse(:),&
       Bsparse(:),Pmu(:),Pnu(:),eigenVals(:),ress(:),&
       alphaR(:),alphaI(:),work(:),EVecs(:,:),VL(:,:),beta(:)
  real*8, allocatable :: dXr(:,:),TmatR(:,:),weightsDMR(:,:),S(:,:),eigenVecs(:,:),&
       iPsi(:,:),jPsi(:,:),H(:,:),B(:,:)
  real*8, allocatable :: Uad(:,:,:),Psi(:,:,:),O(:,:)
  integer LobattoPoints
  integer i,j,k,ifilenodes,ix,sR,sRc,row,mu,nu,inu,jmu,rCountA,numN0A,rCountB,numN0B,probSize, &
       xDim,yDim,ncv,matrixdim,HalfBandWidth,set,&
       NumStates,PsiFlag,CouplingFlag,LegPoints,order,Left,Right,Bottom,Top,xNumPoints,yNumPoints,&
       fparam(64),info,loop,m0,mm,&
       Tstart,Tend,rate
  integer, allocatable :: Hrow(:),Hcol(:),Brow(:),Bcol(:),indexOf(:,:)
  character*64 LegendreFile
  character* 64 filename
  CHARACTER*64 LobattoFile

  read(5,*)
  read(5,*) xNumPoints, yNumPoints
  read(5,*)
  read(5,*)
  read(5,*) LobattoPoints, NumStates
  read(5,*)
  read(5,*)
  read(5,*) V2Depth, Rmax
  
  !open(unit=101,file="eigenVals-Vecs-100x60x60x200.dat")
  !open(unit=200,file="adiabaticPhi.dat")
  !open(unit=201,file="adiabaticPhi-200x60x60x100x100-m12.dat")
  !open(unit=300,file="eigenVals-200x60x60x100x100-m12.dat")
  !open(unit=301,file="eigenVecs-200x60x60x100x100-m12.dat")
  !open(unit=400,file="overlap100x400.dat")
  !open(unit=401,file="overlap200x60x60x100x100-m12.dat")

  do set=1,1
     !if (set==1) filename = "eigenValsCplxAbs-HHLLEven.dat"
     ! naming convention = radial x xpoints x ypoints numchannels x depth
     if (set==1) filename = "Eigenvals.dat"!"eigenVals-100x100x100x20x5-HHHHOdd.dat"
     if (set==2) filename = "eigenVals-100x100x100x20x5-HHHHEven.dat"
     if (set==3) filename = "eigenVals-100x90x90x200x50-m12Even.dat"
     if (set==4) filename = "eigenVals-150x80x80x200x100-m12Even.dat"
     if (set==5) filename = "eigenVals-200x80x80x200x250-m12.dat"
     if (set==6) filename = "eigenVals-200x80x80x200x300-m12.dat"
     open(unit=1,file=filename)

     call system_clock(Tstart)

     sR=LobattoPoints
     sRc=sR-2

     !if(set==1 .OR. set==3 ) then
     !        m1=1d0
     !        m2=m1
     !        m3=m1
     !        m4=m1
     !else 
     m1=1d0
     m2=m1
     m3=1.d0
     m4=m3
     !end if
     mu4=(m1*m2*m3*m4/(m1+m2+m3+m4))**(1.0d0/3.0d0);

     allocate(nodesR(sR),weightsR(sR),dXr(sRc,sR),holder(sR),weightsDMR(sR,sR),TmatR(src,src))

!!!!!!
     !Create the Gauss Lobatto DVR functions X for the R coordinate and the Kmat matrix
!!$     do i=1,sR
!!$        read(ifilenodes,*) n,w
!!$        nodesR(i)=n
!!$        weightsR(i)=w
!!$     end do
     LobattoFile = 'GaussLobatto.dat'
     call GetGaussLobattoFactors(LobattoFile,LobattoPoints,nodesR,weightsR)
     
     call rescale(sR,0d0,Rmax,nodesR)
     do i=1,sRc
        do j=1,sR
           call calcDX(holder,nodesR,nodesR(j),i,sR,weightsR,res)
           dXr(i,j)=res
        end do
     end do
     weightsDMR=0d0
     do i=1,sR
        do j=1,sR
           if(i==j) then
              weightsDMR(i,j)=weightsR(i)*nodesR(i)**2d0
           end if
        end do
     end do
     TmatR=0d0
     TmatR = (0.5d0/mu4)*MATMUL(dXr,MATMUL(weightsDMR,TRANSPOSE(dXr)))
!!!!!!

!!!!!!
     !obtain the adiabatic eigenvalues and eigenvectors from Kelly's code

     NumStates = 5
     LegendreFile='Legendre.dat'
     LegPoints=10

     if(set==1) Shift=-5*40
     if(set==2) Shift=-5*5
     if(set==3) Shift=-5*50
     if(set==4) Shift=-5*100
     if(set==5) Shift=-5*250
     if(set==6) Shift=-5*300
     order=5
     Left=0
     Right=0 
     Bottom=2
     if(set==1) Top=0 !odd parity
     if(set==2) Top=1 !even parity
     alpha=1d0
     m=1d0
     if(set==1)xNumPoints = 80
     if(set==2)xNumPoints = 100
     if(set==2)xNumPoints = 90
     !xNumPoints=80
     xMin=0d0
     xMax=1.570796326794897
     if(set==1)yNumPoints = 80
     if(set==2)yNumPoints = 100
     if(set==3)yNumPoints = 90
     !yNumPoints=80
     yMin=0d0
     yMax =1.570796326794897
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

     print *, NUmStates
     print *, HalfBandWidth, MatrixDim
     allocate(Uad(sRc,numStates,2),Psi(sRc,MatrixDim,numStates),S(HalfBandWidth+1,MatrixDim),&
          indexOf(sRc,NumStates),Pnu(numStates),Pmu(numStates),O(probSize,probSize),&
          iPsi(MatrixDim,numStates),jPsi(MatrixDim,numStates))
     call adiabaticSolver(NumStates,1,0,LegendreFile,LegPoints,Shift,Order,Left,Right,top,Bottom,alpha,m,&
          xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,sRc,RDerivDelt,RFirst,RLast,V2Depth,&
          nodesR(2:sR-1),Uad,Psi,numstates,MatrixDim,S,HalfBandWidth+1,set)

!!!!!
     print *,"About to calculate the overlap matrix..."

!!!!!!
     ! create an array for the collective index
     row=1
     do i=1,sRc
        do j=1,numStates
           indexOf(i,j)=row
           row=row+1
        end do
     end do

!!!!!!!!
     ! Begin calculation of overlap matrix
     O=0d0
     do i=1,sRc
        print *,i
        do j=1,sRc
           iPsi(:,:) = Psi(i,:,:)
           jPsi(:,:) = Psi(j,:,:)
!!$           do k=1,MatrixDim
!!$              do mu=1,numStates
!!$                 iPsi(k,mu)=Psi(i,k,mu)
!!$                 jPsi(k,mu)=Psi(j,k,mu)
!!$              end do
!!$           end do
           call CalcOMatrix(Numstates,HalfBandWidth,MatrixDim,iPsi,jPsi,&
                S,O,sRc,i,j,indexOf,probSize)
        end do
     end do

     allocate(H(probsize,probsize),alphaR(probsize),alphaI(probsize), &
          beta(probsize),VL(probsize,probsize),EVecs(probsize,probsize),work(8*probsize))
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
!!$                 print *, O(inu,jmu)
!!$              endif
!!$           end do
!!$        end do
!!$     end do
!!$     print *, "DIFFS"
!!$     do i=1,sRc
!!$        do nu=1,numStates
!!$           !do j=1,sRc
!!$           !do mu=1,numStates
!!$           inu=indexOf(i,nu)
!!$           jmu=indexOf(i+1,nu)
!!$           print *, O(inu,jmu)-O(jmu,inu)
!!$           !if (O(inu,jmu) .ne. O(jmu,inu)) k=k+1
!!$           !end do
!!$           !end do
!!$        end do
!!$     end do
!!$     print *, k
!!$!!!!!

     m0=100
     write(6,*) "About to diagonalize the Hamiltonian..."
     
     call dgeev('N','V',probsize,H,probsize,alphaR,alphaI,VL,probsize,EVecs,probsize,work,8*probsize,INFO)
     call system_clock(Tend,rate)
     print *,info
     print *,(Tend-Tstart)/rate
     !write(101,*) (Tend-Tstart)/rate, mm, eMin, eMax
     write(1,*) (Tend-Tstart)/rate, numstates, Uad(sRc,1,1)

     do i=1,probsize
        write(1,*) alphaR(i)
        !write(301,*) (EVecs(i,j),j=1,probsize)
     end do

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


     deallocate(nodesR,weightsR,dXr,holder,weightsDMR,TmatR,Uad,Psi,S,&
          H,alphaR,alphaI,beta,VL,EVecs,work,indexOf,Pnu,Pmu,O,iPsi,jPsi)

  end do

end program SVD

subroutine rescale(n,a,b,nodes)
  implicit none
  integer :: n,i
  real*8 :: nodes(*)
  real*8 :: x1,xN,a,b
  x1= nodes(1)
  xN= nodes(n)
  do i=1,n
     if((xN-x1)*(b-a)==0) then
        nodes(i)=0d0
     else 
        nodes(i)=a+(nodes(i)-x1)/(xN-x1)*(b-a)
     end if
  end do
end subroutine rescale

!find derivitave of the nodes
subroutine CalcDX(holder,nodes,node,i,sz,weights,res)
  implicit none
  integer :: i,sz,j,k
  real*8 nodes(sz),holder(sz),weights(sz)
  real*8 val,node,res

  res=0
  do k=1,sz
     if(k /= i) then
        do j=1,sz
           if((nodes(i)-nodes(j))==0d0) then
              holder(j)=0d0
           else 
              holder(j)=(node-nodes(j))/(nodes(i)-nodes(j))
           end if
        end do
        holder(i)=1d0
        holder(k)=1d0
        if((nodes(i)-nodes(k))==0) then
           val=0d0
        else 
           val=PRODUCT(holder)/((nodes(i)-nodes(k)))
        end if
        res = res + val
     end if
  end do
  res=res*(weights(i)**(-.5d0))
end subroutine CalcDX

subroutine CalcOMatrix(NumStates,HalfBandWidth,MatrixDim,iPsi,jPsi,S,O,sRc,i,j,indexOf,probSize)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDIm,i,j,k,R,sRc,nu,mu,&
       indexOf(sRc,numStates),inu,jmu,probSize
  real*8 ddot
  double precision iPsi(MatrixDim,NumStates),jPsi(MatrixDim,numStates),S(HalfBandWidth+1,MatrixDim),&
       O(probSize,probSize),TempPsi(MatrixDim)
  
  do nu=1,numStates
     call dsbmv('U',MatrixDim,Halfbandwidth,1.0d0,S,Halfbandwidth+1,iPsi(:,nu),1,0.0d0,tempPsi,1)
     do mu=1,numStates
        inu=indexOf(i,nu)
        jmu=indexOf(j,mu)
        O(inu,jmu) = ddot(MatrixDim,TempPsi,1,jPsi(1,mu),1)
        !               write(401,*) O(inu,jmu)
     end do
  end do
end subroutine CalcOMatrix
