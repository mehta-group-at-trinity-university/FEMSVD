program SVD
implicit none
        real*8 :: n,w,m,Rmax,res,temp,overlap,tempA,tempB,epsout,&
                Shift,alpha,r0diatom,etaOVERpi,xMin,xMax,yMin,yMax,RDerivDelt,RFirst,RLast,dscale,&
                e,eMin,eMax
        real*8, allocatable :: nodesR(:),weightsR(:),holder(:),Hsparse(:),&
                Bsparse(:),Pmu(:),Pnu(:),eigenVals(:),ress(:),&
                alphaR(:),alphaI(:),work(:),EVecs(:,:),VL(:,:),beta(:)
        real*8, allocatable :: dXr(:,:),TmatR(:,:),weightsDMR(:,:),S(:,:),indexOf(:,:),eigenVecs(:,:),&
                iPsi(:,:),jPsi(:,:),H(:,:),B(:,:)
        real*8, allocatable :: Uad(:,:,:),Psi(:,:,:),O(:,:,:,:)
        integer :: i,j,k,x,sR,sRc,row,mu,nu,inu,jmu,rCountA,numN0A,rCountB,numN0B,probSize,&
                xDim,yDim,ncv,matrixdim,HalfBandWidth,&
                NumStates,PsiFlag,CouplingFlag,LegPoints,order,Left,Right,Bottom,Top,xNumPoints,yNumPoints,&
                fparam(64),info,loop,m0,mm,&
                Tstart,Tend,rate
        integer, allocatable :: Hrow(:),Hcol(:),Brow(:),Bcol(:)
        character*64 LegendreFile

        open(unit=1,file="nodesAndWeights10.dat")
        !open(unit=2,file="nodesAndWeights20.dat")
        !open(unit=3,file="nodesAndWeights30.dat")
        !open(unit=4,file="nodesAndWeights40.dat")
        !open(unit=5,file="nodesAndWeights50.dat")
        !open(unit=6,file="nwTEMP.dat")
        !open(unit=7,file="nodesAndWeights70.dat")
        !open(unit=8,file="nodesAndWeights80.dat")
        !open(unit=9,file="nodesAndWeights90.dat")
        !open(unit=10,file="nodesAndWeights100.dat")

        open(unit=101,file="eigenVals-Vecs-10x10x10.dat")
        !open(unit=200,file="adiabaticPhi.dat")
        open(unit=201,file="adiabaticPhi-10x10x10.dat")
        open(unit=300,file="feastVSdggev-10x10x10.dat")
        !open(unit=400,file="overlap100x400.dat")
        open(unit=401,file="overlap10x100.dat")

        call system_clock(Tstart)

        !Create the Gauss Lobatto DVR functions X for the R coordinate
        x=1
        read(x,*) sR
        Rmax=10d0
        sRc=sR-2
        
        allocate(nodesR(sR),weightsR(sR),dXr(sRc,sR),holder(sR),weightsDMR(sR,sR),TmatR(src,src))

        !!!!!!
        !Create the Gauss Lobatto DVR functions X for the R coordinate and the Kmat matrix
        do i=1,sR
                read(x,*) n,w
                nodesR(i)=n
                weightsR(i)=w
        end do
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
        TmatR = (0.5d0)*MATMUL(dXr,MATMUL(weightsDMR,TRANSPOSE(dXr)))
        !!!!!!
        
        !!!!!!
        !obtain the adiabatic eigenvalues and eigenvectors from Kelly's code
        NumStates = 100
        LegendreFile='LegendreTEMP.dat'
        LegPoints=10
        Shift=-30000d0
        order=5
        Left=0
        Right=0
        Bottom=2
        Top=0
        alpha=1d0
        m=1d0
        xNumPoints=10
        xMin=0d0
        xMax=1.570796326794897
        yNumPoints=10
        yMin=0d0
        yMax =1.570796326794897*2
        RDerivDelt=0.0001d0
        RFirst=nodesR(2)
        RLast=nodesR(sR-1)
        dscale=10000d0
        
        !ncv=2*NumStates
        xDim=xNumPoints+Order-3
        yDim=yNumPoints+Order-3+1
        MatrixDim=xDim*yDim
        HalfBandWidth=yDim*Order+Order
        print *, NUmStates
        print *, HalfBandWidth, MatrixDim
        allocate(Uad(sRc,numStates,2),Psi(sRc,MatrixDim,numStates),S(HalfBandWidth+1,MatrixDim),&
                indexOf(sRc,NumStates),Pnu(numStates),Pmu(numStates),O(sRc,NumStates,sRc,Numstates),&
                iPsi(MatrixDim,numStates),jPsi(MatrixDim,numStates))
        call adiabaticSolver(NumStates,0,0,LegendreFile,LegPoints,Shift,Order,Left,Right,top,Bottom,alpha,m,&
                xNUmPoints,xMin,xMax,yNumPoints,yMin,yMax,sRc,RDerivDelt,RFirst,RLast,dscale,&
                nodesR(2:sR-1),Uad,Psi,numstates,MatrixDim,S,HalfBandWidth+1)
        !!!!!
        !S=0d0
        !do i=1,HalfBandWidth+1
        !        read(201,*) (S(i,j),j=1,MatrixDim)
        !end do
        !do i=1,sRc
        !        read(201,*) (Uad(i,j,1),j=1,numstates)
        !        do j=1,MatrixDim
        !                read(201,*) (Psi(i,j,k),k=1,numstates)
        !        end do
        !end do

        print *,"here"
        probSize=sRc*numStates
        
        eMin=0d0
        do i=1,sRc
                e=Uad(i,1,1)
                if (e<eMin) then
                        eMin=e
                end if
        end do
        eMax=Uad(sRc,1,1)

        !!!!!!
        !creat the complete H matrix in sparse formation
        row=0
        do i=1,sRc
                do j=1,numStates
                        row=row+1
                        indexOf(i,j)=row
                end do
        end do

        !!!!!!!!
        O=0d0
        do i=1,sRc
                print *,i
                do j=1,sRc
                        do k=1,MatrixDim
                                do mu=1,numStates
                                        iPsi(k,mu)=Psi(i,k,mu)
                                        jPsi(k,mu)=Psi(j,k,mu)
                                end do
                        end do
                        call CalcOMatrix(Numstates,HalfBandWidth,MatrixDim,iPsi,jPsi,S,O,sRc,i,j)
                end do
        end do
        !do i=1,sRc
        !        print *, i
        !        do j=1,sRc
        !                do nu=1,numStates
        !                        do mu=1,numStates
        !                                read(401,*)O(i,nu,j,mu)
        !                        end do
        !                end do
        !        end do
        !end do
        !!!!!!!!        

        print *, 'Got S'

        numN0A=0
        rCountA=1
        numN0B=0
        rCountB=1
        allocate(H(probsize,probsize),B(probsize,probsize),alphaR(probsize),alphaI(probsize),beta(probsize),VL(probsize,probsize),EVecs(probsize,probsize),work(8*probsize))
        H=0d0
        B=0d0
        do x=1,2
                print *,x
                do i=1,sRc
                        do nu=1,numStates
                                do j=1,sRc
                                        do mu=1,numStates
                                                inu=indexOf(i,nu)
                                                jmu=indexOf(j,mu)
                                                tempA = 0d0
                                                tempB = 0d0
                                                !Pnu=Psi(i,nu,1:numStates)
                                                !Pmu=Psi(j,mu,1:numStates)
                                                Overlap=O(i,nu,j,mu)
                                                if((overlap /= 0d0).and.(TmatR(i,j) /= 0d0)) then
                                                        tempA=tempA+TmatR(i,j)*overlap
                                                end if
                                                if((i==j).and.(nu==mu).and.(Uad(i,nu,1) /= 0d0)) then
                                                        tempA=tempA+Uad(i,nu,1)
                                                end if
                                                if((i==j).and.(overlap /= 0d0)) then
                                                        tempB=tempB+overlap*nodesR(i+1)**2
                                                end if
                                                if(tempA /= 0) then
                                                        numN0A=numN0A+1
                                                        rCountA=rCountA+1
                                                        if(x==2) then
                                                                Hsparse(numN0A)=tempA
                                                                Hcol(numN0A)=jmu
                                                                Hrow(inu+1)=Hrow(inu+1)+1
                                                        end if
                                                end if
                                                if(tempB /= 0) then
                                                        numN0B=numN0B+1
                                                        rCountB=rCountB+1
                                                        if(x==2) then
                                                                Bsparse(numN0B)=tempB
                                                                Bcol(numN0B)=jmu
                                                                Brow(inu+1)=Brow(inu+1)+1
                                                        end if
                                                end if
                                                if (x==2) then
                                                        H(inu,jmu)=TmatR(i,j)*overlap+Uad(i,nu,1)
                                                        B(inu,jmu)=overlap*nodesR(i+1)**2
                                                end if
                                        end do
                                end do
                        end do
                end do
                if(x==1) then
                        allocate(Hsparse(1:numN0A),Hcol(numN0A),Hrow(probSize+1))
                        Hrow=0
                        Hrow(1)=1
                        numN0A=0
                        rCountA=1
                        allocate(Bsparse(1:numN0B),Bcol(numN0B),Brow(probSize+1))
                        Brow=0
                        Brow(1)=1
                        numN0B=0
                        rCountB=1
                end if
        end do
        do i=2,probSize+1
                Hrow(i)=Hrow(i)+Hrow(i-1)
                Brow(i)=Brow(i)+Brow(i-1)
        end do
        write(*,*) "If",size(Hsparse),"=",Hrow(probSize+1)-1,"=",size(Hcol),"Then Good"
        write(*,*) "If",size(Bsparse),"=",Brow(probSize+1)-1,"=",size(Bcol),"Then Good"
        do i=1,numN0B
                if(Bsparse(i)<=0d0) then
                         !print *, "-B"
                end if
        end do

        m0=1000
        allocate(EigenVals(m0),EigenVecs(probSize,m0),ress(m0*2)) 
        call feastinit(fparam)
        fparam(1)=1
        call dfeast_scsrgv('F',probSize,Hsparse,Hrow,Hcol,Bsparse,Brow,Bcol,&
                fparam,epsout,loop,eMin,eMax,m0,EigenVals,EigenVecs,mm,ress,info)
        
        call dggev('N','V',probsize,H,probsize,B,probsize,alphaR,alphaI,Beta,VL,probsize,EVecs,probsize,work,8*probsize,INFO)
        call system_clock(Tend,rate)
        print *,info
        print *,(Tend-Tstart)/rate
        write(101,*) (Tend-Tstart)/rate, mm, eMin, eMax
        write(300,*) (Tend-Tstart)/rate, mm, eMin, eMax
        
        write(101,*) (EigenVals(j),j=1,mm)
        do i=1,probSize
                write(101,*) (EigenVecs(i,j),j=1,mm)
        end do
        do i=1,probsize
                write(300,*) alphaR(i)/beta(i)
        end do
        
        !close(1)
        !close(2)
        !close(3)
        !close(4)
        !close(5)
        close(6)
        !close(10)
        !close(15)
        !close(100)
        close(101)
        !close(102)
        !close(200)
        close(201)
        close(300)
        !close(400)
        close(401)
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

subroutine CalcOMatrix(NumStates,HalfBandWidth,MatrixDim,iPsi,jPsi,S,O,sRc,i,j)
implicit none
        integer NumStates,HalfBandWidth,MatrixDIm,i,j,k,R,sRc,nu,mu
        real*8 ddot
        double precision iPsi(MatrixDim,NumStates),jPsi(MatrixDim,numStates),S(HalfBandWidth+1,MatrixDim),&
                O(sRc,NumStates,sRc,NumStates),TempPsi(MatrixDim)

        do nu=1,numStates
                call dsbmv('U',MatrixDim,Halfbandwidth,1.0d0,S,Halfbandwidth+1,iPsi(1,nu),1,0.0d0,tempPsi,1)
                do mu=1,numStates
                        O(i,nu,j,mu) = ddot(MatrixDim,TempPsi,1,jPsi(1,mu),1)
                        write(401,*) O(i,nu,j,mu)
                end do
        end do
end subroutine CalcOMatrix
