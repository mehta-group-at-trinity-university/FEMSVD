program ESpacings
use IFPORT
implicit none
        real*8, allocatable :: energies(:),eSpace(:),eNonDg(:)
        real *8 eSAvg,Elevel,epm,rnum, r1,r2
        integer*8 lenn, sz/8/
        integer j,i, eNum,Etot, numNg,eMin,eMax
        external comp
        integer*2 comp
        open(unit=10,file="./eigenVals-100x80x80x50x5-HHLLOdd.dat")
        open(unit=100,file="./eigenValsSorted-100x80x80x50x5-HHLOdd.dat")
        open(unit=200,file="./eigenValsPlotOdd")
        !open(unit=10,file="adiabaticEnergies-100x60x60x200x200-m12.dat")
        !open(unit=11,file="./eSpacings/eSpacings-100x80x80x200x50-m12Even.dat")
        !open(unit=12,file="./eSpacings/eSorted-100x80x80x200x50-m12Even.dat")
        read(10,*) Etot
        print *, Etot
        !Etot = 200
        allocate(energies(Etot),eSpace(Etot-1),eNonDg(Etot))

        do i=1,Etot
                read(10,*) energies(i)
        end do
        lenn = Etot
        !do i=1,Etot
        !        read(10,*) eSAvg, (energies(i),i=1,200)
        !end do
        print *, "read"
        call qsort(energies,lenn,sz,comp)
        print *, energies(1)
        
        do i=1,Etot
                write(100,*) energies(i)
        end do
        
        r1=0
        r2=10
        
        write(200,*) r1, (energies(j),j=1,Etot)
        write(200,*) r2, (energies(j),j=1,Etot)
        
        
        !!!
        !Remove degens
        
        !numNg=1
        !do i=2,Etot
        !        if(energies(i)-energies(i-1) .ge. 0.001d0) then
        !                eNonDg(numNg)=energies(i)
        !                numNg=numNg+1
        !        end if
        !end do
        !print *, Etot, numNg
        !energies=eNonDg
        !Etot=numNg
        !!!
        !50,12-119      100,12-248 
        !Elevel=-119.0d0
        !epm=10d0 
        !eMax=0
        !eMin=0
        !do i=1,Etot
!                write(12,*) energies(i)
        !        if((eMin.eq.0).and.(energies(i) .ge. Elevel-epm)) eMin=i
        !        if((energies(i) .gt. Elevel+epm).and.(eMax.eq.0)) eMax=i-1
        !end do
        !print *, eMax-eMin
        !eSAvg=0d0
        !do i=eMin,eMax-1
        !        eSpace(i)=energies(i+1)-energies(i)
        !        eSAvg=esAvg+eSpace(i)
        !end do
        !print *, "spaced"
        !eSAvg=eSAvg/(eMax-eMin-1)
        !do i=eMin,eMax-1
        !        eSpace(i)=eSpace(i)/eSavg
        !end do
        !call qsort(eSpace,eMax-eMin-1,sz,comp)
        !do i=eMin,eMax-1
!                write(11,*) eSpace(i)
        !end do
        !print *, "printed"
        close(10)
        close(100)
        close(200)
end program ESpacings

     
integer*2 function comp(a,b)
implicit none
        real*8 a,b
        if(a .lt. b) comp = -1
        if(a .eq. b) comp = 0
        if(a .gt. b) comp = 1
        return
end function comp
