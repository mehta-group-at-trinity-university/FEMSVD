program ESpacings
use IFPORT
implicit none
        real*8, allocatable :: energies(:),eSpace(:),eNonDg(:)
        real *8 eSAvg,radi
        integer*8 lenn, sz/8/
        integer j,i, eNum,Etot, numNg,eMin,eMax, rPoints, levels
        external comp
        integer*2 comp
        open(unit=10,file="adiabaticEnergies-150x60x60x200x100-m12.dat")
        !open(unit=10,file="adiabaticEnergies-100x60x60x200x200-m12.dat")
        open(unit=11,file="ADeSpacings44-150x60x60x200x100-m12.dat")
        open(unit=12,file="ADeSorted37-150x60x60x200x100-m12.dat")
        !read(10,*) levels
        levels = 200
        print *, levels
        
        radi=44
        
        !Etot = 200
        allocate(energies(levels),eSpace(levels-1),eNonDg(Etot))
        do i=1,radi-1
                read(10,*)
        end do
        read(10,*) radi, (energies(j),j=1,levels)
        lenn = levels
        !do i=1,Etot
        !        read(10,*) eSAvg, (energies(i),i=1,200)
        !end do
        print *, "read"
        call qsort(energies,lenn,sz,comp)
        print *, energies(1)
        
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
        
        !eMin=0
        !do i=1,Etot
        !        write(12,*) energies(i)
        !        if((eMin.eq.0).and.(energies(i) .ge. -242.5d0)) eMin=i
        !        if((energies(i) .gt. -222.5d0).and.(eMax.eq.0)) eMax=i-1
        !end do
        !print *, eMax-eMin
        eSAvg=0d0
        do i=1,levels-1
                eSpace(i)=energies(i+1)-energies(i)
                eSAvg=esAvg+eSpace(i)
        end do
        print *, "spaced", eSAvg
        eSAvg=eSAvg/(levels-1)
        do i=1,levels-1
                eSpace(i)=eSpace(i)/eSavg
        end do
        call qsort(eSpace,levels-1,sz,comp)
        do i=1,levels-1
                write(11,*) eSpace(i)
        end do
        print *, "printed"
        close(10)
        close(11)
end program ESpacings

     
integer*2 function comp(a,b)
implicit none
        real*8 a,b
        if(a .lt. b) comp = -1
        if(a .eq. b) comp = 0
        if(a .gt. b) comp = 1
        return
end function comp
