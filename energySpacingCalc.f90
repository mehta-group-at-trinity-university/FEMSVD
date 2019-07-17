program ESpacings
use IFPORT
implicit none
        real*8 energies(1000),eSpace(999), eSAvg
        integer*8 lenn/1000/, sz/8/
        integer j,i
        external comp
        integer*2 comp
        open(unit=10,file="eigenVals-Vecs.dat")
        open(unit=11,file="eSpacings.dat")
        read(10,*)
        read(10,*) (energies(j),j=1,1000)
        call qsort(energies,lenn,sz,comp)
        do i=1,999
                eSpace(i)=energies(i+1)-energies(i)
        end do
        eSAvg=0d0
        do i=1,49
                eSAvg=eSAvg+eSpace(i)
        end do
        eSAvg=eSAvg/49
        write(11,*) energies(1)
        do i=1,49
                write(11,*)energies(i),eSpace(i), eSpace(i)/eSAvg
        end do
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
