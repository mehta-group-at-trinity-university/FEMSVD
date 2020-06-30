program ESpacings
use IFPORT
implicit none
        real*8 energies(98,200),R(98), eSAvg
        integer*8 lenn/1000/, sz/8/
        integer j,i
        external comp
        integer*2 comp
        open(unit=10,file="adiabaticEnergies-100x60x50x200x5.dat")
        open(unit=11,file="r,adE.dat")
        do i=1,98
                read(10,*) R(i), (energies(i,j),j=1,200)
        end do
        do i=1,98
                do j=1,200
                       write(11,*) 
                end do
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
