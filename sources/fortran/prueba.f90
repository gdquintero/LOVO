program farrington
    implicit none 

    integer :: samples
    CHARACTER(*), PARAMETER :: fileplace = "/home/gustavo/github/LOVO/data/"


    OPEN(UNIT=10,FILE=fileplace//"seropositives.txt",ACCESS="SEQUENTIAL")
    read(10,*) samples

    print*, samples
end program farrington