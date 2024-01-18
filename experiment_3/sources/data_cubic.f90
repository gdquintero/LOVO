program data_cubic
    implicit none

    real(kind=8) :: seed1,seed2
    real(kind=8) :: a1,b1,a2,b2,t,r,ran1,ran2
    real(kind=8), dimension(4) :: xsol
    integer :: m,i

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    m = 50
    a1 = -0.01d0
    b1 = 0.01d0
    a2 = 5.d0
    b2 = 10.d0

    seed1 = 123456.0d0
    seed2 = 1234.d0
    xsol(:) = (/1.d0,1.d0,-3.d0,1.d0/)

    Open(Unit = 100, File = trim(pwd)//"/../data/cubic.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 200, File = trim(pwd)//"/../data/cubic_latex.txt", ACCESS = "SEQUENTIAL")

    write(100,*) m

    do i = 1, m
        t = (i-1) * 4.d0 / (m-1)
        ran1 = drand(seed1)

        if (ran1 .lt. 0.1d0) then
            ran2 = drand(seed2)

            r = a2 + (b2 - a2) * ran2

            if (ran2 .lt. 0.2d0) then
                write(100,*) t, poly(xsol,t,4) + r
                write(200,10) t, poly(xsol,t,4) + r
            else
                write(100,*) t, poly(xsol,t,4) - r
                write(200,10) t, poly(xsol,t,4) - r
            endif

        else
            r = a1 + (b1 - a1) * ran1
            write(100,*) t, poly(xsol,t,4) + r
            write(200,10) t, poly(xsol,t,4) + r
        endif
    enddo

    ! do i = 1, m
    !     t = (i-1) * 4.d0 / (m-1)
    !     ran1 = drand(seed1)

    !     if (11 .le. i .and. i .le. 20) then
    !         write(100,*) t, 10.d0
    !         write(200,10) t, 10.d0
            
    !     else
    !         r = a1 + (b1 - a1) * ran1
    !         write(100,*) t, poly(xsol,t,4) + r
    !         write(200,10) t, poly(xsol,t,4) + r
    !     endif
        
    ! enddo

    10 format (F5.3,1X,F6.3)

    close(100)
    close(200)

    contains

    function poly(x,t,n)
        implicit none

        real(kind=8) :: poly
        integer :: n
        real(kind=8) :: x(n),t

        poly = x(1) + x(2) * t + x(3) * (t**2) + x(4) * (t**3)

    end function poly

    function drand(ix)

        implicit none
      
        ! This is the random number generator of Schrage:
        !
        ! L. Schrage, A more portable Fortran random number generator, ACM
        ! Transactions on Mathematical Software 5 (1979), 132-138.
      
        ! FUNCTION TYPE
        real(kind=8) :: drand
      
        ! SCALAR ARGUMENT
        real(kind=8), intent(inout) :: ix
      
        ! LOCAL ARRAYS
        real(kind=8) :: a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      
        data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/
      
        xhi= ix/b16
        xhi= xhi - dmod(xhi,1.d0)
        xalo= (ix-xhi*b16)*a
        leftlo= xalo/b16
        leftlo= leftlo - dmod(leftlo,1.d0)
        fhi= xhi*a + leftlo
        k= fhi/b15
        k= k - dmod(k,1.d0)
        ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
        if (ix.lt.0) ix= ix + p
        drand= ix*4.656612875d-10
      
        return
      
      end function drand

end program data_cubic