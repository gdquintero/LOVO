program data_cubic
    implicit none

    real(kind=8) :: seed
    real(kind=8) :: a1,b1,a2,b2,t,r,ran1,ran2
    real(kind=8), dimension(4) :: xsol
    integer :: m,i

    m = 50
    a1 = -0.01d0
    b1 = 0.01d0
    a2 = 0.5d0
    b2 = 1.0d0

    seed = 123456.0d0
    xsol(:) = (/2.0d0,1.0d0,0.1d0,-0.1d0/)

    Open(Unit = 100, File = "data/cubic.txt", ACCESS = "SEQUENTIAL")

    write(100,*) m

    do i = 1, m
        t = (i-1) * 6.d0 / (m-1)
        ran1 = drand(seed)

        if (ran1 .lt. 0.1d0) then
            ran2 = drand(seed)

            r = a2 + (b2 - a2) * ran2

            if (ran2 .lt. 0.3d0) then
                write(100,*) t, poly(xsol,t,4) + r
            else
                write(100,*) t, poly(xsol,t,4) - r
            endif

        else
            r = a1 + (b1 - a1) * ran1
            write(100,*) t, poly(xsol,t,4) + r
        endif
    enddo


    contains

    function poly(x,t,n)
        implicit none

        real(kind=8) :: poly
        integer :: n
        real(kind=8) :: x(n),t

        poly = x(1) + x(2) * (t - 3.d0) + x(3) * ((t - 3.d0)**2) + x(4) * ((t - 3.d0)**3)

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