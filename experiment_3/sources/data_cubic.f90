program data_cubic
    implicit none

    real(kind=8) :: seed1,seed2,a,b,r,ran1,ran2,noise,inf,sup,delta_t,outlier,percent_inf,percent_sup
    real(kind=8), dimension(4) :: xstar
    real(kind=8), allocatable :: t(:),y(:)
    integer :: m,i,allocerr

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    noise = 0.1d0
    percent_sup = 0.5d0
    percent_inf = 0.2d0

    m = 100
    a = -noise
    b = noise
    inf = -1.d0
    sup = 3.d0

    seed1 = 12345678912345.d0
    seed2 = 1234567891234.d0

    xstar(:) = (/1.d0,1.d0,-3.d0,1.d0/)

    Open(Unit = 100, File = trim(pwd)//"/../data/cubic.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 200, File = trim(pwd)//"/../data/cubic_latex.txt", ACCESS = "SEQUENTIAL")

    write(100,*) m

    allocate(t(m),y(m),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
    end if

    delta_t = (sup - inf) / real(m-1)

    do i = 1, m
        t(i) = inf + real(i - 1) * delta_t
        y(i) = poly(xstar,t(i),4)
    enddo

    do i = 1, m
        ran1 = drand(seed1)

        if (ran1 .le. 0.1d0) then
            ran2 = drand(seed2)   
            outlier = maxval(abs(y)) * (percent_inf + (percent_sup - percent_inf) * ran2)      

            if (ran2 .le. 0.2d0) then
                write(100,*) t(i), y(i) + outlier
                write(200,10) t(i), y(i) + outlier
            else
                write(100,*) t(i), y(i) - outlier
                write(200,10) t(i), y(i) - outlier
            endif

        else
            r = a + (b - a) * ran1
            write(100,*) t(i), y(i) + r
            write(200,10) t(i), y(i) + r
        endif
    enddo

    10 format (F6.3,1X,F6.3)

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