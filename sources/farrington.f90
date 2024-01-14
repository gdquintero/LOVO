program farrington
    implicit none

    integer :: samples

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)
    
    Open(Unit = 100, File = trim(pwd)//"/../data/seropositives.txt", Access = "SEQUENTIAL")

    read(100,*) samples

    print*, samples

    contains

    ! !*****************************************************************
    ! !*****************************************************************

    ! subroutine compute_sp(n,x,pdata,res)
    !     implicit none
    !     integer,       intent(in) :: n
    !     real(kind=8),  intent(in) :: x(n)
    !     real(kind=8),  intent(out) :: res

    !     type(pdata_type), intent(inout) :: pdata

    !     integer :: i,kflag

    !     pdata%sp_vector(:) = 0.0d0
    !     kflag = 2
    !     pdata%indices(:) = (/(i, i = 1, pdata%samples)/)

    !     do i = 1, pdata%samples
    !         call fi(n,x,i,pdata,pdata%sp_vector(i))
    !     end do

    !     ! Sorting
    !     call DSORT(pdata%sp_vector,pdata%indices,pdata%samples,kflag)

    !     ! Lovo function
    !     res = sum(pdata%sp_vector(1:pdata%lovo_order))

    ! end subroutine compute_sp

    ! !*****************************************************************
    ! !*****************************************************************

    ! subroutine compute_grad_sp(n,x,pdata,res)
    !     implicit none

    !     integer,       intent(in) :: n
    !     real(kind=8),  intent(in) :: x(n)
    !     real(kind=8),  intent(out) :: res(n)
    !     type(pdata_type), intent(in) :: pdata

    !     real(kind=8) :: gaux1,gaux2,a,b,c,ebt,ti
    !     integer :: i,j

    !     a = x(1)
    !     b = x(2)
    !     c = x(3)

    !     res(:) = 0.0d0

    !     do i = 1, pdata%lovo_order
    !         ti = pdata%t(int(pdata%indices(i)))
    !         ebt = exp(-b * ti)

    !         call model(n,x,int(pdata%indices(i)),pdata,gaux1)

    !         gaux1 = pdata%y(int(pdata%indices(i))) - gaux1
    !         gaux2 = exp((a / b) * ti * ebt + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) - c * ti)

    !         res(1) = res(1) + gaux1 * gaux2 * ((1.0d0 / b**2) * (ebt * (ti * b + 1.0d0) - 1.0d0))

    !         res(2) = res(2) + gaux1 * gaux2 * (ebt * ((-2.0d0 * a * ti / b**2) - ((a * ti**2) / b) &
    !                 - (2.0d0 * a / b**3) + (c / b**2) + (c * ti / b)) + (2.0d0 * a / b**3) - (c / b**2))

    !         res(3) = res(3) + gaux1 * gaux2 * ((1.0d0 / b) * (1.0d0 - ebt) - ti)
    !     enddo

    ! end subroutine compute_grad_sp

    ! !*****************************************************************
    ! !*****************************************************************

    ! subroutine model(n,x,i,pdata,res)
    !     implicit none 

    !     integer,        intent(in) :: n,i
    !     real(kind=8),   intent(in) :: x(n)
    !     real(kind=8),   intent(out) :: res
    !     real(kind=8) :: a,b,c,ti,ebt

    !     type(pdata_type), intent(in) :: pdata

    !     a = x(1)
    !     b = x(2)
    !     c = x(3)
    !     ti = pdata%t(i)
    !     ebt = exp(-b * ti)

    !     res = (a / b) * ti * ebt
    !     res = res + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) 
    !     res = 1.0d0 - exp(res - c * ti)

    ! end subroutine model

    ! !*****************************************************************
    ! !*****************************************************************

    ! subroutine fi(n,x,i,pdata,res)
    !     implicit none

    !     integer,        intent(in) :: n,i
    !     real(kind=8),   intent(in) :: x(n)
    !     real(kind=8),   intent(out) :: res

    !     type(pdata_type), intent(in) :: pdata

    !     call model(n,x,i,pdata,res)
    !     res = res - pdata%y(i)
    !     res = 0.5d0 * (res**2)

    ! end subroutine fi
    
end program farrington