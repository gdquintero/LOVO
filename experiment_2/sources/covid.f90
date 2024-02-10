program main
    use sort

    implicit none

    type :: pdata_type
        integer :: counters(2) = 0
        integer :: samples,inf,sup,lovo_order,dim_Imin,n_train,n_test
        real(kind=8) :: sigma,theta
        real(kind=8), allocatable :: xtrial(:),xk(:),t(:),y(:),data(:,:),indices(:),sp_vector(:),&
        grad_sp(:),gp(:),lbnd(:),ubnd(:),pred(:),re(:)
        integer, allocatable :: outliers(:)
    end type pdata_type

    type(pdata_type), target :: pdata

    integer :: allocerr,i,n

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    n = 4
    
    Open(Unit = 10, File = trim(pwd)//"/../data/cubic.txt", Access = "SEQUENTIAL")

    read(10,*) pdata%samples

    pdata%n_train = pdata%samples - 20
    pdata%n_test = pdata%samples - pdata%n_train

    allocate(pdata%xtrial(n),pdata%xk(n),pdata%t(pdata%n_train),pdata%y(pdata%n_train),&
    pdata%indices(pdata%n_train),pdata%sp_vector(pdata%n_train),pdata%grad_sp(n),pdata%gp(n),&
    pdata%data(2,pdata%samples),pdata%pred(pdata%n_test),pdata%re(pdata%n_test),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
    end if

    do i = 1, pdata%samples
        read(10,*) pdata%data(:,i)
    enddo

    close(10)
  
    pdata%t(1:pdata%n_train) = pdata%data(1,1:pdata%n_train)
    pdata%y(1:pdata%n_train) = pdata%data(2,1:pdata%n_train)

    pdata%inf = 0
    pdata%sup = 10
 
    allocate(pdata%outliers(pdata%n_train*(pdata%sup-pdata%inf+1)),stat=allocerr)
 
    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if
 
    pdata%outliers(:) = 0


    contains

    subroutine lovo_algorithm(n,noutliers,outliers,pdata,fobj)
        implicit none
  
        integer, intent(in) :: n,noutliers
        integer, intent(inout) :: outliers(noutliers)
        real(kind=8), intent(out) :: fobj
        type(pdata_type), intent(inout) :: pdata
  
        real(kind=8) :: sigmin,sigmin_old,epsilon,fxk,fxtrial,alpha,gamma,termination
        integer :: iter_lovo,iter_sub_lovo,max_iter_lovo,max_iter_sub_lovo,k
  
        sigmin = 1.0d-1
        sigmin_old = sigmin
        gamma = 1.d+1
        epsilon = 1.0d-3
        alpha = 1.0d-8
        max_iter_lovo = 100000
        max_iter_sub_lovo = 100
        iter_lovo = 0
        iter_sub_lovo = 0
        pdata%lovo_order = pdata%n_train - noutliers
  
        pdata%xk(1:n) =  -1.d0
  
        call compute_sp(n,pdata%xk,pdata,fxk)      
  
        ! write(*,*) "--------------------------------------------------------"
        ! write(*,10) "#iter","#init","Sp(xstar)","Stop criteria","#Imin"
        ! 10 format (2X,A5,4X,A5,6X,A9,6X,A13,2X,A5)
        ! write(*,*) "--------------------------------------------------------"
  
        do
            iter_lovo = iter_lovo + 1
    
            call compute_grad_sp(n,pdata%xk,pdata,pdata%grad_sp)
    
            termination = norm2(pdata%grad_sp(1:n))
    
            ! write(*,20)  iter_lovo,iter_sub_lovo,fxk,termination,pdata%dim_Imin
            ! 20 format (I6,5X,I4,4X,ES14.6,3X,ES14.6,2X,I2)
    
            if (termination .lt. epsilon) exit
            if (iter_lovo .gt. max_iter_lovo) exit
            
            iter_sub_lovo = 1
            ! pdata%sigma = sigmin_old
            pdata%sigma = sigmin
            k = 1

            do 
                pdata%xtrial(:) = pdata%xk(:) - (1.d0 / pdata%sigma) * pdata%grad_sp(:)

                call compute_sp(n,pdata%xtrial,pdata,fxtrial)

                if (fxtrial .le. (fxk - alpha * norm2(pdata%xtrial(1:n-1) - pdata%xk(1:n-1))**2)) exit
                if (iter_sub_lovo .gt. max_iter_sub_lovo) exit

                k = k + 1

                ! if (k .eq. 2) then
                !     pdata%sigma = sigmin
                ! else
                !     pdata%sigma = gamma * pdata%sigma
                ! endif

                pdata%sigma = gamma * pdata%sigma
                iter_sub_lovo = iter_sub_lovo + 1

            enddo

            sigmin_old = pdata%sigma
  
            fxk = fxtrial
            pdata%xk(:) = pdata%xtrial(:)
            pdata%counters(2) = iter_sub_lovo + pdata%counters(2)
  
        enddo
  
        fobj = fxtrial
        pdata%counters(1) = iter_lovo
  
        ! write(*,*) "--------------------------------------------------------"

  
        outliers(:) = int(pdata%indices(pdata%n_train - noutliers + 1:))
        
  
    end subroutine lovo_algorithm

    !*****************************************************************
    !*****************************************************************
    
    subroutine rmsd(n,o,p,res)
        implicit none

        integer,        intent(in) :: n
        real(kind=8),   intent(in) :: o(n),p(n)
        real(kind=8),   intent(out) :: res
        integer :: i

        res = 0.d0

        do i = 1,n
            res = res + (o(i)-p(i))**2
        enddo

        res = sqrt(res / n)

    end subroutine rmsd

    !*****************************************************************
    !*****************************************************************

    subroutine relative_error(o,p,res)
        implicit none

        real(kind=8),   intent(in) :: o,p
        real(kind=8),   intent(out):: res

        res = abs(p - o) / abs(o)
        ! res = res * 100.d0
    
    end subroutine relative_error

    !*****************************************************************
    !*****************************************************************

    subroutine compute_sp(n,x,pdata,res)
        implicit none
        integer,       intent(in) :: n
        real(kind=8),  intent(in) :: x(n)
        real(kind=8),  intent(out) :: res

        type(pdata_type), intent(inout) :: pdata

        integer :: i,kflag

        pdata%sp_vector(:) = 0.0d0
        kflag = 2
        pdata%indices(:) = (/(i, i = 1, pdata%n_train)/)

        do i = 1, pdata%n_train
            call fi(n,x,i,pdata,pdata%sp_vector(i))
        end do

        ! Sorting
        call DSORT(pdata%sp_vector,pdata%indices,pdata%n_train,kflag)

        ! Lovo function
        res = sum(pdata%sp_vector(1:pdata%lovo_order))

        pdata%dim_Imin = 1

        if ( pdata%sp_vector(pdata%lovo_order) .eq. pdata%sp_vector(pdata%lovo_order + 1) ) then
            pdata%dim_Imin = 2
        endif

    end subroutine compute_sp

    !*****************************************************************
    !*****************************************************************

    subroutine compute_grad_sp(n,x,pdata,res)
        implicit none
  
        integer,       intent(in) :: n
        real(kind=8),  intent(in) :: x(n)
        real(kind=8),  intent(out) :: res(n)
        type(pdata_type), intent(in) :: pdata
  
        real(kind=8) :: gaux,ti,tm
        integer :: i
        
        res(:) = 0.0d0
  
        do i = 1, pdata%lovo_order
            ti = pdata%t(int(pdata%indices(i)))
            tm = pdata%t(pdata%samples)
            call model(n,x,int(pdata%indices(i)),pdata,gaux)
            gaux = gaux - pdata%y(int(pdata%indices(i)))
            
            res(1) = res(1) + gaux * (ti - tm)
            res(2) = res(2) + gaux * ((ti - tm)**2)
            res(3) = res(3) + gaux * ((ti - tm)**3)
        enddo
  
    end subroutine compute_grad_sp

    !*****************************************************************
    !*****************************************************************

    subroutine model(n,x,i,pdata,res)
        implicit none 

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res
        real(kind=8) :: ti,tm

        type(pdata_type), intent(in) :: pdata
   
        ti = pdata%t(i)
        tm = pdata%t(pdata%samples)

        res = pdata%y(pdata%samples) + x(1) * (ti - tm) + &
            x(2) * ((ti - tm)**2) + x(3) * ((ti - tm)**3)

    end subroutine model

    !*****************************************************************
    !*****************************************************************

    subroutine fi(n,x,i,pdata,res)
        implicit none

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res

        type(pdata_type), intent(in) :: pdata

        call model(n,x,i,pdata,res)
        res = res - pdata%y(i)
        res = 0.5d0 * (res**2)

    end subroutine fi
    
end program main