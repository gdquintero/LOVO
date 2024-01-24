program main
    use sort

    implicit none

    type :: pdata_type
        integer :: counters(2) = 0
        integer :: samples,inf,sup,lovo_order,dim_Imin,n_train,n_test
        real(kind=8) :: sigma,theta
        real(kind=8), allocatable :: xstar(:),xtrial(:),xk(:),t(:),y(:),data(:,:),indices(:),sp_vector(:),&
        grad_sp(:),gp(:),lbnd(:),ubnd(:)
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

    allocate(pdata%xstar(n),pdata%xtrial(n),pdata%xk(n),pdata%t(pdata%n_train),pdata%y(pdata%n_train),&
    pdata%indices(pdata%n_train),pdata%sp_vector(pdata%n_train),pdata%grad_sp(n),pdata%gp(n),&
    pdata%data(2,pdata%n_train),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
    end if

    do i = 1, pdata%n_train
        read(10,*) pdata%data(:,i)
    enddo

    close(10)

    pdata%xstar(:) = (/1.d0,1.d0,-3.d0,1.d0/)
  
    pdata%t(:) = pdata%data(1,:)
    pdata%y(:) = pdata%data(2,:)

    pdata%inf = 0
    pdata%sup = 7
 
    allocate(pdata%outliers(pdata%n_train*(pdata%sup-pdata%inf+1)),stat=allocerr)
 
    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if
 
    pdata%outliers(:) = 0

    call mixed_test(n,pdata)

    Open(Unit = 10, File =trim(pwd)//"/../output/outliers.txt", ACCESS = "SEQUENTIAL")

    write(10,100) pdata%sup

    do i = 1, pdata%sup
        write(10,100) pdata%outliers(i)
    enddo

    100 format (I2)

    contains

    subroutine mixed_test(n,pdata)
        implicit none
  
        integer, intent(in) :: n
        integer :: noutliers
        real(kind=8) :: fobj,start,finish
        type(pdata_type), intent(inout) :: pdata

        Open(Unit = 100, File = trim(pwd)//"/../output/solution_cubic.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 200, File = trim(pwd)//"/../output/output_latex.txt", ACCESS = "SEQUENTIAL")
  
        do noutliers = pdata%inf, pdata%sup
            ! write(*,*) "LOVO Algorithm for Measles:"

            call cpu_time(start)
            call lovo_algorithm(n,noutliers,pdata%outliers,pdata,fobj)
            call cpu_time(finish)

            write(100,1000) pdata%xk(1),pdata%xk(2),pdata%xk(3),pdata%xk(4)
            write(200,1100) pdata%xk(1),pdata%xk(2),pdata%xk(3),pdata%xk(4),fobj,norm2(pdata%xk-pdata%xstar),&
            maxval(abs(pdata%xk-pdata%xstar)),pdata%counters(1),pdata%counters(2)
  
            pdata%counters(:) = 0
           
        enddo
  
        Open(Unit = 500, File = trim(pwd)//"/../output/num_mixed_test.txt", ACCESS = "SEQUENTIAL")
        write(500,1200) pdata%inf
        write(500,1200) pdata%sup
  
        1000 format (ES13.6,1X,ES13.6,1X,ES13.6,1X,ES13.6)
        1200 format (I2)
        1100 format (F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,I4,1X,I4)

        close(100)
        close(500)

  
    end subroutine mixed_test

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
  
        write(*,*) "--------------------------------------------------------"
        write(*,10) "#iter","#init","Sp(xstar)","Stop criteria","#Imin"
        10 format (2X,A5,4X,A5,6X,A9,6X,A13,2X,A5)
        write(*,*) "--------------------------------------------------------"
  
        do
            iter_lovo = iter_lovo + 1
    
            call compute_grad_sp(n,pdata%xk,pdata,pdata%grad_sp)
    
            termination = norm2(pdata%grad_sp(1:n))
    
            write(*,20)  iter_lovo,iter_sub_lovo,fxk,termination,pdata%dim_Imin
            20 format (I6,5X,I4,4X,ES14.6,3X,ES14.6,2X,I2)
    
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
  
        write(*,*) "--------------------------------------------------------"

  
        outliers(:) = int(pdata%indices(pdata%n_train - noutliers + 1:))
        
  
    end subroutine lovo_algorithm


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
  
        real(kind=8) :: gaux,ti
        integer :: i
        
        res(:) = 0.0d0
  
        do i = 1, pdata%lovo_order
           ti = pdata%t(int(pdata%indices(i)))
           call model(n,x,int(pdata%indices(i)),pdata,gaux)
           gaux = gaux - pdata%y(int(pdata%indices(i)))
  
           res(1) = res(1) + gaux
           res(2) = res(2) + gaux * ti
           res(3) = res(3) + gaux * (ti**2)
           res(4) = res(4) + gaux * (ti**3)
        enddo
  
    end subroutine compute_grad_sp

    !*****************************************************************
    !*****************************************************************

    subroutine model(n,x,i,pdata,res)
        implicit none 

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res
        real(kind=8) :: ti

        type(pdata_type), intent(in) :: pdata
   
        ti = pdata%t(i)

        res = x(1) + (x(2) * ti) + (x(3) * (ti**2)) + (x(4) * (ti**3))

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