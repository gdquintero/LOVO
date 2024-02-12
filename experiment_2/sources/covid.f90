program main
    use sort

    implicit none

    type :: pdata_type
        integer :: counters(3) = 0
        integer :: lovo_order,n_train,n_test,noutliers,dim_Imin
        real(kind=8) :: sigma,fobj
        real(kind=8), allocatable :: xtrial(:),xk(:),t(:),y(:),t_test(:),y_test(:),indices(:),&
        sp_vector(:),grad_sp(:)
        integer, allocatable :: outliers(:)
    end type pdata_type

    type(pdata_type), target :: pdata

    integer :: allocerr,i,n

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

   ! Number of variables
    n = 3
 
    Open(Unit = 100, File = trim(pwd)//"/../data/covid_train.txt", Access = "SEQUENTIAL")
    Open(Unit = 200, File = trim(pwd)//"/../data/covid_test.txt", Access = "SEQUENTIAL")
    
    read(100,*) pdata%n_train
    read(200,*) pdata%n_test

    pdata%noutliers = 2*int(dble(pdata%n_train) / 7.0d0)

 
    allocate(pdata%t(pdata%n_train),pdata%y(pdata%n_train),pdata%y_test(pdata%n_test),pdata%t_test(pdata%n_test),&
    pdata%xtrial(n),pdata%xk(n),pdata%grad_sp(n),pdata%indices(pdata%n_train),pdata%sp_vector(pdata%n_train),&
    pdata%outliers(pdata%noutliers),stat=allocerr)
 
    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if

    do i = 1, pdata%n_train
        read(100,*) pdata%y(i)
    enddo

    do i = 1, pdata%n_test
        read(200,*) pdata%y_test(i)
    enddo

    close(100)
    close(200)

    pdata%t(1:pdata%n_train) = (/(i, i = 1, pdata%n_train)/)
    pdata%t_test(1:pdata%n_test) = (/(i, i = pdata%n_train + 1, pdata%n_train + pdata%n_test)/)
 
    Open(Unit = 100, File = trim(pwd)//"/../output/solution_covid.txt", ACCESS = "SEQUENTIAL")

    call lovo_algorithm(n,pdata%noutliers,pdata%outliers,pdata,pdata%fobj)

    write(100,10) pdata%xk(1),pdata%xk(2),pdata%xk(3)

    10 format (ES13.6,1X,ES13.6,1X,ES13.6) 
    close(100)
 
    ! allocate(pdata%t(pdata%n_train),pdata%y(pdata%n_train),pdata%indices(pdata%n_train),&
    !          pdata%sp_vector(pdata%n_train),pdata%outliers(pdata%noutliers),stat=allocerr)
 
    ! if ( allocerr .ne. 0 ) then
    !    write(*,*) 'Allocation error.'
    !    stop
    ! end if

    ! Open(Unit = 100, File = trim(pwd)//"/../output/solutions_covid.txt", ACCESS = "SEQUENTIAL")
 
    
    stop

    contains

    subroutine lovo_algorithm(n,noutliers,outliers,pdata,fobj)
        implicit none
  
        integer, intent(in) :: n,noutliers
        integer, intent(inout) :: outliers(noutliers)
        real(kind=8), intent(out) :: fobj
        type(pdata_type), intent(inout) :: pdata
  
        real(kind=8) :: sigmin,epsilon,fxk,fxtrial,alpha,gamma,termination
        integer :: iter_lovo,iter_sub_lovo,max_iter_lovo,max_iter_sub_lovo
  
        sigmin = 1.0d-1
        gamma = 1.0d+1
        epsilon = 1.0d0
        alpha = 1.0d-8
        max_iter_lovo = 1000000
        max_iter_sub_lovo = 100
        iter_lovo = 0
        iter_sub_lovo = 0
        pdata%lovo_order = pdata%n_train - noutliers
  
        pdata%xk(:) = 1.0d-3
        
        call compute_sp(n,pdata%xk,pdata,fxk)      
  
        write(*,*) "--------------------------------------------------------"
        write(*,10) "#iter","#init","Sp(xstar)","Stop criteria","#Imin"
        10 format (2X,A5,4X,A5,6X,A9,6X,A13,2X,A5)
        write(*,*) "--------------------------------------------------------"
  
        do
            iter_lovo = iter_lovo + 1
    
            call compute_grad_sp(n,pdata%xk,pdata,pdata%grad_sp)

            termination = norm2(pdata%grad_sp(:))
            ! termination = maxval(abs(pdata%grad_sp(:)))
            
            write(*,20)  iter_lovo,iter_sub_lovo,fxk,termination,pdata%dim_Imin
            20 format (I8,5X,I4,4X,ES14.6,3X,ES14.6,2X,I2)
    
            if (termination .lt. epsilon) exit
            if (iter_lovo .gt. max_iter_lovo) exit
            
            iter_sub_lovo = 1
            pdata%sigma = sigmin

            do 
                pdata%xtrial(:) = pdata%xk(:) - (1.d0 / pdata%sigma) * pdata%grad_sp(:)

                call compute_sp(n,pdata%xtrial,pdata,fxtrial)

                if (fxtrial .le. (fxk - alpha * norm2(pdata%xtrial(:) - pdata%xk(:))**2)) exit
                if (iter_sub_lovo .gt. max_iter_sub_lovo) exit

                pdata%sigma = max(sigmin,gamma * pdata%sigma)
                ! pdata%sigma = gamma * pdata%sigma
                iter_sub_lovo = iter_sub_lovo + 1

            enddo
  
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
  
        real(kind=8) :: gaux,t,tm
        integer :: i,ix
        
        res(:) = 0.0d0
        tm = pdata%t(pdata%n_train)
  
        do i = 1, pdata%lovo_order
            ix = int(pdata%indices(i))
            t = pdata%t(ix) - tm
        
            call model(n,x,ix,pdata,gaux)

            gaux = gaux - pdata%y(ix)
            
            res(1) = res(1) + gaux * t
            res(2) = res(2) + gaux * (t**2)
            res(3) = res(3) + gaux * (t**3)

            gaux = 0.d0
        enddo
  
    end subroutine compute_grad_sp

    !*****************************************************************
    !*****************************************************************

    subroutine model(n,x,i,pdata,res)
        implicit none 

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res
        real(kind=8) :: t,ym

        type(pdata_type), intent(in) :: pdata
   
        t = pdata%t(i) - pdata%t(pdata%n_train)
        ym = pdata%y(pdata%n_train)

        res = ym + (x(1) * t) + (x(2) * (t**2)) + (x(3) * (t**3))

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