program main
    use sort

    implicit none

    type :: pdata_type
        integer :: counters(2) = 0
        integer :: samples,inf,sup,lovo_order,dim_Imin
        real(kind=8) :: sigma,theta
        real(kind=8), allocatable :: xtrial(:),xk(:),t(:),y(:),data(:,:),indices(:),sp_vector(:),grad_sp(:),&
        gp(:),lbnd(:),ubnd(:),sol_ovo_measles(:,:),sol_ovo_mumps(:,:),sol_ovo_rubella(:,:)
        integer, allocatable :: outliers(:)
    end type pdata_type

    type(pdata_type), target :: pdata

    integer :: allocerr,i,n

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    n = 3
    
    Open(Unit = 10, File = trim(pwd)//"/../data/seropositives.txt", Access = "SEQUENTIAL")

    read(10,*) pdata%samples

    allocate(pdata%xtrial(n),pdata%xk(n),pdata%t(pdata%samples),pdata%y(pdata%samples),pdata%data(5,pdata%samples),&
    pdata%indices(pdata%samples),pdata%sp_vector(pdata%samples),pdata%grad_sp(n),pdata%gp(n),pdata%lbnd(n),&
    pdata%ubnd(n),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
    end if

    do i = 1, pdata%samples
        read(10,*) pdata%data(:,i)
    enddo

    close(10)
  
    pdata%t(:) = pdata%data(1,:)

    pdata%inf = 1
    pdata%sup = 10

    pdata%lbnd(1:n) = 0.0d0
    pdata%ubnd(1:n) = 1.0d+20
 
    allocate(pdata%outliers(3*pdata%samples*(pdata%sup-pdata%inf+1)),pdata%sol_ovo_measles(pdata%sup-pdata%inf+1,3),&
    pdata%sol_ovo_mumps(pdata%sup-pdata%inf+1,3),pdata%sol_ovo_rubella(pdata%sup-pdata%inf+1,3),stat=allocerr)
 
    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if
 
    pdata%outliers(:) = 0

    Open(Unit = 10, File = trim(pwd)//"/../data/sol_ovo_measles.txt", Access = "SEQUENTIAL")
    Open(Unit = 20, File = trim(pwd)//"/../data/sol_ovo_mumps.txt", Access = "SEQUENTIAL")
    Open(Unit = 30, File = trim(pwd)//"/../data/sol_ovo_rubella.txt", Access = "SEQUENTIAL")

    do i = 1, pdata%sup-pdata%inf+1
        read(10,*) pdata%sol_ovo_measles(i,:)
        read(20,*) pdata%sol_ovo_mumps(i,:)
        read(30,*) pdata%sol_ovo_rubella(i,:)
    enddo

    call mixed_test(n,pdata)

    Open(Unit = 10, File =trim(pwd)//"/../output/outliers.txt", ACCESS = "SEQUENTIAL")

    write(10,100) pdata%sup

    do i = 1, 3*pdata%sup
        write(10,100) pdata%outliers(i)
    enddo

    100 format (I2)

    contains

    subroutine mixed_test(n,pdata)
        implicit none
  
        integer, intent(in) :: n
        integer :: noutliers,ind
        real(kind=8) :: fobj,start,finish,lovo_fun_xlovo,ovo_fun_xlovo
        type(pdata_type), intent(inout) :: pdata

        Open(Unit = 100, File = trim(pwd)//"/../output/solutions_mixed_measles.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 200, File = trim(pwd)//"/../output/solutions_mixed_mumps.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 300, File = trim(pwd)//"/../output/solutions_mixed_rubella.txt", ACCESS = "SEQUENTIAL")

        Open(Unit = 110, File = trim(pwd)//"/../output/measles_latex.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 111, File = trim(pwd)//"/../output/measles_latex2.txt", ACCESS = "SEQUENTIAL")

        Open(Unit = 210, File = trim(pwd)//"/../output/mumps_latex.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 211, File = trim(pwd)//"/../output/mumps_latex2.txt", ACCESS = "SEQUENTIAL")

        Open(Unit = 310, File = trim(pwd)//"/../output/rubella_latex.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 311, File = trim(pwd)//"/../output/rubella_latex2.txt", ACCESS = "SEQUENTIAL")
  
        do noutliers = pdata%inf, pdata%sup
            ! write(*,*) "LOVO Algorithm for Measles:"
            ind = 1
            pdata%y(:) = pdata%data(2,:)
            pdata%xk(:) =  (/0.197d0,0.287d0,0.021d0/)

            call cpu_time(start)
            call lovo_algorithm(n,noutliers,pdata%outliers(ind:ind+noutliers-1),pdata,fobj)
            call cpu_time(finish)

            lovo_fun_xlovo = fobj
            ovo_fun_xlovo  = pdata%sp_vector(pdata%lovo_order)
            call compute_sp(n,pdata%sol_ovo_measles(noutliers,:),pdata,fobj) 

            write(100,1000) pdata%xk(1), pdata%xk(2), pdata%xk(3)
            write(110,1010) fobj,pdata%counters(1),pdata%counters(2),finish-start
            write(111,1011) fobj,lovo_fun_xlovo,ovo_fun_xlovo,pdata%sp_vector(pdata%lovo_order)
  
            pdata%counters(:) = 0
        
            ! write(*,*) "LOVO Algorithm for Mumps:"
            ind = ind + noutliers
            pdata%y(:) = pdata%data(3,:)
            pdata%xk(:) = (/0.156d0,0.250d0,0.0d0/)

            call cpu_time(start)
            call lovo_algorithm(n,noutliers,pdata%outliers(ind:ind+noutliers-1),pdata,fobj)
            call cpu_time(finish)

            lovo_fun_xlovo = fobj
            ovo_fun_xlovo  = pdata%sp_vector(pdata%lovo_order)
            call compute_sp(n,pdata%sol_ovo_mumps(noutliers,:),pdata,fobj)

            write(200,1000) pdata%xk(1), pdata%xk(2), pdata%xk(3)
            write(210,1010) fobj,pdata%counters(1),pdata%counters(2),finish-start
            write(211,1011) fobj,lovo_fun_xlovo,ovo_fun_xlovo,pdata%sp_vector(pdata%lovo_order)
  
            pdata%counters(:) = 0
  
            ! write(*,*) "LOVO Algorithm for Rubella:"
            ind = ind + noutliers
            pdata%y(:) = pdata%data(4,:)
            pdata%xk(:) = (/0.0628d0,0.178d0,0.020d0/)
            call cpu_time(start)
            call lovo_algorithm(n,noutliers,pdata%outliers(ind:ind+noutliers-1),pdata,fobj)
            call cpu_time(finish)

            lovo_fun_xlovo = fobj
            ovo_fun_xlovo  = pdata%sp_vector(pdata%lovo_order)
            call compute_sp(n,pdata%sol_ovo_mumps(noutliers,:),pdata,fobj)
            
            write(300,1000) pdata%xk(1), pdata%xk(2), pdata%xk(3)
            write(310,1010) fobj,pdata%counters(1),pdata%counters(2),finish-start
            write(311,1011) fobj,lovo_fun_xlovo,ovo_fun_xlovo,pdata%sp_vector(pdata%lovo_order)
           
        enddo
  
        Open(Unit = 500, File = trim(pwd)//"/../output/num_mixed_test.txt", ACCESS = "SEQUENTIAL")
        write(500,1200) pdata%inf
        write(500,1200) pdata%sup
  
        1000 format (ES12.6,1X,ES12.6,1X,ES12.6)
        1010 format (ES10.3,1X,I4,1X,I4,1X,ES10.3)
        1011 format (ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3)
        1200 format (I2)

        close(100)
        close(200)
        close(300)
        close(500)
        close(110)
        close(111)
        close(210)
        close(211)
        close(310)
        close(311)

  
    end subroutine mixed_test

    subroutine lovo_algorithm(n,noutliers,outliers,pdata,fobj)
        implicit none
  
        integer, intent(in) :: n,noutliers
        integer, intent(inout) :: outliers(noutliers)
        real(kind=8), intent(out) :: fobj
        type(pdata_type), intent(inout) :: pdata
  
        real(kind=8) :: sigmin,epsilon,fxk,fxtrial,alpha,gamma,termination
        integer :: iter_lovo,iter_sub_lovo,max_iter_lovo,max_iter_sub_lovo,i
  
        sigmin = 1.0d-1
        epsilon = 1.0d-3
        alpha = 1.0d-8
        gamma = 5.d0
        max_iter_lovo = 10000
        max_iter_sub_lovo = 100
        iter_lovo = 0
        iter_sub_lovo = 0
        pdata%lovo_order = pdata%samples - noutliers
  
        pdata%theta = 1.d0
  
        ! pdata%xk(1:n) = 1.0d-2
  
        call compute_sp(n,pdata%xk,pdata,fxk)      
  
        ! write(*,*) "--------------------------------------------------------"
        ! write(*,10) "#iter","#init","Sp(xstar)","Stop criteria","#Imin"
        ! 10 format (2X,A5,4X,A5,6X,A9,6X,A13,2X,A5)
        ! write(*,*) "--------------------------------------------------------"
  
        do
            iter_lovo = iter_lovo + 1
    
            call compute_grad_sp(n,pdata%xk,pdata,pdata%grad_sp)
    
            do i = 1, n
                pdata%gp(i) = max(pdata%lbnd(i),min(pdata%xk(i) - pdata%grad_sp(i),pdata%ubnd(i)))
            enddo
    
            termination = norm2(pdata%gp(1:n) - pdata%xk(1:n))
    
            ! write(*,20)  iter_lovo,iter_sub_lovo,fxk,termination,pdata%dim_Imin
            ! 20 format (I6,5X,I4,4X,ES14.6,3X,ES14.6,2X,I2)
    
            if (termination .lt. epsilon) exit
            if (iter_lovo .gt. max_iter_lovo) exit
            
            iter_sub_lovo = 1
            pdata%sigma = sigmin

            do 
                do i = 1, n
                    pdata%xtrial(i) = max(pdata%lbnd(i),min(pdata%xk(i) - (1.d0 / pdata%sigma) * pdata%grad_sp(i),pdata%ubnd(i)))
                enddo

                call compute_sp(n,pdata%xtrial,pdata,fxtrial)

                if (fxtrial .le. (fxk - alpha * norm2(pdata%xtrial(1:n-1) - pdata%xk(1:n-1))**2)) exit
                if (iter_sub_lovo .gt. max_iter_sub_lovo) exit

                pdata%sigma = gamma * pdata%sigma
                iter_sub_lovo = iter_sub_lovo + 1

            enddo
  
            fxk = fxtrial
            pdata%xk(:) = pdata%xtrial(:)
            pdata%counters(2) = iter_sub_lovo + pdata%counters(2)
  
        enddo
  
        fobj = fxtrial
        pdata%counters(1) = iter_lovo
  
        ! write(*,*) "--------------------------------------------------------"

  
        outliers(:) = int(pdata%indices(pdata%samples - noutliers + 1:))
        
  
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
        pdata%indices(:) = (/(i, i = 1, pdata%samples)/)

        do i = 1, pdata%samples
            call fi(n,x,i,pdata,pdata%sp_vector(i))
        end do

        ! Sorting
        call DSORT(pdata%sp_vector,pdata%indices,pdata%samples,kflag)

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

        real(kind=8) :: gaux1,gaux2,a,b,c,ebt,ti
        integer :: i

        a = x(1)
        b = x(2)
        c = x(3)

        res(:) = 0.0d0

        do i = 1, pdata%lovo_order
            ti = pdata%t(int(pdata%indices(i)))
            ebt = exp(-b * ti)

            call model(n,x,int(pdata%indices(i)),pdata,gaux1)

            gaux1 = pdata%y(int(pdata%indices(i))) - gaux1
            gaux2 = exp((a / b) * ti * ebt + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) - c * ti)

            res(1) = res(1) + gaux1 * gaux2 * ((1.0d0 / b**2) * (ebt * (ti * b + 1.0d0) - 1.0d0))

            res(2) = res(2) + gaux1 * gaux2 * (ebt * ((-2.0d0 * a * ti / b**2) - ((a * ti**2) / b) &
                    - (2.0d0 * a / b**3) + (c / b**2) + (c * ti / b)) + (2.0d0 * a / b**3) - (c / b**2))

            res(3) = res(3) + gaux1 * gaux2 * ((1.0d0 / b) * (1.0d0 - ebt) - ti)
        enddo

    end subroutine compute_grad_sp

    !*****************************************************************
    !*****************************************************************

    subroutine model(n,x,i,pdata,res)
        implicit none 

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: res
        real(kind=8) :: a,b,c,ti,ebt

        type(pdata_type), intent(in) :: pdata

        a = x(1)
        b = x(2)
        c = x(3)
        ti = pdata%t(i)
        ebt = exp(-b * ti)

        res = (a / b) * ti * ebt
        res = res + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) 
        res = 1.0d0 - exp(res - c * ti)

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