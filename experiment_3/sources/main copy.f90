program main
    use sort

    implicit none

    type :: pdata_type
        integer :: counters(2) = 0
        integer :: samples,inf,sup,lovo_order,dim_Imin,n_train,n_test
        real(kind=8) :: sigma,theta
        real(kind=8), allocatable :: xstar(:),xtrial(:),xk(:),t(:),y(:),data(:,:),indices(:),sp_vector(:),&
        grad_sp(:),gp(:),lbnd(:),ubnd(:),pred(:),errors(:),hess_sp(:,:),eig_hess_sp(:),aux_mat(:,:),aux_vec(:)
        integer, allocatable :: outliers(:)
        character(len=1) :: JOBZ,UPLO ! lapack variables
        integer :: LDA,LWORK,INFO,NRHS,LDB ! lapack variables
        real(kind=8), allocatable :: WORK(:),IPIV(:) ! lapack variables
    end type pdata_type

    type(pdata_type), target :: pdata

    integer :: allocerr,i,n

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    n = 4

    ! lapack variables
    pdata%JOBZ = 'N'
    pdata%UPLO = 'U'
    pdata%LDA = n
    pdata%LDB = n
    pdata%LWORK = 3*n - 1
    pdata%NRHS = 1
    
    Open(Unit = 10, File = trim(pwd)//"/../data/cubic.txt", Access = "SEQUENTIAL")

    read(10,*) pdata%samples

    pdata%n_train = pdata%samples - 20
    pdata%n_test = pdata%samples - pdata%n_train

    allocate(pdata%xstar(n),pdata%xtrial(n),pdata%xk(n),pdata%t(pdata%n_train),pdata%y(pdata%n_train),&
    pdata%indices(pdata%n_train),pdata%sp_vector(pdata%n_train),pdata%grad_sp(n),pdata%gp(n),&
    pdata%data(2,pdata%samples),pdata%pred(pdata%n_test),pdata%errors(pdata%n_test),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
    end if

    allocate(pdata%hess_sp(n,n),pdata%eig_hess_sp(n),pdata%WORK(pdata%LWORK),pdata%aux_mat(n,n),&
    pdata%aux_vec(n),pdata%IPIV(n),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
    end if

    do i = 1, pdata%samples
        read(10,*) pdata%data(:,i)
    enddo

    close(10)

    pdata%xstar(:) = (/1.d0,1.d0,-3.d0,1.d0/)
  
    pdata%t(1:pdata%n_train) = pdata%data(1,1:pdata%n_train)
    pdata%y(1:pdata%n_train) = pdata%data(2,1:pdata%n_train)

    pdata%inf = 7
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
        integer :: noutliers,i,j
        real(kind=8) :: fobj,start,finish,ti,mean_abs_err,mae_train,p,o,aux
        type(pdata_type), intent(inout) :: pdata

        Open(Unit = 100, File = trim(pwd)//"/../output/solution_cubic.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 200, File = trim(pwd)//"/../output/output_latex.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 300, File = trim(pwd)//"/../output/log_sp.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 400, File = trim(pwd)//"/../output/errors.txt", ACCESS = "SEQUENTIAL")

        
        do noutliers = pdata%inf, pdata%sup
            write(*,*)
            write(*,*) "Outliers: ", noutliers

            call cpu_time(start)
            call lovo_algorithm(n,noutliers,pdata%outliers,pdata,fobj)
            call cpu_time(finish)

            ! print*, pdata%outliers(1:noutliers)

            do i = 1, pdata%n_test
                ti = pdata%data(1,i+pdata%n_train)
                pdata%pred(i) = pdata%xk(1) + (pdata%xk(2) * ti) + (pdata%xk(3) * (ti**2)) + (pdata%xk(4) * (ti**3))
                call absolute_error(pdata%data(2,i+pdata%n_train),pdata%pred(i),pdata%errors(i))
                ! call relative_error(pdata%data(2,i+pdata%n_train),pdata%pred(i),pdata%errors(i))
            enddo

            mean_abs_err = sum(pdata%errors)/pdata%n_test
            mae_train = 0.d0

            do i = 1, pdata%n_train - noutliers
                j = int(pdata%indices(i))
                ti = pdata%data(1,j)
                p = pdata%xk(1) + (pdata%xk(2) * ti) + (pdata%xk(3) * (ti**2)) + (pdata%xk(4) * (ti**3))
                o = pdata%data(2,j)
                call absolute_error(o,p,aux)
                mae_train = mae_train + aux
            enddo

            mae_train = mae_train / (pdata%n_train - noutliers)


            write(100,1000) pdata%xk(1),pdata%xk(2),pdata%xk(3),pdata%xk(4)
            write(200,1100) pdata%xk(1),pdata%xk(2),pdata%xk(3),pdata%xk(4),fobj,norm2(pdata%xk-pdata%xstar),&
            maxval(abs(pdata%xk-pdata%xstar)),mean_abs_err,mae_train,pdata%counters(1),pdata%counters(2)
            write(300,1300) fobj

            write(400,1400) pdata%errors
  
            pdata%counters(:) = 0
           
        enddo
  
        Open(Unit = 500, File = trim(pwd)//"/../output/num_mixed_test.txt", ACCESS = "SEQUENTIAL")
        write(500,1200) pdata%inf
        write(500,1200) pdata%sup
  
        1000 format (ES13.6,1X,ES13.6,1X,ES13.6,1X,ES13.6)
        1200 format (I2)
        1100 format (F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,I4,1X,I4)
        1300 format (ES13.6)
        1400 format (20ES13.6)

        close(100)
        close(200)
        close(300)
        close(400)
        close(500)

  
    end subroutine mixed_test

    subroutine lovo_algorithm(n,noutliers,outliers,pdata,fobj)
        implicit none
  
        integer, intent(in) :: n,noutliers
        integer, intent(inout) :: outliers(noutliers)
        real(kind=8), intent(out) :: fobj
        type(pdata_type), intent(inout) :: pdata
  
        real(kind=8) :: sigmin,epsilon,fxk,fxtrial,alpha,gamma,termination
        integer :: iter_lovo,iter_sub_lovo,max_iter_lovo,max_iter_sub_lovo
  
        sigmin = 1.0d-1
        gamma = 1.d+1
        epsilon = 1.0d-4
        alpha = 1.0d-8
        max_iter_lovo = 100
        max_iter_sub_lovo = 100
        iter_lovo = 0
        iter_sub_lovo = 0
        pdata%lovo_order = pdata%n_train - noutliers
  
        pdata%xk(1:n) =  -1.d0
        ! pdata%xk(:) = 1.0d-1
  
        call compute_sp(n,pdata%xk,pdata,fxk)      
  
        write(*,*) "--------------------------------------------------------"
        write(*,10) "#iter","#init","Sp(xstar)","Stop criteria","#Imin"
        10 format (2X,A5,4X,A5,6X,A9,6X,A13,2X,A5)
        write(*,*) "--------------------------------------------------------"
  
        do
            iter_lovo = iter_lovo + 1
    
            call compute_grad_sp(n,pdata%xk,pdata,pdata%grad_sp)
            call compute_Bkj(n,pdata)
    
            termination = norm2(pdata%grad_sp(1:n))
    
            write(*,20)  iter_lovo,iter_sub_lovo,fxk,termination,pdata%dim_Imin
            20 format (I6,5X,I4,4X,ES14.6,3X,ES14.6,2X,I2)
    
            if (termination .le. epsilon) exit
            if (iter_lovo .ge. max_iter_lovo) exit
            
            iter_sub_lovo = 1
            pdata%sigma = 0.d0

            do                 
                call compute_xtrial(n,pdata)

                ! pdata%xtrial(:) = pdata%xk(:) - (1.d0 / pdata%sigma) * pdata%grad_sp(:)
                call compute_sp(n,pdata%xtrial,pdata,fxtrial)

                if (fxtrial .le. (fxk - alpha * norm2(pdata%xtrial(1:n) - pdata%xk(1:n))**2)) exit
                if (iter_sub_lovo .ge. max_iter_sub_lovo) exit

                pdata%sigma = max(sigmin,gamma * pdata%sigma)
                iter_sub_lovo = iter_sub_lovo + 1

            enddo
  
            fxk = fxtrial
            pdata%xk(:) = pdata%xtrial(:)
            pdata%counters(2) = iter_sub_lovo + pdata%counters(2) + 1
  
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

    subroutine absolute_error(o,p,res)
        implicit none

        real(kind=8),   intent(in) :: o,p
        real(kind=8),   intent(out):: res

        res = abs(p - o)
    
    end subroutine absolute_error

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

    subroutine compute_hess_sp(n,pdata,res)
        implicit none

        integer,       intent(in) :: n
        real(kind=8),  intent(out) :: res(n,n)
        type(pdata_type), intent(in) :: pdata

        real(kind=8) :: ti
        integer :: i

        res(:,:) = 0.0d0

        do i = 1, pdata%lovo_order
            ti = pdata%t(int(pdata%indices(i)))

            res(1,:) = res(1,:) + (/ti**0,ti**1,ti**2,ti**3/)
            res(2,:) = res(2,:) + (/ti**1,ti**2,ti**3,ti**4/) 
            res(3,:) = res(3,:) + (/ti**2,ti**3,ti**4,ti**5/)
            res(4,:) = res(4,:) + (/ti**3,ti**4,ti**5,ti**6/)
            
        enddo
    
    end subroutine compute_hess_sp

    !*****************************************************************
    !*****************************************************************

    subroutine compute_Bkj(n,pdata)
        implicit none

        integer,            intent(in) :: n
        type(pdata_type),   intent(inout) :: pdata
        real(kind=8) :: lambda_min

        call compute_hess_sp(n,pdata,pdata%hess_sp)

        pdata%aux_mat(:,:) = pdata%hess_sp(:,:)

        call dsyev(pdata%JOBZ,pdata%UPLO,n,pdata%aux_mat,pdata%LDA,&
        pdata%eig_hess_sp,pdata%WORK,pdata%LWORK,pdata%INFO)

        lambda_min = minval(pdata%eig_hess_sp)
        call compute_eye(n,pdata%aux_mat)

        pdata%hess_sp(:,:) = pdata%hess_sp(:,:) + &
        max(0.d0,-lambda_min + 1.d-8) * pdata%aux_mat(:,:)
               
    end subroutine compute_Bkj

    !*****************************************************************
    !*****************************************************************

    subroutine compute_xtrial(n,pdata)
        implicit none 

        integer,            intent(in) :: n
        type(pdata_type),   intent(inout) :: pdata

        call compute_eye(n,pdata%aux_mat)

        pdata%aux_mat(:,:) = pdata%hess_sp(:,:) + pdata%sigma * pdata%aux_mat(:,:)

        pdata%aux_vec(:) = matmul(pdata%aux_mat(:,:),pdata%xk(:))
        pdata%aux_vec(:) = pdata%aux_vec(:) - pdata%grad_sp(:)

        call dsysv(pdata%UPLO,n,pdata%NRHS,pdata%aux_mat(:,:),pdata%LDA,pdata%IPIV,&
        pdata%aux_vec(:),pdata%LDB,pdata%WORK,pdata%LWORK,pdata%INFO)

        pdata%xtrial(:) = pdata%aux_vec(:)
    end subroutine compute_xtrial

    !*****************************************************************
    !*****************************************************************

    subroutine compute_eye(n,res)
        implicit none

        integer,        intent(in) :: n
        real(kind=8),   intent(out):: res(n,n)
        integer :: i

        res(:,:) = 0.0d0

        do i = 1, n
            res(i,i) = 1.d0
        enddo
    end subroutine compute_eye

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