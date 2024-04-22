program main
    use sort

    implicit none

    type :: pdata_type
        integer :: counters(2) = 0
        integer :: lovo_order,n_train,n_test,days_test,noutliers,dim_Imin
        real(kind=8) :: sigma,fobj
        real(kind=8), allocatable :: xtrial(:),xk(:),t(:),y(:),t_test(:),y_test(:),indices(:),&
        pred(:),re(:),sp_vector(:),grad_sp(:),hess_sp(:,:),eig_hess_sp(:),aux_mat(:,:),aux_vec(:),&
        test_data(:,:),train_data(:,:)
        integer, allocatable :: outliers(:)
        character(len=1) :: JOBZ,UPLO ! lapack variables
        integer :: LDA,LWORK,INFO,NRHS,LDB ! lapack variables
        real(kind=8), allocatable :: WORK(:),IPIV(:) ! lapack variables
    end type pdata_type

    type(pdata_type), target :: pdata

    integer :: allocerr,n

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    ! Number of variables
    n = 3
    
    ! lapack variables
    pdata%JOBZ = 'N'
    pdata%UPLO = 'U'
    pdata%LDA = n
    pdata%LDB = n
    pdata%LWORK = 3*n - 1
    pdata%NRHS = 1
 
    call single_test()
 
    stop

    contains

    !*****************************************************************
    !*****************************************************************

    subroutine single_test()
        implicit none

        integer :: i

        Open(Unit = 100, File = trim(pwd)//"/../data/covid.txt", Access = "SEQUENTIAL")
    
        read(100,*) pdata%n_train
        read(100,*) pdata%n_test
    
        ! pdata%noutliers = 0*int(dble(pdata%n_train) / 7.0d0)
        pdata%noutliers = 10
    
        allocate(pdata%t(pdata%n_train),pdata%y(pdata%n_train),pdata%y_test(pdata%n_test),pdata%t_test(pdata%n_test),&
        pdata%xtrial(n),pdata%xk(n),pdata%grad_sp(n),pdata%indices(pdata%n_train),stat=allocerr)
     
        if ( allocerr .ne. 0 ) then
           write(*,*) 'Allocation error.'
           stop
        end if
    
        allocate(pdata%sp_vector(pdata%n_train),pdata%outliers(pdata%noutliers),pdata%hess_sp(n,n),pdata%eig_hess_sp(n),&
        pdata%WORK(pdata%LWORK),pdata%aux_mat(n,n),pdata%aux_vec(n),pdata%IPIV(n),stat=allocerr)
     
        if ( allocerr .ne. 0 ) then
           write(*,*) 'Allocation error.'
           stop
        end if
    
        do i = 1, pdata%n_train
            read(100,*) pdata%y(i)
        enddo
    
        do i = 1, pdata%n_test
            read(100,*) pdata%y_test(i)
        enddo

        close(100)
    
        pdata%t(1:pdata%n_train) = (/(i, i = 1, pdata%n_train)/)
        pdata%t_test(1:pdata%n_test) = (/(i, i = pdata%n_train + 1, pdata%n_train + pdata%n_test)/)
     
        Open(Unit = 100, File = trim(pwd)//"/../output/solution_covid.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 200, File = trim(pwd)//"/../output/outliers.txt", ACCESS = "SEQUENTIAL")
    
        call lovo_algorithm(n,pdata%noutliers,pdata%outliers,pdata,.true.,pdata%fobj)

        print*, pdata%counters
    
        write(100,10) pdata%xk(1),pdata%xk(2),pdata%xk(3)
        write(200,20) pdata%noutliers

        do i = 1, pdata%noutliers
            write(200,20) pdata%outliers(i)
        enddo
    
        10 format (ES13.6,1X,ES13.6,1X,ES13.6) 
        20 format (I2)
        30 format (ES10.3)
    
        close(100)
        close(200)

        deallocate(pdata%t,pdata%y,pdata%y_test,pdata%t_test,pdata%xtrial,pdata%xk,pdata%grad_sp,&
        pdata%indices,pdata%sp_vector,pdata%outliers,pdata%hess_sp,pdata%eig_hess_sp,pdata%WORK,&
        pdata%aux_mat,pdata%aux_vec,pdata%IPIV,stat=allocerr)

        if ( allocerr .ne. 0 ) then
            write(*,*) 'Deallocation error.'
            stop
        end if

    end subroutine single_test

    !*****************************************************************
    !*****************************************************************

    subroutine lovo_algorithm(n,noutliers,outliers,pdata,single_type_test,fobj)
        implicit none
        
        logical,        intent(in) :: single_type_test
        integer,        intent(in) :: n,noutliers
        integer,        intent(inout) :: outliers(noutliers)
        real(kind=8),   intent(out) :: fobj
        type(pdata_type), intent(inout) :: pdata
  
        real(kind=8) :: sigmin,epsilon,fxk,fxtrial,alpha,gamma,termination
        integer :: iter_lovo,iter_sub_lovo,max_iter_lovo,max_iter_sub_lovo
  
        sigmin = 1.0d-1
        gamma = 1.0d+1
        epsilon = 1.0d-4
        alpha = 1.0d-8
        max_iter_lovo = 1000
        max_iter_sub_lovo = 100
        iter_lovo = 0
        iter_sub_lovo = 0
        pdata%lovo_order = pdata%n_train - noutliers
  
        pdata%xk(:) = 1.0d-5
        
        call compute_sp(n,pdata%xk,pdata,fxk)  

        if (single_type_test) then
            write(*,*) "--------------------------------------------------------"
            write(*,10) "#iter","#init","Sp(xstar)","Stop criteria","#Imin"
            10 format (2X,A5,4X,A5,6X,A9,6X,A13,2X,A5)
            write(*,*) "--------------------------------------------------------"
        endif
  
        do
            iter_lovo = iter_lovo + 1
    
            call compute_grad_sp(n,pdata%xk,pdata,pdata%grad_sp)
            call compute_Bkj(n,pdata)

            termination = norm2(pdata%grad_sp(:))
            
            if (single_type_test) then
                write(*,20)  iter_lovo,iter_sub_lovo,fxk,termination,pdata%dim_Imin
                20 format (I8,5X,I4,4X,ES14.6,3X,ES14.6,2X,I2)
            endif
    
            if (termination .le. epsilon) exit
            if (iter_lovo .gt. max_iter_lovo) exit
            
            iter_sub_lovo = 1
            pdata%sigma = 0.d0

            do                 
                call compute_xtrial(n,pdata)
                call compute_sp(n,pdata%xtrial,pdata,fxtrial)

                if (fxtrial .le. (fxk - alpha * norm2(pdata%xtrial(:) - pdata%xk(:))**2)) exit
                if (iter_sub_lovo .gt. max_iter_sub_lovo) exit

                pdata%sigma = max(sigmin,gamma * pdata%sigma)
                iter_sub_lovo = iter_sub_lovo + 1

            enddo
  
            fxk = fxtrial
            pdata%xk(:) = pdata%xtrial(:)
            pdata%counters(2) = iter_sub_lovo + pdata%counters(2)
  
        enddo
  
        fobj = fxtrial
        pdata%counters(1) = iter_lovo
  
        if (single_type_test) then
            write(*,*) "--------------------------------------------------------"
        endif

        outliers(:) = int(pdata%indices(pdata%n_train - noutliers + 1:))

    end subroutine lovo_algorithm

    !*****************************************************************
    !*****************************************************************

    subroutine cross_validation(ncv,pdata,covid_data)
        implicit none

        integer, intent(in) :: ncv
        type(pdata_type), intent(inout) :: pdata
        real(kind=8), intent(in) :: covid_data(:)

        integer :: i,init,end
        
        do i = 1, ncv 
            init = 1 + (30 - pdata%n_train) + (i-1) * 40
            end = (30 - pdata%n_train) + (i-1) * 40 + pdata%n_train
            pdata%train_data(i,:) = covid_data(init:end)
            pdata%test_data(i,:) = covid_data(end+1:end + pdata%n_test) 
        enddo
        
    end subroutine cross_validation

    !*****************************************************************
    !*****************************************************************

    subroutine mount_dataset(pdata,covid_data)
        implicit none
  
        type(pdata_type), intent(inout) :: pdata
        real(kind=8), intent(in) :: covid_data(:)
  
        integer :: i
  
        do i = 1, pdata%days_test
           pdata%train_data(i,:) = covid_data(i:i+pdata%n_train-1)
           pdata%test_data(i,:) = covid_data(i+pdata%n_train:i+pdata%n_train+pdata%n_test-1)    
        enddo
  
     end subroutine mount_dataset

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
  
        real(kind=8) :: gaux,ti
        integer :: i
        
        res(:) = 0.0d0
  
        do i = 1, pdata%lovo_order
            ti = pdata%t(int(pdata%indices(i)))
            call model(n,x,int(pdata%indices(i)),pdata,gaux)
            gaux = gaux - pdata%y(int(pdata%indices(i)))
            
            res(1) = res(1) + gaux * (ti - pdata%t(pdata%n_train))
            res(2) = res(2) + gaux * ((ti - pdata%t(pdata%n_train))**2)
            res(3) = res(3) + gaux * ((ti - pdata%t(pdata%n_train))**3)
         enddo
  
    end subroutine compute_grad_sp

    !*****************************************************************
    !*****************************************************************

    subroutine compute_hess_sp(n,pdata,res)
        implicit none

        integer,       intent(in) :: n
        real(kind=8),  intent(out) :: res(n,n)
        type(pdata_type), intent(in) :: pdata

        real(kind=8) :: ti,tm
        integer :: i

        res(:,:) = 0.0d0
        tm = pdata%t(pdata%n_train)

        do i = 1, pdata%lovo_order
            ti = pdata%t(int(pdata%indices(i)))

            res(1,:) = res(1,:) + (/(ti-tm)**2,(ti-tm)**3,(ti-tm)**4/) 
            res(2,:) = res(2,:) + (/(ti-tm)**3,(ti-tm)**4,(ti-tm)**5/) 
            res(3,:) = res(3,:) + (/(ti-tm)**4,(ti-tm)**5,(ti-tm)**6/) 
            
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

        type(pdata_type), intent(in) :: pdata
   
        res = pdata%y(pdata%n_train) + x(1) * (pdata%t(i) - pdata%t(pdata%n_train)) + &
        x(2) * ((pdata%t(i) - pdata%t(pdata%n_train))**2) + x(3) * ((pdata%t(i) - pdata%t(pdata%n_train))**3)

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