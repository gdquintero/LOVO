program main
    use sort

    implicit none

    type :: pdata_type
        integer :: counters(2) = 0,n
        real(kind=8) :: sigma,fobj
        real(kind=8), allocatable :: xtrial(:),xk(:),sp_vector(:),grad_sp(:),hess_sp(:,:),eig_hess_sp(:),aux_mat(:,:),aux_vec(:)
        integer, allocatable :: outliers(:)
        character(len=1) :: JOBZ,UPLO ! lapack variables
        integer :: LDA,LWORK,INFO,NRHS,LDB ! lapack variables
        real(kind=8), allocatable :: WORK(:),IPIV(:) ! lapack variables
    end type pdata_type

    type(pdata_type), target :: pdata

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    ! Number of variables
    pdata%n = 3
    
    ! lapack variables
    pdata%JOBZ = 'N'
    pdata%UPLO = 'U'
    pdata%LDA = pdata%n
    pdata%LDB = pdata%n
    pdata%LWORK = 3*pdata%n - 1
    pdata%NRHS = 1

    call hard_test()
 
    stop

    contains

    !*****************************************************************
    !*****************************************************************
    subroutine hard_test()
        implicit none

        integer :: samples,i,j,k,n_train,n_test,optimal_ntrain,noutliers,allocerr,out_per_ndays,total_test
        real(kind=8) :: fobj,ti,ym,pred,av_err_train,av_err_test,start,finish
        real(kind=8), allocatable :: t(:),t_test(:),covid_data(:),indices(:),sp_vector(:),abs_err(:,:),av_abs_err(:)
        integer, allocatable :: outliers(:)

        Open(Unit = 100, File = trim(pwd)//"/../data/covid_mixed.txt", Access = "SEQUENTIAL")
        Open(Unit = 200, File = trim(pwd)//"/../output/solutions_covid_mixed.txt", Access = "SEQUENTIAL")
        Open(Unit = 300, File = trim(pwd)//"/../output/output_covid_hard.txt", Access = "SEQUENTIAL")
    
        samples = 100
        n_test = 5

        allocate(covid_data(samples),t(25),t_test(n_test),indices(25),sp_vector(25),outliers(10),&
        pdata%hess_sp(pdata%n,pdata%n),pdata%eig_hess_sp(pdata%n),pdata%WORK(pdata%LWORK),&
        pdata%aux_mat(pdata%n,pdata%n),pdata%aux_vec(pdata%n),pdata%IPIV(pdata%n),pdata%xtrial(pdata%n),&
        pdata%xk(pdata%n),pdata%grad_sp(pdata%n),abs_err(21,n_test),av_abs_err(21),stat=allocerr)

        if ( allocerr .ne. 0 ) then
            write(*,*) 'Allocation error.'
            stop
        end if

        t(:) = (/(i, i = 1, 25)/)

        do i = 1, samples
            read(100,*) covid_data(i)
        enddo
    
        close(100)

        out_per_ndays = 1
        total_test = 71
        total_test = 1

        call cpu_time(start)

        do i = 1, total_test
            ym = covid_data(24+i)

            j = 0
        
            do n_train = 5,25,5
                j = j + 1
                noutliers = out_per_ndays * n_train / 5

                t_test = (/(k+1, k = n_train,n_train+4)/)
                indices(:) = (/(k, k = 1, 25)/)

                call lovo_algorithm(t(1:n_train),covid_data(25+i-n_train:24+i),indices(1:n_train),&
                outliers,n_train,noutliers,sp_vector(1:n_train),pdata,.false.,fobj)

                do k = 1, n_test
                    ti = t_test(k) / t(n_train)
                    pred = ym + pdata%xk(1) * (ti - 1.d0) + &
                            pdata%xk(2) * ((ti - 1.d0)**2) + pdata%xk(3) * ((ti - 1.d0)**3)

                    abs_err(j,k) = absolute_error(covid_data(25+i),pred)
                enddo   
                av_abs_err(j) = sum(abs_err(j,:)) / n_test
            enddo    

            optimal_ntrain = 4 + minloc(av_abs_err(1:j),optimal_ntrain)

            av_err_test = av_abs_err(optimal_ntrain-4) 

            indices(:) = (/(k, k = 1, 25)/)
            noutliers = out_per_ndays * optimal_ntrain / 5
            outliers(:) = 0

            call lovo_algorithm(t(1:optimal_ntrain),covid_data(25+i-optimal_ntrain:24+i),&
            indices(1:optimal_ntrain),outliers,optimal_ntrain,noutliers,sp_vector(1:optimal_ntrain),pdata,.false.,fobj)

            do k = 1, optimal_ntrain
                ti = t(k) / t(optimal_ntrain)
                pred = ym + pdata%xk(1) * (ti - 1.d0) + &
                        pdata%xk(2) * ((ti - 1.d0)**2) + pdata%xk(3) * ((ti - 1.d0)**3)

                if (.not. ANY(outliers(1:noutliers) .eq. k)) then
                    av_err_train = av_err_train + absolute_error(covid_data(k),pred)
                endif
            enddo  

            av_err_train = av_err_train / (optimal_ntrain - noutliers)

            write(200,10) pdata%xk(1),pdata%xk(2),pdata%xk(3)
            write(300,20) i,"&",optimal_ntrain,"&",fobj,"&",av_err_train,"&",av_err_test,"\\"

            av_err_train = 0.d0
        enddo

        call cpu_time(finish)

        10 format (ES13.6,1X,ES13.6,1X,ES13.6)
        20 format (I3,1X,A1,1X,I2,1X,A1,1X,ES10.3,1X,A1,1X,ES10.3,1X,A1,1X,ES10.3,1X,A2)

        close(200)
        close(300)

        print*
        print*, "Test Finished!!"
        write(*,30) " Time: ",finish - start

        30 format (A7,1X,F6.4)

    end subroutine hard_test

    !*****************************************************************
    !*****************************************************************

    subroutine lovo_algorithm(t,y,indices,outliers,n_train,noutliers,sp_vector,pdata,single_type_test,fobj)
        implicit none
        
        logical,        intent(in) :: single_type_test
        integer,        intent(in) :: n_train,noutliers
        real(kind=8),   intent(in) :: y(n_train),t(n_train)
        integer,        intent(inout) :: outliers(noutliers)
        real(kind=8),   intent(inout) :: indices(n_train),sp_vector(n_train)
        real(kind=8),   intent(out) :: fobj
        type(pdata_type), intent(inout) :: pdata
  
        real(kind=8) :: sigmin,epsilon,fxk,fxtrial,alpha,gamma,termination
        integer :: iter_lovo,iter_sub_lovo,max_iter_lovo,max_iter_sub_lovo,lovo_order
  
        sigmin = 1.0d-1
        gamma = 1.0d+1
        epsilon = 1.0d-3
        alpha = 1.0d-8
        max_iter_lovo = 1000
        max_iter_sub_lovo = 100
        iter_lovo = 0
        iter_sub_lovo = 0
        lovo_order = n_train - noutliers
  
        pdata%xk(:) = 1.0d-1
        
        call compute_sp(pdata%xk,t,y,indices,sp_vector,pdata%n,n_train,lovo_order,fxk)

        if (single_type_test) then
            write(*,*) "--------------------------------------------------------"
            write(*,10) "#iter","#init","Sp(xstar)","Stop criteria","#Imin"
            10 format (2X,A5,4X,A5,6X,A9,6X,A13,2X,A5)
            write(*,*) "--------------------------------------------------------"
        endif
  
        do
            iter_lovo = iter_lovo + 1
    
            call compute_grad_sp(pdata%xk,t,y,indices,pdata%n,n_train,lovo_order,pdata%grad_sp)
            call compute_Bkj(t,indices,n_train,lovo_order,pdata)

            termination = norm2(pdata%grad_sp(:))
            
            if (single_type_test) then
                write(*,20)  iter_lovo,iter_sub_lovo,fxk,termination
                20 format (I8,5X,I4,4X,ES14.6,3X,ES14.6,2X,I2)
            endif
    
            if (termination .le. epsilon) exit
            if (iter_lovo .gt. max_iter_lovo) exit
            
            iter_sub_lovo = 1
            pdata%sigma = 0.d0

            do                 
                call compute_xtrial(pdata)
                call compute_sp(pdata%xtrial,t,y,indices,sp_vector,pdata%n,n_train,lovo_order,fxtrial)

                if (fxtrial .le. (fxk - alpha * norm2(pdata%xtrial(:) - pdata%xk(:))**2)) exit
                if (iter_sub_lovo .gt. max_iter_sub_lovo) exit

                pdata%sigma = max(sigmin,gamma * pdata%sigma)
                iter_sub_lovo = iter_sub_lovo + 1

            enddo
  
            fxk = fxtrial
            pdata%xk(:) = pdata%xtrial(:)
            pdata%counters(2) = iter_sub_lovo + pdata%counters(2) + 1
  
        enddo
  
        fobj = fxtrial
        pdata%counters(1) = iter_lovo
  
        if (single_type_test) then
            write(*,*) "--------------------------------------------------------"
        endif
        
        outliers(:) = int(indices(n_train - noutliers + 1:))

    end subroutine lovo_algorithm

    !*****************************************************************
    !*****************************************************************

    real(kind=8) function absolute_error(o,p)
        implicit none

        real(kind=8) :: o,p

        absolute_error = abs(p - o)
    
    end function absolute_error

    !*****************************************************************
    !*****************************************************************

    subroutine compute_sp(x,t,y,indices,sp_vector,n,n_train,lovo_order,res)
        implicit none
        integer,        intent(in) :: n,n_train,lovo_order
        real(kind=8),   intent(in) :: x(n),y(n_train),t(n_train)
        real(kind=8),   intent(inout) :: indices(n_train),sp_vector(n_train)
        real(kind=8),   intent(out) :: res

        integer :: i,kflag

        sp_vector(:) = 0.0d0
        kflag = 2
        indices(1:n_train) = (/(i, i = 1,n_train)/)

        do i = 1, n_train
            call fi(x,i,t,y,n,n_train,sp_vector(i))
        end do

        ! Sorting
        call DSORT(sp_vector,indices,n_train,kflag)

        ! Lovo function
        res = sum(sp_vector(1:lovo_order))

    end subroutine compute_sp

    !*****************************************************************
    !*****************************************************************

    subroutine compute_grad_sp(x,t,y,indices,n,n_train,lovo_order,res)
        implicit none
  
        integer,        intent(in) :: n,n_train,lovo_order
        real(kind=8),   intent(in) :: x(n),y(n_train),t(n_train)
        real(kind=8),   intent(inout) :: indices(n_train)
        real(kind=8),   intent(out) :: res(n)
  
        real(kind=8) :: gaux,ti
        integer :: i
        
        res(:) = 0.0d0
  
        do i = 1, lovo_order
            ti = (t(int(indices(i))) / t(n_train)) - 1.d0
            call model(x,int(indices(i)),t,y,n,n_train,gaux)
            gaux = gaux - y(int(indices(i)))
            
            res(1) = res(1) + gaux * (ti)
            res(2) = res(2) + gaux * ((ti)**2)
            res(3) = res(3) + gaux * ((ti)**3)
         enddo
  
    end subroutine compute_grad_sp

    !*****************************************************************
    !*****************************************************************

    subroutine compute_hess_sp(t,indices,n,n_train,lovo_order,res)
        implicit none

        integer,        intent(in) :: n,n_train,lovo_order
        real(kind=8),   intent(in) :: t(n_train),indices(n_train)
        real(kind=8),   intent(out) :: res(n,n)

        real(kind=8) :: ti,tm
        integer :: i

        res(:,:) = 0.0d0
        tm = t(n_train)

        do i = 1, lovo_order
            ti = (t(int(indices(i))) / tm) - 1.d0

            res(1,:) = res(1,:) + (/(ti)**2,(ti)**3,(ti)**4/) 
            res(2,:) = res(2,:) + (/(ti)**3,(ti)**4,(ti)**5/) 
            res(3,:) = res(3,:) + (/(ti)**4,(ti)**5,(ti)**6/) 
            
        enddo
    
    end subroutine compute_hess_sp

    !*****************************************************************
    !*****************************************************************

    subroutine compute_Bkj(t,indices,n_train,lovo_order,pdata)
        implicit none

        integer,            intent(in) :: n_train,lovo_order
        real(kind=8),       intent(in) :: t(n_train)
        real(kind=8),       intent(inout) :: indices(n_train)
        type(pdata_type),   intent(inout) :: pdata
        real(kind=8) :: lambda_min

        call compute_hess_sp(t,indices,pdata%n,n_train,lovo_order,pdata%hess_sp)

        pdata%aux_mat(:,:) = pdata%hess_sp(:,:)

        call dsyev(pdata%JOBZ,pdata%UPLO,pdata%n,pdata%aux_mat,pdata%LDA,&
        pdata%eig_hess_sp,pdata%WORK,pdata%LWORK,pdata%INFO)

        lambda_min = minval(pdata%eig_hess_sp)
        call compute_eye(pdata%n,pdata%aux_mat)

        pdata%hess_sp(:,:) = pdata%hess_sp(:,:) + &
        max(0.d0,-lambda_min + 1.d-8) * pdata%aux_mat(:,:)
               
    end subroutine compute_Bkj

    !*****************************************************************
    !*****************************************************************

    subroutine compute_xtrial(pdata)
        implicit none 

        type(pdata_type),   intent(inout) :: pdata

        call compute_eye(pdata%n,pdata%aux_mat)

        pdata%aux_mat(:,:) = pdata%hess_sp(:,:) + pdata%sigma * pdata%aux_mat(:,:)

        pdata%aux_vec(:) = matmul(pdata%aux_mat(:,:),pdata%xk(:))
        pdata%aux_vec(:) = pdata%aux_vec(:) - pdata%grad_sp(:)

        call dsysv(pdata%UPLO,pdata%n,pdata%NRHS,pdata%aux_mat(:,:),pdata%LDA,pdata%IPIV,&
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

    subroutine fi(x,i,t,y,n,n_train,res)
        implicit none

        integer,        intent(in) :: n,n_train,i
        real(kind=8),   intent(in) :: x(n),y(n_train),t(n_train)
        real(kind=8),   intent(out) :: res

        call model(x,i,t,y,n,n_train,res)
        res = res - y(i)
        res = 0.5d0 * (res**2)

    end subroutine fi

    !*****************************************************************
    !*****************************************************************

    subroutine model(x,i,t,y,n,n_train,res)
        implicit none 

        integer,        intent(in) :: n,n_train,i
        real(kind=8),   intent(in) :: x(n),y(n_train),t(n_train)
        real(kind=8),   intent(out) :: res
        real(kind=8) :: ti

        ti = (t(i) / t(n_train)) - 1.d0

        res = y(n_train) + x(1) * (ti) + &
        x(2) * ((ti)**2) + x(3) * ((ti)**3)

    end subroutine model
    
end program main