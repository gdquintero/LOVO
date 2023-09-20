program table
    implicit none

    real(kind=8) ::  y_true,y_pred,error
    real(kind=8), allocatable :: train_set(:),test_set(:),x(:),accuracy(:,:)
    integer :: i,allocerr,n_train,n_test

    ! Reading data and storing it in the variables t and y
    Open(Unit = 100, File = "output/covid_train.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 200, File = "output/covid_test.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 300, File = "output/solutions_covid_cubic.txt", ACCESS = "SEQUENTIAL")

    ! Set parameters
    read(100,*) n_train
    read(200,*) n_test

    allocate(train_set(n_train),test_set(n_test),x(3),accuracy(n_train-4,n_test),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
    end if

    do i = 1, n_train
        read(100,*) train_set(i)
    enddo

    do i = 1, n_test
        read(200,*) test_set(i)
    enddo

    read(300,*) x

    close(100)
    close(200)

    y_true = test_set(1)
    call cubic_model(x,6,5,train_set(n_train),3,y_pred)

    call percentage_error(y_true,y_pred,error)

    print*, y_true,y_pred,error


    contains

    subroutine cubic_model(x,t,tm,ym,n,res)
        implicit none

        integer,        intent(in) :: n,t,tm
        real(kind=8),   intent(in) :: x(n),ym
        real(kind=8),   intent(out) :: res

        res = ym + x(1) * (t - tm) + x(2) * (t - tm)**2 + x(3) * (t - tm)**3

    end subroutine cubic_model

    subroutine percentage_error(y_true,y_pred,res)
        implicit none
  
        real(kind=8),  intent(in) :: y_true,y_pred
        real(kind=8),  intent(out) :: res
  
        res = (y_pred - y_true) / y_true
        res = abs(res) * 100
  
     end subroutine percentage_error
end program table