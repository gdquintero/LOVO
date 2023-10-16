program table
    implicit none

    real(kind=8) ::  y_true,y_pred,error
    real(kind=8), allocatable :: train_set(:),test_set(:),x(:),accuracy(:,:)
    integer :: i,j,allocerr,n_train,n_test,inf

    ! Reading data and storing it in the variables t and y
    Open(Unit = 100, File = "data/covid_train.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 200, File = "data/covid_test.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 300, File = "output/solutions_covid_cubic.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 400, File = "output/solutions_covid_logistic.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 500, File = "output/table_covid_cubic.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 600, File = "output/inf_sup_covid.txt", ACCESS = "SEQUENTIAL")
    Open(Unit = 700, File = "output/table_covid_logistic.txt", ACCESS = "SEQUENTIAL")

    ! Set parameters
    read(100,*) n_train
    read(200,*) n_test
    read(600,*) inf

    allocate(train_set(n_train),test_set(n_test),x(3),accuracy(n_train - inf + 1,n_test),stat=allocerr)

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

    close(100)
    close(200)
    close(600)

    do i = 1, n_train - inf + 1
        read(300,*) x
        do j = 1, n_test
            y_true = test_set(j)
            call cubic_model(x,inf+i+j-1,inf+i-1,train_set(n_train),3,y_pred)
            call percentage_error(y_true,y_pred,accuracy(i,j))
        enddo
    enddo

    do i = 1, n_train - inf + 1
        write(500,10) inf+i-1,"&",accuracy(i,1),"&",accuracy(i,2),"&",accuracy(i,3),"&",accuracy(i,4),"&",accuracy(i,5),"&",&
        accuracy(i,6),"&",accuracy(i,7),"&",accuracy(i,8),"&",accuracy(i,9),"&",accuracy(i,10),"\\"
    enddo

    close(300)
    close(500)

    do i = 1, n_train - inf + 1
        read(400,*) x
        do j = 1, n_test
            y_true = test_set(j)
            call logistic_model(x,inf+i+j-1,3,y_pred)
            call percentage_error(y_true,y_pred,accuracy(i,j))
        enddo
    enddo

    do i = 1, n_train - inf + 1
        write(700,10) inf+i-1,"&",accuracy(i,1),"&",accuracy(i,2),"&",accuracy(i,3),"&",accuracy(i,4),"&",accuracy(i,5),"&",&
        accuracy(i,6),"&",accuracy(i,7),"&",accuracy(i,8),"&",accuracy(i,9),"&",accuracy(i,10),"\\"
    enddo

    10 format(I2,2X,A1,2X,F8.2,2X,A1,2X,F8.2,2X,A1,2X,F8.2,2X,A1,2X,F8.2,2X,A1,2X,F8.2,2X,A1,2X,F8.2,2X,A1,2X,F8.2,2X,A1,2X,&
            F8.2,2X,A1,2X,F8.2,2X,A1,2X,F8.2,2X,A2)

    close(400)
    close(700)

    contains

    subroutine cubic_model(x,t,tm,ym,n,res)
        implicit none

        integer,        intent(in) :: n,t,tm
        real(kind=8),   intent(in) :: x(n),ym
        real(kind=8),   intent(out) :: res

        res = ym + x(1) * (t - tm) + x(2) * (t - tm)**2 + x(3) * (t - tm)**3

    end subroutine cubic_model

    subroutine logistic_model(x,t,n,res)
        implicit none

        integer,        intent(in) :: n,t
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out):: res

        res = x(1) * x(3) * exp(x(2) * t) / (x(3) + x(1) * (exp(x(2) * t) - 1.0))
    end subroutine logistic_model

    subroutine percentage_error(y_true,y_pred,res)
        implicit none
  
        real(kind=8),  intent(in) :: y_true,y_pred
        real(kind=8),  intent(out) :: res
  
        res = (y_pred - y_true) / y_true
        res = abs(res) * 100
  
     end subroutine percentage_error
end program table