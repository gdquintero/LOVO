program eigen
    implicit none

    character(len=1) :: JOBZ,UPLO
    integer :: N,LDA,LWORK,INFO,i
    real(kind=8), allocatable :: A(:,:),W(:),WORK(:)

    JOBZ = 'N'
    UPLO = 'U'
    N = 3
    LDA = N
    LWORK = 3*N - 1

    allocate(A(n,n),W(n),WORK(LWORK))

    A(1,:) = (/204.d0,-1296.d0,8772.d0/)
    A(2,:) = (/-1296.d0,8772.d0,-61776.d0/)
    A(3,:) = (/8772.d0,-61776.d0,446964.d0/)

    call dsyev(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)

    print*,W
    
end program eigen