program eigen
    implicit none

    character(len=1) :: JOBZ,UPLO
    integer :: N,LDA,LDB,LWORK,INFO,i,NRHS
    real(kind=8) :: dot
    real(kind=8), allocatable :: A(:,:),W(:),WORK(:),IPIV(:),b(:),c(:)

    JOBZ = 'N'
    UPLO = 'U'
    N = 2
    LDA = N
    LDB = N
    LWORK = 3*N - 1
    NRHS = 1

    allocate(A(n,n),W(n),WORK(LWORK),IPIV(n),b(n),c(n))

    A(1,:) = (/1.d0,2.d0/)
    A(2,:) = (/2.d0,3.d0/)

    b(:) = (/2.d0,2.d0/)

    ! call dsyev(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)

    call dsysv (UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO)

    print*, B
end program eigen