program prueba
    implicit none

    integer :: i

    do i = 1, 15
        print*, mod(dble(i),7.0d0)
    enddo

end program prueba