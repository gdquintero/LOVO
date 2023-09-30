program pruebas
    implicit none 

    real, pointer :: p
    real, target :: t1 = 10., t2 = -17.

    p => t1

    WRITE (*,*) 'p, t1, t2 = ', p, t1, t2

    p => t2

    WRITE (*,*) 'p, t1, t2 = ', p, t1, t2

end program pruebas