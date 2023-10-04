program prueba
    implicit none

    CHARACTER(len=255) :: cwd

    call getwd(cwd)

    write(*,*) trim(cwd)

    ! Open(Unit = 100, File = "/output/covid_train.txt", ACCESS = "SEQUENTIAL")

end program prueba