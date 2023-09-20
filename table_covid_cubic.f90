program table
    implicit none

    real(kind=8), allocatable :: y_true(:),y_pred
    real(kind=8) :: error
    integer :: i

    contains

    subroutine mean_squared_error(y_true,y_pred,n,res)
        implicit none
  
        integer,       intent(in) :: n
        real(kind=8),  intent(in) :: y_true(n),y_pred(n)
        real(kind=8),  intent(out) :: res
  
        res = (1.0d0 / n) * norm2(y_pred(1:n) - y_true(1:n))**2
  
     end subroutine mean_squared_error
end program table