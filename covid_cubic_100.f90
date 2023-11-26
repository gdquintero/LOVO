program covid
   use sort
   use bmgencan, only: gencan, genunc
   use iso_c_binding, only: c_ptr, c_loc,c_f_pointer 

   implicit none

   type :: pdata_type
      integer :: counters(3) = 0
      integer :: samples,lovo_order,n_train,n_test,noutliers,days_test
      real(kind=8) :: sigma,theta
      real(kind=8), allocatable :: xtrial(:),xk(:),t(:),y(:),indices(:),sp_vector(:),grad_sp(:),&
                                   gp(:),train_set(:),test_set(:),data_test(:,:),data_train(:,:)
      integer, allocatable :: outliers(:)
   end type pdata_type

   ! LOCAL SCALARS
   logical :: extallowed,hfixstr
   integer :: allocerr,hnnzmax,ierr,istop,iter,maxit,n,nbds
   real(kind=8) :: eps,f,finish,ftarget,gpsupn,start
   real(kind=8), allocatable :: covid_data(:)
   type(pdata_type), target :: pdata
 
   ! LOCAL ARRAYS
   ! character(len=10) :: pname
   logical, allocatable :: lind(:),uind(:)
   real(kind=8), allocatable :: g(:),lbnd(:),ubnd(:),x(:)

   integer :: i,j,k,sam,n_data
   real(kind=8), allocatable :: fobj(:)

   ! Number of variables

   n = 3

   ! Set parameters 
   pdata%n_train  = 30
   pdata%n_test   = 5
   pdata%days_test = 100

   Open(Unit = 100, File = "data/covid.txt", ACCESS = "SEQUENTIAL")
   read(100,*) n_data

   allocate(g(n),lind(n),lbnd(n),uind(n),ubnd(n),x(n),pdata%data_train(100,pdata%n_train),&
            pdata%data_test(100,pdata%n_test),covid_data(n_data),stat=allocerr)

   if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
   end if

   do i = 1, n_data
      read(100,*) covid_data(i)
   enddo

   close(100)

   call mount_dataset(pdata,covid_data)
 
   lbnd(1:n) = - 1.0d+20
   ubnd(1:n) = 1.0d+20
 
   where( lbnd(1:n) .gt. - 1.0d+20 )
      lind(1:n) = .true.
   elsewhere
      lind(1:n) = .false.
   end where
 
   where( ubnd(1:n) .lt.   1.0d+20 )
      uind(1:n) = .true.
   elsewhere
      uind(1:n) = .false.
   end where

   hnnzmax     = 100000
   ftarget     =  - 1.0d+12
   eps         =    1.0d-08
   maxit       =      50000
   hfixstr     =     .true.
   extallowed  =     .true.

   allocate(pdata%train_set(pdata%n_train),pdata%test_set(pdata%n_test),&
   pdata%xtrial(n),pdata%xk(n),pdata%grad_sp(n),pdata%gp(n),stat=allocerr)

   if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
   end if

   pdata%samples = pdata%n_train
   pdata%noutliers = 0*int(dble(pdata%samples) / 7.0d0)

   allocate(pdata%t(pdata%samples),pdata%y(pdata%samples),pdata%indices(pdata%samples),&
            pdata%sp_vector(pdata%samples),pdata%outliers(pdata%noutliers),stat=allocerr)

   if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
   end if

   do k = 1, pdata%days_test

      pdata%train_set(:)   = pdata%data_train(k,:)
      pdata%test_set(:)    = pdata%data_test(k,:)
      pdata%t(:)           = (/(i, i = 1, pdata%samples)/)
      pdata%indices(:)     = (/(i, i = 1, pdata%samples)/)
      pdata%y(:)           = pdata%train_set(:)

      call lovo_algorithm()
      Open(Unit = 100, File = "output/solutions_covid_cubic_100.txt", ACCESS = "SEQUENTIAL")
      write(100,10) pdata%xk(1), pdata%xk(2), pdata%xk(3)

      print*, k * 100/pdata%days_test,"%"

   enddo

   10 format (ES13.6,1X,ES13.6,1X,ES13.6) 
   close(100)

   call cpu_time(start)

   call cpu_time(finish)

   deallocate(pdata%t,pdata%y,pdata%indices,pdata%sp_vector,pdata%outliers,stat=allocerr)
      
   if ( allocerr .ne. 0 ) then
      write(*,*) 'Deallocation error.'
      stop
   end if
   
   deallocate(lind,lbnd,uind,ubnd,x,stat=allocerr)
   
   if ( allocerr .ne. 0 ) then
      write(*,*) 'Deallocation error.'
      stop
   end if

   call export(pdata)
   
   stop
  
   contains

   ! *****************************************************************
   ! LOVO SUBROUTINES
   ! *****************************************************************

   subroutine lovo_algorithm()
      implicit none

      real(kind=8) :: sigmin,epsilon,fxk,alpha,gamma,termination,fxtrial
      integer :: iter_lovo,iter_sub_lovo,max_iter_lovo,max_iter_sub_lovo

      sigmin = 1.0d0
      epsilon = 1.0d-3
      alpha = 1.0d-8
      gamma = 1.0d+1
      max_iter_lovo = 1000
      max_iter_sub_lovo = 100
      iter_lovo = 0
      iter_sub_lovo = 0
      
      pdata%lovo_order = pdata%samples - pdata%noutliers
      
      pdata%theta = 100.d0

      pdata%xk(1:n) = 1.0d-2
      
      call compute_sp(n,pdata%xk,pdata,fxk)

      ! write(*,*)
      ! write(*,99) "Main algorithm with", sam, "previous days"
      ! 99 format (1X,A19,1X,I2,1X,A13)
      ! write(*,*) "--------------------------------------------------"
      ! write(*,10) "#iter","#init","Sp(xstar)","||g(xstar)||"
      ! 10 format (2X,A5,4X,A5,6X,A9,7X,A12)
      ! write(*,*) "--------------------------------------------------"

      do
         iter_lovo = iter_lovo + 1

         call compute_grad_sp(n,x,pdata,pdata%gp)

         ! termination = norm2(pdata%gp(1:n))
         termination = maxval(abs(pdata%gp(1:n)))

         ! write(*,20)  iter_lovo,iter_sub_lovo,fxk,termination
         ! 20 format (I6,5X,I4,4X,ES14.6,3X,ES14.6)

         if (termination .lt. epsilon) exit
         if (iter_lovo .gt. max_iter_lovo) exit

         x(1:n) = pdata%xk(1:n)
         pdata%sigma = sigmin
        
         iter_sub_lovo = 1

         do

            if ( nbds .eq. 0 ) then
               ! write(*,*) 'The problem is unconstrained.'
               call genunc(evalf,evalg,evalh,hnnzmax,hfixstr,n,x,f,g,gpsupn, &
                  ftarget,eps,maxit,extallowed,iter,ierr,istop,stpsub,pdata=c_loc(pdata))
               
            else
               ! write(*,*) 'The problem is a bound constrained problem.'
               call gencan(evalf,evalg,evalh,hnnzmax,hfixstr,n,x,lind,lbnd, &
                  uind,ubnd,f,g,gpsupn,ftarget,eps,maxit,extallowed,iter,ierr, &
                  istop,stpsub,pdata=c_loc(pdata))
            end if

            pdata%xtrial(:) = x(:)

            call compute_sp(n,pdata%xtrial,pdata,fxtrial)

            if (fxtrial .le. (fxk - alpha * norm2(pdata%xtrial(1:n-1) - pdata%xk(1:n-1))**2)) exit
            if (iter_sub_lovo .gt. max_iter_sub_lovo) exit

            pdata%sigma = gamma * pdata%sigma
            iter_sub_lovo = iter_sub_lovo + 1

         enddo

         fxk = fxtrial
         pdata%xk(:) = pdata%xtrial(:)

      enddo

      ! write(*,*) "--------------------------------------------------"

      pdata%outliers(:) = int(pdata%indices(pdata%samples - pdata%noutliers + 1:))

   end subroutine lovo_algorithm

   !*****************************************************************
   !*****************************************************************
   subroutine mount_dataset(pdata,covid_data)
      implicit none

      type(pdata_type), intent(inout) :: pdata
      real(kind=8), intent(in) :: covid_data(:)

      integer :: i

      do i = 1, 100
         pdata%data_train(i,:) = covid_data(i:i+pdata%n_train-1)
         pdata%data_test(i,:) = covid_data(i+pdata%n_train:i+pdata%n_train+pdata%n_test-1)    
      enddo

   end subroutine mount_dataset

   !*****************************************************************
   !*****************************************************************
   subroutine export(pdata)
      implicit none

      type(pdata_type), intent(in) :: pdata

      real(kind=8) ::  y_true,y_pred
      real(kind=8), allocatable :: xsol(:),accuracy(:,:)
      integer :: i,j 


      Open(Unit = 100, File = "output/solutions_covid_cubic_100.txt", ACCESS = "SEQUENTIAL")
      Open(Unit = 200, File = "output/mean_percentage_error.txt", ACCESS = "SEQUENTIAL")

      allocate(xsol(3),accuracy(pdata%days_test,pdata%n_test),stat=allocerr)

      if ( allocerr .ne. 0 ) then
          write(*,*) 'Allocation error.'
          stop
      end if

      do i = 1, pdata%days_test
         read(100,*) xsol

         do j = 1, pdata%n_test
             y_true = pdata%data_test(i,j)
             call cubic_model(xsol,pdata%n_train+j,pdata%n_train,pdata%data_train(i,pdata%n_train),3,y_pred)
             call percentage_error(y_true,y_pred,accuracy(i,j))
         enddo
         write(200,10) sum(accuracy(i,:))/pdata%n_test
         ! write(200,10) pdata%inf+i-1,accuracy(i,:),sum(accuracy(i,:))/pdata%n_test
     enddo

     10 format (F8.2)

     close(100)
     close(200)
     deallocate(xsol)
   end subroutine export

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
      pdata%indices(:) = (/(i, i = 1, pdata%samples)/)

      do i = 1, pdata%samples
         call fi(n,x,i,pdata,pdata%sp_vector(i))
      end do
      
      ! Sorting
      call DSORT(pdata%sp_vector,pdata%indices,pdata%samples,kflag)

      ! Lovo function
      res = sum(pdata%sp_vector(1:pdata%lovo_order))

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
         
         res(1) = res(1) + gaux * (ti - pdata%t(pdata%samples))
         res(2) = res(2) + gaux * ((ti - pdata%t(pdata%samples))**2)
         res(3) = res(3) + gaux * ((ti - pdata%t(pdata%samples))**3)
      enddo

   end subroutine compute_grad_sp

   !*****************************************************************
   !*****************************************************************

   subroutine model(n,x,i,pdata,res)
      implicit none 

      integer,        intent(in) :: n,i
      real(kind=8),   intent(in) :: x(n)
      real(kind=8),   intent(out) :: res

      type(pdata_type), intent(in) :: pdata

      res = pdata%y(pdata%samples) + x(1) * (pdata%t(i) - pdata%t(pdata%samples)) + &
      x(2) * ((pdata%t(i) - pdata%t(pdata%samples))**2) + x(3) * ((pdata%t(i) - pdata%t(pdata%samples))**3)

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

   ! *****************************************************************
   ! GENCAN SUBROUTINES
   ! *****************************************************************

   subroutine stpsub(n,x,gsupn,inhdefstp,stp,ierr,pdataptr)
      use iso_c_binding, only: c_ptr
      implicit none
      integer, intent(in) :: n
      real(kind=8), intent(in) :: gsupn
      real(kind=8), intent(in) :: x(n)
      logical, intent(out) :: inhdefstp,stp
      integer, intent(inout) :: ierr
      type(c_ptr), optional, intent(in) :: pdataptr

      type(pdata_type), pointer :: pdata
      call c_f_pointer(pdataptr, pdata)

      if ( .false. ) write(*,*) ierr

      inhdefstp = .false.
      
      if ( gsupn .le. pdata%theta * maxval(abs(x - pdata%xk))) then
         stp = .true.
      else
         stp = .false.
      endif

   end subroutine stpsub

   ! *****************************************************************
   ! *****************************************************************
 
   subroutine evalf(n,x,f,inform,pdataptr)
     implicit none
 
     ! SCALAR ARGUMENTS
     integer, intent(in) :: n
     integer, intent(inout) :: inform
     real(kind=8), intent(out) :: f
     type(c_ptr), optional, intent(in) :: pdataptr
 
     ! ARRAY ARGUMENTS
     real(kind=8), intent(in) :: x(n)
     
     ! LOCAL SCALARS
   !   integer :: status
     type(pdata_type), pointer :: pdata

     if ( .false. ) write(*,*) inform
     
     call c_f_pointer(pdataptr, pdata)
   !   pdata%counters(1) = pdata%counters(1) + 1
     
   !   call cutest_ufn(status,n,x,f)
   !   if ( status .ne. 0 ) inform = -91

     call compute_sp(n,x,pdata,f)
     call compute_grad_sp(n,x,pdata,pdata%grad_sp)
     f = f + dot_product(pdata%grad_sp(1:n),x(1:n) - pdata%xk(1:n)) + &
         0.5d0 * pdata%sigma * (norm2(x(1:n) - pdata%xk(1:n))**2)
         
   
   end subroutine evalf
 
   ! *****************************************************************
   ! *****************************************************************
 
   subroutine evalg(n,x,g,inform,pdataptr)
 
     implicit none
 
     ! SCALAR ARGUMENTS
     integer, intent(in) :: n
     integer, intent(inout) :: inform
     type(c_ptr), optional, intent(in) :: pdataptr
   
     ! ARRAY ARGUMENTS
     real(kind=8), intent(in) :: x(n)
     real(kind=8), intent(out) :: g(n)
     
     ! LOCAL SCALARS
   !   integer :: status
     type(pdata_type), pointer :: pdata

     if ( .false. ) write(*,*) inform
     
     call c_f_pointer(pdataptr, pdata)
   !   pdata%counters(2) = pdata%counters(2) + 1
     
   !   call cutest_ugr(status,n,x,g)
   !   if ( status .ne. 0 ) inform = -92
 
   !   where ( isnan( g(1:n) ) )
   !      g(1:n) = 0.0d0
   !   end where

     call compute_grad_sp(n,x,pdata,g)
     g = g + pdata%sigma * (x(1:n) - pdata%xk(1:n))
     
   end subroutine evalg
 
   ! *****************************************************************
   ! *****************************************************************
 
   subroutine evalh(n,x,lim,hnnz,hrow,hcol,hval,inform,pdataptr)
 
     implicit none
 
     ! SCALAR ARGUMENTS
     integer, intent(in) :: lim,n
     integer, intent(inout) :: inform
     integer, intent(out) :: hnnz
     type(c_ptr), optional, intent(in) :: pdataptr
 
     ! ARRAY ARGUMENTS
     real(kind=8), intent(in) :: x(n)
     integer, intent(out) :: hcol(lim),hrow(lim)
     real(kind=8), intent(out) :: hval(lim)
   
 
     ! LOCAL SCALARS
   !   integer :: status
     type(pdata_type), pointer :: pdata

     if ( .false. ) write(*,*) inform, x(1:n)
 
     call c_f_pointer(pdataptr, pdata)
   !   pdata%counters(3) = pdata%counters(3) + 1
 
   !   ! Upper triangle
   !   call cutest_ush(status,n,x,hnnz,lim,hval,hrow,hcol)
   !   if ( status .ne. 0 ) inform = -93
     
   !   where ( isnan( hval(1:hnnz) ) )
   !      hval(1:hnnz) = 0.0d0
   !   end where

     hnnz = n

     hrow(1:n) = (/(i, i = 1, n)/)
     hcol(1:n) = (/(i, i = 1, n)/)
     hval(1:n) = pdata%sigma
     
   end subroutine evalh
  
end program covid