program farrington
   use sort
   use bmgencan, only: gencan, genunc
   use iso_c_binding, only: c_ptr, c_loc, c_f_pointer
 
   implicit none
 
   type :: pdata_type
      integer :: counters(3) = 0
      integer :: samples,inf,sup,lovo_order
      real(kind=8) :: sigma,theta
      real(kind=8), allocatable :: xtrial(:),xk(:),t(:),y(:),data(:,:),indices(:),sp_vector(:),grad_sp(:),gp(:)
      integer, allocatable :: outliers(:)
   end type pdata_type
 
   ! LOCAL SCALARS
   logical :: extallowed,hfixstr
   integer :: allocerr,hnnzmax,ierr,istop,iter,maxit,n,nbds,status
   real(kind=8) :: bdsvio,eps,f,finish,ftarget,gpsupn,start 
   type(pdata_type), target :: pdata
 
   ! LOCAL ARRAYS
   character(len=10) :: pname
   logical, allocatable :: lind(:),uind(:)
   real(kind=8), allocatable :: g(:),lbnd(:),ubnd(:),x(:)
   integer :: i

   character(len=128) :: pwd
   call get_environment_variable('PWD',pwd)
 
   n = 3
 
   allocate(g(n),lind(n),lbnd(n),uind(n),ubnd(n),x(n),stat=allocerr)
   
   if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
   end if
 
   lbnd(1:n) = 0.0d0
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

   hnnzmax = 100000
 
   ftarget     =  - 1.0d+12
   eps         =    1.0d-08
   maxit       =      50000
   hfixstr     =     .true.
   extallowed  =     .true.
 
   nbds = count( lind(1:n) ) + count( uind(1:n) )

   ! Reading data and storing it in the variables t and y
   Open(Unit = 10, File = trim(pwd)//"/../data/seropositives.txt", Access = "SEQUENTIAL")

   ! Set parameters
   read(10,*) pdata%samples

   stop

   allocate(pdata%xtrial(n),pdata%xk(n),pdata%t(pdata%samples),pdata%y(pdata%samples),pdata%data(5,pdata%samples),&
   pdata%indices(pdata%samples),pdata%sp_vector(pdata%samples),pdata%grad_sp(n),pdata%gp(n),stat=allocerr)

   if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
   end if

   do i = 1, pdata%samples
      read(10,*) pdata%data(:,i)
   enddo

   pdata%t(:) = pdata%data(1,:)

   close(10)

   pdata%inf = 0
   pdata%sup = 10

   allocate(pdata%outliers(3*pdata%samples*(pdata%sup-pdata%inf+1)),stat=allocerr)

   if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error in main program'
       stop
   end if

   pdata%outliers(:) = 0

   call mixed_test(n)
 
   call export(pdata%sup,pdata%outliers,pdata%samples)

   stop

   contains

   ! *****************************************************************
   ! LOVO SUBROUTINES
   ! *****************************************************************

   subroutine mixed_test(n)
      implicit none

      integer, intent(in) :: n
      integer :: noutliers,ind,it
      real(kind=8) :: fobj


      do noutliers = pdata%inf, pdata%sup
         ! write(*,*) "LOVO Algorithm for Measles:"
         ind = 1
         pdata%y(:) = pdata%data(2,:)
         call cpu_time(start)
         call lovo_algorithm(n,noutliers,pdata%outliers(ind:ind+noutliers-1),fobj,it)
         call cpu_time(finish)
         Open(Unit = 100, File = trim(pwd)//"/../output/solutions_mixed_measles.txt", ACCESS = "SEQUENTIAL")
         Open(Unit = 110, File = trim(pwd)//"/../output/measles_latex.txt", ACCESS = "SEQUENTIAL")
         write(100,1000) pdata%xk(1), pdata%xk(2), pdata%xk(3)
         write(110,1010) fobj,it,pdata%counters(1),pdata%counters(2)

         pdata%counters(:) = 0
      
         ! write(*,*) "LOVO Algorithm for Mumps:"
         ind = ind + noutliers
         pdata%y(:) = pdata%data(3,:)
         call cpu_time(start)
         call lovo_algorithm(n,noutliers,pdata%outliers(ind:ind+noutliers-1),fobj,it)
         call cpu_time(finish)
         Open(Unit = 200, File = trim(pwd)//"/../output/solutions_mixed_mumps.txt", ACCESS = "SEQUENTIAL")
         Open(Unit = 210, File = trim(pwd)//"/../output/mumps_latex.txt", ACCESS = "SEQUENTIAL")
         write(200,1000) pdata%xk(1), pdata%xk(2), pdata%xk(3)
         write(210,1010) fobj,it,pdata%counters(1),pdata%counters(2)

         pdata%counters(:) = 0

         ! write(*,*) "LOVO Algorithm for Rubella:"
         ind = ind + noutliers
         pdata%y(:) = pdata%data(4,:)
         call cpu_time(start)
         call lovo_algorithm(n,noutliers,pdata%outliers(ind:ind+noutliers-1),fobj,it)
         call cpu_time(finish)
         Open(Unit = 300, File = trim(pwd)//"/../output/solutions_mixed_rubella.txt", ACCESS = "SEQUENTIAL")
         Open(Unit = 310, File = trim(pwd)//"/../output/rubella_latex.txt", ACCESS = "SEQUENTIAL")
         write(300,1000) pdata%xk(1), pdata%xk(2), pdata%xk(3)
         write(310,1010) fobj,it,pdata%counters(1),pdata%counters(2)
         
      enddo

      Open(Unit = 500, File = trim(pwd)//"/../output/num_mixed_test.txt", ACCESS = "SEQUENTIAL")
      write(500,1200) pdata%inf
      write(500,1200) pdata%sup

      1000 format (ES12.6,1X,ES12.6,1X,ES12.6)
      1010 format (ES10.3,1X,I4,1X,I6,1X,I6)
      1200 format (I2)
      close(100)
      close(200)
      close(300)
      CLOSE(500)

   end subroutine mixed_test

   !*****************************************************************
   !*****************************************************************

   subroutine lovo_algorithm(n,noutliers,outliers,fobj,it)
      implicit none

      integer, intent(in) :: n,noutliers
      integer, intent(inout) :: outliers(noutliers)
      integer, intent(out) :: it
      real(kind=8), intent(out) :: fobj

      real(kind=8) :: sp,sigmin,epsilon,fxk,fxtrial,alpha,gamma,termination
      integer :: iter_lovo,iter_sub_lovo,max_iter_lovo,max_iter_sub_lovo,i

      sigmin = 1.0d0
      epsilon = 1.0d-3
      alpha = 1.0d-8
      gamma = 1.0d+1
      max_iter_lovo = 10000
      max_iter_sub_lovo = 100
      iter_lovo = 0
      iter_sub_lovo = 0
      pdata%lovo_order = pdata%samples - noutliers

      pdata%theta = 100.d0

      pdata%xk(1:n) = 1.0d-1

      call compute_sp(n,pdata%xk,pdata,fxk)

      ! write(*,*) "--------------------------------------------------"
      ! write(*,10) "#iter","#init","Sp(xstar)","||gp(xstar)||"
      ! 10 format (2X,A5,4X,A5,6X,A9,6X,A13)
      ! write(*,*) "--------------------------------------------------"

      do
         iter_lovo = iter_lovo + 1

         call compute_grad_sp(n,pdata%xk,pdata,pdata%gp)

         do i = 1, n
            pdata%gp(i) = max(lbnd(i),min(pdata%xk(i) - pdata%gp(i),ubnd(i)))
         enddo

         termination = norm2(pdata%gp(1:n) - pdata%xk(1:n))

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

      fobj = fxtrial
      it = iter_lovo

      ! write(*,*) "--------------------------------------------------"

      outliers(:) = int(pdata%indices(pdata%samples - noutliers + 1:))
      

   end subroutine lovo_algorithm

   !*****************************************************************
   !*****************************************************************
   subroutine export(noutliers,outliers,samples)
      implicit none

      integer, intent(in) :: noutliers,outliers(3*samples),samples
      integer :: i

      Open(Unit = 200, File = trim(pwd)//"/../output/outliers.txt", ACCESS = "SEQUENTIAL")

      write(200,210) noutliers

      do i = 1, 3*noutliers
          write(200,210) outliers(i)
      enddo

      210 format (I2)

      close(200)

  end subroutine export

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

      real(kind=8) :: gaux1,gaux2,a,b,c,ebt,ti
      integer :: i,j
      
      a = x(1)
      b = x(2)
      c = x(3)
      
      res(:) = 0.0d0

      do i = 1, pdata%lovo_order
         ti = pdata%t(int(pdata%indices(i)))
         ebt = exp(-b * ti)

         call model(n,x,int(pdata%indices(i)),pdata,gaux1)

         gaux1 = pdata%y(int(pdata%indices(i))) - gaux1
         gaux2 = exp((a / b) * ti * ebt + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) - c * ti)

         res(1) = res(1) + gaux1 * gaux2 * ((1.0d0 / b**2) * (ebt * (ti * b + 1.0d0) - 1.0d0))

         res(2) = res(2) + gaux1 * gaux2 * (ebt * ((-2.0d0 * a * ti / b**2) - ((a * ti**2) / b) &
                  - (2.0d0 * a / b**3) + (c / b**2) + (c * ti / b)) + (2.0d0 * a / b**3) - (c / b**2))
    
         res(3) = res(3) + gaux1 * gaux2 * ((1.0d0 / b) * (1.0d0 - ebt) - ti)
      enddo

   end subroutine compute_grad_sp

   !*****************************************************************
   !*****************************************************************

   subroutine model(n,x,i,pdata,res)
      implicit none 

      integer,        intent(in) :: n,i
      real(kind=8),   intent(in) :: x(n)
      real(kind=8),   intent(out) :: res
      real(kind=8) :: a,b,c,ti,ebt

      type(pdata_type), intent(in) :: pdata

      a = x(1)
      b = x(2)
      c = x(3)
      ti = pdata%t(i)
      ebt = exp(-b * ti)

      res = (a / b) * ti * ebt
      res = res + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) 
      res = 1.0d0 - exp(res - c * ti)

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
     integer :: status
     type(pdata_type), pointer :: pdata
     
     call c_f_pointer(pdataptr, pdata)
     pdata%counters(1) = pdata%counters(1) + 1
     
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
     integer :: status
     type(pdata_type), pointer :: pdata
     
     call c_f_pointer(pdataptr, pdata)
     pdata%counters(2) = pdata%counters(2) + 1
     
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
     integer :: status
     type(pdata_type), pointer :: pdata
 
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
 
 end program farrington
 