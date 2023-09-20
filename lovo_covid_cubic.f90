! *****************************************************************
! *****************************************************************

! This is a main program that calls Algencan to solve a simple
! problem. It is intended to be used as an (incomplete) example of
! usage. Algencan applies to problems of the form
!
! Minimize f(x)
!
! subject to
!
!   heq(x) = 0    heq : R^n \to R^m represents the equality constraints
!   hin(x) <= 0   hin : R^n \to R^p represents the inequality constraints
!   l <= x <= u   l and u \in R^n are the bound constraints.

! *****************************************************************
! *****************************************************************

program algencama
   use sort
   use bmalgencan, only: algencan
   use iso_c_binding, only: c_ptr, c_loc,c_f_pointer

   implicit none

   ! Re-define this type (pdata_type) anyway you want. Algencan
   ! receives a pointer to a 'structure of this type'. Algencan has no
   ! access to the structure. It simple passes the pointer back to the
   ! user defined subroutines evalf, evalg, evalc, evalj, and
   ! evalhl. So, this is a trade safe way (no common blocks) of passing
   ! to the user-provided routines any information related to the
   ! problem. In this example, it is only being used for the user to
   ! count by itself the number of calls to each routine.
   type :: pdata_type
      integer :: counters(5) = 0
   end type pdata_type

   ! LOCAL SCALARS
   logical :: corrin,extallowed,rhoauto,scale
   integer :: allocerr,hlnnzmax,ierr,inform,istop,jnnzmax,m,maxoutit,n,nwcalls,nwtotit,outiter,p,totiter
   real(kind=8) :: bdsvio,csupn,epsfeas,epscompl,epsopt,f,finish,nlpsupn,rhoini,ssupn,start
   type(pdata_type), target :: pdata
   
   ! LOCAL ARRAYS
   logical, allocatable :: lind(:),uind(:)
   real(kind=8), allocatable :: c(:),lbnd(:),ubnd(:),lambda(:),x(:)

   !--> LOVO Algorithm variables <--
   integer :: samples,inf,sup,lovo_order,n_train,n_test,noutliers
   real(kind=8) :: sigma
   real(kind=8), allocatable :: xtrial(:),xk(:),t(:),y(:),data(:,:),indices(:),sp_vector(:),grad_sp(:),&
                                gp(:),train_set(:),test_set(:)
   integer, allocatable :: outliers(:)
   integer :: i

   ! Number of variables

   n = 3

   ! Reading data and storing it in the variables t and y
   Open(Unit = 100, File = "output/covid_train3.txt", ACCESS = "SEQUENTIAL")
   Open(Unit = 200, File = "output/covid_test3.txt", ACCESS = "SEQUENTIAL")

   ! Set parameters
   read(100,*) n_train
   read(200,*) n_test

   allocate(x(n),lind(n),lbnd(n),uind(n),ubnd(n),xtrial(n),xk(n), &
      train_set(n_train),test_set(n_test),grad_sp(n),gp(n),stat=allocerr)

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

   ! Bound constraints

   lind(1:n) = .false.
   lbnd(1:n) = 0.0d0

   uind(1:n) = .false.
   ubnd(1:n) = 0.0d0

   ! Number equality (m) and inequality (p) constraints.
   
   m = 0
   p = 0

   allocate(lambda(m+p),c(m+p),stat=allocerr)

   if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
   end if

   ! Initial guess for the Lagrange multipliers
   
   lambda(1:m+p) = 0.0d0

   ! Number of entries in the JAcobian of the constraints
   
   jnnzmax = 10000

   ! This should be the number of entries in the Hessian of the
   ! Lagrangian. But, in fact, some extra space is need (to store the
   ! Hessian of the Augmented Lagrangian, whose size is hard to
   ! predict, and/or to store the Jacobian of the KKT system). Thus,
   ! declare it as large as possible.
   
   ! hlnnzmax = huge( 1 )
   hlnnzmax = 100000
   ! Feasibility, complementarity, and optimality tolerances
   
   epsfeas  = 1.0d-08
   epscompl = 1.0d-08
   epsopt   = 1.0d-08

   ! Maximum number of outer iterations
   
   maxoutit = 50

   ! rhoauto means that Algencan will automatically set the initial
   ! value of the penalty parameter. If you set rhoauto = .false. then
   ! you must set rhoini below with a meaningful value.
   rhoauto = .true.

   if ( .not. rhoauto ) then
      rhoini = 1.0d-08
   end if

   ! scale = .true. means that you allow Algencan to automatically
   ! scale the constraints. In any case, the feasibility tolerance
   ! (epsfeas) will be always satisfied by the UNSCALED original
   ! constraints.
   scale = .false.

   ! extallowed = .true. means that you allow Gencan (the active-set
   ! method used by Algencan to solve the bound-constrained
   ! subproblems) to perform extrapolations. This strategy may use
   ! extra evaluations of the objective function and the constraints
   ! per iterations; but it uses to provide overal savings. You should
   ! test both choices for the problem at hand.
   extallowed = .true.

   ! extallowed = .true. means that you allow the inertia of the
   ! Jacobian of the KKT system to be corrected during the acceleration
   ! process. You should test both choices for the problem at hand.
   corrin = .false.

   inf = 5
   sup = n_train

   noutliers = 0
   
   do samples = inf, sup
      allocate(t(samples),y(samples),indices(samples),sp_vector(samples),stat=allocerr)

      if ( allocerr .ne. 0 ) then
         write(*,*) 'Allocation error.'
         stop
      end if

      t(1:samples)         = (/(i, i = 1, samples)/)
      indices(1:samples)   = (/(i, i = 1, samples)/)
      y(1:samples)         = train_set(n_train - samples + 1:n_train)

      if (mod(dble(samples),7.d0) .eq. 0) then
         noutliers = noutliers + 1
      endif

      call lovo_algorithm(samples,n,lovo_order,noutliers,t,y,indices,sp_vector,grad_sp,gp)

      Open(Unit = 100, File = "output/solutions_covid_cubic.txt", ACCESS = "SEQUENTIAL")
      write(100,10) xk(1), xk(2), xk(3)

      deallocate(t,y,indices,sp_vector,stat=allocerr)
   
      if ( allocerr .ne. 0 ) then
         write(*,*) 'Deallocation error.'
         stop
      end if

   enddo

   10 format (ES13.6,1X,ES13.6,1X,ES13.6) 
   close(100)

   call cpu_time(start)

   call cpu_time(finish)
   
   deallocate(lind,lbnd,uind,ubnd,x,lambda,c,stat=allocerr)
   
   if ( allocerr .ne. 0 ) then
      write(*,*) 'Deallocation error.'
      stop
   end if
   
   stop
  
   contains

   ! *****************************************************************
   ! LOVO SUBROUTINES
   ! *****************************************************************

   subroutine lovo_algorithm(samples,n,lovo_order,noutliers,t,y,indices,sp_vector,grad_sp,gp)
      implicit none

      integer,       intent(in) :: samples,noutliers,n
      real(kind=8),  intent(in) :: t(samples),y(samples)
      integer,       intent(inout) :: lovo_order
      real(kind=8),  intent(inout) :: indices(samples),sp_vector(samples),grad_sp(n),gp(n)

      real(kind=8) :: sp,sigmin,epsilon,fxk,fxtrial,theta,alpha,gamma,termination
      integer :: iter,iter_sub,max_iter,max_iter_sub,i

      sigmin = 1.0d0
      epsilon = 1.0d-3
      alpha = 1.0d-8
      gamma = 1.0d+1
      max_iter = 1000
      max_iter_sub = 100
      lovo_order = samples - noutliers
      iter = 0
      iter_sub = 0

      xk(1:n) = 1.0d-2

      call compute_sp(samples,lovo_order,n,t,y,xk,indices,sp_vector,fxk)

      write(*,*) "--------------------------------------------------"
      write(*,10) "#iter","#init","Sp(xstar)","||g(xstar)||"
      10 format (2X,A5,4X,A5,6X,A9,7X,A12)
      write(*,*) "--------------------------------------------------"

      do
         iter = iter + 1

         call compute_grad_sp(samples,lovo_order,n,t,y,xk,indices,gp)

         ! do i = 1, n
         !    gp(i) = max(lbnd(i),min(xk(i) - gp(i),ubnd(i)))
         ! enddo

         termination = norm2(gp(1:n))

         write(*,20)  iter,iter_sub,fxk,termination
         20 format (I6,5X,I4,4X,ES14.6,3X,ES14.6)

         if (termination .lt. epsilon) exit
         if (iter .gt. max_iter) exit

         x(1:n) = xk(1:n)
         sigma = sigmin
        
         iter_sub = 1

         do

            call algencan(evalf,evalg,evalc,evalj,evalhl,jnnzmax,hlnnzmax, &
               n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
               scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
               outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))

            xtrial(:) = x(:)

            call compute_sp(samples,lovo_order,n,t,y,xtrial,indices,sp_vector,fxtrial)

            if (fxtrial .le. (fxk - alpha * norm2(xtrial(1:n-1) - xk(1:n-1))**2)) exit
            if (iter_sub .gt. max_iter_sub) exit

            sigma = gamma * sigma
            iter_sub = iter_sub + 1

         enddo

         fxk = fxtrial
         xk(:) = xtrial(:)

      enddo

      write(*,*) "--------------------------------------------------"

   end subroutine lovo_algorithm

   !*****************************************************************
   !*****************************************************************

   subroutine compute_sp(samples,lovo_order,n,t,y,x,indices,sp_vector,res)
      implicit none
      integer,       intent(in) :: samples,n,lovo_order
      real(kind=8),  intent(in) :: x(n),t(samples),y(samples)
      real(kind=8),  intent(inout) :: indices(samples),sp_vector(samples)
      real(kind=8),  intent(out) :: res

      integer :: i,kflag

      sp_vector(:) = 0.0d0
      kflag = 2
      indices(:) = (/(i, i = 1, samples)/)

      do i = 1, samples
         call fi(x,i,n,t,y,samples,sp_vector(i))
      end do
      
      ! Sorting
      call DSORT(sp_vector,indices,samples,kflag)

      ! Lovo function
      res = sum(sp_vector(1:lovo_order))

   end subroutine compute_sp

   !*****************************************************************
   !*****************************************************************

   subroutine compute_grad_sp(samples,lovo_order,n,t,y,x,indices,res)
      implicit none

      integer,       intent(in) :: samples,n,lovo_order
      real(kind=8),  intent(in) :: indices(samples),x(n),t(samples),y(samples)
      real(kind=8),  intent(out) :: res(n)

      real(kind=8) :: gaux,ti
      integer :: i,j
      
      res(:) = 0.0d0

      do i = 1, lovo_order
         ti = t(int(indices(i)))
         call model(x,int(indices(i)),n,t,y,samples,gaux)
         gaux = gaux - y(int(indices(i)))
         
         res(1) = res(1) + gaux * (ti - t(samples))
         res(2) = res(2) + gaux * ((ti - t(samples))**2)
         res(3) = res(3) + gaux * ((ti - t(samples))**3)
      enddo

   end subroutine compute_grad_sp

   !*****************************************************************
   !*****************************************************************

   subroutine model(x,i,n,t,y,samples,res)
      implicit none 

      integer,        intent(in) :: n,i,samples
      real(kind=8),   intent(in) :: x(n),t(samples),y(samples)
      real(kind=8),   intent(out) :: res

      res = y(samples) + x(1) * (t(i) - t(samples)) + &
            x(2) * ((t(i) - t(samples))**2) + x(3) * ((t(i) - t(samples))**3)

   end subroutine model

   !*****************************************************************
   !*****************************************************************

   subroutine fi(x,i,n,t,y,samples,res)
      implicit none

      integer,        intent(in) :: n,i,samples
      real(kind=8),   intent(in) :: x(n),t(samples),y(samples)
      real(kind=8),   intent(out) :: res
      
      call model(x,i,n,t,y,samples,res)
      res = res - y(i)
      res = 0.5d0 * (res**2)

   end subroutine fi

   !*****************************************************************
   !*****************************************************************

   subroutine regularized_taylor(samples,n,lovo_order,sigma,t,y,x,indices,sp_vector,grad_sp,res)

      implicit none

      integer,       intent(in) :: samples,n,lovo_order
      real(kind=8),  intent(in) :: x(n),sigma,t(samples),y(samples)
      real(kind=8),  intent(inout) :: indices(samples),sp_vector(samples),grad_sp(n)
      real(kind=8),  intent(out) :: res

      call compute_sp(samples,lovo_order,n,t,y,x,indices,sp_vector,res)
      call compute_grad_sp(samples,lovo_order,n,t,y,x,indices,grad_sp)
      res = res + dot_product(grad_sp(1:n),x(1:n) - xk(1:n))
      res = res + 0.5d0 * sigma * (norm2(x(1:n) - xk(1:n))**2)

   end subroutine regularized_Taylor

   !*****************************************************************
   !*****************************************************************

   subroutine grad_regularized_taylor(samples,n,lovo_order,sigma,t,y,x,indices,res)
      implicit none

      integer,       intent(in) :: samples,n,lovo_order
      real(kind=8),  intent(in) :: x(n),sigma,t(samples),y(samples)
      real(kind=8),  intent(inout) :: indices(samples)
      real(kind=8),  intent(out) :: res(n)
      
      call compute_grad_sp(samples,lovo_order,n,t,y,x,indices,res)
      res = res + sigma * (x(1:n) - xk(1:n))

   end subroutine grad_regularized_taylor

   ! *****************************************************************
   ! ALGENCAN SUBROUTINES
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

   ! This routine must compute the objective function.
   
   ! LOCAL SCALARS
   type(pdata_type), pointer :: pdata
   
   call c_f_pointer(pdataptr,pdata)
   pdata%counters(1) = pdata%counters(1) + 1
   
   call regularized_taylor(samples,n,lovo_order,sigma,t,y,x,indices,sp_vector,grad_sp,f)
   
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
   
   ! This routine must compute the gradient of the objective
   ! function.
   
   ! LOCAL SCALARS
   type(pdata_type), pointer :: pdata
   
   call c_f_pointer(pdataptr,pdata)
   pdata%counters(2) = pdata%counters(2) + 1
   
   call grad_regularized_taylor(samples,n,lovo_order,sigma,t,y,x,indices,g)

   end subroutine evalg

   ! *****************************************************************
   ! *****************************************************************

   subroutine evalc(n,x,m,p,c,inform,pdataptr)

   implicit none

   ! SCALAR ARGUMENTS
   integer, intent(in) :: m,n,p
   integer, intent(inout) :: inform
   type(c_ptr), optional, intent(in) :: pdataptr

   ! ARRAY ARGUMENTS
   real(kind=8), intent(in) :: x(n)
   real(kind=8), intent(out) :: c(m+p)

   ! This routine must compute all the m+p constraints.
   
   ! LOCAL SCALARS
   type(pdata_type), pointer :: pdata
   
   
   end subroutine evalc

   ! *****************************************************************
   ! *****************************************************************

   subroutine evalj(n,x,m,p,ind,sorted,jsta,jlen,lim,jvar,jval,inform,pdataptr)

   implicit none
   
   ! SCALAR ARGUMENTS
   integer, intent(in) :: lim,m,n,p
   integer, intent(inout) :: inform
   type(c_ptr), optional, intent(in) :: pdataptr

   ! ARRAY ARGUMENTS
   logical, intent(in) :: ind(m+p)
   real(kind=8), intent(in) :: x(n)
   logical, intent(out) :: sorted(m+p)
   integer, intent(out) :: jsta(m+p),jlen(m+p),jvar(lim)
   real(kind=8), intent(out) :: jval(lim)
   
   ! This routine must compute the Jacobian of the constraints. In
   ! fact, only gradients of constraints j such that ind(j) =
   ! .true. need to be computed.
   
   ! LOCAL SCALARS
   integer :: i
   type(pdata_type), pointer :: pdata
   
   
   end subroutine evalj

   ! *****************************************************************
   ! *****************************************************************

   subroutine evalhl(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform,pdataptr)

   implicit none
   
   ! SCALAR ARGUMENTS
   logical, intent(in) :: inclf
   integer, intent(in) :: m,n,lim,p
   integer, intent(out) :: hlnnz
   integer, intent(inout) :: inform
   type(c_ptr), optional, intent(in) :: pdataptr

   ! ARRAY ARGUMENTS
   real(kind=8), intent(in) :: lambda(m+p),x(n)
   integer, intent(out) :: hlrow(lim),hlcol(lim)
   real(kind=8), intent(out) :: hlval(lim)

   ! This routine must compute the Hessian of the Lagrangian. The
   ! Hessian of the objective function must NOT be included if inclf
   ! = .false.
   
   ! LOCAL SCALARS
   type(pdata_type), pointer :: pdata
   
   call c_f_pointer(pdataptr,pdata)
   pdata%counters(5) = pdata%counters(5) + 1

   hlnnz = 0

   ! If .not. inclf then the Hessian of the objective function must not be included
   
   if ( inclf ) then
      hlnnz = n

      if ( hlnnz .gt. lim ) then
          inform = -95
          return
      end if
      
      hlrow(1:n) = (/(i, i = 1, n)/)
      hlcol(1:n) = (/(i, i = 1, n)/)
      hlval(1:n) = sigma

   end if

   ! Note that entries of the Hessian of the Lagrangian can be
   ! repeated. If this is case, them sum of repeated entrances is
   ! considered. This feature simplifies the construction of the
   ! Hessian of the Lagrangian.
   
   if ( hlnnz + 1 .gt. lim ) then
      inform = -95
      return
   end if
   
   end subroutine evalhl
  
end program algencama