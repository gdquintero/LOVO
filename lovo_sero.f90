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
   integer :: samples
   real(kind=8), allocatable :: xtrial(:),xk(:)

   ! Number of variables

   n = 3
   
   allocate(x(n),lind(n),lbnd(n),uind(n),ubnd(n),stat=allocerr)

   if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
   end if

   allocate(xtrial(n),xk(n),stat=allocerr)

   if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
   end if

   ! Initial guess and bound constraints
   
   x(1:n) = 10.0d0

   lind(1:n) = .true.
   lbnd(1:n) = - 100.0d0

   uind(1:n) = .true.
   ubnd(1:n) =   100.0d0

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
   
   jnnzmax = n

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

   call cpu_time(start)

   ! call algencan(evalf,evalg,evalc,evalj,evalhl,jnnzmax,hlnnzmax, &
   !    n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
   !    scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
   !    outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))

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

   subroutine lovo_algorithm()
      implicit none

   end subroutine lovo_algorithm

   !*****************************************************************
   ! MODEL TO BE FITTED TO THE DATA
   !*****************************************************************
   subroutine model(x,i,n,t,samples,res)
      implicit none 

      integer,        intent(in) :: n,i,samples
      real(kind=8),   intent(in) :: x(n-1),t(samples)
      real(kind=8),   intent(out) :: res
      real(kind=8) :: a,b,c,ti,ebt

      a = x(1)
      b = x(2)
      c = x(3)
      ti = t(i)
      ebt = exp(-b * ti)

      res = (a / b) * ti * ebt
      res = res + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) 
      res = 1.0d0 - exp(res - c * ti)

  end subroutine model

   subroutine regularized_Taylor(x,n,ind_train,nuk,sigma,res)

      implicit none

      integer,        intent(in) :: n,nuk,ind_train
      real(kind=8),   intent(in) :: x(n),sigma
      real(kind=8),   intent(out) :: res

      ! call compute_Fmin(x,n,ind_train,res)
      ! call compute_grad_Fi(x,n,nuk,grad_Fi)
      ! res = res + dot_product(grad_Fi,x(1:n) - xk(1:n))
      ! res = res + 0.5d0 * sigma * (norm2(x(1:n) - xk(1:n))**2)

  end subroutine regularized_Taylor

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
   
   f = ( x(1) + 4.0d0 ) ** 4.0d0 + x(2) ** 2.0d0
   
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
   
   g(1) = 4.0d0 * ( x(1) + 4.0d0 ) ** 3.0d0
   g(2) = 2.0d0 * x(2)

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
   
   call c_f_pointer(pdataptr,pdata)
   pdata%counters(3) = pdata%counters(3) + 1
   
   c(1) = x(1) ** 3.0d0 - x(2) - 1.0d0
   
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
   
   call c_f_pointer(pdataptr,pdata)
   pdata%counters(4) = pdata%counters(4) + 1

   ! Only gradients of constraints j such that ind(j) = .true. need
   ! to be computed.
   
   if ( ind(1) ) then
      if ( lim .lt. n ) then
         inform = -94
         return
      end if
      
      jsta(1) = 1
      jlen(1) = n
      
      jvar(1:n) = (/ (i,i=1,n) /)

      jval(1) = 3.0d0 * x(1) ** 2.0d0
      jval(2) = - 1.0d0

      ! Says whether the variables' indices in jvar (related to this
      ! constraint) are in increasing order. In case they are,
      ! Algencan takes advantage of this. Implement sorted gradients
      ! of constraints if you can do this in a natural (cheap)
      ! way. Under no circumnstance use a sorting algorithm. (It is
      ! better to set sorted(1) = .false. in this case.)
      
      sorted(1) = .true.
   end if
   
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
      if ( hlnnz + 2 .gt. lim ) then
         inform = -95
         return
      end if
   
      hlnnz = hlnnz + 1
      
      hlrow(hlnnz) = 1
      hlcol(hlnnz) = 1
      hlval(hlnnz) = 12.0d0 * ( x(1) + 4.0d0 ) ** 2.0d0
   
      hlnnz = hlnnz + 1
      
      hlrow(hlnnz) = 2
      hlcol(hlnnz) = 2
      hlval(hlnnz) = 2.0d0
   end if

   ! Note that entries of the Hessian of the Lagrangian can be
   ! repeated. If this is case, them sum of repeated entrances is
   ! considered. This feature simplifies the construction of the
   ! Hessian of the Lagrangian.
   
   if ( hlnnz + 1 .gt. lim ) then
      inform = -95
      return
   end if
   
   hlnnz = hlnnz + 1
   
   hlrow(hlnnz) = 1
   hlcol(hlnnz) = 1
   hlval(hlnnz) = lambda(1) * 6.0d0 * x(1)
   
   end subroutine evalhl
  
  end program algencama