!   Programmer:  G. Luecke
!   Compile as:  mpiifort prog.f90
!   Execute as:  mpirun -np 1 ./a.out

    program function
    use mpi
    implicit none
    integer,parameter ::  n = 21
    integer  ::  fac, i, ierror

    integer  ::  f !  must declare type for functions if there is no
                   !  interface block or module for the function

    call mpi_init(ierror)
    
   do i = 1, n
      call factorial(i,fac)
      print*,i,'factorial = ', fac, f(i)
   enddo
   call mpi_finalize(ierror)
   end program function


   subroutine factorial(n,ans)
   implicit none
   integer,intent(in)  :: n ! subroutine cannot change n
   integer,intent(out) :: ans ! output value only
   integer  ::  i
   ans = 1
   do i = 2, n
      ans = ans*i
   enddo
   end subroutine factorial

   function f(n)
   implicit none
   integer  ::  f, i
   integer,intent(in)  ::  n
   f = 1 
   do i = 2, n
      f = f*i
   enddo
   end function f

! Output
! [grl@hpc-class fortran]$ mpiifort  prog-6.fcn.sub.factorial.f90 
! [grl@hpc-class fortran]$ ./a.out
!            1 factorial =            1           1
!            2 factorial =            2           2
!            3 factorial =            6           6
!            4 factorial =           24          24
!            5 factorial =          120         120
!            6 factorial =          720         720
!            7 factorial =         5040        5040
!            8 factorial =        40320       40320
!            9 factorial =       362880      362880
!           10 factorial =      3628800     3628800
!           11 factorial =     39916800    39916800
!           12 factorial =    479001600   479001600
!           13 factorial =   1932053504  1932053504
!           14 factorial =   1278945280  1278945280
!           15 factorial =   2004310016  2004310016
!           16 factorial =   2004189184  2004189184
!           17 factorial =   -288522240  -288522240
!           18 factorial =   -898433024  -898433024
!           19 factorial =    109641728   109641728
!           20 factorial =  -2102132736 -2102132736
!           21 factorial =  -1195114496 -1195114496
