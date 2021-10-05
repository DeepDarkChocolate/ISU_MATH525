!    Student's Name: Yonghyun Kwon
!    Hw#3  Prog-9.integration.f90: Serial Integration, due March 4, 2021          

!     Compile as:  ifort -qopenmp program.f90  to use omp_get_wtime()
!     Execute as:  ./a.out

!     Write a program to approximate the value of the integral of
!     f(x) for x in the interval [a,b] using the composite
!     trapezoidal rule of order n. To check answers let
!     f(x) = x*x, a = 0.d0, b = 3.d0 and n = 128, 256 and 512.
!     For the program use "f" where "f" is defined as a
!     "double precision function" after the end of the main program.
!     The values of a and b should be set in your program so they can
!     be easily modified.  You should calculate the error and the error
!     should be in terms of a and b.

!     For the composit trapezoidal rule [a,b] is divided into n
!     subintervals
!     [x(i-1),x(i)], where x(i) = a + i*(b-a)/n, i = 0,1,..., n.        
!     Let h = (b-a)/n  so that x(i) = a + i*h, i = 0,1,..., n 
!     then value of the integral of f(x) from a to b is approximated by
!     the sums of the areas of the trapezoids for each subinterval.
!     Compute and print the average time when  ntrial = 16.

   use omp_lib
   implicit none
   double precision  ::  a = 0.d0, b = 3.d0
   integer,parameter ::  ntrial = 16, n = 512
   integer           ::  i, ktrial
   double precision  :: t1, t2, time(ntrial), psum, f, h, answer

!  Print a and b:
   print*,'a = ',a,' b = ',b

!  Calculate the integral using the composite trapezoidal rule

   do ktrial = 0, ntrial
      t1 = omp_get_wtime() 
      ! Here we use the following formula
      ! int_a^b f(x)dx = [(f(a) + f(b))/2 + sum_{i = 1}^{n-1}f(x_i)] * h
      h = (b - a) / dble(n)
      psum = (f(a) + f(b)) / 2.d0
      do i = 1, n-1
         psum = psum + f(a + dble(i) * h)
      enddo
      psum = psum * h 
                                 
      t2 = omp_get_wtime()
      time(ktrial) = t2 - t1 ! time in seconds
   enddo  
 
   print*,'psum        = ', psum
   answer =  (b*b*b - a*a*a)/3.d0
   print*,'exact answer = ', answer
   print*,'error        = ', psum - answer
   print*,'n            = ', n
   print*,'average time = ', sum(time)/dble(ntrial),' seconds'
   end

   function f(x)
   double precision :: f, x
   f = x*x
   return
   end

!  OUTPUT
!  [yhkwon@hpc-class fortran]$ ifort -qopenmp prog-9.integration.f90
!  [yhkwon@hpc-class fortran]$ ./a.out
!  a =   0.000000000000000E+000  b =    3.00000000000000     
!  psum        =    8.94789897417650
!  exact answer =    9.00000000000000
!  error        =  -5.210102582350373E-002
!  n            =          512
!  average time =   3.129243850708008E-007  seconds
