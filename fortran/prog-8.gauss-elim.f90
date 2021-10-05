!  Due February 25, 2021
!  Student's Name: Yonghyun Kwon   

!  Hw#2:  Solving inear systems with multiple right-hand sides

!  Use omp-lib and the OpenMP timer:  t = omp_get_wtime()
!  Compile as ifort -qopenmp prog-8.gauss-elim.f90
!  Present results for m = 1024 and k = 1024.  You will have to
!  Use  salloc -N 1 -n 16 -t 1:00:0 due to limited memory on hpc-class
!  For program development use m = 512 and k = 512

!  Solve the linear system AX = C with possible multiple right-hand
!  sides using the method described below.  To compute the actual
!  error, we first initalize A and X with a random number generator
!  and then let  C = AX so X will the the exact solution. A is m x m.
!  X is m x k and C is m x k.  Let B be the augmented matrix formed from
!  A
!  and C.  Thus, B = m x (m+k) = m x n, where n = m+k. 

!  Use the following procedure on B
!  1.  Find the row that has the element of largest absolute value in
!  column 1
!      and interchange this row with row 1 when abs(B(1,1)) is not the
!      largest
!      absolute value in column 1.  Next divide the (possible) new row 1
!      by the new
!      B(1,1).  (Thus the updated value of B(1,1) = 1.d0.)  B(1,1) is
!      called
!      the "pivot element".
!  2.  Zero all the elements of B that are ABOVE AND BELOW
!      the pivot element using elementary row operations.
!  3.  Next find the pivot element for the second column by finding the
!      largest element in absolute value in B(2:m,2) and interchanging
!      rows if needed to make B(2,2) the largest in absolute value.
!  4.  Next zero the elements ABOVE AND BELOW B(2,2) usng elementary
!      row operations.
!  5.  Etc.

!  Print the time and gigaflops required to compute the solution X.
!  Calculate and print the max abs of the RESIDUAL = matmul(A,X) - C
!  Calculate and print the max abs of the actual error XE - X where X 
!  is the computed solution.

   use omp_lib
   implicit none
   integer, parameter ::  m =1024, k=1024, n=m+k ! m = 1024, k = 1024
   integer            ::  i, j, kk, index
   double precision   :: A(m,m), X(m,k), XE(m,k), B(m,n), C(m,k)
   double precision   :: RESIDUAL(m,k), t0, t1, time, num_ops, gflops
   double precision   :: D(m,m), E(m,m), temp, tmp(1:n)

!  initialze arrays so the error can be calculated
    call random_number(A)
    B(1:m,1:m) = A
    call random_number(XE) ! XE is the "exact" solution
    C = matmul(A,XE)
    B(1:m,(m+1):(m+k)) = C

    t0 = omp_get_wtime()

!  Transform B(1:m,1:m) to the identity matrix using elementary row ops.
   do i = 1, m
!     Note:  if i=m, then the pivot will be B(m,m), so we can skip
!     everyting
!     except for dividing the m-th row by B(m,m).
      if (i < m) then
!        Find element, B(index,i), of largest absolute value in B(i:m,i)
         index = i 
         do j = i, m            
            if (abs(B(j, i)) > abs(B(index, i))) then
               index = j
            endif                                             
         enddo
!        interchange rows i and index from columns i+1 to n ! columns i to n?
         tmp(i:n) = B(i, i:n)
         B(i, i:n) = B(index, i:n)
         B(index, i:n) = tmp(i:n)       
      endif
!     divide row i by B(i,i) from columns i to n
      temp = B(i, i)
      do j = i, n
         B(i, j) = B(i, j) / temp
      enddo

 !    end of find_pivot(B m, n, i)
      do j = 1, m                                                   
         if (j /= i) then
            temp = B(j, i)
            do kk = i, n
               B(j, kk) = B(j, kk) - temp * B(i, kk)
            enddo                  
            !B(j, i:n) = B(j, i:n) - temp * B(i, i:n)                
         endif                                    
      enddo           
                     
   enddo

   t1 = omp_get_wtime()
   time = t1 - t0 ! time in seconds

!  The number of floating point operations is estimated as:
   num_ops = dble(k)*dble(m)*dble(m)*dble(m)/3.d0
   gflops = (num_ops/time)*1.d-9
   print*,'Performance = ', gflops,' gflops'
   print*,'Time  = ', time,' seconds'

   X(1:m,1:k) = B(1:m,(m+1):(m+k)) ! computed solution of AX = C

   print*,'m = ', m,' k = ', k,' = # of right-hand sides'
   RESIDUAL = matmul(A,X) - C
   print*,'RESIDUAL      = ', maxval(abs(RESIDUAL))
   print*,'maximum ERROR = ', maxval(abs(X-XE))
   end
       
! OUTPUT
! [yhkwon@hpc-class04 fortran]$ ifort -qopenmp prog-8.gauss-elim.f90
! [yhkwon@hpc-class04 fortran]$ ./a.out
!  Performance =    22.9423353367673       gflops
!  Time  =    15.9750030040741       seconds
!  m =         1024  k =         1024  = # of right-hand sides
!  RESIDUAL      =   1.364242052659392E-012
!  maximum ERROR =   4.462652469783279E-012
