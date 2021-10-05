!  Hw9.  Jacobi iteration:  serial, sendrecv, nonblocking,&
!  persistent(extra credit)
!  Due Thursday April 15, 2021 by 11:59 pm.                           
 
!  Name: Yonghyun Kwon
 
!  Write one program with (a) serial Jacobi, (b) sendrecv Jacobi,
!  and (c) nonblocking Jacobi. (d) persistent Jacobi for extra 
!  credit.  Present results for p = 16, 32 and 64.
 
!  Remarks:  Jocobi's iteration that will be presented comes from
!  solving the 2D Laplace equation in a rectangle with values of
!  the solution known on the 4 boundaries.  We shall use the 2D
!  array A to denote approximations to the solution.
!  The domain is a <= x <= b and c <= y <= d.
!  The boundary conditions are 1 on left boundary and 0 elsewhere.
!  x(i) = a + i*(b-a)/(n+1), for i = 0, 1, ...., n+1.
!  y(j) = c + j*(d-c)/(m+1), for j = 0, 1, ...., m+1 
!  Let U(x,y) denote the solution and let
!  A(i,j) = U(x(j),y(i)) the boundary condition mean that
!  A(1:n,0) = 1.d0
!  A(n+1,0:m+1) = 0.d0
!  A(0:n+1,0) = 0.d0
!  A(0:n+1,m+1) = 0.d0
!  Initialize the interior points of A to be zero.
 
  use mpi
  implicit none
  integer, parameter :: n=4096, niter=50, tag=1
  integer :: i, j, m, iter, ierror, left, right, req(4)
  double precision :: t1, t2, time, max_time
  double precision,allocatable :: A(:,:), B(:,:)
  integer, parameter :: dp=mpi_double_precision, comm = mpi_comm_world
  integer :: p, rank, status(mpi_status_size),&
array_of_statuses(mpi_status_size,4)
  call mpi_init(ierror)
  call mpi_comm_size(comm, p, ierror)
  call mpi_comm_rank(comm, rank, ierror)
 
 
! SERIAL JACOBI ITERATION
  if (rank == 0) then
     allocate (A(0:n+1,0:n+1), B(n,n))
!    Initialize A
     A = 0.d0
     A(1:n,0) = 1.d0
 
     t1 = mpi_wtime()
     do iter = 1, niter
       do j = 1, n
         do i = 1, n
           B(i, j) = 0.25d0 * (A(i-1, j) + A(i+1, j) + &
     A(i, j-1) + A(i, j+1) ) 
         enddo
       enddo                                                        
       A(1:n, 1:n) = B(1:n, 1:n)     
     enddo       
     t2 = mpi_wtime()
     time = (t2 - t1)/dble(niter)
     print *,'FOR JACOBI ITERATION with p, n and niter = ', p, n, niter
     print *,'Serial Jacobi Iteration'
     print *,'average time per iteration = ', time,' seconds'
     !print *, A(0:2, 0:2)
     deallocate (A, B)
  end if
  call mpi_barrier(comm, ierror)
 
! INITIALIZATION FOR THE PARALLEL JACOBI ITERATION METHODS
!  Use mpi_proc_null to determine the left and right neighbors when
!  dest or source = mpi_proc_null, then the routine does nothing.
   if (rank == 0) then
      left = mpi_proc_null
   else
      left = rank - 1
   endif
   if (rank == p-1) then
      right = mpi_proc_null
   else
      right = rank + 1
   endif
 
!  Compute the size of the local blocks and allocate A and B
   m = n/p
   if (rank < (n-p*m)) then
      m = m + 1
   endif ! if n = m*p + r, 0 < r < p, then m = n/p + 1 for
         ! rank = 0, 1, ..., r-1.
   allocate (A(0:n+1,0:m+1), B(n,m))
   call mpi_barrier(comm, ierror)
  
! SENDRECV JACOBI ITERATION
!    Initialize A
     A = 0.d0
     if (rank == 0) then
        A(1:n,0) = 1.d0
     endif
     call mpi_barrier(comm, ierror)
     t1 = mpi_wtime()
     do iter = 1, niter
       call mpi_sendrecv(B(1, 1), n, dp, right, 1, &
         A(1, m+1), n, dp, left, 1, &
         comm, mpi_status_ignore, ierror)
       call mpi_sendrecv(B(1, m), n, dp, left, 1, &
         A(1, 0), n, dp, right, 1, &
         comm, mpi_status_ignore, ierror)
       do j = 1, m
         do i = 1, n
           B(i, j) = 0.25d0 * (A(i - 1, j) + A(i + 1, j) + &
             A(i, j - 1) + A(i, j + 1))
         enddo               
       enddo
       A(1:n, 1:m) = B(1:n, 1:m)                 
     enddo        
     t2 = mpi_wtime()
     time = t2 - t1
     call mpi_reduce(time, max_time, 1, dp, mpi_max, 0, comm, ierror)
     if (rank == 0) then
        print *,'  '
        print *,'Jacobi Iteration using mpi_sendrecv'
        print *,'average time per iteration = ', max_time/dble(niter),&
'seconds'
        !print *, A(0:2, 0:2)
     endif
     call mpi_barrier(comm, ierror)
 
 
! NONBLOCKING JACOBI ITERATION
!    Initialize A
     A = 0.d0
     if (rank == 0) then
        A(1:n,0) = 1.d0
     endif
     call mpi_barrier(comm, ierror)
     t1 = mpi_wtime()
     do iter = 1, niter 
       call mpi_irecv(A(1, m+1), n, dp, right, 1, comm, req(1), ierror)
       call mpi_irecv(A(1, 0), n, dp, left, 1, comm, req(2), ierror)
!       Start communication posting receives first so messages can
!       be sent without waiting for a receive to be posted.
       call mpi_isend(B(1, 1), n, dp, left, 1, comm, req(3), ierror)  
       call mpi_isend(B(1, m), n, dp, right, 1, comm, req(4), ierror)               
!       Compute interior
       do j = 1, m
         do i = 1, n
           B(i, j) = 0.25d0 * (A(i - 1, j) + A(i + 1, j) + &
             A(i, j - 1) + A(i, j + 1))
         enddo
       enddo
       A(1:n, 1:m) = B(1:n, 1:m)    
!       Complete communication
       call mpi_waitall(4, req, array_of_statuses, ierror)
     enddo           
     t2 = mpi_wtime()
     time = t2 - t1
     call mpi_reduce(time, max_time, 1, dp, mpi_max, 0, comm, ierror)
     if (rank == 0) then
        print *,'  '
        print *,'Jacobi Iteration using mpi_isends and mpi_irecvs'
        print *,'average time per iteration = ', max_time/dble(niter), &
'seconds'
        !print *, A(0:2, 0:2)
     endif
     call mpi_barrier(comm, ierror)
 
 
! PERSISTENT JACOBI ITERATION - EXTRA CREDIT
!    Initialize A
     A = 0.d0
     if (rank == 0) then
        A(1:n,0) = 1.d0
     endif
!    Create persistent requests
                                                                      
                                                                        
                                                                        
                                                                        
 
     call mpi_barrier(comm, ierror)
     t1 = mpi_wtime()
     do iter = 1, niter
!       Compute boundary columns to allow for fast communication
                           
                                                                 
                                                                   
                
 
!       Start communication posting receives first so messages can
!       be sent without waiting for a receive to be posted.
        !call mpi_startall(4, req, ierror)
 
!       Compute interior
                     
                       
                                                                      
                  
                
                      
                         
                               
                  
                     
!       Complete communication
        !call mpi_waitall(. . .                             
     enddo        
     t2 = mpi_wtime()
     time = t2 - t1
     call mpi_reduce(time, max_time, 1, dp, mpi_max, 0, comm, ierror)
     if (rank == 0) then
        print *,' '
        print *,'Jacobi Iteration using persistent nonblocking routines'
        print *,'average time per iteration = ', max_time/dble(niter), &
 'seconds'
     endif
 
  call mpi_barrier(comm, ierror)
  deallocate(A, B)
  call mpi_finalize(ierror)
  end 
  
!  Results:                       

! FOR JACOBI ITERATION with p, n and niter = 16, 4096, 50
! Serial Jacobi Iteration
! average time per iteration =   5.535970211029053E-002  seconds

! Jacobi Iteration using mpi_sendrecv
! average time per iteration =   1.082376003265381E-002 seconds

! Jacobi Iteration using mpi_isends and mpi_irecvs
! average time per iteration =   1.079095840454102E-002 seconds
! _______________________________________________________________
 
!  FOR JACOBI ITERATION with p, n and niter = 32, 4096, 50
!  Serial Jacobi Iteration
!  average time per iteration =   5.523116111755371E-002  seconds
   
!  Jacobi Iteration using mpi_sendrecv
!  average time per iteration =   5.463256835937500E-003 seconds
   
!  Jacobi Iteration using mpi_isends and mpi_irecvs
!  average time per iteration =   5.466842651367187E-003 seconds
! __________________________________________________________________
 
! FOR JACOBI ITERATION with p, n and niter = 64, 4096, 50
! Serial Jacobi Iteration
! average time per iteration =   5.542436122894287E-002  seconds
   
! Jacobi Iteration using mpi_sendrecv
! average time per iteration =   2.698583602905274E-003 seconds
   
! Jacobi Iteration using mpi_isends and mpi_irecvs
! average time per iteration =   2.713718414306641E-003 seconds

