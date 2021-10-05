
! MPI Homework 8-a: Due April 8, 2021                                          
! MPI Homework 8-b: Due April 8, 2021: Same as 8-a but use mpi_sendrecv
!                   to remove the deadlock.
! MPI Homework 8-c: Due April 8, 2021: Same as 8-a but use nonblocking 
!                   sends and receives to remove the deadlock.
 
!  Name: Yonghyun Kwon
 
! DESCRIPTION FOR HW 8-a:
! Verify that the circular right shift run okay when using mpi_sends    
! and mpi_recvs when n < 65 and deadlocks when n >95 when p = 16.
! Show that the program deadlocks for n = 64 when replacing mpi_send
! with mpi_ssend.  To document this you will have to write on your
! program when your program deadlocks.  Use p = 16 for all experiments.
 
! The Program:                                    
! On each processor initialize A using random_number(A) + dble(rank) and
! double precision  ::  A(n,n), B(n,n), C(n,n), D(n,n)
! Processor of rank i sends A to B on processor i+1, receives
! A from processor i-1 into B, computes C = matmul(A,B) and
! then sends C into D on processor i-1.  When i = p-1 "i+1" becomes
! 0 and when i = 0, "i-1" becomes p-1. After completing the
! above, all non-zero processors send a message to processor 0
! saying that the above tasks have been completed.
! Processor zero will then print:                                  
! "All processors have completed their tasks".
 
  use mpi 
  implicit none
  integer, parameter :: n =  64, ntrial =  128
  integer, parameter :: dp = mpi_double_precision, comm=mpi_comm_world
  double precision   :: t0, t1, time(ntrial), max_time(ntrial)
  double precision   :: x, A(n,n), B(n,n), C(n,n), D(n,n)
  integer            :: i, j, ierror, left, right
  integer            :: p, rank, status(mpi_status_size), ktrial
  call mpi_init(ierror)
  call mpi_comm_size(comm, p, ierror)
  call mpi_comm_rank(comm, rank, ierror)
 
  call random_number(A)
  A = A + dble(rank)
 
! To keep track of the values of rank+1 and rank-1 set
  left = rank-1
  right = rank+1
  if (rank == 0) left = p-1
  if (rank == p-1) right = 0
 
do ktrial = 0, ntrial
   call mpi_barrier(comm, ierror)
   t0 = mpi_wtime()
 
!****************************************************************
! perform calculations
  call mpi_send(A(1, 1), n*n, dp, right, 1, comm, ierror) ! deadblock for ssend
  call mpi_recv(B(1, 1), n*n, dp, left, 1, comm, mpi_status_ignore, &
  ierror)
  C = matmul(A, B)
  call mpi_send(C(1, 1), n*n, dp, left, 1, comm, ierror) ! deadblock for ssend
  call mpi_recv(D(1, 1), n*n, dp, right, 1, comm, mpi_status_ignore, &
  ierror)
  if(rank > 0) then
    call mpi_send(A(1, 1), 1, dp, 0, 1, comm, ierror)
  endif
  if(rank == 0) then
    do i = 1, p-1
      call mpi_recv(B(1, 1), 1, dp, i, 1, comm, mpi_status_ignore, ierror)
    enddo           
  endif
!****************************************************************
 
  t1 = mpi_wtime()
  time(ktrial) = t1 - t0  ! seconds
enddo ! for the ktrial loop
 
  call mpi_reduce(time, max_time, ntrial, dp, mpi_max, 0, comm, ierror)
 
  if (rank == 0) then
      print*,'All tasks have been completed'
      print*,'Times in seconds and n =', n,' and p =', p
      print*,'average time =',sum(max_time(1:ntrial))/dble(ntrial)
  endif
 
  call mpi_finalize(ierror)
  end
 
! OUTPUT

! [yhkwon@hpc-class04 fortran]$ mpirun -np 16 ./a.out
! export TMPDIR=/scratch/yhkwon/27264
! All tasks have been completed
! Times in seconds and n =          64  and p =          16
! average time =  1.383144408464432E-004
 
! Deadblocks when n = 100
! [yhkwon@hpc-class04 fortran]$ mpirun -np 16 ./a.out
! export TMPDIR=/scratch/yhkwon/27264

! Deadblocks for n = 64 when replacing mpi_send with mpi_ssend
! [yhkwon@hpc-class04 fortran]$ mpirun -np 16 ./a.out
! export TMPDIR=/scratch/yhkwon/27264
