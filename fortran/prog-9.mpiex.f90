! Homework 4: Due Tuesday, March 9, 2021      
! Name: Yonghyun Kwon


!  Note: A and B will have different values on diffent processors
!  Initialize A(1:n)  to be 1.d0 on processor 0 and B(1:n) to be 
!  dble(rank) on the rank processor.

! Processor 0 sends A to all other processors.
! Processor i > 0 receives A in A and performs A = A + B and sends
! the new A to processor 0 receiving into B where it is then added to A.  
! When all messages have been received by processor 0, processor 0 
! prints 'Messages received'

  use mpi 
  implicit none
  integer, parameter ::  n = 2, tag = 1
  integer, parameter ::  dp = mpi_double_precision, comm=mpi_comm_world
  double precision   ::  A(n), B(n)
  integer            ::  i, ierror
  integer            ::  p, rank, status(mpi_status_size)
  call mpi_init(ierror) 
  call mpi_comm_size(comm, p, ierror)
  call mpi_comm_rank(comm, rank, ierror)
  
  if(rank == 0) then
    A = 1.d0
    do i = 1, p-1
      call mpi_send(A(1), n, dp, i, tag, comm, ierror) ! dest = i
    enddo
  endif

  if(rank > 0) then
    B = dble(rank)
    call mpi_recv(A(1), n, dp, 0, tag, comm, mpi_status_ignore, ierror)
    A = A + B 
    call mpi_send(A(1), n, dp, 0, tag, comm, ierror)
   endif

  if(rank == 0) then
    do i = 1, p-1
      call mpi_recv(B(1), n, dp, i, tag, comm, mpi_status_ignore, ierror)
      A = A + B
    enddo
    print*, 'Messages received'
  endif    
  
  call mpi_finalize(ierror)
  end

!  Output
!  [yhkwon@hpc-class fortran]$ mpiifort prog-9.mpiex.f90
!  [yhkwon@hpc-class fortran]$ mpirun -np 16 ./a.out
!   Messages received
