!  Homework #6, spring 2021  Due:  March  25
 
!  Name: Yonghyun Kwon

!  Use mpi_bcast for the following:                            
 
!  Initialize A only on process 0 and initialize B on all 
!  processors. From processor 0, broadcast A to all other
!  processors using mpi_bcast
!  Calculate A = A + matmul(A,B) on all processors of rank > 0..
!  Then for all processors of rank > 0 send a message to
!  processor 0 to indicate its computation has completed,
!  print the average time for the calculation.         
 
!   mpiifort prog.f90, mpirun -np 8  ./a.out 
 
  use mpi
  implicit none
  integer,parameter :: n = 512, tag = 1, ntrial = 16 
  integer,parameter :: dp = mpi_double_precision, comm=mpi_comm_world
  double precision  :: A(n,n), B(n,n), t1, t2, time(1:ntrial), average
  integer           :: i, ierror, p, rank, status(mpi_status_size),ktrial
  integer,allocatable :: req(:), array_of_statuses(:, :)
  call mpi_init(ierror)
  call mpi_comm_size(comm, p, ierror)
  call mpi_comm_rank(comm, rank, ierror)
 
! Initialize the array A on each MPI processes
  if(rank==0) call random_number(A)
  B = dble(rank)

  do ktrial = 1, ntrial
    call mpi_barrier(comm, ierror)                             
    t1 = mpi_wtime()                                
    call mpi_bcast(A(1, 1), n*n, dp, 0, comm, ierror) 
    if (rank > 0) then                                                 
      A = A + matmul(A, B)
      call mpi_send(A(1, 1), 1, dp, 0, tag, comm, ierror) !send a message to rk0
    endif
    if (rank == 0) then
      do i = 1, p-1
        call mpi_recv(B(1, 1), 1, dp, i, tag, comm, mpi_status_ignore, ierror)
!receive a message if compuation completed
        !print*, 'i = ', i, 'A(1, 1)', A(1, 1)
      enddo                           
    endif
    t2 = mpi_wtime()                        
    time(ktrial) = t2 - t1                                                        
  enddo ! for the ktrial loop

 if (rank == 0) then
    print*,'n =',n,'ntrial =',ntrial,'p =', p
    average = sum(time(1:ntrial))/dble(ntrial)
    print*,'average time = ', average,'seconds'
 endif
 
  call mpi_finalize(ierror)
  end

! OUTPUT
![yhkwon@hpc-class06 fortran]$ mpiifort prog-5.hw6.bcast.f90
![yhkwon@hpc-class06 fortran]$ mpirun -np 8 ./a.out
!export TMPDIR=/scratch/yhkwon/26691
! n =         512 ntrial =          16 p =           8
! average time =   7.400928437709808E-002 seconds

