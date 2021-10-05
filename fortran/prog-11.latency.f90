! Homework 5:  Due:   Thursday, March 18, 2021

! Name: Yonghyun Kwon                    

! Determine the MPI latency and bandwidth of the student cluster 
! using the ping pong test.  Do timings only on the rank 0
! processor and ping pong between  rank 0 and rank p-1. 
! Divide timings by 2 to get one-way times.
 
! For latency  use ntrial  = 512, n = 1             and p = 16 and 32
! For bandwidth use ntrial = 16,  n = 64*1024*1024  and p = 16 and 32
 
! Latency times are usually measured in microseconds, time*1.d+6
 
! Bandwidth is usually measured in megabytes/second.  Since we are
! passing messages A(1:n) and B(1:n) in double precision so each
! array element has 8 bytes so A and B have 8*n bytes, thus
! bandwidth = dble(8*n)*1.d-6)/time  MBytes/second.
 
  use mpi 
  implicit none
  !integer, parameter ::  n = 1, ntrial = 512, tag = 1 ! For latency
  integer, parameter ::  n = 64*1024*1024, ntrial = 16, tag = 1 ! For bandwidth
  integer, parameter ::  dp = mpi_double_precision, comm=mpi_comm_world
  double precision   ::  t1, t2, time(ntrial), average, bandwidth
  double precision   ::  A(n), B(n)
  integer            ::  i, ierror, ktrial
  integer            ::  p, rank, status(mpi_status_size)
  call mpi_init(ierror) 
  call mpi_comm_size(comm, p, ierror)
  call mpi_comm_rank(comm, rank, ierror)
      
  if (rank==0) call random_number(A)  ! initialize A only on proc 0
 
  do ktrial = 1, ntrial
     if (rank == 0) then
        t1 = mpi_wtime()
        call mpi_send(A(1), n, dp, p-1, tag, comm, ierror)     
        call mpi_recv(A(1), n, dp, p-1, tag, comm, mpi_status_ignore, ierror)                              
        t2 = mpi_wtime()
     endif
     if (rank == p - 1) then                      
        call mpi_recv(B(1), n, dp, 0, tag, comm, mpi_status_ignore, ierror)
        call mpi_send(B(1), n, dp, 0, tag, comm, ierror)
     endif                                                
     time(ktrial) = (t2 - t1) / 2.d0 * 1.d+6  ! compute microseconds

  enddo ! for the ktrial loop

  if (rank == 0) then
     print*,'n = ',n,'ntrial = ',ntrial,'p = ', p
     average = sum(time(1:ntrial))/dble(ntrial)
     if(n == 1)  print*,'average time = ', average,'microseconds'
!    compute bandwidth in MB/second
     average = average*1.d-6 ! time in seconds
     bandwidth = dble(8*n)*1.d-6/average
    if(n > 1) print*,'bandwidth    = ', bandwidth,'MBytes/second'
   endif

  call mpi_finalize(ierror)
  end

! OUTPUT using salloc -N 2 -n 32 -t 2:00:00

! BANDWIDTH TIMES - results show the time bandwidth for p = 16 and 32 are
! about the same.

! [yhkwon@hpc-class04 fortran]$ mpirun -np 16 ./a.out
! n =     67108864 ntrial =           16 p =           16
! bandwidth    =    3527.92090798810      MBytes/second

! [yhkwon@hpc-class04 fortran]$ mpirun -np 32 ./a.out
! n =     67108864 ntrial =           16 p =           32
! bandwidth    =    3328.93344626928      MBytes/second


! LATENCY TIMES - latencies are much smaller within a node

! [yhkwon@hpc-class04 fortran]$ mpirun -np 16 ./a.out
! n =            1 ntrial =          512 p =           16
! average time =   0.898959115147591      microseconds

! [yhkwon@hpc-class04 fortran]$ mpirun -np 32 ./a.out
! n =            1 ntrial =          512 p =           32
! average time =    2.27359123528004      microseconds

