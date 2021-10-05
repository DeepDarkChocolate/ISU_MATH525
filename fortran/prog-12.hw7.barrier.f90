!  Hw#7, spring 2021: barrier  Due: April 1, 2021

!  Name: Yonghyun Kwon

! For this homework, you are to write the barrier operation in several
! different ways and compare the performance of your barriers with
! the performance of mpi_barrier for p = 16, 32 and 64 MPI processes 
! using 1, 2 and 4 nodes.
! To make it easy to compare timings, only compute averages excluding
! the first timing. Use the following as a template for timing.    
 
!!! 1. time the send-recv barrier
!!     call mpi_barrier(comm, ierror)
!!  do ktrial = 1, ntrial !  = 512
!!     call mpi_barrier(comm,ierror)
!!     t0 = mpi_wtime()
!!!    ******************************
!!     call barrier1(p, rank, comm)
!!!    ******************************
!!     t1 = mpi_wtime()
!!     call mpi_barrier(comm, ierror)
!!     time(ktrial) = (t1-t0)*1.d6  !  time in microseconds
!!  enddo
!!
!!  call mpi_reduce(time, max_time, ntrial, dp, mpi_max, 0, comm, ierror)
!!
!!  if (rank == 0) then
!!     print *,'time to call       mpi_barrier is '                         
!!     print*,'average time = ', sum(time(1:ntrial))/dble(ntrial)           
!!  end if
 
 
!!***************************************************************
 
!   Use the following template for each of the barriers you write:
 
!!subroutine barrier1(p, rank, comm) ! using only mpi sends and recvs
!!                                   ! assuming they are not synchronos.
!!! The following program deadlocks when replacing mpi_send with mpi_ssend
!!  use mpi
!!  implicit none
!!  integer :: comm, p, rank, ierror, flag, i
!!  . . .
!!  return
!!  end subroutine barrier1
!!**************************************************************
 
 
!  Note: for each of the barriers that you write have the arguments
!  p, rank and comm, e.g. call barrier1(p,rank,comm) to make this
!  information available to the subroutine.
 
! Write the following barriers
!    barrier1. use only mpi_send's and mpi_recv's assuming mpi_send is not
!              synchronous.  Do not use a central manager.
!    barrier2. use mpi_alltoall.
!    barrier3. use only mpi_bcast.
!    barrier4. use only mpi_gather and mpi_scatter.
!    barrier5. Use only mpi_isend's and mpi_irecv's and do not use a central manager
!    barrier6. use the rank 0 processor as a central manager and use mpi_send     
!              mpi_recv and mpi_bcast.
 
  use mpi
  implicit none
  integer,parameter :: ntrial = 512
  integer :: i,ktrial
  integer, parameter :: dp = mpi_double_precision, comm=mpi_comm_world
  integer            :: ierror, p, rank, status(mpi_status_size)
  double precision   :: t, t0, t1, time(1:ntrial), max_time(1:ntrial)
  call mpi_init(ierror)
  call mpi_comm_size(comm, p, ierror)
  call mpi_comm_rank(comm, rank, ierror)
 
  if (rank == 0) then
      print*,'Times in microseconds, p = ', p
  endif 
 
! time mpi_barrier
  do ktrial = 1, ntrial
     call mpi_barrier(comm,ierror)
     t0 = mpi_wtime()
!    ******************************
     call mpi_barrier(comm,ierror)
!    ******************************
     t1 = mpi_wtime()
     call mpi_barrier(comm, ierror)
     time(ktrial) = (t1-t0)*1.d6  !  time in microseconds
  enddo
 
  call mpi_reduce(time, max_time, ntrial, dp, mpi_max, 0, comm, ierror)
 
  if (rank == 0) then
     print *,'time to call       mpi_barrier is '
     print*,'average time = ', sum(max_time(1:ntrial))/dble(ntrial)
  end if
 
! 1. time the send-recv barrier 
  do ktrial = 1, ntrial
     call mpi_barrier(comm,ierror)
     t0 = mpi_wtime()
!    ******************************
     call barrier1(p, rank, comm)
!    ******************************
     t1 = mpi_wtime()
     call mpi_barrier(comm, ierror)
     time(ktrial) = (t1-t0)*1.d6  !  time in microseconds
  enddo
 
  call mpi_reduce(time, max_time, ntrial, dp, mpi_max, 0, comm, ierror)
 
  if (rank == 0) then
     print *,'time to call the mpi_ send/mpi-recv barrier is '
     print*,'average time = ', sum(max_time(1:ntrial))/dble(ntrial)
  end if
 
 
! 2. time the alltoall barrier
  do ktrial = 1, ntrial
     call mpi_barrier(comm,ierror)
     t0 = mpi_wtime()
!    ******************************
     call barrier2(p, rank, comm)
!    ******************************
     t1 = mpi_wtime()
     call mpi_barrier(comm, ierror)
     time(ktrial) = (t1-t0)*1.d6  !  time in microseconds
  enddo
 
  call mpi_reduce(time, max_time, ntrial, dp, mpi_max, 0, comm, ierror)
 
  if (rank == 0) then
     print *,'time to call the mpi_alltoall barrier is '
     print*,'average time = ', sum(max_time(1:ntrial))/dble(ntrial)
  end if
 
 
! 3. time the bcast barrier 
  do ktrial = 1, ntrial
     call mpi_barrier(comm,ierror)
     t0 = mpi_wtime()
!    ******************************
     call barrier3(p, rank, comm)
!    ******************************
     t1 = mpi_wtime()
     call mpi_barrier(comm, ierror)
     time(ktrial) = (t1-t0)*1.d6  !  time in microseconds
  enddo
 
  call mpi_reduce(time, max_time, ntrial, dp, mpi_max, 0, comm, ierror)
 
  if (rank == 0) then
     print *,'time to call the mpi_bcast barrier is '
     print*,'average time = ', sum(max_time(1:ntrial))/dble(ntrial)
  end if
 
 
! 4. time the gather/scatter barrier
  do ktrial = 1, ntrial
     call mpi_barrier(comm,ierror)
     t0 = mpi_wtime()
!    ******************************
     call barrier4(p, rank, comm)
!    ******************************
     t1 = mpi_wtime()
     call mpi_barrier(comm, ierror)
     time(ktrial) = (t1-t0)*1.d6  !  time in microseconds
  enddo
 
  call mpi_reduce(time, max_time, ntrial, dp, mpi_max, 0, comm, ierror)
 
  if (rank == 0) then
     print *,'time to call the mpi_gath/mpi_scatter barrier is'
     print*,'average time = ', sum(max_time(1:ntrial))/dble(ntrial)
  end if
 
 
! 5. time the isend - irecv barrier
  do ktrial = 1, ntrial
     call mpi_barrier(comm,ierror)
     t0 = mpi_wtime()
!    ******************************
     call barrier5(p, rank, comm)
!    ******************************
     t1 = mpi_wtime()
     call mpi_barrier(comm, ierror)
     time(ktrial) = (t1-t0)*1.d6  !  time in microseconds
  enddo
 
  call mpi_reduce(time, max_time, ntrial, dp, mpi_max, 0, comm, ierror)
 
  if (rank == 0) then
     print *,'time to call the mpi_isend/mpi_irecv barrier is'
     print*,'average time = ', sum(max_time(1:ntrial))/dble(ntrial)
  end if
 
 
! 6. time a central manager (rank=0) barrier
  do ktrial = 1, ntrial
     call mpi_barrier(comm,ierror)
     t0 = mpi_wtime()
!    ******************************
     call barrier6(p, rank, comm)
!    ******************************
     t1 = mpi_wtime()
     call mpi_barrier(comm, ierror)
     time(ktrial) = (t1-t0)*1.d6  !  time in microseconds
  enddo
 
  call mpi_reduce(time, max_time, ntrial, dp, mpi_max, 0, comm, ierror)
 
  if (rank == 0) then
     print *,'time to call c-manager barrier using mpi_send and mpi_bcast is'
     print*,'average time = ', sum(max_time(1:ntrial))/dble(ntrial)
  end if
 
 
  call mpi_finalize(ierror)
  end
 
subroutine barrier1(p, rank, comm) ! using mpi sends and recvs assuming
                                   ! they are not synchonous.
! The following program deadlocks when replacing mpi_send with mpi_ssend
  use mpi
  implicit none
  integer :: comm, p, rank, ierror, i !, . . .
  integer, parameter :: tag = 1
  integer, parameter :: dp = mpi_double_precision
  double precision :: A(1)           
  do i = 0, p - 1  
     if (rank < i) then              
       call mpi_send(A(1), 1, dp, i, tag, comm, ierror)                                 
     endif
     if (rank > i) then
       call mpi_recv(A(1), 1, dp, mpi_any_source, tag, comm, mpi_status_ignore, ierror)     
     endif
  enddo
  return
  end subroutine barrier1
 
 
subroutine barrier2(p, rank, comm) ! using mpi_alltoall
  use mpi
  implicit none
  integer :: p, rank, comm, ierror !, . . .                         
  double precision, allocatable :: sendbuf(:), revcbuf(:)
  integer, parameter :: dp = mpi_double_precision
  allocate(sendbuf(0:p-1), revcbuf(0:p-1))
  call mpi_alltoall(sendbuf(0), 1, dp, revcbuf(0), 1, dp, comm, ierror)
  deallocate(sendbuf, revcbuf)
  return
  end subroutine barrier2
 
 
subroutine barrier3(p, rank, comm) ! using mpi_bcast                 
  use mpi
  implicit none
  integer :: comm, p, rank, ierror, i !, . . .
  integer, parameter :: dp = mpi_double_precision
  double precision :: A(1)   
  do i = 0, p - 1 
    call mpi_bcast(A(1), 1, dp, i, comm, ierror)
  enddo              
  return
  end subroutine barrier3
 
 
subroutine barrier4(p, rank, comm) ! gather/scatter     
  use mpi
  implicit none
  integer :: comm, p, rank, ierror !, . . .     
  double precision :: A(1)
  double precision ,allocatable :: B(:)    
  integer, parameter :: dp = mpi_double_precision
  if(rank == 0) allocate(B(0:p-1))
  call mpi_gather(A(1), 1, dp, B(0), 1, dp, 0, comm, ierror)
  call mpi_scatter(B(0), 1, dp, A(1), 1, dp, 0, comm, ierror)
  if(rank == 0) deallocate(B)                                        
  return
  end subroutine barrier4
 
 
subroutine barrier5(p, rank, comm) ! using isends and irecvs
  use mpi
  implicit none
  integer :: comm, p, rank, ierror, i, count! , . . .    
  integer :: req_arr(1:2*(p-1))
  integer, parameter :: dp = mpi_double_precision
  double precision :: A(1), B(1)  
  count = 0
  do i = 0, p - 1
    if(rank /=  i) then
      count = count + 1
      call mpi_isend(A(1), 1, dp, i, 0, comm, req_arr(2 * count - 1), ierror)
      call mpi_irecv(B(1), 1, dp, i, 0, comm, req_arr(2 * count), ierror)
      call mpi_wait(req_arr(2 * count - 1), mpi_status_ignore, ierror)
      call mpi_wait(req_arr(2 * count), mpi_status_ignore, ierror)
    endif
    !call mpi_waitall(2*(p-1), req_arr, array_of_statuses, ierror)
  enddo
  return
  end subroutine barrier5
 
 
subroutine barrier6(p, rank, comm) ! using rank 0 processor as a central m¡anager
                                   ! and use mpi_send, mpi_recv and mpi_bcast.
  use mpi
  implicit none
  integer :: comm, p, rank, ierror, flag, i
  integer, parameter :: tag = 1
  integer, parameter :: dp = mpi_double_precision
  double precision :: A(1) 
  flag = 0
  if(rank > flag) then                    
    call mpi_send(A(1), 1, dp, 0, tag, comm, ierror)
  endif                                                       
  if(rank == flag) then                            
    do i = 1, p-1
      call mpi_recv(A(1), 1, dp, i, tag, comm, mpi_status_ignore, ierror)                          
    enddo
  endif
  call mpi_bcast(A(1), 1, dp, 0, comm, ierror)
  return
  end subroutine barrier6
 
 
! OUTPUT
![yhkwon@hpc-class67 fortran]$ mpirun -np 16 ./a.out
!export TMPDIR=/scratch/yhkwon/27081
! Times in microseconds, p =           16
! time to call       mpi_barrier is
! average time =    3.56184318661690
! time to call the mpi_ send/mpi-recv barrier is 
! average time =    13.7984752655029
! time to call the mpi_alltoall barrier is 
! average time =    10.6999650597572
! time to call the mpi_bcast barrier is 
! average time =    31.0284085571766
! time to call the mpi_gath/mpi_scatter barrier is
! average time =    11.8883326649666
! time to call the mpi_isend/mpi_irecv barrier is
! average time =    23.8916836678982
! time to call c-manager barrier using mpi_send and mpi_bcast is
! average time =    13.4538859128952
![yhkwon@hpc-class67 fortran]$ mpirun -np 32 ./a.out
!export TMPDIR=/scratch/yhkwon/27081
! Times in microseconds, p =           32
! time to call       mpi_barrier is 
! average time =    45.3409738838673
! time to call the mpi_ send/mpi-recv barrier is 
! average time =    90.3378240764141
! time to call the mpi_alltoall barrier is 
! average time =    66.4447434246540
! time to call the mpi_bcast barrier is 
! average time =    254.192389547825
! time to call the mpi_gath/mpi_scatter barrier is
! average time =    45.2063977718353
! time to call the mpi_isend/mpi_irecv barrier is
! average time =    252.932775765657
! time to call c-manager barrier using mpi_send and mpi_bcast is
! average time =    115.047208964825
![yhkwon@hpc-class67 fortran]$ mpirun -np 64 ./a.out
!export TMPDIR=/scratch/yhkwon/27081
! Times in microseconds, p =           64
! time to call       mpi_barrier is 
! average time =    135.956797748804
! time to call the mpi_ send/mpi-recv barrier is 
! average time =    592.291355133057
! time to call the mpi_alltoall barrier is 
! average time =    221.609603613615
! time to call the mpi_bcast barrier is 
! average time =    1223.60978275537
! time to call the mpi_gath/mpi_scatter barrier is
! average time =    116.388313472271
! time to call the mpi_isend/mpi_irecv barrier is
! average time =    1923.69287833571
! time to call c-manager barrier using mpi_send and mpi_bcast is
! average time =    569.704454392195
