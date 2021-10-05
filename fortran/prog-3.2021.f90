!   Prog-3: Matrix addition: A = B + C with stride 1 and timing
!   Programmer:  G. Luecke

    use mpi
    implicit none
    integer,parameter ::  n = 1024, ntrial = 128
    double precision  ::  A(n,n), B(n,n), C(n,n)
    integer           ::  i, j, ktrial, ierror
    double precision  ::  t0, t1, time(ntrial)
    double precision  ::  min_time, average_time, max_time, mflops, flops
    call mpi_init(ierror)

!   Initialize Data
    call random_number(B)
    call random_number(C)

!   Calculations
    do ktrial = 1, ntrial
       t0 = mpi_wtime()
       do i = 1, n
          do j = 1, n
 !           A(i,j) = B(i,j) + C(i,j) ! bad stride - change loop order
 !           or
             A(j,i) = B(j,i) + C(j,i) ! stride 1
          enddo
       enddo
!!     best way to write is
!!     do j = 1,n
!!        A(1:n,j) = B(1:n,j) + C(1:n,j)
!!     enddo
       t1 = mpi_wtime()
!  OR  time(ktrial) = mpi_wtime() - t0
       time(ktrial) = t1 - t0 ! seconds
     enddo

    print*,' n = ',n
    print*,'number of timings = ', ntrial

    average_time = sum(time)/dble(ntrial)
    min_time     = minval(time)
    max_time     = maxval(time)
    print*,'minimum time =', min_time
    print*,'average time =', average_time
    print*,'maximum time =', max_time
!   number of floating point operations = n*n
    mflops = dble(n*n)/(average_time*1.d6)
    print*,'mflops = ', mflops

    call mpi_finalize(ierror)
    end

! Output
! class fortran]$ mpiifort   prog-3.2021.f90 
! [grl@hpc-class fortran]$ mpirun -np 1 ./a.out
!  n =         1024
! number of timings =          128
! minimum time =  2.403020858764648E-003
! average time =  2.429928630590439E-003
! maximum time =  3.808975219726562E-003
! mflops =    431.525431158532 

