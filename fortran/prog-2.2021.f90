!   Vector addition: z = x + y

!   Compile as:  mpiifort program.f90
!   Execute as:  mpirun -np 1 ./a.out

    use mpi  !  brings in the MPI libraries   
    implicit none
    integer,parameter ::  n = 5
    double precision  ::  x(n), y(n), z(n)
    integer           ::  i, ierror
    call mpi_init(ierror) ! allows one to use MPI routines

!   Initialize Data:
!   Initialize x and y using a random number generator on (0,1)
    call random_number(x) 
    x(1:n) = x(1:n) + 1.d0 ! or x = x + 1.d0
    y(1:n) = x(1:n) + 3.d0 ! or y = y + 3.d0

!   Calculations:
    do i = 1, n
       z(i) = x(i) + y(i)
    enddo
!   OR
!   z(1:n) = x(1:n) + y(1:n)
!   OR
!   z = x + y

!   Print Output
    print*,'z = ', z

    call mpi_finalize(ierror) ! ends MPI communication
    end

!   Copy Output
!   [grl@hpc-class fortran]$ mpiifort  prog-2.2021.f90 
!   [grl@hpc-class fortran]$ mpirun -np 1 ./a.out
!    z =    5.00000078417364        5.05096088551528
!    5.70503232252213     
!      6.33382896304850        6.92611106378931     
