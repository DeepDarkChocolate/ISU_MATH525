!   Compute n facorial
!   Programmer:  G. Luecke
!   Compile as:  mpiifort prog.f90
!   Execute as:  mpirun -np 1 ./a.out

    use mpi
    implicit none
    integer,parameter ::  n = 21
    integer  ::  fac, i, ierror
    call mpi_init(ierror)
    
    fac = 1
    do i = 1, n
       fac = fac*i
       print*,i,'factorial =', fac
    enddo
    
    call mpi_finalize(ierror)
    end

! Output
! [grl@hpc-class fortran]$ mpiifort  prog-5.2021.factorial.f90 
! [grl@hpc-class fortran]$ ./a.out
!            1 factorial =           1
!            2 factorial =           2
!            3 factorial =           6
!            4 factorial =          24
!            5 factorial =         120
!            6 factorial =         720
!            7 factorial =        5040
!            8 factorial =       40320
!            9 factorial =      362880
!           10 factorial =     3628800
!           11 factorial =    39916800
!           12 factorial =   479001600
!           13 factorial =  1932053504
!           14 factorial =  1278945280
!           15 factorial =  2004310016
!           16 factorial =  2004189184
!           17 factorial =  -288522240
!           18 factorial =  -898433024
!           19 factorial =   109641728
!           20 factorial = -2102132736
!           21 factorial = -1195114496
