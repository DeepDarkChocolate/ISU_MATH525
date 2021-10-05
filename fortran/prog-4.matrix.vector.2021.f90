!   Homework 1, due Thursday February 18, 2021
!   Student's Name: Yonghyun Kwon           
!   Matrix Vector: y = y + Ax, compare the performancee of the
!   i-j and j-i versions computing average times and mflops for
!   each version. Use the notation below.
!   Compile as:  mpiifort prog-4.matrix.vector.2021.f90
!   Execute as:  mpirun -np 1 ./a.out

    use mpi
    implicit none
    integer,parameter ::  n = 8*1024, ntrial = 64
    double precision  ::  y(n), A(n,n), x(n)           
    integer           ::  i, j, ktrial, ierror
    double precision  ::  t0, t1, time(ntrial)
    double precision  ::  min_time, average_time, max_time
    double precision  ::  mflops
    call mpi_init(ierror)

!   Initialize Data
    call random_number(y)
    call random_number(A)
    call random_number(x)

!   i-j CALCULATIONS
    do ktrial = 1, ntrial
       t0 = mpi_wtime()
       do i = 1, n
          do j = 1, n
             y(i) = y(i) + A(i, j) * x(j) ! bad stride
          enddo
       enddo                                            
       t1 = mpi_wtime()
       time(ktrial) = t1 - t0 ! seconds 
                                         
     enddo

    print*,' n = ',n
    print*,'number of timings = ', ntrial

    min_time     = minval(time)
    average_time = sum(time)/dble(ntrial)
    max_time     = maxval(time)
    print*,'i-j PERFORMANCE'
    print*,'minimum i-y time =', min_time
    print*,'average i-j time =', average_time
    print*,'maximum i-j time =', max_time
!   number of floating point operations = n * n
    mflops = dble(n * n * 2) / (average_time * 1.d6)                                
    print*,'mflops = ', mflops

!   j-i CALCULATIONS
    do ktrial = 1, ntrial
       t0 = mpi_wtime()
       do j = 1, n
          do i = 1, n
             y(i) = y(i) + A(i, j) * x(j) ! stride 1
          enddo
       enddo
       t1 = mpi_wtime()
       time(ktrial) = t1 - t0 ! seconds                  
                                       
     enddo

    print*,' n = ',n
    print*,'number of timings = ', ntrial

    print*,' '
    print*,'j-i PERFORMANCE'
    min_time     = minval(time)
    average_time = sum(time)/dble(ntrial)
    max_time     = maxval(time)
    print*,'minimum j-i time =', min_time
    print*,'average j-i time =', average_time
    print*,'maximum j-i time =', max_time
    mflops = dble(n * n * 2) / (average_time * 1.d6)                                
    print*,'mflops = ', mflops

    call mpi_finalize(ierror)
    end

! Output
                                                                   

!  n =         8192
!  number of timings =           64
!  i-j PERFORMANCE
!  minimum i-y time =  5.628800392150879E-002
!  average i-j time =  5.642276257276535E-002
!  maximum i-j time =  5.731582641601562E-002
!  mflops =    2378.78689167172
!  n =         8192
!  number of timings =           64

! j-i PERFORMANCE
!  minimum j-i time =  5.628585815429688E-002
!  average j-i time =  5.639070272445679E-002
!  maximum j-i time =  5.720901489257812E-002
!  mflops =    2380.13930515871                                         
