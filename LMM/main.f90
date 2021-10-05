! PPS + Stratified Sampling under CJS setup
! Name: Yonghyun Kwon

MODULE Statbase
CONTAINS
subroutine sample(Ni, n, weights, res)
  implicit none
  integer, intent(in) :: Ni, n
  double precision, intent(in) :: weights(Ni)
  integer, intent(out) :: res(n)
  integer :: i, j, k, tmpval
  double precision :: u1, sumres, sumres2, r8_uniform_01_sample

  sumres = 0.d0
  do i = 1, Ni
    sumres = sumres + weights(i)
  enddo

  i = 1
  do while(i <= n)
    tmpval = 0
    sumres2 = 0.d0
    u1 = r8_uniform_01_sample()
    do j = 1, Ni
       sumres2 = sumres2 + weights(j) / sumres
       if(sumres2 > u1) exit
    enddo
    do k = 1, n
       if(j == res(k)) tmpval = 1
    enddo
    if(tmpval == 0) then
      res(i) = j
      !print*, 'j = ', j
      i = i + 1
    endif
   enddo
end

subroutine sample_SRS(Ni, n, Idx, res)
  implicit none
  integer, intent(in) :: Ni, n, Idx(Ni)
  integer, intent(out) :: res(n)
  integer :: i, j, k, tmpval
  double precision :: u1, sumres2, r8_uniform_01_sample, Nifrac

  i = 1
  Nifrac = 1.d0 / dble(Ni)
  do while(i <= n)
    tmpval = 0
    sumres2 = 0.d0
    u1 = r8_uniform_01_sample()
    do j = 1, Ni
       sumres2 = sumres2 + Nifrac
       if(sumres2 > u1) exit
    enddo
    do k = 1, n
       if(Idx(j) == res(k)) tmpval = 1
    enddo
    if(tmpval == 0) then
      res(i) = Idx(j)
      !print*, 'j = ', j
      !print*, 'Idx(j) = ', Idx(j)
      i = i + 1
    endif
   enddo
end

END MODULE Statbase

use mpi
use Statbase
implicit none
integer, parameter :: n = 200, m = 12, Ni = 5000, Mi = 500
integer, parameter :: B = 500, K = 5000

double precision, parameter :: beta0 = 0.d0, beta1 = 1.d0, &
sigma2a = 1.d0, sigma2e = 1.d0, mua = 0.5d0, w_h = 1.d0 / 3.d0, &
eps_theta = 1.d-3

double precision, parameter :: pi = 3.14159265358979323846264338327950D+00

!double precision :: beta1_ini, mua_ini, sigma2a_ini, sigma2e_ini
double precision :: beta1_t1, mua_t1, sigma2a_t1, sigma2e_t1, etaa_t1, etae_t1

double precision :: theta(4), theta_res(B, 4), theta_final(4)

double precision :: a(Ni), M_vec(Ni), x(Mi, Ni), e(Mi, Ni), pi1(Ni), y(Mi, Ni)
double precision :: pi1_sampled(n), pi2_sampled(m, n)
double precision :: w1_sampled(n), w2_sampled(m, n)
integer :: SampleIdx(n), SampleIdx2(m, n), M1i(n), M2i(n)
double precision:: x_sampled(m, n), y_sampled(m, n), e_sampled(m, n)

double precision :: r8_normal_sample, r8_normal_pdf
double precision :: sigmaa, sigmae, a_tilde_val, a_val, M_val_sum
integer :: i, j, B_p, M_val, m1, m2, simnum
integer :: idx1, idx2, cnt

integer, allocatable :: e_pos(:), e_neg(:)
integer :: e_pos_num, e_neg_num, k1, k2

double precision, parameter :: alpha_qd = 0.0d0, beta_qd = 0.0d0
integer, parameter :: n_qd = 20, kind_qd = 6
double precision :: x_qd(n_qd), w_qd(n_qd)

integer, parameter :: dp = mpi_double_precision, comm = mpi_comm_world
integer :: ierror, p, rank

double precision :: t1, t2, local_time, local_max_time, max_time
double precision :: local_sum_time, sum_time
double precision, allocatable :: local_time_array(:)

sigmaa = sqrt(sigma2a)
sigmae = sqrt(sigma2e)

call mpi_init(ierror)
call mpi_comm_size(comm, p, ierror)
call mpi_comm_rank(comm, rank, ierror)

call advance_state(rank)
B_p = B / p
if(rank < (B - B_p * p)) B_p = B_p + 1
allocate(local_time_array(1:B_p))
call mpi_barrier(comm, ierror)
local_sum_time = 0.d0

do simnum = 1, B_p
!print*, simnum
t1 = mpi_wtime()

do i = 1, Ni
  a_val = r8_normal_sample(mua, sigmaa)
  a(i) = a_val
  a_val = exp(2.5d0 + a_val)
  a_tilde_val = a_val / (1.d0 + a_val)
  M_val = NINT(dble(Mi) * a_tilde_val)
  M_vec(i) = M_val
  do j = 1, M_val
    x(j, i) = r8_normal_sample(0.d0, 1.d0)
    e(j, i) = r8_normal_sample(0.d0, sigmae)
    y(j, i) = x(j, i) * beta1 + a(i) + e(j, i)
  enddo
enddo

M_val_sum = sum(M_vec)
do i = 1, Ni
  pi1(i) = n * M_vec(i) / M_val_sum
enddo

call sample(Ni, n, M_vec, SampleIdx)

do i = 1, n
  pi1_sampled(i) = pi1(SampleIdx(i))
enddo

m1 = NINT(dble(m) * w_h)
m2 = NINT(dble(m) - dble(m) * w_h)

do i = 1, n
  e_pos_num = 0
  do j = 1, M_vec(sampleIdx(i))
    if (e(j, SampleIdx(i)) > 0.d0) e_pos_num = e_pos_num + 1
  enddo
  e_neg_num = M_vec(sampleIdx(i)) - e_pos_num
  allocate(e_pos(1:e_pos_num), e_neg(1:e_neg_num))
  k1 = 0
  k2 = 0
  do j = 1, M_vec(sampleIdx(i))
    if (e(j, SampleIdx(i)) > 0.d0) then
      k1 = k1 + 1
      e_pos(k1) = j
    endif
    if (e(j, SampleIdx(i)) <= 0.d0) then
      k2 = k2 + 1
      e_neg(k2) = j
    endif
  enddo
  M1i(i) = e_pos_num
  M2i(i) = e_neg_num
  call sample_SRS(e_pos_num, m1, e_pos, SampleIdx2(1:m1, i))
  call sample_SRS(e_neg_num, m2, e_neg, SampleIdx2((m1+1):m, i))
  deallocate(e_pos, e_neg)
enddo

do i = 1, n
  idx1 = SampleIdx(i)
  do j = 1, m
    idx2 = SampleIdx2(j, i)
    x_sampled(j, i) = x(idx2, idx1)
    y_sampled(j, i) = y(idx2, idx2)
    e_sampled(j, i) = e(idx2, idx1)
    if(e_sampled(j, i) > 0.d0) pi2_sampled(j, i) = dble(m1) / dble(M1i(i))
    if(e_sampled(j, i) <= 0.d0) pi2_sampled(j, i) = dble(m2) / dble(M2i(i))
  enddo
enddo

w1_sampled = 1.d0 / pi1_sampled
w2_sampled = 1.d0 / pi2_sampled

beta1_t1 = beta1
mua_t1 = mua
sigma2a_t1 = sigma2a
sigma2e_t1 = sigma2e
etaa_t1 = 0.5d0 / sigma2a

call cgqf(n_qd, kind_qd, alpha_qd, beta_qd, mua_t1, etaa_t1, x_qd, w_qd)
!w_qd(1:n) = w_qd(1:n) * sqrt(etaa_t1 / 2.d0) / sqrt(pi)
w_qd = w_qd / sum(w_qd)

call updatebeta(x_sampled, y_sampled, w1_sampled, w2_sampled, &
beta1_t1, sigma2e_t1, mua_t1, sigma2a_t1, n, m, &
m1, m2, M1i, M2i, theta_final)

t2 = mpi_wtime()
local_time = t2 - t1
local_time_array(simnum) = local_time
local_sum_time = local_sum_time + local_time
enddo
call mpi_barrier(comm, ierror)

local_max_time = maxval(local_time_array)

call mpi_reduce(local_max_time, max_time, 1, dp, mpi_max, 0, comm, ierror)
call mpi_reduce(local_sum_time, sum_time, 1, dp, mpi_max, 0, comm, ierror)

!print*, sum(w_qd)

!Test sample ftn
!do i = 1, Ni
!  M_vec(i) = 0.d0
!enddo
!M_vec(1) = 100.d0
!M_vec(2) = 1.d0
!M_vec(3) = 10.d0
!call sample(Ni, 2, M_vec, SampleIdx)

!call sample_SRS(n, m, SampleIdx, SampleIdx2(1:m, 1))


!do i = 1, n_p
!  a_val = r8_normal_sample(mua, sigmaa)
!  a(i) = a_val
!  a_val = exp(2.5d0 + a_val)
!  a_tilde_val = a_val / (1 + a_val)
!  M_val = NINT(dble(Mi) * a_tilde_val)
  !allocate(x(n_p, M_val), y(n_p, M_val))
!enddo

if(rank == 0) then
   print*, 'local_max_time = ', local_max_time, 'seconds'
   print*, 'max_time = ', max_time, 'seconds'
   print*, 'local_sum_time = ', local_sum_time / dble(B_p), 'seconds'
   print*, 'sum_time = ', sum_time, 'seconds'
!  print*, a(1:3)
!  print*, sampleIdx(1:2)
endif

!print*, r8_normal_pdf(sigma2a - 1.d0, sqrt(sigma2a + 1.d0), sigma2a - 1.d0)
!print*, r8_normal_sample(aa, bb)
!print*, x(1:3)



call mpi_finalize(ierror)
end

