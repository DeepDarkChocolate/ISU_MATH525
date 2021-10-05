subroutine updatebeta(x, y, w1, w2, beta1_t1, sigma2e_t1, mua_t1,&
sigma2a_t1, n, m, &
m1, m2, M1i, M2i, theta_final)
  implicit none
  external :: sgesv
  integer, intent(in) :: n, m, m1, m2
  integer, intent(in) :: M1i(n), M2i(n)
  double precision, intent(in) :: x(m, n), y(m, n)
  double precision :: xi(m), yi(m), w1i, w2i(m), pi2mati(m, m)
  double precision, intent(in) :: w1(n), w2(m, n)
  double precision, intent(in) :: beta1_t1, sigma2e_t1, mua_t1, sigma2a_t1
  double precision :: beta1_t2, sigma2e_t2, mua_t2, sigmaa_t2, sigma2a_t2
  double precision :: etaa_t1
  double precision :: theta_final(4)
  double precision :: sigmaa_t1, sigmae_t1
  double precision :: hi_ftn
  double precision :: theta_t1(4), theta_t2(4)
  double precision :: Ev1(1:n), Vv1(1:n)

  double precision, parameter :: pi = 3.14159265358979323846264338327950D+00
  double precision, parameter :: alpha_qd = 0.0d0, beta_qd = 0.0d0
  integer, parameter :: n_qd = 100, kind_qd = 6
  double precision :: x_qd(n_qd), w_qd(n_qd)

  double precision :: denom, num1, num2, tmpval, Ev1i, VV1i
  integer :: i, j, LDA, LDB, INFO
  double precision :: beta_mat(2, 2), beta_vec(2), beta_t2(2)
  Real (kind = 4) :: beta_mat2(2, 2), beta_vec2(2) ! for Lapack use
  integer :: pivot(2)

  sigmaa_t1 = sqrt(sigma2a_t1)
  etaa_t1 = 0.5d0 / sigmaa_t1

  call cgqf(n_qd, kind_qd, alpha_qd, beta_qd, mua_t1, etaa_t1, x_qd, w_qd)
  w_qd = w_qd / sum(w_qd)
  
  beta_mat = 0.d0
  beta_vec = 0.d0
  do i = 1, n
    xi = x(:, i)
    yi = y(:, i)
    w1i = w1(i)
    w2i = w2(:, i)

    pi2mati(1:m1, 1:m1) = dble(m1 * (m1 - 1)) / dble(M1i(i) * (M1i(i) -1))
    pi2mati(1:m1, m1+1:m) = dble(m1 * m2) / dble(M1i(i) * M2i(i))
    pi2mati(m1+1:m, 1:m1) = dble(m1 * m2) / dble(M1i(i) * M2i(i))
    pi2mati(m1+1:m, m1+1:m) = dble(m2 * (m2 - 1)) / dble(M2i(i) * (M2i(i) - 1))
    do j = 1, m1
      pi2mati(j, j) = dble(m1) / dble(M1i(i))
    enddo
    do j = m1+1, m
      pi2mati(j, j) = dble(m2) / dble(M2i(i))
    enddo
    
    num1 = 0.d0
    num2 = 0.d0
    denom = 0.d0
    do j = 1, n_qd
       !if(i == 13 .and. j == 18) then
       tmpval = hi_ftn(x_qd(j), xi, yi, w2i, pi2mati, &
         beta1_t1, sigma2e_t1, mua_t1, sigmaa_t1, m)
       if(i == 13 .and. j == 18) then
         !print*, 'tmpval = ', tmpval
         !print*, 'xi = ', xi
         !print*, 'yi = ', yi
         !print*, 'w2i = ', w2i
         !print*, 'pi2mati = ', pi2mati
         !print*, 'm = ', m
         !print*, 'x_qd = ', x_qd
       endif
       if(isnan(tmpval)) print*, 'i = ', i, 'j = ', j
       denom = denom + tmpval * w_qd(j)
       num1 = num1 + tmpval * x_qd(j) * w_qd(j)
       num2 = num2 + tmpval * x_qd(j) * x_qd(j) * w_qd(j)
       !endif
    enddo
    Ev1i = num1 / denom
    Vv1i = num2 / denom - Ev1i * Ev1i
    Ev1(i) = Ev1i
    Vv1(i) = Vv1i
    
    beta_mat(1, 1) = beta_mat(1, 1) + w1i * sum(w2i)
    beta_mat(2, 1) = beta_mat(2, 1) + w1i * sum(w2i * xi)
    beta_mat(1, 2) = beta_mat(1, 2) + w1i * sum(w2i * xi)
    beta_mat(2, 2) = beta_mat(2, 2) + w1i * sum(w2i * xi * xi)

    beta_vec(1) = beta_vec(1) + w1i * sum(w2i * (yi - Ev1(i)))
    beta_vec(2) = beta_vec(2) + w1i * sum(w2i * (yi - Ev1(i)) * xi)
  enddo

  beta_mat2 = real(beta_mat)
  beta_vec2 = real(beta_vec)
  !beta_mat2(1, 1) = 2.d0
  !beta_mat2(1, 2) = 3.d0
  !beta_mat2(2, 1) = 1.d0
  !beta_mat2(2, 2) = 1.d0
  !beta_vec2(1) = 5.d0
  !beta_vec2(2) = 6.d0
  !print*, 'beta_mat2 = ', beta_mat2
  !print*, 'beta_vec2 = ', beta_vec2
  call SGESV(2, 1, beta_mat2, 2, pivot, beta_vec2, 2, INFO)
  beta_t2(1) = dble(beta_vec2(1))
  beta_t2(2) = dble(beta_vec2(2))
  !print*, sum(Ev1) / dble(n)
  if(info > 0) print*, 'sgesv warning!, info = ', info 
  !print*, 'beta_t2 = ', beta_t2
  !print*, 'beta_vec2 = ', beta_vec2
  !print*, 'beta_mat2 = ', beta_mat2
  sigma2e_t2 = 0.d0
  mua_t2 = 0.d0
  sigma2a_t2 = 0.d0
  do i = 1, n  
    xi = x(:, i)
    yi = y(:, i)
    w1i = w1(i)
    w2i = w2(:, i)
   
    sigma2e_t2 = sigma2e_t2 + w1i * sum( w2i * (-Ev1(i) + yi - beta_t2(1) -&
        beta_t2(2) * xi) ** 2 + Vv1(i) )
    mua_t2 = mua_t2 + w1i * Ev1(i)
    sigma2a_t2 = sigma2a_t2 + w1i * (Ev1(i) ** 2 + Vv1(i))
  enddo
  sigma2e_t2 = sigma2e_t2 / beta_mat(1, 1) 
  mua_t2 = mua_t2 / sum(w1)
  sigma2a_t2 = sigma2a_t2 / sum(w1) - mua_t2 ** 2 

  !print*, 'beta_t2 = ', beta_t2
  !print*, 'sigma2e_t2 = ', sigma2e_t2
  !print*, 'mua_t2 = ', mua_t2
  !print*, 'sigma2a_t2 = ', sigma2a_t2
  theta_final(1) = beta_t2(2)
  theta_final(2) = sigma2e_t2
  theta_final(3) = mua_t2 + beta_t2(1)
  theta_final(4) = sigma2a_t2
  return 
end

function hi_ftn(ai, xi, yi, w2i, pi2mati, &
beta1_t1, sigma2e_t1, beta0_t1, sigmaa_t1, m)
  double precision :: hi_ftn
  integer, intent(in) :: m
  double precision, intent(in) :: ai, xi(m), yi(m), w2i(m), pi2mati(m, m)
  double precision, intent(in) :: beta1_t1, sigma2e_t1, beta0_t1, sigmaa_t1
  double precision :: meanyi, meanxi, pdf1, pdf2
  double precision :: ratioi(m), pi2i(m), pi2mati_tmp(m, m)
  double precision :: pi2i_1m(1, m), pi2i_m1(m, 1), &
                      ratioi_1m(1, m), ratioi_m1(m, 1), vi_tmp(1, 1)
  double precision :: vi, Si, ai_zero, res
  ai_zero = 0.d0
  meanyi = sum(yi) / dble(m)
  meanxi = sum(xi) / dble(m)
  ratioi = (yi(1:m) - meanyi - beta1_t1 * (xi(1:m) - meanxi)) * w2i
  ratioi_1m(1, :) = ratioi
  ratioi_m1(:, 1) = ratioi
  
  pi2i(1:m) = 1.d0 / w2i
  pi2i_1m(1, :) = pi2i
  pi2i_m1(:, 1) = pi2i

  pi2mati_tmp = matmul(pi2i_m1, pi2i_1m) / pi2mati
  vi_tmp = matmul(matmul(ratioi_1m, 1.d0 - pi2mati_tmp), ratioi_m1)
  vi = vi_tmp(1, 1)
  !if(vi < 0) print*, 'vi < 0 warning! vi = ', vi

  Si = sum(w2i * (yi - (beta0_t1 + beta1_t1 * xi + ai)))
  !print*, 'sqrt(vi) = ', sqrt(vi)
  !print*, 'Si = ', Si
  !print*, 'r8_normal_pdf(ai_zero, sigmaa_t1, ai)', r8_normal_pdf(ai_zero, sigmaa_t1, ai)
  !print*, 'r8_normal_pdf(ai_zero, sqrt(vi), Si)', r8_normal_pdf(ai_zero, sqrt(vi), Si)
  pdf1 = r8_normal_pdf(ai_zero, sigmaa_t1, ai)
  pdf2 = r8_normal_pdf(ai_zero, sqrt(vi), Si)
  !hi_ftn = 1.d0
  if(pdf1 <= 0 .or. isnan(pdf1)) then
    pdf1 = 1.d-200
    !hi_ftn = 0.d0
  endif
  if(pdf2 <= 0 .or. isnan(pdf2)) then
    pdf2 = 1.d-200
    !hi_ftn = 0.d0
  endif
  !if(hi_ftn > 0.d0) 
  hi_ftn =  pdf1 * pdf2
  
  return
end
