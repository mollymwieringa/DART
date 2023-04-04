! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module gamma_distribution_mod

use types_mod,               only : r8, PI, missing_r8

use utilities_mod,           only : E_ERR, error_handler

use normal_distribution_mod, only : normal_cdf, inv_cdf

use distribution_params_mod, only : distribution_params_type

use random_seq_mod,          only : random_seq_type, random_uniform

implicit none
private

public :: gamma_cdf,        inv_gamma_cdf,                                  &
          gamma_cdf_params, inv_gamma_cdf_params,                           &
          random_gamma, gamma_pdf, test_gamma, gamma_mn_var_to_shape_scale, &
          gamma_gamma_prod, gamma_shape_scale

character(len=512)          :: errstring
character(len=*), parameter :: source = 'gamma_distribution_mod.f90'

real(r8), parameter :: failed_value = -99.9_r8

contains

!-----------------------------------------------------------------------

subroutine test_gamma

! This routine provides limited tests of the numerics in this module. It begins
! by comparing a handful of cases of the pdf and cdf to results from Matlab. It
! then tests the quality of the inverse cdf for a single shape/scale pair. Failing
! these tests suggests a serious problem. Passing them does not indicate that 
! there are acceptable results for all possible inputs. 

real(r8) :: x, y, inv
real(r8) :: mean, variance, sd, gamma_shape, gamma_scale, max_diff
integer :: i

! Comparative results for a handful of cases from MATLAB21a
real(r8) :: pdf_diff(7), cdf_diff(7)
real(r8) :: mshape(7) = [1.0_r8, 2.0_r8, 3.0_r8, 5.0_r8, 9.0_r8, 7.5_r8, 0.5_r8]
real(r8) :: mscale(7) = [2.0_r8, 2.0_r8, 2.0_r8, 1.0_r8, 0.5_r8, 1.0_r8, 1.0_r8]
real(r8) :: mx(7)     = [1.0_r8, 2.0_r8, 3.0_r8, 4.0_r8, 5.0_r8, 6.0_r8, 7.0_r8]
! Generated by matlab gampdf(mx, mshape, mscale)
real(r8) :: mpdf(7) = [0.303265329856317_r8, 0.183939720585721_r8, 0.125510715083492_r8, &
                       0.195366814813165_r8, 0.225198064298040_r8, 0.151385201555322_r8, &
                       0.000194453010092_r8]
! Generated by matlab gamcdf(mx, mshape, mscale)
real(r8) :: mcdf(7) = [0.393469340287367_r8, 0.264241117657115_r8, 0.191153169461942_r8, &
                       0.371163064820127_r8, 0.667180321249281_r8, 0.320970942909585_r8, &
                       0.999817189367018_r8]

! Compare to matlab
write(*, *) 'Absolute value of differences should be less than 1e-15'
do i = 1, 7
   pdf_diff(i) = gamma_pdf(mx(i), mshape(i), mscale(i)) - mpdf(i)
   cdf_diff(i) = gamma_cdf(mx(i), mshape(i), mscale(i), .true., .false., 0.0_r8, missing_r8) - mcdf(i)
   write(*, *) i, pdf_diff(i), cdf_diff(i)
end do

! Input a mean and variance
mean = 10.0_r8
sd = 1.0_r8
variance = sd**2

! Get shape and scale
gamma_shape = mean**2 / variance
gamma_scale = variance / mean

! Test the inversion of the cdf over +/- 5 standard deviations around mean
max_diff = -1.0_r8
do i = 0, 1000
   x = mean + ((i - 500.0_r8) / 500.0_r8) * 5.0_r8 * sd
   y = gamma_cdf(x, gamma_shape, gamma_scale, .true., .false., 0.0_r8, missing_r8)
   inv = inv_gamma_cdf(y, gamma_shape, gamma_scale, .true., .false., 0.0_r8, missing_r8)
   max_diff = max(abs(x-inv), max_diff)
end do

write(*, *) '----------------------------'
write(*, *) 'max difference in inversion is ', max_diff
write(*, *) 'max difference should be less than 1e-11'

end subroutine test_gamma

!-----------------------------------------------------------------------

function inv_gamma_cdf_params(quantile, p) result(x)

real(r8)                                   :: x
real(r8),                       intent(in) :: quantile
type(distribution_params_type), intent(in) :: p

! Could do error checks for gamma_shape and gamma_scale values here
x = inv_cdf(quantile, gamma_cdf_params, inv_gamma_first_guess_params, p)

end function inv_gamma_cdf_params
!-----------------------------------------------------------------------

function inv_gamma_cdf(quantile, gamma_shape, gamma_scale, &
               bounded_below, bounded_above, lower_bound, upper_bound) result(x)

real(r8)             :: x
real(r8), intent(in) :: quantile
real(r8), intent(in) :: gamma_shape
real(r8), intent(in) :: gamma_scale
logical,  intent(in) :: bounded_below, bounded_above
real(r8), intent(in) :: lower_bound,   upper_bound

! Given a quantile q, finds the value of x for which the gamma cdf
! with shape and scale has approximately this quantile

type(distribution_params_type) :: p

! Load the p type for the generic cdf calls
p%params(1) = gamma_shape; p%params(2) = gamma_scale
p%bounded_below = bounded_below;      p%bounded_above = bounded_above
p%lower_bound   = lower_bound;        p%upper_bound   = upper_bound

x = inv_gamma_cdf_params(quantile, p)

end function inv_gamma_cdf

!---------------------------------------------------------------------------

function gamma_pdf(x, gamma_shape, gamma_scale)

! Returns the probability density of a gamma function with shape and scale
! at the value x

real(r8)             :: gamma_pdf
real(r8), intent(in) :: x, gamma_shape, gamma_scale

! All inputs must be nonnegative
if(x < 0.0_r8 .or. gamma_shape < 0.0_r8 .or. gamma_scale < 0.0_r8) then
   gamma_pdf = failed_value
else
   gamma_pdf = x**(gamma_shape - 1.0_r8) * exp(-x / gamma_scale) / &
      (gamma(gamma_shape) * gamma_scale**gamma_shape)
endif

end function gamma_pdf

!---------------------------------------------------------------------------

function gamma_cdf_params(x, p)

real(r8)                                   :: gamma_cdf_params
real(r8), intent(in)                       :: x
type(distribution_params_type), intent(in) :: p

! A translation routine that is required to use the generic cdf optimization routine
! Extracts the appropriate information from the distribution_params_type that is needed
! for a call to the function gamma_cdf below. 

real(r8) :: gamma_shape, gamma_scale

gamma_shape = p%params(1);     gamma_scale = p%params(2)
gamma_cdf_params = gamma_cdf(x, gamma_shape, gamma_scale, &
                     p%bounded_below, p%bounded_above, p%lower_bound, p%upper_bound)

end function gamma_cdf_params

!---------------------------------------------------------------------------

function gamma_cdf(x, gamma_shape, gamma_scale, bounded_below, bounded_above, lower_bound, upper_bound)

! Returns the cumulative distribution of a gamma function with shape and scale
! at the value x

real(r8) :: gamma_cdf
real(r8), intent(in) :: x, gamma_shape, gamma_scale
logical,  intent(in) :: bounded_below, bounded_above
real(r8), intent(in) :: lower_bound,   upper_bound

! All inputs must be nonnegative
if(x < 0.0_r8 .or. gamma_shape < 0.0_r8 .or. gamma_scale < 0.0_r8) then
   gamma_cdf = failed_value
elseif(x == 0.0_r8) then
   gamma_cdf = 0.0_r8 
else
   ! Use definition as incomplete gamma ratio to gamma
   gamma_cdf = gammad(x / gamma_scale, gamma_shape)
endif

end function gamma_cdf

!---------------------------------------------------------------------------
function gammad (x, p)

implicit none 

real(r8)             :: gammad
real(r8), intent(in) :: x
real(r8), intent(in) :: p

!*****************************************************************************80
!
!! GAMMAD computes the Incomplete Gamma Integral
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by B Shea.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    B Shea,
!    Algorithm AS 239:
!    Chi-squared and Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete 
!    gamma ratio.  0 <= X, and 0 < P.
!
!
!    Output, real ( kind = 8 ) GAMMAD, the value of the incomplete 
!    Gamma integral.
!

real(r8)            :: a, b, c, an, arg, pn(6), rn
real(r8), parameter :: elimit = - 88.0_r8
real(r8), parameter :: oflo   = 1.0e+37_r8
real(r8), parameter :: plimit = 1000.0_r8
real(r8), parameter :: tol    = 1.0e-14_r8
real(r8), parameter :: xbig   = 1.0e+08_r8

! x zero returns zero
if(x == 0.0_r8) then
   gammad = 0.0_r8
elseif(xbig < x) then
   !  If X is large set GAMMAD = 1.
   gammad = 1.0_r8
elseif(plimit < p) then
! If P is large, use a normal approximation.
   pn(1) = 3.0_r8 * sqrt(p) * ((x / p)**(1.0_r8 / 3.0_r8) + &
       1.0_r8 / (9.0_r8 * p) - 1.0_r8)
   gammad = normal_cdf(pn(1), 0.0_r8, 1.0_r8)
elseif(x <= 1.0_r8 .or. x < p) then
!  Use Pearson's series expansion.
!  Original note: (Note that P is not large enough to force overflow in logAM).
   arg = p * log(x) - x - log(gamma(p + 1.0_r8))
   c = 1.0_r8
   gammad = 1.0_r8
   a = p

   do
      a = a + 1.0_r8
      c = c * x / a
      gammad = gammad + c
      if(c <= tol) exit
   end do

   arg = arg + log(gammad)

   if(elimit <= arg) then
     gammad = exp(arg)
   else
      gammad = 0.0_r8
    end if
else 
   !  Use a continued fraction expansion.
   arg = p * log(x) - x - log(gamma(p))
   a = 1.0_r8 - p
   b = a + x + 1.0_r8
   c = 0.0_r8
   pn(1) = 1.0_r8
   pn(2) = x
   pn(3) = x + 1.0_r8
   pn(4) = x * b
   gammad = pn(3) / pn(4)

   do
      a = a + 1.0_r8
      b = b + 2.0_r8
      c = c + 1.0_r8
      an = a * c
      pn(5) = b * pn(3) - an * pn(1)
      pn(6) = b * pn(4) - an * pn(2)

      if (pn(6) /= 0.0_r8) then
         rn = pn(5) / pn(6)
         if(abs(gammad - rn) <= min(tol, tol * rn)) exit
         gammad = rn
      end if

      pn(1) = pn(3)
      pn(2) = pn(4)
      pn(3) = pn(5)
      pn(4) = pn(6)

      ! Re-scale terms in continued fraction if terms are large.
      if (oflo <= abs(pn(5))) pn(1:4) = pn(1:4) / oflo

   end do

   arg = arg + log(gammad)

   if (elimit <= arg) then
      gammad = 1.0_r8 - exp(arg)
   else
      gammad = 1.0_r8
    endif
endif

end function gammad

!---------------------------------------------------------------------------

function random_gamma(r, rshape, rscale)

! Note that this provides same qualitative functionality as a similarly named
! routine in the random_seq_mod that uses a rejection algorithm. However, once
! we have an inverse cdf function for a distribution, it is possible to generate
! random numbers by first getting a draw from a U(0, 1) and then inverting these
! quantiles to get an actual value

type(random_seq_type), intent(inout) :: r
real(r8),              intent(in)    :: rshape
real(r8),              intent(in)    :: rscale
real(r8)                             :: random_gamma

real(r8) :: quantile
if (rshape <= 0.0_r8) then
   write(errstring, *) 'Shape parameter must be positive, was ', rshape
   call error_handler(E_ERR, 'random_gamma', errstring, source)
endif

if (rscale <= 0.0_r8) then
   write(errstring, *) 'Scale parameter (scale=1/rate) must be positive, was ', rscale
   call error_handler(E_ERR, 'random_gamma', errstring, source)
endif

! Draw from U(0, 1) to get a quantile
quantile = random_uniform(r)
! Invert cdf to get a draw from gamma
random_gamma = inv_gamma_cdf(quantile, rshape, rscale, .true., .false., 0.0_r8, missing_r8)

end function random_gamma

!---------------------------------------------------------------------------

subroutine gamma_shape_scale(x, num, gamma_shape, gamma_scale)

integer,  intent(in)  :: num
real(r8), intent(in)  :: x(num)
real(r8), intent(out) :: gamma_shape
real(r8), intent(out) :: gamma_scale

! This subroutine computes a shape and scale from a sample
! It first computes the mean and sd, then converts
! Note that this is NOT the maximum likelihood estimator from the sample 
! and computing that would be an alternative method to get shape and scale

real(r8) :: mean, variance

mean = sum(x) / num
variance  = sum((x - mean)**2) / (num - 1)

call gamma_mn_var_to_shape_scale(mean, variance, gamma_shape, gamma_scale)

end subroutine gamma_shape_scale

!---------------------------------------------------------------------------

subroutine gamma_mn_var_to_shape_scale(mean, variance, gamma_shape, gamma_scale)

real(r8), intent(in)  :: mean, variance
real(r8), intent(out) :: gamma_shape, gamma_scale

gamma_shape = mean**2 / variance
gamma_scale = variance / mean

end subroutine gamma_mn_var_to_shape_scale

!---------------------------------------------------------------------------

subroutine gamma_gamma_prod(prior_shape, prior_scale, like_shape, like_scale, &
   post_shape, post_scale)

real(r8), intent(in)  :: prior_shape, prior_scale, like_shape, like_scale
real(r8), intent(out) :: post_shape, post_scale

! Compute statistics of product of two gammas

post_shape = prior_shape + like_shape - 1
post_scale = prior_scale * like_scale / (prior_scale + like_scale)

end subroutine gamma_gamma_prod

!---------------------------------------------------------------------------

function inv_gamma_first_guess_params(quantile, p)

real(r8)                                   :: inv_gamma_first_guess_params
real(r8), intent(in)                       :: quantile
type(distribution_params_type), intent(in) :: p

! A translation routine that is required to use the generic first_guess for
! the cdf  optimization routine.
! Extracts the appropriate information from the distribution_params_type that is needed
! for a call to the function approx_inv_normal_cdf below (which is nothing).

real(r8) :: gamma_shape, gamma_scale

gamma_shape = p%params(1);     gamma_scale = p%params(2)
inv_gamma_first_guess_params = inv_gamma_first_guess(quantile, gamma_shape, gamma_scale)

end function inv_gamma_first_guess_params

!---------------------------------------------------------------------------

function inv_gamma_first_guess(quantile, gamma_shape, gamma_scale)

real(r8) :: inv_gamma_first_guess
real(r8), intent(in) :: quantile
real(r8), intent(in) :: gamma_shape, gamma_scale

! Need some sort of first guess, should be smarter here
! For starters, take the mean for this shape and scale
inv_gamma_first_guess = gamma_shape * gamma_scale
! Could use info about sd to further refine mean and reduce iterations
!!!sd = sqrt(gamma_shape * gamma_scale**2)

end function inv_gamma_first_guess

!---------------------------------------------------------------------------

end module gamma_distribution_mod
