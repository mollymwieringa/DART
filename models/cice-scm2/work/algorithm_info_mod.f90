! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module algorithm_info_mod

use types_mod, only : r8, i8

use obs_def_mod, only : obs_def_type, get_obs_def_type_of_obs, get_obs_def_error_variance
use obs_kind_mod, only : get_quantity_for_type_of_obs

! Get the QTY definitions that are needed (aka kind)
use obs_kind_mod, only : QTY_SEAICE_VOLUME, QTY_SEAICE_CONCENTR, QTY_SEAICE_SNOWVOLUME, &
      QTY_SEAICE_AGREG_THICKNESS, QTY_SEAICE_AGREG_CONCENTR, QTY_SEAICE_AGREG_FREEBOARD, &
      QTY_SEAICE_VICE01, QTY_SEAICE_VICE02, QTY_SEAICE_VICE03, QTY_SEAICE_VICE04, QTY_SEAICE_VICE05, &
      QTY_SEAICE_VSNO01, QTY_SEAICE_VSNO02, QTY_SEAICE_VSNO03, QTY_SEAICE_VSNO04, QTY_SEAICE_VSNO05, &
      QTY_SEAICE_AICE01, QTY_SEAICE_AICE02, QTY_SEAICE_AICE03, QTY_SEAICE_AICE04, QTY_SEAICE_AICE05

! NOTE: Sadly, the QTY itself is not sufficient for the POWER because there is additional metadata

use assim_model_mod, only : get_state_meta_data
use location_mod, only    : location_type

implicit none
private

! Defining parameter strings for different observation space filters
! For now, retaining backwards compatibility in assim_tools_mod requires using
! these specific integer values and there is no point in using these in assim_tools.
! That will change if backwards compatibility is removed in the future.
integer, parameter :: EAKF               = 1
integer, parameter :: ENKF               = 2
integer, parameter :: UNBOUNDED_RHF      = 8
integer, parameter :: GAMMA_FILTER       = 11
integer, parameter :: BOUNDED_NORMAL_RHF = 101 

! Defining parameter strings for different prior distributions that can be used for probit transform
integer, parameter :: NORMAL_PRIOR            = 1
integer, parameter :: BOUNDED_NORMAL_RH_PRIOR = 2
integer, parameter :: GAMMA_PRIOR             = 3
integer, parameter :: BETA_PRIOR              = 4
integer, parameter :: LOG_NORMAL_PRIOR        = 5
integer, parameter :: UNIFORM_PRIOR           = 6

public :: obs_error_info, probit_dist_info, obs_inc_info, &
          EAKF, ENKF, BOUNDED_NORMAL_RHF, UNBOUNDED_RHF, GAMMA_FILTER, &
          NORMAL_PRIOR, BOUNDED_NORMAL_RH_PRIOR, GAMMA_PRIOR, BETA_PRIOR, LOG_NORMAL_PRIOR, &
          UNIFORM_PRIOR

! Provides routines that give information about details of algorithms for 
! observation error sampling, observation increments, and the transformations
! for regression and inflation in probit space. 
! For now, it is convenient to have these in a single module since several
! users will be developing their own problem specific versions of these
! subroutines. This will avoid constant merge conflicts as other parts of the
! assimilation code are updated.

contains

!-------------------------------------------------------------------------
subroutine obs_error_info(obs_def, error_variance, bounded, bounds)

! Computes information needed to compute error sample for this observation
! This is called by perfect_model_obs when generating noisy obs
type(obs_def_type), intent(in)  :: obs_def
real(r8),           intent(out) :: error_variance
logical,            intent(out) :: bounded(2)
real(r8),           intent(out) :: bounds(2)

integer     :: obs_type, obs_kind
integer(i8) :: state_var_index
type(location_type) :: temp_loc

! Get the kind of the observation
obs_type = get_obs_def_type_of_obs(obs_def)
! If it is negative, it is an identity obs
if(obs_type < 0) then
  state_var_index = -1 * obs_type
  call get_state_meta_data(state_var_index, temp_loc, obs_kind)
else
  obs_kind = get_quantity_for_type_of_obs(obs_type)
endif

! Get the default error variance
error_variance = get_obs_def_error_variance(obs_def)

! Set the observation error details for each type of quantity
SELECT CASE (obs_kind)
    CASE (QTY_SEAICE_AGREG_CONCENTR, &
         QTY_SEAICE_AICE01    , &
         QTY_SEAICE_AICE02    , &
         QTY_SEAICE_AICE03    , &
         QTY_SEAICE_AICE04    , &
         QTY_SEAICE_AICE05    )
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 1.0_r8
    CASE (QTY_SEAICE_VICE01)
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 0.64_r8
    CASE (QTY_SEAICE_VICE02)
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 1.39_r8
    CASE (QTY_SEAICE_VICE03)
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 2.47_r8
    CASE (QTY_SEAICE_VICE04)
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 4.57_r8
     CASE (QTY_SEAICE_VICE05)
      bounded(1) = .true.
      bounded(2) = .false.
      bounds(1) = 0.0_r8
   CASE (QTY_SEAICE_AGREG_THICKNESS, QTY_SEAICE_AGREG_FREEBOARD, &
         QTY_SEAICE_VSNO01    , &
         QTY_SEAICE_VSNO02    , &
         QTY_SEAICE_VSNO03    , &
         QTY_SEAICE_VSNO04    , &
         QTY_SEAICE_VSNO05    )
      bounded(1) = .true.
      bounded(2) = .false.
      bounds(1) = 0.0_r8
   CASE Default
      bounded(:) = .false.
      return
END SELECT

end subroutine obs_error_info


!-------------------------------------------------------------------------


subroutine probit_dist_info(kind, is_state, is_inflation, dist_type, &
   bounded, bounds)

! Computes the details of the probit transform for initial experiments
! with Molly 

integer,  intent(in)  :: kind
logical,  intent(in)  :: is_state      ! True for state variable, false for obs
logical,  intent(in)  :: is_inflation  ! True for inflation transform
integer,  intent(out) :: dist_type
logical,  intent(out) :: bounded(2)
real(r8), intent(out) :: bounds(2)

! Have input information about the kind of the state or observation being transformed
! along with additional logical info that indicates whether this is an observation
! or state variable and about whether the transformation is being done for inflation
! or for regress. 
! Need to select the appropriate transform. At present, options are NORMAL_PRIOR
! which does nothing or BOUNDED_NORMAL_RH_PRIOR. 
! If the BNRH is selected then information about the bounds must also be set.
! The two dimensional logical array 'bounded' is set to false for no bounds and true
! for bounded. the first element of the array is for the lower bound, the second for the upper.
! If bounded is chosen, the corresponding bound value(s) must be set in the two dimensional 
! real array 'bounds'.
! For example, if my_state_kind corresponds to a sea ice fraction then an appropriate choice
! would be:
! bounded(1) = .true.;  bounded(2) = .true.
! bounds(1)  = 0.0_r8;  bounds(2)  = 1.0_r8

! In the long run, may not have to have separate controls for each of the input possibilities
! However, for now these are things that need to be explored for science understanding

if(is_inflation) then
   ! Case for inflation transformation
  SELECT CASE (kind)
    CASE (QTY_SEAICE_CONCENTR  , &
          QTY_SEAICE_AGREG_CONCENTR, &
          QTY_SEAICE_AICE01    , &
          QTY_SEAICE_AICE02    , &
          QTY_SEAICE_AICE03    , &
          QTY_SEAICE_AICE04    , &
          QTY_SEAICE_AICE05    )
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 1.0_r8
    CASE (QTY_SEAICE_VICE01)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 0.64_r8
    CASE (QTY_SEAICE_VICE02)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 1.39_r8
    CASE (QTY_SEAICE_VICE03)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 2.47_r8
    CASE (QTY_SEAICE_VICE04)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 4.57_r8
     CASE (QTY_SEAICE_VICE05)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(1) = .true.
      bounded(2) = .false.
      bounds(1) = 0.0_r8
    CASE (QTY_SEAICE_VOLUME, QTY_SEAICE_SNOWVOLUME, &
          QTY_SEAICE_AGREG_THICKNESS, QTY_SEAICE_AGREG_FREEBOARD, &
          QTY_SEAICE_VSNO02    , &
          QTY_SEAICE_VSNO03    , &
          QTY_SEAICE_VSNO04    , &
          QTY_SEAICE_VSNO05    )
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(1) = .true.
      bounded(2) = .false.
      bounds(1) = 0.0_r8
    CASE Default
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .false.
      return
  END SELECT
elseif(is_state) then
   ! Case for state variable priors
  SELECT CASE (kind)
    CASE (QTY_SEAICE_CONCENTR  , &
          QTY_SEAICE_AGREG_CONCENTR, &
          QTY_SEAICE_AICE01    , &
          QTY_SEAICE_AICE02    , &
          QTY_SEAICE_AICE03    , &
          QTY_SEAICE_AICE04    , &
          QTY_SEAICE_AICE05    )
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 1.0_r8
    CASE (QTY_SEAICE_VICE01)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 0.64_r8
    CASE (QTY_SEAICE_VICE02)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 1.39_r8
    CASE (QTY_SEAICE_VICE03)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 2.47_r8
    CASE (QTY_SEAICE_VICE04)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 4.57_r8
     CASE (QTY_SEAICE_VICE05)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(1) = .true.
      bounded(2) = .false.
      bounds(1) = 0.0_r8
    CASE (QTY_SEAICE_VOLUME, QTY_SEAICE_SNOWVOLUME, &
          QTY_SEAICE_AGREG_THICKNESS, QTY_SEAICE_AGREG_FREEBOARD,&
          QTY_SEAICE_VSNO02    , &
          QTY_SEAICE_VSNO03    , &
          QTY_SEAICE_VSNO04    , &
          QTY_SEAICE_VSNO05    )
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(1) = .true.
      bounded(2) = .false.
      bounds(1) = 0.0_r8
    CASE Default
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .false.
      return
  END SELECT
else
   ! This case is for observation (extended state) priors
  SELECT CASE (kind)
    CASE (QTY_SEAICE_CONCENTR  ,&
          QTY_SEAICE_AGREG_CONCENTR, &
          QTY_SEAICE_AICE01    , &
          QTY_SEAICE_AICE02    , &
          QTY_SEAICE_AICE03    , &
          QTY_SEAICE_AICE04    , &
          QTY_SEAICE_AICE05    )
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 1.0_r8
    CASE (QTY_SEAICE_VICE01)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 0.64_r8
    CASE (QTY_SEAICE_VICE02)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 1.39_r8
    CASE (QTY_SEAICE_VICE03)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 2.47_r8
    CASE (QTY_SEAICE_VICE04)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .true.
      bounds(1) = 0.0_r8
      bounds(2) = 4.57_r8
     CASE (QTY_SEAICE_VICE05)
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(1) = .true.
      bounded(2) = .false.
      bounds(1) = 0.0_r8
    CASE (QTY_SEAICE_VOLUME, QTY_SEAICE_SNOWVOLUME, &
          QTY_SEAICE_AGREG_THICKNESS, QTY_SEAICE_AGREG_FREEBOARD,&
          QTY_SEAICE_VSNO01    , &
          QTY_SEAICE_VSNO02    , &
          QTY_SEAICE_VSNO03    , &
          QTY_SEAICE_VSNO04    , &
          QTY_SEAICE_VSNO05    )
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(1) = .true.
      bounded(2) = .false.
      bounds(1) = 0.0_r8
    CASE Default
      dist_type = BOUNDED_NORMAL_RH_PRIOR
      bounded(:) = .false.
      return
  END SELECT
endif

end subroutine probit_dist_info

!------------------------------------------------------------------------


subroutine obs_inc_info(obs_kind, filter_kind, rectangular_quadrature, gaussian_likelihood_tails, &
   sort_obs_inc, spread_restoration, bounded, bounds)

integer,  intent(in)  :: obs_kind
integer,  intent(inout) :: filter_kind
logical,  intent(inout) :: rectangular_quadrature, gaussian_likelihood_tails
logical,  intent(inout) :: sort_obs_inc
logical,  intent(inout) :: spread_restoration
logical,  intent(inout) :: bounded(2)
real(r8), intent(inout) :: bounds(2)

! The information arguments are all intent (inout). This means that if they are not set
! here, they retain the default values from the assim_tools_mod namelist. Bounds don't exist 
! in that namelist, so default values are set in assim_tools_mod just before the call to here.

! Temporary approach for setting the details of how to assimilate this observation
! This example is designed to reproduce the squared forward operator results from paper


! Set the observation increment details for each type of quantity
SELECT CASE (obs_kind)
    CASE (QTY_SEAICE_AGREG_CONCENTR  , &
         QTY_SEAICE_AICE01    , &
         QTY_SEAICE_AICE02    , &
         QTY_SEAICE_AICE03    , &
         QTY_SEAICE_AICE04    , &
         QTY_SEAICE_AICE05    )
       filter_kind = BOUNDED_NORMAL_RHF
       bounded(:) = .true.
       bounds(1) = 0.0_r8
       bounds(2) = 1.0_r8
      CASE (QTY_SEAICE_VICE01)
        filter_kind = BOUNDED_NORMAL_RHF
        bounded(:) = .true.
        bounds(1) = 0.0_r8
        bounds(2) = 0.64_r8
      CASE (QTY_SEAICE_VICE02)
        filter_kind = BOUNDED_NORMAL_RHF
        bounded(:) = .true.
        bounds(1) = 0.0_r8
        bounds(2) = 1.39_r8
      CASE (QTY_SEAICE_VICE03)
        filter_kind = BOUNDED_NORMAL_RHF
        bounded(:) = .true.
        bounds(1) = 0.0_r8
        bounds(2) = 2.47_r8
      CASE (QTY_SEAICE_VICE04)
        filter_kind = BOUNDED_NORMAL_RHF
        bounded(:) = .true.
        bounds(1) = 0.0_r8
        bounds(2) = 4.57_r8
       CASE (QTY_SEAICE_VICE05)
        filter_kind = BOUNDED_NORMAL_RHF
        bounded(1) = .true.
        bounded(2) = .false.
        bounds(1) = 0.0_r8     
    CASE (QTY_SEAICE_AGREG_THICKNESS, &
         QTY_SEAICE_AGREG_FREEBOARD, &
         QTY_SEAICE_VSNO01    , &
         QTY_SEAICE_VSNO02    , &
         QTY_SEAICE_VSNO03    , &
         QTY_SEAICE_VSNO04    , &
         QTY_SEAICE_VSNO05    ) 
       filter_kind = BOUNDED_NORMAL_RHF
       bounded(1) = .true.
       bounded(2) = .false.
       bounds(1) = 0.0_r8;
   CASE Default
      filter_kind = BOUNDED_NORMAL_RHF
      bounded(:) = .false.
      return
  END SELECT

! Default settings for now for Icepack and tracer model tests
sort_obs_inc = .false.
spread_restoration = .false.

! Only need to set these two for options on old RHF implementation
! rectangular_quadrature = .true.
! gaussian_likelihood_tails = .false.

end subroutine obs_inc_info

!------------------------------------------------------------------------

end module algorithm_info_mod
