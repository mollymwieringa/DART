module fractional_reg_mod

    use types_mod,               only : r8
    use utilities_mod,           only : E_ERR, E_MSG, error_handler
    use assim_tools_mod,         only : update_from_obs_inc
    use probit_transform_mod,    only : transform_to_probit, transform_from_probit
    use distribution_params_mod, only : distribution_params_type

    implicit none
    private

    public :: state_regress_disaggregation, &
              state_regress_relativefrac, &
              state_regress_probit, &
              postprocess

    contains

    ! I don't really want to have to think about redefining distributions. We need the distibutions to 
    ! help define the approach the probit transforms take. Can I load them in from somewhere?

    ! The distributions can be fed in to the functions and used from the distribution_params_mod in the 
    ! mod that calls these subroutines. Here, we just need to be able to define the type for variables that hold 
    ! the parameter information for ensembles of a specific distribution type. 
    !----------------------------------------------------------------------------------!
    
    ! Subroutine to perform regression by disaggregation
    ! Does not care about the distribution of the ensemble
    subroutine state_regress_disaggregation(obs_inc, ens_prior, ens_post, ens_size, nc, &
                                            bounded_above, bounded_below, upper_bound, lower_bound)

        ! Declare routine variables
        integer,   intent(in) :: ens_size, nc
        logical,   intent(in) :: bounded_above, bounded_below
        real(r8),  intent(in) :: upper_bound, lower_bound
        real(r8),  intent(in) :: ens_prior(ens_size, nc)
        real(r8),  intent(in) :: obs_inc(ens_size)     
        real(r8), intent(out) :: ens_post(ens_size, nc) 
        
        ! Declare local variables
        character(len=100)    :: errstring
        integer               :: i, j
        real(r8)              :: weights(ens_size, nc), state_inc(ens_size, nc)

        ! Perform regression by disaggregation
        do i = 1, ens_size
            ! Determine the weights
            weights(i, :) = ens_prior(i, :) / sum(ens_prior(i, :))
            ! Disaggregate the observation increment to state variable increments
            state_inc(i, :) = obs_inc(i) * weights(i, :)
        end do

        ! Calculate the ensemble posterior state
        ens_post = ens_prior + state_inc

        ! Verify that the updated state variables are withing bounds
        ! NOTE: being a bit lazy here and only checking the aggregate. Come back later to
        !      make sure that individual categories are also within bounds if needed.
        do i = 1, ens_size
            if (bounded_below) then 
                if (sum(ens_post(i, :)) < lower_bound) then
                    errstring = "Aggregate of state variables in an ensemble member results in lower bound violation."
                    call error_handler(E_ERR, 'state_regress_disaggregation', trim(errstring))
                end if
            end if
            if (bounded_above) then
                if (sum(ens_post(i, :)) > upper_bound) then
                    errstring = "Aggregate of state variables in an ensemble member results in upper bound violation."
                    call error_handler(E_ERR, 'state_regress_disaggregation', trim(errstring))
                end if
            end if
        end do

    end subroutine state_regress_disaggregation

    ! Subroutine to perform regression with relative fractional amount constraint
    ! Does not care about the distribution of the ensemble
    subroutine state_regress_relativefrac(obs_prior, obs_inc, ens_prior, ens_post, ens_size, nc, &
                                          bounded_above, bounded_below, upper_bound, lower_bound)

        ! This routine assumes that the state variables are fractional amounts that sum to one.
        ! If they are not, the method is invalid. For now, we will leave it to the user to ensure
        ! that the incoming ens_prior meets this requirement.
        ! In the case of representing a quantity whose aggregate can vary between 0 and 1, one 
        ! valid approach would be to calculate and attach a "null" category to the state vector.

        ! Declare routine variables
        integer,   intent(in) :: ens_size, nc
        real(r8),  intent(in) :: obs_prior(ens_size), obs_inc(ens_size)
        real(r8),  intent(in) :: ens_prior(ens_size, nc)
        real(r8), intent(out) :: ens_post(ens_size, nc)
        logical,   intent(in) :: bounded_above, bounded_below
        real(r8),  intent(in) :: upper_bound, lower_bound

        ! Declare local variables
        character(len=100)    :: errstring
        integer               :: i, j, nc_exp
        real(r8)              :: exp_prior(ens_size, nc+1)
        real(r8), allocatable :: z_prior(:, :)
        real(r8), allocatable :: z_post(:, :), z_inc(:, :) 
        real(r8)              :: a_prior(ens_size), a_post(ens_size)
        real(r8)              :: net_a, reg_coef, obs_prior_mean, obs_prior_var


        ! Calculate the fractional amounts
        ! Verify that the sumer of each ensemble member in ens_prior is 1.0
        ! If not, either add a "null category" if less than 1.0 OR squash down to relative
        ! fractional amounts and save the squashing factor.

        nc_exp = nc + 1
        if (all(sum(ens_prior, dim=1) < 1.0_r8)) then
            allocate(z_prior(ens_size, nc_exp))
            allocate(z_post(ens_size, nc_exp))
            allocate(z_inc(ens_size, nc_exp))

            do i = 1, ens_size
                exp_prior(i, 1:nc) = ens_prior(i, :)
                exp_prior(i, nc_exp) = 1.0_r8 - sum(ens_prior(i, :))
                a_prior(i) = 1.0_r8
                z_prior(i, :) = a_prior(i) * exp_prior(i, :)
            end do
        else
            allocate(z_prior(ens_size, nc))
            allocate(z_post(ens_size, nc))
            allocate(z_inc(ens_size, nc))

            do i = 1, ens_size
                if (sum(ens_prior(i,:)) /= 1.0_r8) then
                    a_prior(i) = 1.0_r8 / sum(ens_prior(i,:))
                else
                    a_prior(i) = 1.0_r8
                end if
                z_prior(i,:) = a_prior(i) * ens_prior(i,:)
            end do
        end if

        ! Perform regression in relative fractional space
        net_a = 1.0_r8    ! This is a null value; net_a is not used at this time.
        obs_prior_mean = sum(obs_prior) / ens_size
        obs_prior_var = sum(obs_prior - obs_prior_mean)**2 / (ens_size - 1)
        call update_from_obs_inc(obs_prior, obs_prior_mean, obs_prior_var, & 
                                 obs_inc, z_prior, ens_size, z_inc, &
                                 reg_coef, net_a)

        ! Calculate posterior fractional amounts
        z_post = z_prior + z_inc
        do i = 1, ens_size
            a_post(i) = sum(z_post(i, :))
            ens_post(i, :) = z_post(i, 1:nc) / (a_post(i) * a_prior(i))
        end do
        
        ! ! if necessary, cut out the null category and convert back to state variable space
        ! ens_post = z_post(:, 1:nc) / (a_post * a_prior)

        ! Verify that the updated state variables are within bounds
        do i = 1, ens_size
            ! Check A: verify that all individual state variables are within bounds
            do j = 1, nc
                if (bounded_below) then
                    if (ens_post(i, j) < lower_bound) then
                        errstring = "State variable violates lower bound."
                        call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                    end if
                end if
                if (bounded_above) then
                    if (ens_post(i, j) > upper_bound) then
                        errstring = "State variable violates upper bound."
                        call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                    end if
                end if
            end do
            ! Check B: verify that the sum of all state variables is also within bounds
            if (bounded_below) then
                if (sum(ens_post(i, :)) < lower_bound) then
                    errstring = 'Aggregate of state variables violates lower bound.'
                    call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                end if
            end if
            if (bounded_above) then
                if (sum(ens_post(i, :)) > upper_bound) then
                    errstring = 'Aggregate of state variables violates upper bound.'
                    call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                end if
            end if
        end do

        deallocate(z_prior, z_post, z_inc)

    end subroutine state_regress_relativefrac

    ! Subroutine to perform regression in probit space (DART standard)
    ! Because of the transform in probit space, we do care about the distribution / parameters of the ensemble.
    subroutine state_regress_probit(obs_prior, obs_post, ens_prior, ens_post, ens_size, nc, &
                                    dist_for_obs, dist_for_state, &
                                    bounded_above, bounded_below, upper_bound, lower_bound)

        ! temp_dist_params, and state_dist_params need to be figured out... 
        ! dist_for_obs should be NORMAL_DISTRIBUTION
        ! dist_for_state could be any of NORMAL DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, KDE_DISTRIBUTION

        ! the to_probit_bounded_normal_rhf (etc) called from transform_to_probit calls a function that gathers information about 
        ! the ensemble distribution parameters (p). p needs to go into that function with the distribution assignment, but the rest
        ! of the parameters are filled in within that function. 

        ! p goes in as a type (distribution_params_type) with dimensions for the number of variables in the state / obs list
        !   however, this is otherwise empty to begin

        ! Declare routine variables
        integer,            intent(in) :: ens_size, nc
        integer,            intent(in) :: dist_for_obs, dist_for_state
        logical,            intent(in) :: bounded_above, bounded_below
        real(r8),           intent(in) :: obs_prior(ens_size), obs_post(ens_size)
        real(r8),           intent(in) :: ens_prior(ens_size, nc)
        real(r8),          intent(out) :: ens_post(ens_size, nc)
        real(r8),           intent(in) :: upper_bound, lower_bound

        ! Declare local variables
        character(len=100)             :: errstring
        integer                        :: i, j, ierr
        real(r8)                       :: net_a, reg_coef
        real(r8)                       :: probit_obs_prior_mean, probit_obs_prior_var
        real(r8)                       :: probit_obs_inc(ens_size)
        real(r8)                       :: probit_obs_prior(ens_size), probit_obs_post(ens_size)
        real(r8)                       :: probit_state_inc(ens_size, nc)
        real(r8)                       :: probit_ens_prior(ens_size, nc), probit_ens_post(ens_size, nc)
        type(distribution_params_type) :: obs_dist_params, state_dist_params

        ! 1. Transform the observation space update into probit space and recalculate the increments
        call transform_to_probit(ens_size, obs_prior, dist_for_obs, obs_dist_params, &
                                 probit_obs_prior, .false., bounded_below, bounded_above, &
                                 lower_bound, upper_bound, ierr)
        if (ierr /= 0) then
            errstring = "Error in transform_to_probit for observation prior"
            call error_handler(E_ERR, 'state_regress_probit', trim(errstring))
        end if
        call transform_to_probit(ens_size, obs_post, dist_for_obs, obs_dist_params, &
                                 probit_obs_post, .true., bounded_below, bounded_above, &
                                 lower_bound, upper_bound, ierr)
        if (ierr /= 0) then
            errstring = "Error in transform_to_probit for observation posterior"
            call error_handler(E_ERR, 'state_regress_probit', trim(errstring))
        end if

        probit_obs_inc = probit_obs_post - probit_obs_prior
        probit_obs_prior_mean = sum(probit_obs_prior) / ens_size
        probit_obs_prior_var = sum((probit_obs_prior - probit_obs_prior_mean)**2) / (ens_size - 1)

        ! 2. Transform the state variables to probit space
        call transform_to_probit(ens_size, ens_prior, dist_for_state, state_dist_params, &
                                 probit_ens_prior, .false., bounded_below, bounded_above, &
                                 lower_bound, upper_bound, ierr)
        if (ierr /= 0) then
            errstring = "Error in transform_to_probit for state prior"
            call error_handler(E_ERR, 'state_regress_probit', trim(errstring))
        end if

        ! 3. Call obs_updates ens
        net_a = 1.0_r8   ! This is a null value; net_a is not used at this time.
        call update_from_obs_inc(probit_obs_prior, probit_obs_prior_mean, probit_obs_prior_var, &
                                 probit_obs_inc, probit_ens_prior, ens_size, probit_state_inc, &
                                 reg_coef, net_a)

        probit_ens_post = probit_ens_prior + probit_state_inc

        ! 4. Transform back from probit space to original space
        call transform_from_probit(ens_size, probit_ens_post, state_dist_params, ens_post)

    end subroutine state_regress_probit

    ! Subroutine to perform CICE-style postprocessing after the "regression" step
    subroutine postprocess(ens_post, ens_size, nc, bounded_above, bounded_below, upper_bound, lower_bound)

        ! Declare routine variables
        integer,       intent(in) :: ens_size, nc
        logical,       intent(in) :: bounded_above, bounded_below
        real(r8),   intent(inout) :: ens_post(ens_size, nc) ! aicen
        real(r8),      intent(in) :: upper_bound, lower_bound

        ! Declare local variables 
        real(r8)                  :: agg, agg_temp, squeeze
        integer                   :: i, j
 
        ! Begin process
        do i = 1, ens_size

            if (bounded_above) then
               ens_post(i,:) = min(upper_bound, ens_post(i,:))   ! individual categories must not exceed 1
            end if

            ! calculate aggregates for posterior ensemble
            agg = sum(ens_post(i, :))

            if (bounded_below) then
                ens_post(i, :) = max(lower_bound, ens_post(i,:)) ! individual categories must be non-negative
            end if
            
            ! recalculate aggregate once bounds are enforced
            agg_temp = sum(ens_post(i,:))

            ! Begin squeezing, if necessary
            if (agg <= lower_bound) then
                ens_post(i,:) = lower_bound
            else if (agg > lower_bound) then
                do j = 1, nc
                    ! if both the constrained and unconstrained aggregates are greater than the lower bound
                    if (agg_temp > lower_bound .and. agg > lower_bound) then
                        ! adjust each category to account for any potential "negative" volume
                        ens_post(i,j) = ens_post(i,j) - (agg_temp - agg)*ens_post(i,j)/agg_temp
                    endif
                end do
                
                ! recalculate the aggregate area given the recovered aggregate
                agg = sum(ens_post(i,:))

                ! if the aggregate violates the upper bound
                if (agg > upper_bound) then
                    ! squeeze all the categories down such that the posterior categories sum to the upper bound
                    squeeze = upper_bound / agg
                    ens_post(i,:) = ens_post(i,:) * squeeze
                end if
            end if 
        end do 
    
    end subroutine postprocess

end module fractional_reg_mod