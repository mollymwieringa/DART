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
    subroutine state_regress_disaggregation(obs_prior, obs_post, ens_prior, ens_post, ens_size, nc, &
                                            obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                            state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)

        ! Declare routine variables
        integer,   intent(in) :: ens_size, nc
        logical,   intent(in) :: state_bounded_above, state_bounded_below
        logical,   intent(in) :: obs_bounded_above, obs_bounded_below
        real(r8),  intent(in) :: state_upper_bound, state_lower_bound
        real(r8),  intent(in) :: obs_upper_bound, obs_lower_bound
        real(r8),  intent(in) :: ens_prior(ens_size, nc)
        real(r8),  intent(in) :: obs_prior(ens_size), obs_post(ens_size)    
        real(r8), intent(out) :: ens_post(ens_size, nc) 
        
        ! Declare local variables
        character(len=128)    :: errstring
        integer               :: i, j
        real(r8)              :: weights(ens_size)

        ! Perform regression by disaggregation
        do i = 1, ens_size
            ! We want to be aware of the cases where the obs_prior is close to 0 
            if (obs_prior(i) < 1e-8_r8) then
                ! if obs_posterior is also zero, set weights to zero
                if (obs_post(i) < 1e-8_r8) then
                    ens_post(i, :) = 0.0_r8
                ! else, if the DA process "grew mass", distribute that mass evenly across categories
                ! this is probably fine if we assume the new "mass" is relatively small (no massive jumps)
                ! this is not going to be true for a system like sea ice, where new ice will tend toward the first 
                ! thickness category
                else
                    ens_post(i, :) = obs_post(i) / nc
                end if
            else
                ! calculate weights as a scale factor and then apply the updates to the state
                weights(i) = obs_post(i) / obs_prior(i)
                ens_post(i, :) = weights(i) * ens_prior(i, :)
            end if
        end do

        ! Verify that the updated state variables are within bounds
        do i = 1, ens_size
            ! Check A: verify that all individual state variables are within bounds
            do j = 1, nc
                if (state_bounded_below) then
                    if (ens_post(i, j) < state_lower_bound) then
                        write(errstring, *) "State variable ", j, " in ensemble member ", i, "violates lower bound."
                        write(*, *) 'ensemble at error: ', ens_post
                        call error_handler(E_MSG, 'state_regress_disaggregation', trim(errstring))
                    end if
                end if
                if (state_bounded_above) then
                    if (ens_post(i, j) > state_upper_bound) then
                        write(errstring, *) "State variable ", j, " in ensemble member ", i, "violates upper bound."
                        write(*, *) 'ensemble at error: ', ens_post
                        call error_handler(E_MSG, 'state_regress_disaggregation', trim(errstring))
                    end if
                end if
            end do
            ! Check B: verify that the sum of all state variables is also within bounds
            if (state_bounded_below) then 
                if (sum(ens_post(i, :)) - state_lower_bound < -1e-8) then
                    write(errstring, *) "Aggregate in ensemble member ", i, "results in lower bound violation. Aggregate at error: ", sum(ens_post(i,:))
                    call error_handler(E_MSG, 'state_regress_disaggregation', trim(errstring))
                end if
            end if
            if (state_bounded_above) then
                ! NOTE that this is not totally correct, as the upper bound on the aggregate may be different from the sum of the upper bounds on each category
                !        in some applications. For now (in sea ice, where ub_agg == ub_cat for SIC), this is acceptable. 
                if (sum(ens_post(i, :)) - state_upper_bound > 1e-8) then
                    write(errstring, *) "Aggregate in ensemble member ", i, "results in upper bound violation. Aggregate at error: ", sum(ens_post(i,:))
                    call error_handler(E_MSG, 'state_regress_disaggregation', trim(errstring))
                end if
            end if
        end do

    end subroutine state_regress_disaggregation

    ! Subroutine to perform regression with relative fractional amount constraint
    ! Does not care about the distribution of the ensemble
    subroutine state_regress_relativefrac(obs_prior, obs_post, ens_prior, ens_post, ens_size, nc, &
                                          dist_for_obs, dist_for_state, &
                                          obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                          state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)

        ! This routine assumes that the state variables are fractional amounts that sum to one.
        ! If they are not, the method is invalid. For now, we will leave it to the user to ensure
        ! that the incoming ens_prior meets this requirement.
        ! In the case of representing a quantity whose aggregate can vary between 0 and 1, one 
        ! valid approach would be to calculate and attach a "null" category to the state vector.

        ! Declare routine variables
        integer,   intent(in) :: ens_size, nc
        integer,   intent(in) :: dist_for_obs, dist_for_state
        logical,   intent(in) :: state_bounded_above, state_bounded_below
        logical,   intent(in) :: obs_bounded_above, obs_bounded_below
        real(r8),  intent(in) :: obs_prior(ens_size), obs_post(ens_size)
        real(r8),  intent(in) :: ens_prior(ens_size, nc)
        real(r8), intent(out) :: ens_post(ens_size, nc)
        real(r8),  intent(in) :: state_upper_bound, state_lower_bound
        real(r8),  intent(in) :: obs_upper_bound, obs_lower_bound
        ! Declare local variables
        character(len=100)    :: errstring
        integer               :: i, j, nc_exp
        real(r8)              :: exp_prior(ens_size, nc+1)
        real(r8), allocatable :: xhat_prior(:, :) !, reg_coef(:)
        real(r8), allocatable :: xhat_post(:, :) !, xhat_inc(:, :) 
        real(r8)              :: a_prior(ens_size), a_post(ens_size), xhat_post_j(ens_size)
        ! real(r8)              :: net_a, reg_coef_j, obs_prior_mean, obs_prior_var


        ! net_a = 1.0_r8    ! This is a null value; net_a is not used at this time.
        ! obs_prior_mean = sum(obs_prior) / ens_size
        ! obs_prior_var = sum((obs_prior - obs_prior_mean)**2) / (ens_size - 1)

        ! Calculate the fractional amounts
        ! Verify that the sum of each ensemble member in ens_prior is 1.0
        ! If not, either add a "null category" if less than 1.0 OR squash down to relative
        ! fractional amounts and save the squashing factor.

        nc_exp = nc + 1
        write(*,*) 'Aggregates for ens are ', sum(ens_prior, dim=2)
        if (all(sum(ens_prior, dim=2) < 1.0_r8)) then
            write(*, *) 'Expanding dimensions to include null space...'
            allocate(xhat_prior(ens_size, nc_exp))
            allocate(xhat_post(ens_size, nc_exp))
            ! allocate(xhat_inc(ens_size, nc_exp))
            ! allocate(reg_coef(nc_exp))

            do i = 1, ens_size
                exp_prior(i, 1:nc) = ens_prior(i, :)
                exp_prior(i, nc_exp) = 1.0_r8 - sum(ens_prior(i, :))
                a_prior(i) = 1.0_r8
                xhat_prior(i, :) = a_prior(i) * exp_prior(i, :)
            end do

            do j = 1, nc_exp
                call state_regress_probit(obs_prior, obs_post, xhat_prior(:,j), xhat_post_j, ens_size, 1, &
                                          dist_for_obs, dist_for_state, &
                                          obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                          state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)
                ! call update_from_obs_inc(obs_prior, obs_prior_mean, obs_prior_var, &
                !                          obs_inc, xhat_prior(:,j), ens_size, xhat_inc_j, &
                !                          reg_coef_j, net_a)
                xhat_post(:,j) = xhat_post_j
            end do
        else
            allocate(xhat_prior(ens_size, nc))
            allocate(xhat_post(ens_size, nc))
            ! allocate(xhat_inc(ens_size, nc))
            ! allocate(reg_coef(nc))

            do i = 1, ens_size
                if (sum(ens_prior(i,:)) /= 1.0_r8) then
                    a_prior(i) = 1.0_r8 / sum(ens_prior(i,:))
                else
                    a_prior(i) = 1.0_r8
                end if
                xhat_prior(i,:) = a_prior(i) * ens_prior(i,:)
            end do

            do j = 1, nc
                call state_regress_probit(obs_prior, obs_post, xhat_prior(:,j), xhat_post_j, ens_size, 1, &
                                          dist_for_obs, dist_for_state, &
                                          obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                          state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)
                ! call update_from_obs_inc(obs_prior, obs_prior_mean, obs_prior_var, &
                !                          obs_inc, xhat_prior(:,j), ens_size, xhat_inc_j, &
                !                          reg_coef_j, net_a)
                xhat_post(:,j) = xhat_post_j
            end do
        end if

        write(*,*) 'a from the prior is ', a_prior

        ! ! Perform regression in relative fractional space 
        ! net_a = 1.0_r8    ! This is a null value; net_a is not used at this time.
        ! obs_prior_mean = sum(obs_prior) / ens_size
        ! obs_prior_var = sum((obs_prior - obs_prior_mean)**2) / (ens_size - 1)
        ! call update_from_obs_inc(obs_prior, obs_prior_mean, obs_prior_var, & 
        !                          obs_inc, xhat_prior, ens_size, xhat_inc, &
        !                          reg_coef, net_a)

        ! write(*,*) 'reg_coef is ', reg_coef

        ! Calculate posterior fractional amounts
        ! xhat_post = xhat_prior + xhat_inc
        do i = 1, ens_size
            a_post(i) = sum(xhat_post(i, :))
            ens_post(i, :) = xhat_post(i, 1:nc) / (a_post(i) * a_prior(i))
        end do

        write(*,*) 'a from the posterior is ', a_post

        ! Verify that the updated state variables are within bounds
        do i = 1, ens_size
            ! Check A: verify that all individual state variables are within bounds
            do j = 1, nc
                if (state_bounded_below) then
                    if (ens_post(i, j) < state_lower_bound) then
                        write(errstring, *) "State variable ", j, " in ensemble member ", i, "violates lower bound."
                        write(*, *) 'ensemble at error: ', ens_post
                        call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                    end if
                end if
                if (state_bounded_above) then
                    if (ens_post(i, j) > state_upper_bound) then
                        write(errstring, *) "State variable ", j, " in ensemble member ", i, "violates upper bound."
                        write(*, *) 'ensemble at error: ', ens_post
                        call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                    end if
                end if
            end do
            ! Check B: verify that the sum of all state variables is also within bounds
            if (state_bounded_below) then
                if (sum(ens_post(i, :)) - state_lower_bound < -1e-8) then
                    write(errstring, *) 'Aggregate of state variable in ensemble member ', i, 'violates lower bound. Aggreate at error: ', sum(ens_post(i,:))
                    call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                end if
            end if
            if (state_bounded_above) then
                ! NOTE that this is not totally correct, as the upper bound on the aggregate may be different from the sum of the upper bounds on each category
                !        in some applications. For now (in sea ice, where ub_agg == ub_cat for SIC), this is acceptable.
                if (sum(ens_post(i, :)) - state_upper_bound > 1e-8) then
                    write(errstring, *) 'Aggregate of state variables in ensemble member ', i, 'violates upper bound. Aggregate at error: ', sum(ens_post(i,:))
                    call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                end if
            end if
        end do

        deallocate(xhat_prior, xhat_post)

    end subroutine state_regress_relativefrac

    ! Subroutine to perform regression in probit space (DART standard)
    ! Because of the transform in probit space, we do care about the distribution / parameters of the ensemble.
    subroutine state_regress_probit(obs_prior, obs_post, ens_prior, ens_post, ens_size, nc, &
                                    dist_for_obs, dist_for_state, &
                                    obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                    state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)

        ! Declare routine variables
        integer,            intent(in) :: ens_size, nc
        integer,            intent(in) :: dist_for_obs, dist_for_state
        logical,            intent(in) :: obs_bounded_above, obs_bounded_below
        logical,            intent(in) :: state_bounded_above, state_bounded_below
        real(r8),           intent(in) :: obs_prior(ens_size), obs_post(ens_size)
        real(r8),           intent(in) :: ens_prior(ens_size, nc)
        real(r8),          intent(out) :: ens_post(ens_size, nc)
        real(r8),           intent(in) :: state_upper_bound, state_lower_bound
        real(r8),           intent(in) :: obs_upper_bound, obs_lower_bound
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
        write(*,*) 'Transforming prior to probit...'
        call transform_to_probit(ens_size, obs_prior, dist_for_obs, obs_dist_params, &
                                 probit_obs_prior, .false., obs_bounded_below, obs_bounded_above, &
                                 obs_lower_bound, obs_upper_bound, ierr)
        if (ierr /= 0) then
            ens_post = ens_prior
            write(errstring, *) "Error in transform_to_probit for observation prior, exiting w no update.. "
            call error_handler(E_MSG, 'state_regress_probit', trim(errstring))
        else
            write(*,*) 'Transforming posterior to probit...'
            call transform_to_probit(ens_size, obs_post, dist_for_obs, obs_dist_params, &
                                 probit_obs_post, .true., obs_bounded_below, obs_bounded_above, &
                                 obs_lower_bound, obs_upper_bound, ierr)
            if (ierr /= 0) then
                ens_post = ens_prior
                write(errstring, *) "Error in transform_to_probit for observation posterior, exiting w no update... "
                call error_handler(E_MSG, 'state_regress_probit', trim(errstring))
            else
                probit_obs_inc = probit_obs_post - probit_obs_prior
                probit_obs_prior_mean = sum(probit_obs_prior) / ens_size
                probit_obs_prior_var = sum((probit_obs_prior - probit_obs_prior_mean)**2) / (ens_size - 1)

                ! 2. Transform the state variables to probit space
                write(*,*) 'Transforming prior state to probit...'
                call transform_to_probit(ens_size, ens_prior, dist_for_state, state_dist_params, &
                                         probit_ens_prior, .false., state_bounded_below, state_bounded_above, &
                                         state_lower_bound, state_upper_bound, ierr)
                if (ierr /= 0) then
                    ens_post = ens_prior
                    write(errstring, *) "Error in transform_to_probit for state prior, exiting w no update..."
                    call error_handler(E_MSG, 'state_regress_probit', trim(errstring))
                else
                    ! 3. Call obs_updates ens
                    net_a = 1.0_r8   ! This is a null value; net_a is not used at this time.
                    call update_from_obs_inc(probit_obs_prior, probit_obs_prior_mean, probit_obs_prior_var, &
                                             probit_obs_inc, probit_ens_prior, ens_size, probit_state_inc, &
                                             reg_coef, net_a)

                    probit_ens_post = probit_ens_prior + probit_state_inc

                    ! 4. Transform back from probit space to original space
                    write(*,*) 'Transforming all back from probit...'
                    call transform_from_probit(ens_size, probit_ens_post, state_dist_params, ens_post)
                end if
            end if
        end if

    end subroutine state_regress_probit

    ! Subroutine to perform CICE-style postprocessing after the "regression" step
    subroutine postprocess(ens_post, ens_size, nc, state_bounded_above, state_bounded_below, &
                           state_upper_bound, state_lower_bound, &
                           obs_bounded_above, obs_bounded_below, &
                           obs_upper_bound, obs_lower_bound)

        ! Declare routine variables
        integer,       intent(in) :: ens_size, nc
        logical,       intent(in) :: state_bounded_above, state_bounded_below
        logical,       intent(in) :: obs_bounded_above, obs_bounded_below
        real(r8),   intent(inout) :: ens_post(ens_size, nc) ! aicen
        real(r8),      intent(in) :: state_upper_bound, state_lower_bound
        real(r8),      intent(in) :: obs_upper_bound, obs_lower_bound

        ! Declare local variables 
        real(r8)                  :: agg, agg_temp, squeeze
        integer                   :: i, j
 
        ! Begin process
        do i = 1, ens_size

            if (state_bounded_above) then
               ens_post(i,:) = min(state_upper_bound, ens_post(i,:))   ! individual categories must not exceed 1
            end if

            ! calculate aggregates for posterior ensemble
            agg = sum(ens_post(i, :))

            if (state_bounded_below) then
                ens_post(i, :) = max(state_lower_bound, ens_post(i,:)) ! individual categories must be non-negative
            end if
            
            ! recalculate aggregate once bounds are enforced
            agg_temp = sum(ens_post(i,:))

            ! Begin squeezing, if necessary
            if (agg <= obs_lower_bound) then
                ens_post(i,:) = obs_lower_bound
            else if (agg > obs_lower_bound) then
                do j = 1, nc
                    ! if both the constrained and unconstrained aggregates are greater than the lower bound
                    if (agg_temp > obs_lower_bound .and. agg > obs_lower_bound) then
                        ! adjust each category to account for any potential "negative" volume
                        ens_post(i,j) = ens_post(i,j) - (agg_temp - agg)*ens_post(i,j)/agg_temp
                    endif
                end do
                
                ! recalculate the aggregate area given the recovered aggregate
                agg = sum(ens_post(i,:))

                if (state_bounded_above) then
                !    if the aggregate violates the upper bound
                    ! NOTE that this is not totally correct, as the upper bound on the aggregate may be different from the sum of the upper bounds on each category
                    !        in some applications. For now (in sea ice, where ub_agg == ub_cat for SIC), this is acceptable.
                    if (agg > state_upper_bound) then
                        ! squeeze all the categories down such that the posterior categories sum to the upper bound
                        squeeze = state_upper_bound / agg
                        ens_post(i,:) = ens_post(i,:) * squeeze
                    end if
                end if
            end if 
        end do 
    
    end subroutine postprocess

end module fractional_reg_mod