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
              get_delta_representation_from_state, &
              get_state_from_delta_representation, &
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

        ! This routine assumes that the state variables are fractional amounts whose sum is constrained to be 
        ! to be less than or equal to 1.0. It assumes observations used to update the system are of the 
        ! sum of all state categories, and thus also reflect change in a "null category" that is not explicitly
        ! represented in the state vector. 

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
        integer               :: i, j
        real(r8)              :: mass_check
        real(r8)              :: xhat_prior(ens_size, nc), xhat_post(ens_size, nc)
        real(r8)              :: a_post(ens_size), xhat_post_j(ens_size), Np(ens_size), Na(ens_size), obs_inc(ens_size)

        ! Translate x_prior into xhat_prior (assumes multiplication by a_prior, which is always 1.0)
        xhat_prior = ens_prior

        ! Calculate the prior null category, Np
        Np = 1.0_r8 - sum(ens_prior, dim=2)

        ! Calculate the posterior null category, Na 
        obs_inc = obs_post - obs_prior
        Na = Np - obs_inc

        ! Regress in the space of the fractional amounts
        do j = 1, nc
            call state_regress_probit(obs_prior, obs_post, xhat_prior(:,j), xhat_post_j, ens_size, 1, &
                                      dist_for_obs, dist_for_state, &
                                      obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                      state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)
            xhat_post(:,j) = xhat_post_j
        end do

        ! ! Check that the regression reflects the observation increment across categories
        ! do i = 1, ens_size
        !     ! if the summed change to categories differs from the obs inc by more than 0.1
        !     mass_check = obs_inc(i) - sum(xhat_post(i, :) - xhat_prior(i, :))
        !     if (abs(mass_check) > 0.1_r8) then
        !         write(errstring, *) "Summed change to categories for ensemble member ", i, " differs from obs inc by more than 0.1."
        !         call error_handler(E_MSG, 'state_regress_relativefrac', trim(errstring))
        !         if (mass_check > 0.0_r8) then
        !             write(errstring, *) "regression did not add enough mass to the system. Mass check value: ", mass_check
        !             call error_handler(E_MSG, 'state_regress_relativefrac', trim(errstring))
        !             xhat_post(i, 1) = xhat_post(i, 1) + mass_check * Np(i) / (Np(i) + sum(xhat_prior(i, 2:)))
        !         else if (mass_check < 0.0_r8) then
        !             write(errstring, *) "regression added too much mass to the system. Mass check value: ", mass_check
        !             call error_handler(E_MSG, 'state_regress_relativefrac', trim(errstring))
        !             xhat_post(i, 1) = xhat_post(i, 1) + mass_check * xhat_post(i, 1) / sum(xhat_post(i, :))
        !         else
        !             write(errstring, *) "regression added the correct amount of mass to the system. Mass check value: ", mass_check
        !             call error_handler(E_MSG, 'state_regress_relativefrac', trim(errstring))
        !         end if
        !     end if
        ! end do

        ! Calculate a_post and use it to translate back from xhat_post to ens_post
        open(unit=23, file=trim('alpha_posterior.txt'), status='UNKNOWN', RECL=256)
        do i = 1, ens_size
            if (abs(Na(i) - 1.0_r8) <= 1e-16_r8) then
                write(errstring, *) "Posterior null category for ensemble member ", i, " is 1.0 (to roundoff)"
                call error_handler(E_MSG, 'state_regress_relativefrac', trim(errstring))
                write(23, *) 1e16_r8
                ens_post(i, :) = 0.0_r8
            else if (sum(xhat_post(i, :)) < 1e-16_r8) then
                write(errstring, *) "Posterior xhat values for ensemble member ", i, " sum to zero (to roundoff)"
                call error_handler(E_MSG, 'state_regress_relativefrac', trim(errstring))
                write(23, *) 0.0_r8
                ens_post(i, :) = 0.0_r8
            else
                a_post(i) = sum(xhat_post(i, :) / (1.0_r8 - Na(i)))
                if (a_post(i) < 0.0_r8) then
                    write(errstring, *) "Posterior alpha for ensemble member ", i, " is negative."
                    call error_handler(E_MSG, 'state_regress_relativefrac', trim(errstring))
                    write(errstring, *) "Posterior N value: ", Na(i), " Posterior xhat values: ", sum(xhat_post(i, :))
                    call error_handler(E_MSG, 'state_regress_relativefrac', trim(errstring))
                end if
                write(23, *) a_post(i)
                ens_post(i, :) = xhat_post(i, :) / a_post(i)
            end if
        end do
        close(23)

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
                if (sum(ens_post(i, :)) - state_lower_bound < -1e-8_r8) then
                    write(errstring, *) 'Aggregate of state variable in ensemble member ', i, 'violates lower bound. Aggreate at error: ', sum(ens_post(i,:))
                    call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                end if
            end if
            if (state_bounded_above) then
                ! NOTE that this is not totally correct, as the upper bound on the aggregate may be different from the sum of the upper bounds on each category
                !        in some applications. For now (in sea ice, where ub_agg == ub_cat for SIC), this is acceptable.
                if (sum(ens_post(i, :)) - state_upper_bound > 1e-8_r8) then
                    write(errstring, *) 'Aggregate of state variables in ensemble member ', i, 'violates upper bound. Aggregate at error: ', sum(ens_post(i,:))
                    call error_handler(E_ERR, 'state_regress_relativefrac', trim(errstring))
                end if
            end if
        end do

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

    ! Subroutine to get delta representation (Boden & Grooms, 2025) from state
    subroutine get_delta_representation_from_state(ens_prior_a, ens_prior_v, H_l, H_u, delta_representation, nc, ens_size)

        ! Declare variables
        integer,             intent(in)  :: nc, ens_size
        integer                          :: n
        real(kind=dbl_kind), intent(in)  :: ens_prior_a(ens_size, nc), ens_prior_v(ens_size, nc)
        real(kind=dbl_kind), intent(out) :: delta_representation(ens_size, nc*2),                &
                                            H_l(ens_size, nc), H_u(ens_size, nc)
        real(kind=dbl_kind)              :: h_midpoints(ens_size, nc), deltaH(ens_size, nc),     &
                                            x_r(ens_size, nc), x_l(ens_size, nc),                &
                                            mask_a_zero(ens_size, nc),
        real(kind=dbl_kind)              :: H_bnds(nc)
        real(kind=dbl_kind)              :: cc1, cc2, cc3, x1, rncat


        ! Step 1: calculate the midpoint thickness from the prior
        mask_a_zero = (ens_prior_a > 0)
        h_midpoints = ens_prior_v / ens_prior_a
        h_midpoints(.not. mask_a_zero) = 0

        ! Step 2: calculate the ITD bounds
        get_ITD_bounds(h_midpoints, H_l, H_u, nc, ens_size)
        deltaH = H_u - H_l

        ! Step 4: calculate the delta representation components for each category and ensemble member
        x_r = ens_prior_a * (H_u - h_midpoints)/deltaH
        x_l = ens_prior_a * (h_midpoints - H_l)/deltaH
        x_r(.not. mask_a_zero) = 0
        x_l(.not. mask_a_zero) = 0

        ! Step 5. interleave the left and right components into the delta representation array
        do n = 1, nc
            delta_representation(:, 2*n-1) = x_l(:, n)
            delta_representation(:, 2*n)   = x_r(:, n)
        end do

        ! Step 6. verify that the sum of the delta representation is the same as the sum of the original category fractions
        ! in each ensemble member
        if (any(abs(sum(delta_representation, dim=2) - sum(ens_prior_a, dim=2)) > 1e-12)) then
            print *, "Warning: Delta representation verification failed."
        endif

    end subroutine get_delta_representation_from_state
    
    ! Subroutine to get state back from delta representation
    subroutine get_state_from_delta_representation(ens_post_a, ens_post_v, H_l, H_u, delta_representation, nc, ens_size)

        ! Declare variables
        integer,             intent(in)  :: nc, ens_size
        integer                          :: n
        real(kind=dbl_kind), intent(in)  :: delta_representation(ens_size, nc*2), &
                                            H_l(ens_size, nc), H_u(ens_size, nc)
        real(kind=dbl_kind), intent(out) :: ens_post_a(ens_size, nc), ens_post_v(ens_size, nc)

        ! Step 1: extract the left and right components from the delta representation array
        ! Note that the v posterior calculation assumes that the upper bound on the thickest category stays the same
        ! as in the prior. 
        do n = 1, nc
            ens_post_a(:, n) = delta_representation(:, 2*n-1) + delta_representation(:, 2*n)
            ens_post_v(:, n) = delta_representation(:, 2*n-1) * H_l(:, n) + delta_representation(:, 2*n) * H_u(:, n)
        end do

    end subroutine get_state_from_delta_representation

    ! Subroutine to get the category bounds based on nc
    subroutine get_ITD_bounds(h_midpoints, H_l, H_u, nc, ens_size)

        ! Declare variables
        integer,             intent(in)  :: n, nc, ens_size
        real(kind=dbl_kind), intent(in)  :: h_midpoints(ens_size, nc)
        real(kind=dbl_kind), intent(out) :: H_l(ens_size, nc), H_u(ens_size, nc)
        real(kind=dbl_kind)              :: H_bnds(nc)
        real(kind=dbl_kind)              :: cc1, cc2, cc3, x1, rncat
        
        ! Calculate category bounds
        rncat = real(nc, kind=dbl_kind)
        cc1 = 3.0_dbl_kind / rncat
        cc2 = 15.0_dbl_kind * cc1
        cc3 = 3.0_dbl_kind
        H_bnds(1) = 0.0_dbl_kind
        do n = 2, nc
            x1 = real(n-1, kind=dbl_kind) / rncat
            H_bnds(n) = H_bnds(n-1) + cc1 + cc2 * (1.0_dbl_kind + tanh(cc3*(x1-1.0_dbl_kind)))
        end do 

        ! Calculate the ensemble array of upper and lower bounds
        H_l(:, :) = H_bnds
        H_u(:, 1:nc-1) = H_bnds(2:nc)
        H_u(:, nc) = 3 * h_midpoints(:, nc) - 2*H_bnds(nc)

    end subroutine get_ITD_bounds

end module fractional_reg_mod
