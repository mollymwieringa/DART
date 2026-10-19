program state_regression

    use types_mod,               only : r8
    use utilities_mod,           only : initialize_utilities, finalize_utilities,   &
                                        check_namelist_read, find_namelist_in_file
    use netcdf_utilities_mod,    only : nc_create_file, nc_close_file, nc_define_dimension, &
                                        nc_define_double_variable, nc_add_global_attribute, &
                                        nc_end_define_mode, nc_put_variable
    use distribution_params_mod, only : NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, KDE_DISTRIBUTION
    use fractional_reg_mod,      only : state_regress_disaggregation, state_regress_relativefrac,   &
                                        state_regress_probit, get_delta_representation_from_state,  &
                                        get_state_from_delta_representation, rebalance_volume, postprocess

    ! Declare variables
    integer               :: ens_size, nc, nc_exp, nc2, kcatbound
    integer               :: iunit, io, i, j, istatus, ncid
    integer               :: dist_for_obs, dist_for_state
    logical               :: state_bounded_below, state_bounded_above, fexists
    logical               :: obs_bounded_below, obs_bounded_above
    logical               :: delta_representation, regression_check
    character(len=128)    :: obs_dist, regression_type, filename
    real(r8)              :: start_time, end_time, inflate
    real(r8)              :: state_lower_bound, state_upper_bound
    real(r8)              :: obs_lower_bound, obs_upper_bound
    real(r8), allocatable :: obs_prior(:), obs_post(:), obs_inc(:), &
                             ens_prior_a(:,:), ens_prior_v(:,:), ens_prior_dr(:,:), &
                             ens_post_a(:,:), ens_post_v(:,:), ens_post_dr(:,:), &
                             a_post(:), xhat_prior(:,:), xhat_post(:,:), &
                             H_l(:,:), H_u(:,:), h_midpoints(:,:)
    real(r8), allocatable :: asprd(:), fsprd(:), mean(:), inflate_cat(:)

    ! Handle namelist reading 
    namelist / state_regression_nml / obs_dist, obs_bounded_below, obs_bounded_above, &
                                      obs_lower_bound, obs_upper_bound, &
                                      state_bounded_below, state_bounded_above, &
                                      state_lower_bound, state_upper_bound, &
                                      regression_type, regression_check, &
                                      delta_representation, ens_size, inflate, nc, &
                                      kcatbound
    
    call initialize_utilities()

    call find_namelist_in_file("input.nml", "state_regression_nml", iunit)
    read(iunit, nml = state_regression_nml, iostat = io)
    call check_namelist_read(iunit, io, "state_regression_nml")

    ! Allocate variables of interest
    nc2 = nc * 2
    allocate(obs_prior(ens_size))
    allocate(obs_post(ens_size))
    allocate(obs_inc(ens_size))
    allocate(ens_prior_a(ens_size, nc))
    allocate(ens_post_a(ens_size, nc))
    allocate(ens_prior_v(ens_size, nc))
    allocate(ens_post_v(ens_size, nc))
    allocate(ens_prior_dr(ens_size, nc2))
    allocate(ens_post_dr(ens_size, nc2))
    allocate(xhat_prior(ens_size, nc))
    allocate(xhat_post(ens_size, nc))
    allocate(H_l(ens_size, nc))
    allocate(H_u(ens_size, nc))
    allocate(h_midpoints(ens_size, nc))
    allocate(a_post(ens_size))
    allocate(asprd(nc))
    allocate(fsprd(nc))
    allocate(mean(nc))
    allocate(inflate_cat(nc))
    
    ! -------------------------------------------------------------------------------------------------
    ! BEGIN 
    ! -------------------------------------------------------------------------------------------------
    ! Get the state prior ensemble information
    open(unit=23, file='prior_ensemble_area.txt', status='OLD')
    do i = 1, ens_size
        read(23, *) ens_prior_a(i,:)
    end do
    close(23)

    open(unit=23, file='prior_ensemble_vol.txt', status='OLD')
    do i = 1, ens_size
        read(23, *) ens_prior_v(i,:)
    end do
    close(23)

    ! calculate the ens_prior spread across ensemble member for each state variable
    do j = 1, nc
        mean(j) = sum(ens_prior_a(:,j)) / real(ens_size, r8)
        fsprd(j) = sqrt(sum((ens_prior_a(:,j) - mean(j))**2) / real(ens_size-1, r8))
    end do

    ! get the obs prior ensemble information
    open(unit=23, file='obs_prior_ensemble.txt', status='OLD')
    do i = 1, ens_size
        read(23, *) obs_prior(i)
    end do
    close(23)

    ! Get the obs increments from the prior ensemble and observation information
    if (obs_dist == 'bnrh') then
        dist_for_obs = BOUNDED_NORMAL_RH_DISTRIBUTION
        dist_for_state = BOUNDED_NORMAL_RH_DISTRIBUTION
        ! get the obs posterior information
        open(unit=23, file='obs_post_ensemble_bnrhf.txt', status='OLD')
        do i = 1, ens_size
            read(23, *) obs_post(i)
        end do
        close(23)

        open(unit=23, file='obs_increments_bnrhf.txt', status='OLD')
    else if (obs_dist == 'kde') then
        dist_for_obs = KDE_DISTRIBUTION
        dist_for_state = KDE_DISTRIBUTION
        ! get the obs posterior information
        open(unit=23, file='obs_post_ensemble_kde.txt', status='OLD')
        do i = 1, ens_size
            read(23, *) obs_post(i)
        end do
        close(23)
        open(unit=23, file='obs_increments_kde.txt', status='OLD')
    else 
        dist_for_obs = NORMAL_DISTRIBUTION
        dist_for_state = NORMAL_DISTRIBUTION
        ! get the obs posterior information
        open(unit=23, file='obs_post_ensemble_norm.txt', status='OLD')
        do i = 1, ens_size
            read(23, *) obs_post(i)
        end do
        close(23)
        print *, 'Supplied observation distribution type is: ', obs_dist, '. By default, normal increments assumed.'
        print *, 'Please ensure this is correct if obs_dist is not normal in the namelist.'
        open(unit=23, file='obs_increments_norm.txt', status='OLD')
    end if

    print *, 'Reading observation increments from obs adjustments with a ', obs_dist, ' distribution in the QCEF.'
    do i = 1, ens_size
        read(23, *) obs_inc(i)
    end do
    close(23)

    ! -------------------------------------------------------------------------------------------------
    ! PERFORM REGRESSION METHODS
    ! -------------------------------------------------------------------------------------------------
    ! Initialize netcdf file for output
    filename = trim('state_regression_output.nc')
    ncid = nc_create_file(filename, 'state_regression')
    call nc_define_dimension(ncid, 'member', ens_size, 'state_regression')
    call nc_define_dimension(ncid, 'nc', nc, 'state_regression')
    call nc_add_global_attribute(ncid, 'regression_type', trim(regression_type), 'state_regression')
    call nc_add_global_attribute(ncid, 'obs_dist', trim(obs_dist), 'state_regression')
    call nc_define_double_variable(ncid, 'state_prior_area', 'member', 'nc', 'state_regression')
    call nc_define_double_variable(ncid, 'state_prior_volume', 'member', 'nc', 'state_regression')
    call nc_put_variable(ncid, 'state_prior_area', ens_prior_a, 'state_regression')
    call nc_put_variable(ncid, 'state_prior_volume', ens_prior_v, 'state_regression')

    if (regression_type == 'disaggregation') then
        ! --- Disaggregation ------------------------------------------------------------------------------
        call cpu_time(start_time)
        call state_regress_disaggregation(obs_prior, obs_post, ens_prior_a, ens_post_a, ens_size, nc, &
                                          obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                          state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)
        call rebalance_volume(ens_prior_a, ens_prior_v, ens_post_a, ens_post_v, nc, ens_size, kcatbound)
        call cpu_time(end_time)
        print *, 'Time for regression by disaggregation: ', end_time - start_time

    else if (regression_type == 'relativefrac') then 
        ! --- Relative Fractional -------------------------------------------------------------------------
        call cpu_time(start_time)

        if (delta_representation) then
            call nc_define_dimension(ncid, 'nc2', nc2, 'state_regression')
            call nc_define_double_variable(ncid, 'state_prior_dr', 'member', 'nc2', 'state_regression')
            call nc_define_double_variable(ncid, 'state_post_dr', 'member', 'nc2', 'state_regression')
            call nc_define_double_variable(ncid, 'H_l', 'member', 'nc', 'state_regression')
            call nc_define_double_variable(ncid, 'H_u', 'member', 'nc', 'state_regression')
            
            ! Get the delta_representation from the state
            call get_delta_representation_from_state(ens_prior_a, ens_prior_v, H_l, H_u, ens_prior_dr, &
                                                     h_midpoints, nc, ens_size, kcatbound)
            call state_regress_relativefrac(obs_prior, obs_post, ens_prior_dr, ens_post_dr, ens_size, nc2, &
                                            a_post, xhat_prior, xhat_post, dist_for_obs, dist_for_state, &
                                            obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                            state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound, &
                                            regression_check)
            call get_state_from_delta_representation(ens_post_a, ens_post_v, H_l, H_u, ens_post_dr, nc, ens_size)

            ! write the delta representation to the netcdf file
            call nc_put_variable(ncid, 'state_prior_dr', ens_prior_dr, 'state_regression')
            call nc_put_variable(ncid, 'state_post_dr', ens_post_dr, 'state_regression')
            call nc_put_variable(ncid, 'H_l', H_l, 'state_regression')
            call nc_put_variable(ncid, 'H_u', H_u, 'state_regression')
        else
            call state_regress_relativefrac(obs_prior, obs_post, ens_prior_a, ens_post_a, ens_size, nc, &
                                            a_post, xhat_prior, xhat_post, dist_for_obs, dist_for_state, &
                                            obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                            state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound, &
                                            regression_check)
            call rebalance_volume(ens_prior_a, ens_prior_v, ens_post_a, ens_post_v, nc, ens_size, kcatbound)
        end if 
        call nc_define_double_variable(ncid, 'alpha', 'member', 'state_regression')
        call nc_put_variable(ncid, 'alpha', a_post, 'state_regression')
        call cpu_time(end_time)
        print *, 'Time for regression by relative fractional amount: ', end_time - start_time

    else if (regression_type == 'probit_postprocess') then
        ! --- Probit + postprocessing ---------------------------------------------------------------------
        call cpu_time(start_time)
        call state_regress_probit(obs_prior, obs_post, ens_prior_a, ens_post_a, ens_size, nc, &
                                  dist_for_obs, dist_for_state, obs_bounded_above, obs_bounded_below, &
                                  obs_upper_bound, obs_lower_bound, state_bounded_above, & 
                                  state_bounded_below, state_upper_bound, state_lower_bound)
        call cpu_time(end_time)
        print *, 'Time for regression by probit (DART default): ', end_time - start_time
        
        call nc_define_double_variable(ncid, 'state_post_a_raw', 'member', 'nc', 'state_regression')
        call nc_put_variable(ncid, 'state_post_a_raw', ens_post_a, 'state_regression')

        call cpu_time(start_time)
        call postprocess(ens_post_a, ens_size, nc, state_bounded_above, state_bounded_below, &
                         state_upper_bound, state_lower_bound, &
                         obs_bounded_above, obs_bounded_below, &
                         obs_upper_bound, obs_lower_bound)
        call rebalance_volume(ens_prior_a, ens_prior_v, ens_post_a, ens_post_v, nc, ens_size, kcatbound)
        call cpu_time(end_time)
        print *, 'Time for probit postprocessing: ', end_time - start_time

    else if (regression_type == 'linear_postprocess') then
        ! --- Linear regression + postprocessing ----------------------------------------------------------
        call cpu_time(start_time)
        call state_regress_probit(obs_prior, obs_post, ens_prior_a, ens_post_a, ens_size, nc, &
                                  NORMAL_DISTRIBUTION, NORMAL_DISTRIBUTION, &
                                  .false., .false., obs_upper_bound, obs_lower_bound, &
                                  .false., .false., state_upper_bound, state_lower_bound)
        call cpu_time(end_time)
        print *, 'Time for regression by linear (DART default): ', end_time - start_time
        
        call nc_define_double_variable(ncid, 'state_post_a_raw', 'member', 'nc', 'state_regression')
        call nc_put_variable(ncid, 'state_post_a_raw', ens_post_a, 'state_regression')

        call cpu_time(start_time)
        call postprocess(ens_post_a, ens_size, nc, state_bounded_above, state_bounded_below, &
                         state_upper_bound, state_lower_bound, &
                         obs_bounded_above, obs_bounded_below, &
                         obs_upper_bound, obs_lower_bound)
        call rebalance_volume(ens_prior_a, ens_prior_v, ens_post_a, ens_post_v, nc, ens_size, kcatbound)
        call cpu_time(end_time)
        print *, 'Time for linear regression postprocessing: ', end_time - start_time
    else
        write(*,*) 'Regression type not recognized. Leaving this to crash...'
    end if

    !--------------------------------------------------------------------------------------------------
    ! Do RTPS posterior inflation
    !--------------------------------------------------------------------------------------------------
    do j = 1, nc
        mean(j) = sum(ens_post_a(:,j)) / real(ens_size, r8)
        asprd(j) = sqrt(sum((ens_post_a(:,j) - mean(j))**2) / real(ens_size-1, r8))
    end do

    do j = 1, nc
        if (asprd(j) .gt. 0.0_r8 .and. fsprd(j) .gt. 0.0_r8) then
            inflate_cat(j) = 1.0_r8 + inflate * ((fsprd(j)-asprd(j)) / asprd(j))
        else
            inflate_cat(j) = 1.0_r8
        end if
        ens_post_a(:,j) = ens_post_a(:,j) * inflate_cat(j) + mean(j) * (1.0_r8 - inflate_cat(j))
    end do

    !--------------------------------------------------------------------------------------------------
    ! Write the final posterior ensemble and inflation factors to file
    !--------------------------------------------------------------------------------------------------
    call nc_define_double_variable(ncid, 'state_post_a', 'member', 'nc', 'state_regression')
    call nc_define_double_variable(ncid, 'state_post_v', 'member', 'nc', 'state_regression')
    call nc_define_double_variable(ncid, 'inflation', 'nc', 'state_regression')
    call nc_put_variable(ncid, 'state_post_a', ens_post_a, 'state_regression')
    call nc_put_variable(ncid, 'state_post_v', ens_post_v, 'state_regression')
    call nc_put_variable(ncid, 'inflation', inflate_cat, 'state_regression')
    call nc_end_define_mode(ncid, 'state_regression')
    call nc_close_file(ncid, 'state_regression')

    ! -------------------------------------------------------------------------------------------------
    ! END
    ! -------------------------------------------------------------------------------------------------

    call finalize_utilities

end program state_regression
