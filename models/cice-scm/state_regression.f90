program state_regression

    use types_mod,               only : r8
    use utilities_mod,           only : initialize_utilities, finalize_utilities,   &
                                        check_namelist_read, find_namelist_in_file
    use distribution_params_mod, only : NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, KDE_DISTRIBUTION
    use fractional_reg_mod,      only : state_regress_disaggregation, state_regress_relativefrac,   &
                                        state_regress_probit, postprocess

    ! Declare variables
    integer               :: ens_size, nc, nc_exp
    integer               :: iunit, io, i, j, istatus
    integer               :: dist_for_obs, dist_for_state
    logical               :: state_bounded_below, state_bounded_above, fexists
    logical               :: obs_bounded_below, obs_bounded_above
    character(len=128)    :: obs_dist, regression_type
    real(r8)              :: start_time, end_time
    real(r8)              :: state_lower_bound, state_upper_bound
    real(r8)              :: obs_lower_bound, obs_upper_bound
    real(r8), allocatable :: obs_prior(:), obs_post(:), obs_inc(:), ens_prior(:,:), ens_post(:,:), ens_post_temp(:)

    ! Handle namelist reading 
    namelist / state_regression_nml / obs_dist, obs_bounded_below, obs_bounded_above, &
                                      obs_lower_bound, obs_upper_bound, &
                                      state_bounded_below, state_bounded_above, &
                                      state_lower_bound, state_upper_bound, &
                                      regression_type, ens_size, nc
    
    call initialize_utilities()

    call find_namelist_in_file("input.nml", "state_regression_nml", iunit)
    read(iunit, nml = state_regression_nml, iostat = io)
    call check_namelist_read(iunit, io, "state_regression_nml")

    ! Allocate variables of interest
    allocate(obs_prior(ens_size))
    allocate(obs_post(ens_size))
    allocate(obs_inc(ens_size))
    allocate(ens_prior(ens_size, nc))
    allocate(ens_post(ens_size, nc))
    allocate(ens_post_temp(ens_size))
    
    ! -------------------------------------------------------------------------------------------------
    ! BEGIN 
    ! -------------------------------------------------------------------------------------------------
    ! Get the state prior ensemble information
    open(unit=23, file='prior_ensemble.txt', status='OLD')
    do i = 1, ens_size
        read(23, *) ens_prior(i,:)
    end do
    close(23)

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

    if (regression_type == 'disaggregation') then
        ! --- Disaggregation ------------------------------------------------------------------------------
        call cpu_time(start_time)
        call state_regress_disaggregation(obs_prior, obs_post, ens_prior, ens_post, ens_size, nc, &
                                          obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                          state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)
        call cpu_time(end_time)
        print *, 'Time for regression by disaggregation: ', end_time - start_time

        open(unit=23, file=trim('ens_post_disaggregation_')//trim(obs_dist)//trim('.txt'), status='UNKNOWN', RECL=256)
        do i = 1, ens_size
            write(23, *) ens_post(i,:)
        end do
        close(23)
    else if (regression_type == 'relativefrac') then 
        ! --- Relative Fractional -------------------------------------------------------------------------
        call cpu_time(start_time)
        ! We need some approach here to make sure that everything sums to one in the prior
        ! and that nc goes up if necessary
        call state_regress_relativefrac(obs_prior, obs_post, ens_prior, ens_post, ens_size, nc, &
                                        dist_for_obs, dist_for_state, &
                                        obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                        state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)
        call cpu_time(end_time)
        print *, 'Time for regression by relative fractional amount: ', end_time - start_time

        inquire(file='alpha_prior.txt', exist=fexists)
        if (fexists) then
            istatus = rename('alpha_prior.txt', trim('alpha_prior_')//trim(obs_dist)//trim('.txt'))
        endif
        istatus = rename('alpha_posterior.txt', trim('alpha_posterior_')//trim(obs_dist)//trim('.txt'))
        
        open(unit=23, file=trim('ens_post_relativefrac_')//trim(obs_dist)//trim('.txt'), status='UNKNOWN', RECL=256)
        do i = 1, ens_size
            write(23, *) ens_post(i,:)
        end do
        close(23)
    else if (regression_type == 'probit_postprocess') then
        ! --- Probit + postprocessing ---------------------------------------------------------------------
        call cpu_time(start_time)
        do j = 1, nc
            call state_regress_probit(obs_prior, obs_post, ens_prior(:,j), ens_post_temp, ens_size, 1, &
                                      dist_for_obs, dist_for_state, &
                                      obs_bounded_above, obs_bounded_below, obs_upper_bound, obs_lower_bound, &
                                      state_bounded_above, state_bounded_below, state_upper_bound, state_lower_bound)
                ! call update_from_obs_inc(obs_prior, obs_prior_mean, obs_prior_var, &
                !                          obs_inc, xhat_prior(:,j), ens_size, xhat_inc_j, &
                !                          reg_coef_j, net_a)
            ens_post(:,j) = ens_post_temp
        end do
        call cpu_time(end_time)
        print *, 'Time for regression by probit (DART default): ', end_time - start_time
        
        open(unit=23, file=trim('ens_post_probit_raw_')//trim(obs_dist)//trim('.txt'), status='UNKNOWN', RECL=256)
        do i = 1, ens_size
            write(23, *) ens_post(i, :)
        end do
        close(23)

        call cpu_time(start_time)
        call postprocess(ens_post, ens_size, nc, state_bounded_above, state_bounded_below, &
                         state_upper_bound, state_lower_bound, &
                         obs_bounded_above, obs_bounded_below, &
                         obs_upper_bound, obs_lower_bound)
        call cpu_time(end_time)
        print *, 'Time for probit postprocessing: ', end_time - start_time

        open(unit=23, file=trim('ens_post_probit_postprocess_')//trim(obs_dist)//trim('.txt'), status='UNKNOWN', RECL=256)
        do i = 1, ens_size
            write(23, *) ens_post(i, :)
        end do
        close(23)

    else if (regression_type == 'linear_postprocess') then
        ! --- Linear regression + postprocessing ----------------------------------------------------------
        call cpu_time(start_time)
        ! nc_exp = nc + 1
        do j = 1, nc
            call state_regress_probit(obs_prior, obs_post, ens_prior(:,j), ens_post_temp, ens_size, 1, &
                                      NORMAL_DISTRIBUTION, NORMAL_DISTRIBUTION, &
                                      .false., .false., obs_upper_bound, obs_lower_bound, &
                                      .false., .false., state_upper_bound, state_lower_bound)
            ens_post(:,j) = ens_post_temp
        end do
        call cpu_time(end_time)
        print *, 'Time for regression by linear (DART default): ', end_time - start_time
        
        open(unit=23, file=trim('ens_post_linear_raw_')//trim(obs_dist)//trim('.txt'), status='UNKNOWN', RECL=256)
        do i = 1, ens_size
            write(23, *) ens_post(i, :)
        end do
        close(23)

        call cpu_time(start_time)
        call postprocess(ens_post, ens_size, nc, state_bounded_above, state_bounded_below, &
                        state_upper_bound, state_lower_bound, &
                        obs_bounded_above, obs_bounded_below, &
                        obs_upper_bound, obs_lower_bound)
        call cpu_time(end_time)
        print *, 'Time for linear regression postprocessing: ', end_time - start_time

        open(unit=23, file=trim('ens_post_linear_postprocess_')//trim(obs_dist)//trim('.txt'), status='UNKNOWN', RECL=256)
        do i = 1, ens_size
            write(23, *) ens_post(i, :)
        end do
        close(23)
    else
        write(*,*) 'Regression type not recognized. Leaving this to crash...'
    end if

    ! -------------------------------------------------------------------------------------------------
    ! END
    ! -------------------------------------------------------------------------------------------------

    call finalize_utilities

end program state_regression