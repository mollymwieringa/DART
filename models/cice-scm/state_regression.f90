program state_regression

    use types_mod,               only : r8
    use utilities_mod,           only : initialize_utilities, finalize_utilities,   &
                                        check_namelist_read, find_namelist_in_file
    use distribution_params_mod, only : NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION
    use fractional_reg_mod,      only : state_regress_disaggregation, state_regress_relativefrac,   &
                                        state_regress_probit, postprocess

    ! Declare variables
    integer               :: ens_size, nc
    integer               :: iunit, io, i
    integer               :: dist_for_obs, dist_for_state
    logical               :: bounded_below, bounded_above
    character(len=128)    :: obs_dist
    real(r8)              :: start_time, end_time
    real(r8)              :: lower_bound, upper_bound
    real(r8), allocatable :: obs_prior(:), obs_post(:), obs_inc(:), ens_prior(:,:), ens_post(:,:)

    ! Handle namelist reading 
    namelist / state_regression_nml / obs_dist, bounded_below, bounded_above, &
                                           lower_bound, upper_bound, ens_size, nc
    
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
    
    ! -------------------------------------------------------------------------------------------------
    ! BEGIN 
    ! -------------------------------------------------------------------------------------------------
    ! Get the state prior ensemble information
    open(unit=23, file='prior_ensemble.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    do i = 1, ens_size
        read(23, REC=i) ens_prior(i,:)
    end do
    close(23)

    ! get the obs prior ensemble information
    open(unit=23, file='obs_prior_ensemble.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    do i = 1, ens_size
        read(23, REC=i) obs_prior(i)
    end do
    close(23)

    ! get the obs posterior information
    open(unit=23, file='obs_post_ensemble.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    do i = 1, ens_size
        read(23, REC=i) obs_post(i)
    end do
    close(23)

    ! Get the obs increments from the prior ensemble and observation information
    if (obs_dist == 'bnrh') then
        dist_for_obs = BOUNDED_NORMAL_RH_DISTRIBUTION
        dist_for_state = BOUNDED_NORMAL_RH_DISTRIBUTION
        open(unit=23, file='obs_increments_bnrh.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    else if (obs_dist == 'kde') then
        dist_for_obs = KDE_DISTRIBUTION
        dist_for_state = KDE_DISTRIBUTION
        open(unit=23, file='obs_increments_kde.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    else 
        dist_for_obs = NORMAL_DISTRIBUTION
        dist_for_state = NORMAL_DISTRIBUTION
        print *, 'Supplied observation distribution type is: ', obs_dist, '. By default, normal increments assumed.'
        print *, 'Please ensure this is correct if obs_dist is not normal in the namelist.'
        open(unit=23, file='obs_increments_norm.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    end if

    print *, 'Reading observation increments from obs adjustments with a ', obs_dist, ' distribution in the QCEF.'
    do i = 1, ens_size
        read(23, REC=i) obs_inc(i)
    end do
    close(23)

    ! -------------------------------------------------------------------------------------------------
    ! PERFORM REGRESSION METHODS
    ! -------------------------------------------------------------------------------------------------
    ! --- Disaggregation ------------------------------------------------------------------------------
    call cpu_time(start_time)
    call state_regress_disaggregation(obs_inc, ens_prior, ens_post, ens_size, nc, &
                                      bounded_above, bounded_below, upper_bound, lower_bound)
    call cpu_time(end_time)
    print *, 'Time for regression by disaggregation: ', end_time - start_time

    open(unit=23, file='ens_post_disaggregation.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    write(23, REC=1) ens_post
    close(23)

    ! --- Relative Fractional -------------------------------------------------------------------------
    call cpu_time(start_time)
    ! We need some approach here to make sure that everything sums to one in the prior
    ! and that nc goes up if necessary
    call state_regress_relativefrac(obs_prior, obs_inc, ens_prior, ens_post, ens_size, nc, &
                                    bounded_above, bounded_below, upper_bound, lower_bound)
    call cpu_time(end_time)
    print *, 'Time for regression by relative fractional amount: ', end_time - start_time

    open(unit=23, file='ens_post_relativefrac.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    write(23, REC=1) ens_post
    close(23)

    ! --- Probit + postprocessing ---------------------------------------------------------------------
    call cpu_time(start_time)
    call state_regress_probit(obs_prior, obs_post, ens_prior, ens_post, ens_size, nc, &
                              dist_for_obs, dist_for_state, &
                              bounded_above, bounded_below, upper_bound, lower_bound)
    call cpu_time(end_time)
    print *, 'Time for regression by probit (DART default): ', end_time - start_time

    open(unit=23, file='ens_post_probit_raw.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    write(23, REC=1) ens_post
    close(23)

    call cpu_time(start_time)
    call postprocess(ens_post, ens_size, nc, bounded_above, bounded_below, upper_bound, lower_bound)
    call cpu_time(end_time)
    print *, 'Time for probit postprocessing: ', end_time - start_time

    open(unit=23, file='ens_post_probit_postprocess.txt', access='DIRECT', form='UNFORMATTED', status='UNKNOWN', RECL=8*ens_size)
    write(23, REC=1) ens_post
    close(23)

    ! -------------------------------------------------------------------------------------------------
    ! END
    ! -------------------------------------------------------------------------------------------------

    call finalize_utilities

end program state_regression