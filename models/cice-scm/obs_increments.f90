program obs_increments
   use types_mod,             only : r8
   use utilities_mod,         only : initialize_utilities, finalize_utilities,                &
                                     check_namelist_read, find_namelist_in_file
   use kde_distribution_mod,  only : separate_ensemble, likelihood_function, kde_cdf,         &
                                     inv_kde_cdf, obs_dist_types
   use random_seq_mod,        only : random_seq_type, random_gaussian, init_random_seq,       &
                                     random_uniform
   use       sort_mod,        only : index_sort
   use bnrh_distribution_mod, only : inv_bnrh_cdf, bnrh_cdf, inv_bnrh_cdf_like
   use assim_tools_mod,       only : obs_increment_kde, obs_increment_bounded_norm_rhf,       &
                                     obs_increment_eakf

   logical  :: bounded_below, bounded_above
   real(r8) :: lower_bound, upper_bound
   integer  :: obs_dist_type
   real(r8) :: y, obs_var
   integer  :: ens_size
   real(r8), allocatable :: ens(:), obs_inc(:), likelihood(:)
   real(r8) :: like_sum, prior_mean, prior_var, net_a=0._r8
   real(r8) :: start_time, end_time

   integer  :: iunit, io, i

   namelist / obs_increments_nml / bounded_below, bounded_above, lower_bound, upper_bound, &
                                        ens_size

   ! All tests involve a normal likelihood
   obs_dist_type = obs_dist_types%normal

   call initialize_utilities()

   call find_namelist_in_file("input.nml", "obs_increments_nml", iunit)
   read(iunit, nml = obs_increments_nml, iostat = io)
   call check_namelist_read(iunit, io, "obs_increments_nml")

   allocate(ens(ens_size))
   allocate(obs_inc(ens_size))
   allocate(likelihood(ens_size))

   ! Read prior ensemble and observation information
   open(unit=23,file='obs_prior_ensemble.txt',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8*ens_size)
   read(23,REC=1) ens
   close(23)

   open(unit=23,file='obs_info.txt',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8)
   read(23,REC=1) y
   read(23,REC=2) obs_var
   close(23)

   ! Get obs increments from kde+quadrature filter; write out analysis ensemble
   call cpu_time(start_time)
   call obs_increment_kde(ens, ens_size, y, obs_var, &
                          bounded_below, bounded_above, lower_bound, upper_bound, obs_inc)
   call cpu_time(end_time)
   print *, 'Time for KDE: ', end_time - start_time

   open(unit=23,file='obs_post_ensemble_kde.txt',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8*ens_size)
   write(23,REC=1) ens + obs_inc
   close(23)

   open(unit=23,file='obs_increments_kde.txt',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8*ens_size)
   write(23,REC=1) obs_inc
   close(23)

   ! Get obs increments from the bnrh filter; write out analysis ensemble
   call cpu_time(start_time)
   do i = 1, ens_size
      likelihood(i) = exp( -0.5_r8 * (ens(i) - y)**2 / obs_var )
   end do

   ! Normalize the likelihood here
   like_sum = sum(likelihood)
   ! If likelihood underflow, assume flat likelihood, so no increments
   if(like_sum <= 0.0_r8) then
      obs_inc = 0.0_r8
   else
      likelihood = likelihood / like_sum
      prior_mean = sum(ens) / real(ens_size, r8)
      prior_var  = sum((ens - prior_mean)**2) / real(ens_size - 1, r8)
      call obs_increment_bounded_norm_rhf(ens, likelihood, ens_size, prior_var, &
         obs_inc, bounded_below, bounded_above, lower_bound, upper_bound)
   endif
   call cpu_time(end_time)
   print *, 'Time for BNRHF: ', end_time - start_time

   open(unit=23,file='obs_post_ensemble_bnrhf.txt',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8*ens_size)
   write(23,REC=1) ens + obs_inc
   close(23)

   open(unit=23,file='obs_increments_bnrhf.txt',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8*ens_size)
   write(23,REC=1) obs_inc
   close(23)

   ! For the unbounded cases call the EAKF
   if(.not. (bounded_below .or. bounded_above)) then
      call cpu_time(start_time)
      prior_mean = sum(ens) / real(ens_size, r8)
      prior_var  = sum((ens - prior_mean)**2) / real(ens_size - 1, r8)
      call obs_increment_eakf(ens, ens_size, prior_mean, prior_var, &
      y, obs_var, obs_inc, net_a)
      call cpu_time(end_time)
      print *, 'Time for EAKF: ', end_time - start_time
      open(unit=23,file='obs_post_ensemble_norm.txt',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8*ens_size)
      write(23,REC=1) ens + obs_inc
      close(23)

      open(unit=23,file='obs_increments_norm.txt',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8*ens_size)
      write(23,REC=1) obs_inc
      close(23)
   endif

   call finalize_utilities()

end program obs_increments
