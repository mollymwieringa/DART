program obs_increments
   use types_mod,             only : r8
   use utilities_mod,         only : initialize_utilities, finalize_utilities,                &
                                     check_namelist_read, find_namelist_in_file
   use netcdf_utilities_mod,  only : nc_open_file_readonly, nc_get_variable, &
                                     nc_create_file, nc_close_file, nc_define_dimension, &
                                     nc_define_double_variable, nc_add_global_attribute, &
                                     nc_end_define_mode, nc_put_variable
   use kde_distribution_mod,  only : obs_dist_types
   use random_seq_mod,        only : random_seq_type, random_gaussian, init_random_seq,       &
                                     random_uniform
   use       sort_mod,        only : index_sort
   use bnrh_distribution_mod, only : inv_bnrh_cdf, bnrh_cdf, inv_bnrh_cdf_like
   use assim_tools_mod,       only : obs_increment_bounded_norm_rhf, obs_increment_eakf

   logical  :: bounded_below, bounded_above
   real(r8) :: lower_bound, upper_bound
   integer  :: obs_dist_type
   real(r8) :: y, obs_var
   integer  :: ens_size, nc
   character(len=256) :: filename
   character(len=20)  :: state_dims(2)
   real(r8), allocatable :: ens(:), state_ens_a(:,:), state_ens_v(:,:), obs_inc(:), likelihood(:)
   real(r8) :: like_sum, prior_mean, prior_var, net_a=0._r8
   real(r8) :: start_time, end_time

   integer  :: iunit, io, i, ncid

   namelist / obs_increments_nml / bounded_below, bounded_above, &
                                   lower_bound, upper_bound, &
                                   ens_size, nc

   ! All tests involve a normal likelihood
   obs_dist_type = obs_dist_types%normal

   call initialize_utilities()

   call find_namelist_in_file("input.nml", "obs_increments_nml", iunit)
   read(iunit, nml = obs_increments_nml, iostat = io)
   call check_namelist_read(iunit, io, "obs_increments_nml")

   allocate(ens(ens_size))
   allocate(obs_inc(ens_size))
   allocate(likelihood(ens_size))
   allocate(state_ens_a(ens_size, nc))
   allocate(state_ens_v(ens_size, nc))

   ! Read prior ensemble and observation information
   open(unit=23, file='prior_ens_a.txt', status='OLD')
   do i=1, ens_size
      read(23,*) state_ens_a(i, :)
   end do
   close(23)

   open(unit=23, file='prior_ens_v.txt', status='OLD')
   do i=1, ens_size
      read(23,*) state_ens_v(i, :)
   end do
   close(23)
   
   ! calculate ens by summing state_ens over the second dimension (nc)
   ens = sum(state_ens_a, dim=2)

   open(unit=23,file='obs_info.txt', status='OLD')
   read(23, *) y
   read(23, *) obs_var
   print *, 'obs info: ', y, obs_var
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

     ! For the unbounded cases call the EAKF
   if(.not. (bounded_below .or. bounded_above)) then
      call cpu_time(start_time)
      prior_mean = sum(ens) / real(ens_size, r8)
      prior_var  = sum((ens - prior_mean)**2) / real(ens_size - 1, r8)
      call obs_increment_eakf(ens, ens_size, prior_mean, prior_var, &
      y, obs_var, obs_inc, net_a)
      call cpu_time(end_time)
      print *, 'Time for EAKF: ', end_time - start_time
   endif

   ! write out the analysis ensemble and increments to a netcdf file 
   filename = trim('obs_increments_output.nc')
   state_dims(1) = 'member'
   state_dims(2) = 'nc'
   ncid = nc_create_file(filename, 'create output file for obs increments')
   call nc_define_dimension(ncid, 'member', ens_size, 'create output file for obs increments')
   call nc_define_dimension(ncid, 'nc', nc, 'create_output file for obs increments')
   call nc_define_double_variable(ncid, 'obs_prior', 'member', 'create output file for obs increments')
   call nc_define_double_variable(ncid, 'obs_post', 'member', 'create output file for obs increments')
   call nc_define_double_variable(ncid, 'obs_incs', 'member', 'create output file for obs increments')
   call nc_define_double_variable(ncid, 'state_prior_a', state_dims, 'create output file for obs increments')
   call nc_define_double_variable(ncid, 'state_prior_v', state_dims, 'create output file for obs increments')
   call nc_end_define_mode(ncid, 'create output file for obs increments')
   call nc_put_variable(ncid, 'obs_prior', ens, 'create output file for obs increments')
   call nc_put_variable(ncid, 'obs_post', ens + obs_inc, 'create output file for obs increments')
   call nc_put_variable(ncid, 'obs_incs', obs_inc, 'create output file for obs increments')
   call nc_put_variable(ncid, 'state_prior_a', state_ens_a, 'create output file for obs increments')
   call nc_put_variable(ncid, 'state_prior_v', state_ens_v, 'create output file for obs increments')
   call nc_close_file(ncid, 'close output file for obs increments')

   call finalize_utilities()

end program obs_increments
