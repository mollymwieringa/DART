! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_cice

!----------------------------------------------------------------------
! purpose: implement a 'partition function' to modify the cice state 
!          to be consistent with the states from assimilation
!
! method: Read in restart (restart with prior) and out restart (restart 
!         with posterior) written by DART after filter. 
!
! author: M Wieringa (2023) based on C Bitz (2016) and C Riedel (2022)
!----------------------------------------------------------------------

use  types_mod, only : r8
use  utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             file_exist, error_handler, E_ERR, E_MSG, to_upper
use ice_postprocessing_mod, only : cice_rebalancing, area_simple_squeeze, &
                                   volume_simple_squeeze, get_3d_variable, &
                                   write_3d_variable
use  netcdf_utilities_mod, only : nc_check
use  netcdf 

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! SET NAMELIST AND ALLOCATE VARIABLES    
!------------------------------------------------------------------
character(len=256) :: dart_to_cice_input_file = 'post_filter_restart.nc'
character(len=256) :: original_cice_restart_file = 'pre_filter_restart.nc'
character(len=256) :: postprocessed_output_file = 'postprocessed_restart.nc'
character(len=128) :: balance_method = 'simple_squeeze'
character(len=128) :: postprocess = 'cice'

namelist /dart_to_cice_nml/ dart_to_cice_input_file,    &
                            original_cice_restart_file, &
                            postprocessed_output_file,  &
                            balance_method,             &
                            postprocess

! general variable iniatlization
character(len=512) :: string1, string2, msgstring
character(len=128) :: method
character(len=3)   :: nchar

integer :: iunit, io, ncid, dimid, l, n, i, j, Ncat, Nx, Ny
real(r8), allocatable :: aicen_original(:,:,:), vicen_original(:,:,:), vsnon_original(:,:,:)
real(r8), allocatable :: aicen(:,:,:), vicen(:,:,:), vsnon(:,:,:), Tsfcn(:,:,:)
real(r8), allocatable :: qice001(:,:,:), qice002(:,:,:), qice003(:,:,:), qice004(:,:,:), qice005(:,:,:), qice006(:,:,:), qice007(:,:,:), qice008(:,:,:)
real(r8), allocatable :: sice001(:,:,:), sice002(:,:,:), sice003(:,:,:), sice004(:,:,:), sice005(:,:,:), sice006(:,:,:), sice007(:,:,:), sice008(:,:,:)
real(r8), allocatable :: qsno001(:,:,:), qsno002(:,:,:), qsno003(:,:,:)

!------------------------------------------------------------------
! INIALIZE AND PERFORM CHECKS ON FILES                
!------------------------------------------------------------------
call initialize_utilities(progname='dart_to_cice')

call find_namelist_in_file("input.nml", "dart_to_cice_nml", iunit)
read(iunit, nml = dart_to_cice_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cice_nml")

method = balance_method
call to_upper(method)

write(string1,*) 'converting DART output file "'// &
                 &trim(dart_to_cice_input_file)//'" to one CICE will like'
write(string2,*) 'using the "'//trim(balance_method)//'" method.'
call error_handler(E_MSG,'dart_to_cice',string1,text2=string2)

if ( .not. file_exist(dart_to_cice_input_file) ) then
   write(string1,*) 'cannot open "', trim(dart_to_cice_input_file),'" for updating.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(dart_to_cice_input_file))
endif

if ( .not. file_exist(original_cice_restart_file) ) then
   write(string1,*) 'cannot open "', trim(original_cice_restart_file),'" for reading.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(original_cice_restart_file))
endif

!------------------------------------------------------------------
! READ VARIABLES FROM RESTART FILES               
!------------------------------------------------------------------
! Read the pre-assim variables
call nc_check( nf90_open(trim(original_cice_restart_file), NF90_NOWRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(original_cice_restart_file)//'"')

! get dimension information
call nc_check(nf90_inq_dimid(ncid, "ncat", dimid), &
              'dart_to_cice', 'inquire ncat dimid from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Ncat), &
              'dart_to_cice', 'inquire ncat from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inq_dimid(ncid,"ni",dimid), &
               'dart_to_cice', 'inquire ni dimid from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inquire_dimension(ncid,dimid,len=Nx),&
               'dart_to_cice', 'inquire ni from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inq_dimid(ncid,"nj",dimid), &
               'dart_to_cice', 'inquire nj dimid from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inquire_dimension(ncid,dimid,len=Ny),&
               'dart_to_cice', 'inquire nj from "'//trim(original_cice_restart_file)//'"')

call get_3d_variable(ncid, 'aicen', aicen_original, original_cice_restart_file)
call get_3d_variable(ncid, 'vicen', vicen_original, original_cice_restart_file)
call get_3d_variable(ncid, 'vsnon', vsnon_original, original_cice_restart_file)

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(original_cice_restart_file))

! Read the post-assim variables 
call nc_check( nf90_open(trim(dart_to_cice_input_file), NF90_WRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(dart_to_cice_input_file)//'"')

! get the key restart variables post-assimilation (allocated in routine)
call get_3d_variable(ncid, 'aicen', aicen, dart_to_cice_input_file)
call get_3d_variable(ncid, 'vicen', vicen, dart_to_cice_input_file)
call get_3d_variable(ncid, 'vsnon', vsnon, dart_to_cice_input_file)
call get_3d_variable(ncid, 'Tsfcn', Tsfcn, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice001', sice001, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice002', sice002, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice003', sice003, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice004', sice004, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice005', sice005, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice006', sice006, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice007', sice007, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice008', sice008, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice001', qice001, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice002', qice002, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice003', qice003, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice004', qice004, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice005', qice005, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice006', qice006, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice007', qice007, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice008', qice008, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qsno001', qsno001, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qsno002', qsno002, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qsno003', qsno003, dart_to_cice_input_file)

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(dart_to_cice_input_file))

! write(*,*) 'aicen dimensions are ', shape(aicen)
!------------------------------------------------------------------
! PERFORM POSTPROCESSING 
!------------------------------------------------------------------
if (postprocess == 'cice') then
   write(*,*) 'calling cice postprocessing...'
   call cice_rebalancing(qice001, qice002,     &
                         qice003, qice004,     &
                         qice005, qice006,     &
                         qice007, qice008,     &
                         sice001, sice002,     &
                         sice003, sice004,     &
                         sice005, sice006,     &
                         sice007, sice008,     &
                         qsno001, qsno002,     &
                         qsno003, aicen,       &
                         vicen, vsnon,         &
                         aicen_original,       &
                         vicen_original,       &
                         vsnon_original,       &
                         Tsfcn,                &
                         Ncat, Nx, Ny)
   write(*,*) 'cice postprocessing function completed...'

else if (postprocess == 'aice') then 
   write(*,*) 'calling aice postprocessing...'
   call area_simple_squeeze(qice001, qice002,     &
                            qice003, qice004,     &
                            qice005, qice006,     &
                            qice007, qice008,     &
                            sice001, sice002,     &
                            sice003, sice004,     &
                            sice005, sice006,     &
                            sice007, sice008,     &
                            qsno001, qsno002,     &
                            qsno003, aicen,       &
                            vicen, vsnon,         &
                            aicen_original,       &
                            vicen_original,       &
                            vsnon_original,       &
                            Tsfcn,                &
                            Ncat, Nx, Ny)
   write(*,*) 'aice postprocessing function completed...'
else if (postprocess == 'vice') then
   write(*,*) 'calling vice postprocessing...'
   call volume_simple_squeeze(qice001, qice002,     &
                              qice003, qice004,     &
                              qice005, qice006,     &
                              qice007, qice008,     &
                              sice001, sice002,     &
                              sice003, sice004,     &
                              sice005, sice006,     &
                              sice007, sice008,     &
                              qsno001, qsno002,     &
                              qsno003, aicen,       &
                              vicen, vsnon,         &
                              aicen_original,       &
                              vicen_original,       &
                              vsnon_original,       &
                              Tsfcn,                &
                              Ncat, Nx, Ny)
   write(*,*) 'vice postprocessing function completed...'
else
   write(*,*) 'No valid postprocessing method called. No adjustments will be made.'
end if

!------------------------------------------------------------------
! WRITE VARIABLES TO RESTART FILE
!------------------------------------------------------------------
call nc_check( nf90_open(trim(postprocessed_output_file), NF90_WRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(postprocessed_output_file)//'"')

call write_3d_variable(ncid, 'aicen', aicen, postprocessed_output_file)
call write_3d_variable(ncid, 'vicen', vicen, postprocessed_output_file)
call write_3d_variable(ncid, 'vsnon', vsnon, postprocessed_output_file)
call write_3d_variable(ncid, 'Tsfcn', Tsfcn, postprocessed_output_file)
call write_3d_variable(ncid, 'qice001', qice001, postprocessed_output_file)
call write_3d_variable(ncid, 'qice002', qice002, postprocessed_output_file)
call write_3d_variable(ncid, 'qice003', qice003, postprocessed_output_file)
call write_3d_variable(ncid, 'qice004', qice004, postprocessed_output_file)
call write_3d_variable(ncid, 'qice005', qice005, postprocessed_output_file)     
call write_3d_variable(ncid, 'qice006', qice006, postprocessed_output_file)
call write_3d_variable(ncid, 'qice007', qice007, postprocessed_output_file)
call write_3d_variable(ncid, 'qice008', qice008, postprocessed_output_file)
call write_3d_variable(ncid, 'sice001', sice001, postprocessed_output_file)
call write_3d_variable(ncid, 'sice002', sice002, postprocessed_output_file)
call write_3d_variable(ncid, 'sice003', sice003, postprocessed_output_file)
call write_3d_variable(ncid, 'sice004', sice004, postprocessed_output_file)
call write_3d_variable(ncid, 'sice005', sice005, postprocessed_output_file)
call write_3d_variable(ncid, 'sice006', sice006, postprocessed_output_file)     
call write_3d_variable(ncid, 'sice007', sice007, postprocessed_output_file)
call write_3d_variable(ncid, 'sice008', sice008, postprocessed_output_file)
call write_3d_variable(ncid, 'qsno001', qsno001, postprocessed_output_file)
call write_3d_variable(ncid, 'qsno002', qsno002, postprocessed_output_file)
call write_3d_variable(ncid, 'qsno003', qsno003, postprocessed_output_file)

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(postprocessed_output_file))


!------------------------------------------------------------------
! DEALLOCATE AND FINALIZE                
!------------------------------------------------------------------
deallocate(aicen, vicen, vsnon, Tsfcn, aicen_original, vicen_original, vsnon_original)
deallocate(qice001, qice002, qice003, qice004, qice005, qice006, qice007, qice008)
deallocate(sice001, sice002, sice003, sice004, sice005, sice006, sice007, sice008)
deallocate(qsno001, qsno002, qsno003)

call finalize_utilities('dart_to_cice')

!------------------------------------------------------------------
contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! FUNCTIONS       
!------------------------------------------------------------------
function enthalpy_mush(zTin, zSin) result(zqin)

   ! enthalpy of mush from mush temperature and bulk salinity

   real(r8), intent(in) :: &
        zTin, & ! ice layer temperature (C)
        zSin    ! ice layer bulk salinity (ppt)

   real(r8) :: &
        zqin    ! ice layer enthalpy (J m-3) 

   real(r8) :: &
        phi     ! ice liquid fraction 

! from shr_const_mod.F90
   real(r8),parameter :: SHR_CONST_CPSW  = 3.996e3_R8   ! specific heat of sea water ~ J/kg/K
   real(R8),parameter :: SHR_CONST_CPICE = 2.11727e3_R8 ! specific heat of fresh ice ~ J/kg/K
   real(R8),parameter :: SHR_CONST_RHOSW = 1.026e3_R8   ! density of sea water ~ kg/m^3
   real(R8),parameter :: SHR_CONST_RHOICE= 0.917e3_R8   ! density of ice        ~ kg/m^3
   real(R8),parameter :: SHR_CONST_LATICE= 3.337e5_R8   ! latent heat of fusion ~ J/kg


! from cice/src/drivers/cesm/ice_constants.F90
   real(r8) :: cp_ocn, cp_ice, rhoi, rhow, Lfresh

   cp_ice    = SHR_CONST_CPICE  ! specific heat of fresh ice (J/kg/K)
   cp_ocn    = SHR_CONST_CPSW   ! specific heat of ocn    (J/kg/K)
   rhoi      = SHR_CONST_RHOICE ! density of ice (kg/m^3)
   rhow      = SHR_CONST_RHOSW  ! density of seawater (kg/m^3)
   Lfresh    = SHR_CONST_LATICE ! latent heat of melting of fresh ice (J/kg)

   phi = liquid_fraction(zTin, zSin)

   zqin = phi * (cp_ocn * rhow - cp_ice * rhoi) * zTin + &
          rhoi * cp_ice * zTin - (1._r8 - phi) * rhoi * Lfresh

 end function enthalpy_mush

 function liquid_fraction(zTin, zSin) result(phi)

   ! liquid fraction of mush from mush temperature and bulk salinity

   real(r8), intent(in) :: &
        zTin, & ! ice layer temperature (C)
        zSin    ! ice layer bulk salinity (ppt)

   real(r8) :: &
        phi , & ! liquid fraction
        Sbr     ! brine salinity (ppt)

   real (r8), parameter :: puny = 1.0e-11_r8 ! cice/src/drivers/cesm/ice_constants.F90

   Sbr = max(liquidus_brine_salinity_mush(zTin),puny)
   phi = zSin / max(Sbr, zSin)

 end function liquid_fraction

 function snow_enthaply(Ti) result(qsno)
   real(r8), intent(in) :: Ti

   real(r8),parameter :: rhos = 330.0_r8, &
                       Lfresh = 2.835e6_r8 - 2.501e6_r8, &
                       cp_ice = 2106._r8
   real(r8) :: qsno

   qsno = -rhos*(Lfresh - cp_ice*min(0.0_r8,Ti))
end function snow_enthaply

function liquidus_brine_salinity_mush(zTin) result(Sbr)

   ! liquidus relation: equilibrium brine salinity as function of temperature
   ! based on empirical data from Assur (1958)

   real(r8), intent(in) :: &
        zTin         ! ice layer temperature (C)

   real(r8) :: &
        Sbr          ! ice brine salinity (ppt)

   real(r8) :: &
        t_high   , & ! mask for high temperature liquidus region
        lsubzero     ! mask for sub-zero temperatures

   !constant numbers from ice_constants.F90
   real(r8), parameter :: &
        c1      = 1.0_r8 , &
        c1000   = 1000_r8

   ! liquidus relation - higher temperature region
   real(r8), parameter :: &
        az1_liq = -18.48_r8 ,&
        bz1_liq =   0.0_r8

   ! liquidus relation - lower temperature region
   real(r8), parameter :: &
        az2_liq = -10.3085_r8,  &
        bz2_liq =  62.4_r8

   ! liquidus break
   real(r8), parameter :: &
        Tb_liq = -7.6362968855167352_r8
        
   ! basic liquidus relation constants
   real(r8), parameter :: &
        az1p_liq = az1_liq / c1000, &
        bz1p_liq = bz1_liq / c1000, &
        az2p_liq = az2_liq / c1000, &
        bz2p_liq = bz2_liq / c1000

   ! temperature to brine salinity
   real(r8), parameter :: &
      J1_liq = bz1_liq / az1_liq         , &
      K1_liq = c1 / c1000                , &
      L1_liq = (c1 + bz1p_liq) / az1_liq , &
      J2_liq = bz2_liq  / az2_liq        , &
      K2_liq = c1 / c1000                , &
      L2_liq = (c1 + bz2p_liq) / az2_liq

   t_high   = merge(1._r8, 0._r8, (zTin > Tb_liq))
   lsubzero = merge(1._r8, 0._r8, (zTin <= 1._r8))

   Sbr = ((zTin + J1_liq) / (K1_liq * zTin + L1_liq)) * t_high + &
         ((zTin + J2_liq) / (K2_liq * zTin + L2_liq)) * (1._r8 - t_high)

   Sbr = Sbr * lsubzero

end function liquidus_brine_salinity_mush

function liquidus_temperature_mush(Sbr) result(zTin)

   ! liquidus relation: equilibrium temperature as function of brine salinity
   ! based on empirical data from Assur (1958)

   real(r8), intent(in) :: &
        Sbr    ! ice brine salinity (ppt)

   real(r8) :: &
        zTin   ! ice layer temperature (C)

   real(r8) :: &
        t_high ! mask for high temperature liquidus region

   ! liquidus break
   real(r8), parameter :: &
      Sb_liq =  123.66702800276086_r8    ! salinity of liquidus break

   ! constant numbers from ice_constants.F90
   real(r8), parameter :: &
        c1      = 1.0_r8 , &
        c1000   = 1000_r8

   ! liquidus relation - higher temperature region
   real(r8), parameter :: &
        az1_liq = -18.48_r8 ,&
        bz1_liq =   0.0_r8

   ! liquidus relation - lower temperature region
   real(r8), parameter :: &
        az2_liq = -10.3085_r8,  &
        bz2_liq =  62.4_r8

   ! basic liquidus relation constants
   real(r8), parameter :: &
        az1p_liq = az1_liq / c1000, &
        bz1p_liq = bz1_liq / c1000, &
        az2p_liq = az2_liq / c1000, &
        bz2p_liq = bz2_liq / c1000

 ! brine salinity to temperature
   real(r8), parameter :: &
      M1_liq = az1_liq            , &
      N1_liq = -az1p_liq          , &
      O1_liq = -bz1_liq / az1_liq , &
      M2_liq = az2_liq            , &
      N2_liq = -az2p_liq          , &
      O2_liq = -bz2_liq / az2_liq

   t_high = merge(1._r8, 0._r8, (Sbr <= Sb_liq))

   zTin = ((Sbr / (M1_liq + N1_liq * Sbr)) + O1_liq) * t_high + &
         ((Sbr / (M2_liq + N2_liq * Sbr)) + O2_liq) * (1._r8 - t_high)

end function liquidus_temperature_mush

!------------------------------------------------------------------
! END             
!------------------------------------------------------------------
end program dart_to_cice

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
