! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module ice_postprocessing_mod

use types_mod, only : r8
use netcdf

implicit none

! general variable iniatlization
character(len=512) :: string1, string2, msgstring
character(len=15)  :: varname



! -----------------------------------------------------------------------------
contains
! -----------------------------------------------------------------------------

subroutine get_3d_variable(ncid, varname, var, filename)

    integer,               intent(in)  :: ncid
    character(len=*),      intent(in)  :: varname
    real(r8), allocatable, intent(out) :: var(:,:,:)
    character(len=*),      intent(in)  :: filename
 
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimLengths
    integer                               :: ndims
    character(len=NF90_MAX_NAME)          :: dimName
    
    write(msgstring,*) trim(varname)//' '//trim(filename)
    
    io = nf90_inq_varid(ncid, trim(varname), VarID)
    call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(msgstring))
    
    io = nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ndims)
    call nc_check(io, 'dart_to_cice', 'inquire_variable '//trim(msgstring))
    
    if (ndims /= 3) then
       write(string2,*) 'expected 3 dimension, got ', ndims
       call error_handler(E_ERR,'dart_to_cice',msgstring,text2=string2)
    endif
    
    dimLengths = 1
    DimensionLoop : do i = 1,ndims
    
       write(string1,'(''inquire dimension'',i2,A)') i,trim(msgstring)
       io = nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimLengths(i))
       call nc_check(io, 'dart_to_cice', string1)
    
    enddo DimensionLoop
    
    allocate( var(dimLengths(1), dimLengths(2), dimLengths(3)) )
    
    call nc_check(nf90_get_var(ncid, VarID, var), 'dart_to_cice', &
             'get_var '//trim(msgstring))
    
 end subroutine get_3d_variable
! ----------------------------------------------------------------------------- 

 subroutine area_simple_squeeze(qice001, qice002,     &
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
 
    real(r8), allocatable, intent(inout), dimension(:,:,:) :: &
                                        aicen, vicen, vsnon,  &
                                        qice001, qice002, qice003, qice004, &
                                        qice005, qice006, qice007, qice008, &
                                        sice001, sice002, sice003, sice004, &
                                        sice005, sice006, sice007, sice008, &
                                        qsno001, qsno002, qsno003,   &
                                        Tsfcn
    real(r8), allocatable, intent(in), dimension(:,:,:) :: &
                                        aicen_original,    &
                                        vicen_original,    &
                                        vsnon_original  
    integer, intent(in) :: Ncat, Nx, Ny
 
    real(r8), allocatable, dimension(:,:,:) :: hicen_original, hsnon_original, aicen_temp
    real(r8), allocatable, dimension(:,:) :: aice, aice_temp
    real(r8), allocatable, dimension(:) :: hin_max, hcat_midpoint
    real(r8) :: squeeze, cc1, cc2, x1, Si0new, Ti, qsno_hold, qi0new
    real(r8), parameter :: Tsmelt = 0._r8,        &
                           cc3 = 3._r8,           &
                           c1 = 1._r8,            &
                           phi_init = 0.75_r8,    &
                           dSin0_frazil = 3.0_r8, &
                           sss = 34.7_r8
 
    ! calculate dependent variables 
    cc1 = cc3/real(Ncat,kind=r8)
    cc2 = 15.0_r8*cc1
    Si0new = sss - dSin0_frazil              
     
    ! calculate bounds and midpoints of thickness distribution
    allocate(hin_max(Ncat))
    allocate(hcat_midpoint(Ncat))
    hin_max(0) = 0.0_r8
    do n = 1, Ncat
       x1 = real(n-1,kind=r8) / real(Ncat,kind=r8)
       hin_max(n) = hin_max(n-1) &
                     + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
       hcat_midpoint(n)=0.5_r8*(hin_max(n-1)+hin_max(n))
    enddo
    
    ! allocate(sice001(Nx,Ny,Ncat),sice002(Nx,Ny,Ncat),sice003(Nx,Ny,Ncat),sice004(Nx,Ny,Ncat), &
    !          sice005(Nx,Ny,Ncat),sice006(Nx,Ny,Ncat),sice007(Nx,Ny,Ncat),sice008(Nx,Ny,Ncat), &
    !          qice001(Nx,Ny,Ncat),qice002(Nx,Ny,Ncat),qice003(Nx,Ny,Ncat),qice004(Nx,Ny,Ncat), &
    !          qice005(Nx,Ny,Ncat),qice006(Nx,Ny,Ncat),qice007(Nx,Ny,Ncat),qice008(Nx,Ny,Ncat), &
    !          qsno001(Nx,Ny,Ncat),qsno002(Nx,Ny,Ncat),qsno003(Nx,Ny,Ncat))
    ! allocate(aicen(Nx,Ny,Ncat),vicen(Nx,Ny,Ncat),vsnon(Nx,Ny,Ncat),Tsfcn(Nx,Ny,Ncat))
    ! Begin process 
    sice001  = max(0.0_r8, sice001)  ! salinities must be non-negative
    sice002  = max(0.0_r8, sice002)  ! salinities must be non-negative
    sice003  = max(0.0_r8, sice003)  ! salinities must be non-negative
    sice004  = max(0.0_r8, sice004)  ! salinities must be non-negative
    sice005  = max(0.0_r8, sice005)  ! salinities must be non-negative
    sice006  = max(0.0_r8, sice006)  ! salinities must be non-negative
    sice007  = max(0.0_r8, sice007)  ! salinities must be non-negative
    sice008  = max(0.0_r8, sice008)  ! salinities must be non-negative
    qice001  = min(0.0_r8, qice001)  ! enthalpies (ice) must be non-positive
    qice002  = min(0.0_r8, qice002)  ! enthalpies (ice) must be non-positive
    qice003  = min(0.0_r8, qice003)  ! enthalpies (ice) must be non-positive
    qice004  = min(0.0_r8, qice004)  ! enthalpies (ice) must be non-positive
    qice005  = min(0.0_r8, qice005)  ! enthalpies (ice) must be non-positive
    qice006  = min(0.0_r8, qice006)  ! enthalpies (ice) must be non-positive
    qice007  = min(0.0_r8, qice007)  ! enthalpies (ice) must be non-positive
    qice008  = min(0.0_r8, qice008)  ! enthalpies (ice) must be non-positive
    qsno001  = min(0.0_r8, qsno001)  ! enthalphies (snow) must be non-positive
    qsno002  = min(0.0_r8, qsno002)  ! enthalphies (snow) must be non-positive
    qsno003  = min(0.0_r8, qsno003)  ! enthalphies (snow) must be non-positive
    aicen = min(1.0_r8,aicen)  ! concentrations must not exceed 1 
    Tsfcn = min(Tsmelt,Tsfcn)  ! ice/snow surface must not exceed melting
 
    ! calculate aice, which might be negative or >1 at this point
    allocate(aice(Nx,Ny))
    aice = aicen(:,:,1)
    do n = 2, Ncat  
       aice = aice+aicen(:,:,n)
    enddo
    
    ! set negative aicen to zero
    aicen = max(0.0_r8,aicen)   ! concentrations must be non-negative
    vicen = max(0.0_r8,vicen)   ! same for volumes (ice)
    vsnon = max(0.0_r8,vsnon)   ! same for volumes (snow)
 
    ! reclaculate aice, now it should be non-negative\
    allocate(aice_temp(Nx,Ny))
    aice_temp = aicen(:,:,1)
    do n = 2, Ncat
       aice_temp = aice_temp + aicen(:,:,n)
    enddo
  
    ! if aice <0, then set every category to 0
    do j = 1, Ny
       do i = 1, Nx
          if (aice(i,j)<0._r8) then
             aicen(i,j,:) = 0._r8
          endif
       enddo
    enddo
 
    ! shift negative concentration values 
    do n=1, Ncat
       do j=1, Ny
          do i=1, Nx
             if (aice_temp(i,j) > 0._r8 .and. aice(i,j)>0._r8) then
                aicen(i,j,n) = aicen(i,j,n) - (aice_temp(i,j)-aice(i,j))*aicen(i,j,n)/aice_temp(i,j)
             endif  
          enddo
       enddo
    enddo
 
    ! now squeeze aicen 
     do j = 1, Ny
       do i = 1, Nx
          if (aice(i,j) > 1.0_r8) then
             squeeze        = 1.0_r8 / aice(i,j)
             aicen(i,j,:)   = aicen(i,j,:)*squeeze
          endif
       enddo
    enddo
 
    ! update vsnon and vicen using conserved category thickness values
    allocate(aicen_temp(Nx,Ny,Ncat))
    aicen_temp = aicen_original
    where(aicen_temp==0) aicen_temp = -999
 
    allocate(hicen_original(Nx,Ny,Ncat))
    allocate(hsnon_original(Nx,Ny,Ncat))
    do n=1,Ncat
       do j=1,Ny
          do i=1,Nx
             hicen_original(i,j,n) = vicen_original(i,j,n)/aicen_temp(i,j,n)
             hsnon_original(i,j,n) = vsnon_original(i,j,n)/aicen_temp(i,j,n)
          end do
       end do
    end do
 
    where(hicen_original < 0)  hicen_original = 0.0_r8
    where(hsnon_original < 0)  hsnon_original = 0.0_r8
 
    vicen  = aicen*hicen_original
    vsnon  = aicen*hsnon_original
 
   ! consider special cases
    do n = 1, Ncat
       do j = 1, Ny
          do i = 1, Nx
          ! If there is no ice post-adjustment... 
             if (aicen(i,j,n)==0._r8 ) then
                vicen(i,j,n)   = 0._r8
                vsnon(i,j,n) = 0._r8
                sice001(i,j,n) = 0._r8
                sice002(i,j,n) = 0._r8
                sice003(i,j,n) = 0._r8
                sice004(i,j,n) = 0._r8
                sice005(i,j,n) = 0._r8
                sice006(i,j,n) = 0._r8
                sice007(i,j,n) = 0._r8
                sice008(i,j,n) = 0._r8
                qice001(i,j,n) = 0._r8
                qice002(i,j,n) = 0._r8
                qice003(i,j,n) = 0._r8
                qice004(i,j,n) = 0._r8
                qice005(i,j,n) = 0._r8
                qice006(i,j,n) = 0._r8
                qice007(i,j,n) = 0._r8
                qice008(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
                Tsfcn(i,j,n)   = -1.836_r8
             ! If the adjustment introduced new ice.. 
             else if (aicen(i,j,n)>0._r8 .and. aicen_original(i,j,n)==0._r8) then
                ! allow no snow volume or enthalpy
                vsnon(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
 
                ! require ice volume for thickness = category boundary midpoint
                vicen(i,j,n) =  aicen(i,j,n) * hcat_midpoint(n)
 
                ! salinity of mushy ice, see add_new_ice in ice_therm_itd.F90
                Si0new = sss - dSin0_frazil ! given our choice of sss
                sice001(i,j,n) = Si0new
                sice002(i,j,n) = Si0new
                sice003(i,j,n) = Si0new
                sice004(i,j,n) = Si0new
                sice005(i,j,n) = Si0new
                sice006(i,j,n) = Si0new
                sice007(i,j,n) = Si0new
                sice008(i,j,n) = Si0new
 
                ! temperature and enthalpy
                Ti          = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
                qi0new      = enthalpy_mush(Ti, Si0new)
                qice001(i,j,n) = qi0new
                qice002(i,j,n) = qi0new
                qice003(i,j,n) = qi0new
                qice004(i,j,n) = qi0new
                qice005(i,j,n) = qi0new
                qice006(i,j,n) = qi0new
                qice007(i,j,n) = qi0new
                qice008(i,j,n) = qi0new
                Tsfcn(i,j,n)  = Ti
             endif
          enddo
       enddo
    enddo
 
    deallocate(aice, aice_temp, aicen_temp, hicen_original, hsnon_original)
    deallocate(hin_max, hcat_midpoint)
 
 end subroutine area_simple_squeeze
! -----------------------------------------------------------------------------
 
 subroutine volume_simple_squeeze(qice001, qice002,     &
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
                                  aicen_original,      &
                                  vicen_original,      &
                                  vsnon_original,      &
                                  Tsfcn,               &
                                  Ncat, Nx, Ny)
 
    real(r8), allocatable, intent(inout), dimension(:,:,:) :: &
                                        aicen, vicen, vsnon,  &
                          qice001, qice002, qice003, qice004, &
                          qice005, qice006, qice007, qice008, &
                          sice001, sice002, sice003, sice004, &
                          sice005, sice006, sice007, sice008, &
                                   qsno001, qsno002, qsno003, &
                                                        Tsfcn
    real(r8), allocatable, intent(in), dimension(:,:,:) :: &
                                           aicen_original, &
                                           vicen_original, &
                                           vsnon_original   
    integer, intent(in) :: Ncat, Nx, Ny
 
    real(r8), allocatable, dimension(:,:,:) :: hicen_original, hsnon_original, aicen_temp
    real(r8), allocatable, dimension(:,:) :: aice, vice, vice_temp
    real(r8), allocatable, dimension(:) :: hin_max, hcat_midpoint
    real(r8) :: squeeze, cc1, cc2, x1, Si0new, Ti, qsno_hold, qi0new
    real(r8), parameter :: Tsmelt = 0._r8,        &
                           cc3 = 3._r8,           &
                           c1 = 1._r8,            &
                           phi_init = 0.75_r8,    &
                           dSin0_frazil = 3.0_r8, &
                           sss = 34.7_r8
 
   ! calculate dependent variables 
     cc1 = cc3/real(Ncat,kind=r8)
     cc2 = 15.0_r8*cc1
     Si0new = sss - dSin0_frazil              
 
    ! calculate bounds and midpoints of thickness distribution
    allocate(hin_max(Ncat))
    allocate(hcat_midpoint(Ncat))
    hin_max(0) = 0.0_r8
    do n = 1, Ncat
       x1 = real(n-1,kind=r8) / real(Ncat,kind=r8)
       hin_max(n) = hin_max(n-1) &
                     + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
       hcat_midpoint(n)=0.5_r8*(hin_max(n-1)+hin_max(n))
    enddo
 
    ! allocate(sice001(Nx,Ny,Ncat),sice002(Nx,Ny,Ncat),sice003(Nx,Ny,Ncat),sice004(Nx,Ny,Ncat), &
    !          sice005(Nx,Ny,Ncat),sice006(Nx,Ny,Ncat),sice007(Nx,Ny,Ncat),sice008(Nx,Ny,Ncat), &
    !          qice001(Nx,Ny,Ncat),qice002(Nx,Ny,Ncat),qice003(Nx,Ny,Ncat),qice004(Nx,Ny,Ncat), &
    !          qice005(Nx,Ny,Ncat),qice006(Nx,Ny,Ncat),qice007(Nx,Ny,Ncat),qice008(Nx,Ny,Ncat), &
    !          qsno001(Nx,Ny,Ncat),qsno002(Nx,Ny,Ncat),qsno003(Nx,Ny,Ncat))
    ! allocate(aicen(Nx,Ny,Ncat),vicen(Nx,Ny,Ncat),vsnon(Nx,Ny,Ncat),Tsfcn(Nx,Ny,Ncat))
 
    ! Begin process 
    sice001  = max(0.0_r8, sice001)  ! salinities must be non-negative
    sice002  = max(0.0_r8, sice002)  ! salinities must be non-negative
    sice003  = max(0.0_r8, sice003)  ! salinities must be non-negative
    sice004  = max(0.0_r8, sice004)  ! salinities must be non-negative
    sice005  = max(0.0_r8, sice005)  ! salinities must be non-negative
    sice006  = max(0.0_r8, sice006)  ! salinities must be non-negative
    sice007  = max(0.0_r8, sice007)  ! salinities must be non-negative
    sice008  = max(0.0_r8, sice008)  ! salinities must be non-negative
    qice001  = min(0.0_r8, qice001)  ! enthalpies (ice) must be non-positive
    qice002  = min(0.0_r8, qice002)  ! enthalpies (ice) must be non-positive
    qice003  = min(0.0_r8, qice003)  ! enthalpies (ice) must be non-positive
    qice004  = min(0.0_r8, qice004)  ! enthalpies (ice) must be non-positive
    qice005  = min(0.0_r8, qice005)  ! enthalpies (ice) must be non-positive
    qice006  = min(0.0_r8, qice006)  ! enthalpies (ice) must be non-positive
    qice007  = min(0.0_r8, qice007)  ! enthalpies (ice) must be non-positive
    qice008  = min(0.0_r8, qice008)  ! enthalpies (ice) must be non-positive
    qsno001  = min(0.0_r8, qsno001)  ! enthalphies (snow) must be non-positive
    qsno002  = min(0.0_r8, qsno002)  ! enthalphies (snow) must be non-positive
    qsno003  = min(0.0_r8, qsno003)  ! enthalphies (snow) must be non-positive    ! aicen = min(1.0_r8,aicen)  ! concentrations must not exceed 1 
    Tsfcn = min(Tsmelt,Tsfcn)  ! ice/snow surface must not exceed melting
 
    ! calculate aice, which might be negative or >1 at this point
    allocate(vice(Nx,Ny))
    vice = vicen(:,:,1)
    do n = 2, Ncat  
       vice = vice+vicen(:,:,n)
    enddo
 
    ! set negative aicen to zero
    vicen = max(0.0_r8,vicen)   ! same for volumes (ice)
 
    ! reclaculate aice, now it should be non-negative
    allocate(vice_temp(Nx,Ny))
    vice_temp = vicen(:,:,1)
    do n = 2, Ncat
       vice_temp = vice_temp + vicen(:,:,n)
    enddo
 
    ! if vice <0, then set every category to 0
    do j = 1, Ny
       do i = 1, Nx
          if (vice(i,j)<0._r8) then
             vicen(i,j,:) = 0._r8
          endif
       enddo
    enddo
 
    ! shift negative post-adjustment volume values 
    do n=1, Ncat
       do j=1, Ny
          do i=1, Nx
             if (vice_temp(i,j) > 0._r8 .and. vice(i,j)>0._r8) then
                vicen(i,j,n) = vicen(i,j,n) - (vice_temp(i,j)-vice(i,j))*vicen(i,j,n)/vice_temp(i,j)
             endif  
          enddo
       enddo
    enddo
 
    ! calculate orignal caterogy thickness values
    allocate(aicen_temp(Nx,Ny,Ncat))
    aicen_temp = aicen_original
    where(aicen_temp==0) aicen_temp = -999
 
    allocate(hicen_original(Nx,Ny,Ncat))
    allocate(hsnon_original(Nx,Ny,Ncat))
    do n=1,Ncat
       do j=1,Ny
          do i=1,Nx
             hicen_original(i,j,n) = vicen_original(i,j,n)/aicen_temp(i,j,n)
             hsnon_original(i,j,n) = vsnon_original(i,j,n)/aicen_temp(i,j,n)
          end do
       end do
    end do
 
    where(hicen_original < 0)  hicen_original = 0.0_r8
    where(hsnon_original < 0)  hsnon_original = 0.0_r8
   
    ! calculate the area implied by original category thickness and updated volume
    aicen = vicen/hicen_original
    allocate(aice(Nx,Ny))
    aice = aicen(:,:,1)
    do n = 2, Ncat  
       aice = aice+aicen(:,:,n)
    enddo
 
    ! now squeeze aicen implied by original category thickness and updated volume
    do j = 1, Ny
       do i = 1, Nx
          if (aice(i,j) > 1.0_r8) then
             squeeze = 1.0_r8/aice(i,j)
             aicen(i,j,:) = aicen(i,j,:)*squeeze
          endif
       enddo
    enddo
 
    ! recalculate volume and snow volume with squeezed vicen
    vicen = aicen*hicen_original
    vsnon = aicen*hsnon_original
 
    ! consider special cases
    do n = 1, Ncat
       do j = 1, Ny
          do i = 1, Nx
          ! If there is no ice post-adjustment... 
             if (vicen(i,j,n)==0._r8 ) then
                aicen(i,j,n)   = 0._r8
                vsnon(i,j,n)   = 0._r8
                sice001(i,j,n) = 0._r8
                sice002(i,j,n) = 0._r8
                sice003(i,j,n) = 0._r8
                sice004(i,j,n) = 0._r8
                sice005(i,j,n) = 0._r8
                sice006(i,j,n) = 0._r8
                sice007(i,j,n) = 0._r8
                sice008(i,j,n) = 0._r8
                qice001(i,j,n) = 0._r8
                qice002(i,j,n) = 0._r8
                qice003(i,j,n) = 0._r8
                qice004(i,j,n) = 0._r8
                qice005(i,j,n) = 0._r8
                qice006(i,j,n) = 0._r8
                qice007(i,j,n) = 0._r8
                qice008(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
                Tsfcn(i,j,n)   = -1.836_r8
             ! If the adjustment introduced new ice.. 
             else if (aicen(i,j,n)>0._r8 .and. aicen_original(i,j,n)==0._r8) then
                ! allow no snow volume or enthalpy
                vsnon(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
 
                ! require ice volume for thickness = category boundary midpoint
                aicen(i,j,n) =  vicen(i,j,n)/hcat_midpoint(n)
 
                ! salinity of mushy ice, see add_new_ice in ice_therm_itd.F90
                Si0new = sss - dSin0_frazil ! given our choice of sss
                sice001(i,j,n) = Si0new
                sice002(i,j,n) = Si0new
                sice003(i,j,n) = Si0new
                sice004(i,j,n) = Si0new
                sice005(i,j,n) = Si0new
                sice006(i,j,n) = Si0new
                sice007(i,j,n) = Si0new
                sice008(i,j,n) = Si0new
 
                ! temperature and enthalpy
                Ti          = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
                qi0new      = enthalpy_mush(Ti, Si0new)
                qice001(i,j,n) = qi0new
                qice002(i,j,n) = qi0new
                qice003(i,j,n) = qi0new
                qice004(i,j,n) = qi0new
                qice005(i,j,n) = qi0new
                qice006(i,j,n) = qi0new
                qice007(i,j,n) = qi0new
                qice008(i,j,n) = qi0new
                Tsfcn(i,j,n)  = Ti
             endif
          enddo
       enddo
    enddo
 
    deallocate(aice, vice, vice_temp, aicen_temp, hicen_original, hsnon_original)
    deallocate(hin_max, hcat_midpoint)
 
 end subroutine volume_simple_squeeze
! -----------------------------------------------------------------------------
 
 subroutine cice_rebalancing(qice001, qice002,     &
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
 
    real(r8), allocatable, intent(inout), dimension(:,:,:) :: &
                                        aicen, vicen, vsnon,  &
                          qice001, qice002, qice003, qice004, &
                          qice005, qice006, qice007, qice008, &
                          sice001, sice002, sice003, sice004, &
                          sice005, sice006, sice007, sice008, &
                                   qsno001, qsno002, qsno003, &
                                                       Tsfcn
    real(r8), allocatable, intent(in), dimension(:,:,:) :: &
                                           aicen_original, &
                                           vicen_original, &
                                           vsnon_original   
    integer, intent(in) :: Ncat, Nx, Ny
 
    real(r8), allocatable :: aice(:,:), vice(:,:), vsno(:,:), aice_temp(:,:), vice_temp(:,:), vsno_temp(:,:) 
    real(r8), allocatable :: hin_max(:), hcat_midpoint(:)
    real(r8) :: squeeze, cc1, cc2, x1, Si0new, Ti, qsno_hold, qi0new, hicen
    real(r8), parameter :: Tsmelt = 0._r8,        &
                           cc3 = 3._r8,           &
                           c1 = 1._r8,            &
                           phi_init = 0.75_r8,    &
                           dSin0_frazil = 3.0_r8, &
                           sss = 34.7_r8
 
    ! calculate dependent variables 
    cc1 = cc3/real(Ncat,kind=r8)
    cc2 = 15.0_r8*cc1
    Si0new = sss - dSin0_frazil              
 
    ! calculate bounds and midpoints of thickness distribution
    allocate(hin_max(Ncat))
    allocate(hcat_midpoint(Ncat))
    hin_max(0) = 0.0_r8
    do n = 1, Ncat
       x1 = real(n-1,kind=r8) / real(Ncat,kind=r8)
       hin_max(n) = hin_max(n-1) &
                   + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
       hcat_midpoint(n)=0.5_r8*(hin_max(n-1)+hin_max(n))
    enddo
 
    ! allocate(sice001(Nx,Ny,Ncat),sice002(Nx,Ny,Ncat),sice003(Nx,Ny,Ncat),sice004(Nx,Ny,Ncat), &
    !          sice005(Nx,Ny,Ncat),sice006(Nx,Ny,Ncat),sice007(Nx,Ny,Ncat),sice008(Nx,Ny,Ncat), &
    !          qice001(Nx,Ny,Ncat),qice002(Nx,Ny,Ncat),qice003(Nx,Ny,Ncat),qice004(Nx,Ny,Ncat), &
    !          qice005(Nx,Ny,Ncat),qice006(Nx,Ny,Ncat),qice007(Nx,Ny,Ncat),qice008(Nx,Ny,Ncat), &
    !          qsno001(Nx,Ny,Ncat),qsno002(Nx,Ny,Ncat),qsno003(Nx,Ny,Ncat))
    ! allocate(aicen(Nx,Ny,Ncat),vicen(Nx,Ny,Ncat),vsnon(Nx,Ny,Ncat),Tsfcn(Nx,Ny,Ncat))
 
    ! Begin process 
    write(*,*) 'beginning cice postprocessing process...'
    sice001  = max(0.0_r8, sice001)  ! salinities must be non-negative
    sice002  = max(0.0_r8, sice002)  ! salinities must be non-negative
    sice003  = max(0.0_r8, sice003)  ! salinities must be non-negative
    sice004  = max(0.0_r8, sice004)  ! salinities must be non-negative
    sice005  = max(0.0_r8, sice005)  ! salinities must be non-negative
    sice006  = max(0.0_r8, sice006)  ! salinities must be non-negative
    sice007  = max(0.0_r8, sice007)  ! salinities must be non-negative
    sice008  = max(0.0_r8, sice008)  ! salinities must be non-negative
    qice001  = min(0.0_r8, qice001)  ! enthalpies (ice) must be non-positive
    qice002  = min(0.0_r8, qice002)  ! enthalpies (ice) must be non-positive
    qice003  = min(0.0_r8, qice003)  ! enthalpies (ice) must be non-positive
    qice004  = min(0.0_r8, qice004)  ! enthalpies (ice) must be non-positive
    qice005  = min(0.0_r8, qice005)  ! enthalpies (ice) must be non-positive
    qice006  = min(0.0_r8, qice006)  ! enthalpies (ice) must be non-positive
    qice007  = min(0.0_r8, qice007)  ! enthalpies (ice) must be non-positive
    qice008  = min(0.0_r8, qice008)  ! enthalpies (ice) must be non-positive
    qsno001  = min(0.0_r8, qsno001)  ! enthalphies (snow) must be non-positive
    qsno002  = min(0.0_r8, qsno002)  ! enthalphies (snow) must be non-positive
    qsno003  = min(0.0_r8, qsno003)  ! enthalphies (snow) must be non-positive
    Tsfcn = min(Tsmelt,Tsfcn)  ! ice/snow surface must not exceed melting
    aicen = min(1.0_r8,aicen)  ! concentrations must not exceed 1 
    
    ! calculate aggregates for post-adjustment category variables 
    allocate(aice(Nx,Ny))
    allocate(vice(Nx,Ny))
    allocate(vsno(Nx,Ny))
 
    write(*,*) 'calculating aggregates...'
    aice = aicen(:,:,1)
    vice = vicen(:,:,1)
    vsno = vsnon(:,:,1)
    do n = 2, Ncat  
       aice = aice+aicen(:,:,n)
       vice = vice+vicen(:,:,n)
       vsno = vsno+vsnon(:,:,n)
    enddo
 
    ! impose bounds on categories
    write(*,*) 'imposing bounds on categories...'
    aicen = max(0.0_r8,aicen) ! concentration must be non-negative
    vicen = max(0.0_r8,vicen) ! volumes (ice) must be non-negative
    vsnon = max(0.0_r8,vsnon) ! volumes (snow) must be non-negative
 
    ! re-calculate aggregates once bounds are enforced
    allocate(aice_temp(Nx,Ny))
    allocate(vice_temp(Nx,Ny))
    allocate(vsno_temp(Nx,Ny))
 
    write(*,*) 'recalculating aggregates...'
    aice_temp = aicen(:,:,1)
    vice_temp = vicen(:,:,1)
    vsno_temp = vsnon(:,:,1)
    do n = 2, Ncat  
       aice_temp = aice_temp+aicen(:,:,n)
       vice_temp = vice_temp+vicen(:,:,n)
       vsno_temp = vsno_temp+vsnon(:,:,n)
    enddo
    
    write(*,*) 'begin squeezing...'
    do j = 1, Ny
       do i = 1, Nx
       ! if the post-adjustment concentartion was 0 or less than 0, remove all ice 
          if (aice(i,j) <= 0.0_r8) then
             aicen(i,j,:) = 0.0_r8
             vicen(i,j,:) = 0.0_r8
             vsnon(i,j,:) = 0.0_r8
          else if (aice(i,j) > 0.0_r8) then
             do n=1,Ncat
                if (aice_temp(i,j) > 0.0_r8 .and. aice(i,j) > 0.0_r8) then
                   aicen(i,j,n) = aicen(i,j,n) - (aice_temp(i,j) - aice(i,j))*aicen(i,j,n)/aice_temp(i,j) 
                endif
                if (vice_temp(i,j) > 0.0_r8 .and. vice(i,j) > 0.0_r8) then
                   vicen(i,j,n) = vicen(i,j,n) - (vice_temp(i,j) - vice(i,j))*vicen(i,j,n)/vice_temp(i,j)
                endif
                if (vsno_temp(i,j) > 0.0_r8 .and. vsno(i,j) > 0.0_r8) then
                   vsnon(i,j,n) = vsnon(i,j,n) - (vsno_temp(i,j) - vsno(i,j))*vsnon(i,j,n)/vsno_temp(i,j)
                endif
                if (aicen(i,j,n) > 0.0_r8) then
                   hicen = vicen(i,j,n)/aicen(i,j,n)
                   if (n == Ncat) then
                      if (hicen < hin_max(n-1)) then
                          aicen(i,j,n) = vicen(i,j,n)/hin_max(n-1)
                      endif
                   else
                      if (hicen > hin_max(n) .or. hicen < hin_max(n-1)) then
                         aicen(i,j,n) = vicen(i,j,n)/hcat_midpoint(n)
                      endif
                   endif
                else
                   vicen(i,j,n) = 0.0_r8
                   vsnon(i,j,n) = 0.0_r8
                endif
             enddo
 
             ! recalculate the aggregate area
             aice = aicen(:,:,1)
             do n = 2, Ncat  
                aice = aice+aicen(:,:,n)
             enddo
 
             ! If the post-adjustment concentration is greater than 1, squeeze it down
             if (aice(i,j) > 1.0_r8) then
                squeeze = 1.0_r8/aice(i,j)
                aicen(i,j,:) = aicen(i,j,:)*squeeze
             endif        
          endif
 
          !! if ice exists in both the post-adjustment and post-bounds variables, 
          !! shift the post-adjustment negative values of each category variable 
          !do n=1,Ncat
          !   if (aice_temp(i,j) > 0.0_r8 .and. aice(i,j) > 0.0_r8) then
          !      aicen(i,j,n) = aicen(i,j,n) - (aice_temp(i,j) - aice(i,j))*aicen(i,j,n)/aice_temp(i,j)
          !   endif
          !   if (vice_temp(i,j) > 0.0_r8 .and. vice(i,j) > 0.0_r8) then
          !      vicen(i,j,n) = vicen(i,j,n) - (vice_temp(i,j) - vice(i,j))*vicen(i,j,n)/vice_temp(i,j)
          !   endif
          !   if (vsno_temp(i,j) > 0.0_r8 .and. vsno(i,j) > 0.0_r8) then
          !      vsnon(i,j,n) = vsnon(i,j,n) - (vsno_temp(i,j) - vsno(i,j))*vsnon(i,j,n)/vsno_temp(i,j)
          !   endif
          !enddo
 
          !! If the post-adjustment concentration is greater than 1, squeeze it down
          !if (aice(i,j) > 1.0_r8) then
          !   squeeze = 1.0_r8/aice(i,j)
          !   aicen(i,j,:) = aicen(i,j,:)*squeeze
          !endif
 
          ! Adjust the volume, snow, salinities and enthalphies to be consistent with the squeezed concentrations
          do n=1,Ncat
             ! if the adjustment and the original category both have ice in them... 
             if (aicen(i,j,n) > 0.0_r8 .and. aicen_original(i,j,n) > 0.0_r8) then
                ! calculate the volume corresponding to the area and midpoint thickness, if there's no volume
                if (vicen(i,j,n) == 0.0_r8) vicen(i,j,n) = aicen(i,j,n)*hcat_midpoint(n)
                ! calculate the enthalphy required to accomodate any new snow in the category
                if (vsnon(i,j,n) > 0.0_r8 .and. vsnon_original(i,j,n) == 0.0_r8) then
                   Ti = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
                   qsno_hold = snow_enthaply(Ti)
                   qsno001(i,j,n) = qsno_hold
                   qsno002(i,j,n) = qsno_hold
                   qsno003(i,j,n) = qsno_hold
                endif
                ! if the adjustment doesn't have ice but the original does...
             else if (aicen(i,j,n) == 0.0_r8 .and. aicen_original(i,j,n) > 0.0_r8) then
                vicen(i,j,n) = 0.0_r8
                sice001(i,j,n) = 0._r8
                sice002(i,j,n) = 0._r8
                sice003(i,j,n) = 0._r8
                sice004(i,j,n) = 0._r8
                sice005(i,j,n) = 0._r8
                sice006(i,j,n) = 0._r8
                sice007(i,j,n) = 0._r8
                sice008(i,j,n) = 0._r8
                qice001(i,j,n) = 0._r8
                qice002(i,j,n) = 0._r8
                qice003(i,j,n) = 0._r8
                qice004(i,j,n) = 0._r8
                qice005(i,j,n) = 0._r8
                qice006(i,j,n) = 0._r8
                qice007(i,j,n) = 0._r8
                qice008(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
                vsnon(i,j,n) = 0.0_r8
                Tsfcn(i,j,n) = -1.836_r8
             ! if the adjustment has ice but the original doesn't... 
             else if (aicen(i,j,n)>0.0_r8 .and. aicen_original(i,j,n) == 0.0_r8) then
                if (vicen(i,j,n) == 0.0_r8) vicen(i,j,n) =  aicen(i,j,n) * hcat_midpoint(n)
                sice001(i,j,n) = Si0new
                sice002(i,j,n) = Si0new
                sice003(i,j,n) = Si0new
                sice004(i,j,n) = Si0new
                sice005(i,j,n) = Si0new
                sice006(i,j,n) = Si0new
                sice007(i,j,n) = Si0new
                sice008(i,j,n) = Si0new
                Ti = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
                qi0new = enthalpy_mush(Ti, Si0new)
                qice001(i,j,n) = qi0new
                qice002(i,j,n) = qi0new
                qice003(i,j,n) = qi0new
                qice004(i,j,n) = qi0new
                qice005(i,j,n) = qi0new
                qice006(i,j,n) = qi0new
                qice007(i,j,n) = qi0new
                qice008(i,j,n) = qi0new
 
                if (vsnon(i,j,n) == 0.0_r8 .and. vsnon_original(i,j,n) > 0.0_r8) then
                   qsno001(i,j,n) = 0._r8
                   qsno002(i,j,n) = 0._r8
                   qsno003(i,j,n) = 0._r8
                else if (vsnon(i,j,n) > 0.0_r8 .and. vsnon_original(i,j,n) == 0.0_r8) then
                   qsno_hold = snow_enthaply(Ti)
                   qsno001(i,j,n) = qsno_hold
                   qsno002(i,j,n) = qsno_hold
                   qsno003(i,j,n) = qsno_hold
                endif
                Tsfcn(i,j,n) = Ti
             ! If neither the adjustment nor the original category have ice in them... 
             else if (aicen(i,j,n) == 0.0_r8) then
                vicen(i,j,n) = 0.0_r8
                vsnon(i,j,n) = 0.0_r8
             endif
          enddo
       enddo
    enddo
 
    deallocate(aice, vice, vsno, aice_temp, vice_temp, vsno_temp)
    deallocate(hin_max, hcat_midpoint)
 
 write(*,*) 'finishing squeezing and reinitalizing associated variables...'
 
 end subroutine cice_rebalancing
 
 !------------------------------------------------------------------------

end module ice_postprocessing_mod