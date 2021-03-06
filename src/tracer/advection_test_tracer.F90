module advection_test_tracer
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, 2002                                           *
!*                                                                     *
!*    This file contains an example of the code that is needed to set  *
!*  up and use a set (in this case eleven) of dynamically passive      *
!*  tracers.  These tracers dye the inflowing water or water initially *
!*  within a range of latitudes or water initially in a range of       *
!*  depths.                                                            *
!*                                                                     *
!*    A single subroutine is called from within each file to register  *
!*  each of the tracers for reinitialization and advection and to      *
!*  register the subroutine that initializes the tracers and set up    *
!*  their output and the subroutine that does any tracer physics or    *
!*  chemistry along with diapycnal mixing (included here because some  *
!*  tracers may float or swim vertically or dye diapycnal processes).  *
!*                                                                     *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, tr                                    *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_diag_to_Z, only : register_Z_tracer, diag_to_Z_CS
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data, slasher, vardesc
use MOM_restart, only : register_restart_field, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type, get_time
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_registry, only : add_tracer_diagnostics, add_tracer_OBC_values
use MOM_tracer_registry, only : tracer_vertdiff
use MOM_variables, only : surface, ocean_OBC_type

use coupler_util, only : set_coupler_values, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_advection_test_tracer, initialize_advection_test_tracer
public advection_test_tracer_surface_state, advection_test_tracer_end
public advection_test_tracer_column_physics

! ntr is the number of tracers in this module.
integer, parameter :: ntr = 11

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: advection_test_tracer_CS ; private
  logical :: coupled_tracers = .false.  ! These tracers are not offered to the
                                        ! coupler.
  character(len = 200) :: tracer_IC_file ! The full path to the IC file, or " "
                                   ! to initialize internally.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL()
  real, pointer :: tr(:,:,:,:) => NULL()   ! The array of tracers used in this
                                           ! subroutine, in g m-3?
  real, pointer :: tr_aux(:,:,:,:) => NULL() ! The masked tracer concentration
                                             ! for output, in g m-3.
  type(p3d), dimension(NTR) :: &
    tr_adx, &! Tracer zonal advective fluxes in g m-3 m3 s-1.
    tr_ady, &! Tracer meridional advective fluxes in g m-3 m3 s-1.
    tr_dfx, &! Tracer zonal diffusive fluxes in g m-3 m3 s-1.
    tr_dfy   ! Tracer meridional diffusive fluxes in g m-3 m3 s-1.
  real :: land_val(NTR) = -1.0 ! The value of tr used where land is masked out.
  logical :: mask_tracers  ! If true, tracers are masked out in massless layers.
  logical :: use_sponge

  real :: x_origin, x_width ! Parameters describing the test functions
  real :: y_origin, y_width ! Parameters describing the test functions

  integer, dimension(NTR) :: ind_tr ! Indices returned by aof_set_coupler_flux
             ! if it is used and the surface tracer concentrations are to be
             ! provided to the coupler.

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer, dimension(NTR) :: id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1
  integer, dimension(NTR) :: id_tr_dfx = -1, id_tr_dfy = -1

  type(vardesc) :: tr_desc(NTR)
end type advection_test_tracer_CS

contains

function register_advection_test_tracer(G, param_file, CS, diag, tr_Reg, &
                                      restart_CS)
  type(ocean_grid_type),   intent(in) :: G
  type(param_file_type),   intent(in) :: param_file
  type(advection_test_tracer_CS),    pointer    :: CS
  type(diag_ctrl), target, intent(in) :: diag
  type(tracer_registry_type),  pointer    :: tr_Reg
  type(MOM_restart_CS),   pointer    :: restart_CS
! This subroutine is used to register tracer fields and subroutines
! to be used with MOM.
! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  tr_Reg - A pointer that is set to point to the control structure
!                  for the tracer advection and diffusion module.
!  (in)      restart_CS - A pointer to the restart control structure.
  character(len=80)  :: name, longname
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "advection_test_tracer" ! This module's name.
  character(len=200) :: inputdir
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_advection_test_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_advection_test_tracer called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

  CS%diag => diag
  call get_param(param_file, mod, "ADVECTION_TEST_X_ORIGIN", CS%x_origin, &
        "The x-coorindate of the center of the test-functions.\n", default=0.)
  call get_param(param_file, mod, "ADVECTION_TEST_Y_ORIGIN", CS%y_origin, &
        "The y-coorindate of the center of the test-functions.\n", default=0.)
  call get_param(param_file, mod, "ADVECTION_TEST_X_WIDTH", CS%x_width, &
        "The x-width of the test-functions.\n", default=0.)
  call get_param(param_file, mod, "ADVECTION_TEST_Y_WIDTH", CS%y_width, &
        "The y-width of the test-functions.\n", default=0.)
  call get_param(param_file, mod, "ADVECTION_TEST_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial \n"//&
                 "conditions for the tracers, or blank to initialize \n"//&
                 "them internally.", default=" ")

  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    CS%tracer_IC_file = trim(slasher(inputdir))//trim(CS%tracer_IC_file)
    call log_param(param_file, mod, "INPUTDIR/ADVECTION_TEST_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mod, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified from MOM_initialization.F90.", default=.false.)

  call get_param(param_file, mod, "MASK_TRACERS_IN_MASSLESS_LAYERS", CS%mask_tracers, &
                 "If true, tracers will be masked out in massless layers. \n", &
                 default=.false.)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR)) ; CS%tr(:,:,:,:) = 0.0
  if (CS%mask_tracers) then
    allocate(CS%tr_aux(isd:ied,jsd:jed,nz,NTR)) ; CS%tr_aux(:,:,:,:) = 0.0
  endif

  do m=1,NTR
    CS%tr_desc(m) = vardesc("tr","Tracer",'h','L','s',"kg kg-1")
    if (m < 10) then ; write(name,'("tr",I1.1)') m
    else ; write(name,'("tr",I2.2)') m ; endif
    write(longname,'("Concentration of Tracer ",I2.2)') m
    CS%tr_desc(m)%name = name
    CS%tr_desc(m)%longname = longname
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for the restart file.
    call register_restart_field(tr_ptr, CS%tr_desc(m), .true., restart_CS)
    ! Register the tracer for horizontal advection & diffusion.
    call register_tracer(tr_ptr, CS%tr_desc(m)%name, param_file, tr_Reg)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(CS%tr_desc(m)%name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_advection_test_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  register_advection_test_tracer = .true.
end function register_advection_test_tracer

subroutine initialize_advection_test_tracer(restart, day, G, h, OBC, CS, sponge_CSp, &
                                  diag_to_Z_CSp)
  logical,                            intent(in) :: restart
  type(time_type), target,            intent(in) :: day
  type(ocean_grid_type),              intent(in) :: G
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h
  type(ocean_OBC_type),               pointer    :: OBC
  type(advection_test_tracer_CS),               pointer    :: CS
  type(sponge_CS),                    pointer    :: sponge_CSp
  type(diag_to_Z_CS),                 pointer    :: diag_to_Z_CSp
!   This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.

! Arguments: restart - .true. if the fields have already been read from
!                     a restart file.
!  (in)      day - Time of the start of the run.
!  (in)      G - The ocean's grid structure.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in/out)  CS - The control structure returned by a previous call to
!                 register_advection_test_tracer.
!  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
!                         they are in use.  Otherwise this may be unassociated.
!  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
!                            in depth space.
  real, allocatable :: temp(:,:,:)
  real, pointer, dimension(:,:,:) :: &
    OBC_tr1_u => NULL(), & ! These arrays should be allocated and set to
    OBC_tr1_v => NULL()    ! specify the values of tracer 1 that should come
                           ! in through u- and v- points through the open
                           ! boundary conditions, in the same units as tr.
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for tracer fluxes, usually
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  real :: PI     ! 3.1415926... calculated as 4*atan(1)
  real :: tr_y   ! Initial zonally uniform tracer concentrations.
  real :: dist2  ! The distance squared from a line, in m2.
  real :: h_neglect         ! A thickness that is so small it is usually lost
                            ! in roundoff and can be neglected, in m.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB
  real :: tmpx, tmpy, locx, locy

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  h_neglect = G%H_subroundoff

  if (.not.restart) then
    do m=1,NTR
      do k=1,nz ; do j=js,je ; do i=is,ie
        CS%tr(i,j,k,m) = 0.0
      enddo ; enddo ; enddo
      k=1 ! Square wave
      do j=js,je ; do i=is,ie
        if (abs(G%geoLonT(i,j)-CS%x_origin)<0.5*CS%x_width .and. &
            abs(G%geoLatT(i,j)-CS%y_origin)<0.5*CS%y_width) CS%tr(i,j,k,m) = 1.0
      enddo ; enddo
      k=2 ! Triangle wave
      do j=js,je ; do i=is,ie
        locx=abs(G%geoLonT(i,j)-CS%x_origin)/CS%x_width
        locy=abs(G%geoLatT(i,j)-CS%y_origin)/CS%y_width
        CS%tr(i,j,k,m) = max(0.0, 1.0-locx)*max(0.0, 1.0-locy)
      enddo ; enddo
      k=3 ! Cosine bell
      do j=js,je ; do i=is,ie
        locx=min(1.0, abs(G%geoLonT(i,j)-CS%x_origin)/CS%x_width)*(acos(0.0)*2.)
        locy=min(1.0, abs(G%geoLatT(i,j)-CS%y_origin)/CS%y_width)*(acos(0.0)*2.)
        CS%tr(i,j,k,m) = (1.0+cos(locx))*(1.0+cos(locy))*0.25
      enddo ; enddo
      k=4 ! Cylinder
      do j=js,je ; do i=is,ie
        locx=abs(G%geoLonT(i,j)-CS%x_origin)/CS%x_width
        locy=abs(G%geoLatT(i,j)-CS%y_origin)/CS%y_width
        if (locx**2+locy**2<=1.0) CS%tr(i,j,k,m) = 1.0
      enddo ; enddo
      k=5 ! Cut cylinder
      do j=js,je ; do i=is,ie
        locx=(G%geoLonT(i,j)-CS%x_origin)/CS%x_width
        locy=(G%geoLatT(i,j)-CS%y_origin)/CS%y_width
        if (locx**2+locy**2<=1.0) CS%tr(i,j,k,m) = 1.0
        if (locx>0.0.and.abs(locy)<0.2) CS%tr(i,j,k,m) = 0.0
      enddo ; enddo
    enddo
  endif ! restart

  ! This needs to be changed if the units of tracer are changed above.
  if (G%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
  else ; flux_units = "kg s-1" ; endif

  do m=1,NTR
    ! Register the tracer for the restart file.
    name = CS%tr_desc(m)%name ; longname = CS%tr_desc(m)%longname
    units = CS%tr_desc(m)%units
    CS%id_tracer(m) = register_diag_field("ocean_model", trim(name), CS%diag%axesTL, &
        day, trim(longname) , trim(units))
    CS%id_tr_adx(m) = register_diag_field("ocean_model", trim(name)//"_adx", &
        CS%diag%axesCuL, day, trim(longname)//" advective zonal flux" , &
        trim(flux_units))
    CS%id_tr_ady(m) = register_diag_field("ocean_model", trim(name)//"_ady", &
        CS%diag%axesCvL, day, trim(longname)//" advective meridional flux" , &
        trim(flux_units))
    CS%id_tr_dfx(m) = register_diag_field("ocean_model", trim(name)//"_dfx", &
        CS%diag%axesCuL, day, trim(longname)//" diffusive zonal flux" , &
        trim(flux_units))
    CS%id_tr_dfy(m) = register_diag_field("ocean_model", trim(name)//"_dfy", &
        CS%diag%axesCvL, day, trim(longname)//" diffusive zonal flux" , &
        trim(flux_units))
    if (CS%id_tr_adx(m) > 0) call safe_alloc_ptr(CS%tr_adx(m)%p,IsdB,IedB,jsd,jed,nz)
    if (CS%id_tr_ady(m) > 0) call safe_alloc_ptr(CS%tr_ady(m)%p,isd,ied,JsdB,JedB,nz)
    if (CS%id_tr_dfx(m) > 0) call safe_alloc_ptr(CS%tr_dfx(m)%p,IsdB,IedB,jsd,jed,nz)
    if (CS%id_tr_dfy(m) > 0) call safe_alloc_ptr(CS%tr_dfy(m)%p,isd,ied,JsdB,JedB,nz)

!    Register the tracer for horizontal advection & diffusion.
    if ((CS%id_tr_adx(m) > 0) .or. (CS%id_tr_ady(m) > 0) .or. &
        (CS%id_tr_dfx(m) > 0) .or. (CS%id_tr_dfy(m) > 0)) &
      call add_tracer_diagnostics(name, CS%tr_Reg, CS%tr_adx(m)%p, &
                                  CS%tr_ady(m)%p,CS%tr_dfx(m)%p,CS%tr_dfy(m)%p)

    call register_Z_tracer(CS%tr(:,:,:,m), trim(name), longname, units, &
                           day, G, diag_to_Z_CSp)
  enddo

end subroutine initialize_advection_test_tracer


subroutine advection_test_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, CS)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                         intent(in) :: fluxes
  real,                                  intent(in) :: dt
  type(ocean_grid_type),                 intent(in) :: G
  type(advection_test_tracer_CS),        pointer    :: CS
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      h_new -  Layer thickness after entrainment, in m or kg m-2.
!  (in)      ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m or kg m-2.
!  (in)      eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m or kg m-2.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_advection_test_tracer.
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

  real :: b1(SZI_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))  ! tridiagonal solver.
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) return

! do m=1,NTR
!   call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G)
! enddo

  if (CS%mask_tracers) then
    do m = 1,NTR ; if (CS%id_tracer(m) > 0) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        if (h_new(i,j,k) < 1.1*G%Angstrom) then
          CS%tr_aux(i,j,k,m) = CS%land_val(m)
        else
          CS%tr_aux(i,j,k,m) = CS%tr(i,j,k,m)
        endif
      enddo ; enddo ; enddo
    endif ; enddo
  endif

  do m=1,NTR
    if (CS%mask_tracers) then
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr_aux(:,:,:,m),CS%diag)
    else
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr(:,:,:,m),CS%diag)
    endif
    if (CS%id_tr_adx(m)>0) &
      call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(m)>0) &
      call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfx(m)>0) &
      call post_data(CS%id_tr_dfx(m),CS%tr_dfx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfy(m)>0) &
      call post_data(CS%id_tr_dfy(m),CS%tr_dfy(m)%p(:,:,:),CS%diag)
  enddo
end subroutine advection_test_tracer_column_physics

subroutine advection_test_tracer_surface_state(state, h, G, CS)
  type(surface),                         intent(inout) :: state
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h
  type(ocean_grid_type),                 intent(in)    :: G
  type(advection_test_tracer_CS),        pointer       :: CS
!   This particular tracer package does not report anything back to the coupler.
! The code that is here is just a rough guide for packages that would.
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_advection_test_tracer.
  integer :: i, j, m, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  
  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,ntr
      !   This call loads the surface vlues into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_values(CS%tr(:,:,1,1), state%tr_fields, CS%ind_tr(m), &
                              ind_csurf, is, ie, js, je)
    enddo
  endif
end subroutine advection_test_tracer_surface_state

subroutine advection_test_tracer_end(CS)
  type(advection_test_tracer_CS), pointer :: CS
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    if (associated(CS%tr_aux)) deallocate(CS%tr_aux)
    do m=1,NTR
      if (associated(CS%tr_adx(m)%p)) deallocate(CS%tr_adx(m)%p)
      if (associated(CS%tr_ady(m)%p)) deallocate(CS%tr_ady(m)%p)
      if (associated(CS%tr_dfx(m)%p)) deallocate(CS%tr_dfx(m)%p)
      if (associated(CS%tr_dfy(m)%p)) deallocate(CS%tr_dfy(m)%p)
    enddo

    deallocate(CS)
  endif
end subroutine advection_test_tracer_end

end module advection_test_tracer
