module user_initialization
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
!*  By Robert Hallberg, April 1994 - June 2002                         *
!*                                                                     *
!*    This subroutine initializes the fields for the simulations.      *
!*  The one argument passed to initialize, Time, is set to the         *
!*  current time of the simulation.  The fields which are initialized  *
!*  here are:                                                          *
!*    u - Zonal velocity in m s-1.                                     *
!*    v - Meridional velocity in m s-1.                                *
!*    h - Layer thickness in m.  (Must be positive.)                   *
!*    G%bathyT - Basin depth in m.  (Must be positive.)                *
!*    G%CoriolisBu - The Coriolis parameter, in s-1.                   *
!*    G%g_prime - The reduced gravity at each interface, in m s-2.     *
!*    G%Rlay - Layer potential density (coordinate variable) in kg m-3.*
!*  If ENABLE_THERMODYNAMICS is defined:                               *
!*    T - Temperature in C.                                            *
!*    S - Salinity in psu.                                             *
!*  If BULKMIXEDLAYER is defined:                                      *
!*    Rml - Mixed layer and buffer layer potential densities in        *
!*          units of kg m-3.                                           *
!*  If SPONGE is defined:                                              *
!*    A series of subroutine calls are made to set up the damping      *
!*    rates and reference profiles for all variables that are damped   *
!*    in the sponge.                                                   *
!*  Any user provided tracer code is also first linked through this    *
!*  subroutine.                                                        *
!*                                                                     *
!*    Forcing-related fields (taux, tauy, buoy, ustar, etc.) are set   *
!*  in MOM_surface_forcing.F90.                                        *
!*                                                                     *
!*    These variables are all set in the set of subroutines (in this   *
!*  file) USER_initialize_bottom_depth, USER_initialize_thickness,     *
!*  USER_initialize_velocity,  USER_initialize_temperature_salinity,   *
!*  USER_initialize_mixed_layer_density, USER_initialize_sponges,      *
!*  USER_set_coord, and USER_set_ref_profile.                          *
!*                                                                     *
!*    The names of these subroutines should be self-explanatory. They  *
!*  start with "USER_" to indicate that they will likely have to be    *
!*  modified for each simulation to set the initial conditions and     *
!*  boundary conditions.  Most of these take two arguments: an integer *
!*  argument specifying whether the fields are to be calculated        *
!*  internally or read from a NetCDF file; and a string giving the     *
!*  path to that file.  If the field is initialized internally, the    *
!*  path is ignored.                                                   *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v, tauy                                  *
!*    j    x ^ x ^ x   At >:  u, taux                                  *
!*    j    > o > o >   At o:  h, bathyT, buoy, tr, T, S, Rml, ustar    *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type, add_tracer_OBC_values
use MOM_variables, only : thermo_var_ptrs, ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_variables, only : OBC_FLATHER_E, OBC_FLATHER_W, OBC_FLATHER_N, OBC_FLATHER_S
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
implicit none ; private

#include <MOM_memory.h>

public USER_set_coord, USER_initialize_topography, USER_initialize_thickness
public USER_initialize_velocity, USER_init_temperature_salinity
public USER_init_mixed_layer_density, USER_initialize_sponges
public USER_set_Open_Bdry_Conds, USER_set_rotation
! Joe
public MOM_read_topodrag, read_topodrag, init_topodrag, topo_drag
public set_rotation_fath,USER_param
! endJoe

logical :: first_call = .true.

contains

subroutine USER_set_coord(Rlay, g_prime, G, param_file, eqn_of_state)
  real, dimension(:), intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(EOS_type),        pointer    :: eqn_of_state

  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_set_coord: " // &
   "Unmodified user routine called - you must edit the routine to use it")
  Rlay(:) = 0.0
  g_prime(:) = 0.0
  
  if (first_call) call write_user_log(param_file)

end subroutine USER_set_coord

subroutine USER_initialize_topography(D, G, param_file)
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_initialize_topography: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  D(:,:) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_topography

!---Joe------------------
subroutine USER_initialize_thickness(h, G, param_file,T)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: h
  type(ocean_grid_type),                  intent(in)  :: G
  type(param_file_type),                  intent(in)  :: param_file
  real, intent(in), dimension(NIMEM_,NJMEM_, NKMEM_)  :: T
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  This subroutine initializes the layer thicknesses.
  character(len=40)  :: mod = "USER_initialize_thickness" ! This subroutine's name.
  real :: e0(SZK_(G)+1)   ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  integer :: i, j, k, is, ie, js, je, nz
  real, save :: tdepth

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (G%max_depth<=0.) call MOM_error(FATAL,"USER_initialize_thickness: "// &
      "MAXIMUM_DEPTH has a non-sensical value! Was it set?")

  do k=1,nz
    e0(K) = -G%max_depth * real(k-1) / real(nz)
  enddo
  e0(nz) = G%tdepth !comment this line to use uniform layer depths (negative)
                    !tdepth is defined in /core/mOM_grid.F90
                    !It is TOP_LAYER_DEPTH found in MOM_input

  do j=js,je ; do i=is,ie
!    This sets the initial thickness (in m) of the layers.  The      !
!  thicknesses are set to insure that: 1.  each layer is at least an !
!  Angstrom thick, and 2.  the interfaces are where they should be   !
!  based on the resting depths and interface height perturbations,   !
!  as long at this doesn't interfere with 1.                         !
    eta1D(nz+1) = -1.0*G%bathyT(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + G%Angstrom_z)) then
        eta1D(K) = eta1D(K+1) + G%Angstrom_z
        h(i,j,k) = G%Angstrom_z
      else
        h(i,j,k) = eta1D(K) - eta1D(K+1)
      endif
    enddo
  enddo ; enddo

  if (first_call) call write_user_log(param_file)
end subroutine USER_initialize_thickness
! endJoe
! -----------------------------------------------------------------------------
! Joe---Get upper layer depth and save it (see beginning of this module)-------
subroutine USER_param(G,param_file,tdepth)
  type(ocean_grid_type),                  intent(in)  :: G
  type(param_file_type),                  intent(in)  :: param_file
  real :: tdepth
! Arguments:
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.i
! (out)     tdepth - depth of upper layer (2-layer). It is
!                    saved (see the beginning of this module)
!  This subroutine reads in the initial layer depth for 2-layer runs
!  From MOM_input.

!  character(len=40)  :: mod = "USER_param" ! This subroutine's name.
  
! call get_param(param_file, mod, "TOP_LAYER_DEPTH",tdepth, &
!                 "Depth of upper layer (2-layer)", default=0.0)

!  if (first_call) call write_user_log(param_file)

end subroutine USER_param
! endJoe------------------------------------------

subroutine USER_initialize_velocity(u, v, G, param_file)
  real, dimension(NIMEMB_, NJMEM_, NKMEM_), intent(out) :: u
  real, dimension(NIMEM_, NJMEMB_, NKMEM_), intent(out) :: v
  type(ocean_grid_type),                    intent(in)  :: G
  type(param_file_type),                    intent(in)  :: param_file
! Joe
  integer i,j,k
! endJoe
 
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_initialize_velocity: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  u(:,:,1) = 0.0
  v(:,:,1) = 0.0

!!    This initializes the velocities to 0--from Brian's C code 
!    for (k=0;k<=NZ-1;k++) {
!      for (j=Y0;j<=ny;j++) {
!        for (i=X0;i<=nx;i++) {
!          u[0][k][j][i] = 0.0; v[0][k][j][i] = 0.0;
!       

 if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_velocity

subroutine USER_init_temperature_salinity(T, S, G, param_file, eqn_of_state)
  real, dimension(NIMEM_, NJMEM_, NKMEM_), intent(out) :: T, S
  type(ocean_grid_type),                   intent(in)  :: G
  type(param_file_type),                   intent(in)  :: param_file
  type(EOS_type),                          pointer     :: eqn_of_state
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_init_temperature_salinity: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  T(:,:,1) = 0.0
  S(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_init_temperature_salinity

subroutine USER_init_mixed_layer_density(Rml, G, param_file, use_temperature, &
                                         eqn_of_state, T, S, P_Ref)
  real, dimension(NIMEM_, NJMEM_, NKMEM_),       intent(out) :: Rml
  type(ocean_grid_type),                         intent(in)  :: G
  type(param_file_type),                         intent(in)  :: param_file
  logical,                                       intent(in)  :: use_temperature
  type(EOS_type),                      optional, pointer     :: eqn_of_state
  real, dimension(NIMEM_, NJMEM_, NKMEM_), optional, intent(in) :: T, S
  real,                                optional, intent(in)  :: P_Ref
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_init_mixed_layer_density: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  Rml(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_init_mixed_layer_density

subroutine USER_initialize_sponges(G, use_temperature, tv, param_file, CSp, h)
  type(ocean_grid_type), intent(in) :: G
  logical,               intent(in) :: use_temperature
  type(thermo_var_ptrs), intent(in) :: tv
  type(param_file_type), intent(in) :: param_file
  type(sponge_CS),       pointer    :: CSp
  real, dimension(NIMEM_, NJMEM_, NKMEM_), intent(in) :: h
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_initialize_sponges: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_sponges

subroutine USER_set_Open_Bdry_Conds(OBC, tv, G, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC
  type(thermo_var_ptrs),      intent(in) :: tv
  type(ocean_grid_type),      intent(in) :: G
  type(param_file_type),      intent(in) :: param_file
  type(tracer_registry_type), pointer    :: tr_Reg
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_set_Open_Bdry_Conds: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_Open_Bdry_Conds

subroutine USER_set_rotation(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_set_rotation: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_rotation

subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "user_initialization" ! This module's name.

  call log_version(param_file, mod, version)
  first_call = .false.

end subroutine write_user_log
!-----------------------------------------------------------------------
! Joe
subroutine set_rotation_fath(fath, G, param_file)
  type(ocean_grid_type),                        intent(in)  :: G
  real, dimension(G%isd:G%ied,G%jsd:G%jed), intent(out) :: fath
  type(param_file_type),                        intent(in)  :: param_file
! Arguments: fath          - Coriolis parameter at h points in s^-1
!     (in)   G          - grid type
!     (in)   param_file - parameter file type

! This subroutine sets up the Coriolis parameter for a sphere
  character(len=30) :: mod = "set_rotation_fath" ! This subroutine's name.
  integer :: I, J
  real    :: PI, omega

!  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, "set_rotation_fath", "OMEGA", omega, &
                 "The rotation rate of earth at h points.", units="s-1", &
                 default=7.2921e-5)
  PI = 4.0*atan(1.0)

!write(6,*) G%isd, G%ied

  do I=G%isd,G%ied ; do J=G%jsd,G%jed
    fath(I,J) = ( 2.0 * omega ) * sin( ( PI * G%geoLatT(I,J) ) / 180.)
  enddo ; enddo

!  call callTree_leave(trim(mod)//'()')
end subroutine set_rotation_fath
! -----------------------------------------------------------------------------
! endJoe
! Joe
subroutine MOM_read_topodrag(t11, t12,t21,t22,hkmin,hkmax,ss,dragmask, &
                                       G, PF)
  real, dimension(NIMEM_,NJMEM_), intent(out) :: t11
  real, dimension(NIMEM_,NJMEM_), intent(out) :: t12
  real, dimension(NIMEM_,NJMEM_), intent(out) :: t21
  real, dimension(NIMEM_,NJMEM_), intent(out) :: t22
  real, dimension(NIMEM_,NJMEM_), intent(out) :: hkmin
  real, dimension(NIMEM_,NJMEM_), intent(out) :: hkmax
  real, dimension(NIMEM_,NJMEM_), intent(out) :: ss
  real, dimension(NIMEM_,NJMEM_), intent(out) :: dragmask
  type(ocean_grid_type),          intent(in)  :: G
  type(param_file_type),          intent(in)  :: PF
! Arguments: t11  - tensor component t11
!  (in)      G  - The ocean's grid structure.
!  (in)      PF - A structure indicating the open file to parse for
!                         model parameter values.

  character(len=40)  :: mod = "MOM_read_topodrag" ! This subroutine's name.
!  character(len=200) :: config

   call read_topodrag(t11,t12,t21,t22,hkmin,hkmax,ss,dragmask, G, PF)

end subroutine MOM_read_topodrag
! endJoe
!!----------------------------------------------------------------------------
!! Joe
! -----------------------------------------------------------------------------
subroutine read_topodrag(t11,t12,t21,t22,hkmin,hkmax,ss,dragmask, &
                                       G, param_file )
  real, dimension(NIMEM_,NJMEM_), intent(out) :: t11
  real, dimension(NIMEM_,NJMEM_), intent(out) :: t12
  real, dimension(NIMEM_,NJMEM_), intent(out) :: t21
  real, dimension(NIMEM_,NJMEM_), intent(out) :: t22
  real, dimension(NIMEM_,NJMEM_), intent(out) :: hkmin
  real, dimension(NIMEM_,NJMEM_), intent(out) :: hkmax
  real, dimension(NIMEM_,NJMEM_), intent(out) :: dragmask
  real, dimension(NIMEM_,NJMEM_), intent(out) :: ss
  type(ocean_grid_type),          intent(in)  :: G
  type(param_file_type),          intent(in)  :: param_file
! Arguments: t11          - tensor component t11
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine reads the tensor component t11 from a file and puts it into t11(:,:).
  character(len=200) :: filename, topo_file, inputdir ! Strings for file/path
  character(len=200) :: topo_varname                  ! Variable name in file
  character(len=40)  :: mod = "read_topodrag" ! This subroutine's name.

!  call callTree_enter(trim(mod)//"(), user_initialization.F90")

  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

!  filename = trim(inputdir)//trim(topo_file)
  filename = trim(inputdir)//trim("t11_86to82.nc")
!  call log_param(param_file, mod, "INPUTDIR/TOPO_FILE", filename)
!  call log_param(param_file, mod, "INPUTDIR", filename)
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topodrag: Unable to open "//trim(filename))
!  call read_data(filename,trim(topo_varname),D,domain=G%Domain%mpp_domain)
  call read_data(filename,trim('T11'),t11,domain=G%Domain%mpp_domain)
!  call callTree_leave(trim(mod)//'()')
!!
  filename = trim(inputdir)//trim("t12_86to82.nc")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topodrag: Unable to open "//trim(filename))
  call read_data(filename,trim('T12'),t12,domain=G%Domain%mpp_domain) 
!!
  filename = trim(inputdir)//trim("t21_86to82.nc")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topodrag: Unable to open "//trim(filename))
  call read_data(filename,trim('T21'),t21,domain=G%Domain%mpp_domain)
 !!
  filename = trim(inputdir)//trim("t22_86to82.nc")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topodrag: Unable to open "//trim(filename))
  call read_data(filename,trim('T22'),t22,domain=G%Domain%mpp_domain)
 !!
  filename = trim(inputdir)//trim("hkmin_86to82.nc")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topodrag: Unable to open "//trim(filename))
  call read_data(filename,trim('del2chimin'),hkmin,domain=G%Domain%mpp_domain)
 !!
  filename = trim(inputdir)//trim("hkmax_86to82.nc")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topodrag: Unable to open "//trim(filename))
  call read_data(filename,trim('del2chimax'),hkmax,domain=G%Domain%mpp_domain)

 !!
  filename = trim(inputdir)//trim("n_matrix_levitus98_86to82.nc")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topodrag: Unable to open "//trim(filename))
  call read_data(filename,trim('N'),ss,domain=G%Domain%mpp_domain)
 !!
  filename = trim(inputdir)//trim("dragmask100m_86to82.nc")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topodrag: Unable to open "//trim(filename))
  call read_data(filename,trim('mask'),dragmask,domain=G%Domain%mpp_domain)
end subroutine read_topodrag
!!--------------------------------------------------------------------------
!endJoe

!!--------------------------------------------------------------------   
subroutine init_topodrag(fath, t11, t12, t21, t22, ss, & 
                           hkmin, hkmax, dragmask, dragfac, & 
                           ssharmonic, G,param_file)
			   
   type(ocean_grid_type), intent(in)  :: G
   real, dimension(NIMEM_,NJMEM_) :: fath        !/* The Coriolis parameter in s-1
   real, dimension(NIMEM_,NJMEM_) :: dragfac     !/* topo drag factor */
   real, dimension(NIMEM_,NJMEM_) :: dragmask    !
   real, dimension(NIMEM_,NJMEM_) :: t11   !
   real, dimension(NIMEM_,NJMEM_) :: t12   !
   real, dimension(NIMEM_,NJMEM_) :: t21   !
   real, dimension(NIMEM_,NJMEM_) :: t22   !
   real, dimension(NIMEM_,NJMEM_) :: hkmin  !
   real, dimension(NIMEM_,NJMEM_) :: hkmax  !
   real, dimension(NIMEM_,NJMEM_) :: ss !
   real, dimension(NIMEM_,NJMEM_) :: ssharmonic !
 
   type(param_file_type), intent(in)  :: param_file

   integer i,j
   real frmin, frmax, frclp, taulin, taup, taun, gterm
   real dummy1, frcrit, anonlin, beta, gamma1, omega1
   real rnew

   real PI
  character(len=30) :: mod = "init_topodrag" ! This subroutine's name.
   PI = 4.0*atan(1.0)
!   
   frcrit=1.0
   anonlin=10.0
   rnew = 1.0/(3.8*3600.0*24.0) ! (3.8day)^(-1) decay rate as in Ansong et al 2015
   beta=1.0
   gamma1=0.35
!   omega1=1.405189*0.0001
   
   call get_param(param_file, mod, "TIDE_M2_FREQ", omega1, &
                   "Frequency of the M2 tidal constituent. \n"//&
                   "This is only used if TIDES and TIDE_M2 \n"//&
                   " are true.", units="s-1", default=1.405198e-4)
!!    call set_rotation_fath(fath, G, param_file)
  
   do j= G%jsd,G%jed    !1,1345
    do i= G%isd,G%ied   !1,2880

   !  calculate simple harmonic equivalent of static stability 

      if (omega1 > abs(fath(i,j))) then
        ssharmonic(i,j)=ss(i,j)* &
         ( (omega1*omega1-fath(i,j)*fath(i,j))**0.5)*(1.0/omega1)
      endif
      
      if (omega1 <= abs(fath(i,j))) then
          ssharmonic(i,j)=0.0
      endif

   enddo
  enddo 
!!------------------  
   do j=G%jsd,G%jed
    do i=G%isd,G%ied 

   !  calculate min/max Froude numbers

       if (omega1 > abs(fath(i,j))) then
      
        frmin = 0.5*ssharmonic(i,j)*hkmin(i,j)* &
                   omega1/(omega1*omega1-fath(i,j)*fath(i,j))
        frmax = 0.5*ssharmonic(i,j)*hkmax(i,j)* &
                   omega1/(omega1*omega1-fath(i,j)*fath(i,j))
         
  !      if (frmax < frmin+1.0e-4) then
  !        frmax = frmin+1.0e-4
  !      endif
        if (frmax >= frmin+1.0e-4) then
            frmax = frmax
        else
           frmax = frmin+1.0e-4
        endif
   
        if (frmin >= frcrit) then
          dummy1 = frmin
        else
          dummy1 = frcrit
        endif
        
        if (frmax <= dummy1) then
          frclp = frmax
        else
        frclp = dummy1
        endif
        
        !!get total drag in linear limit  
        
        taulin =(frmax**(gamma1+2.0)-frmin**(gamma1+2.0))/(gamma1+2.0)
        gterm  =(frmax**(gamma1-beta)-frclp**(gamma1-beta))/(gamma1-beta)
        
        !!get propagating and nonpropagating parts of total drag 
        
        taup =  ( frclp**(gamma1+2.0)-frmin**(gamma1+2.0) )/(gamma1+2.0)+gterm 
        taun = anonlin*(1.0/(beta+1.0))*((1.0/(gamma1+1.0))*( frmax**(gamma1+1.0)-frclp**(gamma1+1.0) )-gterm)
        
       ! dragfac(i,j) = dragmask(i,j)*ssharmonic(i,j)*(taup+taun)/taulin
       dragfac(i,j) = dragmask(i,j)*rnew
        
      endif
      
      if (omega1 <= abs(fath(i,j)) ) then
        dragfac(i,j) = 0.0 
      endif

   enddo
  enddo  
!-------------    
  end subroutine init_topodrag
!--------------------------------------------------------------------------
!! Joe
 subroutine topo_drag(u, v, h, dt, G, t11, t12, t21, t22, dragfac)

  type(ocean_grid_type), intent(in)  :: G
  real dt
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), target, intent(inout) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), target, intent(inout) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_), target,  intent(inout) :: h

  real, dimension(NIMEM_,NJMEM_), intent(in) :: dragfac     !/* topo drag factor */
  real, dimension(NIMEM_,NJMEM_), intent(in) :: t11   !
  real, dimension(NIMEM_,NJMEM_), intent(in) :: t12   !
  real, dimension(NIMEM_,NJMEM_), intent(in) :: t21   !
  real, dimension(NIMEM_,NJMEM_), intent(in) :: t22   !
     
  real, dimension(SZI_(G),SZJ_(G)) :: um 
  real, dimension(SZI_(G),SZJ_(G)) :: vm
  real, dimension(SZI_(G),SZJ_(G)) :: deltaum
  real, dimension(SZI_(G),SZJ_(G)) :: deltavm
     
  real alin, d11,d21, d12, d22, rnorm

  real dummy2, dummy3

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: Isp, Iep, Jsp, Jep
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  Isp = G%Isd; Iep =G%Ied; Jsp = G%Jsd; Jep = G%Jed

  !!alin=4.0
  alin = G%alinj !defined in /core/MOM_grid.F90,
                 !value of DRAG_STRENGTH in MOM_input

   do j=Jsp+1,Jep; do I=Isp+1,Iep    ! do j=js,je; do I=Isq,Ieq
  ! interpolate velocities to mass points: */  
    um(i,j)= 0.5*(u(i-1,j,nz)+u(i,j,nz))
    vm(i,j)= 0.5*(v(i,j-1,nz)+v(i,j,nz))
   enddo; enddo
   
   ! update velocity  */

   do j= G%jsd,G%jed !js,je
    do i= G%isd,G%ied  !is,ie
          
       if (h(i,j,nz)>= 1.0 ) then
         dummy2 = h(i,j,nz)
       else
         dummy2 = 1.0
       endif
       
!       dummy3 = alin*dragfac(i,j)*dt/dummy2
       dummy3 = -alin*dragfac(i,j)*dt !negative added for d11, etc 

       d11 = dummy3 !t11(i,j)*dummy3
       d21 = 0.0    !t21(i,j)*dummy3
       d12 = 0.0    !t12(i,j)*dummy3
       d22 = dummy3 !t22(i,j)*dummy3
       rnorm = (1.0-d11)*(1.0-d22)-d21*d12
       deltaum(i,j) = -um(i,j)+(um(i,j)*(1.0-d22)+vm(i,j)*d21)/rnorm
       deltavm(i,j) = -vm(i,j)+(vm(i,j)*(1.0-d11)+um(i,j)*d12)/rnorm
       
     enddo
   enddo
   
   ! interpolate delta velocities back to the appropriate grids */
!   write(6,*) SIZE(deltaum,1), SIZE(deltaum,2)
     
   do J= G%jsd,G%jed-1 !Jsp,Jep-1 
    do I=G%isd,G%ied-1 !Isp,Iep-1 
 
    u(i,j,nz) =  u(i,j,nz)+0.5*(deltaum(i+1,j)+deltaum(i,j))
    v(i,j,nz) =  v(i,j,nz)+0.5*(deltavm(i,j+1)+deltavm(i,j))

    enddo
   enddo 
  !-----------
  end subroutine topo_drag
!!endJoe
end module user_initialization
