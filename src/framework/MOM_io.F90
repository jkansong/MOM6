module MOM_io

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
!*  By R. Hallberg July 1999 - June 2002                               *
!*                                                                     *
!*   This file contains a number of subroutines that manipulate        *
!*  NetCDF files and handle input and output of fields.  These         *
!*  subroutines, along with their purpose, are:                        *
!*                                                                     *
!*   create_file: create a new file and set up structures that are     *
!*       needed for subsequent output and write out the coordinates.   *
!*   reopen_file: reopen an existing file for writing and set up       *
!*       structures that are needed for subsequent output.             *
!*   open_input_file: open the indicated file for reading only.        *
!*   close_file: close an open file.                                   *
!*   synch_file: flush the buffers, completing all pending output.     *
!*                                                                     *
!*   write_field: write a field to an open file.                       *
!*   write_time: write a value of the time axis to an open file.       *
!*   read_field: read a field from an open file.                       *
!*   read_time: read a time from an open file.                         *
!*                                                                     *
!*   name_output_file: provide a name for an output file based on a    *
!*       name root and the time of the output.                         *
!*   find_input_file: find a file that has been previously written by  *
!*       MOM and named by name_output_file and open it for reading.    *
!*                                                                     *
!*   handle_error: write an error code and quit.                       *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**


use MOM_error_handler, only : MOM_error, NOTE, FATAL, WARNING
use MOM_domains, only : MOM_domain_type
use MOM_file_parser, only : log_version, param_file_type
use MOM_string_functions, only : lowercase
use MOM_grid, only : ocean_grid_type

use fms_mod, only : write_version_number, open_namelist_file, check_nml_error
use fms_io_mod, only : file_exist, field_size, read_data
use fms_io_mod, only : field_exists => field_exist, io_infra_end=>fms_io_exit
use mpp_domains_mod, only : domain1d, mpp_get_domain_components
use mpp_domains_mod, only : CENTER, CORNER, NORTH_FACE=>NORTH, EAST_FACE=>EAST
use mpp_io_mod, only : open_file => mpp_open, close_file => mpp_close
use mpp_io_mod, only : mpp_write_meta, write_field => mpp_write, mpp_get_info
use mpp_io_mod, only : mpp_get_atts, mpp_get_axes, mpp_get_axis_data, axistype
use mpp_io_mod, only : mpp_get_fields, fieldtype, axistype, flush_file => mpp_flush
use mpp_io_mod, only : APPEND_FILE=>MPP_APPEND, ASCII_FILE=>MPP_ASCII
use mpp_io_mod, only : MULTIPLE=>MPP_MULTI, NETCDF_FILE=>MPP_NETCDF
use mpp_io_mod, only : OVERWRITE_FILE=>MPP_OVERWR, READONLY_FILE=>MPP_RDONLY
use mpp_io_mod, only : SINGLE_FILE=>MPP_SINGLE, WRITEONLY_FILE=>MPP_WRONLY
use mpp_io_mod, only : MPP_APPEND, MPP_MULTI, MPP_OVERWR, MPP_NETCDF, MPP_RDONLY
use mpp_io_mod, only : get_file_info=>mpp_get_info, get_file_atts=>mpp_get_atts
use mpp_io_mod, only : get_file_fields=>mpp_get_fields, get_file_times=>mpp_get_times
use mpp_io_mod, only : read_field=>mpp_read, io_infra_init=>mpp_io_init

use netcdf

implicit none ; private

public :: close_file, create_file, field_exists, field_size, fieldtype
public :: file_exists, flush_file, get_file_info, get_file_atts, get_file_fields
public :: get_file_times, open_file, read_axis_data, read_data, read_field
public :: num_timelevels, MOM_read_data
public :: reopen_file, slasher, write_field, write_version_number, MOM_io_init
public :: open_namelist_file, check_nml_error, io_infra_init, io_infra_end
public :: APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
public :: READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE

type, public :: vardesc
  character(len=64)  :: name     ! The variable name in a NetCDF file.
  character(len=240) :: longname ! The long name of the variable.
  character(len=8)   :: hor_grid ! The hor. grid:  u, v, h, q, Cu, Cv, T, Bu, or 1.
  character(len=8)   :: z_grid   ! The vert. grid:  L, i, or 1.
  character(len=8)   :: t_grid   ! The time description: s, p, or 1.
  character(len=48)  :: units    ! The dimensions of the variable.
end type vardesc

interface file_exists
  module procedure file_exist
  module procedure MOM_file_exists
end interface

interface MOM_read_data
  module procedure MOM_read_data_3d
  module procedure MOM_read_data_2d
  module procedure MOM_read_data_1d
end interface


contains

subroutine create_file(unit, filename, vars, novars, G, fields, threading, &
                       timeunit)
  integer,               intent(out)   :: unit
  character(len=*),      intent(in)    :: filename
  type(vardesc),         intent(in)    :: vars(:)
  integer,               intent(in)    :: novars
  type(ocean_grid_type), intent(in)    :: G
  type(fieldtype),       intent(inout) :: fields(:)
  integer, optional,     intent(in)    :: threading
  real, optional,        intent(in)    :: timeunit
!   create_file creates a new NetCDF file.  It also sets up
! structures that describe this file and the variables that will
! later be written to this file.
! Arguments: unit - The unit id of an open file or -1 on a
!                   nonwriting PE with single file output.
!  (in)      filename - The full path to the file to create.
!  (in)      vars - An array of structures describing the fields
!                   that will be written to this file.
!  (in)      novars - The number of fields that will be written to
!                     this file.
!  (in)      G - The ocean's grid structure.
!  (out)     fields - An array of fieldtypes for each variable.
!  (in,opt)  threading - SINGLE_FILE or MULTIPLE, optional.
!  (in,opt)  timeunit - The length, in seconds, of the units for time. The
!                       default value is 86400.0, for 1 day.

  logical :: use_lath, use_lonh, use_latq, use_lonq, use_time
  logical :: use_layer, use_int, use_periodic
  type(axistype) :: axis_lath, axis_latq, axis_lonh, axis_lonq
  type(axistype) :: axis_layer, axis_int, axis_time, axis_periodic
  type(axistype) :: axes(4)
  type(domain1d) :: x_domain, y_domain
  integer :: numaxes, pack, thread, k, nz
  integer :: var_periods, num_periods=0
  real :: layer_val(G%ks:G%ke), interface_val(G%ks:G%ke+1)
  real, dimension(:), allocatable :: period_val
  character(len=40) :: time_units
  character(len=8) :: t_grid, t_grid_read

  nz = G%ke

  use_lath = .false. ; use_lonh = .false.
  use_latq = .false. ; use_lonq = .false.
  use_time = .false. ; use_periodic = .false.
  use_layer = .false. ; use_int = .false.

  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  if ((thread == SINGLE_FILE) .or. .not.G%Domain%use_io_layout) then
    call open_file(unit, filename, MPP_OVERWR, MPP_NETCDF, threading=thread)
  else
    call open_file(unit, filename, MPP_OVERWR, MPP_NETCDF, domain=G%Domain%mpp_domain)
  endif

  interface_val(1:nz+1) = G%GV%sInterface(1:nz+1)
  layer_val(1:nz) = G%GV%sLayer(1:nz)

  call mpp_get_domain_components(G%Domain%mpp_domain, x_domain, y_domain)

! Define the coordinates.
  do k=1,novars
    select case (vars(k)%hor_grid)
      case ('h') ; use_lath = .true. ; use_lonh = .true.
      case ('q') ; use_latq = .true. ; use_lonq = .true.
      case ('u') ; use_lath = .true. ; use_lonq = .true.
      case ('v') ; use_latq = .true. ; use_lonh = .true.
      case ('T')  ; use_lath = .true. ; use_lonh = .true.
      case ('Bu') ; use_latq = .true. ; use_lonq = .true.
      case ('Cu') ; use_lath = .true. ; use_lonq = .true.
      case ('Cv') ; use_latq = .true. ; use_lonh = .true.
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized hor_grid "//trim(vars(k)%hor_grid))
    end select
    select case (vars(k)%z_grid)
      case ('L') ; use_layer = .true.
      case ('i') ; use_int = .true.
      case ('1') ! Do nothing.
      case default
        call MOM_error(FATAL, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized z_grid "//trim(vars(k)%z_grid))
    end select
    t_grid = adjustl(vars(k)%t_grid)
    select case (t_grid(1:1))
      case ('s', 'a', 'm') ; use_time = .true.
      case ('p') ; use_periodic = .true.
        if (len_trim(t_grid(2:8)) <= 0) call MOM_error(FATAL, &
          "MOM_io create_file: No periodic axis length was specified in "//&
          trim(vars(k)%t_grid) // " in the periodic axes of variable "//&
          trim(vars(k)%name)//" in file "//trim(filename))
        var_periods = -9999999
        t_grid_read = adjustl(t_grid(2:8))
        read(t_grid_read,*) var_periods
        if (var_periods == -9999999) call MOM_error(FATAL, &
          "MOM_io create_file: Failed to read the number of periods from "//&
          trim(vars(k)%t_grid) // " in the periodic axes of variable "//&
          trim(vars(k)%name)//" in file "//trim(filename))
        if (var_periods < 1) call MOM_error(FATAL, "MOM_io create_file: "//&
           "variable "//trim(vars(k)%name)//" in file "//trim(filename)//&
           " uses a periodic time axis, and must have a positive "//&
           "value for the number of periods in "//vars(k)%t_grid )
        if ((num_periods > 0) .and. (var_periods /= num_periods)) &
          call MOM_error(FATAL, "MOM_io create_file: "//&
            "Only one value of the number of periods can be used in the "//&
            "create_file call for file "//trim(filename)//".  The second is "//&
            "variable "//trim(vars(k)%name)//" with t_grid "//vars(k)%t_grid )
  
        num_periods = var_periods
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized t_grid "//trim(vars(k)%t_grid))
    end select
  end do
! Specify all optional arguments to mpp_write_meta: name, units, longname, cartesian, calendar, sense, domain, data, min)
! Otherwise if optional arguments are added to mpp_write_meta the compiler may (and in case of GNU is) get confused abd crash
  if (use_lath) &
    call mpp_write_meta(unit, axis_lath, name="lath", units=G%y_axis_units, longname="Latitude", &
                   cartesian='Y', domain = y_domain, data=G%gridLatT(G%jsg:G%jeg))

  if (use_lonh) &
    call mpp_write_meta(unit, axis_lonh, name="lonh", units=G%x_axis_units, longname="Longitude", &
                   cartesian='X', domain = x_domain, data=G%gridLonT(G%isg:G%ieg))

  if (use_latq) &
    call mpp_write_meta(unit, axis_latq, name="latq", units=G%y_axis_units, longname="Latitude", &
                   cartesian='Y', domain = y_domain, data=G%gridLatB(G%JsgB:G%JegB))

  if (use_lonq) &
    call mpp_write_meta(unit, axis_lonq, name="lonq", units=G%x_axis_units, longname="Longitude", &
                   cartesian='X', domain = x_domain, data=G%gridLonB(G%IsgB:G%IegB))

  if (use_layer) &
    call mpp_write_meta(unit, axis_layer, name="Layer", units=trim(G%GV%zAxisUnits), &
          longname="Layer "//trim(G%GV%zAxisLongName), cartesian='Z', sense=1, data=layer_val)

  if (use_int) &
    call mpp_write_meta(unit, axis_int, name="Interface", units=trim(G%GV%zAxisUnits), &
          longname="Interface "//trim(G%GV%zAxisLongName), cartesian='Z', sense=1, data=interface_val)

  if (use_time) then ; if (present(timeunit)) then
    ! Set appropriate units, depending on the value.
    if (timeunit < 0.0) then
      time_units = "days" ! The default value.
    else if ((timeunit >= 0.99) .and. (timeunit < 1.01)) then
      time_units = "seconds"
    else if ((timeunit >= 3599.0) .and. (timeunit < 3601.0)) then
      time_units = "hours"
    else if ((timeunit >= 86399.0) .and. (timeunit < 86401.0)) then
      time_units = "days"
    else if ((timeunit >= 3.0e7) .and. (timeunit < 3.2e7)) then
      time_units = "years"
    else
      write(time_units,'(es8.2," s")') timeunit
    endif

    call mpp_write_meta(unit, axis_time, name="Time", units=time_units, longname="Time", cartesian='T')
  else
    call mpp_write_meta(unit, axis_time, name="Time", units="days", longname="Time",cartesian= 'T')
  endif ; endif

  if (use_periodic) then
    if (num_periods <= 1) call MOM_error(FATAL, "MOM_io create_file: "//&
      "num_periods for file "//trim(filename)//" must be at least 1.")
    ! Define a periodic axis with unit labels.
    allocate(period_val(num_periods))
    do k=1,num_periods ; period_val(k) = real(k) ; enddo
    call mpp_write_meta(unit, axis_periodic, name="Period", units="nondimensional", &
          longname="Periods for cyclical varaiables", cartesian= 't', data=period_val)
    deallocate(period_val)
  endif

  do k=1,novars
    numaxes = 0
    select case (vars(k)%hor_grid)
      case ('h') ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath
      case ('q') ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq
      case ('u') ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath
      case ('v') ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq
      case ('T')  ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath
      case ('Bu') ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq
      case ('Cu') ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath
      case ('Cv') ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized hor_grid "//trim(vars(k)%hor_grid))
    end select
    select case (vars(k)%z_grid)
      case ('L') ; numaxes = numaxes+1 ; axes(numaxes) = axis_layer
      case ('i') ; numaxes = numaxes+1 ; axes(numaxes) = axis_int
      case ('1') ! Do nothing.
      case default
        call MOM_error(FATAL, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized z_grid "//trim(vars(k)%z_grid))
    end select
    t_grid = adjustl(vars(k)%t_grid)
    select case (t_grid(1:1))
      case ('s', 'a', 'm') ; numaxes = numaxes+1 ; axes(numaxes) = axis_time
      case ('p') ; numaxes = numaxes+1 ; axes(numaxes) = axis_periodic
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized t_grid "//trim(vars(k)%t_grid))
    end select
    pack = 1

    call mpp_write_meta(unit, fields(k), axes(1:numaxes), vars(k)%name, vars(k)%units, &
           vars(k)%longname, pack = pack)
  enddo

  if (use_lath) call write_field(unit, axis_lath)
  if (use_latq) call write_field(unit, axis_latq)
  if (use_lonh) call write_field(unit, axis_lonh)
  if (use_lonq) call write_field(unit, axis_lonq)
  if (use_layer) call write_field(unit, axis_layer)
  if (use_int) call write_field(unit, axis_int)
  if (use_periodic) call write_field(unit, axis_periodic)

end subroutine create_file

subroutine reopen_file(unit, filename, vars, novars, G, fields, threading, timeunit)
  integer,               intent(out)   :: unit
  character(len=*),      intent(in)    :: filename
  type(vardesc),         intent(in)    :: vars(:)
  integer,               intent(in)    :: novars
  type(ocean_grid_type), intent(in)    :: G
  type(fieldtype),       intent(inout) :: fields(:)
  integer, optional,     intent(in)    :: threading
  real, optional,        intent(in)    :: timeunit
!   reopen_file opens an existing NetCDF file for output.  If it
! does not find the file, a new file is created.  It also sets up
! structures that describe this file and the variables that will
! later be written to this file.
! Arguments: unit - The unit id of an open file or -1 on a
!                   nonwriting PE with single file output.
!  (in)      filename - The full path to the file to create.
!  (in)      vars - An array of structures describing the fields
!                   that will be written to this file.
!  (in)      novars - The number of fields that will be written to
!                     this file.
!  (in)      G - The ocean's grid structure.
!  (out)     fields - An array of fieldtypes for each variable.
!  (in)      threading - SINGLE_FILE or MULTIPLE, optional.
!  (in,opt)  timeunit - The length, in seconds, of the units for time. The
!                       default value is 86400.0, for 1 day.
  character(len=200) :: check_name, mesg
  integer :: length, ndim, nvar, natt, ntime, thread
  logical :: exists

  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  check_name = filename
  length = len(trim(check_name))
  if (check_name(length-2:length) /= ".nc") check_name = trim(check_name)//".nc"
  if (thread /= SINGLE_FILE) check_name = trim(check_name)//".0000"

  inquire(file=check_name,EXIST=exists)

  if (.not.exists) then
    call create_file(unit, filename, vars, novars, G, fields, threading, timeunit)
  else
    if ((thread == SINGLE_FILE) .or. .not.G%Domain%use_io_layout) then
      call open_file(unit, filename, MPP_APPEND, MPP_NETCDF, threading=thread)
    else
      call open_file(unit, filename, MPP_APPEND, MPP_NETCDF, domain=G%Domain%mpp_domain)
    endif
    if (unit < 0) return

    call mpp_get_info(unit, ndim, nvar, natt, ntime)

    if (nvar /= novars) then
      write (mesg,*) "Reopening file ",trim(filename)," with ",novars,&
                     " variables instead of ",nvar,"."
      call MOM_error(FATAL,"MOM_io: "//mesg)
    endif

    call mpp_get_fields(unit,fields(1:nvar))

    ! Check the field names...
!    do i=1,nvar
!      call mpp_get_field_atts(fields(i),name)
!      !if (trim(name) /= trim(vars%name) then
!      !write (mesg,'("Reopening file ",a," variable ",a," is called ",a,".")',&
!      !    filename,vars%name,name);
!      !call MOM_error(NOTE,"MOM_io: "//mesg)
!    enddo
  endif

end subroutine reopen_file

subroutine read_axis_data(filename, axis_name, var)
  character(len=*),   intent(in)  :: filename, axis_name
  real, dimension(:), intent(out) :: var

  integer :: i,len,unit, ndim, nvar, natt, ntime
  logical :: axis_found
  type(axistype), allocatable :: axes(:)
  type(axistype) :: time_axis
  character(len=32) :: name, units

  call open_file(unit, trim(filename), action=MPP_RDONLY, form=MPP_NETCDF, &
                 threading=MPP_MULTI, fileset=SINGLE_FILE)
!Find the number of variables (nvar) in this file
  call mpp_get_info(unit, ndim, nvar, natt, ntime)
! -------------------------------------------------------------------
! Allocate space for the number of axes in the data file.
! -------------------------------------------------------------------
  allocate(axes(ndim))
  call mpp_get_axes(unit, axes, time_axis)

  axis_found = .false.
  do i = 1, ndim
    call mpp_get_atts(axes(i), name=name,len=len,units=units)
    if (name == axis_name) then
      axis_found = .true.
      call mpp_get_axis_data(axes(i),var)
      exit
    endif
  enddo

  if (.not.axis_found) call MOM_error(FATAL, "MOM_io read_axis_data: "//&
    "Unable to find axis "//trim(axis_name)//" in file "//trim(filename))

  deallocate(axes)

end subroutine read_axis_data

function num_timelevels(filename, varname, min_dims) result(n_time)
  character(len=*),  intent(in) :: filename, varname
  integer, optional, intent(in) :: min_dims
  integer :: n_time

!  This subroutine determines how many time levels a variable has.

! Arguments: filename - The name of the file to read.
!  (in)      varname - The name of the variable whose number of time levels
!                      are to be returned.
!  (in, opt) min_dims - The minimum number of dimensions a variable must have
!                       if it has a time dimension.  If the variable has 1 less
!                       dimension than this, 0 is returned.
!  (result)  n_time - The number of time levels varname has in file filename.

  logical :: found
  character(len=200) :: msg
  character(len=nf90_max_name) :: name
  integer :: ncid, nvars, status, varid, ndims, n
  integer, allocatable :: varids(:)
  integer, dimension(nf90_max_var_dims) :: dimids

  n_time = -1
  found = .false.

  status = NF90_OPEN(filename, NF90_NOWRITE, ncid)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " Difficulties opening "//trim(filename)//" - "//&
        trim(NF90_STRERROR(status)))
    return
  endif

  status = NF90_INQUIRE(ncid, nVariables=nvars)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " Difficulties getting the number of variables in file "//&
        trim(filename)//" - "//trim(NF90_STRERROR(status)))
    return
  endif

  if (nvars < 1) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " There appear not to be any variables in "//trim(filename))
    return
  endif


  allocate(varids(nvars))

  status = nf90_inq_varids(ncid, nvars, varids)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " Difficulties getting the variable IDs in file "//&
        trim(filename)//" - "//trim(NF90_STRERROR(status)))
    deallocate(varids) ; return
  endif
  
  do n = 1,nvars
    status = nf90_inquire_variable(ncid, varids(n), name=name)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING,"num_timelevels: "//&
          " Difficulties getting a variable name in file "//&
          trim(filename)//" - "//trim(NF90_STRERROR(status)))
    endif

    if (trim(lowercase(name)) == trim(lowercase(varname))) then
      if (found) then
        call MOM_error(WARNING,"num_timelevels: "//&
          " Two variables match the case-insensitive name "//trim(varname)//&
          " in file "//trim(filename)//" - "//trim(NF90_STRERROR(status)))
      else
        varid = varids(n) ; found = .true.
      endif
    endif
  enddo    

  deallocate(varids)

  if (.not.found) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " variable "//trim(varname)//" was not found in file "//&
        trim(filename))
    return
  endif

  status = nf90_inquire_variable(ncid, varid, ndims = ndims)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
      trim(NF90_STRERROR(status))//" Getting number of dimensions of "//&
      trim(varname)//" in "//trim(filename))
    return
  endif

  if (present(min_dims)) then
    if (ndims < min_dims-1) then
      write(msg, '(I3)') min_dims
      call MOM_error(WARNING, "num_timelevels: variable "//trim(varname)//&
        " in file "//trim(filename)//" has fewer than min_dims = "//trim(msg)//&
        " dimensions.")
    elseif (ndims == min_dims - 1) then
      n_time = 0 ; return
    endif
  endif

  status = nf90_inquire_variable(ncid, varid, dimids = dimids(1:ndims))
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
      trim(NF90_STRERROR(status))//" Getting last dimension ID for "//&
      trim(varname)//" in "//trim(filename))
    return
  endif

  status = nf90_Inquire_Dimension(ncid, dimids(ndims), len=n_time)
  if (status /= NF90_NOERR) call MOM_error(WARNING,"num_timelevels: "//&
      trim(NF90_STRERROR(status))//" Getting number of time levels of "//&
      trim(varname)//" in "//trim(filename))

  return

end function num_timelevels


function slasher(dir)
  character(len=*), intent(in) :: dir
  character(len=len(dir)) :: slasher

  if (len_trim(dir) == 0) then
    if (len(dir) < 2) call MOM_error(FATAL, &
        "Argument to MOM_io slasher must be at least two characters long.")
    slasher = "./"
  elseif (dir(len_trim(dir):len_trim(dir)) == '/') then
    slasher = trim(dir)
  else
    if (len_trim(dir) == len(dir)) call MOM_error(FATAL, &
        "Argument too short for MOM_io slasher to add needed slash to "//dir)
    slasher = trim(dir)//"/"
  endif
end function slasher

function MOM_file_exists(file_name, MOM_Domain)
  character(len=*),       intent(in) :: file_name
  type(MOM_domain_type),  intent(in) :: MOM_domain

!   This function uses the fms_io function file_exist to determine whether
! a named file (or its decomposed variant) exists.

  logical :: MOM_file_exists

  MOM_file_exists = file_exist(file_name, MOM_Domain%mpp_domain)

end function MOM_file_exists

subroutine MOM_read_data_1d(filename, fieldname, data)
  character(len=*),                 intent(in)    :: filename, fieldname
  real, dimension(:),               intent(inout) :: data ! 1 dimensional data

!   This function uses the fms_io function read_data to read 1-D
!  data field named "fieldname" from file "filename".

  call read_data(filename, fieldname, data)

end subroutine MOM_read_data_1d

subroutine MOM_read_data_2d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position)
  character(len=*),                 intent(in)    :: filename, fieldname
  real, dimension(:,:),             intent(inout) :: data ! 2 dimensional data
  type(MOM_domain_type),            intent(in)    :: MOM_Domain
  integer,                optional, intent(in)    :: timelevel, position

!   This function uses the fms_io function read_data to read a distributed
! 2-D data field named "fieldname" from file "filename".  Valid values for
! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.

  call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=position)

end subroutine MOM_read_data_2d

subroutine MOM_read_data_3d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position)
  character(len=*),                 intent(in)    :: filename, fieldname
  real, dimension(:,:,:),           intent(inout) :: data ! 2 dimensional data    
  type(MOM_domain_type),            intent(in)    :: MOM_Domain
  integer,                optional, intent(in)    :: timelevel, position

!   This function uses the fms_io function read_data to read a distributed
! 2-D data field named "fieldname" from file "filename".  Valid values for
! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.

  call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=position)

end subroutine MOM_read_data_3d

subroutine MOM_io_init(param_file)
  type(param_file_type), intent(in) :: param_file
! Argument:  param_file - A structure indicating the open file to parse for
!                         model parameter values.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_io" ! This module's name.

  call log_version(param_file, mod, version)

end subroutine MOM_io_init

end module MOM_io
