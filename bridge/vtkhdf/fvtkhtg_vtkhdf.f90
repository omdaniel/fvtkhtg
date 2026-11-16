!> fvtkhtg_vtkhdf: optional bridge from an fvtkhtg grid to the fvtkhdf
!> HyperTreeGrid writer (`vtkhdf_htg_file`).
!>
!> This is the only fvtkhtg code that depends on fvtkhdf (and through it on
!> HDF5 and, in MPI builds, MPI). The core `fvtkhtg` library stays pure
!> Fortran; fvtkhdf knows nothing about fvtkhtg. The bridge only calls public
!> methods of both sides:
!>
!>   fvtkhtg                       bridge                      fvtkhdf
!>   grid%export_breadth_first --> write_topology ----------> write_mesh_topology
!>   bf%gather(field, values) ---> write_cell_data ---------> write_cell_data
!>
!> The bridge never calls `create`, `start_time_step`, `finalize_time_step` or
!> `close`; the caller owns the file and its time stepping. Every routine is
!> collective in MPI builds of fvtkhdf: all ranks call it, in the same order,
!> with grids that agree on dimensions, branch factor, root indexing and the
!> names, order, types and component counts of their CellData arrays. Each
!> rank passes its own trees (no tree may be split between ranks).
!>
!> Errors detected on the fvtkhtg side are returned in `status`; the fvtkhdf
!> writer reports its own errors by aborting.
!>
!> Requires an fvtkhdf whose HyperTreeGrid writer accepts the interface array
!> names on write_mesh. The whole path is exercised through real files, read
!> back with VTK's vtkHDFReader, by test/vtkhdf.
module fvtkhtg_vtkhdf
   use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64
   use fvtkhtg, only: htg_grid_t, htg_breadth_first_t, htg_field_array_t, &
      htg_status_t, htg_set_ok, htg_fail, id_k, i32_k, &
      HTG_TYPE_REAL64, HTG_TYPE_REAL32, HTG_TYPE_INT64, HTG_TYPE_INT32, HTG_TYPE_BIT, &
      HTG_ERR_INVALID_ARG, HTG_ERR_INVALID_STATE, HTG_ERR_INTERFACE
   use vtkhdf_htg_file_type, only: vtkhdf_htg_file, vtkhdf_cell_data_handle
   implicit none
   private

   public :: htg_vtkhdf_write_mesh
   public :: htg_vtkhdf_write_geometry
   public :: htg_vtkhdf_write_topology
   public :: htg_vtkhdf_write_cell_data
   public :: htg_vtkhdf_write_all_cell_data
   public :: htg_vtkhdf_register_temporal_cell_data
   public :: htg_vtkhdf_write_temporal_cell_data
   public :: htg_vtkhdf_write_static

contains

   !> Define the root grid: point dimensions, branch factor, root indexing
   !> and, when the grid has a material interface, the names of its normals
   !> and intercepts arrays. Call once per file.
   !>
   !> The interface names are stored as file attributes that VTK matches
   !> against cell array names, so they must be names the writer stores
   !> unchanged: non-blank, without '/', '.' or spaces.
   subroutine htg_vtkhdf_write_mesh(file, grid, status)
      type(vtkhdf_htg_file), intent(inout) :: file
      type(htg_grid_t), intent(in) :: grid
      type(htg_status_t), intent(out), optional :: status

      character(len=:), allocatable :: normals, intercepts

      if (.not. grid%is_dimensions_set()) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_vtkhdf_write_mesh: grid dimensions are not set', status)
         return
      end if

      if (grid%get_has_interface()) then
         normals = grid%get_interface_normals_name()
         intercepts = grid%get_interface_intercepts_name()
         if (.not. (stored_unchanged(normals) .and. stored_unchanged(intercepts))) then
            call htg_fail(HTG_ERR_INTERFACE, &
               'htg_vtkhdf_write_mesh: interface array names must be non-blank and ' // &
               'contain no "/", "." or space', status)
            return
         end if
         call file%write_mesh(int(grid%get_dimensions()), &
            branch_factor=int(grid%get_branch_factor()), &
            transposed_root_indexing=grid%get_transposed_root_indexing(), &
            interface_normals_name=normals, interface_intercepts_name=intercepts)
      else
         call file%write_mesh(int(grid%get_dimensions()), &
            branch_factor=int(grid%get_branch_factor()), &
            transposed_root_indexing=grid%get_transposed_root_indexing())
      end if
      call htg_set_ok(status)
   end subroutine htg_vtkhdf_write_mesh

   !> Would fvtkhdf store this cell array name as it is?
   pure logical function stored_unchanged(name)
      character(len=*), intent(in) :: name
      stored_unchanged = (len_trim(name) > 0) .and. (scan(trim(name), '/. ') == 0)
   end function stored_unchanged

   !> Write the root grid's coordinate vectors (only rank zero's are stored).
   subroutine htg_vtkhdf_write_geometry(file, grid, status)
      type(vtkhdf_htg_file), intent(inout) :: file
      type(htg_grid_t), intent(in) :: grid
      type(htg_status_t), intent(out), optional :: status

      if (.not. (grid%has_x_coords() .and. grid%has_y_coords() .and. grid%has_z_coords())) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_vtkhdf_write_geometry: grid coordinates are incomplete', status)
         return
      end if
      call file%write_mesh_geometry(grid%x_coordinates, grid%y_coordinates, grid%z_coordinates)
      call htg_set_ok(status)
   end subroutine htg_vtkhdf_write_geometry

   !> Export the grid's trees and write them. `bf` is returned so cell data
   !> can be gathered in the same order. Pass `force_mask = .true.` on every
   !> rank when any rank's grid has a mask: fvtkhdf requires all ranks to
   !> supply a mask or none.
   subroutine htg_vtkhdf_write_topology(file, grid, bf, status, force_mask)
      type(vtkhdf_htg_file), intent(inout) :: file
      type(htg_grid_t), intent(in) :: grid
      type(htg_breadth_first_t), intent(out) :: bf
      type(htg_status_t), intent(out), optional :: status
      logical, intent(in), optional :: force_mask

      type(htg_status_t) :: st

      call grid%export_breadth_first(bf, st, force_mask=force_mask)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      ! bf%mask is unallocated when no mask is exported, which Fortran passes
      ! as an absent optional argument.
      call file%write_mesh_topology(bf%tree_ids, bf%xtree, bf%refined, bf%mask)
      call htg_set_ok(status)
   end subroutine htg_vtkhdf_write_topology

   !> Write one CellData array as a static cell dataset (modes STATIC,
   !> FIXED_MESH and MOVING_MESH). `name` defaults to the array's own name.
   subroutine htg_vtkhdf_write_cell_data(file, bf, field, status, name, attribute)
      type(vtkhdf_htg_file), intent(inout) :: file
      type(htg_breadth_first_t), intent(in) :: bf
      type(htg_field_array_t), intent(in) :: field
      type(htg_status_t), intent(out), optional :: status
      character(len=*), intent(in), optional :: name
      character(len=*), intent(in), optional :: attribute

      character(len=:), allocatable :: dataset_name

      if (present(name)) then
         dataset_name = name
      else
         dataset_name = field%get_name()
      end if
      call write_field(file, bf, field, status, dataset_name=dataset_name, attribute=attribute)
   end subroutine htg_vtkhdf_write_cell_data

   !> Write every CellData array of the grid as a static cell dataset. The
   !> active scalars array, if any, gets the VTK attribute role 'Scalars'.
   subroutine htg_vtkhdf_write_all_cell_data(file, grid, bf, status)
      type(vtkhdf_htg_file), intent(inout) :: file
      type(htg_grid_t), intent(in) :: grid
      type(htg_breadth_first_t), intent(in) :: bf
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st
      integer(i32_k) :: i
      character(len=:), allocatable :: active

      active = grid%cell_data%get_active_scalars()
      do i = 0_i32_k, grid%cell_data%get_number_of_arrays() - 1_i32_k
         associate (field => grid%cell_data%arrays(i)%item)
            if (len(active) > 0 .and. field%get_name() == active) then
               call htg_vtkhdf_write_cell_data(file, bf, field, st, attribute='Scalars')
            else
               call htg_vtkhdf_write_cell_data(file, bf, field, st)
            end if
         end associate
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
      end do
      call htg_set_ok(status)
   end subroutine htg_vtkhdf_write_all_cell_data

   !> Register a temporal cell dataset whose kind and component count match
   !> a CellData array. Call before the file's first start_time_step. BIT
   !> arrays are registered as int8.
   subroutine htg_vtkhdf_register_temporal_cell_data(file, field, handle, status, name, attribute)
      type(vtkhdf_htg_file), intent(inout) :: file
      type(htg_field_array_t), intent(in) :: field
      type(vtkhdf_cell_data_handle), intent(out) :: handle
      type(htg_status_t), intent(out), optional :: status
      character(len=*), intent(in), optional :: name
      character(len=*), intent(in), optional :: attribute

      character(len=:), allocatable :: dataset_name
      integer :: nc

      if (.not. field%is_initialized()) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_vtkhdf_register_temporal_cell_data: field array is not initialized', status)
         return
      end if
      if (present(name)) then
         dataset_name = name
      else
         dataset_name = field%get_name()
      end if
      nc = int(field%get_number_of_components())

      select case (field%get_data_type())
      case (HTG_TYPE_REAL64)
         if (nc == 1) then
            handle = file%register_temporal_cell_data(dataset_name, 0.0_real64, attribute)
         else
            handle = file%register_temporal_cell_data(dataset_name, spread(0.0_real64, 1, nc), attribute)
         end if
      case (HTG_TYPE_REAL32)
         if (nc == 1) then
            handle = file%register_temporal_cell_data(dataset_name, 0.0_real32, attribute)
         else
            handle = file%register_temporal_cell_data(dataset_name, spread(0.0_real32, 1, nc), attribute)
         end if
      case (HTG_TYPE_INT64)
         if (nc == 1) then
            handle = file%register_temporal_cell_data(dataset_name, 0_int64, attribute)
         else
            handle = file%register_temporal_cell_data(dataset_name, spread(0_int64, 1, nc), attribute)
         end if
      case (HTG_TYPE_INT32)
         if (nc == 1) then
            handle = file%register_temporal_cell_data(dataset_name, 0_int32, attribute)
         else
            handle = file%register_temporal_cell_data(dataset_name, spread(0_int32, 1, nc), attribute)
         end if
      case (HTG_TYPE_BIT)
         if (nc == 1) then
            handle = file%register_temporal_cell_data(dataset_name, 0_int8, attribute)
         else
            handle = file%register_temporal_cell_data(dataset_name, spread(0_int8, 1, nc), attribute)
         end if
      case default
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_vtkhdf_register_temporal_cell_data: unsupported data type', status)
         return
      end select
      call htg_set_ok(status)
   end subroutine htg_vtkhdf_register_temporal_cell_data

   !> Write one CellData array to a registered temporal cell dataset for the
   !> current time step.
   subroutine htg_vtkhdf_write_temporal_cell_data(file, handle, bf, field, status)
      type(vtkhdf_htg_file), intent(inout) :: file
      type(vtkhdf_cell_data_handle), intent(in) :: handle
      type(htg_breadth_first_t), intent(in) :: bf
      type(htg_field_array_t), intent(in) :: field
      type(htg_status_t), intent(out), optional :: status

      call write_field(file, bf, field, status, handle=handle)
   end subroutine htg_vtkhdf_write_temporal_cell_data

   !> Write a complete static HyperTreeGrid: validate the grid, then mesh,
   !> geometry, topology and every CellData array. The file must have been
   !> created in VTKHDF_STATIC (or FIXED_MESH / MOVING_MESH) mode.
   subroutine htg_vtkhdf_write_static(file, grid, status, force_mask, validate)
      type(vtkhdf_htg_file), intent(inout) :: file
      type(htg_grid_t), intent(in) :: grid
      type(htg_status_t), intent(out), optional :: status
      logical, intent(in), optional :: force_mask
      logical, intent(in), optional :: validate

      type(htg_breadth_first_t) :: bf
      type(htg_status_t) :: st
      logical :: do_validate

      do_validate = .true.
      if (present(validate)) do_validate = validate
      if (do_validate) then
         call grid%validate(st)
         if (.not. st%ok()) then
            call htg_fail(st%code, 'htg_vtkhdf_write_static: ' // trim(st%message), status)
            return
         end if
      end if

      call htg_vtkhdf_write_mesh(file, grid, st)
      if (st%ok()) call htg_vtkhdf_write_geometry(file, grid, st)
      if (st%ok()) call htg_vtkhdf_write_topology(file, grid, bf, st, force_mask=force_mask)
      if (st%ok()) call htg_vtkhdf_write_all_cell_data(file, grid, bf, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      call htg_set_ok(status)
   end subroutine htg_vtkhdf_write_static

   ! --------------------------------------------------------------------------
   ! Gather a field into export order and hand it to the writer, either as a
   ! new static dataset (dataset_name) or to a temporal handle.
   ! --------------------------------------------------------------------------
   subroutine write_field(file, bf, field, status, dataset_name, attribute, handle)
      type(vtkhdf_htg_file), intent(inout) :: file
      type(htg_breadth_first_t), intent(in) :: bf
      type(htg_field_array_t), intent(in) :: field
      type(htg_status_t), intent(out), optional :: status
      character(len=*), intent(in), optional :: dataset_name
      character(len=*), intent(in), optional :: attribute
      type(vtkhdf_cell_data_handle), intent(in), optional :: handle

      type(htg_status_t) :: st
      logical :: scalar
      real(real64), allocatable :: r64_1(:), r64_2(:, :)
      real(real32), allocatable :: r32_1(:), r32_2(:, :)
      integer(int64), allocatable :: i64_1(:), i64_2(:, :)
      integer(int32), allocatable :: i32_1(:), i32_2(:, :)
      integer(int8), allocatable :: i8_1(:), i8_2(:, :)

      if (.not. field%is_initialized()) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'fvtkhtg_vtkhdf: field array is not initialized', status)
         return
      end if
      scalar = (field%get_number_of_components() == 1_i32_k)

      select case (field%get_data_type())
      case (HTG_TYPE_REAL64)
         if (scalar) then
            call bf%gather(field, r64_1, st)
            if (st%ok()) call put_r64_1(r64_1)
         else
            call bf%gather(field, r64_2, st)
            if (st%ok()) call put_r64_2(r64_2)
         end if
      case (HTG_TYPE_REAL32)
         if (scalar) then
            call bf%gather(field, r32_1, st)
            if (st%ok()) call put_r32_1(r32_1)
         else
            call bf%gather(field, r32_2, st)
            if (st%ok()) call put_r32_2(r32_2)
         end if
      case (HTG_TYPE_INT64)
         if (scalar) then
            call bf%gather(field, i64_1, st)
            if (st%ok()) call put_i64_1(i64_1)
         else
            call bf%gather(field, i64_2, st)
            if (st%ok()) call put_i64_2(i64_2)
         end if
      case (HTG_TYPE_INT32)
         if (scalar) then
            call bf%gather(field, i32_1, st)
            if (st%ok()) call put_i32_1(i32_1)
         else
            call bf%gather(field, i32_2, st)
            if (st%ok()) call put_i32_2(i32_2)
         end if
      case (HTG_TYPE_BIT)
         if (scalar) then
            call bf%gather(field, i8_1, st)
            if (st%ok()) call put_i8_1(i8_1)
         else
            call bf%gather(field, i8_2, st)
            if (st%ok()) call put_i8_2(i8_2)
         end if
      case default
         call st%set(HTG_ERR_INVALID_ARG, 'fvtkhtg_vtkhdf: unsupported field data type')
      end select

      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      call htg_set_ok(status)

   contains

      ! One small writer call per kind and rank keeps generic resolution on
      ! the fvtkhdf side unambiguous.
      subroutine put_r64_1(v)
         real(real64), intent(in) :: v(:)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_r64_1
      subroutine put_r64_2(v)
         real(real64), intent(in) :: v(:, :)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_r64_2
      subroutine put_r32_1(v)
         real(real32), intent(in) :: v(:)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_r32_1
      subroutine put_r32_2(v)
         real(real32), intent(in) :: v(:, :)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_r32_2
      subroutine put_i64_1(v)
         integer(int64), intent(in) :: v(:)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_i64_1
      subroutine put_i64_2(v)
         integer(int64), intent(in) :: v(:, :)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_i64_2
      subroutine put_i32_1(v)
         integer(int32), intent(in) :: v(:)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_i32_1
      subroutine put_i32_2(v)
         integer(int32), intent(in) :: v(:, :)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_i32_2
      subroutine put_i8_1(v)
         integer(int8), intent(in) :: v(:)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_i8_1
      subroutine put_i8_2(v)
         integer(int8), intent(in) :: v(:, :)
         if (present(handle)) then
            call file%write_cell_data(handle, v)
         else
            call file%write_cell_data(dataset_name, v, attribute)
         end if
      end subroutine put_i8_2

   end subroutine write_field

end module fvtkhtg_vtkhdf
