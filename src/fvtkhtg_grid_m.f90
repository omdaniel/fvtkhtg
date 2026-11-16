module fvtkhtg_grid_m
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use fvtkhtg_kinds_m, only: i32_k, i64_k, id_k, r32_k, r64_k, &
      htg_status_t, HTG_ERR_INVALID_ARG, HTG_ERR_OUT_OF_BOUNDS, &
      HTG_ERR_INVALID_STATE, HTG_ERR_INDEXING, HTG_ERR_NOT_FOUND, &
      HTG_INVALID_ID, HTG_INDEX_UNSET, HTG_TYPE_REAL64, HTG_TYPE_REAL32, &
      HTG_TYPE_INT64, HTG_TYPE_INT32, htg_set_ok, htg_fail
   use fvtkhtg_bit_array_m, only: htg_bit_array_t
   use fvtkhtg_cell_data_m, only: htg_cell_data_t, htg_field_array_t
   use fvtkhtg_tree_m, only: htg_tree_t
   use fvtkhtg_breadth_first_m, only: htg_breadth_first_t
   implicit none
   private

   integer(i32_k), parameter, public :: HTG_AXIS_X = 0_i32_k
   integer(i32_k), parameter, public :: HTG_AXIS_Y = 1_i32_k
   integer(i32_k), parameter, public :: HTG_AXIS_Z = 2_i32_k
   integer(i32_k), parameter, public :: HTG_INVALID_ORIENTATION = -1_i32_k

   integer(id_k), parameter :: MIN_TREE_SLOT_GROWTH = 8_id_k

   ! Boxed sparse tree entry sorted by 0-based tree_id.
   ! Using an allocatable item moved via move_alloc preserves heap addresses
   ! of existing htg_tree_t instances across container growth and insertions.
   type, public :: htg_tree_slot_t
      integer(id_k) :: tree_id = HTG_INVALID_ID
      type(htg_tree_t), allocatable :: tree
   end type htg_tree_slot_t

   type, public :: htg_grid_t
      logical :: dimensions_set = .false.
      integer(i32_k) :: dimensions(0:2) = [0_i32_k, 0_i32_k, 0_i32_k]
      integer(i32_k) :: cell_dims(0:2) = [0_i32_k, 0_i32_k, 0_i32_k]
      integer(i32_k) :: dimension = 0_i32_k
      integer(i32_k) :: orientation = HTG_INVALID_ORIENTATION
      integer(i32_k) :: axes(0:2) = [-1_i32_k, -1_i32_k, -1_i32_k]
      integer(i32_k) :: branch_factor = 2_i32_k
      integer(i32_k) :: number_of_children = 0_i32_k
      logical :: transposed_root_indexing = .false.
      character(len=:), allocatable :: mode_squeeze
      integer(i32_k) :: depth_limiter = huge(0_i32_k)

      logical :: has_x_coordinates = .false.
      logical :: has_y_coordinates = .false.
      logical :: has_z_coordinates = .false.
      logical :: explicit_x_coordinates = .false.
      logical :: explicit_y_coordinates = .false.
      logical :: explicit_z_coordinates = .false.
      real(r64_k), allocatable :: x_coordinates(:)
      real(r64_k), allocatable :: y_coordinates(:)
      real(r64_k), allocatable :: z_coordinates(:)

      integer(id_k) :: num_trees = 0_id_k
      integer(id_k) :: tree_slot_capacity = 0_id_k
      type(htg_tree_slot_t), allocatable :: tree_slots(:)

      logical :: mask_active = .false.
      type(htg_bit_array_t) :: mask
      logical :: pure_mask_valid = .false.
      type(htg_bit_array_t) :: pure_mask

      type(htg_cell_data_t) :: cell_data
      logical :: has_interface = .false.
      character(len=:), allocatable :: interface_normals_name
      character(len=:), allocatable :: interface_intercepts_name
   contains
      procedure :: initialize => grid_initialize
      procedure :: clear => grid_clear

      procedure, private :: grid_set_dimensions_arr_i32
      procedure, private :: grid_set_dimensions_arr_i64
      procedure, private :: grid_set_dimensions_ijk_i32
      procedure, private :: grid_set_dimensions_ijk_i64
      generic :: set_dimensions => grid_set_dimensions_arr_i32, &
         grid_set_dimensions_arr_i64, grid_set_dimensions_ijk_i32, &
         grid_set_dimensions_ijk_i64
      generic :: set_dimensions_ijk => grid_set_dimensions_ijk_i32, &
         grid_set_dimensions_ijk_i64

      procedure :: is_dimensions_set => grid_is_dimensions_set
      procedure :: get_dimensions => grid_get_dimensions
      procedure :: get_cell_dims => grid_get_cell_dims
      procedure :: get_grid_sizes => grid_get_cell_dims
      procedure :: get_dimension => grid_get_dimension
      procedure :: get_orientation => grid_get_orientation
      procedure :: get_axes => grid_get_axes
      procedure :: get_1d_axis => grid_get_1d_axis
      procedure :: get_2d_axes => grid_get_2d_axes

      procedure :: set_mode_squeeze => grid_set_mode_squeeze
      procedure :: get_mode_squeeze => grid_get_mode_squeeze
      procedure :: set_depth_limiter => grid_set_depth_limiter
      procedure :: get_depth_limiter => grid_get_depth_limiter

      procedure :: set_branch_factor => grid_set_branch_factor
      procedure :: get_branch_factor => grid_get_branch_factor
      procedure :: get_number_of_children => grid_get_number_of_children

      procedure :: set_transposed_root_indexing => grid_set_transposed_root_indexing
      procedure :: get_transposed_root_indexing => grid_get_transposed_root_indexing

      procedure :: set_x_coordinates => grid_set_x_coordinates
      procedure :: set_y_coordinates => grid_set_y_coordinates
      procedure :: set_z_coordinates => grid_set_z_coordinates
      procedure :: has_x_coords => grid_has_x_coords
      procedure :: has_y_coords => grid_has_y_coords
      procedure :: has_z_coords => grid_has_z_coords
      procedure :: get_x_coordinates_view => grid_get_x_coordinates_view
      procedure :: get_y_coordinates_view => grid_get_y_coordinates_view
      procedure :: get_z_coordinates_view => grid_get_z_coordinates_view

      procedure, private :: grid_get_index_from_l0_i32
      procedure, private :: grid_get_index_from_l0_i64
      generic :: get_index_from_level_zero_coordinates => &
         grid_get_index_from_l0_i32, grid_get_index_from_l0_i64

      procedure, private :: grid_get_l0_from_index_i32
      procedure, private :: grid_get_l0_from_index_i64
      generic :: get_level_zero_coordinates_from_index => &
         grid_get_l0_from_index_i32, grid_get_l0_from_index_i64

      procedure :: get_active_axes => grid_get_active_axes
      procedure :: get_child_position => grid_get_child_position
      procedure :: get_level_zero_origin_and_size_from_index => grid_get_l0_origin_and_size
      procedure :: get_level_zero_origin_from_index => grid_get_l0_origin

      procedure :: get_max_number_of_trees => grid_get_max_number_of_trees
      procedure :: get_number_of_non_empty_trees => grid_get_number_of_non_empty_trees
      procedure :: has_tree => grid_has_tree
      procedure :: get_tree => grid_get_tree
      procedure :: set_tree => grid_set_tree
      procedure :: remove_tree => grid_remove_tree
      procedure :: get_tree_by_pos => grid_get_tree_by_pos

      procedure :: get_number_of_vertices => grid_get_number_of_vertices
      procedure :: get_number_of_cells => grid_get_number_of_vertices
      procedure :: get_number_of_leaves => grid_get_number_of_leaves
      procedure :: get_number_of_levels => grid_get_number_of_levels
      procedure :: get_max_global_index => grid_get_max_global_index
      procedure :: get_global_node_index_max => grid_get_global_node_index_max
      procedure :: initialize_local_index_node => grid_initialize_local_index_node
      procedure :: shift_all_trees_implicit_global_index => grid_initialize_local_index_node

      procedure :: set_mask => grid_set_mask
      procedure :: clear_mask => grid_clear_mask
      procedure :: has_mask => grid_has_mask
      procedure :: get_mask => grid_get_mask
      procedure :: set_cell_mask => grid_set_cell_mask
      procedure :: is_cell_masked => grid_is_cell_masked
      procedure :: invalidate_pure_mask => grid_invalidate_pure_mask
      procedure :: is_pure_mask_valid => grid_is_pure_mask_valid
      procedure :: compute_pure_mask_into => grid_compute_pure_mask_into
      procedure :: get_pure_mask => grid_get_pure_mask

      procedure :: get_cell_data => grid_get_cell_data
      procedure :: copy_cell_attributes => grid_copy_cell_attributes
      procedure :: clear_cell_attributes => grid_clear_cell_attributes
      procedure :: set_has_interface => grid_set_has_interface
      procedure :: get_has_interface => grid_get_has_interface
      procedure :: set_interface_normals_name => grid_set_interface_normals_name
      procedure :: get_interface_normals_name => grid_get_interface_normals_name
      procedure :: set_interface_intercepts_name => grid_set_interface_intercepts_name
      procedure :: get_interface_intercepts_name => grid_get_interface_intercepts_name

      procedure :: validate => grid_validate
      procedure :: export_breadth_first => grid_export_breadth_first
   end type htg_grid_t

   public :: htg_compute_squeezed_geometry

   interface
      module subroutine grid_validate(self, status)
         class(htg_grid_t), intent(in) :: self
         type(htg_status_t), intent(out), optional :: status
      end subroutine grid_validate

      !> Export the forest in breadth-first order for a serializer (see
      !> fvtkhtg_breadth_first_m for the index convention). Implemented in
      !> submodule fvtkhtg_grid_export_s.
      module subroutine grid_export_breadth_first(self, bf, status, force_mask)
         class(htg_grid_t), intent(in) :: self
         type(htg_breadth_first_t), intent(out) :: bf
         type(htg_status_t), intent(out), optional :: status
         logical, intent(in), optional :: force_mask
      end subroutine grid_export_breadth_first
   end interface

contains

   pure subroutine htg_compute_squeezed_geometry(dims, branch_factor, &
         cell_dims, dimension, orientation, axes, num_children)
      integer(i32_k), intent(in) :: dims(0:2)
      integer(i32_k), intent(in) :: branch_factor
      integer(i32_k), intent(out) :: cell_dims(0:2)
      integer(i32_k), intent(out) :: dimension
      integer(i32_k), intent(out) :: orientation
      integer(i32_k), intent(out) :: axes(0:2)
      integer(i32_k), intent(out) :: num_children

      integer(i32_k) :: a

      dimension = 0_i32_k
      orientation = HTG_INVALID_ORIENTATION
      axes = [-1_i32_k, -1_i32_k, -1_i32_k]

      do a = 0_i32_k, 2_i32_k
         if (dims(a) <= 1_i32_k) then
            cell_dims(a) = 1_i32_k
         else
            cell_dims(a) = dims(a) - 1_i32_k
            axes(dimension) = a
            dimension = dimension + 1_i32_k
         end if
      end do

      select case (dimension)
      case (0_i32_k)
         orientation = HTG_INVALID_ORIENTATION
         axes = [-1_i32_k, -1_i32_k, -1_i32_k]
         num_children = 1_i32_k
      case (1_i32_k)
         orientation = axes(0)
         axes = [orientation, -1_i32_k, -1_i32_k]
         num_children = branch_factor
      case (2_i32_k)
         do a = 0_i32_k, 2_i32_k
            if (dims(a) <= 1_i32_k) then
               orientation = a
               exit
            end if
         end do
         ! Cyclic right-handed plane axes matching VTK v9.7.2:
         ! YZ (orient=0) -> (1, 2); XZ (orient=1) -> (2, 0); XY (orient=2) -> (0, 1)
         axes(0) = mod(orientation + 1_i32_k, 3_i32_k)
         axes(1) = mod(orientation + 2_i32_k, 3_i32_k)
         axes(2) = -1_i32_k
         num_children = branch_factor * branch_factor
      case (3_i32_k)
         orientation = HTG_INVALID_ORIENTATION
         axes = [0_i32_k, 1_i32_k, 2_i32_k]
         num_children = branch_factor * branch_factor * branch_factor
      end select
   end subroutine htg_compute_squeezed_geometry

   pure subroutine grid_clear_trees_only(self)
      class(htg_grid_t), intent(inout) :: self
      integer(id_k) :: i

      if (allocated(self%tree_slots)) then
         do i = 0_id_k, self%num_trees - 1_id_k
            if (allocated(self%tree_slots(i)%tree)) then
               call self%tree_slots(i)%tree%clear()
               deallocate(self%tree_slots(i)%tree)
            end if
         end do
         deallocate(self%tree_slots)
      end if
      self%num_trees = 0_id_k
      self%tree_slot_capacity = 0_id_k
      self%pure_mask_valid = .false.
   end subroutine grid_clear_trees_only

   pure subroutine grid_clear(self)
      class(htg_grid_t), intent(inout) :: self

      call grid_clear_trees_only(self)

      self%dimensions_set = .false.
      self%dimensions = [0_i32_k, 0_i32_k, 0_i32_k]
      self%cell_dims = [0_i32_k, 0_i32_k, 0_i32_k]
      self%dimension = 0_i32_k
      self%orientation = HTG_INVALID_ORIENTATION
      self%axes = [-1_i32_k, -1_i32_k, -1_i32_k]
      self%branch_factor = 2_i32_k
      self%number_of_children = 0_i32_k
      self%transposed_root_indexing = .false.
      if (allocated(self%mode_squeeze)) deallocate(self%mode_squeeze)
      self%depth_limiter = huge(0_i32_k)

      self%has_x_coordinates = .false.
      self%has_y_coordinates = .false.
      self%has_z_coordinates = .false.
      self%explicit_x_coordinates = .false.
      self%explicit_y_coordinates = .false.
      self%explicit_z_coordinates = .false.
      if (allocated(self%x_coordinates)) deallocate(self%x_coordinates)
      if (allocated(self%y_coordinates)) deallocate(self%y_coordinates)
      if (allocated(self%z_coordinates)) deallocate(self%z_coordinates)

      self%mask_active = .false.
      call self%mask%clear()
      self%pure_mask_valid = .false.
      call self%pure_mask%clear()

      call self%cell_data%clear()
      self%has_interface = .false.
      if (allocated(self%interface_normals_name)) deallocate(self%interface_normals_name)
      if (allocated(self%interface_intercepts_name)) deallocate(self%interface_intercepts_name)
   end subroutine grid_clear

   pure subroutine grid_initialize(self)
      class(htg_grid_t), intent(inout) :: self
      call self%clear()
   end subroutine grid_initialize

   pure subroutine sync_squeezed_default_coord(dim_len, coords, has_flag, explicit_flag)
      integer(i32_k), intent(in) :: dim_len
      real(r64_k), allocatable, intent(inout) :: coords(:)
      logical, intent(inout) :: has_flag
      logical, intent(inout) :: explicit_flag

      if (dim_len == 1_i32_k) then
         if (.not. has_flag) then
            allocate(coords(0:0))
            coords(0) = 0.0_r64_k
            has_flag = .true.
            explicit_flag = .false.
         end if
      else
         if (has_flag .and. .not. explicit_flag) then
            if (allocated(coords)) deallocate(coords)
            has_flag = .false.
         end if
      end if
   end subroutine sync_squeezed_default_coord

   pure subroutine grid_set_dimensions_ijk_i32(self, ni, nj, nk, status)
      class(htg_grid_t), intent(inout) :: self
      integer(i32_k), intent(in) :: ni, nj, nk
      type(htg_status_t), intent(out), optional :: status

      integer(i32_k) :: d(0:2), cdims(0:2), ndim, orient, ax(0:2), nchild

      d = [ni, nj, nk]
      if (any(d < 1_i32_k)) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_dimensions: point dimensions must all be >= 1', status)
         return
      end if

      if (self%dimensions_set .and. all(self%dimensions == d)) then
         call htg_set_ok(status)
         return
      end if

      if (self%num_trees > 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%set_dimensions: cannot change dimensions while trees exist', &
            status)
         return
      end if

      call htg_compute_squeezed_geometry(d, self%branch_factor, &
         cdims, ndim, orient, ax, nchild)

      self%dimensions = d
      self%cell_dims = cdims
      self%dimension = ndim
      self%orientation = orient
      self%axes = ax
      self%number_of_children = nchild
      self%dimensions_set = .true.
      self%pure_mask_valid = .false.

      call sync_squeezed_default_coord(d(0), self%x_coordinates, &
         self%has_x_coordinates, self%explicit_x_coordinates)
      call sync_squeezed_default_coord(d(1), self%y_coordinates, &
         self%has_y_coordinates, self%explicit_y_coordinates)
      call sync_squeezed_default_coord(d(2), self%z_coordinates, &
         self%has_z_coordinates, self%explicit_z_coordinates)

      call htg_set_ok(status)
   end subroutine grid_set_dimensions_ijk_i32

   pure subroutine grid_set_dimensions_ijk_i64(self, ni, nj, nk, status)
      class(htg_grid_t), intent(inout) :: self
      integer(i64_k), intent(in) :: ni, nj, nk
      type(htg_status_t), intent(out), optional :: status

      if (ni < 1_i64_k .or. nj < 1_i64_k .or. nk < 1_i64_k .or. &
          ni > int(huge(0_i32_k), i64_k) .or. &
          nj > int(huge(0_i32_k), i64_k) .or. &
          nk > int(huge(0_i32_k), i64_k)) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_dimensions: dimension out of valid int32 range', status)
         return
      end if

      call self%grid_set_dimensions_ijk_i32( &
         int(ni, i32_k), int(nj, i32_k), int(nk, i32_k), status)
   end subroutine grid_set_dimensions_ijk_i64

   pure subroutine grid_set_dimensions_arr_i32(self, dims, status)
      class(htg_grid_t), intent(inout) :: self
      integer(i32_k), intent(in) :: dims(0:)
      type(htg_status_t), intent(out), optional :: status

      if (size(dims) /= 3) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_dimensions: dims array must have exactly 3 elements', status)
         return
      end if
      call self%grid_set_dimensions_ijk_i32(dims(0), dims(1), dims(2), status)
   end subroutine grid_set_dimensions_arr_i32

   pure subroutine grid_set_dimensions_arr_i64(self, dims, status)
      class(htg_grid_t), intent(inout) :: self
      integer(i64_k), intent(in) :: dims(0:)
      type(htg_status_t), intent(out), optional :: status

      if (size(dims) /= 3) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_dimensions: dims array must have exactly 3 elements', status)
         return
      end if
      call self%grid_set_dimensions_ijk_i64(dims(0), dims(1), dims(2), status)
   end subroutine grid_set_dimensions_arr_i64

   pure logical function grid_is_dimensions_set(self) result(flag)
      class(htg_grid_t), intent(in) :: self
      flag = self%dimensions_set
   end function grid_is_dimensions_set

   pure function grid_get_dimensions(self) result(dims)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k) :: dims(0:2)
      dims = self%dimensions
   end function grid_get_dimensions

   pure function grid_get_cell_dims(self) result(cdims)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k) :: cdims(0:2)
      cdims = self%cell_dims
   end function grid_get_cell_dims

   pure function grid_get_dimension(self) result(ndim)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k) :: ndim
      ndim = self%dimension
   end function grid_get_dimension

   pure function grid_get_orientation(self) result(orient)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k) :: orient
      orient = self%orientation
   end function grid_get_orientation

   pure function grid_get_axes(self) result(ax)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k) :: ax(0:2)
      ax = self%axes
   end function grid_get_axes

   pure subroutine grid_get_1d_axis(self, axis, status)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k), intent(out) :: axis
      type(htg_status_t), intent(out), optional :: status

      axis = -1_i32_k
      if (.not. self%dimensions_set .or. self%dimension /= 1_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%get_1d_axis: grid dimension is not 1', status)
         return
      end if
      axis = self%axes(0)
      call htg_set_ok(status)
   end subroutine grid_get_1d_axis

   pure subroutine grid_get_2d_axes(self, axis1, axis2, status)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k), intent(out) :: axis1, axis2
      type(htg_status_t), intent(out), optional :: status

      axis1 = -1_i32_k
      axis2 = -1_i32_k
      if (.not. self%dimensions_set .or. self%dimension /= 2_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%get_2d_axes: grid dimension is not 2', status)
         return
      end if
      axis1 = self%axes(0)
      axis2 = self%axes(1)
      call htg_set_ok(status)
   end subroutine grid_get_2d_axes

   pure subroutine grid_set_mode_squeeze(self, mode)
      class(htg_grid_t), intent(inout) :: self
      character(len=*), intent(in) :: mode
      self%mode_squeeze = trim(mode)
   end subroutine grid_set_mode_squeeze

   pure function grid_get_mode_squeeze(self) result(mode)
      class(htg_grid_t), intent(in) :: self
      character(len=:), allocatable :: mode
      if (allocated(self%mode_squeeze)) then
         mode = self%mode_squeeze
      else
         mode = ''
      end if
   end function grid_get_mode_squeeze

   pure subroutine grid_set_depth_limiter(self, limit, status)
      class(htg_grid_t), intent(inout) :: self
      integer(i32_k), intent(in) :: limit
      type(htg_status_t), intent(out), optional :: status

      if (limit < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_depth_limiter: depth_limiter must be >= 0', status)
         return
      end if
      if (self%depth_limiter /= limit) then
         self%depth_limiter = limit
         self%pure_mask_valid = .false.
      end if
      call htg_set_ok(status)
   end subroutine grid_set_depth_limiter

   pure function grid_get_depth_limiter(self) result(limit)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k) :: limit
      limit = self%depth_limiter
   end function grid_get_depth_limiter

   pure subroutine grid_set_branch_factor(self, factor, status)
      class(htg_grid_t), intent(inout) :: self
      integer(i32_k), intent(in) :: factor
      type(htg_status_t), intent(out), optional :: status

      integer(i32_k) :: cdims(0:2), ndim, orient, ax(0:2), nchild

      if (factor /= 2_i32_k .and. factor /= 3_i32_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_branch_factor: branch_factor must be 2 or 3', status)
         return
      end if

      if (self%branch_factor == factor) then
         call htg_set_ok(status)
         return
      end if

      if (self%num_trees > 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%set_branch_factor: cannot change while trees exist', &
            status)
         return
      end if

      self%branch_factor = factor
      if (self%dimensions_set) then
         call htg_compute_squeezed_geometry(self%dimensions, self%branch_factor, &
            cdims, ndim, orient, ax, nchild)
         self%number_of_children = nchild
      end if
      self%pure_mask_valid = .false.
      call htg_set_ok(status)
   end subroutine grid_set_branch_factor

   pure function grid_get_branch_factor(self) result(factor)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k) :: factor
      factor = self%branch_factor
   end function grid_get_branch_factor

   pure function grid_get_number_of_children(self) result(nchild)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k) :: nchild
      nchild = self%number_of_children
   end function grid_get_number_of_children

   pure subroutine grid_set_transposed_root_indexing(self, flag, status)
      class(htg_grid_t), intent(inout) :: self
      logical, intent(in) :: flag
      type(htg_status_t), intent(out), optional :: status

      if (self%transposed_root_indexing .eqv. flag) then
         call htg_set_ok(status)
         return
      end if

      if (self%num_trees > 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%set_transposed_root_indexing: cannot change after trees exist', status)
         return
      end if

      self%transposed_root_indexing = flag
      call htg_set_ok(status)
   end subroutine grid_set_transposed_root_indexing

   pure logical function grid_get_transposed_root_indexing(self) result(flag)
      class(htg_grid_t), intent(in) :: self
      flag = self%transposed_root_indexing
   end function grid_get_transposed_root_indexing

   pure subroutine validate_and_copy_coords(coords, expected_len, check_len, &
         axis_label, out_coords, out_has, status)
      real(r64_k), intent(in) :: coords(0:)
      integer(i32_k), intent(in) :: expected_len
      logical, intent(in) :: check_len
      character(len=*), intent(in) :: axis_label
      real(r64_k), allocatable, intent(inout) :: out_coords(:)
      logical, intent(inout) :: out_has
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: n, i
      logical :: increasing

      n = size(coords, kind=id_k)
      if (n < 1_id_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_' // axis_label // &
            '_coordinates: coordinate array cannot be empty', status)
         return
      end if

      if (check_len .and. n /= int(expected_len, id_k)) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_' // axis_label // &
            '_coordinates: length does not match dimensions', status)
         return
      end if

      do i = 0_id_k, n - 1_id_k
         if (.not. ieee_is_finite(coords(i))) then
            call htg_fail(HTG_ERR_INVALID_ARG, &
               'htg_grid_t%set_' // axis_label // &
               '_coordinates: non-finite coordinate value', status)
            return
         end if
      end do

      if (n >= 2_id_k) then
         if (.not. (coords(1) > coords(0)) .and. .not. (coords(1) < coords(0))) then
            call htg_fail(HTG_ERR_INVALID_ARG, &
               'htg_grid_t%set_' // axis_label // &
               '_coordinates: duplicate coordinates (not strictly monotonic)', status)
            return
         end if
         increasing = (coords(1) > coords(0))
         do i = 1_id_k, n - 1_id_k
            if (increasing) then
               if (.not. (coords(i) > coords(i - 1_id_k))) then
                  call htg_fail(HTG_ERR_INVALID_ARG, &
                     'htg_grid_t%set_' // axis_label // &
                     '_coordinates: coordinates must be strictly monotonic', status)
                  return
               end if
            else
               if (.not. (coords(i) < coords(i - 1_id_k))) then
                  call htg_fail(HTG_ERR_INVALID_ARG, &
                     'htg_grid_t%set_' // axis_label // &
                     '_coordinates: coordinates must be strictly monotonic', status)
                  return
               end if
            end if
         end do
      end if

      if (allocated(out_coords)) deallocate(out_coords)
      allocate(out_coords(0_id_k : n - 1_id_k))
      out_coords(0_id_k : n - 1_id_k) = coords(0_id_k : n - 1_id_k)
      out_has = .true.
      call htg_set_ok(status)
   end subroutine validate_and_copy_coords

   pure subroutine grid_set_x_coordinates(self, coords, status)
      class(htg_grid_t), intent(inout) :: self
      real(r64_k), intent(in) :: coords(0:)
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st

      call validate_and_copy_coords(coords, self%dimensions(0), &
         self%dimensions_set, 'x', self%x_coordinates, self%has_x_coordinates, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      self%explicit_x_coordinates = .true.
      call htg_set_ok(status)
   end subroutine grid_set_x_coordinates

   pure subroutine grid_set_y_coordinates(self, coords, status)
      class(htg_grid_t), intent(inout) :: self
      real(r64_k), intent(in) :: coords(0:)
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st

      call validate_and_copy_coords(coords, self%dimensions(1), &
         self%dimensions_set, 'y', self%y_coordinates, self%has_y_coordinates, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      self%explicit_y_coordinates = .true.
      call htg_set_ok(status)
   end subroutine grid_set_y_coordinates

   pure subroutine grid_set_z_coordinates(self, coords, status)
      class(htg_grid_t), intent(inout) :: self
      real(r64_k), intent(in) :: coords(0:)
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st

      call validate_and_copy_coords(coords, self%dimensions(2), &
         self%dimensions_set, 'z', self%z_coordinates, self%has_z_coordinates, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      self%explicit_z_coordinates = .true.
      call htg_set_ok(status)
   end subroutine grid_set_z_coordinates

   pure logical function grid_has_x_coords(self) result(flag)
      class(htg_grid_t), intent(in) :: self
      flag = self%has_x_coordinates .and. allocated(self%x_coordinates)
   end function grid_has_x_coords

   pure logical function grid_has_y_coords(self) result(flag)
      class(htg_grid_t), intent(in) :: self
      flag = self%has_y_coordinates .and. allocated(self%y_coordinates)
   end function grid_has_y_coords

   pure logical function grid_has_z_coords(self) result(flag)
      class(htg_grid_t), intent(in) :: self
      flag = self%has_z_coordinates .and. allocated(self%z_coordinates)
   end function grid_has_z_coords

   subroutine grid_get_x_coordinates_view(self, ptr)
      class(htg_grid_t), target, intent(inout) :: self
      real(r64_k), pointer, intent(out) :: ptr(:)
      integer(id_k) :: n

      if (.not. self%has_x_coordinates .or. .not. allocated(self%x_coordinates)) then
         ptr => null()
      else
         n = size(self%x_coordinates, kind=id_k)
         ptr(0_id_k : n - 1_id_k) => self%x_coordinates(0_id_k : n - 1_id_k)
      end if
   end subroutine grid_get_x_coordinates_view

   subroutine grid_get_y_coordinates_view(self, ptr)
      class(htg_grid_t), target, intent(inout) :: self
      real(r64_k), pointer, intent(out) :: ptr(:)
      integer(id_k) :: n

      if (.not. self%has_y_coordinates .or. .not. allocated(self%y_coordinates)) then
         ptr => null()
      else
         n = size(self%y_coordinates, kind=id_k)
         ptr(0_id_k : n - 1_id_k) => self%y_coordinates(0_id_k : n - 1_id_k)
      end if
   end subroutine grid_get_y_coordinates_view

   subroutine grid_get_z_coordinates_view(self, ptr)
      class(htg_grid_t), target, intent(inout) :: self
      real(r64_k), pointer, intent(out) :: ptr(:)
      integer(id_k) :: n

      if (.not. self%has_z_coordinates .or. .not. allocated(self%z_coordinates)) then
         ptr => null()
      else
         n = size(self%z_coordinates, kind=id_k)
         ptr(0_id_k : n - 1_id_k) => self%z_coordinates(0_id_k : n - 1_id_k)
      end if
   end subroutine grid_get_z_coordinates_view

   pure subroutine grid_get_index_from_l0_i64(self, tree_id, i, j, k, status)
      class(htg_grid_t), intent(in) :: self
      integer(id_k), intent(out) :: tree_id
      integer(id_k), intent(in) :: i, j, k
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: ncx, ncy, ncz

      tree_id = HTG_INVALID_ID
      if (.not. self%dimensions_set) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%get_index_from_level_zero_coordinates: dimensions not set', status)
         return
      end if

      ncx = int(self%cell_dims(0), id_k)
      ncy = int(self%cell_dims(1), id_k)
      ncz = int(self%cell_dims(2), id_k)

      if (i < 0_id_k .or. i >= ncx .or. &
          j < 0_id_k .or. j >= ncy .or. &
          k < 0_id_k .or. k >= ncz) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%get_index_from_level_zero_coordinates: (i,j,k) out of bounds', &
            status)
         return
      end if

      if (.not. self%transposed_root_indexing) then
         tree_id = i + j * ncx + k * ncx * ncy
      else
         tree_id = k + j * ncz + i * ncz * ncy
      end if

      call htg_set_ok(status)
   end subroutine grid_get_index_from_l0_i64

   pure subroutine grid_get_index_from_l0_i32(self, tree_id, i, j, k, status)
      class(htg_grid_t), intent(in) :: self
      integer(id_k), intent(out) :: tree_id
      integer(i32_k), intent(in) :: i, j, k
      type(htg_status_t), intent(out), optional :: status

      call self%grid_get_index_from_l0_i64(tree_id, &
         int(i, id_k), int(j, id_k), int(k, id_k), status)
   end subroutine grid_get_index_from_l0_i32

   pure subroutine grid_get_l0_from_index_i64(self, tree_id, i, j, k, status)
      class(htg_grid_t), intent(in) :: self
      integer(id_k), intent(in) :: tree_id
      integer(id_k), intent(out) :: i, j, k
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: ncx, ncy, ncz, max_trees, rem

      i = -1_id_k
      j = -1_id_k
      k = -1_id_k

      if (.not. self%dimensions_set) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%get_level_zero_coordinates_from_index: dimensions not set', status)
         return
      end if

      max_trees = self%get_max_number_of_trees()
      if (tree_id < 0_id_k .or. tree_id >= max_trees) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%get_level_zero_coordinates_from_index: tree_id out of bounds', status)
         return
      end if

      ncx = int(self%cell_dims(0), id_k)
      ncy = int(self%cell_dims(1), id_k)
      ncz = int(self%cell_dims(2), id_k)

      if (.not. self%transposed_root_indexing) then
         i = mod(tree_id, ncx)
         rem = tree_id / ncx
         j = mod(rem, ncy)
         k = rem / ncy
      else
         k = mod(tree_id, ncz)
         rem = tree_id / ncz
         j = mod(rem, ncy)
         i = rem / ncy
      end if

      call htg_set_ok(status)
   end subroutine grid_get_l0_from_index_i64

   pure subroutine grid_get_l0_from_index_i32(self, tree_id, i, j, k, status)
      class(htg_grid_t), intent(in) :: self
      integer(id_k), intent(in) :: tree_id
      integer(i32_k), intent(out) :: i, j, k
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: i64, j64, k64
      type(htg_status_t) :: st

      i = -1_i32_k
      j = -1_i32_k
      k = -1_i32_k
      call self%grid_get_l0_from_index_i64(tree_id, i64, j64, k64, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      i = int(i64, i32_k)
      j = int(j64, i32_k)
      k = int(k64, i32_k)
      call htg_set_ok(status)
   end subroutine grid_get_l0_from_index_i32

   ! --- Cell Geometry Support ---

   !> The grid's non-squeezed axes in ascending order (x before y before z):
   !> active(0 : nactive - 1), with -1 in the unused entries. This is the
   !> order child indices are built in: digit p of a child index (base
   !> branch_factor) is the child's position along active(p). It differs from
   !> get_axes() for a grid in the xz plane, whose axes VTK lists as (z, x)
   !> while still numbering children x first.
   pure subroutine grid_get_active_axes(self, active, nactive)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k), intent(out) :: active(0:2)
      integer(i32_k), intent(out) :: nactive

      integer(i32_k) :: a

      active = -1_i32_k
      nactive = 0_i32_k
      if (.not. self%dimensions_set) return
      do a = 0_i32_k, 2_i32_k
         if (self%dimensions(a) > 1_i32_k) then
            active(nactive) = a
            nactive = nactive + 1_i32_k
         end if
      end do
   end subroutine grid_get_active_axes

   !> Position of child `ichild` within its parent along physical axis `axis`:
   !> 0 .. branch_factor - 1, counted from the low-coordinate side. Returns -1
   !> if the axis is squeezed or an argument is out of range.
   pure function grid_get_child_position(self, ichild, axis) result(pos)
      class(htg_grid_t), intent(in) :: self
      integer(i32_k), intent(in) :: ichild
      integer(i32_k), intent(in) :: axis
      integer(i32_k) :: pos

      integer(i32_k) :: active(0:2), nactive, p, stride

      pos = -1_i32_k
      if (ichild < 0_i32_k .or. ichild >= self%number_of_children) return
      call self%get_active_axes(active, nactive)
      stride = 1_i32_k
      do p = 0_i32_k, nactive - 1_i32_k
         if (active(p) == axis) then
            pos = mod(ichild / stride, self%branch_factor)
            return
         end if
         stride = stride * self%branch_factor
      end do
   end function grid_get_child_position

   !> Origin (low corner) and edge lengths of root cell `tree_id`, matching
   !> vtkHyperTreeGrid::GetLevelZeroOriginAndSizeFromIndex: the origin is the
   !> coordinate at the root's (i, j, k); the size is the distance to the next
   !> coordinate, and 0 along a squeezed axis. With descending coordinates
   !> the size is negative, as in VTK.
   pure subroutine grid_get_l0_origin_and_size(self, tree_id, origin, cell_size, status)
      class(htg_grid_t), intent(in) :: self
      integer(id_k), intent(in) :: tree_id
      real(r64_k), intent(out) :: origin(0:2)
      real(r64_k), intent(out) :: cell_size(0:2)
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: i, j, k
      type(htg_status_t) :: st

      origin = 0.0_r64_k
      cell_size = 0.0_r64_k
      call self%grid_get_l0_from_index_i64(tree_id, i, j, k, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      if (.not. coords_match(self%x_coordinates, self%has_x_coordinates, self%dimensions(0)) .or. &
          .not. coords_match(self%y_coordinates, self%has_y_coordinates, self%dimensions(1)) .or. &
          .not. coords_match(self%z_coordinates, self%has_z_coordinates, self%dimensions(2))) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%get_level_zero_origin_and_size_from_index: coordinates are missing ' // &
            'or do not match the dimensions', status)
         return
      end if

      origin(0) = self%x_coordinates(i)
      origin(1) = self%y_coordinates(j)
      origin(2) = self%z_coordinates(k)
      if (self%dimensions(0) > 1_i32_k) cell_size(0) = self%x_coordinates(i + 1_id_k) - origin(0)
      if (self%dimensions(1) > 1_i32_k) cell_size(1) = self%y_coordinates(j + 1_id_k) - origin(1)
      if (self%dimensions(2) > 1_i32_k) cell_size(2) = self%z_coordinates(k + 1_id_k) - origin(2)
      call htg_set_ok(status)
   end subroutine grid_get_l0_origin_and_size

   pure subroutine grid_get_l0_origin(self, tree_id, origin, status)
      class(htg_grid_t), intent(in) :: self
      integer(id_k), intent(in) :: tree_id
      real(r64_k), intent(out) :: origin(0:2)
      type(htg_status_t), intent(out), optional :: status

      real(r64_k) :: cell_size(0:2)

      call self%get_level_zero_origin_and_size_from_index(tree_id, origin, cell_size, status)
   end subroutine grid_get_l0_origin

   pure logical function coords_match(coords, has_flag, dim_len) result(ok)
      real(r64_k), allocatable, intent(in) :: coords(:)
      logical, intent(in) :: has_flag
      integer(i32_k), intent(in) :: dim_len

      ok = .false.
      if (.not. has_flag .or. .not. allocated(coords)) return
      ok = (lbound(coords, 1) == 0) .and. (size(coords, kind=id_k) == int(dim_len, id_k))
   end function coords_match

   ! --- Sparse Forest Management ---

   pure function grid_get_max_number_of_trees(self) result(max_trees)
      class(htg_grid_t), intent(in) :: self
      integer(id_k) :: max_trees

      if (.not. self%dimensions_set) then
         max_trees = 0_id_k
      else
         max_trees = int(self%cell_dims(0), id_k) * &
                     int(self%cell_dims(1), id_k) * &
                     int(self%cell_dims(2), id_k)
      end if
   end function grid_get_max_number_of_trees

   pure function grid_get_number_of_non_empty_trees(self) result(ntrees)
      class(htg_grid_t), intent(in) :: self
      integer(id_k) :: ntrees
      ntrees = self%num_trees
   end function grid_get_number_of_non_empty_trees

   ! Binary search returning 0-based insertion/match position in [0, num_trees]
   pure subroutine lower_bound_tree_slot(self, tree_id, pos, found)
      class(htg_grid_t), intent(in) :: self
      integer(id_k), intent(in) :: tree_id
      integer(id_k), intent(out) :: pos
      logical, intent(out) :: found

      integer(id_k) :: lo, hi, mid

      found = .false.
      if (self%num_trees <= 0_id_k .or. .not. allocated(self%tree_slots)) then
         pos = 0_id_k
         return
      end if

      ! Fast path for monotonic ascending append
      if (tree_id > self%tree_slots(self%num_trees - 1_id_k)%tree_id) then
         pos = self%num_trees
         return
      end if

      lo = 0_id_k
      hi = self%num_trees
      do while (lo < hi)
         mid = lo + (hi - lo) / 2_id_k
         if (self%tree_slots(mid)%tree_id < tree_id) then
            lo = mid + 1_id_k
         else
            hi = mid
         end if
      end do
      pos = lo
      if (pos < self%num_trees) then
         found = (self%tree_slots(pos)%tree_id == tree_id)
      end if
   end subroutine lower_bound_tree_slot

   pure logical function grid_has_tree(self, tree_id) result(flag)
      class(htg_grid_t), intent(in) :: self
      integer(id_k), intent(in) :: tree_id

      integer(id_k) :: pos
      logical :: found

      flag = .false.
      if (.not. self%dimensions_set) return
      if (tree_id < 0_id_k .or. tree_id >= self%get_max_number_of_trees()) return
      call lower_bound_tree_slot(self, tree_id, pos, found)
      if (found) then
         flag = allocated(self%tree_slots(pos)%tree)
      end if
   end function grid_has_tree

   subroutine ensure_tree_slot_capacity(self, min_cap)
      class(htg_grid_t), target, intent(inout) :: self
      integer(id_k), intent(in) :: min_cap

      integer(id_k) :: new_cap, i
      type(htg_tree_slot_t), allocatable, target :: tmp(:)

      if (min_cap <= self%tree_slot_capacity .and. allocated(self%tree_slots)) return

      new_cap = max(min_cap, max(MIN_TREE_SLOT_GROWTH, 2_id_k * self%tree_slot_capacity))
      allocate(tmp(0_id_k : new_cap - 1_id_k))
      if (allocated(self%tree_slots)) then
         do i = 0_id_k, self%num_trees - 1_id_k
            tmp(i)%tree_id = self%tree_slots(i)%tree_id
            call move_alloc(self%tree_slots(i)%tree, tmp(i)%tree)
         end do
      end if
      call move_alloc(tmp, self%tree_slots)
      self%tree_slot_capacity = new_cap
   end subroutine ensure_tree_slot_capacity

   subroutine grid_get_tree(self, tree_id, tree_ptr, create, status)
      class(htg_grid_t), target, intent(inout) :: self
      integer(id_k), intent(in) :: tree_id
      type(htg_tree_t), pointer, intent(out) :: tree_ptr
      logical, intent(in), optional :: create
      type(htg_status_t), intent(out), optional :: status

      logical :: do_create, found
      integer(id_k) :: pos, i
      type(htg_tree_t), allocatable, target :: new_tree
      type(htg_status_t) :: init_st

      tree_ptr => null()
      if (.not. self%dimensions_set) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%get_tree: grid dimensions must be set first', status)
         return
      end if

      if (tree_id < 0_id_k .or. tree_id >= self%get_max_number_of_trees()) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%get_tree: tree_id out of bounds', status)
         return
      end if

      do_create = .false.
      if (present(create)) do_create = create

      call lower_bound_tree_slot(self, tree_id, pos, found)
      if (found) then
         tree_ptr => self%tree_slots(pos)%tree
         call htg_set_ok(status)
         return
      end if

      if (.not. do_create) then
         tree_ptr => null()
         call htg_set_ok(status)
         return
      end if

      if (self%dimension == 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%get_tree: cannot instantiate HyperTree on 0D (1,1,1) grid', &
            status)
         return
      end if

      allocate(new_tree)
      call new_tree%initialize(self%branch_factor, self%dimension, tree_id, init_st)
      if (.not. init_st%ok()) then
         deallocate(new_tree)
         call htg_fail(init_st%code, trim(init_st%message), status)
         return
      end if

      call ensure_tree_slot_capacity(self, self%num_trees + 1_id_k)
      do i = self%num_trees - 1_id_k, pos, -1_id_k
         self%tree_slots(i + 1_id_k)%tree_id = self%tree_slots(i)%tree_id
         call move_alloc(self%tree_slots(i)%tree, self%tree_slots(i + 1_id_k)%tree)
      end do

      self%tree_slots(pos)%tree_id = tree_id
      call move_alloc(new_tree, self%tree_slots(pos)%tree)
      self%num_trees = self%num_trees + 1_id_k
      self%pure_mask_valid = .false.

      tree_ptr => self%tree_slots(pos)%tree
      call htg_set_ok(status)
   end subroutine grid_get_tree

   subroutine grid_set_tree(self, tree_id, tree, status)
      class(htg_grid_t), target, intent(inout) :: self
      integer(id_k), intent(in) :: tree_id
      type(htg_tree_t), intent(in) :: tree
      type(htg_status_t), intent(out), optional :: status

      logical :: found
      integer(id_k) :: pos, i, max_gid
      type(htg_tree_t), allocatable, target :: copy_tree
      type(htg_status_t) :: cd_st

      if (.not. self%dimensions_set) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%set_tree: grid dimensions must be set first', status)
         return
      end if

      if (tree_id < 0_id_k .or. tree_id >= self%get_max_number_of_trees()) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%set_tree: tree_id out of bounds', status)
         return
      end if

      if (.not. tree%is_initialized()) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_tree: supplied tree is not initialized', status)
         return
      end if

      if (tree%get_branch_factor() /= self%branch_factor .or. &
          tree%get_dimension() /= self%dimension) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%set_tree: tree (branch_factor, dimension) mismatch with grid', status)
         return
      end if

      allocate(copy_tree, source=tree)
      copy_tree%tree_index = tree_id

      call lower_bound_tree_slot(self, tree_id, pos, found)
      if (found) then
         if (allocated(self%tree_slots(pos)%tree)) then
            call self%tree_slots(pos)%tree%clear()
            deallocate(self%tree_slots(pos)%tree)
         end if
         call move_alloc(copy_tree, self%tree_slots(pos)%tree)
      else
         call ensure_tree_slot_capacity(self, self%num_trees + 1_id_k)
         do i = self%num_trees - 1_id_k, pos, -1_id_k
            self%tree_slots(i + 1_id_k)%tree_id = self%tree_slots(i)%tree_id
            call move_alloc(self%tree_slots(i)%tree, self%tree_slots(i + 1_id_k)%tree)
         end do
         self%tree_slots(pos)%tree_id = tree_id
         call move_alloc(copy_tree, self%tree_slots(pos)%tree)
         self%num_trees = self%num_trees + 1_id_k
      end if

      if (self%cell_data%get_number_of_arrays() > 0_i32_k) then
         max_gid = self%tree_slots(pos)%tree%get_max_global_index()
         if (max_gid >= 0_id_k) then
            if (max_gid == huge(0_id_k)) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_grid_t%set_tree: max_global_index == huge(id_k) overflows tuple count', &
                  status)
               return
            end if
            call self%cell_data%ensure_tuples(max_gid + 1_id_k, cd_st)
            if (.not. cd_st%ok()) then
               call htg_fail(cd_st%code, trim(cd_st%message), status)
               return
            end if
         end if
      end if

      self%pure_mask_valid = .false.
      call htg_set_ok(status)
   end subroutine grid_set_tree

   subroutine grid_remove_tree(self, tree_id, status)
      class(htg_grid_t), target, intent(inout) :: self
      integer(id_k), intent(in) :: tree_id
      type(htg_status_t), intent(out), optional :: status

      logical :: found
      integer(id_k) :: pos, i

      if (.not. self%dimensions_set) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_grid_t%remove_tree: grid dimensions must be set first', status)
         return
      end if

      if (tree_id < 0_id_k .or. tree_id >= self%get_max_number_of_trees()) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%remove_tree: tree_id out of bounds', status)
         return
      end if

      call lower_bound_tree_slot(self, tree_id, pos, found)
      if (.not. found) then
         call htg_fail(HTG_ERR_NOT_FOUND, &
            'htg_grid_t%remove_tree: tree not present at tree_id', status)
         return
      end if

      if (allocated(self%tree_slots(pos)%tree)) then
         call self%tree_slots(pos)%tree%clear()
         deallocate(self%tree_slots(pos)%tree)
      end if

      do i = pos, self%num_trees - 2_id_k
         self%tree_slots(i)%tree_id = self%tree_slots(i + 1_id_k)%tree_id
         call move_alloc(self%tree_slots(i + 1_id_k)%tree, self%tree_slots(i)%tree)
      end do
      self%tree_slots(self%num_trees - 1_id_k)%tree_id = HTG_INVALID_ID
      self%num_trees = self%num_trees - 1_id_k
      self%pure_mask_valid = .false.

      call htg_set_ok(status)
   end subroutine grid_remove_tree

   subroutine grid_get_tree_by_pos(self, pos, tree_id, tree_ptr, status)
      class(htg_grid_t), target, intent(inout) :: self
      integer(id_k), intent(in) :: pos
      integer(id_k), intent(out) :: tree_id
      type(htg_tree_t), pointer, intent(out) :: tree_ptr
      type(htg_status_t), intent(out), optional :: status

      tree_id = HTG_INVALID_ID
      tree_ptr => null()

      if (pos < 0_id_k .or. pos >= self%num_trees .or. .not. allocated(self%tree_slots)) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%get_tree_by_pos: pos out of range [0, num_non_empty_trees-1]', status)
         return
      end if

      tree_id = self%tree_slots(pos)%tree_id
      tree_ptr => self%tree_slots(pos)%tree
      call htg_set_ok(status)
   end subroutine grid_get_tree_by_pos

   ! --- Forest Aggregations ---

   pure function grid_get_number_of_vertices(self) result(total)
      class(htg_grid_t), intent(in) :: self
      integer(id_k) :: total
      integer(id_k) :: i

      total = 0_id_k
      if (.not. allocated(self%tree_slots)) return
      do i = 0_id_k, self%num_trees - 1_id_k
         if (allocated(self%tree_slots(i)%tree)) then
            total = total + self%tree_slots(i)%tree%get_number_of_vertices()
         end if
      end do
   end function grid_get_number_of_vertices

   pure function grid_get_number_of_leaves(self) result(total)
      class(htg_grid_t), intent(in) :: self
      integer(id_k) :: total
      integer(id_k) :: i

      total = 0_id_k
      if (.not. allocated(self%tree_slots)) return
      do i = 0_id_k, self%num_trees - 1_id_k
         if (allocated(self%tree_slots(i)%tree)) then
            total = total + self%tree_slots(i)%tree%get_number_of_leaves()
         end if
      end do
   end function grid_get_number_of_leaves

   pure function grid_get_number_of_levels(self) result(max_lvl)
      class(htg_grid_t), intent(in) :: self
      integer(id_k) :: max_lvl
      integer(id_k) :: i

      max_lvl = 0_id_k
      if (.not. allocated(self%tree_slots)) return
      do i = 0_id_k, self%num_trees - 1_id_k
         if (allocated(self%tree_slots(i)%tree)) then
            max_lvl = max(max_lvl, self%tree_slots(i)%tree%get_number_of_levels())
         end if
      end do
   end function grid_get_number_of_levels

   function grid_get_max_global_index(self) result(max_gid)
      class(htg_grid_t), intent(in) :: self
      integer(id_k) :: max_gid
      integer(id_k) :: i, t_max

      max_gid = HTG_INVALID_ID
      if (.not. allocated(self%tree_slots)) return
      do i = 0_id_k, self%num_trees - 1_id_k
         if (allocated(self%tree_slots(i)%tree)) then
            t_max = self%tree_slots(i)%tree%get_max_global_index()
            if (t_max > max_gid) max_gid = t_max
         end if
      end do
   end function grid_get_max_global_index

   function grid_get_global_node_index_max(self) result(max_gid)
      class(htg_grid_t), intent(in) :: self
      integer(id_k) :: max_gid

      max_gid = self%get_max_global_index()
      if (max_gid < 0_id_k) max_gid = 0_id_k
   end function grid_get_global_node_index_max

   subroutine grid_initialize_local_index_node(self, status)
      class(htg_grid_t), intent(inout) :: self
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: i, offset
      type(htg_status_t) :: st

      offset = 0_id_k
      if (allocated(self%tree_slots)) then
         do i = 0_id_k, self%num_trees - 1_id_k
            if (allocated(self%tree_slots(i)%tree)) then
               call self%tree_slots(i)%tree%set_global_index_start(offset, st)
               if (.not. st%ok()) then
                  call htg_fail(st%code, trim(st%message), status)
                  return
               end if
               offset = offset + self%tree_slots(i)%tree%get_number_of_vertices()
            end if
         end do
      end if

      if (self%cell_data%get_number_of_arrays() > 0_i32_k .and. offset > 0_id_k) then
         call self%cell_data%ensure_tuples(offset, st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
      end if

      self%pure_mask_valid = .false.
      call htg_set_ok(status)
   end subroutine grid_initialize_local_index_node

   ! --- Masking & PureMask ---

   subroutine grid_set_mask(self, mask_array, status)
      class(htg_grid_t), intent(inout) :: self
      type(htg_bit_array_t), intent(in) :: mask_array
      type(htg_status_t), intent(out), optional :: status

      self%mask = mask_array
      self%mask_active = .true.
      self%pure_mask_valid = .false.
      call htg_set_ok(status)
   end subroutine grid_set_mask

   pure subroutine grid_clear_mask(self)
      class(htg_grid_t), intent(inout) :: self
      call self%mask%clear()
      self%mask_active = .false.
      self%pure_mask_valid = .false.
   end subroutine grid_clear_mask

   pure logical function grid_has_mask(self) result(flag)
      class(htg_grid_t), intent(in) :: self
      flag = self%mask_active .and. (self%mask%get_number_of_values() > 0_id_k)
   end function grid_has_mask

   subroutine grid_get_mask(self, mask_ptr)
      class(htg_grid_t), target, intent(inout) :: self
      type(htg_bit_array_t), pointer, intent(out) :: mask_ptr

      if (self%mask_active) then
         mask_ptr => self%mask
      else
         mask_ptr => null()
      end if
   end subroutine grid_get_mask

   subroutine grid_set_cell_mask(self, global_id, state, status)
      class(htg_grid_t), intent(inout) :: self
      integer(id_k), intent(in) :: global_id
      logical, intent(in) :: state
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st

      if (global_id < 0_id_k) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%set_cell_mask: global_id must be >= 0', status)
         return
      end if

      call self%mask%insert_value(global_id, state, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      self%mask_active = .true.
      self%pure_mask_valid = .false.
      call htg_set_ok(status)
   end subroutine grid_set_cell_mask

   logical function grid_is_cell_masked(self, global_id, status) result(masked)
      class(htg_grid_t), intent(in) :: self
      integer(id_k), intent(in) :: global_id
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st

      masked = .false.
      if (global_id < 0_id_k) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%is_cell_masked: global_id must be >= 0', status)
         return
      end if

      if (self%mask_active .and. global_id < self%mask%get_number_of_values()) then
         masked = self%mask%get_value(global_id, st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
      end if
      call htg_set_ok(status)
   end function grid_is_cell_masked

   pure subroutine grid_invalidate_pure_mask(self)
      class(htg_grid_t), intent(inout) :: self
      self%pure_mask_valid = .false.
   end subroutine grid_invalidate_pure_mask

   pure logical function grid_is_pure_mask_valid(self) result(flag)
      class(htg_grid_t), intent(in) :: self
      flag = self%pure_mask_valid
   end function grid_is_pure_mask_valid

   ! Helper checking whether an unmasked leaf cell is impure due to material
   ! interface intercepts (VTK v9.7.2: Intercepts[gid][2] < 2.0)
   logical function leaf_has_active_interface(intercepts_arr, gid) result(impure)
      type(htg_field_array_t), intent(in) :: intercepts_arr
      integer(id_k), intent(in) :: gid

      real(r64_k) :: val64
      real(r32_k) :: val32
      integer(i64_k) :: vali64
      integer(i32_k) :: vali32
      type(htg_status_t) :: st

      impure = .false.
      if (intercepts_arr%get_number_of_components() /= 3_i32_k) return
      if (gid < 0_id_k .or. gid >= intercepts_arr%get_number_of_tuples()) return

      select case (intercepts_arr%get_data_type())
      case (HTG_TYPE_REAL64)
         val64 = intercepts_arr%get_component_r64(gid, 2_i32_k, st)
         if (st%ok()) impure = (val64 < 2.0_r64_k)
      case (HTG_TYPE_REAL32)
         val32 = intercepts_arr%get_component_r32(gid, 2_i32_k, st)
         if (st%ok()) impure = (val32 < 2.0_r32_k)
      case (HTG_TYPE_INT64)
         vali64 = intercepts_arr%get_component_i64(gid, 2_i32_k, st)
         if (st%ok()) impure = (vali64 < 2_i64_k)
      case (HTG_TYPE_INT32)
         vali32 = intercepts_arr%get_component_i32(gid, 2_i32_k, st)
         if (st%ok()) impure = (vali32 < 2_i32_k)
      end select
   end function leaf_has_active_interface

   subroutine grid_compute_pure_mask_into(self, out_pure_mask, status)
      class(htg_grid_t), intent(in) :: self
      type(htg_bit_array_t), intent(out) :: out_pure_mask
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: max_gid, num_bits, tpos, nv, v, gid, elder, c, child_gid
      integer(i32_k) :: k, ichild, sp, new_cap, idx_arr, normals_idx, intercepts_idx
      integer(i32_k) :: d
      logical :: is_leaf_v, is_masked_v, impure, use_interface
      type(htg_status_t) :: st
      integer(id_k), allocatable :: stk_v(:)
      logical, allocatable :: stk_post(:)
      integer(id_k), allocatable :: tmp_v(:)
      logical, allocatable :: tmp_post(:)
      integer(i32_k), allocatable :: stk_d(:), tmp_d(:)

      call out_pure_mask%clear()

      ! Empty grid: VTK v9.7.2 GetGlobalNodeIndexMax() == 0 -> 1 bit [0] = false
      if (self%num_trees <= 0_id_k) then
         call out_pure_mask%initialize(1_id_k, .false., st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
         call htg_set_ok(status)
         return
      end if

      do tpos = 0_id_k, self%num_trees - 1_id_k
         if (allocated(self%tree_slots(tpos)%tree)) then
            if (self%tree_slots(tpos)%tree%get_indexing_mode() == HTG_INDEX_UNSET) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_grid_t%get_pure_mask: global indices must be assigned first', &
                  status)
               return
            end if
         end if
      end do

      max_gid = self%get_max_global_index()
      if (max_gid < 0_id_k) then
         call htg_fail(HTG_ERR_INDEXING, &
            'htg_grid_t%get_pure_mask: global indices must be assigned first', status)
         return
      end if

      if (max_gid == huge(0_id_k)) then
         call htg_fail(HTG_ERR_INDEXING, &
            'htg_grid_t%get_pure_mask: max_global_index == huge(id_k) overflows bit count', &
            status)
         return
      end if

      num_bits = max_gid + 1_id_k
      call out_pure_mask%initialize(num_bits, .false., st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      ! VTK v9.7.2: enable interface leaf check iff HasInterface is true AND both
      ! InterfaceNormals and InterfaceIntercepts exist in CellData with 3 components.
      normals_idx = -1_i32_k
      intercepts_idx = -1_i32_k
      if (self%has_interface .and. &
          allocated(self%interface_normals_name) .and. &
          allocated(self%interface_intercepts_name)) then
         if (len_trim(self%interface_normals_name) > 0 .and. &
             len_trim(self%interface_intercepts_name) > 0 .and. &
             allocated(self%cell_data%arrays)) then
            do idx_arr = 0_i32_k, self%cell_data%num_arrays - 1_i32_k
               if (allocated(self%cell_data%arrays(idx_arr)%item)) then
                  if (self%cell_data%arrays(idx_arr)%item%get_number_of_components() &
                      == 3_i32_k) then
                     if (self%cell_data%arrays(idx_arr)%item%get_name() == &
                         trim(self%interface_normals_name)) then
                        normals_idx = idx_arr
                     end if
                     if (self%cell_data%arrays(idx_arr)%item%get_name() == &
                         trim(self%interface_intercepts_name)) then
                        intercepts_idx = idx_arr
                     end if
                  end if
               end if
            end do
         end if
      end if
      use_interface = (normals_idx >= 0_i32_k .and. intercepts_idx >= 0_i32_k)

      allocate(stk_v(0:63))
      allocate(stk_post(0:63))
      allocate(stk_d(0:63))

      do tpos = 0_id_k, self%num_trees - 1_id_k
         if (.not. allocated(self%tree_slots(tpos)%tree)) cycle
         associate (tree => self%tree_slots(tpos)%tree)
            nv = tree%get_number_of_vertices()
            if (nv > 0_id_k) then
               k = tree%get_number_of_children()

               sp = 0_i32_k
               stk_v(0) = 0_id_k
               stk_post(0) = .false.
               stk_d(0) = 0_i32_k

               do while (sp >= 0_i32_k)
                  v = stk_v(sp)
                  if (.not. stk_post(sp)) then
                     d = stk_d(sp)
                     sp = sp - 1_i32_k
                     gid = tree%get_global_index_from_local(v, st)
                     if (.not. st%ok()) then
                        call htg_fail(st%code, trim(st%message), status)
                        return
                     end if

                     is_masked_v = self%is_cell_masked(gid, st)
                     if (.not. st%ok()) then
                        call htg_fail(st%code, trim(st%message), status)
                        return
                     end if

                     ! Masked coarse or leaf short-circuits immediately (no ToChild)
                     if (is_masked_v) then
                        call out_pure_mask%set_value(gid, .true., st)
                        cycle
                     end if

                     ! VTK v9.7.2 cursor IsLeaf: a cell at the depth limiter is a leaf,
                     ! so its descendants are never visited and keep PureMask == false
                     if (d >= self%depth_limiter) then
                        is_leaf_v = .true.
                     else
                        is_leaf_v = tree%is_leaf(v, st)
                        if (.not. st%ok()) then
                           call htg_fail(st%code, trim(st%message), status)
                           return
                        end if
                     end if

                     if (is_leaf_v) then
                        if (use_interface) then
                           impure = leaf_has_active_interface( &
                              self%cell_data%arrays(intercepts_idx)%item, gid)
                        else
                           impure = .false.
                        end if
                        call out_pure_mask%set_value(gid, impure, st)
                     else
                        elder = tree%get_elder_child_index(v, st)
                        if (.not. st%ok()) then
                           call htg_fail(st%code, trim(st%message), status)
                           return
                        end if

                        ! Ensure stack capacity for 1 post-order parent + k children
                        if (sp + 1_i32_k + k > ubound(stk_v, 1)) then
                           new_cap = max(2_i32_k * size(stk_v, kind=i32_k), &
                              sp + 2_i32_k + k)
                           allocate(tmp_v(0 : new_cap - 1_i32_k))
                           allocate(tmp_post(0 : new_cap - 1_i32_k))
                           allocate(tmp_d(0 : new_cap - 1_i32_k))
                           if (sp >= 0_i32_k) then
                              tmp_v(0:sp) = stk_v(0:sp)
                              tmp_post(0:sp) = stk_post(0:sp)
                              tmp_d(0:sp) = stk_d(0:sp)
                           end if
                           call move_alloc(tmp_v, stk_v)
                           call move_alloc(tmp_post, stk_post)
                           call move_alloc(tmp_d, stk_d)
                        end if

                        sp = sp + 1_i32_k
                        stk_v(sp) = v
                        stk_post(sp) = .true.
                        stk_d(sp) = d

                        ! Push all k children so every child is visited (|=)
                        do ichild = k - 1_i32_k, 0_i32_k, -1_i32_k
                           sp = sp + 1_i32_k
                           stk_v(sp) = elder + int(ichild, id_k)
                           stk_post(sp) = .false.
                           stk_d(sp) = d + 1_i32_k
                        end do
                     end if
                  else
                     sp = sp - 1_i32_k
                     gid = tree%get_global_index_from_local(v, st)
                     elder = tree%get_elder_child_index(v, st)
                     impure = .false.
                     do ichild = 0_i32_k, k - 1_i32_k
                        c = elder + int(ichild, id_k)
                        child_gid = tree%get_global_index_from_local(c, st)
                        if (out_pure_mask%get_value(child_gid, st)) then
                           impure = .true.
                           exit
                        end if
                     end do
                     call out_pure_mask%set_value(gid, impure, st)
                  end if
               end do
            end if
         end associate
      end do

      call htg_set_ok(status)
   end subroutine grid_compute_pure_mask_into

   subroutine grid_get_pure_mask(self, pure_mask_ptr, status)
      class(htg_grid_t), target, intent(inout) :: self
      type(htg_bit_array_t), pointer, intent(out) :: pure_mask_ptr
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st

      pure_mask_ptr => null()
      if (.not. self%pure_mask_valid) then
         call self%compute_pure_mask_into(self%pure_mask, st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
         self%pure_mask_valid = .true.
      end if

      pure_mask_ptr => self%pure_mask
      call htg_set_ok(status)
   end subroutine grid_get_pure_mask

   ! --- CellData & Material Interface Metadata ---

   subroutine grid_get_cell_data(self, ptr)
      class(htg_grid_t), target, intent(inout) :: self
      type(htg_cell_data_t), pointer, intent(out) :: ptr
      ptr => self%cell_data
   end subroutine grid_get_cell_data

   !> Copy the cell attributes of global id `from_id` onto `to_id`: every
   !> CellData tuple (in arrays long enough to hold both) and the mask bit.
   !> `from_id` keeps its values. Use it to move a cell to another global id.
   subroutine grid_copy_cell_attributes(self, from_id, to_id, status)
      class(htg_grid_t), intent(inout) :: self
      integer(id_k), intent(in) :: from_id, to_id
      type(htg_status_t), intent(out), optional :: status

      logical :: masked
      type(htg_status_t) :: st

      if (from_id < 0_id_k .or. to_id < 0_id_k) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%copy_cell_attributes: global ids must be >= 0', status)
         return
      end if

      call self%cell_data%copy_tuple(to_id, from_id)

      if (self%mask_active) then
         masked = .false.
         if (from_id < self%mask%get_number_of_values()) masked = self%mask%get_value(from_id, st)
         if (to_id < self%mask%get_number_of_values()) then
            call self%mask%set_value(to_id, masked, st)
         else if (masked) then
            call self%mask%insert_value(to_id, .true., st)
         end if
      end if
      self%pure_mask_valid = .false.
      call htg_set_ok(status)
   end subroutine grid_copy_cell_attributes

   !> Reset the cell attributes of a global id that no cell uses any more:
   !> its CellData tuples become zero and its mask bit is cleared, so the id
   !> can be reused for a new cell.
   subroutine grid_clear_cell_attributes(self, global_id, status)
      class(htg_grid_t), intent(inout) :: self
      integer(id_k), intent(in) :: global_id
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st

      if (global_id < 0_id_k) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_grid_t%clear_cell_attributes: global id must be >= 0', status)
         return
      end if

      call self%cell_data%clear_tuple(global_id)
      if (self%mask_active .and. global_id < self%mask%get_number_of_values()) then
         call self%mask%set_value(global_id, .false., st)
      end if
      self%pure_mask_valid = .false.
      call htg_set_ok(status)
   end subroutine grid_clear_cell_attributes

   pure subroutine grid_set_has_interface(self, flag)
      class(htg_grid_t), intent(inout) :: self
      logical, intent(in) :: flag
      if (self%has_interface .neqv. flag) then
         self%has_interface = flag
         self%pure_mask_valid = .false.
      end if
   end subroutine grid_set_has_interface

   pure logical function grid_get_has_interface(self) result(flag)
      class(htg_grid_t), intent(in) :: self
      flag = self%has_interface
   end function grid_get_has_interface

   pure subroutine grid_set_interface_normals_name(self, name)
      class(htg_grid_t), intent(inout) :: self
      character(len=*), intent(in) :: name
      self%interface_normals_name = trim(name)
      self%pure_mask_valid = .false.
   end subroutine grid_set_interface_normals_name

   pure function grid_get_interface_normals_name(self) result(name)
      class(htg_grid_t), intent(in) :: self
      character(len=:), allocatable :: name
      if (allocated(self%interface_normals_name)) then
         name = self%interface_normals_name
      else
         name = ''
      end if
   end function grid_get_interface_normals_name

   pure subroutine grid_set_interface_intercepts_name(self, name)
      class(htg_grid_t), intent(inout) :: self
      character(len=*), intent(in) :: name
      self%interface_intercepts_name = trim(name)
      self%pure_mask_valid = .false.
   end subroutine grid_set_interface_intercepts_name

   pure function grid_get_interface_intercepts_name(self) result(name)
      class(htg_grid_t), intent(in) :: self
      character(len=:), allocatable :: name
      if (allocated(self%interface_intercepts_name)) then
         name = self%interface_intercepts_name
      else
         name = ''
      end if
   end function grid_get_interface_intercepts_name

end module fvtkhtg_grid_m
