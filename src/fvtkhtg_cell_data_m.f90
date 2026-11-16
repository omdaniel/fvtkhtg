module fvtkhtg_cell_data_m
   use fvtkhtg_kinds_m, only: byte_k, i32_k, i64_k, id_k, r32_k, r64_k, &
      htg_status_t, HTG_ERR_INVALID_ARG, HTG_ERR_INVALID_STATE, &
      HTG_ERR_NOT_FOUND, HTG_ERR_OUT_OF_BOUNDS, HTG_TYPE_BIT, &
      HTG_TYPE_INT32, HTG_TYPE_INT64, HTG_TYPE_REAL32, HTG_TYPE_REAL64, &
      htg_fail, htg_set_ok
   use fvtkhtg_bit_array_m, only: htg_bit_array_t
   implicit none
   private

   integer(id_k), parameter :: MIN_TUPLE_GROWTH = 8_id_k
   integer(i32_k), parameter :: MIN_ARRAY_GROWTH = 4_i32_k

   ! Named multi-component 0-based cell field array
   type, public :: htg_field_array_t
      character(len=:), allocatable :: name
      integer(i32_k) :: data_type = HTG_TYPE_REAL64
      integer(i32_k) :: num_components = 1_i32_k
      integer(id_k) :: num_tuples = 0_id_k
      integer(id_k) :: tuple_capacity = 0_id_k
      real(r64_k), allocatable :: r64_values(:, :)
      real(r32_k), allocatable :: r32_values(:, :)
      integer(i64_k), allocatable :: i64_values(:, :)
      integer(i32_k), allocatable :: i32_values(:, :)
      type(htg_bit_array_t) :: bit_values
   contains
      procedure :: initialize => field_initialize
      procedure :: clear => field_clear
      procedure :: is_initialized => field_is_initialized
      procedure :: get_name => field_get_name
      procedure :: get_data_type => field_get_data_type
      procedure :: get_number_of_components => field_get_num_components
      procedure :: get_number_of_tuples => field_get_num_tuples
      procedure :: get_tuple_capacity => field_get_tuple_capacity
      procedure :: ensure_tuples => field_ensure_tuples
      procedure :: set_number_of_tuples => field_set_number_of_tuples
      procedure :: copy_tuple => field_copy_tuple
      procedure :: clear_tuple => field_clear_tuple

      procedure :: set_component_r64 => field_set_component_r64
      procedure :: get_component_r64 => field_get_component_r64
      procedure :: set_component_r32 => field_set_component_r32
      procedure :: get_component_r32 => field_get_component_r32
      procedure :: set_component_i64 => field_set_component_i64
      procedure :: get_component_i64 => field_get_component_i64
      procedure :: set_component_i32 => field_set_component_i32
      procedure :: get_component_i32 => field_get_component_i32
      procedure :: set_component_bit => field_set_component_bit
      procedure :: get_component_bit => field_get_component_bit

      procedure :: set_scalar_r64 => field_set_scalar_r64
      procedure :: get_scalar_r64 => field_get_scalar_r64
      procedure :: set_scalar_r32 => field_set_scalar_r32
      procedure :: get_scalar_r32 => field_get_scalar_r32
      procedure :: set_scalar_i64 => field_set_scalar_i64
      procedure :: get_scalar_i64 => field_get_scalar_i64
      procedure :: set_scalar_i32 => field_set_scalar_i32
      procedure :: get_scalar_i32 => field_get_scalar_i32
      procedure :: set_scalar_bit => field_set_scalar_bit
      procedure :: get_scalar_bit => field_get_scalar_bit

      procedure :: set_tuple_r64 => field_set_tuple_r64
      procedure :: get_tuple_r64 => field_get_tuple_r64
      procedure :: set_tuple_r32 => field_set_tuple_r32
      procedure :: get_tuple_r32 => field_get_tuple_r32
      procedure :: set_tuple_i64 => field_set_tuple_i64
      procedure :: get_tuple_i64 => field_get_tuple_i64
      procedure :: set_tuple_i32 => field_set_tuple_i32
      procedure :: get_tuple_i32 => field_get_tuple_i32
      procedure :: set_tuple_bit => field_set_tuple_bit
      procedure :: get_tuple_bit => field_get_tuple_bit

      procedure :: get_view_r64 => field_get_view_r64
      procedure :: get_view_r32 => field_get_view_r32
      procedure :: get_view_i64 => field_get_view_i64
      procedure :: get_view_i32 => field_get_view_i32
      procedure :: get_view_bit => field_get_view_bit
      procedure :: get_view_bit_bytes => field_get_view_bit_bytes
   end type htg_field_array_t

   ! Boxed field slot preserving heap address across container reallocation
   type, public :: htg_field_box_t
      type(htg_field_array_t), allocatable :: item
   end type htg_field_box_t

   ! Dynamic collection of named cell data arrays indexed by 0-based global ID
   type, public :: htg_cell_data_t
      integer(i32_k) :: num_arrays = 0_i32_k
      integer(i32_k) :: array_capacity = 0_i32_k
      integer(id_k) :: num_tuples = 0_id_k
      type(htg_field_box_t), allocatable :: arrays(:)
      character(len=:), allocatable :: active_scalars
   contains
      procedure :: clear => cell_data_clear
      procedure :: get_number_of_arrays => cell_data_get_num_arrays
      procedure :: get_number_of_tuples => cell_data_get_num_tuples
      procedure :: has_array => cell_data_has_array
      procedure :: add_array => cell_data_add_array
      procedure :: remove_array => cell_data_remove_array
      procedure :: get_array => cell_data_get_array
      procedure :: get_array_by_index => cell_data_get_array_by_index
      procedure :: set_active_scalars => cell_data_set_active_scalars
      procedure :: get_active_scalars => cell_data_get_active_scalars
      procedure :: has_active_scalars => cell_data_has_active_scalars
      procedure :: ensure_tuples => cell_data_ensure_tuples
      procedure :: set_number_of_tuples => cell_data_set_number_of_tuples
      procedure :: copy_tuple => cell_data_copy_tuple
      procedure :: clear_tuple => cell_data_clear_tuple
   end type htg_cell_data_t

contains

   pure logical function is_valid_data_type(dtype) result(valid)
      integer(i32_k), intent(in) :: dtype
      select case (dtype)
      case (HTG_TYPE_REAL64, HTG_TYPE_REAL32, HTG_TYPE_INT64, HTG_TYPE_INT32, HTG_TYPE_BIT)
         valid = .true.
      case default
         valid = .false.
      end select
   end function is_valid_data_type

   pure subroutine field_clear(self)
      class(htg_field_array_t), intent(inout) :: self
      if (allocated(self%name)) deallocate(self%name)
      if (allocated(self%r64_values)) deallocate(self%r64_values)
      if (allocated(self%r32_values)) deallocate(self%r32_values)
      if (allocated(self%i64_values)) deallocate(self%i64_values)
      if (allocated(self%i32_values)) deallocate(self%i32_values)
      call self%bit_values%clear()
      self%data_type = HTG_TYPE_REAL64
      self%num_components = 1_i32_k
      self%num_tuples = 0_id_k
      self%tuple_capacity = 0_id_k
   end subroutine field_clear

   pure logical function field_is_initialized(self) result(flag)
      class(htg_field_array_t), intent(in) :: self
      flag = allocated(self%name)
   end function field_is_initialized

   pure function field_get_name(self) result(res)
      class(htg_field_array_t), intent(in) :: self
      character(len=:), allocatable :: res
      if (allocated(self%name)) then
         res = self%name
      else
         res = ''
      end if
   end function field_get_name

   pure function field_get_data_type(self) result(dtype)
      class(htg_field_array_t), intent(in) :: self
      integer(i32_k) :: dtype
      dtype = self%data_type
   end function field_get_data_type

   pure function field_get_num_components(self) result(nc)
      class(htg_field_array_t), intent(in) :: self
      integer(i32_k) :: nc
      nc = self%num_components
   end function field_get_num_components

   pure function field_get_num_tuples(self) result(nt)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k) :: nt
      nt = self%num_tuples
   end function field_get_num_tuples

   pure function field_get_tuple_capacity(self) result(cap)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k) :: cap
      cap = self%tuple_capacity
   end function field_get_tuple_capacity

   pure subroutine field_initialize(self, name, data_type, num_components, num_tuples, status)
      class(htg_field_array_t), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer(i32_k), intent(in) :: data_type
      integer(i32_k), intent(in), optional :: num_components
      integer(id_k), intent(in), optional :: num_tuples
      type(htg_status_t), intent(out), optional :: status

      integer(i32_k) :: nc
      integer(id_k) :: nt
      type(htg_status_t) :: bit_st

      nc = 1_i32_k
      if (present(num_components)) nc = num_components

      nt = 0_id_k
      if (present(num_tuples)) nt = num_tuples

      if (len_trim(name) == 0) then
         call self%clear()
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_field_array_t%initialize: array name must be non-empty', status)
         return
      end if

      if (.not. is_valid_data_type(data_type)) then
         call self%clear()
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_field_array_t%initialize: unrecognized data_type', status)
         return
      end if

      if (nc < 1_i32_k) then
         call self%clear()
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_field_array_t%initialize: num_components must be >= 1', status)
         return
      end if

      if (nt < 0_id_k) then
         call self%clear()
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_field_array_t%initialize: num_tuples must be >= 0', status)
         return
      end if

      call self%clear()
      self%name = trim(name)
      self%data_type = data_type
      self%num_components = nc
      self%num_tuples = nt
      self%tuple_capacity = nt

      if (nt > 0_id_k) then
         select case (data_type)
         case (HTG_TYPE_REAL64)
            allocate(self%r64_values(0_i32_k : nc - 1_i32_k, 0_id_k : nt - 1_id_k))
            self%r64_values = 0.0_r64_k
         case (HTG_TYPE_REAL32)
            allocate(self%r32_values(0_i32_k : nc - 1_i32_k, 0_id_k : nt - 1_id_k))
            self%r32_values = 0.0_r32_k
         case (HTG_TYPE_INT64)
            allocate(self%i64_values(0_i32_k : nc - 1_i32_k, 0_id_k : nt - 1_id_k))
            self%i64_values = 0_i64_k
         case (HTG_TYPE_INT32)
            allocate(self%i32_values(0_i32_k : nc - 1_i32_k, 0_id_k : nt - 1_id_k))
            self%i32_values = 0_i32_k
         case (HTG_TYPE_BIT)
            call self%bit_values%initialize(int(nc, id_k) * nt, .false., bit_st)
            if (.not. bit_st%ok()) then
               call self%clear()
               call htg_fail(bit_st%code, trim(bit_st%message), status)
               return
            end if
            self%tuple_capacity = (self%bit_values%byte_capacity * 8_id_k) / int(nc, id_k)
         end select
      end if

      call htg_set_ok(status)
   end subroutine field_initialize

   pure subroutine field_ensure_tuples(self, min_tuples, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: min_tuples
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: new_cap, old_nt
      integer(i32_k) :: nc
      real(r64_k), allocatable :: tmp_r64(:, :)
      real(r32_k), allocatable :: tmp_r32(:, :)
      integer(i64_k), allocatable :: tmp_i64(:, :)
      integer(i32_k), allocatable :: tmp_i32(:, :)
      type(htg_status_t) :: bit_st

      if (.not. allocated(self%name)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%ensure_tuples: field array is not initialized', status)
         return
      end if

      if (min_tuples < 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_field_array_t%ensure_tuples: min_tuples must be >= 0', status)
         return
      end if

      if (min_tuples <= self%num_tuples) then
         call htg_set_ok(status)
         return
      end if

      nc = self%num_components
      old_nt = self%num_tuples

      if (self%data_type == HTG_TYPE_BIT) then
         call self%bit_values%resize(int(nc, id_k) * min_tuples, bit_st)
         if (.not. bit_st%ok()) then
            call htg_fail(bit_st%code, trim(bit_st%message), status)
            return
         end if
         self%num_tuples = min_tuples
         self%tuple_capacity = max(self%tuple_capacity, &
            (self%bit_values%byte_capacity * 8_id_k) / int(nc, id_k))
         call htg_set_ok(status)
         return
      end if

      if (min_tuples > self%tuple_capacity) then
         new_cap = max(min_tuples, max(MIN_TUPLE_GROWTH, 2_id_k * self%tuple_capacity))
         select case (self%data_type)
         case (HTG_TYPE_REAL64)
            allocate(tmp_r64(0_i32_k : nc - 1_i32_k, 0_id_k : new_cap - 1_id_k))
            tmp_r64 = 0.0_r64_k
            if (old_nt > 0_id_k .and. allocated(self%r64_values)) then
               tmp_r64(:, 0_id_k : old_nt - 1_id_k) = self%r64_values(:, 0_id_k : old_nt - 1_id_k)
            end if
            call move_alloc(tmp_r64, self%r64_values)
         case (HTG_TYPE_REAL32)
            allocate(tmp_r32(0_i32_k : nc - 1_i32_k, 0_id_k : new_cap - 1_id_k))
            tmp_r32 = 0.0_r32_k
            if (old_nt > 0_id_k .and. allocated(self%r32_values)) then
               tmp_r32(:, 0_id_k : old_nt - 1_id_k) = self%r32_values(:, 0_id_k : old_nt - 1_id_k)
            end if
            call move_alloc(tmp_r32, self%r32_values)
         case (HTG_TYPE_INT64)
            allocate(tmp_i64(0_i32_k : nc - 1_i32_k, 0_id_k : new_cap - 1_id_k))
            tmp_i64 = 0_i64_k
            if (old_nt > 0_id_k .and. allocated(self%i64_values)) then
               tmp_i64(:, 0_id_k : old_nt - 1_id_k) = self%i64_values(:, 0_id_k : old_nt - 1_id_k)
            end if
            call move_alloc(tmp_i64, self%i64_values)
         case (HTG_TYPE_INT32)
            allocate(tmp_i32(0_i32_k : nc - 1_i32_k, 0_id_k : new_cap - 1_id_k))
            tmp_i32 = 0_i32_k
            if (old_nt > 0_id_k .and. allocated(self%i32_values)) then
               tmp_i32(:, 0_id_k : old_nt - 1_id_k) = self%i32_values(:, 0_id_k : old_nt - 1_id_k)
            end if
            call move_alloc(tmp_i32, self%i32_values)
         end select
         self%tuple_capacity = new_cap
      else
         select case (self%data_type)
         case (HTG_TYPE_REAL64)
            self%r64_values(:, old_nt : min_tuples - 1_id_k) = 0.0_r64_k
         case (HTG_TYPE_REAL32)
            self%r32_values(:, old_nt : min_tuples - 1_id_k) = 0.0_r32_k
         case (HTG_TYPE_INT64)
            self%i64_values(:, old_nt : min_tuples - 1_id_k) = 0_i64_k
         case (HTG_TYPE_INT32)
            self%i32_values(:, old_nt : min_tuples - 1_id_k) = 0_i32_k
         end select
      end if

      self%num_tuples = min_tuples
      call htg_set_ok(status)
   end subroutine field_ensure_tuples

   pure subroutine field_set_number_of_tuples(self, num_tuples, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: num_tuples
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: bit_st

      if (.not. allocated(self%name)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%set_number_of_tuples: field array is not initialized', status)
         return
      end if

      if (num_tuples < 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_field_array_t%set_number_of_tuples: num_tuples must be >= 0', status)
         return
      end if

      if (num_tuples > self%num_tuples) then
         call self%ensure_tuples(num_tuples, status)
         return
      else if (num_tuples < self%num_tuples) then
         select case (self%data_type)
         case (HTG_TYPE_REAL64)
            if (allocated(self%r64_values)) then
               self%r64_values(:, num_tuples : self%num_tuples - 1_id_k) = 0.0_r64_k
            end if
         case (HTG_TYPE_REAL32)
            if (allocated(self%r32_values)) then
               self%r32_values(:, num_tuples : self%num_tuples - 1_id_k) = 0.0_r32_k
            end if
         case (HTG_TYPE_INT64)
            if (allocated(self%i64_values)) then
               self%i64_values(:, num_tuples : self%num_tuples - 1_id_k) = 0_i64_k
            end if
         case (HTG_TYPE_INT32)
            if (allocated(self%i32_values)) then
               self%i32_values(:, num_tuples : self%num_tuples - 1_id_k) = 0_i32_k
            end if
         case (HTG_TYPE_BIT)
            call self%bit_values%resize(int(self%num_components, id_k) * num_tuples, bit_st)
            if (.not. bit_st%ok()) then
               call htg_fail(bit_st%code, trim(bit_st%message), status)
               return
            end if
         end select
         self%num_tuples = num_tuples
      end if

      call htg_set_ok(status)
   end subroutine field_set_number_of_tuples

   !> Copy every component of tuple src_id onto tuple dst_id.
   pure subroutine field_copy_tuple(self, dst_id, src_id, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: dst_id, src_id
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: nc, c
      logical :: bit
      type(htg_status_t) :: st

      if (.not. allocated(self%name)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%copy_tuple: field array is not initialized', status)
         return
      end if
      if (dst_id < 0_id_k .or. dst_id >= self%num_tuples .or. &
          src_id < 0_id_k .or. src_id >= self%num_tuples) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_field_array_t%copy_tuple: tuple id out of bounds', status)
         return
      end if

      select case (self%data_type)
      case (HTG_TYPE_REAL64)
         self%r64_values(:, dst_id) = self%r64_values(:, src_id)
      case (HTG_TYPE_REAL32)
         self%r32_values(:, dst_id) = self%r32_values(:, src_id)
      case (HTG_TYPE_INT64)
         self%i64_values(:, dst_id) = self%i64_values(:, src_id)
      case (HTG_TYPE_INT32)
         self%i32_values(:, dst_id) = self%i32_values(:, src_id)
      case (HTG_TYPE_BIT)
         nc = int(self%num_components, id_k)
         do c = 0_id_k, nc - 1_id_k
            bit = btest(self%bit_values%bytes((src_id * nc + c) / 8_id_k), &
               7 - int(mod(src_id * nc + c, 8_id_k)))
            call self%bit_values%set_value(dst_id * nc + c, bit, st)
         end do
      end select
      call htg_set_ok(status)
   end subroutine field_copy_tuple

   !> Reset every component of a tuple to zero (.false. for BIT arrays).
   pure subroutine field_clear_tuple(self, tuple_id, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: nc, c
      type(htg_status_t) :: st

      if (.not. allocated(self%name)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%clear_tuple: field array is not initialized', status)
         return
      end if
      if (tuple_id < 0_id_k .or. tuple_id >= self%num_tuples) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_field_array_t%clear_tuple: tuple id out of bounds', status)
         return
      end if

      select case (self%data_type)
      case (HTG_TYPE_REAL64)
         self%r64_values(:, tuple_id) = 0.0_r64_k
      case (HTG_TYPE_REAL32)
         self%r32_values(:, tuple_id) = 0.0_r32_k
      case (HTG_TYPE_INT64)
         self%i64_values(:, tuple_id) = 0_i64_k
      case (HTG_TYPE_INT32)
         self%i32_values(:, tuple_id) = 0_i32_k
      case (HTG_TYPE_BIT)
         nc = int(self%num_components, id_k)
         do c = 0_id_k, nc - 1_id_k
            call self%bit_values%set_value(tuple_id * nc + c, .false., st)
         end do
      end select
      call htg_set_ok(status)
   end subroutine field_clear_tuple

   pure subroutine check_component_access(self, expected_type, tuple_id, comp_id, ok, status)
      class(htg_field_array_t), intent(in) :: self
      integer(i32_k), intent(in) :: expected_type
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      logical, intent(out) :: ok
      type(htg_status_t), intent(out), optional :: status

      ok = .false.
      if (.not. allocated(self%name)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t: field array is not initialized', status)
         return
      end if
      if (self%data_type /= expected_type) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t: mismatched field data_type', status)
         return
      end if
      if (tuple_id < 0_id_k .or. tuple_id >= self%num_tuples) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_field_array_t: tuple_id out of bounds', status)
         return
      end if
      if (comp_id < 0_i32_k .or. comp_id >= self%num_components) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_field_array_t: comp_id out of bounds', status)
         return
      end if
      ok = .true.
   end subroutine check_component_access

   pure subroutine check_tuple_access(self, expected_type, tuple_id, nvals, ok, status)
      class(htg_field_array_t), intent(in) :: self
      integer(i32_k), intent(in) :: expected_type
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: nvals
      logical, intent(out) :: ok
      type(htg_status_t), intent(out), optional :: status

      ok = .false.
      if (.not. allocated(self%name)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t: field array is not initialized', status)
         return
      end if
      if (self%data_type /= expected_type) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t: mismatched field data_type', status)
         return
      end if
      if (tuple_id < 0_id_k .or. tuple_id >= self%num_tuples) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_field_array_t: tuple_id out of bounds', status)
         return
      end if
      if (nvals /= self%num_components) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_field_array_t: tuple slice length does not match num_components', status)
         return
      end if
      ok = .true.
   end subroutine check_tuple_access

   ! --- Component getters and setters ---

   pure subroutine field_set_component_r64(self, tuple_id, comp_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      real(r64_k), intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      call check_component_access(self, HTG_TYPE_REAL64, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      self%r64_values(comp_id, tuple_id) = val
      call htg_set_ok(status)
   end subroutine field_set_component_r64

   real(r64_k) function field_get_component_r64(self, tuple_id, comp_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      val = 0.0_r64_k
      call check_component_access(self, HTG_TYPE_REAL64, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      val = self%r64_values(comp_id, tuple_id)
      call htg_set_ok(status)
   end function field_get_component_r64

   pure subroutine field_set_component_r32(self, tuple_id, comp_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      real(r32_k), intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      call check_component_access(self, HTG_TYPE_REAL32, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      self%r32_values(comp_id, tuple_id) = val
      call htg_set_ok(status)
   end subroutine field_set_component_r32

   real(r32_k) function field_get_component_r32(self, tuple_id, comp_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      val = 0.0_r32_k
      call check_component_access(self, HTG_TYPE_REAL32, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      val = self%r32_values(comp_id, tuple_id)
      call htg_set_ok(status)
   end function field_get_component_r32

   pure subroutine field_set_component_i64(self, tuple_id, comp_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      integer(i64_k), intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      call check_component_access(self, HTG_TYPE_INT64, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      self%i64_values(comp_id, tuple_id) = val
      call htg_set_ok(status)
   end subroutine field_set_component_i64

   integer(i64_k) function field_get_component_i64(self, tuple_id, comp_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      val = 0_i64_k
      call check_component_access(self, HTG_TYPE_INT64, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      val = self%i64_values(comp_id, tuple_id)
      call htg_set_ok(status)
   end function field_get_component_i64

   pure subroutine field_set_component_i32(self, tuple_id, comp_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      integer(i32_k), intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      call check_component_access(self, HTG_TYPE_INT32, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      self%i32_values(comp_id, tuple_id) = val
      call htg_set_ok(status)
   end subroutine field_set_component_i32

   integer(i32_k) function field_get_component_i32(self, tuple_id, comp_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      val = 0_i32_k
      call check_component_access(self, HTG_TYPE_INT32, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      val = self%i32_values(comp_id, tuple_id)
      call htg_set_ok(status)
   end function field_get_component_i32

   pure subroutine field_set_component_bit(self, tuple_id, comp_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      logical, intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      logical :: ok
      integer(id_k) :: bit_id

      call check_component_access(self, HTG_TYPE_BIT, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      bit_id = tuple_id * int(self%num_components, id_k) + int(comp_id, id_k)
      call self%bit_values%set_value(bit_id, val, status)
   end subroutine field_set_component_bit

   logical function field_get_component_bit(self, tuple_id, comp_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: comp_id
      type(htg_status_t), intent(out), optional :: status
      logical :: ok
      integer(id_k) :: bit_id

      val = .false.
      call check_component_access(self, HTG_TYPE_BIT, tuple_id, comp_id, ok, status)
      if (.not. ok) return
      bit_id = tuple_id * int(self%num_components, id_k) + int(comp_id, id_k)
      val = self%bit_values%get_value(bit_id, status)
   end function field_get_component_bit

   ! --- Scalar convenience getters and setters (component 0) ---

   pure subroutine field_set_scalar_r64(self, tuple_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      real(r64_k), intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      call self%set_component_r64(tuple_id, 0_i32_k, val, status)
   end subroutine field_set_scalar_r64

   real(r64_k) function field_get_scalar_r64(self, tuple_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      type(htg_status_t), intent(out), optional :: status
      val = self%get_component_r64(tuple_id, 0_i32_k, status)
   end function field_get_scalar_r64

   pure subroutine field_set_scalar_r32(self, tuple_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      real(r32_k), intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      call self%set_component_r32(tuple_id, 0_i32_k, val, status)
   end subroutine field_set_scalar_r32

   real(r32_k) function field_get_scalar_r32(self, tuple_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      type(htg_status_t), intent(out), optional :: status
      val = self%get_component_r32(tuple_id, 0_i32_k, status)
   end function field_get_scalar_r32

   pure subroutine field_set_scalar_i64(self, tuple_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i64_k), intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      call self%set_component_i64(tuple_id, 0_i32_k, val, status)
   end subroutine field_set_scalar_i64

   integer(i64_k) function field_get_scalar_i64(self, tuple_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      type(htg_status_t), intent(out), optional :: status
      val = self%get_component_i64(tuple_id, 0_i32_k, status)
   end function field_get_scalar_i64

   pure subroutine field_set_scalar_i32(self, tuple_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      call self%set_component_i32(tuple_id, 0_i32_k, val, status)
   end subroutine field_set_scalar_i32

   integer(i32_k) function field_get_scalar_i32(self, tuple_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      type(htg_status_t), intent(out), optional :: status
      val = self%get_component_i32(tuple_id, 0_i32_k, status)
   end function field_get_scalar_i32

   pure subroutine field_set_scalar_bit(self, tuple_id, val, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      logical, intent(in) :: val
      type(htg_status_t), intent(out), optional :: status
      call self%set_component_bit(tuple_id, 0_i32_k, val, status)
   end subroutine field_set_scalar_bit

   logical function field_get_scalar_bit(self, tuple_id, status) result(val)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      type(htg_status_t), intent(out), optional :: status
      val = self%get_component_bit(tuple_id, 0_i32_k, status)
   end function field_get_scalar_bit

   ! --- Tuple getters and setters (explicit 0: lower bound) ---

   pure subroutine field_set_tuple_r64(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      real(r64_k), intent(in) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      call check_tuple_access(self, HTG_TYPE_REAL64, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      self%r64_values(0_i32_k : self%num_components - 1_i32_k, tuple_id) = &
         vals(0_i32_k : self%num_components - 1_i32_k)
      call htg_set_ok(status)
   end subroutine field_set_tuple_r64

   pure subroutine field_get_tuple_r64(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      real(r64_k), intent(out) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      vals = 0.0_r64_k
      call check_tuple_access(self, HTG_TYPE_REAL64, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      vals(0_i32_k : self%num_components - 1_i32_k) = &
         self%r64_values(0_i32_k : self%num_components - 1_i32_k, tuple_id)
      call htg_set_ok(status)
   end subroutine field_get_tuple_r64

   pure subroutine field_set_tuple_r32(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      real(r32_k), intent(in) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      call check_tuple_access(self, HTG_TYPE_REAL32, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      self%r32_values(0_i32_k : self%num_components - 1_i32_k, tuple_id) = &
         vals(0_i32_k : self%num_components - 1_i32_k)
      call htg_set_ok(status)
   end subroutine field_set_tuple_r32

   pure subroutine field_get_tuple_r32(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      real(r32_k), intent(out) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      vals = 0.0_r32_k
      call check_tuple_access(self, HTG_TYPE_REAL32, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      vals(0_i32_k : self%num_components - 1_i32_k) = &
         self%r32_values(0_i32_k : self%num_components - 1_i32_k, tuple_id)
      call htg_set_ok(status)
   end subroutine field_get_tuple_r32

   pure subroutine field_set_tuple_i64(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i64_k), intent(in) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      call check_tuple_access(self, HTG_TYPE_INT64, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      self%i64_values(0_i32_k : self%num_components - 1_i32_k, tuple_id) = &
         vals(0_i32_k : self%num_components - 1_i32_k)
      call htg_set_ok(status)
   end subroutine field_set_tuple_i64

   pure subroutine field_get_tuple_i64(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i64_k), intent(out) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      vals = 0_i64_k
      call check_tuple_access(self, HTG_TYPE_INT64, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      vals(0_i32_k : self%num_components - 1_i32_k) = &
         self%i64_values(0_i32_k : self%num_components - 1_i32_k, tuple_id)
      call htg_set_ok(status)
   end subroutine field_get_tuple_i64

   pure subroutine field_set_tuple_i32(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(in) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      call check_tuple_access(self, HTG_TYPE_INT32, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      self%i32_values(0_i32_k : self%num_components - 1_i32_k, tuple_id) = &
         vals(0_i32_k : self%num_components - 1_i32_k)
      call htg_set_ok(status)
   end subroutine field_set_tuple_i32

   pure subroutine field_get_tuple_i32(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      integer(i32_k), intent(out) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok

      vals = 0_i32_k
      call check_tuple_access(self, HTG_TYPE_INT32, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      vals(0_i32_k : self%num_components - 1_i32_k) = &
         self%i32_values(0_i32_k : self%num_components - 1_i32_k, tuple_id)
      call htg_set_ok(status)
   end subroutine field_get_tuple_i32

   pure subroutine field_set_tuple_bit(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id
      logical, intent(in) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok
      integer(i32_k) :: c
      integer(id_k) :: base_bit
      type(htg_status_t) :: bit_st

      call check_tuple_access(self, HTG_TYPE_BIT, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      base_bit = tuple_id * int(self%num_components, id_k)
      do c = 0_i32_k, self%num_components - 1_i32_k
         call self%bit_values%set_value(base_bit + int(c, id_k), vals(c), bit_st)
         if (.not. bit_st%ok()) then
            call htg_fail(bit_st%code, trim(bit_st%message), status)
            return
         end if
      end do
      call htg_set_ok(status)
   end subroutine field_set_tuple_bit

   subroutine field_get_tuple_bit(self, tuple_id, vals, status)
      class(htg_field_array_t), intent(in) :: self
      integer(id_k), intent(in) :: tuple_id
      logical, intent(out) :: vals(0:)
      type(htg_status_t), intent(out), optional :: status
      logical :: ok
      integer(i32_k) :: c
      integer(id_k) :: base_bit
      type(htg_status_t) :: bit_st

      vals = .false.
      call check_tuple_access(self, HTG_TYPE_BIT, tuple_id, size(vals, kind=i32_k), ok, status)
      if (.not. ok) return
      base_bit = tuple_id * int(self%num_components, id_k)
      do c = 0_i32_k, self%num_components - 1_i32_k
         vals(c) = self%bit_values%get_value(base_bit + int(c, id_k), bit_st)
         if (.not. bit_st%ok()) then
            call htg_fail(bit_st%code, trim(bit_st%message), status)
            return
         end if
      end do
      call htg_set_ok(status)
   end subroutine field_get_tuple_bit

   ! --- Zero-copy 0-based pointer views ---

   subroutine field_get_view_r64(self, ptr, status)
      class(htg_field_array_t), target, intent(inout) :: self
      real(r64_k), pointer, intent(out) :: ptr(:, :)
      type(htg_status_t), intent(out), optional :: status

      ptr => null()
      if (.not. allocated(self%name) .or. self%data_type /= HTG_TYPE_REAL64) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%get_view_r64: array is not initialized as REAL64', status)
         return
      end if
      if (self%num_tuples > 0_id_k .and. allocated(self%r64_values)) then
         ptr(0_i32_k:, 0_id_k:) => self%r64_values(:, 0_id_k : self%num_tuples - 1_id_k)
      end if
      call htg_set_ok(status)
   end subroutine field_get_view_r64

   subroutine field_get_view_r32(self, ptr, status)
      class(htg_field_array_t), target, intent(inout) :: self
      real(r32_k), pointer, intent(out) :: ptr(:, :)
      type(htg_status_t), intent(out), optional :: status

      ptr => null()
      if (.not. allocated(self%name) .or. self%data_type /= HTG_TYPE_REAL32) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%get_view_r32: array is not initialized as REAL32', status)
         return
      end if
      if (self%num_tuples > 0_id_k .and. allocated(self%r32_values)) then
         ptr(0_i32_k:, 0_id_k:) => self%r32_values(:, 0_id_k : self%num_tuples - 1_id_k)
      end if
      call htg_set_ok(status)
   end subroutine field_get_view_r32

   subroutine field_get_view_i64(self, ptr, status)
      class(htg_field_array_t), target, intent(inout) :: self
      integer(i64_k), pointer, intent(out) :: ptr(:, :)
      type(htg_status_t), intent(out), optional :: status

      ptr => null()
      if (.not. allocated(self%name) .or. self%data_type /= HTG_TYPE_INT64) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%get_view_i64: array is not initialized as INT64', status)
         return
      end if
      if (self%num_tuples > 0_id_k .and. allocated(self%i64_values)) then
         ptr(0_i32_k:, 0_id_k:) => self%i64_values(:, 0_id_k : self%num_tuples - 1_id_k)
      end if
      call htg_set_ok(status)
   end subroutine field_get_view_i64

   subroutine field_get_view_i32(self, ptr, status)
      class(htg_field_array_t), target, intent(inout) :: self
      integer(i32_k), pointer, intent(out) :: ptr(:, :)
      type(htg_status_t), intent(out), optional :: status

      ptr => null()
      if (.not. allocated(self%name) .or. self%data_type /= HTG_TYPE_INT32) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%get_view_i32: array is not initialized as INT32', status)
         return
      end if
      if (self%num_tuples > 0_id_k .and. allocated(self%i32_values)) then
         ptr(0_i32_k:, 0_id_k:) => self%i32_values(:, 0_id_k : self%num_tuples - 1_id_k)
      end if
      call htg_set_ok(status)
   end subroutine field_get_view_i32

   subroutine field_get_view_bit(self, ptr, status)
      class(htg_field_array_t), target, intent(inout) :: self
      type(htg_bit_array_t), pointer, intent(out) :: ptr
      type(htg_status_t), intent(out), optional :: status

      ptr => null()
      if (.not. allocated(self%name) .or. self%data_type /= HTG_TYPE_BIT) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%get_view_bit: array is not initialized as BIT', status)
         return
      end if
      ptr => self%bit_values
      call htg_set_ok(status)
   end subroutine field_get_view_bit

   subroutine field_get_view_bit_bytes(self, ptr, status)
      class(htg_field_array_t), target, intent(inout) :: self
      integer(byte_k), pointer, intent(out) :: ptr(:)
      type(htg_status_t), intent(out), optional :: status

      ptr => null()
      if (.not. allocated(self%name) .or. self%data_type /= HTG_TYPE_BIT) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_field_array_t%get_view_bit_bytes: array is not initialized as BIT', status)
         return
      end if
      call self%bit_values%get_raw_bytes(ptr)
      call htg_set_ok(status)
   end subroutine field_get_view_bit_bytes

   ! --- htg_cell_data_t methods ---

   pure function find_array_index(self, name) result(idx)
      class(htg_cell_data_t), intent(in) :: self
      character(len=*), intent(in) :: name
      integer(i32_k) :: idx
      integer(i32_k) :: i

      idx = -1_i32_k
      if (len_trim(name) == 0 .or. .not. allocated(self%arrays)) return
      do i = 0_i32_k, self%num_arrays - 1_i32_k
         if (allocated(self%arrays(i)%item)) then
            if (allocated(self%arrays(i)%item%name)) then
               if (self%arrays(i)%item%name == trim(name)) then
                  idx = i
                  return
               end if
            end if
         end if
      end do
   end function find_array_index

   pure subroutine cell_data_clear(self)
      class(htg_cell_data_t), intent(inout) :: self
      integer(i32_k) :: i

      if (allocated(self%arrays)) then
         do i = 0_i32_k, self%num_arrays - 1_i32_k
            if (allocated(self%arrays(i)%item)) then
               call self%arrays(i)%item%clear()
               deallocate(self%arrays(i)%item)
            end if
         end do
         deallocate(self%arrays)
      end if
      if (allocated(self%active_scalars)) deallocate(self%active_scalars)
      self%num_arrays = 0_i32_k
      self%array_capacity = 0_i32_k
      self%num_tuples = 0_id_k
   end subroutine cell_data_clear

   pure function cell_data_get_num_arrays(self) result(n)
      class(htg_cell_data_t), intent(in) :: self
      integer(i32_k) :: n
      n = self%num_arrays
   end function cell_data_get_num_arrays

   pure function cell_data_get_num_tuples(self) result(nt)
      class(htg_cell_data_t), intent(in) :: self
      integer(id_k) :: nt
      nt = self%num_tuples
   end function cell_data_get_num_tuples

   pure logical function cell_data_has_array(self, name) result(found)
      class(htg_cell_data_t), intent(in) :: self
      character(len=*), intent(in) :: name
      found = (find_array_index(self, name) >= 0_i32_k)
   end function cell_data_has_array

   subroutine cell_data_add_array(self, name, data_type, num_components, num_tuples, status)
      class(htg_cell_data_t), target, intent(inout) :: self
      character(len=*), intent(in) :: name
      integer(i32_k), intent(in) :: data_type
      integer(i32_k), intent(in), optional :: num_components
      integer(id_k), intent(in), optional :: num_tuples
      type(htg_status_t), intent(out), optional :: status

      integer(i32_k) :: nc, new_cap, i
      integer(id_k) :: nt
      type(htg_field_array_t), allocatable, target :: candidate
      type(htg_field_box_t), allocatable, target :: tmp(:)
      type(htg_status_t) :: init_st

      if (len_trim(name) == 0) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_cell_data_t%add_array: array name must be non-empty', status)
         return
      end if

      if (self%has_array(name)) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_cell_data_t%add_array: duplicate array name', status)
         return
      end if

      nc = 1_i32_k
      if (present(num_components)) nc = num_components

      if (present(num_tuples)) then
         nt = num_tuples
      else
         nt = self%num_tuples
      end if

      allocate(candidate)
      call candidate%initialize(name, data_type, nc, nt, init_st)
      if (.not. init_st%ok()) then
         deallocate(candidate)
         call htg_fail(init_st%code, trim(init_st%message), status)
         return
      end if

      if (self%num_arrays >= self%array_capacity) then
         new_cap = max(self%num_arrays + 1_i32_k, &
            max(MIN_ARRAY_GROWTH, 2_i32_k * self%array_capacity))
         allocate(tmp(0_i32_k : new_cap - 1_i32_k))
         if (allocated(self%arrays)) then
            do i = 0_i32_k, self%num_arrays - 1_i32_k
               call move_alloc(self%arrays(i)%item, tmp(i)%item)
            end do
         end if
         call move_alloc(tmp, self%arrays)
         self%array_capacity = new_cap
      end if

      call move_alloc(candidate, self%arrays(self%num_arrays)%item)
      self%num_arrays = self%num_arrays + 1_i32_k
      self%num_tuples = max(self%num_tuples, nt)
      call htg_set_ok(status)
   end subroutine cell_data_add_array

   subroutine cell_data_remove_array(self, name, status)
      class(htg_cell_data_t), target, intent(inout) :: self
      character(len=*), intent(in) :: name
      type(htg_status_t), intent(out), optional :: status

      integer(i32_k) :: idx, i

      idx = find_array_index(self, name)
      if (idx < 0_i32_k) then
         call htg_fail(HTG_ERR_NOT_FOUND, &
            'htg_cell_data_t%remove_array: array not found', status)
         return
      end if

      call self%arrays(idx)%item%clear()
      deallocate(self%arrays(idx)%item)
      do i = idx, self%num_arrays - 2_i32_k
         call move_alloc(self%arrays(i + 1_i32_k)%item, self%arrays(i)%item)
      end do
      self%num_arrays = self%num_arrays - 1_i32_k

      if (allocated(self%active_scalars)) then
         if (self%active_scalars == trim(name)) deallocate(self%active_scalars)
      end if

      self%num_tuples = 0_id_k
      do i = 0_i32_k, self%num_arrays - 1_i32_k
         if (allocated(self%arrays(i)%item)) then
            self%num_tuples = max(self%num_tuples, self%arrays(i)%item%num_tuples)
         end if
      end do
      call htg_set_ok(status)
   end subroutine cell_data_remove_array

   subroutine cell_data_get_array(self, name, ptr, status)
      class(htg_cell_data_t), target, intent(inout) :: self
      character(len=*), intent(in) :: name
      type(htg_field_array_t), pointer, intent(out) :: ptr
      type(htg_status_t), intent(out), optional :: status

      integer(i32_k) :: idx

      ptr => null()
      idx = find_array_index(self, name)
      if (idx < 0_i32_k) then
         call htg_fail(HTG_ERR_NOT_FOUND, &
            'htg_cell_data_t%get_array: array not found', status)
         return
      end if

      ptr => self%arrays(idx)%item
      call htg_set_ok(status)
   end subroutine cell_data_get_array

   subroutine cell_data_get_array_by_index(self, idx, ptr, status)
      class(htg_cell_data_t), target, intent(inout) :: self
      integer(i32_k), intent(in) :: idx
      type(htg_field_array_t), pointer, intent(out) :: ptr
      type(htg_status_t), intent(out), optional :: status

      ptr => null()
      if (idx < 0_i32_k .or. idx >= self%num_arrays) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_cell_data_t%get_array_by_index: index out of bounds', status)
         return
      end if

      ptr => self%arrays(idx)%item
      call htg_set_ok(status)
   end subroutine cell_data_get_array_by_index

   pure subroutine cell_data_set_active_scalars(self, name, status)
      class(htg_cell_data_t), intent(inout) :: self
      character(len=*), intent(in) :: name
      type(htg_status_t), intent(out), optional :: status

      if (len_trim(name) == 0) then
         if (allocated(self%active_scalars)) deallocate(self%active_scalars)
         call htg_set_ok(status)
         return
      end if

      if (.not. self%has_array(name)) then
         call htg_fail(HTG_ERR_NOT_FOUND, &
            'htg_cell_data_t%set_active_scalars: array not found', status)
         return
      end if

      self%active_scalars = trim(name)
      call htg_set_ok(status)
   end subroutine cell_data_set_active_scalars

   pure function cell_data_get_active_scalars(self) result(res)
      class(htg_cell_data_t), intent(in) :: self
      character(len=:), allocatable :: res
      if (allocated(self%active_scalars)) then
         res = self%active_scalars
      else
         res = ''
      end if
   end function cell_data_get_active_scalars

   pure logical function cell_data_has_active_scalars(self) result(flag)
      class(htg_cell_data_t), intent(in) :: self
      flag = allocated(self%active_scalars)
   end function cell_data_has_active_scalars

   pure subroutine cell_data_ensure_tuples(self, min_tuples, status)
      class(htg_cell_data_t), intent(inout) :: self
      integer(id_k), intent(in) :: min_tuples
      type(htg_status_t), intent(out), optional :: status

      integer(i32_k) :: i
      type(htg_status_t) :: arr_st

      if (min_tuples < 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_cell_data_t%ensure_tuples: min_tuples must be >= 0', status)
         return
      end if

      do i = 0_i32_k, self%num_arrays - 1_i32_k
         call self%arrays(i)%item%ensure_tuples(min_tuples, arr_st)
         if (.not. arr_st%ok()) then
            call htg_fail(arr_st%code, trim(arr_st%message), status)
            return
         end if
      end do

      self%num_tuples = max(self%num_tuples, min_tuples)
      call htg_set_ok(status)
   end subroutine cell_data_ensure_tuples

   !> Copy tuple src_id onto tuple dst_id in every array that holds both.
   !> Arrays too short to hold either tuple are left untouched.
   pure subroutine cell_data_copy_tuple(self, dst_id, src_id)
      class(htg_cell_data_t), intent(inout) :: self
      integer(id_k), intent(in) :: dst_id, src_id

      integer(i32_k) :: i

      if (dst_id < 0_id_k .or. src_id < 0_id_k) return
      do i = 0_i32_k, self%num_arrays - 1_i32_k
         if (dst_id < self%arrays(i)%item%num_tuples .and. &
             src_id < self%arrays(i)%item%num_tuples) then
            call self%arrays(i)%item%copy_tuple(dst_id, src_id)
         end if
      end do
   end subroutine cell_data_copy_tuple

   !> Reset a tuple to zero in every array that holds it.
   pure subroutine cell_data_clear_tuple(self, tuple_id)
      class(htg_cell_data_t), intent(inout) :: self
      integer(id_k), intent(in) :: tuple_id

      integer(i32_k) :: i

      if (tuple_id < 0_id_k) return
      do i = 0_i32_k, self%num_arrays - 1_i32_k
         if (tuple_id < self%arrays(i)%item%num_tuples) then
            call self%arrays(i)%item%clear_tuple(tuple_id)
         end if
      end do
   end subroutine cell_data_clear_tuple

   pure subroutine cell_data_set_number_of_tuples(self, num_tuples, status)
      class(htg_cell_data_t), intent(inout) :: self
      integer(id_k), intent(in) :: num_tuples
      type(htg_status_t), intent(out), optional :: status

      integer(i32_k) :: i
      type(htg_status_t) :: arr_st

      if (num_tuples < 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_cell_data_t%set_number_of_tuples: num_tuples must be >= 0', status)
         return
      end if

      do i = 0_i32_k, self%num_arrays - 1_i32_k
         call self%arrays(i)%item%set_number_of_tuples(num_tuples, arr_st)
         if (.not. arr_st%ok()) then
            call htg_fail(arr_st%code, trim(arr_st%message), status)
            return
         end if
      end do

      self%num_tuples = num_tuples
      call htg_set_ok(status)
   end subroutine cell_data_set_number_of_tuples

end module fvtkhtg_cell_data_m
