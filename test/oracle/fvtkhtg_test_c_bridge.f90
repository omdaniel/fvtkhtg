! Test-only Fortran 2018 `iso_c_binding` bridge for C++ VTK 9.7+ differential
! and LLDB lockstep verification.
!
! IMPORTANT: `libfvtkhtg` in `src/` is 100% pure Fortran 2018 with ZERO C/C++
! bindings. This module lives strictly inside `test/oracle/` and is compiled
! only when `-DFVTKHTG_BUILD_ORACLE_TESTS=ON`.
module fvtkhtg_test_c_bridge_m
   use, intrinsic :: iso_c_binding, only: &
      c_ptr, c_null_ptr, c_associated, c_loc, c_f_pointer, &
      c_bool, c_char, c_null_char, c_int8_t, c_int32_t, c_int64_t, &
      c_float, c_double
   use fvtkhtg, only: &
      id_k, i64_k, i32_k, byte_k, r64_k, r32_k, &
      HTG_OK, HTG_ERR_INVALID_ARG, HTG_ERR_INVALID_STATE, HTG_INVALID_ID, &
      htg_status_t, htg_bit_array_t, htg_field_array_t, htg_cell_data_t, &
      htg_tree_t, htg_grid_t, htg_cursor_t
   implicit none
   private

   type :: i64_vec_box_t
      integer(id_k), allocatable :: data(:)
   end type i64_vec_box_t

contains

   pure function c_chars_to_f_str(c_chars, n) result(f_str)
      character(kind=c_char), intent(in) :: c_chars(*)
      integer(c_int32_t), intent(in) :: n
      character(len=max(0, int(n))) :: f_str
      integer :: i

      do i = 1, max(0, int(n))
         f_str(i:i) = c_chars(i)
      end do
   end function c_chars_to_f_str

   subroutine write_opt_i32(ptr, val)
      type(c_ptr), value :: ptr
      integer(c_int32_t), intent(in) :: val
      integer(c_int32_t), pointer :: fp

      if (c_associated(ptr)) then
         call c_f_pointer(ptr, fp)
         fp = val
      end if
   end subroutine write_opt_i32

   subroutine write_opt_i64(ptr, val)
      type(c_ptr), value :: ptr
      integer(c_int64_t), intent(in) :: val
      integer(c_int64_t), pointer :: fp

      if (c_associated(ptr)) then
         call c_f_pointer(ptr, fp)
         fp = val
      end if
   end subroutine write_opt_i64

   ! =========================================================================
   ! 1. BitArray (`htg_bit_array_t`)
   ! =========================================================================

   function ftn_htg_bit_array_create(num_bits, init_val, status_out) &
         result(res) bind(C, name="ftn_htg_bit_array_create")
      integer(c_int64_t), value :: num_bits
      logical(c_bool), value :: init_val
      type(c_ptr), value :: status_out
      type(c_ptr) :: res

      type(htg_bit_array_t), pointer :: bptr
      type(htg_status_t) :: st

      allocate(bptr)
      call bptr%initialize(int(num_bits, id_k), logical(init_val), st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      if (.not. st%ok()) then
         deallocate(bptr)
         res = c_null_ptr
      else
         res = c_loc(bptr)
      end if
   end function ftn_htg_bit_array_create

   subroutine ftn_htg_bit_array_destroy(bit_arr) &
         bind(C, name="ftn_htg_bit_array_destroy")
      type(c_ptr), value :: bit_arr
      type(htg_bit_array_t), pointer :: bptr

      if (.not. c_associated(bit_arr)) return
      call c_f_pointer(bit_arr, bptr)
      call bptr%clear()
      deallocate(bptr)
   end subroutine ftn_htg_bit_array_destroy

   function ftn_htg_bit_array_resize(bit_arr, num_bits) &
         result(code) bind(C, name="ftn_htg_bit_array_resize")
      type(c_ptr), value :: bit_arr
      integer(c_int64_t), value :: num_bits
      integer(c_int32_t) :: code
      type(htg_bit_array_t), pointer :: bptr
      type(htg_status_t) :: st

      call c_f_pointer(bit_arr, bptr)
      call bptr%resize(int(num_bits, id_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_bit_array_resize

   function ftn_htg_bit_array_get_number_of_values(bit_arr) &
         result(nbits) bind(C, name="ftn_htg_bit_array_get_number_of_values")
      type(c_ptr), value :: bit_arr
      integer(c_int64_t) :: nbits
      type(htg_bit_array_t), pointer :: bptr

      call c_f_pointer(bit_arr, bptr)
      nbits = int(bptr%get_number_of_values(), c_int64_t)
   end function ftn_htg_bit_array_get_number_of_values

   function ftn_htg_bit_array_get_number_of_bytes(bit_arr) &
         result(nbytes) bind(C, name="ftn_htg_bit_array_get_number_of_bytes")
      type(c_ptr), value :: bit_arr
      integer(c_int64_t) :: nbytes
      type(htg_bit_array_t), pointer :: bptr

      call c_f_pointer(bit_arr, bptr)
      nbytes = int(bptr%get_number_of_bytes(), c_int64_t)
   end function ftn_htg_bit_array_get_number_of_bytes

   function ftn_htg_bit_array_set_value(bit_arr, id, val) &
         result(code) bind(C, name="ftn_htg_bit_array_set_value")
      type(c_ptr), value :: bit_arr
      integer(c_int64_t), value :: id
      logical(c_bool), value :: val
      integer(c_int32_t) :: code
      type(htg_bit_array_t), pointer :: bptr
      type(htg_status_t) :: st

      call c_f_pointer(bit_arr, bptr)
      call bptr%set_value(int(id, id_k), logical(val), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_bit_array_set_value

   function ftn_htg_bit_array_insert_value(bit_arr, id, val) &
         result(code) bind(C, name="ftn_htg_bit_array_insert_value")
      type(c_ptr), value :: bit_arr
      integer(c_int64_t), value :: id
      logical(c_bool), value :: val
      integer(c_int32_t) :: code
      type(htg_bit_array_t), pointer :: bptr
      type(htg_status_t) :: st

      call c_f_pointer(bit_arr, bptr)
      call bptr%insert_value(int(id, id_k), logical(val), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_bit_array_insert_value

   function ftn_htg_bit_array_insert_next_value(bit_arr, val, inserted_id_out) &
         result(code) bind(C, name="ftn_htg_bit_array_insert_next_value")
      type(c_ptr), value :: bit_arr
      logical(c_bool), value :: val
      type(c_ptr), value :: inserted_id_out
      integer(c_int32_t) :: code
      type(htg_bit_array_t), pointer :: bptr
      integer(id_k) :: ins_id
      type(htg_status_t) :: st

      call c_f_pointer(bit_arr, bptr)
      call bptr%insert_next_value(logical(val), ins_id, st)
      call write_opt_i64(inserted_id_out, int(ins_id, c_int64_t))
      code = int(st%code, c_int32_t)
   end function ftn_htg_bit_array_insert_next_value

   function ftn_htg_bit_array_get_value(bit_arr, id, status_out) &
         result(res) bind(C, name="ftn_htg_bit_array_get_value")
      type(c_ptr), value :: bit_arr
      integer(c_int64_t), value :: id
      type(c_ptr), value :: status_out
      logical(c_bool) :: res
      type(htg_bit_array_t), pointer :: bptr
      logical :: fval
      type(htg_status_t) :: st

      call c_f_pointer(bit_arr, bptr)
      fval = bptr%get_value(int(id, id_k), st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      res = logical(fval, kind=c_bool)
   end function ftn_htg_bit_array_get_value

   subroutine ftn_htg_bit_array_fill(bit_arr, val) &
         bind(C, name="ftn_htg_bit_array_fill")
      type(c_ptr), value :: bit_arr
      logical(c_bool), value :: val
      type(htg_bit_array_t), pointer :: bptr

      call c_f_pointer(bit_arr, bptr)
      call bptr%fill(logical(val))
   end subroutine ftn_htg_bit_array_fill

   subroutine ftn_htg_bit_array_get_raw_bytes(bit_arr, ptr_out, nbytes_out) &
         bind(C, name="ftn_htg_bit_array_get_raw_bytes")
      type(c_ptr), value :: bit_arr
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: nbytes_out
      type(htg_bit_array_t), pointer :: bptr
      integer(byte_k), pointer :: raw(:)

      call c_f_pointer(bit_arr, bptr)
      call bptr%get_raw_bytes(raw)
      if (associated(raw) .and. size(raw) > 0) then
         ptr_out = c_loc(raw(0_id_k))
         nbytes_out = int(size(raw, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         nbytes_out = 0_c_int64_t
      end if
   end subroutine ftn_htg_bit_array_get_raw_bytes

   ! =========================================================================
   ! 2. Grid (`htg_grid_t`) Lifecycle, Dimensions, Coordinates & Root Indexing
   ! =========================================================================

   function ftn_htg_grid_create() result(res) bind(C, name="ftn_htg_grid_create")
      type(c_ptr) :: res
      type(htg_grid_t), pointer :: gptr

      allocate(gptr)
      call gptr%initialize()
      res = c_loc(gptr)
   end function ftn_htg_grid_create

   subroutine ftn_htg_grid_destroy(grid) bind(C, name="ftn_htg_grid_destroy")
      type(c_ptr), value :: grid
      type(htg_grid_t), pointer :: gptr

      if (.not. c_associated(grid)) return
      call c_f_pointer(grid, gptr)
      call gptr%clear()
      deallocate(gptr)
   end subroutine ftn_htg_grid_destroy

   subroutine ftn_htg_grid_initialize(grid) bind(C, name="ftn_htg_grid_initialize")
      type(c_ptr), value :: grid
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      call gptr%initialize()
   end subroutine ftn_htg_grid_initialize

   function ftn_htg_grid_set_dimensions(grid, ni, nj, nk) &
         result(code) bind(C, name="ftn_htg_grid_set_dimensions")
      type(c_ptr), value :: grid
      integer(c_int32_t), value :: ni, nj, nk
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%set_dimensions(int(ni, i32_k), int(nj, i32_k), int(nk, i32_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_set_dimensions

   subroutine ftn_htg_grid_get_dimensions(grid, dims_out) &
         bind(C, name="ftn_htg_grid_get_dimensions")
      type(c_ptr), value :: grid
      integer(c_int32_t), intent(out) :: dims_out(0:2)
      type(htg_grid_t), pointer :: gptr
      integer(i32_k) :: d(0:2)

      call c_f_pointer(grid, gptr)
      d = gptr%get_dimensions()
      dims_out(0:2) = int(d(0:2), c_int32_t)
   end subroutine ftn_htg_grid_get_dimensions

   subroutine ftn_htg_grid_get_cell_dims(grid, cell_dims_out) &
         bind(C, name="ftn_htg_grid_get_cell_dims")
      type(c_ptr), value :: grid
      integer(c_int32_t), intent(out) :: cell_dims_out(0:2)
      type(htg_grid_t), pointer :: gptr
      integer(i32_k) :: cd(0:2)

      call c_f_pointer(grid, gptr)
      cd = gptr%get_cell_dims()
      cell_dims_out(0:2) = int(cd(0:2), c_int32_t)
   end subroutine ftn_htg_grid_get_cell_dims

   function ftn_htg_grid_get_dimension(grid) &
         result(ndim) bind(C, name="ftn_htg_grid_get_dimension")
      type(c_ptr), value :: grid
      integer(c_int32_t) :: ndim
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      ndim = int(gptr%get_dimension(), c_int32_t)
   end function ftn_htg_grid_get_dimension

   function ftn_htg_grid_get_orientation(grid) &
         result(orient) bind(C, name="ftn_htg_grid_get_orientation")
      type(c_ptr), value :: grid
      integer(c_int32_t) :: orient
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      orient = int(gptr%get_orientation(), c_int32_t)
   end function ftn_htg_grid_get_orientation

   subroutine ftn_htg_grid_get_axes(grid, axes_out, naxes_out) &
         bind(C, name="ftn_htg_grid_get_axes")
      type(c_ptr), value :: grid
      integer(c_int32_t), intent(out) :: axes_out(0:2)
      integer(c_int32_t), intent(out) :: naxes_out
      type(htg_grid_t), pointer :: gptr
      integer(i32_k) :: ax(0:2)

      call c_f_pointer(grid, gptr)
      ax = gptr%get_axes()
      axes_out(0:2) = int(ax(0:2), c_int32_t)
      naxes_out = int(gptr%get_dimension(), c_int32_t)
   end subroutine ftn_htg_grid_get_axes

   function ftn_htg_grid_set_branch_factor(grid, factor) &
         result(code) bind(C, name="ftn_htg_grid_set_branch_factor")
      type(c_ptr), value :: grid
      integer(c_int32_t), value :: factor
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%set_branch_factor(int(factor, i32_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_set_branch_factor

   function ftn_htg_grid_set_depth_limiter(grid, limit) &
         result(code) bind(C, name="ftn_htg_grid_set_depth_limiter")
      type(c_ptr), value :: grid
      integer(c_int32_t), value :: limit
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%set_depth_limiter(int(limit, i32_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_set_depth_limiter

   function ftn_htg_grid_get_depth_limiter(grid) &
         result(limit) bind(C, name="ftn_htg_grid_get_depth_limiter")
      type(c_ptr), value :: grid
      integer(c_int32_t) :: limit
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      limit = int(gptr%get_depth_limiter(), c_int32_t)
   end function ftn_htg_grid_get_depth_limiter

   function ftn_htg_grid_get_branch_factor(grid) &
         result(factor) bind(C, name="ftn_htg_grid_get_branch_factor")
      type(c_ptr), value :: grid
      integer(c_int32_t) :: factor
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      factor = int(gptr%get_branch_factor(), c_int32_t)
   end function ftn_htg_grid_get_branch_factor

   function ftn_htg_grid_get_number_of_children(grid) &
         result(nchild) bind(C, name="ftn_htg_grid_get_number_of_children")
      type(c_ptr), value :: grid
      integer(c_int32_t) :: nchild
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      nchild = int(gptr%get_number_of_children(), c_int32_t)
   end function ftn_htg_grid_get_number_of_children

   subroutine ftn_htg_grid_set_transposed_root_indexing(grid, flag) &
         bind(C, name="ftn_htg_grid_set_transposed_root_indexing")
      type(c_ptr), value :: grid
      logical(c_bool), value :: flag
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%set_transposed_root_indexing(logical(flag), st)
   end subroutine ftn_htg_grid_set_transposed_root_indexing

   function ftn_htg_grid_get_transposed_root_indexing(grid) &
         result(flag) bind(C, name="ftn_htg_grid_get_transposed_root_indexing")
      type(c_ptr), value :: grid
      logical(c_bool) :: flag
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      flag = logical(gptr%get_transposed_root_indexing(), kind=c_bool)
   end function ftn_htg_grid_get_transposed_root_indexing

   function ftn_htg_grid_set_x_coordinates(grid, coords_ptr, n) &
         result(code) bind(C, name="ftn_htg_grid_set_x_coordinates")
      type(c_ptr), value :: grid
      type(c_ptr), value :: coords_ptr
      integer(c_int64_t), value :: n
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      real(r64_k), pointer :: fcoords(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      if (.not. c_associated(coords_ptr) .or. n <= 0_c_int64_t) then
         code = int(HTG_ERR_INVALID_ARG, c_int32_t)
         return
      end if
      call c_f_pointer(coords_ptr, fcoords, [n])
      call gptr%set_x_coordinates(fcoords(1:n), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_set_x_coordinates

   function ftn_htg_grid_set_y_coordinates(grid, coords_ptr, n) &
         result(code) bind(C, name="ftn_htg_grid_set_y_coordinates")
      type(c_ptr), value :: grid
      type(c_ptr), value :: coords_ptr
      integer(c_int64_t), value :: n
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      real(r64_k), pointer :: fcoords(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      if (.not. c_associated(coords_ptr) .or. n <= 0_c_int64_t) then
         code = int(HTG_ERR_INVALID_ARG, c_int32_t)
         return
      end if
      call c_f_pointer(coords_ptr, fcoords, [n])
      call gptr%set_y_coordinates(fcoords(1:n), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_set_y_coordinates

   function ftn_htg_grid_set_z_coordinates(grid, coords_ptr, n) &
         result(code) bind(C, name="ftn_htg_grid_set_z_coordinates")
      type(c_ptr), value :: grid
      type(c_ptr), value :: coords_ptr
      integer(c_int64_t), value :: n
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      real(r64_k), pointer :: fcoords(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      if (.not. c_associated(coords_ptr) .or. n <= 0_c_int64_t) then
         code = int(HTG_ERR_INVALID_ARG, c_int32_t)
         return
      end if
      call c_f_pointer(coords_ptr, fcoords, [n])
      call gptr%set_z_coordinates(fcoords(1:n), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_set_z_coordinates

   subroutine ftn_htg_grid_get_x_coordinates_view(grid, ptr_out, len_out) &
         bind(C, name="ftn_htg_grid_get_x_coordinates_view")
      type(c_ptr), value :: grid
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: len_out
      type(htg_grid_t), pointer :: gptr
      real(r64_k), pointer :: view(:)

      call c_f_pointer(grid, gptr)
      call gptr%get_x_coordinates_view(view)
      if (associated(view) .and. size(view) > 0) then
         ptr_out = c_loc(view(0_id_k))
         len_out = int(size(view, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         len_out = 0_c_int64_t
      end if
   end subroutine ftn_htg_grid_get_x_coordinates_view

   subroutine ftn_htg_grid_get_y_coordinates_view(grid, ptr_out, len_out) &
         bind(C, name="ftn_htg_grid_get_y_coordinates_view")
      type(c_ptr), value :: grid
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: len_out
      type(htg_grid_t), pointer :: gptr
      real(r64_k), pointer :: view(:)

      call c_f_pointer(grid, gptr)
      call gptr%get_y_coordinates_view(view)
      if (associated(view) .and. size(view) > 0) then
         ptr_out = c_loc(view(0_id_k))
         len_out = int(size(view, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         len_out = 0_c_int64_t
      end if
   end subroutine ftn_htg_grid_get_y_coordinates_view

   subroutine ftn_htg_grid_get_z_coordinates_view(grid, ptr_out, len_out) &
         bind(C, name="ftn_htg_grid_get_z_coordinates_view")
      type(c_ptr), value :: grid
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: len_out
      type(htg_grid_t), pointer :: gptr
      real(r64_k), pointer :: view(:)

      call c_f_pointer(grid, gptr)
      call gptr%get_z_coordinates_view(view)
      if (associated(view) .and. size(view) > 0) then
         ptr_out = c_loc(view(0_id_k))
         len_out = int(size(view, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         len_out = 0_c_int64_t
      end if
   end subroutine ftn_htg_grid_get_z_coordinates_view

   function ftn_htg_grid_get_index_from_level_zero_coordinates( &
         grid, i, j, k, tree_id_out) result(code) &
         bind(C, name="ftn_htg_grid_get_index_from_level_zero_coordinates")
      type(c_ptr), value :: grid
      integer(c_int32_t), value :: i, j, k
      integer(c_int64_t), intent(out) :: tree_id_out
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      integer(id_k) :: tid
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_index_from_level_zero_coordinates( &
         tid, int(i, i32_k), int(j, i32_k), int(k, i32_k), st)
      tree_id_out = int(tid, c_int64_t)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_get_index_from_level_zero_coordinates

   function ftn_htg_grid_get_level_zero_coordinates_from_index( &
         grid, tree_id, i_out, j_out, k_out) result(code) &
         bind(C, name="ftn_htg_grid_get_level_zero_coordinates_from_index")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: tree_id
      integer(c_int32_t), intent(out) :: i_out, j_out, k_out
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      integer(i32_k) :: fi, fj, fk
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_level_zero_coordinates_from_index( &
         int(tree_id, id_k), fi, fj, fk, st)
      i_out = int(fi, c_int32_t)
      j_out = int(fj, c_int32_t)
      k_out = int(fk, c_int32_t)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_get_level_zero_coordinates_from_index

   ! =========================================================================
   ! 3. Sparse Forest & Aggregations
   ! =========================================================================

   function ftn_htg_grid_get_max_number_of_trees(grid) &
         result(res) bind(C, name="ftn_htg_grid_get_max_number_of_trees")
      type(c_ptr), value :: grid
      integer(c_int64_t) :: res
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      res = int(gptr%get_max_number_of_trees(), c_int64_t)
   end function ftn_htg_grid_get_max_number_of_trees

   function ftn_htg_grid_get_number_of_non_empty_trees(grid) &
         result(res) bind(C, name="ftn_htg_grid_get_number_of_non_empty_trees")
      type(c_ptr), value :: grid
      integer(c_int64_t) :: res
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      res = int(gptr%get_number_of_non_empty_trees(), c_int64_t)
   end function ftn_htg_grid_get_number_of_non_empty_trees

   function ftn_htg_grid_get_number_of_vertices(grid) &
         result(res) bind(C, name="ftn_htg_grid_get_number_of_vertices")
      type(c_ptr), value :: grid
      integer(c_int64_t) :: res
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      res = int(gptr%get_number_of_vertices(), c_int64_t)
   end function ftn_htg_grid_get_number_of_vertices

   function ftn_htg_grid_get_number_of_leaves(grid) &
         result(res) bind(C, name="ftn_htg_grid_get_number_of_leaves")
      type(c_ptr), value :: grid
      integer(c_int64_t) :: res
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      res = int(gptr%get_number_of_leaves(), c_int64_t)
   end function ftn_htg_grid_get_number_of_leaves

   function ftn_htg_grid_get_number_of_levels(grid) &
         result(res) bind(C, name="ftn_htg_grid_get_number_of_levels")
      type(c_ptr), value :: grid
      integer(c_int64_t) :: res
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      res = int(gptr%get_number_of_levels(), c_int64_t)
   end function ftn_htg_grid_get_number_of_levels

   function ftn_htg_grid_get_global_node_index_max(grid) &
         result(res) bind(C, name="ftn_htg_grid_get_global_node_index_max")
      type(c_ptr), value :: grid
      integer(c_int64_t) :: res
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      res = int(gptr%get_global_node_index_max(), c_int64_t)
   end function ftn_htg_grid_get_global_node_index_max

   function ftn_htg_grid_initialize_local_index_node(grid) &
         result(code) bind(C, name="ftn_htg_grid_initialize_local_index_node")
      type(c_ptr), value :: grid
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%initialize_local_index_node(st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_initialize_local_index_node

   function ftn_htg_grid_has_tree(grid, tree_id) &
         result(flag) bind(C, name="ftn_htg_grid_has_tree")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: tree_id
      logical(c_bool) :: flag
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      flag = logical(gptr%has_tree(int(tree_id, id_k)), kind=c_bool)
   end function ftn_htg_grid_has_tree

   function ftn_htg_grid_get_tree(grid, tree_id, create, status_out) &
         result(res) bind(C, name="ftn_htg_grid_get_tree")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: tree_id
      logical(c_bool), value :: create
      type(c_ptr), value :: status_out
      type(c_ptr) :: res
      type(htg_grid_t), pointer :: gptr
      type(htg_tree_t), pointer :: tptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_tree(int(tree_id, id_k), tptr, create=logical(create), status=st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      if (associated(tptr)) then
         res = c_loc(tptr)
      else
         res = c_null_ptr
      end if
   end function ftn_htg_grid_get_tree

   function ftn_htg_grid_remove_tree(grid, tree_id) &
         result(code) bind(C, name="ftn_htg_grid_remove_tree")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: tree_id
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%remove_tree(int(tree_id, id_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_remove_tree

   function ftn_htg_grid_get_tree_by_pos(grid, pos, tree_id_out, tree_ptr_out) &
         result(code) bind(C, name="ftn_htg_grid_get_tree_by_pos")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: pos
      integer(c_int64_t), intent(out) :: tree_id_out
      type(c_ptr), intent(out) :: tree_ptr_out
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_tree_t), pointer :: tptr
      integer(id_k) :: tid
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_tree_by_pos(int(pos, id_k), tid, tptr, st)
      tree_id_out = int(tid, c_int64_t)
      if (associated(tptr)) then
         tree_ptr_out = c_loc(tptr)
      else
         tree_ptr_out = c_null_ptr
      end if
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_get_tree_by_pos

   ! =========================================================================
   ! 4. Tree (`htg_tree_t`) Topology, Global Indexing, BF Descriptors & Views
   ! =========================================================================

   function ftn_htg_tree_get_tree_index(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_tree_index")
      type(c_ptr), value :: tree
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_tree_index(), c_int64_t)
   end function ftn_htg_tree_get_tree_index

   function ftn_htg_tree_get_branch_factor(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_branch_factor")
      type(c_ptr), value :: tree
      integer(c_int32_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_branch_factor(), c_int32_t)
   end function ftn_htg_tree_get_branch_factor

   function ftn_htg_tree_get_dimension(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_dimension")
      type(c_ptr), value :: tree
      integer(c_int32_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_dimension(), c_int32_t)
   end function ftn_htg_tree_get_dimension

   function ftn_htg_tree_get_number_of_children(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_number_of_children")
      type(c_ptr), value :: tree
      integer(c_int32_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_number_of_children(), c_int32_t)
   end function ftn_htg_tree_get_number_of_children

   function ftn_htg_tree_get_number_of_levels(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_number_of_levels")
      type(c_ptr), value :: tree
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_number_of_levels(), c_int64_t)
   end function ftn_htg_tree_get_number_of_levels

   function ftn_htg_tree_get_number_of_vertices(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_number_of_vertices")
      type(c_ptr), value :: tree
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_number_of_vertices(), c_int64_t)
   end function ftn_htg_tree_get_number_of_vertices

   function ftn_htg_tree_get_number_of_nodes(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_number_of_nodes")
      type(c_ptr), value :: tree
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_number_of_nodes(), c_int64_t)
   end function ftn_htg_tree_get_number_of_nodes

   function ftn_htg_tree_get_number_of_leaves(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_number_of_leaves")
      type(c_ptr), value :: tree
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_number_of_leaves(), c_int64_t)
   end function ftn_htg_tree_get_number_of_leaves

   function ftn_htg_tree_get_number_of_vertices_at_depth(tree, depth, status_out) &
         result(res) bind(C, name="ftn_htg_tree_get_number_of_vertices_at_depth")
      type(c_ptr), value :: tree
      integer(c_int32_t), value :: depth
      type(c_ptr), value :: status_out
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_number_of_vertices_at_depth(int(depth, i32_k), st), c_int64_t)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
   end function ftn_htg_tree_get_number_of_vertices_at_depth

   function ftn_htg_tree_is_leaf(tree, local_id, status_out) &
         result(flag) bind(C, name="ftn_htg_tree_is_leaf")
      type(c_ptr), value :: tree
      integer(c_int64_t), value :: local_id
      type(c_ptr), value :: status_out
      logical(c_bool) :: flag
      type(htg_tree_t), pointer :: tptr
      logical :: fval
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      fval = tptr%is_leaf(int(local_id, id_k), st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      flag = logical(fval, kind=c_bool)
   end function ftn_htg_tree_is_leaf

   function ftn_htg_tree_is_terminal_node(tree, local_id, status_out) &
         result(flag) bind(C, name="ftn_htg_tree_is_terminal_node")
      type(c_ptr), value :: tree
      integer(c_int64_t), value :: local_id
      type(c_ptr), value :: status_out
      logical(c_bool) :: flag
      type(htg_tree_t), pointer :: tptr
      logical :: fval
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      fval = tptr%is_terminal_node(int(local_id, id_k), st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      flag = logical(fval, kind=c_bool)
   end function ftn_htg_tree_is_terminal_node

   function ftn_htg_tree_get_elder_child_index(tree, local_id, status_out) &
         result(res) bind(C, name="ftn_htg_tree_get_elder_child_index")
      type(c_ptr), value :: tree
      integer(c_int64_t), value :: local_id
      type(c_ptr), value :: status_out
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_elder_child_index(int(local_id, id_k), st), c_int64_t)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
   end function ftn_htg_tree_get_elder_child_index

   function ftn_htg_tree_subdivide_leaf(tree, local_id, level) &
         result(code) bind(C, name="ftn_htg_tree_subdivide_leaf")
      type(c_ptr), value :: tree
      integer(c_int64_t), value :: local_id
      integer(c_int32_t), value :: level
      integer(c_int32_t) :: code
      type(htg_tree_t), pointer :: tptr
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      call tptr%subdivide_leaf(int(local_id, id_k), int(level, i32_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_tree_subdivide_leaf

   function ftn_htg_tree_set_global_index_start(tree, start_id) &
         result(code) bind(C, name="ftn_htg_tree_set_global_index_start")
      type(c_ptr), value :: tree
      integer(c_int64_t), value :: start_id
      integer(c_int32_t) :: code
      type(htg_tree_t), pointer :: tptr
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      call tptr%set_global_index_start(int(start_id, id_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_tree_set_global_index_start

   function ftn_htg_tree_get_global_index_start(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_global_index_start")
      type(c_ptr), value :: tree
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_global_index_start(), c_int64_t)
   end function ftn_htg_tree_get_global_index_start

   function ftn_htg_tree_set_global_index_from_local(tree, local_id, global_id) &
         result(code) bind(C, name="ftn_htg_tree_set_global_index_from_local")
      type(c_ptr), value :: tree
      integer(c_int64_t), value :: local_id, global_id
      integer(c_int32_t) :: code
      type(htg_tree_t), pointer :: tptr
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      call tptr%set_global_index_from_local(int(local_id, id_k), int(global_id, id_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_tree_set_global_index_from_local

   function ftn_htg_tree_get_global_index_from_local(tree, local_id, status_out) &
         result(res) bind(C, name="ftn_htg_tree_get_global_index_from_local")
      type(c_ptr), value :: tree
      integer(c_int64_t), value :: local_id
      type(c_ptr), value :: status_out
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_global_index_from_local(int(local_id, id_k), st), c_int64_t)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
   end function ftn_htg_tree_get_global_index_from_local

   function ftn_htg_tree_get_indexing_mode(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_indexing_mode")
      type(c_ptr), value :: tree
      integer(c_int32_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_indexing_mode(), c_int32_t)
   end function ftn_htg_tree_get_indexing_mode

   function ftn_htg_tree_get_max_global_index(tree) &
         result(res) bind(C, name="ftn_htg_tree_get_max_global_index")
      type(c_ptr), value :: tree
      integer(c_int64_t) :: res
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(tree, tptr)
      res = int(tptr%get_max_global_index(), c_int64_t)
   end function ftn_htg_tree_get_max_global_index

   subroutine ftn_htg_tree_get_parent_to_elder_child_view(tree, ptr_out, len_out) &
         bind(C, name="ftn_htg_tree_get_parent_to_elder_child_view")
      type(c_ptr), value :: tree
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: len_out
      type(htg_tree_t), pointer :: tptr
      integer(id_k), pointer :: view(:)

      call c_f_pointer(tree, tptr)
      call tptr%get_parent_to_elder_child_view(view)
      if (associated(view) .and. size(view) > 0) then
         ptr_out = c_loc(view(0_id_k))
         len_out = int(size(view, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         len_out = 0_c_int64_t
      end if
   end subroutine ftn_htg_tree_get_parent_to_elder_child_view

   subroutine ftn_htg_tree_get_vertices_per_depth_view(tree, ptr_out, len_out) &
         bind(C, name="ftn_htg_tree_get_vertices_per_depth_view")
      type(c_ptr), value :: tree
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: len_out
      type(htg_tree_t), pointer :: tptr
      integer(id_k), pointer :: view(:)

      call c_f_pointer(tree, tptr)
      call tptr%get_vertices_per_depth_view(view)
      if (associated(view) .and. size(view) > 0) then
         ptr_out = c_loc(view(0_id_k))
         len_out = int(size(view, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         len_out = 0_c_int64_t
      end if
   end subroutine ftn_htg_tree_get_vertices_per_depth_view

   subroutine ftn_htg_tree_get_global_index_table_view(tree, ptr_out, len_out) &
         bind(C, name="ftn_htg_tree_get_global_index_table_view")
      type(c_ptr), value :: tree
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: len_out
      type(htg_tree_t), pointer :: tptr
      integer(id_k), pointer :: view(:)

      call c_f_pointer(tree, tptr)
      call tptr%get_global_index_table_view(view)
      if (associated(view) .and. size(view) > 0) then
         ptr_out = c_loc(view(0_id_k))
         len_out = int(size(view, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         len_out = 0_c_int64_t
      end if
   end subroutine ftn_htg_tree_get_global_index_table_view

   function ftn_htg_tree_build_from_bf_descriptor( &
         tree, descriptor_bit_arr, num_bits, start_index) result(code) &
         bind(C, name="ftn_htg_tree_build_from_bf_descriptor")
      type(c_ptr), value :: tree
      type(c_ptr), value :: descriptor_bit_arr
      integer(c_int64_t), value :: num_bits
      integer(c_int64_t), value :: start_index
      integer(c_int32_t) :: code
      type(htg_tree_t), pointer :: tptr
      type(htg_bit_array_t), pointer :: bptr
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      call c_f_pointer(descriptor_bit_arr, bptr)
      call tptr%build_from_breadth_first_order_descriptor( &
         bptr, int(num_bits, id_k), start_index=int(start_index, id_k), status=st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_tree_build_from_bf_descriptor

   function ftn_htg_i64_vec_create() result(res) bind(C, name="ftn_htg_i64_vec_create")
      type(c_ptr) :: res
      type(i64_vec_box_t), pointer :: vptr

      allocate(vptr)
      res = c_loc(vptr)
   end function ftn_htg_i64_vec_create

   subroutine ftn_htg_i64_vec_destroy(vec) bind(C, name="ftn_htg_i64_vec_destroy")
      type(c_ptr), value :: vec
      type(i64_vec_box_t), pointer :: vptr

      if (.not. c_associated(vec)) return
      call c_f_pointer(vec, vptr)
      if (allocated(vptr%data)) deallocate(vptr%data)
      deallocate(vptr)
   end subroutine ftn_htg_i64_vec_destroy

   subroutine ftn_htg_i64_vec_get_view(vec, ptr_out, len_out) &
         bind(C, name="ftn_htg_i64_vec_get_view")
      type(c_ptr), value :: vec
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: len_out
      type(i64_vec_box_t), pointer :: vptr
      integer(id_k) :: lb

      call c_f_pointer(vec, vptr)
      if (allocated(vptr%data) .and. size(vptr%data) > 0) then
         lb = lbound(vptr%data, 1, kind=id_k)
         ptr_out = c_loc(vptr%data(lb))
         len_out = int(size(vptr%data, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         len_out = 0_c_int64_t
      end if
   end subroutine ftn_htg_i64_vec_get_view

   function ftn_htg_tree_compute_bf_descriptor_append( &
         tree, inout_descriptor_bit_arr, inout_vpd_vec, inout_bfid_vec, &
         opt_input_mask_bit_arr, has_depth_limiter, depth_limiter) result(code) &
         bind(C, name="ftn_htg_tree_compute_bf_descriptor_append")
      type(c_ptr), value :: tree
      type(c_ptr), value :: inout_descriptor_bit_arr
      type(c_ptr), value :: inout_vpd_vec
      type(c_ptr), value :: inout_bfid_vec
      type(c_ptr), value :: opt_input_mask_bit_arr
      logical(c_bool), value :: has_depth_limiter
      integer(c_int32_t), value :: depth_limiter
      integer(c_int32_t) :: code

      type(htg_tree_t), pointer :: tptr
      type(htg_bit_array_t), pointer :: desc_ptr, mask_ptr
      type(i64_vec_box_t), pointer :: vpd_ptr, bfid_ptr
      type(htg_status_t) :: st

      call c_f_pointer(tree, tptr)
      call c_f_pointer(inout_descriptor_bit_arr, desc_ptr)
      call c_f_pointer(inout_vpd_vec, vpd_ptr)
      call c_f_pointer(inout_bfid_vec, bfid_ptr)

      if (c_associated(opt_input_mask_bit_arr)) then
         call c_f_pointer(opt_input_mask_bit_arr, mask_ptr)
         if (logical(has_depth_limiter)) then
            call tptr%compute_breadth_first_order_descriptor( &
               desc_ptr, vpd_ptr%data, bfid_ptr%data, &
               input_mask=mask_ptr, depth_limiter=int(depth_limiter, i32_k), status=st)
         else
            call tptr%compute_breadth_first_order_descriptor( &
               desc_ptr, vpd_ptr%data, bfid_ptr%data, &
               input_mask=mask_ptr, status=st)
         end if
      else
         if (logical(has_depth_limiter)) then
            call tptr%compute_breadth_first_order_descriptor( &
               desc_ptr, vpd_ptr%data, bfid_ptr%data, &
               depth_limiter=int(depth_limiter, i32_k), status=st)
         else
            call tptr%compute_breadth_first_order_descriptor( &
               desc_ptr, vpd_ptr%data, bfid_ptr%data, status=st)
         end if
      end if

      code = int(st%code, c_int32_t)
   end function ftn_htg_tree_compute_bf_descriptor_append

   ! =========================================================================
   ! 5. Build-Time NonOrientedCursor (`htg_cursor_t`)
   ! =========================================================================

   function ftn_htg_cursor_create(grid, tree_id, create, status_out) &
         result(res) bind(C, name="ftn_htg_cursor_create")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: tree_id
      logical(c_bool), value :: create
      type(c_ptr), value :: status_out
      type(c_ptr) :: res
      type(htg_grid_t), pointer :: gptr
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      allocate(cptr)
      call cptr%initialize(gptr, int(tree_id, id_k), create=logical(create), status=st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      if (.not. st%ok()) then
         deallocate(cptr)
         res = c_null_ptr
      else
         res = c_loc(cptr)
      end if
   end function ftn_htg_cursor_create

   subroutine ftn_htg_cursor_destroy(cursor) bind(C, name="ftn_htg_cursor_destroy")
      type(c_ptr), value :: cursor
      type(htg_cursor_t), pointer :: cptr

      if (.not. c_associated(cursor)) return
      call c_f_pointer(cursor, cptr)
      call cptr%reset()
      deallocate(cptr)
   end subroutine ftn_htg_cursor_destroy

   function ftn_htg_cursor_initialize(cursor, grid, tree_id, create) &
         result(code) bind(C, name="ftn_htg_cursor_initialize")
      type(c_ptr), value :: cursor
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: tree_id
      logical(c_bool), value :: create
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      call c_f_pointer(grid, gptr)
      call cptr%initialize(gptr, int(tree_id, id_k), create=logical(create), status=st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_initialize

   function ftn_htg_cursor_has_tree(cursor) &
         result(flag) bind(C, name="ftn_htg_cursor_has_tree")
      type(c_ptr), value :: cursor
      logical(c_bool) :: flag
      type(htg_cursor_t), pointer :: cptr

      call c_f_pointer(cursor, cptr)
      flag = logical(cptr%has_tree(), kind=c_bool)
   end function ftn_htg_cursor_has_tree

   function ftn_htg_cursor_get_tree(cursor) &
         result(res) bind(C, name="ftn_htg_cursor_get_tree")
      type(c_ptr), value :: cursor
      type(c_ptr) :: res
      type(htg_cursor_t), pointer :: cptr
      type(htg_tree_t), pointer :: tptr

      call c_f_pointer(cursor, cptr)
      call cptr%get_tree(tptr)
      if (associated(tptr)) then
         res = c_loc(tptr)
      else
         res = c_null_ptr
      end if
   end function ftn_htg_cursor_get_tree

   function ftn_htg_cursor_get_vertex_id(cursor) &
         result(res) bind(C, name="ftn_htg_cursor_get_vertex_id")
      type(c_ptr), value :: cursor
      integer(c_int64_t) :: res
      type(htg_cursor_t), pointer :: cptr

      call c_f_pointer(cursor, cptr)
      res = int(cptr%get_vertex_id(), c_int64_t)
   end function ftn_htg_cursor_get_vertex_id

   function ftn_htg_cursor_get_global_node_index(cursor, status_out) &
         result(res) bind(C, name="ftn_htg_cursor_get_global_node_index")
      type(c_ptr), value :: cursor
      type(c_ptr), value :: status_out
      integer(c_int64_t) :: res
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      res = int(cptr%get_global_node_index(st), c_int64_t)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
   end function ftn_htg_cursor_get_global_node_index

   function ftn_htg_cursor_get_dimension(cursor) &
         result(res) bind(C, name="ftn_htg_cursor_get_dimension")
      type(c_ptr), value :: cursor
      integer(c_int32_t) :: res
      type(htg_cursor_t), pointer :: cptr

      call c_f_pointer(cursor, cptr)
      res = int(cptr%get_dimension(), c_int32_t)
   end function ftn_htg_cursor_get_dimension

   function ftn_htg_cursor_get_number_of_children(cursor) &
         result(res) bind(C, name="ftn_htg_cursor_get_number_of_children")
      type(c_ptr), value :: cursor
      integer(c_int32_t) :: res
      type(htg_cursor_t), pointer :: cptr

      call c_f_pointer(cursor, cptr)
      res = int(cptr%get_number_of_children(), c_int32_t)
   end function ftn_htg_cursor_get_number_of_children

   function ftn_htg_cursor_get_level(cursor) &
         result(res) bind(C, name="ftn_htg_cursor_get_level")
      type(c_ptr), value :: cursor
      integer(c_int32_t) :: res
      type(htg_cursor_t), pointer :: cptr

      call c_f_pointer(cursor, cptr)
      res = int(cptr%get_level(), c_int32_t)
   end function ftn_htg_cursor_get_level

   function ftn_htg_cursor_is_leaf(cursor, status_out) &
         result(flag) bind(C, name="ftn_htg_cursor_is_leaf")
      type(c_ptr), value :: cursor
      type(c_ptr), value :: status_out
      logical(c_bool) :: flag
      type(htg_cursor_t), pointer :: cptr
      logical :: fval
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      fval = cptr%is_leaf(st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      flag = logical(fval, kind=c_bool)
   end function ftn_htg_cursor_is_leaf

   function ftn_htg_cursor_is_terminal_node(cursor, status_out) &
         result(flag) bind(C, name="ftn_htg_cursor_is_terminal_node")
      type(c_ptr), value :: cursor
      type(c_ptr), value :: status_out
      logical(c_bool) :: flag
      type(htg_cursor_t), pointer :: cptr
      logical :: fval
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      fval = cptr%is_terminal_node(st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      flag = logical(fval, kind=c_bool)
   end function ftn_htg_cursor_is_terminal_node

   function ftn_htg_cursor_is_root(cursor) &
         result(flag) bind(C, name="ftn_htg_cursor_is_root")
      type(c_ptr), value :: cursor
      logical(c_bool) :: flag
      type(htg_cursor_t), pointer :: cptr

      call c_f_pointer(cursor, cptr)
      flag = logical(cptr%is_root(), kind=c_bool)
   end function ftn_htg_cursor_is_root

   function ftn_htg_cursor_to_root(cursor) &
         result(code) bind(C, name="ftn_htg_cursor_to_root")
      type(c_ptr), value :: cursor
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      call cptr%to_root(st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_to_root

   function ftn_htg_cursor_to_parent(cursor) &
         result(code) bind(C, name="ftn_htg_cursor_to_parent")
      type(c_ptr), value :: cursor
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      call cptr%to_parent(st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_to_parent

   function ftn_htg_cursor_to_child(cursor, ichild) &
         result(code) bind(C, name="ftn_htg_cursor_to_child")
      type(c_ptr), value :: cursor
      integer(c_int32_t), value :: ichild
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      call cptr%to_child(int(ichild, i32_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_to_child

   function ftn_htg_cursor_subdivide_leaf(cursor) &
         result(code) bind(C, name="ftn_htg_cursor_subdivide_leaf")
      type(c_ptr), value :: cursor
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      call cptr%subdivide_leaf(st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_subdivide_leaf

   function ftn_htg_cursor_set_global_index_start(cursor, offset) &
         result(code) bind(C, name="ftn_htg_cursor_set_global_index_start")
      type(c_ptr), value :: cursor
      integer(c_int64_t), value :: offset
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      call cptr%set_global_index_start(int(offset, id_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_set_global_index_start

   function ftn_htg_cursor_set_global_index_from_local(cursor, global_id) &
         result(code) bind(C, name="ftn_htg_cursor_set_global_index_from_local")
      type(c_ptr), value :: cursor
      integer(c_int64_t), value :: global_id
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      call cptr%set_global_index_from_local(int(global_id, id_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_set_global_index_from_local

   function ftn_htg_cursor_set_mask(cursor, state) &
         result(code) bind(C, name="ftn_htg_cursor_set_mask")
      type(c_ptr), value :: cursor
      logical(c_bool), value :: state
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      call cptr%set_mask(logical(state), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_set_mask

   function ftn_htg_cursor_is_masked(cursor, status_out) &
         result(flag) bind(C, name="ftn_htg_cursor_is_masked")
      type(c_ptr), value :: cursor
      type(c_ptr), value :: status_out
      logical(c_bool) :: flag
      type(htg_cursor_t), pointer :: cptr
      logical :: fval
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      fval = cptr%is_masked(st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      flag = logical(fval, kind=c_bool)
   end function ftn_htg_cursor_is_masked

   ! =========================================================================
   ! 6. Grid Mask, PureMask, CellData, Material Interfaces & Validate
   ! =========================================================================

   function ftn_htg_grid_set_mask(grid, bit_arr) &
         result(code) bind(C, name="ftn_htg_grid_set_mask")
      type(c_ptr), value :: grid
      type(c_ptr), value :: bit_arr
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_bit_array_t), pointer :: bptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      if (.not. c_associated(bit_arr)) then
         call gptr%clear_mask()
         code = int(HTG_OK, c_int32_t)
      else
         call c_f_pointer(bit_arr, bptr)
         call gptr%set_mask(bptr, st)
         code = int(st%code, c_int32_t)
      end if
   end function ftn_htg_grid_set_mask

   function ftn_htg_grid_clear_mask(grid) &
         result(code) bind(C, name="ftn_htg_grid_clear_mask")
      type(c_ptr), value :: grid
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      call gptr%clear_mask()
      code = int(HTG_OK, c_int32_t)
   end function ftn_htg_grid_clear_mask

   function ftn_htg_grid_has_mask(grid) &
         result(flag) bind(C, name="ftn_htg_grid_has_mask")
      type(c_ptr), value :: grid
      logical(c_bool) :: flag
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      flag = logical(gptr%has_mask(), kind=c_bool)
   end function ftn_htg_grid_has_mask

   subroutine ftn_htg_grid_get_mask_raw_bytes(grid, ptr_out, nbits_out, nbytes_out) &
         bind(C, name="ftn_htg_grid_get_mask_raw_bytes")
      type(c_ptr), value :: grid
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: nbits_out, nbytes_out
      type(htg_grid_t), pointer :: gptr
      type(htg_bit_array_t), pointer :: mptr
      integer(byte_k), pointer :: raw(:)

      call c_f_pointer(grid, gptr)
      call gptr%get_mask(mptr)
      if (.not. associated(mptr)) then
         ptr_out = c_null_ptr
         nbits_out = 0_c_int64_t
         nbytes_out = 0_c_int64_t
         return
      end if
      nbits_out = int(mptr%get_number_of_values(), c_int64_t)
      call mptr%get_raw_bytes(raw)
      if (associated(raw) .and. size(raw) > 0) then
         ptr_out = c_loc(raw(0_id_k))
         nbytes_out = int(size(raw, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         nbytes_out = 0_c_int64_t
      end if
   end subroutine ftn_htg_grid_get_mask_raw_bytes

   function ftn_htg_grid_set_cell_mask(grid, global_id, state) &
         result(code) bind(C, name="ftn_htg_grid_set_cell_mask")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: global_id
      logical(c_bool), value :: state
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%set_cell_mask(int(global_id, id_k), logical(state), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_set_cell_mask

   function ftn_htg_grid_is_cell_masked(grid, global_id, status_out) &
         result(flag) bind(C, name="ftn_htg_grid_is_cell_masked")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: global_id
      type(c_ptr), value :: status_out
      logical(c_bool) :: flag
      type(htg_grid_t), pointer :: gptr
      logical :: fval
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      fval = gptr%is_cell_masked(int(global_id, id_k), st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      flag = logical(fval, kind=c_bool)
   end function ftn_htg_grid_is_cell_masked

   function ftn_htg_grid_get_pure_mask_raw_bytes(grid, ptr_out, nbits_out, nbytes_out) &
         result(code) bind(C, name="ftn_htg_grid_get_pure_mask_raw_bytes")
      type(c_ptr), value :: grid
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int64_t), intent(out) :: nbits_out, nbytes_out
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_bit_array_t), pointer :: pmptr
      integer(byte_k), pointer :: raw(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_pure_mask(pmptr, st)
      code = int(st%code, c_int32_t)
      if (.not. st%ok() .or. .not. associated(pmptr)) then
         ptr_out = c_null_ptr
         nbits_out = 0_c_int64_t
         nbytes_out = 0_c_int64_t
         return
      end if
      nbits_out = int(pmptr%get_number_of_values(), c_int64_t)
      call pmptr%get_raw_bytes(raw)
      if (associated(raw) .and. size(raw) > 0) then
         ptr_out = c_loc(raw(0_id_k))
         nbytes_out = int(size(raw, kind=id_k), c_int64_t)
      else
         ptr_out = c_null_ptr
         nbytes_out = 0_c_int64_t
      end if
   end function ftn_htg_grid_get_pure_mask_raw_bytes

   function ftn_htg_grid_get_pure_mask_bit(grid, global_id, status_out) &
         result(flag) bind(C, name="ftn_htg_grid_get_pure_mask_bit")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: global_id
      type(c_ptr), value :: status_out
      logical(c_bool) :: flag
      type(htg_grid_t), pointer :: gptr
      type(htg_bit_array_t), pointer :: pmptr
      logical :: fval
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_pure_mask(pmptr, st)
      if (.not. st%ok() .or. .not. associated(pmptr)) then
         call write_opt_i32(status_out, int(st%code, c_int32_t))
         flag = logical(.false., kind=c_bool)
         return
      end if
      fval = pmptr%get_value(int(global_id, id_k), st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      flag = logical(fval, kind=c_bool)
   end function ftn_htg_grid_get_pure_mask_bit

   subroutine ftn_htg_grid_set_has_interface(grid, flag) &
         bind(C, name="ftn_htg_grid_set_has_interface")
      type(c_ptr), value :: grid
      logical(c_bool), value :: flag
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      call gptr%set_has_interface(logical(flag))
   end subroutine ftn_htg_grid_set_has_interface

   function ftn_htg_grid_get_has_interface(grid) &
         result(flag) bind(C, name="ftn_htg_grid_get_has_interface")
      type(c_ptr), value :: grid
      logical(c_bool) :: flag
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      flag = logical(gptr%get_has_interface(), kind=c_bool)
   end function ftn_htg_grid_get_has_interface

   subroutine ftn_htg_grid_set_interface_normals_name(grid, name, name_len) &
         bind(C, name="ftn_htg_grid_set_interface_normals_name")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      call gptr%set_interface_normals_name(c_chars_to_f_str(name, name_len))
   end subroutine ftn_htg_grid_set_interface_normals_name

   subroutine ftn_htg_grid_set_interface_intercepts_name(grid, name, name_len) &
         bind(C, name="ftn_htg_grid_set_interface_intercepts_name")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      type(htg_grid_t), pointer :: gptr

      call c_f_pointer(grid, gptr)
      call gptr%set_interface_intercepts_name(c_chars_to_f_str(name, name_len))
   end subroutine ftn_htg_grid_set_interface_intercepts_name

   function ftn_htg_grid_cell_data_add_array( &
         grid, name, name_len, data_type, num_components, num_tuples) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_add_array")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int32_t), value :: data_type
      integer(c_int32_t), value :: num_components
      integer(c_int64_t), value :: num_tuples
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%add_array(c_chars_to_f_str(name, name_len), &
         int(data_type, i32_k), int(num_components, i32_k), int(num_tuples, id_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_add_array

   function ftn_htg_grid_cell_data_ensure_tuples(grid, min_tuples) &
         result(code) bind(C, name="ftn_htg_grid_cell_data_ensure_tuples")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: min_tuples
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%ensure_tuples(int(min_tuples, id_k), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_ensure_tuples

   function ftn_htg_grid_cell_data_set_tuple_r64( &
         grid, name, name_len, tuple_id, vals_ptr, ncomp) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_set_tuple_r64")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      type(c_ptr), value :: vals_ptr
      integer(c_int32_t), value :: ncomp
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      real(r64_k), pointer :: fvals(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      call c_f_pointer(vals_ptr, fvals, [ncomp])
      call fptr%set_tuple_r64(int(tuple_id, id_k), fvals(1:ncomp), st)
      call gptr%invalidate_pure_mask()
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_set_tuple_r64

   function ftn_htg_grid_cell_data_get_tuple_r64( &
         grid, name, name_len, tuple_id, vals_out_ptr, ncomp) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_get_tuple_r64")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      type(c_ptr), value :: vals_out_ptr
      integer(c_int32_t), value :: ncomp
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      real(r64_k), pointer :: fvals(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      call c_f_pointer(vals_out_ptr, fvals, [ncomp])
      call fptr%get_tuple_r64(int(tuple_id, id_k), fvals(1:ncomp), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_get_tuple_r64

   function ftn_htg_grid_cell_data_set_tuple_r32( &
         grid, name, name_len, tuple_id, vals_ptr, ncomp) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_set_tuple_r32")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      type(c_ptr), value :: vals_ptr
      integer(c_int32_t), value :: ncomp
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      real(r32_k), pointer :: fvals(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      call c_f_pointer(vals_ptr, fvals, [ncomp])
      call fptr%set_tuple_r32(int(tuple_id, id_k), fvals(1:ncomp), st)
      call gptr%invalidate_pure_mask()
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_set_tuple_r32

   function ftn_htg_grid_cell_data_get_tuple_r32( &
         grid, name, name_len, tuple_id, vals_out_ptr, ncomp) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_get_tuple_r32")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      type(c_ptr), value :: vals_out_ptr
      integer(c_int32_t), value :: ncomp
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      real(r32_k), pointer :: fvals(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      call c_f_pointer(vals_out_ptr, fvals, [ncomp])
      call fptr%get_tuple_r32(int(tuple_id, id_k), fvals(1:ncomp), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_get_tuple_r32

   function ftn_htg_grid_cell_data_set_tuple_i64( &
         grid, name, name_len, tuple_id, vals_ptr, ncomp) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_set_tuple_i64")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      type(c_ptr), value :: vals_ptr
      integer(c_int32_t), value :: ncomp
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      integer(i64_k), pointer :: fvals(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      call c_f_pointer(vals_ptr, fvals, [ncomp])
      call fptr%set_tuple_i64(int(tuple_id, id_k), fvals(1:ncomp), st)
      call gptr%invalidate_pure_mask()
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_set_tuple_i64

   function ftn_htg_grid_cell_data_get_tuple_i64( &
         grid, name, name_len, tuple_id, vals_out_ptr, ncomp) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_get_tuple_i64")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      type(c_ptr), value :: vals_out_ptr
      integer(c_int32_t), value :: ncomp
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      integer(i64_k), pointer :: fvals(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      call c_f_pointer(vals_out_ptr, fvals, [ncomp])
      call fptr%get_tuple_i64(int(tuple_id, id_k), fvals(1:ncomp), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_get_tuple_i64

   function ftn_htg_grid_cell_data_set_tuple_i32( &
         grid, name, name_len, tuple_id, vals_ptr, ncomp) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_set_tuple_i32")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      type(c_ptr), value :: vals_ptr
      integer(c_int32_t), value :: ncomp
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      integer(i32_k), pointer :: fvals(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      call c_f_pointer(vals_ptr, fvals, [ncomp])
      call fptr%set_tuple_i32(int(tuple_id, id_k), fvals(1:ncomp), st)
      call gptr%invalidate_pure_mask()
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_set_tuple_i32

   function ftn_htg_grid_cell_data_get_tuple_i32( &
         grid, name, name_len, tuple_id, vals_out_ptr, ncomp) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_get_tuple_i32")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      type(c_ptr), value :: vals_out_ptr
      integer(c_int32_t), value :: ncomp
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      integer(i32_k), pointer :: fvals(:)
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      call c_f_pointer(vals_out_ptr, fvals, [ncomp])
      call fptr%get_tuple_i32(int(tuple_id, id_k), fvals(1:ncomp), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_get_tuple_i32

   function ftn_htg_grid_cell_data_set_component_bit( &
         grid, name, name_len, tuple_id, comp_id, val) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_set_component_bit")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      integer(c_int32_t), value :: comp_id
      logical(c_bool), value :: val
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      call fptr%set_component_bit(int(tuple_id, id_k), int(comp_id, i32_k), logical(val), st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_set_component_bit

   function ftn_htg_grid_cell_data_get_component_bit( &
         grid, name, name_len, tuple_id, comp_id, status_out) result(flag) &
         bind(C, name="ftn_htg_grid_cell_data_get_component_bit")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      integer(c_int64_t), value :: tuple_id
      integer(c_int32_t), value :: comp_id
      type(c_ptr), value :: status_out
      logical(c_bool) :: flag
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      logical :: fval
      type(htg_status_t) :: st

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         call write_opt_i32(status_out, int(st%code, c_int32_t))
         flag = logical(.false., kind=c_bool)
         return
      end if
      fval = fptr%get_component_bit(int(tuple_id, id_k), int(comp_id, i32_k), st)
      call write_opt_i32(status_out, int(st%code, c_int32_t))
      flag = logical(fval, kind=c_bool)
   end function ftn_htg_grid_cell_data_get_component_bit

   function ftn_htg_grid_cell_data_get_view_r64( &
         grid, name, name_len, ptr_out, ncomp_out, ntuples_out) result(code) &
         bind(C, name="ftn_htg_grid_cell_data_get_view_r64")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t), value :: name_len
      type(c_ptr), intent(out) :: ptr_out
      integer(c_int32_t), intent(out) :: ncomp_out
      integer(c_int64_t), intent(out) :: ntuples_out
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_cell_data_t), pointer :: cdptr
      type(htg_field_array_t), pointer :: fptr
      real(r64_k), pointer :: view(:, :)
      type(htg_status_t) :: st

      ptr_out = c_null_ptr
      ncomp_out = 0_c_int32_t
      ntuples_out = 0_c_int64_t

      call c_f_pointer(grid, gptr)
      call gptr%get_cell_data(cdptr)
      call cdptr%get_array(c_chars_to_f_str(name, name_len), fptr, st)
      if (.not. st%ok() .or. .not. associated(fptr)) then
         code = int(st%code, c_int32_t)
         return
      end if
      ncomp_out = int(fptr%get_number_of_components(), c_int32_t)
      ntuples_out = int(fptr%get_number_of_tuples(), c_int64_t)
      call fptr%get_view_r64(view, st)
      if (st%ok() .and. associated(view) .and. size(view) > 0) then
         ptr_out = c_loc(view(0_i32_k, 0_id_k))
      end if
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_cell_data_get_view_r64

   function ftn_htg_grid_validate(grid, err_msg_out, err_msg_cap) &
         result(code) bind(C, name="ftn_htg_grid_validate")
      type(c_ptr), value :: grid
      character(kind=c_char), intent(out) :: err_msg_out(*)
      integer(c_int32_t), value :: err_msg_cap
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st
      integer :: i, mlen, copy_len

      call c_f_pointer(grid, gptr)
      call gptr%validate(st)
      code = int(st%code, c_int32_t)

      if (err_msg_cap > 0_c_int32_t) then
         if (.not. st%ok()) then
            mlen = len_trim(st%message)
            copy_len = min(mlen, int(err_msg_cap) - 1)
            do i = 1, copy_len
               err_msg_out(i) = st%message(i:i)
            end do
            err_msg_out(copy_len + 1) = c_null_char
         else
            err_msg_out(1) = c_null_char
         end if
      end if
   end function ftn_htg_grid_validate

   ! =========================================================================
   ! Cell geometry, face neighbours & coarsening
   ! =========================================================================

   function ftn_htg_grid_get_level_zero_origin_and_size(grid, tree_id, origin_out, size_out) &
         result(code) bind(C, name="ftn_htg_grid_get_level_zero_origin_and_size")
      type(c_ptr), value :: grid
      integer(c_int64_t), value :: tree_id
      real(c_double), intent(out) :: origin_out(3)
      real(c_double), intent(out) :: size_out(3)
      integer(c_int32_t) :: code
      type(htg_grid_t), pointer :: gptr
      type(htg_status_t) :: st
      real(r64_k) :: origin(0:2), cell_size(0:2)

      call c_f_pointer(grid, gptr)
      call gptr%get_level_zero_origin_and_size_from_index(int(tree_id, id_k), origin, cell_size, st)
      origin_out = origin
      size_out = cell_size
      code = int(st%code, c_int32_t)
   end function ftn_htg_grid_get_level_zero_origin_and_size

   function ftn_htg_cursor_get_origin(cursor, origin_out) &
         result(code) bind(C, name="ftn_htg_cursor_get_origin")
      type(c_ptr), value :: cursor
      real(c_double), intent(out) :: origin_out(3)
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st
      real(r64_k) :: origin(0:2)

      call c_f_pointer(cursor, cptr)
      call cptr%get_origin(origin, st)
      origin_out = origin
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_get_origin

   function ftn_htg_cursor_get_size(cursor, size_out) &
         result(code) bind(C, name="ftn_htg_cursor_get_size")
      type(c_ptr), value :: cursor
      real(c_double), intent(out) :: size_out(3)
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st
      real(r64_k) :: cell_size(0:2)

      call c_f_pointer(cursor, cptr)
      call cptr%get_size(cell_size, st)
      size_out = cell_size
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_get_size

   function ftn_htg_cursor_get_bounds(cursor, bounds_out) &
         result(code) bind(C, name="ftn_htg_cursor_get_bounds")
      type(c_ptr), value :: cursor
      real(c_double), intent(out) :: bounds_out(6)
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st
      real(r64_k) :: bounds(0:5)

      call c_f_pointer(cursor, cptr)
      call cptr%get_bounds(bounds, st)
      bounds_out = bounds
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_get_bounds

   function ftn_htg_cursor_get_point(cursor, point_out) &
         result(code) bind(C, name="ftn_htg_cursor_get_point")
      type(c_ptr), value :: cursor
      real(c_double), intent(out) :: point_out(3)
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st
      real(r64_k) :: point(0:2)

      call c_f_pointer(cursor, cptr)
      call cptr%get_point(point, st)
      point_out = point
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_get_point

   function ftn_htg_cursor_get_child_index(cursor) &
         result(res) bind(C, name="ftn_htg_cursor_get_child_index")
      type(c_ptr), value :: cursor
      integer(c_int32_t) :: res
      type(htg_cursor_t), pointer :: cptr

      call c_f_pointer(cursor, cptr)
      res = int(cptr%get_child_index(), c_int32_t)
   end function ftn_htg_cursor_get_child_index

   !> Face neighbour of the cursor's cell along physical `axis` in `direction`.
   !> When there is one, reports its tree, local and global id (-1 if unset),
   !> level, leaf and mask state and bounds.
   function ftn_htg_cursor_get_face_neighbor(cursor, axis, direction, has_tree_out, &
         tree_id_out, vertex_id_out, global_id_out, level_out, is_leaf_out, is_masked_out, &
         bounds_out) result(code) bind(C, name="ftn_htg_cursor_get_face_neighbor")
      type(c_ptr), value :: cursor
      integer(c_int32_t), value :: axis
      integer(c_int32_t), value :: direction
      logical(c_bool), intent(out) :: has_tree_out
      integer(c_int64_t), intent(out) :: tree_id_out
      integer(c_int64_t), intent(out) :: vertex_id_out
      integer(c_int64_t), intent(out) :: global_id_out
      integer(c_int32_t), intent(out) :: level_out
      logical(c_bool), intent(out) :: is_leaf_out
      logical(c_bool), intent(out) :: is_masked_out
      real(c_double), intent(out) :: bounds_out(6)
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_cursor_t) :: nb
      type(htg_status_t) :: st, st2
      real(r64_k) :: bounds(0:5)
      integer(id_k) :: gid
      logical :: flag

      call c_f_pointer(cursor, cptr)
      has_tree_out = .false._c_bool
      tree_id_out = int(HTG_INVALID_ID, c_int64_t)
      vertex_id_out = int(HTG_INVALID_ID, c_int64_t)
      global_id_out = int(HTG_INVALID_ID, c_int64_t)
      level_out = -1_c_int32_t
      is_leaf_out = .false._c_bool
      is_masked_out = .false._c_bool
      bounds_out = 0.0_c_double

      call cptr%get_face_neighbor(int(axis, i32_k), int(direction, i32_k), nb, st)
      code = int(st%code, c_int32_t)
      if (.not. st%ok()) return

      tree_id_out = int(nb%get_tree_id(), c_int64_t)
      if (.not. nb%has_tree()) return
      has_tree_out = .true._c_bool
      vertex_id_out = int(nb%get_vertex_id(), c_int64_t)
      level_out = int(nb%get_level(), c_int32_t)
      gid = nb%get_global_node_index(st2)
      if (st2%ok()) global_id_out = int(gid, c_int64_t)
      flag = nb%is_leaf(st2)
      is_leaf_out = logical(flag, kind=c_bool)
      flag = nb%is_masked(st2)
      is_masked_out = logical(flag, kind=c_bool)
      call nb%get_bounds(bounds, st2)
      bounds_out = bounds
   end function ftn_htg_cursor_get_face_neighbor

   function ftn_htg_cursor_coarsen_node(cursor) &
         result(code) bind(C, name="ftn_htg_cursor_coarsen_node")
      type(c_ptr), value :: cursor
      integer(c_int32_t) :: code
      type(htg_cursor_t), pointer :: cptr
      type(htg_status_t) :: st

      call c_f_pointer(cursor, cptr)
      call cptr%coarsen_node(st)
      code = int(st%code, c_int32_t)
   end function ftn_htg_cursor_coarsen_node

end module fvtkhtg_test_c_bridge_m
