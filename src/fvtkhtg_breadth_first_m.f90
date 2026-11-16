!> Breadth-first export of a forest: the seam between fvtkhtg (construction,
!> refinement, validation) and a serializer such as the fvtkhdf
!> `vtkhdf_htg_file` writer.
!>
!> INDEX CONVENTION AT THE SEAM
!> Everything inside fvtkhtg is 0-based, matching VTK. `htg_breadth_first_t`
!> is the one deliberate exception: its arrays follow the fvtkhdf writer's
!> Fortran conventions, so they can be passed to
!> `write_mesh_topology(tree_ids, xtree, refined [, mask])` unchanged.
!>   - arrays have the default lower bound 1;
!>   - `tree_ids` holds 1-based root-cell numbers (fvtkhtg tree_id + 1);
!>   - `xtree` holds 1-based starting positions of each tree's cells;
!>   - `global_ids` holds fvtkhtg's own 0-based global IDs (CellData tuple
!>     indices) and is only ever used to gather cell data.
module fvtkhtg_breadth_first_m
   use fvtkhtg_kinds_m, only: byte_k, i32_k, i64_k, id_k, r32_k, r64_k, &
      htg_status_t, HTG_ERR_INVALID_ARG, HTG_ERR_INVALID_STATE, &
      HTG_ERR_OUT_OF_BOUNDS, HTG_TYPE_BIT, HTG_TYPE_INT32, HTG_TYPE_INT64, &
      HTG_TYPE_REAL32, HTG_TYPE_REAL64, htg_fail, htg_set_ok
   use fvtkhtg_cell_data_m, only: htg_field_array_t
   implicit none
   private

   !> Per-cell breadth-first description of the trees of one grid (one MPI
   !> rank's partition). Cells are listed tree by tree in ascending tree_id
   !> order; within a tree, level by level, with the children of each refined
   !> cell contiguous, in child order, and in the order of their parents.
   type, public :: htg_breadth_first_t
      integer(id_k), allocatable :: tree_ids(:)    !< (ntree)   1-based root-cell numbers
      integer(id_k), allocatable :: xtree(:)       !< (ntree+1) 1-based start of each tree's cells
      logical, allocatable :: refined(:)           !< (ncell)   .true. for a refined cell
      logical, allocatable :: mask(:)              !< (ncell)   allocated only when a mask is exported
      integer(id_k), allocatable :: global_ids(:)  !< (ncell)   0-based fvtkhtg global ID of each cell
      !> .true. when global_ids(n) == n - 1 for every cell, i.e. CellData is
      !> already stored in export order and needs no permutation.
      logical :: identity_order = .false.
   contains
      procedure :: number_of_trees => bf_number_of_trees
      procedure :: number_of_cells => bf_number_of_cells
      procedure :: clear => bf_clear
      procedure, private :: bf_gather_r64_1
      procedure, private :: bf_gather_r64_2
      procedure, private :: bf_gather_r32_1
      procedure, private :: bf_gather_r32_2
      procedure, private :: bf_gather_i64_1
      procedure, private :: bf_gather_i64_2
      procedure, private :: bf_gather_i32_1
      procedure, private :: bf_gather_i32_2
      procedure, private :: bf_gather_bit_1
      procedure, private :: bf_gather_bit_2
      !> Gather a CellData array into export order: rank-1 `values(ncell)`
      !> for single-component arrays, rank-2 `values(ncomp, ncell)` for any
      !> array. BIT arrays gather into `integer(byte_k)` 0/1 values, the
      !> fvtkhdf writer having no bit type.
      generic :: gather => bf_gather_r64_1, bf_gather_r64_2, &
         bf_gather_r32_1, bf_gather_r32_2, &
         bf_gather_i64_1, bf_gather_i64_2, &
         bf_gather_i32_1, bf_gather_i32_2, &
         bf_gather_bit_1, bf_gather_bit_2
   end type htg_breadth_first_t

contains

   pure function bf_number_of_trees(self) result(n)
      class(htg_breadth_first_t), intent(in) :: self
      integer(id_k) :: n
      n = 0_id_k
      if (allocated(self%tree_ids)) n = size(self%tree_ids, kind=id_k)
   end function bf_number_of_trees

   pure function bf_number_of_cells(self) result(n)
      class(htg_breadth_first_t), intent(in) :: self
      integer(id_k) :: n
      n = 0_id_k
      if (allocated(self%refined)) n = size(self%refined, kind=id_k)
   end function bf_number_of_cells

   pure subroutine bf_clear(self)
      class(htg_breadth_first_t), intent(inout) :: self
      if (allocated(self%tree_ids)) deallocate(self%tree_ids)
      if (allocated(self%xtree)) deallocate(self%xtree)
      if (allocated(self%refined)) deallocate(self%refined)
      if (allocated(self%mask)) deallocate(self%mask)
      if (allocated(self%global_ids)) deallocate(self%global_ids)
      self%identity_order = .false.
   end subroutine bf_clear

   !> Common preconditions of every gather. On failure the caller's output
   !> is left unallocated.
   logical function gather_ok(self, field, dtype, scalar_only, status) result(ok)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      integer(i32_k), intent(in) :: dtype
      logical, intent(in) :: scalar_only
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: ncell

      ok = .false.
      if (.not. field%is_initialized()) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_breadth_first_t%gather: field array is not initialized', status)
         return
      end if
      if (field%get_data_type() /= dtype) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_breadth_first_t%gather: output kind does not match field "' // &
            field%get_name() // '" data type', status)
         return
      end if
      if (scalar_only .and. field%get_number_of_components() /= 1_i32_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_breadth_first_t%gather: rank-1 output requires a single-component field; "' // &
            field%get_name() // '" has more', status)
         return
      end if
      ncell = self%number_of_cells()
      if (ncell > 0_id_k) then
         if (.not. allocated(self%global_ids)) then
            call htg_fail(HTG_ERR_INVALID_STATE, &
               'htg_breadth_first_t%gather: export has no global_ids', status)
            return
         end if
         if (minval(self%global_ids) < 0_id_k .or. &
             maxval(self%global_ids) >= field%get_number_of_tuples()) then
            call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
               'htg_breadth_first_t%gather: field "' // field%get_name() // &
               '" has fewer tuples than the exported global IDs require', status)
            return
         end if
      end if
      ok = .true.
   end function gather_ok

   subroutine bf_gather_r64_1(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      real(r64_k), allocatable, intent(out) :: values(:)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n

      if (.not. gather_ok(self, field, HTG_TYPE_REAL64, .true., status)) return
      allocate(values(self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         values(n) = field%r64_values(0, self%global_ids(n))
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_r64_1

   subroutine bf_gather_r64_2(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      real(r64_k), allocatable, intent(out) :: values(:, :)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n

      if (.not. gather_ok(self, field, HTG_TYPE_REAL64, .false., status)) return
      allocate(values(field%get_number_of_components(), self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         values(:, n) = field%r64_values(:, self%global_ids(n))
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_r64_2

   subroutine bf_gather_r32_1(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      real(r32_k), allocatable, intent(out) :: values(:)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n

      if (.not. gather_ok(self, field, HTG_TYPE_REAL32, .true., status)) return
      allocate(values(self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         values(n) = field%r32_values(0, self%global_ids(n))
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_r32_1

   subroutine bf_gather_r32_2(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      real(r32_k), allocatable, intent(out) :: values(:, :)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n

      if (.not. gather_ok(self, field, HTG_TYPE_REAL32, .false., status)) return
      allocate(values(field%get_number_of_components(), self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         values(:, n) = field%r32_values(:, self%global_ids(n))
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_r32_2

   subroutine bf_gather_i64_1(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      integer(i64_k), allocatable, intent(out) :: values(:)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n

      if (.not. gather_ok(self, field, HTG_TYPE_INT64, .true., status)) return
      allocate(values(self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         values(n) = field%i64_values(0, self%global_ids(n))
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_i64_1

   subroutine bf_gather_i64_2(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      integer(i64_k), allocatable, intent(out) :: values(:, :)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n

      if (.not. gather_ok(self, field, HTG_TYPE_INT64, .false., status)) return
      allocate(values(field%get_number_of_components(), self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         values(:, n) = field%i64_values(:, self%global_ids(n))
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_i64_2

   subroutine bf_gather_i32_1(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      integer(i32_k), allocatable, intent(out) :: values(:)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n

      if (.not. gather_ok(self, field, HTG_TYPE_INT32, .true., status)) return
      allocate(values(self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         values(n) = field%i32_values(0, self%global_ids(n))
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_i32_1

   subroutine bf_gather_i32_2(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      integer(i32_k), allocatable, intent(out) :: values(:, :)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n

      if (.not. gather_ok(self, field, HTG_TYPE_INT32, .false., status)) return
      allocate(values(field%get_number_of_components(), self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         values(:, n) = field%i32_values(:, self%global_ids(n))
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_i32_2

   subroutine bf_gather_bit_1(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      integer(byte_k), allocatable, intent(out) :: values(:)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n
      logical :: bit

      if (.not. gather_ok(self, field, HTG_TYPE_BIT, .true., status)) return
      allocate(values(self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         bit = field%bit_values%get_value(self%global_ids(n))
         values(n) = merge(1_byte_k, 0_byte_k, bit)
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_bit_1

   subroutine bf_gather_bit_2(self, field, values, status)
      class(htg_breadth_first_t), intent(in) :: self
      type(htg_field_array_t), intent(in) :: field
      integer(byte_k), allocatable, intent(out) :: values(:, :)
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: n, nc, c
      logical :: bit

      nc = int(field%get_number_of_components(), id_k)
      if (.not. gather_ok(self, field, HTG_TYPE_BIT, .false., status)) return
      allocate(values(nc, self%number_of_cells()))
      do n = 1_id_k, self%number_of_cells()
         do c = 0_id_k, nc - 1_id_k
            bit = field%bit_values%get_value(self%global_ids(n) * nc + c)
            values(c + 1_id_k, n) = merge(1_byte_k, 0_byte_k, bit)
         end do
      end do
      call htg_set_ok(status)
   end subroutine bf_gather_bit_2

end module fvtkhtg_breadth_first_m
