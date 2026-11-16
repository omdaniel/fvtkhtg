!> Example 06: Adapting a Mesh with Cell Geometry, Face Neighbours & Coarsening
!>
!> Demonstrates the three operations a solver needs to adapt a mesh in place:
!>   - `cursor%get_bounds` / `get_point`: refine the cells that a circle
!>     passes through, decided from each cell's own geometry
!>   - `cursor%get_face_neighbor`: enforce 2:1 balance (no leaf more than one
!>     level finer than a face neighbour), refining THROUGH the neighbour
!>     cursor, across root-cell boundaries as well as inside a tree
!>   - `cursor%coarsen_node`: when the circle moves, merge the cells it has
!>     left; CellData follows the cells that are renumbered
!>
!> Global ids are implicit with a fixed range per tree (`tree_id * CAPACITY`),
!> so trees can grow and shrink without ever colliding.
program example_06_adaptive_refinement
   use fvtkhtg, only: &
      htg_grid_t, htg_cursor_t, htg_cell_data_t, htg_field_array_t, &
      htg_breadth_first_t, htg_status_t, new_non_oriented_cursor, &
      HTG_TYPE_REAL64, HTG_AXIS_X, HTG_AXIS_Y, id_k, i32_k, r64_k
   implicit none

   integer(i32_k), parameter :: NX = 4_i32_k, NY = 2_i32_k   ! root cells
   integer(i32_k), parameter :: MAX_LEVEL = 4_i32_k
   integer(id_k), parameter :: CAPACITY = 512_id_k           ! >= cells of a full tree (341)
   real(r64_k), parameter :: RADIUS = 0.6_r64_k

   type(htg_grid_t), target :: grid
   type(htg_cursor_t) :: cursor
   type(htg_cell_data_t), pointer :: cell_data => null()
   type(htg_field_array_t), pointer :: distance => null()
   type(htg_breadth_first_t) :: bf
   type(htg_status_t) :: status
   real(r64_k) :: centre(0:1)
   real(r64_k), allocatable :: distance_bf(:)
   integer(id_k) :: tree_id, ntrees
   integer(i32_k) :: i
   integer :: nrefined, nbalanced, ncoarsened, step

   ! 1. Root grid and one unrefined tree per root cell
   call grid%initialize()
   call grid%set_dimensions(NX + 1_i32_k, NY + 1_i32_k, 1_i32_k, status)
   call grid%set_x_coordinates([(real(i, r64_k), i = 0_i32_k, NX)], status)
   call grid%set_y_coordinates([(real(i, r64_k), i = 0_i32_k, NY)], status)
   call check_ok('grid setup')
   ntrees = grid%get_max_number_of_trees()
   do tree_id = 0_id_k, ntrees - 1_id_k
      call new_non_oriented_cursor(grid, tree_id, cursor, create=.true., status=status)
      call cursor%set_global_index_start(tree_id * CAPACITY, status)
   end do
   call grid%get_cell_data(cell_data)
   call cell_data%add_array('Distance', HTG_TYPE_REAL64, num_components=1_i32_k, &
      num_tuples=ntrees * CAPACITY, status=status)
   call cell_data%get_array('Distance', distance, status)
   call check_ok('cell data setup')

   ! 2. Follow the circle as it moves from left to right
   do step = 0, 2
      centre = [0.9_r64_k + 1.1_r64_k * real(step, r64_k), 1.0_r64_k]

      ! Merge cells the circle has left (their parent keeps its own value)
      ncoarsened = 0
      do tree_id = 0_id_k, ntrees - 1_id_k
         call new_non_oriented_cursor(grid, tree_id, cursor, create=.false., status=status)
         call coarsen_where_unneeded(cursor, ncoarsened)
      end do

      ! Refine the cells the circle now passes through
      nrefined = 0
      do tree_id = 0_id_k, ntrees - 1_id_k
         call new_non_oriented_cursor(grid, tree_id, cursor, create=.false., status=status)
         call refine_where_needed(cursor, nrefined)
      end do

      ! Restore 2:1 balance, sweeping until nothing changes
      nbalanced = 0
      do
         i = 0_i32_k
         do tree_id = 0_id_k, ntrees - 1_id_k
            call new_non_oriented_cursor(grid, tree_id, cursor, create=.false., status=status)
            call balance(cursor, i)
         end do
         if (i == 0_i32_k) exit
         nbalanced = nbalanced + int(i)
      end do

      ! Sample the field on every cell, then hand the grid to a serializer
      do tree_id = 0_id_k, ntrees - 1_id_k
         call new_non_oriented_cursor(grid, tree_id, cursor, create=.false., status=status)
         call sample_field(cursor)
      end do
      call grid%validate(status)
      call check_ok('validate')
      if (.not. is_balanced()) error stop 'grid is not 2:1 balanced'
      call grid%export_breadth_first(bf, status)
      call bf%gather(distance, distance_bf, status)
      call check_ok('export')

      print '(a,i0,a,f4.1,a)', 'Step ', step, ': circle centre x = ', centre(0), ''
      print '(a,i0,a,i0,a,i0)', '   coarsened ', ncoarsened, ' cells, refined ', nrefined, &
         ', refined for balance ', nbalanced
      print '(a,i0,a,i0,a,i0,a)', '   ', grid%get_number_of_cells(), ' cells, ', &
         grid%get_number_of_leaves(), ' leaves, ', grid%get_number_of_levels(), ' levels; 2:1 balanced'
      print '(a,i0,a,f7.4)', '   exported ', bf%number_of_cells(), &
         ' cells; smallest |distance to circle| = ', minval(abs(distance_bf))
   end do
   print '(a)', 'Example 06 completed successfully.'

contains

   subroutine check_ok(what)
      character(len=*), intent(in) :: what
      if (.not. status%ok()) error stop what // ': ' // trim(status%message)
   end subroutine check_ok

   !> Does the circle pass through the cell under the cursor? Decided from
   !> the cell's bounds: the circle crosses a box when the box's nearest point
   !> is inside the circle and its farthest point is outside.
   logical function circle_crosses(c) result(crosses)
      type(htg_cursor_t), intent(in) :: c
      real(r64_k) :: b(0:5), nearest, farthest, lo, hi, d
      type(htg_status_t) :: st
      integer :: a

      call c%get_bounds(b, st)
      nearest = 0.0_r64_k
      farthest = 0.0_r64_k
      do a = 0, 1
         lo = b(2*a) - centre(a)
         hi = b(2*a + 1) - centre(a)
         d = max(lo, -hi, 0.0_r64_k)            ! distance from the centre to the box along this axis
         nearest = nearest + d*d
         farthest = farthest + max(lo*lo, hi*hi)
      end do
      crosses = (nearest < RADIUS*RADIUS) .and. (farthest > RADIUS*RADIUS)
   end function circle_crosses

   recursive subroutine refine_where_needed(c, count)
      type(htg_cursor_t), intent(inout) :: c
      integer, intent(inout) :: count
      type(htg_status_t) :: st
      integer(i32_k) :: ic
      logical :: leaf

      leaf = c%is_leaf(st)
      if (leaf) then
         if (c%get_level() >= MAX_LEVEL) return
         if (.not. circle_crosses(c)) return
         call c%subdivide_leaf(st)
         count = count + 1
      end if
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, st)
         call refine_where_needed(c, count)
         call c%to_parent(st)
      end do
   end subroutine refine_where_needed

   !> Bottom-up: a refined cell whose children are all leaves is merged back
   !> into one leaf when the circle no longer passes through it.
   recursive subroutine coarsen_where_unneeded(c, count)
      type(htg_cursor_t), intent(inout) :: c
      integer, intent(inout) :: count
      type(htg_status_t) :: st
      integer(i32_k) :: ic
      logical :: leaf, terminal, needed

      leaf = c%is_leaf(st)
      if (leaf) return
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, st)
         call coarsen_where_unneeded(c, count)
         call c%to_parent(st)
      end do
      terminal = c%is_terminal_node(st)
      needed = circle_crosses(c)
      if (terminal .and. .not. needed) then
         call c%coarsen_node(st)       ! the cursor stays on this cell, now a leaf
         if (.not. st%ok()) error stop 'coarsen_node: ' // trim(st%message)
         count = count + 1
      end if
   end subroutine coarsen_where_unneeded

   !> Refine any leaf that is two or more levels coarser than a leaf next to
   !> it. Looked at from the finer side: the face neighbour of a leaf is a
   !> leaf of a level at least 2 lower.
   recursive subroutine balance(c, count)
      type(htg_cursor_t), intent(inout) :: c
      integer(i32_k), intent(inout) :: count
      type(htg_cursor_t) :: neighbor
      type(htg_status_t) :: st
      integer(i32_k) :: ic, axis, direction
      logical :: leaf, neighbor_is_leaf

      leaf = c%is_leaf(st)
      if (.not. leaf) then
         do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
            call c%to_child(ic, st)
            call balance(c, count)
            call c%to_parent(st)
         end do
         return
      end if
      do axis = HTG_AXIS_X, HTG_AXIS_Y
         do direction = -1_i32_k, 1_i32_k, 2_i32_k
            call c%get_face_neighbor(axis, direction, neighbor, st)
            if (.not. neighbor%has_tree()) cycle            ! edge of the grid
            neighbor_is_leaf = neighbor%is_leaf(st)
            if (neighbor_is_leaf .and. neighbor%get_level() < c%get_level() - 1_i32_k) then
               call neighbor%subdivide_leaf(st)             ! refine through the neighbour cursor
               count = count + 1_i32_k
            end if
         end do
      end do
   end subroutine balance

   logical function is_balanced() result(ok)
      integer(id_k) :: t
      type(htg_cursor_t) :: c
      type(htg_status_t) :: st

      ok = .true.
      do t = 0_id_k, ntrees - 1_id_k
         call new_non_oriented_cursor(grid, t, c, create=.false., status=st)
         call check_balance(c, ok)
      end do
   end function is_balanced

   recursive subroutine check_balance(c, ok)
      type(htg_cursor_t), intent(inout) :: c
      logical, intent(inout) :: ok
      type(htg_cursor_t) :: neighbor
      type(htg_status_t) :: st
      integer(i32_k) :: ic, axis, direction
      logical :: leaf, neighbor_is_leaf

      leaf = c%is_leaf(st)
      if (.not. leaf) then
         do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
            call c%to_child(ic, st)
            call check_balance(c, ok)
            call c%to_parent(st)
         end do
         return
      end if
      do axis = HTG_AXIS_X, HTG_AXIS_Y
         do direction = -1_i32_k, 1_i32_k, 2_i32_k
            call c%get_face_neighbor(axis, direction, neighbor, st)
            if (.not. neighbor%has_tree()) cycle
            neighbor_is_leaf = neighbor%is_leaf(st)
            if (neighbor_is_leaf .and. neighbor%get_level() < c%get_level() - 1_i32_k) ok = .false.
         end do
      end do
   end subroutine check_balance

   !> Signed distance from each cell's centre to the circle
   recursive subroutine sample_field(c)
      type(htg_cursor_t), intent(inout) :: c
      type(htg_status_t) :: st
      real(r64_k) :: p(0:2)
      integer(i32_k) :: ic
      logical :: leaf

      call c%get_point(p, st)
      call distance%set_scalar_r64(c%get_global_node_index(), &
         norm2(p(0:1) - centre) - RADIUS, st)
      leaf = c%is_leaf(st)
      if (leaf) return
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, st)
         call sample_field(c)
         call c%to_parent(st)
      end do
   end subroutine sample_field

end program example_06_adaptive_refinement
