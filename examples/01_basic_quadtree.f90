!> Example 01: Building and Validating a 2D Rectilinear Quadtree Grid
!>
!> Demonstrates:
!>   - Declaring `type(htg_grid_t), target :: grid` for compile-time cursor safety
!>   - Configuring 2D rectilinear XY geometry (`set_branch_factor`, `set_dimensions`,
!>     and 0-based monotonic `set_x/y/z_coordinates`)
!>   - Creating a build-time cursor via `new_non_oriented_cursor`
!>   - Recursive and selective AMR refinement (`subdivide_leaf`, `to_child`, `to_parent`)
!>   - Assigning contiguous implicit global IDs (`initialize_local_index_node`)
!>   - Querying forest-wide cell/leaf/level counts and running `grid%validate(status)`
program example_01_basic_quadtree
   use fvtkhtg, only: &
      htg_grid_t, htg_cursor_t, htg_status_t, &
      new_non_oriented_cursor, id_k, i32_k, r64_k
   implicit none

   type(htg_grid_t), target :: grid
   type(htg_cursor_t) :: cursor
   type(htg_status_t) :: status
   real(r64_k) :: x_coords(0:2), y_coords(0:2), z_coords(0:0)
   integer(i32_k) :: axis1, axis2

   ! 1. Initialize grid and configure a 2x2 forest of binary quadtrees in XY (dz = 1)
   call grid%initialize()
   call grid%set_branch_factor(2_i32_k, status)
   call check_ok(status, 'set_branch_factor')

   call grid%set_dimensions(3_i32_k, 3_i32_k, 1_i32_k, status)
   call check_ok(status, 'set_dimensions')

   x_coords = [0.0_r64_k, 0.5_r64_k, 1.0_r64_k]
   y_coords = [0.0_r64_k, 0.5_r64_k, 1.0_r64_k]
   z_coords = [0.0_r64_k]

   call grid%set_x_coordinates(x_coords, status)
   call check_ok(status, 'set_x_coordinates')
   call grid%set_y_coordinates(y_coords, status)
   call check_ok(status, 'set_y_coordinates')
   call grid%set_z_coordinates(z_coords, status)
   call check_ok(status, 'set_z_coordinates')

   call grid%get_2d_axes(axis1, axis2, status)
   call check_ok(status, 'get_2d_axes')

   ! 2. Tree 0: recursive diagonal refinement along child 3 up to depth 3 (4 levels)
   call new_non_oriented_cursor(grid, 0_id_k, cursor, create=.true., status=status)
   call check_ok(status, 'new_non_oriented_cursor(tree 0)')
   call refine_corner_recursively(cursor, 3_i32_k, 3_i32_k)
   if (.not. cursor%is_root()) error stop 'Expected cursor back at root of tree 0'

   ! 3. Tree 1: selective refinement on children 0 and 2 (leaving tree 2 empty/sparse)
   call new_non_oriented_cursor(grid, 1_id_k, cursor, create=.true., status=status)
   call check_ok(status, 'new_non_oriented_cursor(tree 1)')
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'tree 1 subdivide root')

   call cursor%to_child(0_i32_k, status)
   call check_ok(status, 'tree 1 to_child(0)')
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'tree 1 subdivide child 0')
   if (.not. cursor%is_terminal_node(status)) then
      error stop 'Expected refined child 0 to be a terminal coarse node'
   end if
   call cursor%to_parent(status)
   call check_ok(status, 'tree 1 to_parent')

   call cursor%to_child(2_i32_k, status)
   call check_ok(status, 'tree 1 to_child(2)')
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'tree 1 subdivide child 2')
   call cursor%to_root(status)
   call check_ok(status, 'tree 1 to_root')

   ! 4. Tree 3: single unrefined root leaf cell
   call new_non_oriented_cursor(grid, 3_id_k, cursor, create=.true., status=status)
   call check_ok(status, 'new_non_oriented_cursor(tree 3)')

   ! 5. Assign contiguous implicit global indices [0 .. N_cells - 1] across active trees
   call grid%initialize_local_index_node(status)
   call check_ok(status, 'initialize_local_index_node')

   ! 6. Run whole-grid structural invariant validation
   call grid%validate(status)
   call check_ok(status, 'grid%validate')

   ! 7. Verify expected topology and indexing counts
   if (grid%get_dimension() /= 2_i32_k) error stop 'Expected dimension == 2'
   if (grid%get_number_of_children() /= 4_i32_k) error stop 'Expected 4 children'
   if (grid%get_max_number_of_trees() /= 4_id_k) error stop 'Expected max 4 trees'
   if (grid%get_number_of_non_empty_trees() /= 3_id_k) error stop 'Expected 3 active trees'
   if (grid%get_number_of_cells() /= 27_id_k) error stop 'Expected 27 total vertices'
   if (grid%get_number_of_leaves() /= 21_id_k) error stop 'Expected 21 leaf cells'
   if (grid%get_number_of_levels() /= 4_id_k) error stop 'Expected 4 levels'
   if (grid%get_max_global_index() /= 26_id_k) error stop 'Expected max global ID == 26'

   write(*, '(a)') '=== Example 01: Basic 2D Quadtree ==='
   write(*, '(a,i0,a,i0,a,i0,a)') '  Active plane axes : (', axis1, ', ', axis2, &
      '), children/node = ', grid%get_number_of_children(), ''
   write(*, '(a,i0,a,i0)') '  Non-empty trees   : ', &
      grid%get_number_of_non_empty_trees(), ' / ', grid%get_max_number_of_trees()
   write(*, '(a,i0)') '  Total cells       : ', grid%get_number_of_cells()
   write(*, '(a,i0)') '  Leaf cells        : ', grid%get_number_of_leaves()
   write(*, '(a,i0)') '  Max levels        : ', grid%get_number_of_levels()
   write(*, '(a,i0)') '  Max global ID     : ', grid%get_max_global_index()
   write(*, '(a)') '  Validation        : PASSED'

contains

   subroutine check_ok(st, context)
      type(htg_status_t), intent(in) :: st
      character(len=*), intent(in) :: context
      if (.not. st%ok()) then
         error stop context // ': ' // trim(st%message)
      end if
   end subroutine check_ok

   recursive subroutine refine_corner_recursively(cur, target_child, max_depth)
      type(htg_cursor_t), intent(inout) :: cur
      integer(i32_k), intent(in) :: target_child, max_depth
      type(htg_status_t) :: st

      if (cur%get_level() >= max_depth) return

      call cur%subdivide_leaf(st)
      call check_ok(st, 'refine_corner_recursively: subdivide_leaf')

      call cur%to_child(target_child, st)
      call check_ok(st, 'refine_corner_recursively: to_child')

      call refine_corner_recursively(cur, target_child, max_depth)

      call cur%to_parent(st)
      call check_ok(st, 'refine_corner_recursively: to_parent')
   end subroutine refine_corner_recursively

end program example_01_basic_quadtree
