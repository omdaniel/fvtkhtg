!> Example 03: Cell Masking, Material Interfaces & Lazy PureMask Caching
!>
!> Demonstrates:
!>   - Masking leaf and coarse cells (`cursor%set_mask`, `is_masked`, `grid%has_mask`)
!>     and VTK v9.7+ precondition enforcement (`HTG_ERR_MASK` on masked subdivide/descend)
!>   - Configuring Material Interface metadata (`set_has_interface`,
!>     `set_interface_normals_name`, `set_interface_intercepts_name`) with 3-component
!>     `Normals` and `Intercepts` arrays in `CellData`
!>   - Querying lazy `PureMask` (`get_pure_mask`, where `.false.` = pure and
!>     `.true.` = masked or mixed/interfaced) and observing automatic cache invalidation
program example_03_mask_puremask_interface
   use fvtkhtg, only: &
      htg_grid_t, htg_cursor_t, htg_cell_data_t, htg_field_array_t, &
      htg_bit_array_t, htg_status_t, new_non_oriented_cursor, &
      HTG_TYPE_REAL64, HTG_ERR_MASK, id_k, i32_k, r64_k
   implicit none

   type(htg_grid_t), target :: grid
   type(htg_cursor_t) :: cursor
   type(htg_cell_data_t), pointer :: cell_data => null()
   type(htg_field_array_t), pointer :: normals => null()
   type(htg_field_array_t), pointer :: intercepts => null()
   type(htg_bit_array_t), pointer :: pure_mask => null()
   type(htg_status_t) :: status, expected_err

   real(r64_k) :: x(0:2), y(0:1), z(0:0), nrm(0:2), icp(0:2)
   integer(id_k) :: gid, n_cells
   logical :: is_masked_now

   ! 1. Build a 2-tree 2D XY quadtree grid (dimensions 3 x 2 x 1 -> 2 x 1 root cells)
   call grid%initialize()
   call grid%set_branch_factor(2_i32_k, status)
   call check_ok(status, 'set_branch_factor')
   call grid%set_dimensions(3_i32_k, 2_i32_k, 1_i32_k, status)
   call check_ok(status, 'set_dimensions')

   x = [0.0_r64_k, 1.0_r64_k, 2.0_r64_k]
   y = [0.0_r64_k, 1.0_r64_k]
   z = [0.0_r64_k]
   call grid%set_x_coordinates(x, status)
   call check_ok(status, 'set_x_coordinates')
   call grid%set_y_coordinates(y, status)
   call check_ok(status, 'set_y_coordinates')
   call grid%set_z_coordinates(z, status)
   call check_ok(status, 'set_z_coordinates')

   ! Tree 0: subdivide root (0 -> 1..4), child 0 (1 -> 5..8), and child 1 (2 -> 9..12)
   call new_non_oriented_cursor(grid, 0_id_k, cursor, create=.true., status=status)
   call check_ok(status, 'cursor tree 0')
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'subdivide root')

   call cursor%to_child(0_i32_k, status)
   call check_ok(status, 'to_child(0)')
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'subdivide child 0')
   call cursor%to_parent(status)
   call check_ok(status, 'to_parent from child 0')

   call cursor%to_child(1_i32_k, status)
   call check_ok(status, 'to_child(1)')
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'subdivide child 1')
   call cursor%to_root(status)
   call check_ok(status, 'to_root')

   ! Tree 1: unrefined single-leaf root (global ID 13 after indexing)
   call new_non_oriented_cursor(grid, 1_id_k, cursor, create=.true., status=status)
   call check_ok(status, 'cursor tree 1')

   call grid%initialize_local_index_node(status)
   call check_ok(status, 'initialize_local_index_node')
   n_cells = grid%get_number_of_cells() ! 14 total vertices (0..13)

   ! 2. Attach 3-component Material Interface Normals and Intercepts in CellData
   call grid%get_cell_data(cell_data)
   call cell_data%add_array('Normals', HTG_TYPE_REAL64, &
      num_components=3_i32_k, num_tuples=n_cells, status=status)
   call check_ok(status, 'add_array Normals')
   call cell_data%add_array('Intercepts', HTG_TYPE_REAL64, &
      num_components=3_i32_k, num_tuples=n_cells, status=status)
   call check_ok(status, 'add_array Intercepts')

   call cell_data%get_array('Normals', normals, status)
   call check_ok(status, 'get_array Normals')
   call cell_data%get_array('Intercepts', intercepts, status)
   call check_ok(status, 'get_array Intercepts')

   ! Default all cells to pure material sentinel: Intercepts(2) = 2.0 (>= 2.0 is pure)
   nrm = [0.0_r64_k, 0.0_r64_k, 0.0_r64_k]
   icp = [0.0_r64_k, 0.0_r64_k, 2.0_r64_k]
   do gid = 0_id_k, n_cells - 1_id_k
      call normals%set_tuple_r64(gid, nrm, status)
      call check_ok(status, 'init Normals')
      call intercepts%set_tuple_r64(gid, icp, status)
      call check_ok(status, 'init Intercepts')
   end do

   ! Mark leaf gid = 5 (under coarse parent gid = 1) with an active material interface:
   ! In VTK v9.7+, Intercepts(2) < 2.0 (-1, 0, or 1) indicates a mixed/interfaced leaf
   nrm = [0.7071067811865476_r64_k, 0.7071067811865476_r64_k, 0.0_r64_k]
   icp = [-0.25_r64_k, 0.0_r64_k, 0.0_r64_k]
   call normals%set_tuple_r64(5_id_k, nrm, status)
   call check_ok(status, 'set interface Normals(5)')
   call intercepts%set_tuple_r64(5_id_k, icp, status)
   call check_ok(status, 'set interface Intercepts(5)')

   call grid%set_has_interface(.true.)
   call grid%set_interface_normals_name('Normals')
   call grid%set_interface_intercepts_name('Intercepts')

   ! 3. Mask coarse node gid = 2 (child 1 of Tree 0) and leaf gid = 3 (child 2 of Tree 0)
   call new_non_oriented_cursor(grid, 0_id_k, cursor, create=.false., status=status)
   call check_ok(status, 're-open cursor tree 0')

   call cursor%to_child(1_i32_k, status)
   call check_ok(status, 'to_child(1) for masking')
   call cursor%set_mask(.true., status)
   call check_ok(status, 'set_mask on coarse child 1')
   if (.not. grid%has_mask()) error stop 'Expected grid mask to be active'
   is_masked_now = cursor%is_masked(status)
   call check_ok(status, 'is_masked on coarse child 1')
   if (.not. is_masked_now) error stop 'Expected coarse child 1 to be masked'

   ! Descending into a masked coarse node is rejected with HTG_ERR_MASK
   call cursor%to_child(0_i32_k, expected_err)
   if (expected_err%code /= HTG_ERR_MASK) then
      error stop 'Expected HTG_ERR_MASK when descending into masked coarse cell'
   end if
   call cursor%to_parent(status)
   call check_ok(status, 'to_parent after masked coarse')

   call cursor%to_child(2_i32_k, status)
   call check_ok(status, 'to_child(2) for masking')
   call cursor%set_mask(.true., status)
   call check_ok(status, 'set_mask on leaf child 2')

   ! Subdividing a masked leaf is rejected with HTG_ERR_MASK
   call cursor%subdivide_leaf(expected_err)
   if (expected_err%code /= HTG_ERR_MASK) then
      error stop 'Expected HTG_ERR_MASK when subdividing masked leaf cell'
   end if

   ! 4. Compute and inspect lazy PureMask (0/.false. = pure, 1/.true. = impure)
   if (grid%is_pure_mask_valid()) error stop 'Expected PureMask cache initially invalid'
   call grid%get_pure_mask(pure_mask, status)
   call check_ok(status, 'get_pure_mask')
   if (.not. grid%is_pure_mask_valid()) error stop 'Expected PureMask cache valid'

   ! Verify PureMask bottom-up semantics:
   !   gid 5 (interfaced leaf)       -> .true.  (impure)
   !   gid 6 (pure sibling leaf)     -> .false. (pure)
   !   gid 1 (coarse parent of 5..8) -> .true.  (impure propagated from child 5)
   !   gid 2 (masked coarse node)    -> .true.  (short-circuits; child 9 stays .false.)
   !   gid 3 (masked leaf)           -> .true.  (impure)
   !   gid 4 (unmasked pure leaf)    -> .false. (pure)
   !   gid 13 (Tree 1 pure root)     -> .false. (pure)
   if (.not. pure_mask%get_value(5_id_k, status)) error stop 'Expected gid 5 impure'
   if (pure_mask%get_value(6_id_k, status)) error stop 'Expected gid 6 pure'
   if (.not. pure_mask%get_value(1_id_k, status)) error stop 'Expected gid 1 impure'
   if (.not. pure_mask%get_value(2_id_k, status)) error stop 'Expected gid 2 impure'
   if (pure_mask%get_value(9_id_k, status)) error stop 'Expected unvisited gid 9 false'
   if (.not. pure_mask%get_value(3_id_k, status)) error stop 'Expected gid 3 impure'
   if (pure_mask%get_value(4_id_k, status)) error stop 'Expected gid 4 pure'
   if (pure_mask%get_value(13_id_k, status)) error stop 'Expected Tree 1 root pure'

   ! 5. Mutating mask or topology automatically invalidates cached PureMask
   call new_non_oriented_cursor(grid, 1_id_k, cursor, create=.false., status=status)
   call check_ok(status, 'cursor tree 1 for mask update')
   call cursor%set_mask(.true., status)
   call check_ok(status, 'mask tree 1 root')
   if (grid%is_pure_mask_valid()) then
      error stop 'Expected PureMask cache invalidated after cursor%set_mask'
   end if

   call grid%get_pure_mask(pure_mask, status)
   call check_ok(status, 'recompute get_pure_mask')
   if (.not. pure_mask%get_value(13_id_k, status)) then
      error stop 'Expected Tree 1 root impure after masking'
   end if

   call grid%validate(status)
   call check_ok(status, 'grid%validate')

   write(*, '(a)') '=== Example 03: Masking, Material Interfaces & PureMask ==='
   write(*, '(a,l1,a,a,a,a)') '  Interface active  : ', grid%get_has_interface(), &
      ' (', grid%get_interface_normals_name(), ', ', &
      grid%get_interface_intercepts_name() // ')'
   write(*, '(a,i0)') '  PureMask bits     : ', pure_mask%get_number_of_values()
   write(*, '(a,l1)') '  gid=5 (interface) : impure = ', pure_mask%get_value(5_id_k)
   write(*, '(a,l1)') '  gid=6 (pure leaf) : impure = ', pure_mask%get_value(6_id_k)
   write(*, '(a,l1)') '  gid=1 (ancestor)  : impure = ', pure_mask%get_value(1_id_k)
   write(*, '(a,l1)') '  gid=2 (masked)    : impure = ', pure_mask%get_value(2_id_k)
   write(*, '(a)')    '  Validation        : PASSED'

contains

   subroutine check_ok(st, context)
      type(htg_status_t), intent(in) :: st
      character(len=*), intent(in) :: context
      if (.not. st%ok()) then
         error stop context // ': ' // trim(st%message)
      end if
   end subroutine check_ok

end program example_03_mask_puremask_interface
