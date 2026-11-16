!> Example 05: Breadth-First Export for a Serializer (the fvtkhtg -> fvtkhdf seam)
!>
!> Demonstrates:
!>   - `grid%export_breadth_first(bf)`: the per-cell breadth-first description
!>     that the fvtkhdf writer's `write_mesh_topology(tree_ids, xtree, refined [, mask])`
!>     takes unchanged
!>   - `bf%gather(field, values)`: CellData permuted from fvtkhtg global-ID order
!>     into export order, ready for `write_cell_data(name, values)`
!>   - the index convention at the seam: fvtkhtg is 0-based everywhere; the
!>     export's `tree_ids` and `xtree` are 1-based, as the fvtkhdf writer expects
!>
!> No HDF5 or fvtkhdf is needed to run this example. With fvtkhdf, the arrays
!> printed below are passed straight to the writer:
!>
!>    call file%write_mesh(grid%get_dimensions(), grid%get_branch_factor())
!>    call file%write_mesh_geometry(x, y, z)
!>    call file%write_mesh_topology(bf%tree_ids, bf%xtree, bf%refined, bf%mask)
!>    call file%write_cell_data('Level', level)
!>
!> or, with the optional `fvtkhtg_vtkhdf` bridge library, in one call:
!>
!>    call htg_vtkhdf_write_static(file, grid, status)
program example_05_breadth_first_export
   use fvtkhtg, only: &
      htg_grid_t, htg_cursor_t, htg_cell_data_t, htg_field_array_t, &
      htg_breadth_first_t, htg_status_t, new_non_oriented_cursor, &
      HTG_TYPE_INT32, id_k, i32_k, r64_k
   implicit none

   type(htg_grid_t), target :: grid
   type(htg_cursor_t) :: cursor
   type(htg_cell_data_t), pointer :: cell_data => null()
   type(htg_field_array_t), pointer :: level_field => null()
   type(htg_breadth_first_t) :: bf
   type(htg_status_t) :: status
   integer(i32_k), allocatable :: level(:)
   integer(id_k) :: gid

   ! 1. A 2 x 1 quadtree forest. In tree 0, child 3 is refined before child 0,
   !    so fvtkhtg's storage order differs from breadth-first order.
   call grid%initialize()
   call grid%set_dimensions(3_i32_k, 2_i32_k, 1_i32_k, status)
   call check_ok(status, 'set_dimensions')
   call grid%set_x_coordinates([0.0_r64_k, 1.0_r64_k, 2.0_r64_k], status)
   call grid%set_y_coordinates([0.0_r64_k, 1.0_r64_k], status)

   call new_non_oriented_cursor(grid, 0_id_k, cursor, create=.true., status=status)
   call check_ok(status, 'open tree 0')
   call cursor%subdivide_leaf(status)
   call cursor%to_child(3_i32_k, status)
   call cursor%subdivide_leaf(status)
   call cursor%to_root(status)
   call cursor%to_child(0_i32_k, status)
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'refine tree 0')

   call new_non_oriented_cursor(grid, 1_id_k, cursor, create=.true., status=status)
   call check_ok(status, 'open tree 1')

   ! 2. Global indices and a cell field holding each cell's refinement level
   call grid%initialize_local_index_node(status)
   call grid%get_cell_data(cell_data)
   call cell_data%add_array('Level', HTG_TYPE_INT32, num_components=1_i32_k, &
      num_tuples=grid%get_number_of_cells(), status=status)
   call cell_data%get_array('Level', level_field, status)
   call set_levels(cursor, grid, level_field)

   call grid%validate(status)
   call check_ok(status, 'validate')

   ! 3. Export and gather
   call grid%export_breadth_first(bf, status)
   call check_ok(status, 'export_breadth_first')
   call bf%gather(level_field, level, status)
   call check_ok(status, 'gather Level')

   print '(a,i0,a,i0,a)', 'Exported ', bf%number_of_trees(), ' trees, ', &
      bf%number_of_cells(), ' cells'
   print '(a,*(1x,i0))', '  tree_ids (1-based root cells):', bf%tree_ids
   print '(a,*(1x,i0))', '  xtree    (1-based offsets)   :', bf%xtree
   print '(a,*(1x,l1))', '  refined                      :', bf%refined
   print '(a,*(1x,i0))', '  global_ids (0-based, fvtkhtg):', bf%global_ids
   print '(a,*(1x,i0))', '  Level gathered (export order):', level
   print '(a,l1)', '  identity_order (no gather needed): ', bf%identity_order

   ! The export is in breadth-first order: levels never decrease within a tree
   do gid = 2_id_k, bf%xtree(2) - 1_id_k
      if (level(gid) < level(gid - 1_id_k)) error stop 'export is not breadth-first'
   end do
   print '(a)', 'Example 05 completed successfully.'

contains

   subroutine check_ok(st, what)
      type(htg_status_t), intent(in) :: st
      character(len=*), intent(in) :: what
      if (.not. st%ok()) error stop what // ': ' // trim(st%message)
   end subroutine check_ok

   !> Walk tree 0 depth-first and store each cell's level in the field
   subroutine set_levels(c, g, field)
      type(htg_cursor_t), intent(inout) :: c
      type(htg_grid_t), target, intent(inout) :: g
      type(htg_field_array_t), intent(inout) :: field
      type(htg_status_t) :: st
      integer(i32_k) :: i, j

      call new_non_oriented_cursor(g, 0_id_k, c, create=.false., status=st)
      call field%set_scalar_i32(c%get_global_node_index(), 0_i32_k, st)
      do i = 0_i32_k, 3_i32_k
         call c%to_root(st)
         call c%to_child(i, st)
         call field%set_scalar_i32(c%get_global_node_index(), 1_i32_k, st)
         if (.not. c%is_leaf()) then
            do j = 0_i32_k, 3_i32_k
               call c%to_child(j, st)
               call field%set_scalar_i32(c%get_global_node_index(), 2_i32_k, st)
               call c%to_parent(st)
            end do
         end if
      end do
      call new_non_oriented_cursor(g, 1_id_k, c, create=.false., status=st)
      call field%set_scalar_i32(c%get_global_node_index(), 0_i32_k, st)
   end subroutine set_levels

end program example_05_breadth_first_export
