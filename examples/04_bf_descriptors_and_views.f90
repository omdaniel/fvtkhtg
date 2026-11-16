!> Example 04: Breadth-First Order Descriptors & Zero-Copy Serialization Views
!>
!> Demonstrates:
!>   - Building trees from concatenated Breadth-First Order bit descriptors
!>     (`build_from_breadth_first_order_descriptor`), including VTK v9.7+'s
!>     canonical `number_of_bits = 0` single-leaf root convention
!>   - Exporting Breadth-First Order descriptors (`compute_breadth_first_order_descriptor`)
!>     both losslessly across a forest and pruned with `input_mask` + `depth_limiter`
!>   - Extracting zero-copy `0:`-bounded pointer views (`get_parent_to_elder_child_view`,
!>     `get_vertices_per_depth_view`, `get_raw_bytes` MSB-first packed bytes,
!>     coordinate views, and `CellData` 2D views) for external file/MPI serializers
program example_04_bf_descriptors_and_views
   use fvtkhtg, only: &
      htg_grid_t, htg_tree_t, htg_bit_array_t, htg_cell_data_t, htg_field_array_t, &
      htg_status_t, HTG_TYPE_REAL64, HTG_TYPE_INT32, HTG_TREE_LEAF_SENTINEL, &
      id_k, i32_k, byte_k, r64_k
   implicit none

   type(htg_grid_t), target :: grid
   type(htg_tree_t), pointer :: tree_ptr => null()
   type(htg_bit_array_t), target :: in_desc, out_desc, pruned_desc, prune_mask
   type(htg_cell_data_t), pointer :: cell_data => null()
   type(htg_field_array_t), pointer :: pressure => null()
   type(htg_field_array_t), pointer :: mat_id => null()
   type(htg_status_t) :: status

   integer(id_k), allocatable :: out_vpd(:), out_id_map(:)
   integer(id_k), allocatable :: pruned_vpd(:), pruned_id_map(:)

   real(r64_k), pointer :: x_view(:) => null()
   real(r64_k), pointer :: y_view(:) => null()
   real(r64_k), pointer :: z_view(:) => null()
   integer(id_k), pointer :: elder_view(:) => null()
   integer(id_k), pointer :: vpd_view(:) => null()
   integer(byte_k), pointer :: raw_bytes(:) => null()
   real(r64_k), pointer :: p_view(:, :) => null()
   integer(i32_k), pointer :: m_view(:, :) => null()

   real(r64_k) :: x(0:3), y(0:1), z(0:0)
   logical :: bits(0:5)
   integer(id_k) :: pos, tid, gid, n_cells
   integer(i32_k) :: unsigned_byte0

   ! 1. Configure a 3 x 1 2D quadtree grid (point dimensions 4 x 2 x 1)
   call grid%initialize()
   call grid%set_branch_factor(2_i32_k, status)
   call check_ok(status, 'set_branch_factor')
   call grid%set_dimensions(4_i32_k, 2_i32_k, 1_i32_k, status)
   call check_ok(status, 'set_dimensions')

   x = [0.0_r64_k, 1.0_r64_k, 2.0_r64_k, 3.0_r64_k]
   y = [0.0_r64_k, 1.0_r64_k]
   z = [0.0_r64_k]
   call grid%set_x_coordinates(x, status)
   call check_ok(status, 'set_x_coordinates')
   call grid%set_y_coordinates(y, status)
   call check_ok(status, 'set_y_coordinates')
   call grid%set_z_coordinates(z, status)
   call check_ok(status, 'set_z_coordinates')

   ! 2. Pack a concatenated 6-bit Breadth-First Order forest descriptor stream:
   !    - Tree 0 (5 bits at start_index = 0): `1 | 1 0 0 1`
   !        depth 0: [1] -> root refined into 4 children (1..4)
   !        depth 1: [1, 0, 0, 1] -> children 0 and 3 refined (+8 leaves at depth 2)
   !        depth 2: all 8 leaves omitted (deepest level is never encoded in VTK!)
   !    - Tree 1 (0 bits at start_index = 5): ``
   !        VTK v9.7+ canonical single-leaf root uses `number_of_bits = 0`!
   !    - Tree 2 (1 bit  at start_index = 5): `1`
   !        depth 0: [1] -> root refined into 4 leaves at depth 1
   bits = [.true., .true., .false., .false., .true., .true.]
   call in_desc%initialize(6_id_k, .false., status)
   call check_ok(status, 'in_desc%initialize')
   do gid = 0_id_k, 5_id_k
      call in_desc%set_value(gid, bits(gid), status)
      call check_ok(status, 'in_desc%set_value')
   end do

   call grid%get_tree(0_id_k, tree_ptr, create=.true., status=status)
   call check_ok(status, 'get_tree(0)')
   call tree_ptr%build_from_breadth_first_order_descriptor( &
      in_desc, 5_id_k, start_index=0_id_k, status=status)
   call check_ok(status, 'build tree 0 from 5-bit slice')

   call grid%get_tree(1_id_k, tree_ptr, create=.true., status=status)
   call check_ok(status, 'get_tree(1)')
   call tree_ptr%build_from_breadth_first_order_descriptor( &
      in_desc, 0_id_k, start_index=5_id_k, status=status)
   call check_ok(status, 'build tree 1 from 0-bit single-leaf slice')
   if (tree_ptr%get_number_of_levels() /= 1_id_k .or. &
       tree_ptr%get_number_of_vertices() /= 1_id_k) then
      error stop 'Expected 1-level, 1-vertex single-leaf root for 0-bit descriptor'
   end if

   call grid%get_tree(2_id_k, tree_ptr, create=.true., status=status)
   call check_ok(status, 'get_tree(2)')
   call tree_ptr%build_from_breadth_first_order_descriptor( &
      in_desc, 1_id_k, start_index=5_id_k, status=status)
   call check_ok(status, 'build tree 2 from 1-bit slice')

   ! 3. Assign implicit global indices [0..18] and populate CellData arrays
   call grid%initialize_local_index_node(status)
   call check_ok(status, 'initialize_local_index_node')
   n_cells = grid%get_number_of_cells() ! 13 + 1 + 5 = 19 cells

   call grid%get_cell_data(cell_data)
   call cell_data%add_array('Pressure', HTG_TYPE_REAL64, &
      num_components=1_i32_k, num_tuples=n_cells, status=status)
   call check_ok(status, 'add_array Pressure')
   call cell_data%add_array('MaterialId', HTG_TYPE_INT32, &
      num_components=1_i32_k, num_tuples=n_cells, status=status)
   call check_ok(status, 'add_array MaterialId')

   call cell_data%get_array('Pressure', pressure, status)
   call check_ok(status, 'get_array Pressure')
   call cell_data%get_array('MaterialId', mat_id, status)
   call check_ok(status, 'get_array MaterialId')

   do gid = 0_id_k, n_cells - 1_id_k
      call pressure%set_scalar_r64(gid, 101.325_r64_k + real(gid, r64_k), status)
      call check_ok(status, 'set Pressure')
      call mat_id%set_scalar_i32(gid, int(mod(gid, 3_id_k), i32_k), status)
      call check_ok(status, 'set MaterialId')
   end do

   ! 4. Export concatenated Breadth-First Order descriptor across all active trees
   do pos = 0_id_k, grid%get_number_of_non_empty_trees() - 1_id_k
      call grid%get_tree_by_pos(pos, tid, tree_ptr, status)
      call check_ok(status, 'get_tree_by_pos')
      call tree_ptr%compute_breadth_first_order_descriptor( &
         out_desc, out_vpd, out_id_map, status=status)
      call check_ok(status, 'compute_breadth_first_order_descriptor')
   end do

   if (out_desc%get_number_of_values() /= 6_id_k) then
      error stop 'Expected 6 exported descriptor bits across the 3 trees'
   end if
   if (lbound(out_vpd, 1) /= 0 .or. size(out_vpd) /= 6) then
      error stop 'Expected 6 depth entries (3 + 1 + 2) with 0: lower bound'
   end if
   if (lbound(out_id_map, 1) /= 0 .or. size(out_id_map) /= 19) then
      error stop 'Expected 19 mapped global IDs with 0: lower bound'
   end if

   ! Demonstrate pruned export on Tree 0 with input_mask (masking coarse gid = 4)
   ! and depth_limiter = 2 (pruning gid = 4's 4 children -> 9 visited instead of 13)
   call prune_mask%initialize(n_cells, .false., status)
   call check_ok(status, 'prune_mask%initialize')
   call prune_mask%set_value(4_id_k, .true., status)
   call check_ok(status, 'prune_mask%set_value(4)')

   call grid%get_tree(0_id_k, tree_ptr, create=.false., status=status)
   call check_ok(status, 'get_tree(0) for pruned export')
   call tree_ptr%compute_breadth_first_order_descriptor( &
      pruned_desc, pruned_vpd, pruned_id_map, &
      input_mask=prune_mask, depth_limiter=2_i32_k, status=status)
   call check_ok(status, 'pruned compute_breadth_first_order_descriptor')
   if (size(pruned_id_map) /= 9 .or. pruned_vpd(2) /= 4_id_k) then
      error stop 'Expected masked coarse gid=4 subtree to be pruned to 9 vertices'
   end if

   ! 5. Extract Zero-Copy 0-Based Pointer Views for External Serialization
   call grid%get_x_coordinates_view(x_view)
   call grid%get_y_coordinates_view(y_view)
   call grid%get_z_coordinates_view(z_view)
   call tree_ptr%get_parent_to_elder_child_view(elder_view)
   call tree_ptr%get_vertices_per_depth_view(vpd_view)
   call out_desc%get_raw_bytes(raw_bytes)
   call pressure%get_view_r64(p_view, status)
   call check_ok(status, 'pressure%get_view_r64')
   call mat_id%get_view_i32(m_view, status)
   call check_ok(status, 'mat_id%get_view_i32')

   if (lbound(x_view, 1) /= 0 .or. lbound(y_view, 1) /= 0 .or. &
       lbound(z_view, 1) /= 0) then
      error stop 'Expected 0-based coordinate pointer views'
   end if
   if (lbound(elder_view, 1) /= 0 .or. elder_view(0) /= 1_id_k .or. &
       elder_view(2) /= HTG_TREE_LEAF_SENTINEL) then
      error stop 'Unexpected Tree 0 parent_to_elder_child view contents'
   end if
   if (lbound(vpd_view, 1) /= 0 .or. vpd_view(2) /= 8_id_k) then
      error stop 'Unexpected Tree 0 vertices_per_depth view contents'
   end if

   ! MSB-first packed bits `110011` + 2 zero slack bits = 0b11001100 = 0xCC (204)
   if (lbound(raw_bytes, 1) /= 0 .or. size(raw_bytes) /= 1) then
      error stop 'Expected 1 packed MSB-first byte with 0: lower bound'
   end if
   unsigned_byte0 = iand(int(raw_bytes(0), i32_k), 255_i32_k)
   if (unsigned_byte0 /= 204_i32_k) then
      error stop 'Expected MSB-first packed byte 0xCC (204) for bits 110011'
   end if

   if (lbound(p_view, 1) /= 0 .or. lbound(p_view, 2) /= 0 .or. &
       lbound(m_view, 1) /= 0 .or. lbound(m_view, 2) /= 0) then
      error stop 'Expected (0:nc-1, 0:nt-1) bounds on 2D CellData views'
   end if

   call verify_zero_based_dummy_view(x_view)

   call grid%validate(status)
   call check_ok(status, 'grid%validate')

   write(*, '(a)') '=== Example 04: Breadth-First Descriptors & Zero-Copy Views ==='
   write(*, '(a,i0,a,z2.2)') '  Forest BF bits    : ', &
      out_desc%get_number_of_values(), ', MSB-first byte[0] = 0x', unsigned_byte0
   write(*, '(a,i0,a,i0)') '  Tree 0 vertices   : full = ', &
      tree_ptr%get_number_of_vertices(), ', mask-pruned = ', size(pruned_id_map)
   write(*, '(a,i0,a,i0,a)') '  Pressure 2D view  : (0:', ubound(p_view, 1), &
      ', 0:', ubound(p_view, 2), ')'
   write(*, '(a)') '  Validation        : PASSED'

contains

   subroutine check_ok(st, context)
      type(htg_status_t), intent(in) :: st
      character(len=*), intent(in) :: context
      if (.not. st%ok()) then
         error stop context // ': ' // trim(st%message)
      end if
   end subroutine check_ok

   ! CRITICAL FORTRAN RULE: Always declare `0:` on assumed-shape dummy arguments
   ! receiving fvtkhtg views; otherwise Fortran resets the dummy's lower bound to 1!
   subroutine verify_zero_based_dummy_view(coords)
      real(r64_k), intent(in) :: coords(0:)
      if (lbound(coords, 1) /= 0 .or. ubound(coords, 1) /= 3) then
         error stop 'Expected 0:3 bounds inside dummy view consumer'
      end if
      if (abs(coords(0) - 0.0_r64_k) > 1.0e-12_r64_k .or. &
          abs(coords(3) - 3.0_r64_k) > 1.0e-12_r64_k) then
         error stop 'Unexpected coordinate values in zero-copy view'
      end if
   end subroutine verify_zero_based_dummy_view

end program example_04_bf_descriptors_and_views
