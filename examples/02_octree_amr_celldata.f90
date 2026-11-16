!> Example 02: 3D Sparse Multi-Tree Octree AMR, Global Indexing & Dynamic CellData
!>
!> Demonstrates:
!>   - 3D sparse multi-tree octree setup (`branch_factor = 2`, `k = 8` children)
!>   - Level-zero `(i, j, k) <-> tree_id` mapping in both default KJI (`i` fastest)
!>     and transposed IJK (`transposed_root_indexing = .true.`, `k` fastest) modes
!>   - Combining implicit (`set_global_index_start`) and explicit
!>     (`set_global_index_from_local`) global node indexing across trees
!>   - Attaching dynamic `CellData` scalar (`HTG_TYPE_REAL64`) and 3-component vector
!>     (`HTG_TYPE_REAL64` and `HTG_TYPE_INT32`) arrays that grow automatically as
!>     leaves subdivide
program example_02_octree_amr_celldata
   use fvtkhtg, only: &
      htg_grid_t, htg_cursor_t, htg_cell_data_t, htg_field_array_t, htg_status_t, &
      new_non_oriented_cursor, HTG_TYPE_REAL64, HTG_TYPE_INT32, &
      id_k, i32_k, r64_k
   implicit none

   type(htg_grid_t), target :: grid, grid_ijk
   type(htg_cursor_t) :: cursor
   type(htg_cell_data_t), pointer :: cell_data => null()
   type(htg_field_array_t), pointer :: density => null()
   type(htg_field_array_t), pointer :: velocity => null()
   type(htg_field_array_t), pointer :: cell_ijk => null()
   type(htg_status_t) :: status

   real(r64_k) :: x(0:3), y(0:2), z(0:2), vel_vec(0:2), vel_out(0:2)
   integer(i32_k) :: tag_vec(0:2), tag_out(0:2)
   integer(id_k) :: tid_kji, tid_ijk, tid_a, tid_b, gid, base_explicit
   integer(i32_k) :: ri, rj, rk, ichild

   ! 1. Compare KJI (default) vs IJK (transposed_root_indexing) level-zero mapping
   !    on a 3 x 2 x 2 cell grid (point dimensions 4 x 3 x 3 -> 12 tree slots)
   call grid_ijk%initialize()
   call grid_ijk%set_dimensions(4_i32_k, 3_i32_k, 3_i32_k, status)
   call check_ok(status, 'grid_ijk set_dimensions')
   call grid_ijk%set_transposed_root_indexing(.true., status)
   call check_ok(status, 'grid_ijk set_transposed_root_indexing')

   call grid%initialize()
   call grid%set_branch_factor(2_i32_k, status)
   call check_ok(status, 'set_branch_factor')
   call grid%set_dimensions(4_i32_k, 3_i32_k, 3_i32_k, status)
   call check_ok(status, 'set_dimensions')

   x = [0.0_r64_k, 1.0_r64_k, 2.0_r64_k, 3.0_r64_k]
   y = [0.0_r64_k, 1.0_r64_k, 2.0_r64_k]
   z = [0.0_r64_k, 1.0_r64_k, 2.0_r64_k]
   call grid%set_x_coordinates(x, status)
   call check_ok(status, 'set_x_coordinates')
   call grid%set_y_coordinates(y, status)
   call check_ok(status, 'set_y_coordinates')
   call grid%set_z_coordinates(z, status)
   call check_ok(status, 'set_z_coordinates')

   ! Verify root coordinate formulas for cell (i=2, j=1, k=1):
   !   KJI (default): tree_id = i + j*nx + k*nx*ny = 2 + 1*3 + 1*6 = 11
   !   IJK (transp.): tree_id = k + j*nz + i*nz*ny = 1 + 1*2 + 2*4 = 11 (corner)
   ! For (i=2, j=0, k=1): KJI = 2 + 0 + 6 = 8; IJK = 1 + 0 + 8 = 9
   call grid%get_index_from_level_zero_coordinates(tid_kji, 2_i32_k, 0_i32_k, 1_i32_k, status)
   call check_ok(status, 'KJI get_index_from_level_zero_coordinates')
   call grid_ijk%get_index_from_level_zero_coordinates(tid_ijk, 2_i32_k, 0_i32_k, 1_i32_k, status)
   call check_ok(status, 'IJK get_index_from_level_zero_coordinates')
   if (tid_kji /= 8_id_k .or. tid_ijk /= 9_id_k) then
      error stop 'Unexpected KJI/IJK root index mapping'
   end if

   call grid%get_level_zero_coordinates_from_index(tid_kji, ri, rj, rk, status)
   call check_ok(status, 'get_level_zero_coordinates_from_index')
   if (ri /= 2_i32_k .or. rj /= 0_i32_k .or. rk /= 1_i32_k) then
      error stop 'Level-zero inverse coordinate mismatch'
   end if

   ! 2. Attach dynamic CellData arrays BEFORE building trees so we can observe
   !    automatic tuple growth during cursor%subdivide_leaf and explicit assignment
   call grid%get_cell_data(cell_data)
   call cell_data%add_array('Density', HTG_TYPE_REAL64, num_components=1_i32_k, status=status)
   call check_ok(status, 'add_array Density')
   call cell_data%set_active_scalars('Density', status)
   call check_ok(status, 'set_active_scalars')

   call cell_data%add_array('Velocity', HTG_TYPE_REAL64, num_components=3_i32_k, status=status)
   call check_ok(status, 'add_array Velocity')

   call cell_data%add_array('CellTag', HTG_TYPE_INT32, num_components=3_i32_k, status=status)
   call check_ok(status, 'add_array CellTag')

   call cell_data%get_array('Density', density, status)
   call check_ok(status, 'get_array Density')
   call cell_data%get_array('Velocity', velocity, status)
   call check_ok(status, 'get_array Velocity')
   call cell_data%get_array('CellTag', cell_ijk, status)
   call check_ok(status, 'get_array CellTag')

   ! 3. Sparse Tree A at (0, 0, 0) using IMPLICIT global indexing starting at 0
   call grid%get_index_from_level_zero_coordinates(tid_a, 0_i32_k, 0_i32_k, 0_i32_k, status)
   call check_ok(status, 'tree A index')
   call new_non_oriented_cursor(grid, tid_a, cursor, create=.true., status=status)
   call check_ok(status, 'cursor tree A')

   call cursor%set_global_index_start(0_id_k, status)
   call check_ok(status, 'tree A set_global_index_start(0)')
   if (cell_data%get_number_of_tuples() /= 1_id_k) then
      error stop 'Expected CellData tuples == 1 after root implicit start'
   end if

   ! Subdivide root -> +8 children (9 vertices -> CellData auto-grows to 9 tuples)
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'tree A subdivide root')
   if (cell_data%get_number_of_tuples() /= 9_id_k) then
      error stop 'Expected CellData auto-growth to 9 tuples'
   end if

   ! Subdivide octant 7 -> +8 children (17 vertices -> CellData auto-grows to 17 tuples)
   call cursor%to_child(7_i32_k, status)
   call check_ok(status, 'tree A to_child(7)')
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'tree A subdivide child 7')
   if (cell_data%get_number_of_tuples() /= 17_id_k) then
      error stop 'Expected CellData auto-growth to 17 tuples'
   end if

   ! Populate CellData on fine octant children of octant 7
   do ichild = 0_i32_k, 7_i32_k
      call cursor%to_child(ichild, status)
      call check_ok(status, 'tree A fine to_child')
      gid = cursor%get_global_node_index(status)
      call check_ok(status, 'tree A fine gid')

      call density%set_scalar_r64(gid, 1.25_r64_k + 0.1_r64_k * real(ichild, r64_k), status)
      call check_ok(status, 'set_scalar_r64 Density')

      vel_vec = [real(ichild, r64_k), -0.5_r64_k * real(ichild, r64_k), 2.0_r64_k]
      call velocity%set_tuple_r64(gid, vel_vec, status)
      call check_ok(status, 'set_tuple_r64 Velocity')

      tag_vec = [0_i32_k, cursor%get_level(), ichild]
      call cell_ijk%set_tuple_i32(gid, tag_vec, status)
      call check_ok(status, 'set_tuple_i32 CellTag')

      call cursor%to_parent(status)
      call check_ok(status, 'tree A fine to_parent')
   end do

   ! 4. Sparse Tree B at (2, 0, 1) using EXPLICIT permuted global IDs [100..108]
   !    (leaving 10 of the 12 root slots unallocated in the sparse forest)
   tid_b = tid_kji
   call new_non_oriented_cursor(grid, tid_b, cursor, create=.true., status=status)
   call check_ok(status, 'cursor tree B')
   call cursor%subdivide_leaf(status)
   call check_ok(status, 'tree B subdivide root')

   base_explicit = 100_id_k
   call cursor%set_global_index_from_local(base_explicit + 8_id_k, status)
   call check_ok(status, 'tree B explicit root gid = 108')

   do ichild = 0_i32_k, 7_i32_k
      call cursor%to_child(ichild, status)
      call check_ok(status, 'tree B to_child')
      gid = base_explicit + int(7_i32_k - ichild, id_k)
      call cursor%set_global_index_from_local(gid, status)
      call check_ok(status, 'tree B explicit child gid')

      call density%set_scalar_r64(gid, 10.0_r64_k + real(ichild, r64_k), status)
      call check_ok(status, 'tree B Density')
      vel_vec = [0.0_r64_k, 1.0_r64_k, real(ichild, r64_k)]
      call velocity%set_tuple_r64(gid, vel_vec, status)
      call check_ok(status, 'tree B Velocity')
      tag_vec = [int(tid_b, i32_k), cursor%get_level(), ichild]
      call cell_ijk%set_tuple_i32(gid, tag_vec, status)
      call check_ok(status, 'tree B CellTag')

      call cursor%to_parent(status)
      call check_ok(status, 'tree B to_parent')
   end do

   ! 5. Validate forest topology, disjoint implicit/explicit global IDs, and CellData
   call grid%validate(status)
   call check_ok(status, 'grid%validate')

   if (grid%get_number_of_non_empty_trees() /= 2_id_k) error stop 'Expected 2 sparse trees'
   if (grid%get_number_of_cells() /= 26_id_k) error stop 'Expected 26 total vertices'
   if (grid%get_max_global_index() /= 108_id_k) error stop 'Expected max global ID == 108'
   if (cell_data%get_number_of_tuples() /= 109_id_k) error stop 'Expected 109 CellData tuples'

   call velocity%get_tuple_r64(16_id_k, vel_out, status)
   call check_ok(status, 'get_tuple_r64(16)')
   call cell_ijk%get_tuple_i32(100_id_k, tag_out, status)
   call check_ok(status, 'get_tuple_i32(100)')
   if (abs(vel_out(0) - 7.0_r64_k) > 1.0e-12_r64_k .or. tag_out(2) /= 7_i32_k) then
      error stop 'CellData round-trip verification failed'
   end if

   write(*, '(a)') '=== Example 02: 3D Sparse Octree AMR & Dynamic CellData ==='
   write(*, '(a,i0,a,i0)') '  Root (2,0,1) ID   : KJI = ', tid_kji, ', IJK = ', tid_ijk
   write(*, '(a,i0,a,i0)') '  Sparse octrees    : ', &
      grid%get_number_of_non_empty_trees(), ' / ', grid%get_max_number_of_trees()
   write(*, '(a,i0,a,i0)') '  Cells / Leaves    : ', &
      grid%get_number_of_cells(), ' / ', grid%get_number_of_leaves()
   write(*, '(a,i0)') '  CellData arrays   : ', cell_data%get_number_of_arrays()
   write(*, '(a,i0)') '  CellData tuples   : ', cell_data%get_number_of_tuples()
   write(*, '(a,a)')  '  Active scalars    : ', cell_data%get_active_scalars()
   write(*, '(a)')    '  Validation        : PASSED'

contains

   subroutine check_ok(st, context)
      type(htg_status_t), intent(in) :: st
      character(len=*), intent(in) :: context
      if (.not. st%ok()) then
         error stop context // ': ' // trim(st%message)
      end if
   end subroutine check_ok

end program example_02_octree_amr_celldata
