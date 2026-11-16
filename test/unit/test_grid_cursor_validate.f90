program test_grid_cursor_validate
   use fvtkhtg
   implicit none

   call test_dimensions_and_coordinates()
   call test_root_indexing_modes()
   call test_sparse_forest_and_cursor_stability()
   call test_deep_cursor_amr_and_preconditions()
   call test_mask_pure_mask_and_interfaces()
   call test_pure_mask_honors_depth_limiter()
   call test_validation_corruption_suite()

   print '(a)', 'ALL test_grid_cursor_validate TESTS PASSED'

contains

   subroutine check(cond, msg)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) then
         error stop 'FAIL: ' // msg
      end if
   end subroutine check

   subroutine build_linspace(n, x0, dx, arr)
      integer(i32_k), intent(in) :: n
      real(r64_k), intent(in) :: x0, dx
      real(r64_k), allocatable, intent(out) :: arr(:)
      integer(i32_k) :: i

      allocate(arr(0 : n - 1_i32_k))
      do i = 0_i32_k, n - 1_i32_k
         arr(i) = x0 + real(i, r64_k) * dx
      end do
   end subroutine build_linspace

   ! ---------------------------------------------------------------------------
   ! 1. 1D (X, Y, Z), 2D (XY, YZ, XZ), and 3D rectilinear grids & coordinates
   ! ---------------------------------------------------------------------------
   subroutine test_dimensions_and_coordinates()
      type(htg_grid_t), target :: g
      type(htg_status_t) :: st
      real(r64_k), allocatable :: cx(:), cy(:), cz(:)
      real(r64_k), pointer :: vx(:), vy(:), vz(:)
      type(htg_tree_t), pointer :: t_null
      type(htg_bit_array_t), pointer :: pm0
      integer(i32_k) :: d(0:2), cd(0:2), ax(0:2), a1, a2

      ! 0D degenerate (1, 1, 1) -> dimension=0, cell_dims=[1,1,1], orient=-1, axes=[-1,-1,-1]
      call g%clear()
      call g%set_dimensions_ijk(1_i32_k, 1_i32_k, 1_i32_k, st)
      call check(st%ok(), '0D set_dimensions(1,1,1)')
      call check(g%get_dimension() == 0_i32_k, '0D dimension == 0')
      call check(g%get_orientation() == HTG_INVALID_ORIENTATION, '0D orientation == -1')
      ax = g%get_axes()
      call check(all(ax == [-1_i32_k, -1_i32_k, -1_i32_k]), '0D axes == [-1,-1,-1]')
      cd = g%get_cell_dims()
      call check(all(cd == [1_i32_k, 1_i32_k, 1_i32_k]), '0D cell_dims == [1,1,1]')
      call g%get_pure_mask(pm0, st)
      call check(st%ok(), '0D empty-grid PureMask status')
      call check(pm0%get_number_of_values() == 1_id_k, 'empty-grid PureMask has 1 bit')
      call check(.not. pm0%get_value(0_id_k), 'empty-grid PureMask(0) == false')
      call g%get_tree(0_id_k, t_null, create=.true., status=st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INVALID_STATE, &
         'reject tree creation on 0D (1,1,1) grid')
      call g%validate(st)
      call check(st%ok(), '0D empty grid validate: ' // trim(st%message))

      ! 1D X (5, 1, 1), branch factor 2 -> nchild = 2
      call g%clear()
      call g%set_branch_factor(2_i32_k, st)
      call check(st%ok(), '1D-X set_branch_factor(2)')
      call g%set_dimensions_ijk(5_i32_k, 1_i32_k, 1_i32_k, st)
      call check(st%ok(), '1D-X set_dimensions')
      call check(g%get_dimension() == 1_i32_k, '1D-X dimension == 1')
      call check(g%get_orientation() == HTG_AXIS_X, '1D-X orientation == X')
      call g%get_1d_axis(a1, st)
      call check(st%ok() .and. a1 == HTG_AXIS_X, '1D-X get_1d_axis == X')
      ax = g%get_axes()
      call check(ax(0) == HTG_AXIS_X, '1D-X axes(0) == X')
      cd = g%get_cell_dims()
      call check(all(cd == [4_i32_k, 1_i32_k, 1_i32_k]), '1D-X cell_dims')
      call check(g%get_number_of_children() == 2_i32_k, '1D-X nchildren == 2')
      call build_linspace(5_i32_k, 0.0_r64_k, 0.25_r64_k, cx)
      call g%set_x_coordinates(cx, st)
      call check(st%ok(), '1D-X set_x_coordinates')
      call g%validate(st)
      call check(st%ok(), '1D-X validate: ' // trim(st%message))

      ! 1D Y (1, 4, 1), branch factor 3 -> nchild = 3, descending coordinates
      call g%clear()
      call g%set_branch_factor(3_i32_k, st)
      call g%set_dimensions([1_i32_k, 4_i32_k, 1_i32_k], st)
      call check(st%ok(), '1D-Y set_dimensions')
      call check(g%get_dimension() == 1_i32_k, '1D-Y dimension == 1')
      call check(g%get_orientation() == HTG_AXIS_Y, '1D-Y orientation == Y')
      ax = g%get_axes()
      call check(ax(0) == HTG_AXIS_Y, '1D-Y axes(0) == Y')
      call check(g%get_number_of_children() == 3_i32_k, '1D-Y nchildren == 3')
      call build_linspace(4_i32_k, 10.0_r64_k, -2.0_r64_k, cy)
      call g%set_y_coordinates(cy, st)
      call check(st%ok(), '1D-Y descending coords accepted')
      call g%validate(st)
      call check(st%ok(), '1D-Y validate')

      ! 1D Z (1, 1, 6), branch factor 2 -> nchild = 2
      call g%clear()
      call g%set_dimensions(1_i32_k, 1_i32_k, 6_i32_k, st)
      call check(st%ok(), '1D-Z set_dimensions')
      call check(g%get_dimension() == 1_i32_k, '1D-Z dimension == 1')
      call check(g%get_orientation() == HTG_AXIS_Z, '1D-Z orientation == Z')
      call build_linspace(6_i32_k, -1.0_r64_k, 0.5_r64_k, cz)
      call g%set_z_coordinates(cz, st)
      call check(st%ok(), '1D-Z set_z_coordinates')
      call g%validate(st)
      call check(st%ok(), '1D-Z validate')

      ! 2D XY (4, 3, 1), f=2 -> nchild=4, orientation=2 (Z normal)
      call g%clear()
      call g%set_dimensions(4_i32_k, 3_i32_k, 1_i32_k, st)
      call check(g%get_dimension() == 2_i32_k, '2D-XY dimension == 2')
      call check(g%get_orientation() == HTG_AXIS_Z, '2D-XY orientation == Z')
      ax = g%get_axes()
      call check(ax(0) == HTG_AXIS_X .and. ax(1) == HTG_AXIS_Y, '2D-XY axes')
      call check(g%get_number_of_children() == 4_i32_k, '2D-XY nchildren == 4')
      call build_linspace(4_i32_k, 0.0_r64_k, 1.0_r64_k, cx)
      call build_linspace(3_i32_k, 0.0_r64_k, 2.0_r64_k, cy)
      call g%set_x_coordinates(cx, st)
      call g%set_y_coordinates(cy, st)
      call g%validate(st)
      call check(st%ok(), '2D-XY validate')

      ! 2D YZ (1, 4, 5), f=3 -> nchild=9, orientation=0 (X normal)
      call g%clear()
      call g%set_branch_factor(3_i32_k, st)
      call g%set_dimensions(1_i32_k, 4_i32_k, 5_i32_k, st)
      call check(g%get_dimension() == 2_i32_k, '2D-YZ dimension == 2')
      call check(g%get_orientation() == HTG_AXIS_X, '2D-YZ orientation == X')
      ax = g%get_axes()
      call check(ax(0) == HTG_AXIS_Y .and. ax(1) == HTG_AXIS_Z, '2D-YZ axes')
      call check(g%get_number_of_children() == 9_i32_k, '2D-YZ nchildren == 9')
      call build_linspace(4_i32_k, 0.0_r64_k, 1.0_r64_k, cy)
      call build_linspace(5_i32_k, 0.0_r64_k, 1.0_r64_k, cz)
      call g%set_y_coordinates(cy, st)
      call g%set_z_coordinates(cz, st)
      call g%validate(st)
      call check(st%ok(), '2D-YZ validate')

      ! 2D XZ (3, 1, 4), f=2 -> nchild=4, orientation=1 (Y normal), cyclic axes=(2,0)
      call g%clear()
      call g%set_branch_factor(2_i32_k, st)
      call g%set_dimensions(3_i32_k, 1_i32_k, 4_i32_k, st)
      call check(g%get_dimension() == 2_i32_k, '2D-XZ dimension == 2')
      call check(g%get_orientation() == HTG_AXIS_Y, '2D-XZ orientation == Y')
      ax = g%get_axes()
      call check(ax(0) == HTG_AXIS_Z .and. ax(1) == HTG_AXIS_X, '2D-XZ cyclic axes (2,0)')
      call g%get_2d_axes(a1, a2, st)
      call check(st%ok() .and. a1 == HTG_AXIS_Z .and. a2 == HTG_AXIS_X, &
         '2D-XZ get_2d_axes == (Z, X)')
      call build_linspace(3_i32_k, 0.0_r64_k, 1.0_r64_k, cx)
      call build_linspace(4_i32_k, 0.0_r64_k, 1.0_r64_k, cz)
      call g%set_x_coordinates(cx, st)
      call g%set_z_coordinates(cz, st)
      call g%validate(st)
      call check(st%ok(), '2D-XZ validate')

      ! 3D XYZ (3, 4, 5), f=3 -> nchild=27
      call g%clear()
      call g%set_branch_factor(3_i32_k, st)
      call g%set_dimensions(3_i32_k, 4_i32_k, 5_i32_k, st)
      d = g%get_dimensions()
      cd = g%get_cell_dims()
      call check(all(d == [3_i32_k, 4_i32_k, 5_i32_k]), '3D dimensions')
      call check(all(cd == [2_i32_k, 3_i32_k, 4_i32_k]), '3D cell_dims')
      call check(g%get_dimension() == 3_i32_k, '3D dimension == 3')
      call check(g%get_number_of_children() == 27_i32_k, '3D nchildren == 27')
      call check(g%get_max_number_of_trees() == 24_id_k, '3D max trees == 24')
      call build_linspace(3_i32_k, 0.0_r64_k, 1.0_r64_k, cx)
      call build_linspace(4_i32_k, 1.0_r64_k, 0.5_r64_k, cy)
      call build_linspace(5_i32_k, -2.0_r64_k, 1.0_r64_k, cz)
      call g%set_x_coordinates(cx, st)
      call g%set_y_coordinates(cy, st)
      call g%set_z_coordinates(cz, st)
      call g%get_x_coordinates_view(vx)
      call g%get_y_coordinates_view(vy)
      call g%get_z_coordinates_view(vz)
      call check(associated(vx) .and. lbound(vx, 1) == 0 .and. ubound(vx, 1) == 2, 'vx 0-based bounds')
      call check(associated(vy) .and. lbound(vy, 1) == 0 .and. ubound(vy, 1) == 3, 'vy 0-based bounds')
      call check(associated(vz) .and. lbound(vz, 1) == 0 .and. ubound(vz, 1) == 4, 'vz 0-based bounds')
      call g%validate(st)
      call check(st%ok(), '3D validate')

      ! Negative guard tests on dimensions, branch factor, and coordinates
      call g%set_branch_factor(4_i32_k, st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INVALID_ARG, &
         'reject branch_factor=4')
      call g%set_dimensions(0_i32_k, 3_i32_k, 3_i32_k, st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INVALID_ARG, &
         'reject zero dimension')
      call g%set_x_coordinates([0.0_r64_k, 1.0_r64_k], st)
      call check(.not. st%ok(), 'reject mismatched x_coordinates length')
      call g%set_x_coordinates([0.0_r64_k, 1.0_r64_k, 1.0_r64_k], st)
      call check(.not. st%ok(), 'reject duplicate x_coordinates')
      call g%set_x_coordinates([0.0_r64_k, 2.0_r64_k, 1.5_r64_k], st)
      call check(.not. st%ok(), 'reject non-monotonic x_coordinates')
   end subroutine test_dimensions_and_coordinates

   ! ---------------------------------------------------------------------------
   ! 2. Root indexing permutations (TransposedRootIndexing .false. vs .true.)
   ! ---------------------------------------------------------------------------
   subroutine test_root_indexing_modes()
      type(htg_grid_t) :: g
      type(htg_status_t) :: st
      integer(i32_k) :: i, j, k, io, jo, ko
      integer(id_k) :: tid, expected_tid

      call g%set_dimensions(4_i32_k, 3_i32_k, 5_i32_k, st)
      call check(st%ok(), 'root indexing set_dimensions (cell_dims = 3,2,4)')

      ! Standard mode (.false.): tid = i + j*3 + k*3*2
      call g%set_transposed_root_indexing(.false., st)
      call check(.not. g%get_transposed_root_indexing(), 'transposed is false')
      do k = 0_i32_k, 3_i32_k
         do j = 0_i32_k, 1_i32_k
            do i = 0_i32_k, 2_i32_k
               expected_tid = int(i + j * 3_i32_k + k * 6_i32_k, id_k)
               call g%get_index_from_level_zero_coordinates(tid, i, j, k, st)
               call check(st%ok() .and. tid == expected_tid, 'std l0 -> index')
               call g%get_level_zero_coordinates_from_index(tid, io, jo, ko, st)
               call check(st%ok() .and. io == i .and. jo == j .and. ko == k, 'std index -> l0')
            end do
         end do
      end do

      ! Transposed mode (.true.): tid = k + j*4 + i*4*2
      call g%set_transposed_root_indexing(.true., st)
      call check(g%get_transposed_root_indexing(), 'transposed is true')
      do i = 0_i32_k, 2_i32_k
         do j = 0_i32_k, 1_i32_k
            do k = 0_i32_k, 3_i32_k
               expected_tid = int(k + j * 4_i32_k + i * 8_i32_k, id_k)
               call g%get_index_from_level_zero_coordinates(tid, i, j, k, st)
               call check(st%ok() .and. tid == expected_tid, 'transposed l0 -> index')
               call g%get_level_zero_coordinates_from_index(tid, io, jo, ko, st)
               call check(st%ok() .and. io == i .and. jo == j .and. ko == k, 'transposed index -> l0')
            end do
         end do
      end do

      ! Out-of-bounds checks
      call g%get_index_from_level_zero_coordinates(tid, 3_i32_k, 0_i32_k, 0_i32_k, st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_OUT_OF_BOUNDS, 'l0 i out of bounds')
      call g%get_level_zero_coordinates_from_index(24_id_k, io, jo, ko, st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_OUT_OF_BOUNDS, 'l0 tid out of bounds')
   end subroutine test_root_indexing_modes

   ! ---------------------------------------------------------------------------
   ! 3. Sparse 500^3 (125M-cell) root grid with 5+ trees & cursor stability
   ! ---------------------------------------------------------------------------
   subroutine test_sparse_forest_and_cursor_stability()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: cur_anchor, cur_empty, cur_tmp
      type(htg_tree_t), pointer :: t_ptr, t_anchor
      type(htg_status_t) :: st
      real(r64_k), allocatable :: c501(:)
      integer(id_k), parameter :: IDS_5(5) = [ &
         50000000_id_k, 17_id_k, 124999999_id_k, 0_id_k, 25000000_id_k]
      integer(id_k), parameter :: SORTED_5(5) = [ &
         0_id_k, 17_id_k, 25000000_id_k, 50000000_id_k, 124999999_id_k]
      integer(id_k) :: pos, tid, extra_id
      integer(i32_k) :: idx

      call build_linspace(501_i32_k, 0.0_r64_k, 1.0_r64_k, c501)
      call g%set_dimensions(501_i32_k, 501_i32_k, 501_i32_k, st)
      call check(st%ok(), '500^3 set_dimensions')
      call check(g%get_max_number_of_trees() == 125000000_id_k, '500^3 max_number_of_trees == 125M')
      call g%set_x_coordinates(c501, st)
      call g%set_y_coordinates(c501, st)
      call g%set_z_coordinates(c501, st)

      ! Attach anchor cursor to tree 50,000,000 FIRST
      call htg_new_non_oriented_cursor(g, IDS_5(1), cur_anchor, create=.true., status=st)
      call check(st%ok() .and. cur_anchor%has_tree(), 'anchor cursor created')
      call cur_anchor%get_tree(t_anchor)
      call cur_anchor%subdivide_leaf(st)
      call check(st%ok(), 'anchor root subdivide')

      ! Insert remaining 4 sparse trees out of order
      do idx = 2, 5
         call htg_new_non_oriented_cursor(g, IDS_5(idx), cur_tmp, create=.true., status=st)
         call check(st%ok(), 'insert sparse tree')
      end do
      call check(g%get_number_of_non_empty_trees() == 5_id_k, '5 active trees in 125M grid')

      ! Verify deterministic ascending iteration over the 5 trees
      do pos = 0_id_k, 4_id_k
         call g%get_tree_by_pos(pos, tid, t_ptr, st)
         call check(st%ok() .and. tid == SORTED_5(pos + 1), 'deterministic ascending tree_id order')
      end do

      ! Query an uninstantiated root with create=.false.
      call new_non_oriented_cursor(g, 999999_id_k, cur_empty, create=.false., status=st)
      call check(st%ok() .and. .not. cur_empty%has_tree(), 'empty root cursor has_tree == false')
      call check(cur_empty%is_leaf(st), 'empty root cursor is_leaf == true')
      call check(cur_empty%get_vertex_id() == HTG_INVALID_ID, 'empty root vertex_id == -1')

      ! Insert 20 additional trees both before and after 50,000,000 to force multiple
      ! slot shifts and capacity reallocations while cur_anchor remains live!
      do idx = 1, 20
         extra_id = int(idx * 1000, id_k)
         call new_non_oriented_cursor(g, extra_id, cur_tmp, create=.true., status=st)
         call check(st%ok(), 'insert extra tree forcing slot growth/shift')
      end do
      call check(g%get_number_of_non_empty_trees() == 25_id_k, '25 active trees')

      ! Verify cur_anchor and t_anchor did NOT dangle and still mutate the exact tree in g!
      call g%get_tree(50000000_id_k, t_ptr, create=.false., status=st)
      call check(associated(t_ptr, t_anchor), 'heap tree address preserved across reallocations')
      call cur_anchor%to_child(3_i32_k, st)
      call check(st%ok() .and. cur_anchor%get_level() == 1_i32_k, 'cur_anchor to_child(3) after reallocs')
      call cur_anchor%subdivide_leaf(st)
      call check(st%ok(), 'cur_anchor subdivide child 3 after reallocs')
      call check(t_ptr%get_number_of_vertices() == 17_id_k, 'tree 50M has 1 + 8 + 8 = 17 vertices')

      ! Remove the 20 extra trees, initialize implicit global indices, and validate
      do idx = 1, 20
         extra_id = int(idx * 1000, id_k)
         call g%remove_tree(extra_id, st)
         call check(st%ok(), 'remove extra tree')
      end do
      call check(g%get_number_of_non_empty_trees() == 5_id_k, 'back to 5 active trees')
      call g%initialize_local_index_node(st)
      call check(st%ok(), 'initialize_local_index_node on 5 sparse trees')
      call check(g%get_number_of_vertices() == 21_id_k, 'total vertices = 17 + 4*1 = 21')
      call check(g%get_number_of_leaves() == 19_id_k, 'total leaves = 15 + 4*1 = 19')
      call check(g%get_number_of_levels() == 3_id_k, 'max levels == 3')

      call g%validate(st)
      call check(st%ok(), '500^3 sparse grid validate: ' // trim(st%message))
   end subroutine test_sparse_forest_and_cursor_stability

   ! ---------------------------------------------------------------------------
   ! 4. Deep cursor AMR construction (12+ levels) & cursor guard checks
   ! ---------------------------------------------------------------------------
   subroutine test_deep_cursor_amr_and_preconditions()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: cur
      type(htg_tree_t), pointer :: t_ptr
      type(htg_status_t) :: st
      real(r64_k), allocatable :: cx(:), cy(:)
      integer(i32_k) :: depth
      integer(id_k) :: nv_before

      call g%set_dimensions(3_i32_k, 3_i32_k, 1_i32_k, st)
      call build_linspace(3_i32_k, 0.0_r64_k, 1.0_r64_k, cx)
      call build_linspace(3_i32_k, 0.0_r64_k, 1.0_r64_k, cy)
      call g%set_x_coordinates(cx, st)
      call g%set_y_coordinates(cy, st)

      call new_non_oriented_cursor(g, 0_id_k, cur, create=.true., status=st)
      call cur%set_global_index_start(0_id_k, st)
      call check(cur%is_root(), 'initial cursor at root')
      call check(cur%is_leaf(), 'initial cursor at leaf')

      ! Precondition checks at root leaf
      call cur%to_parent(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INVALID_STATE, &
         'reject to_parent at root')
      call cur%to_child(0_i32_k, st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_TOPOLOGY, &
         'reject to_child on leaf')

      ! Subdivide 12 levels deep (expands cursor entry stack past INITIAL_ENTRY_CAPACITY=8)
      do depth = 0_i32_k, 11_i32_k
         call check(cur%get_level() == depth, 'cursor depth check')
         call cur%subdivide_leaf(st)
         call check(st%ok(), 'deep subdivide')
         call check(.not. cur%is_leaf(), 'newly split node is not leaf')
         call check(cur%is_terminal_node(), 'newly split node is terminal')
         call cur%to_child(mod(depth, 4_i32_k), st)
         call check(st%ok(), 'descend to child')
      end do
      call check(cur%get_level() == 12_i32_k, 'reached level 12')
      call check(g%get_number_of_levels() == 13_id_k, '13 levels total (0..12)')

      ! Walk back up one level: parent of level-12 leaves is terminal
      call cur%to_parent(st)
      call check(st%ok() .and. cur%get_level() == 11_i32_k, 'stepped to level 11')
      call check(cur%is_terminal_node(), 'level 11 node is terminal')

      ! VTK v9.7.2: cursor%subdivide_leaf on coarse node is a silent no-op,
      ! whereas direct tree%subdivide_leaf returns HTG_ERR_TOPOLOGY!
      nv_before = g%get_number_of_vertices()
      call cur%subdivide_leaf(st)
      call check(st%ok(), 'cursor%subdivide_leaf on coarse node is no-op')
      call check(g%get_number_of_vertices() == nv_before, &
         'coarse no-op left vertex count unchanged')
      call cur%get_tree(t_ptr)
      call t_ptr%subdivide_leaf(cur%get_vertex_id(), cur%get_level(), st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_TOPOLOGY, &
         'direct tree%subdivide_leaf rejects coarse node')

      call cur%to_child(4_i32_k, st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_OUT_OF_BOUNDS, &
         'reject to_child(4) on quadtree')

      ! Masked coarse node blocks ToChild (VTK v9.7.2 assert(!IsMasked))
      call cur%set_mask(.true., st)
      call check(st%ok(), 'mask coarse node at level 11')
      call cur%to_child(0_i32_k, st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_MASK, &
         'reject to_child on masked coarse node')
      call cur%set_mask(.false., st)
      call g%clear_mask()

      call cur%to_parent(st)
      call check(st%ok() .and. cur%get_level() == 10_i32_k, 'stepped to level 10')
      call check(.not. cur%is_terminal_node(), 'level 10 node is not terminal')

      ! Lowering depth_limiter below current cursor level (level=10 > depth_limiter=4):
      ! VTK v9.7.2: is_leaf returns true at >= depth_limiter, and both to_child and
      ! subdivide_leaf reject level > depth_limiter with HTG_ERR_INVALID_STATE!
      call g%set_depth_limiter(4_i32_k, st)
      call check(st%ok(), 'set_depth_limiter(4) while cursor at level 10')
      call check(cur%is_leaf(), 'cursor at level > depth_limiter reports is_leaf == true')
      call cur%to_child(0_i32_k, st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INVALID_STATE, &
         'to_child rejects cursor level > depth_limiter')
      call cur%subdivide_leaf(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INVALID_STATE, &
         'subdivide_leaf rejects cursor level > depth_limiter')
      call g%set_depth_limiter(huge(0_i32_k), st)

      ! DepthLimiter clamping check at exact equality (level == depth_limiter)
      call cur%to_root(st)
      call check(st%ok() .and. cur%is_root() .and. cur%get_level() == 0_i32_k, &
         'to_root resets to level 0')
      call g%set_depth_limiter(1_i32_k, st)
      call cur%to_child(0_i32_k, st)
      call check(st%ok() .and. cur%get_level() == 1_i32_k, 'stepped to depth_limiter=1')
      call check(cur%is_leaf(), 'cursor treats node at depth_limiter as leaf')
      call cur%to_child(1_i32_k, st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_TOPOLOGY, &
         'to_child blocked at depth_limiter')
      call g%set_depth_limiter(huge(0_i32_k), st)

      call g%validate(st)
      call check(st%ok(), 'deep AMR grid validate')
   end subroutine test_deep_cursor_amr_and_preconditions

   ! ---------------------------------------------------------------------------
   ! 5. Mask, PureMask invalidation/caching, CellData growth & Material Interfaces
   ! ---------------------------------------------------------------------------
   subroutine test_mask_pure_mask_and_interfaces()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: cur
      type(htg_cell_data_t), pointer :: cd
      type(htg_field_array_t), pointer :: arr_temp, arr_norm, arr_inter
      type(htg_bit_array_t) :: empty_mask
      type(htg_bit_array_t), pointer :: pm, m_ptr
      type(htg_status_t) :: st
      real(r64_k), allocatable :: cx(:), cy(:)
      integer(id_k) :: gid_root, gid_c0, gid_c1, gid_c2, gid_c3

      call g%set_dimensions(2_i32_k, 2_i32_k, 1_i32_k, st)
      call build_linspace(2_i32_k, 0.0_r64_k, 1.0_r64_k, cx)
      call build_linspace(2_i32_k, 0.0_r64_k, 1.0_r64_k, cy)
      call g%set_x_coordinates(cx, st)
      call g%set_y_coordinates(cy, st)

      ! VTK v9.7.2 HasMask() parity on 0-tuple Mask attached before global indexing:
      ! get_mask() is associated, has_mask() is false, and cursor subdivide/to_child
      ! work without global indices assigned!
      call empty_mask%initialize(0_id_k, .false., st)
      call g%set_mask(empty_mask, st)
      call check(st%ok(), 'set_mask(empty_mask) ok')
      call check(.not. g%has_mask(), 'has_mask() is false when mask has 0 tuples')
      call g%get_mask(m_ptr)
      call check(associated(m_ptr), 'get_mask() still returns attached empty mask')
      call new_non_oriented_cursor(g, 0_id_k, cur, create=.true., status=st)
      call cur%subdivide_leaf(st)
      call check(st%ok(), 'subdivide_leaf succeeds with 0-tuple mask before global index')
      call cur%to_child(0_i32_k, st)
      call check(st%ok(), 'to_child succeeds with 0-tuple mask before global index')
      call htg_validate_grid(g, st, allow_unset_global_index=.true.)
      call check(st%ok(), 'validate allows 0-tuple mask before global index')
      call g%remove_tree(0_id_k, st)
      call g%clear_mask()

      call g%get_cell_data(cd)
      call cd%add_array('Temperature', HTG_TYPE_REAL64, 1_i32_k, 0_id_k, st)
      call check(st%ok(), 'add Temperature array')
      call cd%add_array('Normals', HTG_TYPE_REAL64, 3_i32_k, 0_id_k, st)
      call check(st%ok(), 'add Normals array')
      call cd%add_array('Intercepts', HTG_TYPE_REAL64, 3_i32_k, 0_id_k, st)
      call check(st%ok(), 'add Intercepts array')
      call cd%set_active_scalars('Temperature', st)
      call check(st%ok(), 'set_active_scalars')

      call new_non_oriented_cursor(g, 0_id_k, cur, create=.true., status=st)
      call cur%set_global_index_start(0_id_k, st)
      call check(cd%get_number_of_tuples() == 1_id_k, 'CellData grew to 1 tuple')

      ! Subdivide root -> 5 vertices total (0, 1, 2, 3, 4); CellData must auto-grow to 5 tuples!
      call cur%subdivide_leaf(st)
      call check(cd%get_number_of_tuples() == 5_id_k, 'CellData auto-grew to 5 tuples on subdivide')

      gid_root = cur%get_global_node_index(st)
      call cur%to_child(0_i32_k, st)
      gid_c0 = cur%get_global_node_index(st)
      call cur%to_parent(st)
      call cur%to_child(1_i32_k, st)
      gid_c1 = cur%get_global_node_index(st)
      ! Subdivide child 1 -> vertices 5,6,7,8 (9 vertices total)
      call cur%subdivide_leaf(st)
      call check(cd%get_number_of_tuples() == 9_id_k, 'CellData auto-grew to 9 tuples')
      call cur%to_parent(st)
      call cur%to_child(2_i32_k, st)
      gid_c2 = cur%get_global_node_index(st)
      call cur%to_parent(st)
      call cur%to_child(3_i32_k, st)
      gid_c3 = cur%get_global_node_index(st)

      ! Initially no mask & no interface -> all PureMask bits should be .false.
      call g%get_pure_mask(pm, st)
      call check(st%ok() .and. g%is_pure_mask_valid(), 'initial PureMask valid')
      call check(.not. pm%get_value(gid_root), 'initial root is pure (bit=false)')

      ! Mask leaf child 3 (gid_c3): invalidates PureMask!
      call cur%set_mask(.true., st)
      call check(st%ok(), 'child 3 set_mask ok')
      call check(cur%is_masked(), 'child 3 masked')
      call check(.not. g%is_pure_mask_valid(), 'set_mask invalidated PureMask')

      ! Cannot subdivide a masked leaf!
      call cur%subdivide_leaf(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_MASK, &
         'reject subdividing masked leaf')

      ! Recompute PureMask: root (gid 0) and child 3 (gid 4) are impure (.true.);
      ! child 0 (1), coarse child 1 (2) + leaves (5..8), and child 2 (3) are pure (.false.)
      call g%get_pure_mask(pm, st)
      call check(st%ok() .and. g%is_pure_mask_valid(), 'recomputed PureMask valid')
      call check(pm%get_value(gid_root), 'root impure due to masked child 3')
      call check(pm%get_value(gid_c3), 'masked child 3 impure')
      call check(.not. pm%get_value(gid_c0), 'unmasked leaf child 0 pure')
      call check(.not. pm%get_value(gid_c1), 'coarse child 1 with pure children is pure')
      call check(.not. pm%get_value(gid_c2), 'unmasked leaf child 2 pure')

      ! VTK v9.7.2: Enable Material Interfaces with Intercepts[gid][2] < 2.0 rule
      ! (-1, 0, 1 = mixed interface; 2.0 = pure cell)
      call cd%get_array('Normals', arr_norm, st)
      call cd%get_array('Intercepts', arr_inter, st)
      call cd%get_array('Temperature', arr_temp, st)
      call arr_temp%set_scalar_r64(5_id_k, 300.0_r64_k, st)
      call arr_norm%set_tuple_r64(5_id_k, [1.0_r64_k, 0.0_r64_k, 0.0_r64_k], st)

      block
         integer(id_k) :: idx_t
         do idx_t = 0_id_k, 8_id_k
            call arr_inter%set_tuple_r64(idx_t, [0.0_r64_k, 0.0_r64_k, 2.0_r64_k], st)
         end do
      end block
      call arr_inter%set_tuple_r64(5_id_k, [-0.1_r64_k, 0.1_r64_k, 1.0_r64_k], st)

      call g%set_interface_normals_name('Normals')
      call g%set_interface_intercepts_name('Intercepts')
      call g%set_has_interface(.true.)
      call check(.not. g%is_pure_mask_valid(), 'set_has_interface invalidated PureMask')

      call g%get_pure_mask(pm, st)
      call check(st%ok(), 'PureMask with material interface')
      call check(pm%get_value(5_id_k), 'interface leaf 5 (<2.0) marked impure in PureMask')
      call check(pm%get_value(gid_c1), 'coarse gid_c1 marked impure via interface child 5')
      call check(.not. pm%get_value(6_id_k), 'pure sibling leaf 6 (2.0) remains false')

      ! Masked coarse node short-circuit test (VTK v9.7.2 RecursivelyInitializePureMask):
      ! Masking coarse gid_c1 marks gid_c1 = true WITHOUT visiting child 5, so child 5's
      ! PureMask bit stays at its initial 0 (false)!
      call g%set_cell_mask(gid_c1, .true., st)
      call g%get_pure_mask(pm, st)
      call check(st%ok(), 'PureMask with masked coarse gid_c1')
      call check(pm%get_value(gid_c1), 'masked coarse gid_c1 is true in PureMask')
      call check(.not. pm%get_value(5_id_k), &
         'unvisited child 5 under masked coarse gid_c1 stays false in PureMask')
      call g%set_cell_mask(gid_c1, .false., st)
      call g%get_pure_mask(pm, st)
      call check(pm%get_value(5_id_k), 'unmasked coarse gid_c1 visits child 5 again')

      call g%validate(st)
      call check(st%ok(), &
         'grid with mask + PureMask + interfaces validates: ' // trim(st%message))
   end subroutine test_mask_pure_mask_and_interfaces

   ! ---------------------------------------------------------------------------
   ! 5b. PureMask stops at the depth limiter (VTK v9.7.2 cursor IsLeaf semantics).
   !     Expected bits verified against C++ vtkHyperTreeGrid::GetPureMask().
   ! ---------------------------------------------------------------------------
   subroutine test_pure_mask_honors_depth_limiter()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_bit_array_t), pointer :: pm
      type(htg_status_t) :: st
      real(r64_k), allocatable :: cx(:)
      integer(id_k) :: i

      call g%set_dimensions(2_i32_k, 2_i32_k, 1_i32_k, st)
      call build_linspace(2_i32_k, 0.0_r64_k, 1.0_r64_k, cx)
      call g%set_x_coordinates(cx, st)
      call g%set_y_coordinates(cx, st)

      ! Root (gid 0) -> children gids 1..4; child 0 (gid 1) -> grandchildren gids 5..8
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=st)
      call c%subdivide_leaf(st)
      call c%to_child(0_i32_k, st)
      call c%subdivide_leaf(st)
      call g%initialize_local_index_node(st)
      call g%set_cell_mask(5_id_k, .true., st)
      call g%set_cell_mask(8_id_k, .false., st)

      ! Unlimited: masked grandchild 5 propagates to gid 1 and the root
      call g%get_pure_mask(pm, st)
      call check(st%ok() .and. pm%get_number_of_values() == 9_id_k, 'unlimited PureMask size')
      do i = 0_id_k, 8_id_k
         call check(pm%get_value(i) .eqv. (i == 0_id_k .or. i == 1_id_k .or. i == 5_id_k), &
            'unlimited PureMask bits == 110001000')
      end do

      ! depth_limiter = 1: depth-1 cells are leaves, grandchildren are never visited
      call g%set_depth_limiter(1_i32_k, st)
      call check(.not. g%is_pure_mask_valid(), 'set_depth_limiter invalidated PureMask')
      call g%get_pure_mask(pm, st)
      call check(st%ok() .and. pm%get_number_of_values() == 9_id_k, 'limited PureMask size')
      do i = 0_id_k, 8_id_k
         call check(.not. pm%get_value(i), 'depth_limiter=1 PureMask bits == 000000000')
      end do
      call g%validate(st)
      call check(st%ok(), 'depth-limited PureMask grid validates: ' // trim(st%message))

      ! depth_limiter = 0: masked root is still reported
      call g%set_cell_mask(0_id_k, .true., st)
      call g%set_depth_limiter(0_i32_k, st)
      call g%get_pure_mask(pm, st)
      call check(st%ok(), 'depth_limiter=0 PureMask computed')
      call check(pm%get_value(0_id_k), 'depth_limiter=0 PureMask marks the masked root')
      call check(.not. pm%get_value(1_id_k), 'depth_limiter=0 PureMask leaves child 0 unset')
   end subroutine test_pure_mask_honors_depth_limiter

   ! ---------------------------------------------------------------------------
   ! 6. Deliberate internal corruption suite proving htg_validate_grid catches all violations
   ! ---------------------------------------------------------------------------
   subroutine make_baseline_valid_grid(g)
      type(htg_grid_t), target, intent(out) :: g
      type(htg_cursor_t) :: c0, c1
      type(htg_cell_data_t), pointer :: cd
      type(htg_bit_array_t), pointer :: pm
      type(htg_status_t) :: st
      real(r64_k), allocatable :: cx(:), cy(:)

      call g%set_dimensions(3_i32_k, 2_i32_k, 1_i32_k, st)
      call build_linspace(3_i32_k, 0.0_r64_k, 1.0_r64_k, cx)
      call build_linspace(2_i32_k, 0.0_r64_k, 1.0_r64_k, cy)
      call g%set_x_coordinates(cx, st)
      call g%set_y_coordinates(cy, st)

      ! Tree 0: implicit global indexing [0..8] (root + child 0 subdivided -> 9 vertices)
      call new_non_oriented_cursor(g, 0_id_k, c0, create=.true., status=st)
      call c0%set_global_index_start(0_id_k, st)
      call c0%subdivide_leaf(st)
      call c0%to_child(0_i32_k, st)
      call c0%subdivide_leaf(st)

      ! Tree 1: explicit global indexing [100, 101, 102, 103, 104] (5 vertices)
      call new_non_oriented_cursor(g, 1_id_k, c1, create=.true., status=st)
      call c1%subdivide_leaf(st)
      call c1%set_global_index_from_local(100_id_k, st)
      block
         integer(i32_k) :: ch
         do ch = 0_i32_k, 3_i32_k
            call c1%to_root(st)
            call c1%to_child(ch, st)
            call c1%set_global_index_from_local(101_id_k + int(ch, id_k), st)
         end do
      end block

      call g%get_cell_data(cd)
      call cd%add_array('Density', HTG_TYPE_REAL64, 1_i32_k, 105_id_k, st)
      call cd%add_array('Normals', HTG_TYPE_REAL64, 3_i32_k, 105_id_k, st)
      call cd%add_array('Intercepts', HTG_TYPE_REAL64, 3_i32_k, 105_id_k, st)
      call g%set_interface_normals_name('Normals')
      call g%set_interface_intercepts_name('Intercepts')
      call g%set_has_interface(.true.)
      call g%set_cell_mask(4_id_k, .true., st)
      call g%get_pure_mask(pm, st)

      call g%validate(st)
      call check(st%ok(), 'baseline grid must be valid: ' // trim(st%message))
   end subroutine make_baseline_valid_grid

   subroutine test_validation_corruption_suite()
      type(htg_grid_t), target :: g
      type(htg_status_t) :: st

      ! (a) Corrupted squeezed geometry (orientation mismatch)
      call make_baseline_valid_grid(g)
      g%orientation = 0_i32_k
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INVALID_STATE, 'catch corrupted orientation')

      ! (b) Corrupted non-monotonic coordinate injected directly into x_coordinates
      call make_baseline_valid_grid(g)
      g%x_coordinates(2) = 0.5_r64_k
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INVALID_ARG, 'catch corrupted non-monotonic coord')

      ! (c) Unsorted sparse forest tree_slots
      call make_baseline_valid_grid(g)
      g%tree_slots(0)%tree_id = 1_id_k
      g%tree_slots(0)%tree%tree_index = 1_id_k
      call g%validate(st)
      call check(.not. st%ok(), 'catch unsorted/duplicate tree_slots')

      ! (d) Violated k-ary tree counting identity (number_of_leaves off by +1)
      call make_baseline_valid_grid(g)
      g%tree_slots(0)%tree%number_of_leaves = g%tree_slots(0)%tree%number_of_leaves + 1_id_k
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_TOPOLOGY, 'catch violated leaf counting identity')

      ! (e) Violated per-depth sum / count mismatch
      call make_baseline_valid_grid(g)
      g%tree_slots(0)%tree%number_of_vertices_per_depth(1) = 3_id_k
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_TOPOLOGY, 'catch per-depth sum mismatch')

      ! (f) Violated topological bijection: two coarse parents sharing the same elder child block (1..4),
      ! leaving vertices 5..8 orphaned!
      call make_baseline_valid_grid(g)
      g%tree_slots(0)%tree%parent_to_elder_child(1) = 1_id_k
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_TOPOLOGY, 'catch cycle / shared child / orphan block')

      ! (g) Out-of-bounds elder child index in parent_to_elder_child
      call make_baseline_valid_grid(g)
      g%tree_slots(0)%tree%parent_to_elder_child(1) = 6_id_k
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_TOPOLOGY, 'catch out-of-bounds elder block')

      ! (h) Unassigned global index on active tree
      call make_baseline_valid_grid(g)
      g%tree_slots(1)%tree%global_index_table(3) = HTG_INVALID_ID
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INDEXING, 'catch incomplete explicit global index')

      ! (i) Duplicate explicit global ID within a tree
      call make_baseline_valid_grid(g)
      g%tree_slots(1)%tree%global_index_table(3) = 101_id_k
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INDEXING, 'catch duplicate explicit global ID')

      ! (j) Cross-mode collision: explicit global ID in Tree 1 colliding with implicit interval [0..8] of Tree 0
      call make_baseline_valid_grid(g)
      g%tree_slots(1)%tree%global_index_table(2) = 7_id_k
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INDEXING, 'catch explicit vs implicit global ID collision')

      ! (k) Overlapping implicit global index intervals across two implicit trees
      call make_baseline_valid_grid(g)
      g%tree_slots(1)%tree%indexing_mode = HTG_INDEX_IMPLICIT
      deallocate(g%tree_slots(1)%tree%global_index_table)
      g%tree_slots(1)%tree%global_index_start = 8_id_k
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INDEXING, 'catch overlapping implicit intervals [0..8] & [8..12]')

      ! (l) Stale / corrupted cached PureMask bit
      call make_baseline_valid_grid(g)
      call g%pure_mask%set_value(0_id_k, .false., st)
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_MASK, 'catch stale/corrupted cached pure_mask')

      ! (m) Non-zero trailing slack bit in packed Mask byte
      call make_baseline_valid_grid(g)
      g%mask%bytes(0) = ibset(g%mask%bytes(0), 0)
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_MASK, 'catch non-zero trailing slack bit in mask')

      ! (n) Insufficient CellData tuple count for max_global_index
      call make_baseline_valid_grid(g)
      call g%cell_data%set_number_of_tuples(50_id_k, st)
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_OUT_OF_BOUNDS, 'catch insufficient cell_data tuples')

      ! (o) Missing Material Interface normals array
      call make_baseline_valid_grid(g)
      call g%cell_data%remove_array('Normals', st)
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INTERFACE, 'catch missing interface normals array')

      ! (p) Wrong component count on Material Interface normals array (2 instead of 3)
      call make_baseline_valid_grid(g)
      call g%cell_data%remove_array('Normals', st)
      call g%cell_data%add_array('Normals', HTG_TYPE_REAL64, 2_i32_k, 105_id_k, st)
      call g%validate(st)
      call check(.not. st%ok() .and. st%code == HTG_ERR_INTERFACE, &
         'catch wrong component count on Normals')

      ! (q) Explicit global ID == huge(id_k) overflow guard in validate and get_pure_mask
      block
         type(htg_bit_array_t), pointer :: pm_ovf
         call make_baseline_valid_grid(g)
         g%tree_slots(1)%tree%global_index_table(4) = huge(0_id_k)
         call g%invalidate_pure_mask()
         call g%validate(st)
         call check(.not. st%ok() .and. st%code == HTG_ERR_OUT_OF_BOUNDS, &
            'catch huge(id_k) global ID in validate without int64 overflow')
         call g%get_pure_mask(pm_ovf, st)
         call check(.not. st%ok() .and. st%code == HTG_ERR_INDEXING, &
            'catch huge(id_k) global ID in get_pure_mask without int64 overflow')
      end block
   end subroutine test_validation_corruption_suite

end program test_grid_cursor_validate
