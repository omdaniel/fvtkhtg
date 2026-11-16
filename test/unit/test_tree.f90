program test_tree
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use fvtkhtg_kinds_m, only: i32_k, id_k, htg_status_t, &
      HTG_ERR_INDEXING, HTG_ERR_INVALID_ARG, HTG_ERR_INVALID_STATE, &
      HTG_ERR_MASK, HTG_ERR_OUT_OF_BOUNDS, HTG_ERR_TOPOLOGY, &
      HTG_INDEX_EXPLICIT, HTG_INDEX_IMPLICIT, HTG_INDEX_UNSET, &
      HTG_INVALID_ID
   use fvtkhtg_bit_array_m, only: htg_bit_array_t
   use fvtkhtg_tree_m, only: htg_tree_t, HTG_TREE_LEAF_SENTINEL
   implicit none

   call test_uninitialized_guards()
   call test_all_branch_dimension_combinations()
   call test_deep_asymmetric_subdivision_and_terminal_nodes()
   call test_global_indexing_state_machine()
   call test_bf_descriptor_single_leaf_and_balanced()
   call test_bf_descriptor_asymmetric_and_forest_concat()
   call test_bf_descriptor_mask_and_depth_limiter()
   call test_bf_descriptor_malformed_rejection()
   call test_vtk_9_7_2_audit_exact_vectors()

   write (output_unit, '(a)') 'PASS: all fvtkhtg_tree_m unit tests succeeded.'

contains

   subroutine check(cond, msg)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) then
         write (error_unit, '(a,a)') 'FAIL: ', msg
         error stop 1
      end if
   end subroutine check

   subroutine test_uninitialized_guards()
      type(htg_tree_t), target :: tree
      type(htg_status_t) :: st
      type(htg_bit_array_t) :: desc
      integer(id_k), pointer :: p_elder(:), p_depth(:), p_gid(:)
      integer(id_k), allocatable :: nv_depth(:), bf_map(:)
      integer(id_k) :: val
      logical :: flag

      call check(.not. tree%is_initialized(), 'default tree must be uninitialized')
      call tree%get_parent_to_elder_child_view(p_elder)
      call tree%get_vertices_per_depth_view(p_depth)
      call tree%get_global_index_table_view(p_gid)
      call check(.not. associated(p_elder), 'uninit elder view null')
      call check(.not. associated(p_depth), 'uninit depth view null')
      call check(.not. associated(p_gid), 'uninit global view null')

      flag = tree%is_leaf(0_id_k, st)
      call check(.not. flag .and. st%code == HTG_ERR_INVALID_STATE, 'uninit is_leaf guard')

      flag = tree%is_terminal_node(0_id_k, st)
      call check(.not. flag .and. st%code == HTG_ERR_INVALID_STATE, 'uninit is_terminal_node guard')

      val = tree%get_elder_child_index(0_id_k, st)
      call check(val == HTG_TREE_LEAF_SENTINEL .and. st%code == HTG_ERR_INVALID_STATE, &
         'uninit get_elder_child_index guard')

      call tree%subdivide_leaf(0_id_k, 0_i32_k, st)
      call check(st%code == HTG_ERR_INVALID_STATE, 'uninit subdivide_leaf guard')

      call tree%set_global_index_start(0_id_k, st)
      call check(st%code == HTG_ERR_INVALID_STATE, 'uninit set_global_index_start guard')

      call tree%set_global_index_from_local(0_id_k, 0_id_k, st)
      call check(st%code == HTG_ERR_INVALID_STATE, 'uninit set_global_index_from_local guard')

      val = tree%get_global_index_from_local(0_id_k, st)
      call check(val == HTG_INVALID_ID .and. st%code == HTG_ERR_INVALID_STATE, &
         'uninit get_global_index_from_local guard')

      val = tree%get_max_global_index(st)
      call check(val == HTG_INVALID_ID .and. st%code == HTG_ERR_INVALID_STATE, &
         'uninit get_max_global_index guard')

      call desc%initialize(0_id_k)
      call tree%build_from_breadth_first_order_descriptor(desc, 0_id_k, status=st)
      call check(st%code == HTG_ERR_INVALID_STATE, 'uninit build_from_bf guard')

      call tree%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map, status=st)
      call check(st%code == HTG_ERR_INVALID_STATE, 'uninit compute_bf guard')

      ! Invalid initialize arguments
      call tree%initialize(1_i32_k, 3_i32_k, status=st)
      call check(st%code == HTG_ERR_INVALID_ARG .and. .not. tree%is_initialized(), 'bad bf=1')

      call tree%initialize(4_i32_k, 2_i32_k, status=st)
      call check(st%code == HTG_ERR_INVALID_ARG .and. .not. tree%is_initialized(), 'bad bf=4')

      call tree%initialize(2_i32_k, 0_i32_k, status=st)
      call check(st%code == HTG_ERR_INVALID_ARG .and. .not. tree%is_initialized(), 'bad dim=0')

      call tree%initialize(3_i32_k, 4_i32_k, status=st)
      call check(st%code == HTG_ERR_INVALID_ARG .and. .not. tree%is_initialized(), 'bad dim=4')

      call tree%initialize(2_i32_k, 3_i32_k, -2_id_k, st)
      call check(st%code == HTG_ERR_INVALID_ARG .and. .not. tree%is_initialized(), 'bad tree_id=-2')
   end subroutine test_uninitialized_guards

   subroutine test_all_branch_dimension_combinations()
      integer(i32_k), parameter :: bfs(6) = [2, 2, 2, 3, 3, 3]
      integer(i32_k), parameter :: dims(6) = [1, 2, 3, 1, 2, 3]
      integer(i32_k), parameter :: expected_k(6) = [2, 4, 8, 3, 9, 27]

      type(htg_tree_t), target :: tree
      type(htg_status_t) :: st
      integer(id_k), pointer :: p_elder(:), p_depth(:)
      integer(i32_k) :: icase, bf, d, k
      integer(id_k) :: last_child, elder

      do icase = 1, 6
         bf = bfs(icase)
         d = dims(icase)
         k = expected_k(icase)

         call tree%initialize(bf, d, int(icase - 1, id_k), st)
         call check(st%ok(), 'initialize ok for valid (bf, d)')
         call check(tree%is_initialized(), 'is_initialized true')
         call check(tree%get_tree_index() == int(icase - 1, id_k), 'tree_index match')
         call check(tree%get_branch_factor() == bf, 'branch_factor match')
         call check(tree%get_dimension() == d, 'dimension match')
         call check(tree%get_number_of_children() == k, 'number_of_children match')
         call check(tree%get_number_of_levels() == 1_id_k, 'root levels == 1')
         call check(tree%get_number_of_vertices() == 1_id_k, 'root vertices == 1')
         call check(tree%get_number_of_nodes() == 0_id_k, 'root nodes == 0')
         call check(tree%get_number_of_leaves() == 1_id_k, 'root leaves == 1')
         call check(tree%get_number_of_vertices_at_depth(0_i32_k) == 1_id_k, 'depth 0 == 1')
         call check(tree%is_leaf(0_id_k), 'root is leaf initially')
         call check(.not. tree%is_terminal_node(0_id_k), 'root leaf is not terminal coarse node')

         call tree%get_parent_to_elder_child_view(p_elder)
         call tree%get_vertices_per_depth_view(p_depth)
         call check(associated(p_elder), 'p_elder associated')
         call check(lbound(p_elder, 1) == 0 .and. ubound(p_elder, 1) == 0, 'p_elder bounds 0:0')
         call check(p_elder(0) == HTG_TREE_LEAF_SENTINEL, 'root elder sentinel')
         call check(lbound(p_depth, 1) == 0 .and. ubound(p_depth, 1) == 0, 'p_depth bounds 0:0')
         call check(p_depth(0) == 1_id_k, 'p_depth(0) == 1')

         ! Subdivide root -> k children at indices 1 .. k
         call tree%subdivide_leaf(0_id_k, 0_i32_k, st)
         call check(st%ok(), 'subdivide root ok')
         call check(.not. tree%is_leaf(0_id_k), 'root now coarse')
         call check(tree%is_terminal_node(0_id_k), 'root now terminal coarse node')
         elder = tree%get_elder_child_index(0_id_k, st)
         call check(st%ok() .and. elder == 1_id_k, 'root elder == 1')
         call check(tree%get_number_of_levels() == 2_id_k, 'levels == 2')
         call check(tree%get_number_of_nodes() == 1_id_k, 'nodes == 1')
         call check(tree%get_number_of_leaves() == int(k, id_k), 'leaves == k')
         call check(tree%get_number_of_vertices() == 1_id_k + int(k, id_k), 'vertices == 1 + k')
         call check(tree%get_number_of_vertices_at_depth(1_id_k) == int(k, id_k), 'depth 1 == k')

         ! Subdivide last child (index k) at level 1
         last_child = int(k, id_k)
         call tree%subdivide_leaf(last_child, 1_id_k, st)
         call check(st%ok(), 'subdivide last child ok')
         call check(.not. tree%is_terminal_node(0_id_k), 'root no longer terminal after child split')
         call check(tree%is_terminal_node(last_child), 'split child is now terminal node')
         call check(tree%get_number_of_levels() == 3_id_k, 'levels == 3')
         call check(tree%get_number_of_nodes() == 2_id_k, 'nodes == 2')
         call check(tree%get_number_of_leaves() == 1_id_k + 2_id_k * int(k - 1, id_k), &
            'k-ary leaf identity')
         call check(tree%get_number_of_vertices() == 1_id_k + 2_id_k * int(k, id_k), &
            'k-ary vertex identity')

         call tree%get_parent_to_elder_child_view(p_elder)
         call check(lbound(p_elder, 1) == 0 .and. ubound(p_elder, 1) == int(last_child), &
            'p_elder active_len == last_child + 1')
         call check(p_elder(0) == 1_id_k, 'p_elder(0) == 1')
         if (k > 1) then
            call check(all(p_elder(1 : last_child - 1_id_k) == HTG_TREE_LEAF_SENTINEL), &
               'intermediate leaf slots hold sentinel')
         end if
         call check(p_elder(last_child) == 1_id_k + int(k, id_k), 'p_elder(last_child) == 1 + k')

         call tree%clear()
         call check(.not. tree%is_initialized(), 'clear resets initialized')
      end do
   end subroutine test_all_branch_dimension_combinations

   subroutine test_deep_asymmetric_subdivision_and_terminal_nodes()
      type(htg_tree_t), target :: tree
      type(htg_status_t) :: st
      integer(id_k), pointer :: p_elder(:), p_depth(:)
      integer(id_k) :: cur_vid, elder, sum_depths
      integer(i32_k) :: depth

      ! Octree (2, 3) refined down a diagonal chain to depth 10 (11 levels)
      call tree%initialize(2_i32_k, 3_i32_k, 42_id_k)
      cur_vid = 0_id_k
      do depth = 0_i32_k, 9_i32_k
         call tree%subdivide_leaf(cur_vid, depth)
         elder = tree%get_elder_child_index(cur_vid)
         ! Pick varying child offsets along chain
         cur_vid = elder + int(mod(depth * 3_i32_k + 2_i32_k, 8_i32_k), id_k)
      end do

      call check(tree%get_number_of_levels() == 11_id_k, 'deep octree levels == 11')
      call check(tree%get_number_of_nodes() == 10_id_k, 'deep octree nodes == 10')
      call check(tree%get_number_of_leaves() == 1_id_k + 10_id_k * 7_id_k, 'deep octree leaves == 71')
      call check(tree%get_number_of_vertices() == 1_id_k + 10_id_k * 8_id_k, &
         'deep octree vertices == 81')

      ! Also refine a sibling branch at depth 1 -> depth 2 to test non-increasing level split
      call tree%subdivide_leaf(1_id_k, 1_i32_k, st)
      call check(st%ok(), 'refine sibling at depth 1 ok')
      call check(tree%get_number_of_levels() == 11_id_k, 'levels unchanged at 11')
      call check(tree%get_number_of_vertices_at_depth(2_i32_k) == 16_id_k, 'depth 2 now 16 vertices')
      call check(tree%is_terminal_node(1_id_k), 'newly split depth-1 sibling is terminal')
      call check(.not. tree%is_terminal_node(0_id_k), 'root is not terminal')

      call tree%get_vertices_per_depth_view(p_depth)
      call tree%get_parent_to_elder_child_view(p_elder)
      call check(lbound(p_depth, 1) == 0 .and. ubound(p_depth, 1) == 10, 'p_depth bounds 0:10')
      sum_depths = sum(p_depth(0:10))
      call check(sum_depths == tree%get_number_of_vertices(), 'sum of per-depth counts == n_verts')

      ! Negative subdivision & query checks
      call tree%subdivide_leaf(0_id_k, 0_i32_k, st)
      call check(st%code == HTG_ERR_TOPOLOGY, 'reject subdividing coarse root')

      call tree%subdivide_leaf(2_id_k, 0_i32_k, st)
      call check(st%code == HTG_ERR_INVALID_ARG, 'reject non-root at level 0')

      call tree%subdivide_leaf(2_id_k, 11_i32_k, st)
      call check(st%code == HTG_ERR_INVALID_ARG, 'reject level >= number_of_levels')

      call tree%subdivide_leaf(-1_id_k, 1_i32_k, st)
      call check(st%code == HTG_ERR_OUT_OF_BOUNDS, 'reject negative local_id')

      call tree%subdivide_leaf(999_id_k, 1_i32_k, st)
      call check(st%code == HTG_ERR_OUT_OF_BOUNDS, 'reject local_id >= number_of_vertices')

      elder = tree%get_elder_child_index(2_id_k, st)
      call check(st%code == HTG_ERR_TOPOLOGY .and. elder == HTG_TREE_LEAF_SENTINEL, &
         'get_elder_child_index on leaf returns sentinel and topology error')
   end subroutine test_deep_asymmetric_subdivision_and_terminal_nodes

   subroutine test_global_indexing_state_machine()
      type(htg_tree_t), target :: tree
      type(htg_status_t) :: st
      integer(id_k), pointer :: p_gid(:)
      integer(id_k) :: gid, i

      ! 1. Unset mode checks
      call tree%initialize(2_i32_k, 2_i32_k, 0_id_k)
      call tree%subdivide_leaf(0_id_k, 0_i32_k)
      call check(tree%get_indexing_mode() == HTG_INDEX_UNSET, 'initially HTG_INDEX_UNSET')
      call check(tree%get_global_index_start() == HTG_INVALID_ID, 'unset start == -1')
      call check(tree%get_max_global_index() == HTG_INVALID_ID, 'unset max_global_index == -1')
      gid = tree%get_max_global_index(st)
      call check(gid == HTG_INVALID_ID .and. st%code == HTG_ERR_INDEXING, &
         'unset max_global_index status == HTG_ERR_INDEXING')
      gid = tree%get_global_index_from_local(0_id_k, st)
      call check(gid == HTG_INVALID_ID .and. st%code == HTG_ERR_INDEXING, &
         'unset get_global_index_from_local fails with HTG_ERR_INDEXING')
      call tree%get_global_index_table_view(p_gid)
      call check(.not. associated(p_gid), 'unset global table view is null')

      ! 2. Implicit mode & mode-mixing rejection
      call tree%set_global_index_start(-5_id_k, st)
      call check(st%code == HTG_ERR_INDEXING .and. tree%get_indexing_mode() == HTG_INDEX_UNSET, &
         'negative global_index_start rejected')

      call tree%set_global_index_start(100_id_k, st)
      call check(st%ok() .and. tree%get_indexing_mode() == HTG_INDEX_IMPLICIT, &
         'switched to HTG_INDEX_IMPLICIT')
      call check(tree%get_global_index_start() == 100_id_k, 'start == 100')
      do i = 0_id_k, 4_id_k
         call check(tree%get_global_index_from_local(i) == 100_id_k + i, 'implicit offset formula')
      end do
      call check(tree%get_max_global_index() == 104_id_k, 'implicit max == 104')

      ! Re-setting implicit start (grid shift) allowed
      call tree%set_global_index_start(250_id_k, st)
      call check(st%ok(), 're-setting implicit start status ok')
      call check(tree%get_global_index_from_local(3_id_k) == 253_id_k, &
         're-setting implicit start works')

      ! Mixing explicit into implicit rejected
      call tree%set_global_index_from_local(0_id_k, 10_id_k, st)
      call check(st%code == HTG_ERR_INDEXING, 'reject explicit call in implicit mode')
      call tree%get_global_index_table_view(p_gid)
      call check(.not. associated(p_gid), 'implicit global table view is null')

      ! 3. Explicit mode & dynamic growth during subdivision
      call tree%initialize(2_i32_k, 2_i32_k, 1_id_k)
      call tree%subdivide_leaf(0_id_k, 0_i32_k) ! vertices 0..4

      call tree%set_global_index_from_local(0_id_k, -1_id_k, st)
      call check(st%code == HTG_ERR_INDEXING .and. tree%get_indexing_mode() == HTG_INDEX_UNSET, &
         'reject negative global_id')

      call tree%set_global_index_from_local(5_id_k, 10_id_k, st)
      call check(st%code == HTG_ERR_OUT_OF_BOUNDS, 'reject out-of-bounds local_id')

      call tree%set_global_index_from_local(0_id_k, 50_id_k, st)
      call check(st%ok() .and. tree%get_indexing_mode() == HTG_INDEX_EXPLICIT, &
         'switched to HTG_INDEX_EXPLICIT')

      ! Mixing implicit into explicit rejected
      call tree%set_global_index_start(0_id_k, st)
      call check(st%code == HTG_ERR_INDEXING, 'reject implicit start in explicit mode')

      ! Unassigned slot in explicit table fails cleanly
      gid = tree%get_global_index_from_local(2_id_k, st)
      call check(gid == HTG_INVALID_ID .and. st%code == HTG_ERR_INDEXING, &
         'unassigned explicit slot fails with HTG_ERR_INDEXING')

      do i = 1_id_k, 4_id_k
         call tree%set_global_index_from_local(i, 50_id_k + 10_id_k * i)
      end do
      call check(tree%get_max_global_index() == 90_id_k, 'explicit max_global_index == 90')

      ! Subdivide multiple times in explicit mode to trigger geometric capacity growth
      call tree%subdivide_leaf(1_id_k, 1_i32_k) ! adds 5..8
      call tree%subdivide_leaf(5_id_k, 2_i32_k) ! adds 9..12
      call tree%subdivide_leaf(9_id_k, 3_i32_k) ! adds 13..16
      call tree%subdivide_leaf(13_id_k, 4_i32_k) ! adds 17..20 (exceeds INITIAL_GLOBAL_CAPACITY=16)

      call tree%get_global_index_table_view(p_gid)
      call check(associated(p_gid), 'explicit global table view associated')
      call check(lbound(p_gid, 1) == 0 .and. ubound(p_gid, 1) == 20, &
         'explicit global table view bounds 0:20')
      call check(p_gid(0) == 50_id_k .and. p_gid(4) == 90_id_k, &
         'existing explicit entries preserved across growth')
      call check(all(p_gid(5:20) == HTG_INVALID_ID), &
         'newly subdivided explicit slots initialized to HTG_INVALID_ID')
   end subroutine test_global_indexing_state_machine

   subroutine test_bf_descriptor_single_leaf_and_balanced()
      integer(i32_k), parameter :: bfs(6) = [2, 2, 2, 3, 3, 3]
      integer(i32_k), parameter :: dims(6) = [1, 2, 3, 1, 2, 3]

      type(htg_tree_t) :: src_tree, dst_tree
      type(htg_bit_array_t) :: desc
      integer(id_k), allocatable :: nv_depth(:), bf_map(:)
      integer(i32_k) :: icase, bf, d, k
      integer(id_k) :: c

      ! 1. Single-leaf root round-trip
      call src_tree%initialize(2_i32_k, 3_i32_k, 0_id_k)
      call src_tree%set_global_index_start(77_id_k)
      call src_tree%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map)

      call check(desc%get_number_of_values() == 0_id_k, '1-leaf tree descriptor has 0 bits')
      call check(lbound(nv_depth, 1) == 0 .and. size(nv_depth) == 1, '1-leaf nv_depth bounds 0:0')
      call check(nv_depth(0) == 1_id_k, '1-leaf nv_depth(0) == 1')
      call check(lbound(bf_map, 1) == 0 .and. size(bf_map) == 1, '1-leaf bf_map bounds 0:0')
      call check(bf_map(0) == 77_id_k, '1-leaf bf_map(0) == 77')

      call dst_tree%initialize(2_i32_k, 3_i32_k, 0_id_k)
      call dst_tree%subdivide_leaf(0_id_k, 0_i32_k) ! dirty it first
      call dst_tree%build_from_breadth_first_order_descriptor(desc, 0_id_k)
      call check(dst_tree%get_number_of_levels() == 1_id_k, 'rebuilt 1-leaf levels == 1')
      call check(dst_tree%get_number_of_vertices() == 1_id_k, 'rebuilt 1-leaf vertices == 1')
      call check(dst_tree%get_number_of_nodes() == 0_id_k, 'rebuilt 1-leaf nodes == 0')
      call check(dst_tree%get_number_of_leaves() == 1_id_k, 'rebuilt 1-leaf leaves == 1')
      call check(dst_tree%is_leaf(0_id_k), 'rebuilt 1-leaf root is leaf')

      ! 2. Full 3-level balanced tree round-trip for all 6 (bf, dim) pairs
      do icase = 1, 6
         bf = bfs(icase)
         d = dims(icase)
         k = bf ** d

         call desc%clear()
         if (allocated(nv_depth)) deallocate(nv_depth)
         if (allocated(bf_map)) deallocate(bf_map)

         call src_tree%initialize(bf, d, 0_id_k)
         call src_tree%subdivide_leaf(0_id_k, 0_i32_k)
         do c = 1_id_k, int(k, id_k)
            call src_tree%subdivide_leaf(c, 1_i32_k)
         end do
         call src_tree%set_global_index_start(1000_id_k)

         call src_tree%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map)
         call check(desc%get_number_of_values() == 1_id_k + int(k, id_k), &
            'balanced 3-level descriptor bit count == 1 + k')
         call check(size(nv_depth) == 3 .and. nv_depth(0) == 1_id_k .and. &
            nv_depth(1) == int(k, id_k) .and. nv_depth(2) == int(k * k, id_k), &
            'balanced 3-level per-depth counts')
         call check(size(bf_map) == 1 + k + k * k, 'balanced 3-level bf_map size')

         call dst_tree%initialize(bf, d, 0_id_k)
         call dst_tree%build_from_breadth_first_order_descriptor(desc, desc%get_number_of_values())
         call check(dst_tree%get_number_of_levels() == 3_id_k, 'rebuilt balanced levels == 3')
         call check(dst_tree%get_number_of_nodes() == 1_id_k + int(k, id_k), &
            'rebuilt balanced nodes == 1 + k')
         call check(dst_tree%get_number_of_vertices() == src_tree%get_number_of_vertices(), &
            'rebuilt balanced vertices match')
         call check(dst_tree%get_number_of_leaves() == int(k * k, id_k), &
            'rebuilt balanced leaves == k^2')
         do c = 1_id_k, int(k, id_k)
            call check(dst_tree%is_terminal_node(c), 'depth-1 nodes are terminal in 3-level tree')
         end do
      end do
   end subroutine test_bf_descriptor_single_leaf_and_balanced

   subroutine test_bf_descriptor_asymmetric_and_forest_concat()
      type(htg_tree_t) :: t0, t1, r0, r1
      type(htg_bit_array_t) :: desc, desc_recheck
      integer(id_k), allocatable :: nv_depth(:), bf_map(:), nv_recheck(:), map_recheck(:)
      integer(id_k) :: bits_t0, bits_t1, i

      ! Build t0 with out-of-BFS-order subdivisions:
      ! Quadtree (f=2, d=2, k=4):
      ! Split 0 -> children 1,2,3,4 (depth 1)
      ! Split 3 -> children 5,6,7,8 (depth 2)
      ! Split 6 -> children 9,10,11,12 (depth 3)
      ! Split 1 -> children 13,14,15,16 (depth 2) -- subdivided AFTER depth-2 child 6!
      call t0%initialize(2_i32_k, 2_i32_k, 0_id_k)
      call t0%subdivide_leaf(0_id_k, 0_i32_k)
      call t0%subdivide_leaf(3_id_k, 1_i32_k)
      call t0%subdivide_leaf(6_id_k, 2_i32_k)
      call t0%subdivide_leaf(1_id_k, 1_i32_k)
      call t0%set_global_index_start(100_id_k)

      ! Build t1 (quadtree: root 0 split into 1..4, child 4 split into 5..8)
      call t1%initialize(2_i32_k, 2_i32_k, 1_id_k)
      call t1%subdivide_leaf(0_id_k, 0_i32_k)
      call t1%subdivide_leaf(4_id_k, 1_i32_k)
      call t1%set_global_index_start(500_id_k)

      ! Concatenate t0 then t1 into shared descriptor, nv_depth, bf_map
      call t0%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map)
      bits_t0 = desc%get_number_of_values()
      ! t0 has depths: d0=1, d1=4, d2=8, d3=4 -> descriptor encodes d0..d2 = 1 + 4 + 8 = 13 bits
      call check(bits_t0 == 13_id_k, 't0 descriptor bits == 13')
      call check(all(nv_depth(0:3) == [1_id_k, 4_id_k, 8_id_k, 4_id_k]), 't0 nv_depth exact')
      ! Check BFS ordering of depth 2 in bf_map:
      ! Depth 0 (slot 0): local 0 -> 100
      ! Depth 1 (slots 1..4): locals 1,2,3,4 -> 101,102,103,104
      ! Depth 2 (slots 5..12): children of 1 (locals 13..16 -> 113..116) BEFORE children of 3 (locals 5..8 -> 105..108)!
      ! Depth 3 (slots 13..16): children of 6 (locals 9..12 -> 109..112)
      call check(all(bf_map(5:8) == [113_id_k, 114_id_k, 115_id_k, 116_id_k]), &
         'BFS map orders children of local 1 before children of local 3')
      call check(all(bf_map(9:12) == [105_id_k, 106_id_k, 107_id_k, 108_id_k]), &
         'BFS map places children of local 3 second at depth 2')
      call check(all(bf_map(13:16) == [109_id_k, 110_id_k, 111_id_k, 112_id_k]), &
         'BFS map places depth 3 children of local 6 at end')

      call t1%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map)
      bits_t1 = desc%get_number_of_values() - bits_t0
      call check(bits_t1 == 5_id_k, 't1 descriptor bits == 5 (d0=1 + d1=4)')
      call check(lbound(nv_depth, 1) == 0 .and. size(nv_depth) == 7, &
         'concatenated nv_depth has 7 entries with 0: lower bound')
      call check(lbound(bf_map, 1) == 0 .and. size(bf_map) == 17 + 9, &
         'concatenated bf_map has 26 entries with 0: lower bound')

      ! Rebuild r0 from [0, bits_t0) and r1 from [bits_t0, bits_t0 + bits_t1)
      call r0%initialize(2_i32_k, 2_i32_k, 0_id_k)
      call r0%build_from_breadth_first_order_descriptor(desc, bits_t0, 0_id_k)
      call r0%set_global_index_start(100_id_k)

      call r1%initialize(2_i32_k, 2_i32_k, 1_id_k)
      call r1%build_from_breadth_first_order_descriptor(desc, bits_t1, bits_t0)
      call r1%set_global_index_start(500_id_k)

      call check(r0%get_number_of_levels() == 4_id_k .and. &
         r0%get_number_of_vertices() == 17_id_k .and. r0%get_number_of_nodes() == 4_id_k, &
         'r0 topology counts match t0')
      call check(r1%get_number_of_levels() == 3_id_k .and. &
         r1%get_number_of_vertices() == 9_id_k .and. r1%get_number_of_nodes() == 2_id_k, &
         'r1 topology counts match t1')

      ! Re-computing descriptors from r0 + r1 must match desc & nv_depth bit-for-bit!
      call r0%compute_breadth_first_order_descriptor(desc_recheck, nv_recheck, map_recheck)
      call r1%compute_breadth_first_order_descriptor(desc_recheck, nv_recheck, map_recheck)
      call check(desc_recheck%get_number_of_values() == desc%get_number_of_values(), &
         'round-trip concatenated descriptor length identical')
      do i = 0_id_k, desc%get_number_of_values() - 1_id_k
         call check(desc_recheck%get_value(i) .eqv. desc%get_value(i), &
            'round-trip concatenated descriptor bit identical')
      end do
      call check(all(nv_recheck == nv_depth), 'round-trip concatenated nv_depth identical')
   end subroutine test_bf_descriptor_asymmetric_and_forest_concat

   subroutine test_bf_descriptor_mask_and_depth_limiter()
      type(htg_tree_t) :: tree
      type(htg_bit_array_t) :: mask, desc
      integer(id_k), allocatable :: nv_depth(:), bf_map(:)

      ! Bintree (f=2, d=1, k=2):
      ! 0 -> [1, 2] (depth 1)
      ! 2 -> [3, 4] (depth 2)
      ! 3 -> [5, 6] (depth 3)
      call tree%initialize(2_i32_k, 1_i32_k, 0_id_k)
      call tree%subdivide_leaf(0_id_k, 0_i32_k)
      call tree%subdivide_leaf(2_id_k, 1_i32_k)
      call tree%subdivide_leaf(3_id_k, 2_i32_k)
      call tree%set_global_index_start(10_id_k) ! global IDs 10..16

      ! Case 1: Mask coarse node local 2 (global 12) -> prunes [3,4,5,6], reduces effective levels from 4 to 2!
      call mask%initialize(20_id_k, .false.)
      call mask%set_value(12_id_k, .true.)

      call tree%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map, input_mask=mask)
      call check(size(nv_depth) == 2, 'masking deep branch squeezes trailing empty depths to 2')
      call check(all(nv_depth(0:1) == [1_id_k, 2_id_k]), 'masked nv_depth == [1, 2]')
      call check(desc%get_number_of_values() == 1_id_k, 'masked descriptor has 1 bit')
      call check(desc%get_value(0_id_k), 'masked descriptor bit 0 is true')
      call check(all(bf_map(0:2) == [10_id_k, 11_id_k, 12_id_k]), &
         'masked coarse node 12 still listed as leaf at depth 1, descendants omitted')

      ! Case 2: Mask root itself (global 10) -> reduces effective levels to 1, 0 descriptor bits!
      call desc%clear()
      deallocate(nv_depth, bf_map)
      call mask%set_value(10_id_k, .true.)
      call tree%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map, input_mask=mask)
      call check(size(nv_depth) == 1 .and. nv_depth(0) == 1_id_k, 'masked root nv_depth == [1]')
      call check(desc%get_number_of_values() == 0_id_k, 'masked root descriptor has 0 bits')
      call check(size(bf_map) == 1 .and. bf_map(0) == 10_id_k, 'masked root bf_map == [10]')

      ! Case 3: depth_limiter = 2 (unmasked) -> visits depths 0, 1, 2; prunes depth 3 ([5, 6])
      call desc%clear()
      deallocate(nv_depth, bf_map)
      call tree%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map, depth_limiter=2_i32_k)
      call check(size(nv_depth) == 3, 'depth_limiter=2 gives 3 levels (depths 0, 1, 2)')
      call check(all(nv_depth(0:2) == [1_id_k, 2_id_k, 2_id_k]), 'depth_limiter=2 nv_depth')
      call check(desc%get_number_of_values() == 3_id_k, 'depth_limiter=2 descriptor bits == 3')
      call check(desc%get_value(0_id_k), 'depth_limiter=2 descriptor bit 0 == 1')
      call check(.not. desc%get_value(1_id_k), 'depth_limiter=2 descriptor bit 1 == 0')
      call check(desc%get_value(2_id_k), 'depth_limiter=2 descriptor bit 2 == 1')
      call check(all(bf_map(0:4) == [10_id_k, 11_id_k, 12_id_k, 13_id_k, 14_id_k]), &
         'depth_limiter=2 bf_map')

      ! Case 4: depth_limiter = 0 -> root only
      call desc%clear()
      deallocate(nv_depth, bf_map)
      call tree%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map, depth_limiter=0_i32_k)
      call check(size(nv_depth) == 1 .and. desc%get_number_of_values() == 0_id_k, &
         'depth_limiter=0 collapses to single root leaf')
   end subroutine test_bf_descriptor_mask_and_depth_limiter

   subroutine test_bf_descriptor_malformed_rejection()
      type(htg_tree_t) :: tree
      type(htg_bit_array_t) :: desc
      type(htg_status_t) :: st
      integer(id_k), allocatable :: nv_depth(:), bf_map(:)

      call tree%initialize(2_i32_k, 3_i32_k, 0_id_k) ! octree k=8
      call tree%subdivide_leaf(0_id_k, 0_i32_k) ! 9 vertices, 2 levels

      ! 1. Negative bit count or start_index
      call desc%initialize(16_id_k, .false.)
      call tree%build_from_breadth_first_order_descriptor(desc, -1_id_k, 0_id_k, st)
      call check(st%code == HTG_ERR_INVALID_ARG, 'reject negative number_of_bits')
      call tree%build_from_breadth_first_order_descriptor(desc, 1_id_k, -1_id_k, st)
      call check(st%code == HTG_ERR_INVALID_ARG, 'reject negative start_index')

      ! 2. Out-of-bounds slice
      call tree%build_from_breadth_first_order_descriptor(desc, 10_id_k, 8_id_k, st)
      call check(st%code == HTG_ERR_OUT_OF_BOUNDS, 'reject out-of-bounds slice')

      ! 3. Truncated level: root=1 (expects 8 bits at depth 1), but pass number_of_bits = 6 (only 5 depth-1 bits)
      call desc%fill(.false.)
      call desc%set_value(0_id_k, .true.)
      call desc%set_value(2_id_k, .true.)
      call tree%build_from_breadth_first_order_descriptor(desc, 6_id_k, 0_id_k, st)
      call check(st%code == HTG_ERR_TOPOLOGY, 'reject truncated breadth-first level')
      call check(tree%get_number_of_vertices() == 9_id_k, &
         'failed build leaves existing tree state intact')

      ! 4. Terminal all-zero leaf depth encoded: [0] with number_of_bits = 1
      call desc%fill(.false.)
      call tree%build_from_breadth_first_order_descriptor(desc, 1_id_k, 0_id_k, st)
      call check(st%code == HTG_ERR_TOPOLOGY, 'reject terminal leaf depth [0] with 1 bit')

      ! 5. Terminal all-zero depth 1 encoded: [1, 0,0,0,0,0,0,0,0] with number_of_bits = 9
      call desc%fill(.false.)
      call desc%set_value(0_id_k, .true.)
      call tree%build_from_breadth_first_order_descriptor(desc, 9_id_k, 0_id_k, st)
      call check(st%code == HTG_ERR_TOPOLOGY, 'reject terminal all-zero depth-1 bits')

      ! 6. Trailing orphan bits after all-zero depth: [1, 0,0,0,0,0,0,0,0, 1] with 10 bits
      call desc%set_value(9_id_k, .true.)
      call tree%build_from_breadth_first_order_descriptor(desc, 10_id_k, 0_id_k, st)
      call check(st%code == HTG_ERR_TOPOLOGY, 'reject trailing orphan bits')

      ! 7. compute_breadth_first_order_descriptor error guards
      call tree%compute_breadth_first_order_descriptor(desc, nv_depth, bf_map, status=st)
      call check(st%code == HTG_ERR_INDEXING, 'reject compute_bf when global indexing unset')

      call tree%set_global_index_start(0_id_k)
      call tree%compute_breadth_first_order_descriptor( &
         desc, nv_depth, bf_map, depth_limiter=-1_i32_k, status=st)
      call check(st%code == HTG_ERR_INVALID_ARG, 'reject negative depth_limiter')

      call desc%initialize(4_id_k, .false.) ! too short for global IDs 0..8
      call tree%compute_breadth_first_order_descriptor( &
         desc, nv_depth, bf_map, input_mask=desc, status=st)
      call check(st%code == HTG_ERR_MASK, 'reject undersized input_mask')
   end subroutine test_bf_descriptor_malformed_rejection

   subroutine test_vtk_9_7_2_audit_exact_vectors()
      type(htg_tree_t), target :: qtree, btree, rebuilt
      type(htg_bit_array_t) :: desc
      integer(id_k), pointer :: p_elder(:), p_depth(:)
      integer(id_k), allocatable :: vpd(:), id_map(:)
      logical, parameter :: expected_bits(0:6) = &
         [.true., .true., .true., .false., .false., .true., .false.]
      integer(id_k) :: i

      ! Audit 3.2: Quadtree (f=2, d=2, k=4) lazy ParentToElderChild sizing
      call qtree%initialize(2_i32_k, 2_i32_k, 0_id_k)
      call qtree%get_parent_to_elder_child_view(p_elder)
      call check(size(p_elder) == 1 .and. p_elder(0) == 4294967295_id_k, &
         'audit 3.2 step 1: [4294967295]')

      call qtree%subdivide_leaf(0_id_k, 0_i32_k)
      call qtree%get_parent_to_elder_child_view(p_elder)
      call check(size(p_elder) == 1 .and. p_elder(0) == 1_id_k, &
         'audit 3.2 step 2: [1]')

      call qtree%subdivide_leaf(3_id_k, 1_i32_k)
      call qtree%get_parent_to_elder_child_view(p_elder)
      call check(size(p_elder) == 4, 'audit 3.2 step 3: size 4')
      call check(all(p_elder(0:3) == [1_id_k, 4294967295_id_k, 4294967295_id_k, 5_id_k]), &
         'audit 3.2 step 3: [1, 4294967295, 4294967295, 5]')
      call check(.not. qtree%is_terminal_node(1_id_k), &
         'audit 3.2 IsTerminalNode quirk fix: intermediate leaf 1 is NOT terminal node')
      call check(qtree%is_terminal_node(3_id_k), 'audit 3.2 coarse node 3 IS terminal node')

      ! Audit 5.2: Out-of-order bintree (GlobalIndexStart=100, ParentToElderChild=[1, 5, 3, 7])
      call btree%initialize(2_i32_k, 1_i32_k, 0_id_k)
      call btree%subdivide_leaf(0_id_k, 0_i32_k) ! 0 -> [1, 2]
      call btree%subdivide_leaf(2_id_k, 1_i32_k) ! 2 -> [3, 4]
      call btree%subdivide_leaf(1_id_k, 1_i32_k) ! 1 -> [5, 6]
      call btree%subdivide_leaf(3_id_k, 2_i32_k) ! 3 -> [7, 8]
      call btree%set_global_index_start(100_id_k)

      call btree%get_parent_to_elder_child_view(p_elder)
      call check(size(p_elder) == 4 .and. all(p_elder(0:3) == [1_id_k, 5_id_k, 3_id_k, 7_id_k]), &
         'audit 5.2 ParentToElderChild == [1, 5, 3, 7]')

      call btree%compute_breadth_first_order_descriptor(desc, vpd, id_map)
      call check(size(vpd) == 4 .and. all(vpd(0:3) == [1_id_k, 2_id_k, 4_id_k, 2_id_k]), &
         'audit 5.2 vpd == [1, 2, 4, 2]')
      call check(desc%get_number_of_values() == 7_id_k, 'audit 5.2 descriptor 7 bits')
      do i = 0_id_k, 6_id_k
         call check(desc%get_value(i) .eqv. expected_bits(i), 'audit 5.2 descriptor 1110010')
      end do
      call check(size(id_map) == 9 .and. all(id_map(0:8) == &
         [100_id_k, 101_id_k, 102_id_k, 105_id_k, 106_id_k, 103_id_k, 104_id_k, 107_id_k, 108_id_k]), &
         'audit 5.2 breadthFirstIdMap exact match')

      ! Rebuild from audit 5.2 descriptor and check ParentToElderChild size == max(1, numberOfBits) == 7
      call rebuilt%initialize(2_i32_k, 1_i32_k, 0_id_k)
      call rebuilt%build_from_breadth_first_order_descriptor(desc, 7_id_k, 0_id_k)
      call rebuilt%get_parent_to_elder_child_view(p_elder)
      call rebuilt%get_vertices_per_depth_view(p_depth)
      call check(size(p_elder) == 7, 'audit 5.2 rebuilt ParentToElderChild size == 7')
      call check(all(p_elder(0:6) == &
         [1_id_k, 3_id_k, 5_id_k, 4294967295_id_k, 4294967295_id_k, 7_id_k, 4294967295_id_k]), &
         'audit 5.2 rebuilt ParentToElderChild exact match')
      call check(all(p_depth(0:3) == [1_id_k, 2_id_k, 4_id_k, 2_id_k]), &
         'audit 5.2 rebuilt vertices_per_depth exact match')

      ! Audit 5.2 depthLimiter = 1
      call desc%clear()
      deallocate(vpd, id_map)
      call btree%compute_breadth_first_order_descriptor(desc, vpd, id_map, depth_limiter=1_i32_k)
      call check(size(vpd) == 2 .and. all(vpd(0:1) == [1_id_k, 2_id_k]), &
         'audit 5.2 depthLimiter=1 vpd == [1, 2]')
      call check(size(id_map) == 3 .and. all(id_map(0:2) == [100_id_k, 101_id_k, 102_id_k]), &
         'audit 5.2 depthLimiter=1 idMap == [100, 101, 102]')
      call check(desc%get_number_of_values() == 1_id_k, 'audit 5.2 depthLimiter=1 descriptor len 1')
      call check(desc%get_value(0_id_k), 'audit 5.2 depthLimiter=1 descriptor == [1]')
   end subroutine test_vtk_9_7_2_audit_exact_vectors

end program test_tree
