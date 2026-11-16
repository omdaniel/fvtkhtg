module fvtkhtg_validate_m
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use fvtkhtg_kinds_m, only: i32_k, id_k, r64_k, htg_status_t, &
      HTG_ERR_INVALID_ARG, HTG_ERR_OUT_OF_BOUNDS, HTG_ERR_INVALID_STATE, &
      HTG_ERR_TOPOLOGY, HTG_ERR_INDEXING, HTG_ERR_MASK, HTG_ERR_INTERFACE, &
      HTG_INVALID_ID, HTG_LEAF_ELDER_SENTINEL, &
      HTG_INDEX_UNSET, HTG_INDEX_IMPLICIT, HTG_INDEX_EXPLICIT, &
      HTG_TYPE_BIT, htg_set_ok, htg_fail
   use fvtkhtg_bit_array_m, only: htg_bit_array_t
   use fvtkhtg_tree_m, only: htg_tree_t, HTG_TREE_LEAF_SENTINEL
   use fvtkhtg_grid_m, only: htg_grid_t, htg_compute_squeezed_geometry
   implicit none
   private

   public :: htg_validate_tree
   public :: htg_validate_grid

contains

   ! --- In-place O(N log N) heapsort helpers (zero recursion, O(1) space) ---

   pure subroutine sift_down_i64(arr, root_in, end_idx)
      integer(id_k), intent(inout) :: arr(0:)
      integer(id_k), intent(in) :: root_in, end_idx

      integer(id_k) :: root, child, swap_idx, tmp

      root = root_in
      do while (2_id_k * root + 1_id_k <= end_idx)
         child = 2_id_k * root + 1_id_k
         swap_idx = root
         if (arr(swap_idx) < arr(child)) swap_idx = child
         if (child + 1_id_k <= end_idx) then
            if (arr(swap_idx) < arr(child + 1_id_k)) swap_idx = child + 1_id_k
         end if
         if (swap_idx == root) exit
         tmp = arr(root)
         arr(root) = arr(swap_idx)
         arr(swap_idx) = tmp
         root = swap_idx
      end do
   end subroutine sift_down_i64

   pure subroutine heapsort_i64(arr, n)
      integer(id_k), intent(inout) :: arr(0:)
      integer(id_k), intent(in) :: n

      integer(id_k) :: start_idx, end_idx, tmp

      if (n <= 1_id_k) return
      do start_idx = (n - 2_id_k) / 2_id_k, 0_id_k, -1_id_k
         call sift_down_i64(arr, start_idx, n - 1_id_k)
      end do
      do end_idx = n - 1_id_k, 1_id_k, -1_id_k
         tmp = arr(end_idx)
         arr(end_idx) = arr(0_id_k)
         arr(0_id_k) = tmp
         call sift_down_i64(arr, 0_id_k, end_idx - 1_id_k)
      end do
   end subroutine heapsort_i64

   pure subroutine sift_down_intervals(lo_arr, hi_arr, root_in, end_idx)
      integer(id_k), intent(inout) :: lo_arr(0:), hi_arr(0:)
      integer(id_k), intent(in) :: root_in, end_idx

      integer(id_k) :: root, child, swap_idx, tmp_lo, tmp_hi

      root = root_in
      do while (2_id_k * root + 1_id_k <= end_idx)
         child = 2_id_k * root + 1_id_k
         swap_idx = root
         if (lo_arr(swap_idx) < lo_arr(child) .or. &
             (lo_arr(swap_idx) == lo_arr(child) .and. hi_arr(swap_idx) < hi_arr(child))) then
            swap_idx = child
         end if
         if (child + 1_id_k <= end_idx) then
            if (lo_arr(swap_idx) < lo_arr(child + 1_id_k) .or. &
                (lo_arr(swap_idx) == lo_arr(child + 1_id_k) .and. &
                 hi_arr(swap_idx) < hi_arr(child + 1_id_k))) then
               swap_idx = child + 1_id_k
            end if
         end if
         if (swap_idx == root) exit
         tmp_lo = lo_arr(root)
         tmp_hi = hi_arr(root)
         lo_arr(root) = lo_arr(swap_idx)
         hi_arr(root) = hi_arr(swap_idx)
         lo_arr(swap_idx) = tmp_lo
         hi_arr(swap_idx) = tmp_hi
         root = swap_idx
      end do
   end subroutine sift_down_intervals

   pure subroutine heapsort_intervals(lo_arr, hi_arr, n)
      integer(id_k), intent(inout) :: lo_arr(0:), hi_arr(0:)
      integer(id_k), intent(in) :: n

      integer(id_k) :: start_idx, end_idx, tmp_lo, tmp_hi

      if (n <= 1_id_k) return
      do start_idx = (n - 2_id_k) / 2_id_k, 0_id_k, -1_id_k
         call sift_down_intervals(lo_arr, hi_arr, start_idx, n - 1_id_k)
      end do
      do end_idx = n - 1_id_k, 1_id_k, -1_id_k
         tmp_lo = lo_arr(end_idx)
         tmp_hi = hi_arr(end_idx)
         lo_arr(end_idx) = lo_arr(0_id_k)
         hi_arr(end_idx) = hi_arr(0_id_k)
         lo_arr(0_id_k) = tmp_lo
         hi_arr(0_id_k) = tmp_hi
         call sift_down_intervals(lo_arr, hi_arr, 0_id_k, end_idx - 1_id_k)
      end do
   end subroutine heapsort_intervals

   ! --- BitArray internal integrity check ---

   pure subroutine validate_bit_array_internal(bits, err_code, label, status)
      type(htg_bit_array_t), intent(in) :: bits
      integer(i32_k), intent(in) :: err_code
      character(len=*), intent(in) :: label
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: active_bytes
      integer(i32_k) :: rem, pos

      if (bits%num_bits < 0_id_k .or. bits%byte_capacity < 0_id_k) then
         call htg_fail(err_code, label // ': negative num_bits or byte_capacity', status)
         return
      end if

      if (allocated(bits%bytes)) then
         if (lbound(bits%bytes, 1) /= 0) then
            call htg_fail(err_code, label // ': bytes lower bound is not 0', status)
            return
         end if
         if (size(bits%bytes, kind=id_k) /= bits%byte_capacity) then
            call htg_fail(err_code, label // ': bytes size does not match byte_capacity', status)
            return
         end if
      end if

      if (bits%num_bits > 0_id_k) then
         active_bytes = (bits%num_bits + 7_id_k) / 8_id_k
         if (.not. allocated(bits%bytes) .or. active_bytes > bits%byte_capacity) then
            call htg_fail(err_code, label // ': insufficient allocated bytes for num_bits', status)
            return
         end if
         rem = int(mod(bits%num_bits, 8_id_k), i32_k)
         if (rem /= 0_i32_k) then
            do pos = 0_i32_k, 7_i32_k - rem
               if (btest(bits%bytes(active_bytes - 1_id_k), pos)) then
                  call htg_fail(err_code, label // ': non-zero trailing slack bit in final byte', status)
                  return
               end if
            end do
         end if
      end if

      call htg_set_ok(status)
   end subroutine validate_bit_array_internal

   ! --- Coordinate monotonicity & length audit ---

   pure subroutine validate_axis_coordinates(coords, has_flag, expected_len, axis_name, status)
      real(r64_k), allocatable, intent(in) :: coords(:)
      logical, intent(in) :: has_flag
      integer(i32_k), intent(in) :: expected_len
      character(len=*), intent(in) :: axis_name
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: n, i
      logical :: increasing

      if (.not. has_flag) then
         if (allocated(coords)) then
            call htg_fail(HTG_ERR_INVALID_STATE, &
               'htg_validate_grid: ' // axis_name // ' coordinates allocated while has_flag is false', status)
            return
         end if
         if (expected_len > 1_i32_k) then
            call htg_fail(HTG_ERR_INVALID_STATE, &
               'htg_validate_grid: missing ' // axis_name // ' coordinates for active axis', status)
            return
         end if
         call htg_set_ok(status)
         return
      end if

      if (.not. allocated(coords)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_validate_grid: ' // axis_name // ' coordinates not allocated despite has_flag', status)
         return
      end if

      if (lbound(coords, 1) /= 0) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_validate_grid: ' // axis_name // ' coordinates lower bound must be 0', status)
         return
      end if

      n = size(coords, kind=id_k)
      if (n /= int(expected_len, id_k)) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_validate_grid: ' // axis_name // ' coordinate length mismatch with dimensions', status)
         return
      end if

      do i = 0_id_k, n - 1_id_k
         if (.not. ieee_is_finite(coords(i))) then
            call htg_fail(HTG_ERR_INVALID_ARG, &
               'htg_validate_grid: non-finite value in ' // axis_name // ' coordinates', status)
            return
         end if
      end do

      if (n >= 2_id_k) then
         if (.not. (coords(1) > coords(0)) .and. .not. (coords(1) < coords(0))) then
            call htg_fail(HTG_ERR_INVALID_ARG, &
               'htg_validate_grid: duplicate entries in ' // axis_name // ' coordinates', status)
            return
         end if
         increasing = (coords(1) > coords(0))
         do i = 1_id_k, n - 1_id_k
            if (increasing) then
               if (.not. (coords(i) > coords(i - 1_id_k))) then
                  call htg_fail(HTG_ERR_INVALID_ARG, &
                     'htg_validate_grid: ' // axis_name // ' coordinates are not strictly monotonic', status)
                  return
               end if
            else
               if (.not. (coords(i) < coords(i - 1_id_k))) then
                  call htg_fail(HTG_ERR_INVALID_ARG, &
                     'htg_validate_grid: ' // axis_name // ' coordinates are not strictly monotonic', status)
                  return
               end if
            end if
         end do
      end if

      call htg_set_ok(status)
   end subroutine validate_axis_coordinates

   ! --- Comprehensive Single-Tree Structural Invariant Audit ---

   subroutine htg_validate_tree(tree, status, require_global_index)
      type(htg_tree_t), intent(in) :: tree
      type(htg_status_t), intent(out), optional :: status
      logical, intent(in), optional :: require_global_index

      logical :: req_gid, is_leaf_v
      integer(i32_k) :: bf, ndim, k, expected_k
      integer(id_k) :: nlevels, nv, nnodes, nleaves, k64
      integer(id_k) :: lvl, cnt, prev_cnt, depth_sum
      integer(id_k) :: sp, new_cap, v, d, elder, c, child
      integer(id_k) :: obs_nodes, obs_leaves, ptec_hi, raw_elder
      logical, allocatable :: visited(:)
      integer(id_k), allocatable :: obs_per_depth(:)
      integer(id_k), allocatable :: stk_v(:), stk_d(:), tmp_v(:), tmp_d(:)
      integer(id_k), allocatable :: sorted_gids(:)
      type(htg_status_t) :: st

      req_gid = .false.
      if (present(require_global_index)) req_gid = require_global_index

      if (.not. tree%is_initialized()) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_validate_tree: tree is not initialized', status)
         return
      end if

      bf = tree%branch_factor
      ndim = tree%dimension
      k = tree%number_of_children

      if (bf /= 2_i32_k .and. bf /= 3_i32_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_validate_tree: invalid branch_factor (expected 2 or 3)', status)
         return
      end if

      if (ndim < 1_i32_k .or. ndim > 3_i32_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_validate_tree: invalid dimension (expected 1, 2, or 3)', status)
         return
      end if

      expected_k = bf ** ndim
      if (k /= expected_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: number_of_children /= branch_factor ** dimension', status)
         return
      end if
      k64 = int(k, id_k)

      nlevels = tree%number_of_levels
      nv = tree%number_of_vertices
      nnodes = tree%number_of_nodes
      nleaves = tree%number_of_leaves

      if (nlevels < 1_id_k .or. nv < 1_id_k .or. nnodes < 0_id_k .or. nleaves < 1_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: non-positive vertex/level/leaf counter or negative node counter', status)
         return
      end if

      ! k-ary tree counting identities:
      ! N_v = 1 + k * N_n,  N_l = 1 + (k - 1) * N_n,  N_v = N_n + N_l
      if (nnodes > (huge(0_id_k) - 1_id_k) / k64) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: counter overflow in k * number_of_nodes', status)
         return
      end if

      if (nv /= 1_id_k + k64 * nnodes) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: violated identity N_vertices == 1 + k * N_nodes', status)
         return
      end if

      if (nleaves /= 1_id_k + (k64 - 1_id_k) * nnodes) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: violated identity N_leaves == 1 + (k - 1) * N_nodes', status)
         return
      end if

      if (nv /= nnodes + nleaves) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: violated identity N_vertices == N_nodes + N_leaves', status)
         return
      end if

      if (nnodes == 0_id_k .and. nlevels /= 1_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: single-leaf root tree must have number_of_levels == 1', status)
         return
      end if

      if (nnodes > 0_id_k .and. nlevels < 2_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: refined tree with N_nodes > 0 must have number_of_levels >= 2', status)
         return
      end if

      ! Per-depth array bounds & sum checks
      if (.not. allocated(tree%number_of_vertices_per_depth)) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: number_of_vertices_per_depth is not allocated', status)
         return
      end if

      if (lbound(tree%number_of_vertices_per_depth, 1) /= 0) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: number_of_vertices_per_depth lower bound must be 0', status)
         return
      end if

      if (ubound(tree%number_of_vertices_per_depth, 1, kind=id_k) < nlevels - 1_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: number_of_vertices_per_depth shorter than number_of_levels', status)
         return
      end if

      if (tree%number_of_vertices_per_depth(0) /= 1_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: number_of_vertices_per_depth(0) must be 1', status)
         return
      end if

      depth_sum = 0_id_k
      do lvl = 0_id_k, nlevels - 1_id_k
         cnt = tree%number_of_vertices_per_depth(lvl)
         if (cnt <= 0_id_k) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_validate_tree: active depth level has <= 0 vertices', status)
            return
         end if
         if (lvl >= 1_id_k) then
            prev_cnt = tree%number_of_vertices_per_depth(lvl - 1_id_k)
            if ((cnt + k64 - 1_id_k) / k64 > prev_cnt) then
               call htg_fail(HTG_ERR_TOPOLOGY, &
                  'htg_validate_tree: N_depth(lvl) exceeds k * N_depth(lvl - 1)', status)
               return
            end if
         end if
         depth_sum = depth_sum + cnt
      end do

      if (depth_sum /= nv) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: sum(number_of_vertices_per_depth) /= number_of_vertices', status)
         return
      end if

      ! ParentToElderChild raw container audit
      if (tree%elder_active_len < 0_id_k .or. tree%elder_active_len > nv) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: elder_active_len out of range [0, number_of_vertices]', status)
         return
      end if

      if (allocated(tree%parent_to_elder_child)) then
         if (lbound(tree%parent_to_elder_child, 1) /= 0) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_validate_tree: parent_to_elder_child lower bound must be 0', status)
            return
         end if
         if (size(tree%parent_to_elder_child, kind=id_k) < tree%elder_active_len) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_validate_tree: parent_to_elder_child smaller than elder_active_len', status)
            return
         end if
      else if (nnodes > 0_id_k .or. tree%elder_active_len > 0_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: parent_to_elder_child unallocated despite active entries', status)
         return
      end if

      if (allocated(tree%parent_to_elder_child) .and. tree%elder_active_len > 0_id_k) then
         ptec_hi = tree%elder_active_len - 1_id_k
         obs_nodes = 0_id_k
         do v = 0_id_k, ptec_hi
            raw_elder = tree%parent_to_elder_child(v)
            if (raw_elder /= HTG_TREE_LEAF_SENTINEL .and. &
                raw_elder /= HTG_LEAF_ELDER_SENTINEL) then
               if (nv == 1_id_k .or. raw_elder < 1_id_k .or. &
                   raw_elder + k64 - 1_id_k >= nv) then
                  call htg_fail(HTG_ERR_TOPOLOGY, &
                     'htg_validate_tree: parent_to_elder_child entry out of range', &
                     status)
                  return
               end if
               obs_nodes = obs_nodes + 1_id_k
            end if
         end do
         if (obs_nodes /= nnodes) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_validate_tree: coarse entries in parent_to_elder_child /= number_of_nodes', &
               status)
            return
         end if
      else if (nnodes > 0_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: elder_active_len == 0 while number_of_nodes > 0', status)
         return
      end if

      ! Full iterative DFS topological traversal from root (vertex 0, depth 0)
      allocate(visited(0_id_k : nv - 1_id_k))
      visited = .false.
      allocate(obs_per_depth(0_id_k : nlevels - 1_id_k))
      obs_per_depth = 0_id_k

      allocate(stk_v(0_id_k : 63_id_k))
      allocate(stk_d(0_id_k : 63_id_k))

      visited(0) = .true.
      stk_v(0) = 0_id_k
      stk_d(0) = 0_id_k
      sp = 0_id_k
      obs_nodes = 0_id_k
      obs_leaves = 0_id_k

      do while (sp >= 0_id_k)
         v = stk_v(sp)
         d = stk_d(sp)
         sp = sp - 1_id_k

         if (d < 0_id_k .or. d >= nlevels) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_validate_tree: vertex depth exceeds number_of_levels - 1', status)
            return
         end if
         obs_per_depth(d) = obs_per_depth(d) + 1_id_k

         is_leaf_v = tree%is_leaf(v, st)
         if (.not. st%ok()) then
            call htg_fail(HTG_ERR_TOPOLOGY, trim(st%message), status)
            return
         end if

         if (is_leaf_v) then
            obs_leaves = obs_leaves + 1_id_k
         else
            obs_nodes = obs_nodes + 1_id_k
            if (d + 1_id_k >= nlevels) then
               call htg_fail(HTG_ERR_TOPOLOGY, &
                  'htg_validate_tree: coarse node at deepest level has children beyond number_of_levels', status)
               return
            end if

            elder = tree%get_elder_child_index(v, st)
            if (.not. st%ok()) then
               call htg_fail(HTG_ERR_TOPOLOGY, trim(st%message), status)
               return
            end if

            if (elder < 1_id_k .or. elder + k64 - 1_id_k >= nv) then
               call htg_fail(HTG_ERR_TOPOLOGY, &
                  'htg_validate_tree: elder child block outside [1, number_of_vertices - 1]', status)
               return
            end if

            if (sp + k64 > ubound(stk_v, 1, kind=id_k)) then
               new_cap = max(2_id_k * size(stk_v, kind=id_k), sp + k64 + 1_id_k)
               allocate(tmp_v(0_id_k : new_cap - 1_id_k))
               allocate(tmp_d(0_id_k : new_cap - 1_id_k))
               if (sp >= 0_id_k) then
                  tmp_v(0_id_k : sp) = stk_v(0_id_k : sp)
                  tmp_d(0_id_k : sp) = stk_d(0_id_k : sp)
               end if
               call move_alloc(tmp_v, stk_v)
               call move_alloc(tmp_d, stk_d)
            end if

            do c = k64 - 1_id_k, 0_id_k, -1_id_k
               child = elder + c
               if (visited(child)) then
                  call htg_fail(HTG_ERR_TOPOLOGY, &
                     'htg_validate_tree: cycle or shared child detected in parent_to_elder_child', status)
                  return
               end if
               visited(child) = .true.
               sp = sp + 1_id_k
               stk_v(sp) = child
               stk_d(sp) = d + 1_id_k
            end do

            ! Sibling blocks tile 1 .. nv-1 in steps of k, and the parent table
            ! (used by get_parent_index, get_level and coarsen_node) must agree
            ! with parent_to_elder_child.
            if (mod(elder - 1_id_k, k64) /= 0_id_k) then
               call htg_fail(HTG_ERR_TOPOLOGY, &
                  'htg_validate_tree: child block does not start at 1 + m * number_of_children', status)
               return
            end if
            if (.not. allocated(tree%block_parent)) then
               call htg_fail(HTG_ERR_TOPOLOGY, &
                  'htg_validate_tree: block_parent is not allocated on a refined tree', status)
               return
            end if
            if (lbound(tree%block_parent, 1) /= 0 .or. &
                ubound(tree%block_parent, 1, kind=id_k) < (elder - 1_id_k) / k64) then
               call htg_fail(HTG_ERR_TOPOLOGY, &
                  'htg_validate_tree: block_parent bounds do not cover every sibling block', status)
               return
            end if
            if (tree%block_parent((elder - 1_id_k) / k64) /= v) then
               call htg_fail(HTG_ERR_TOPOLOGY, &
                  'htg_validate_tree: block_parent disagrees with parent_to_elder_child', status)
               return
            end if
         end if
      end do

      if (.not. all(visited)) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: orphan / unreachable vertices detected', status)
         return
      end if

      if (obs_nodes /= nnodes .or. obs_leaves /= nleaves) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_tree: observed DFS node/leaf counts disagree with stored counters', status)
         return
      end if

      do lvl = 0_id_k, nlevels - 1_id_k
         if (obs_per_depth(lvl) /= tree%number_of_vertices_per_depth(lvl)) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_validate_tree: observed per-depth vertex count disagrees with number_of_vertices_per_depth', status)
            return
         end if
      end do

      ! Global indexing audit within the tree
      select case (tree%indexing_mode)
      case (HTG_INDEX_UNSET)
         if (req_gid) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_validate_tree: global indexing is unset on active tree', status)
            return
         end if
         if (tree%global_index_start /= HTG_INVALID_ID) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_validate_tree: global_index_start set while indexing_mode is UNSET', status)
            return
         end if

      case (HTG_INDEX_IMPLICIT)
         if (tree%global_index_start < 0_id_k) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_validate_tree: negative global_index_start in IMPLICIT mode', status)
            return
         end if
         if (tree%global_index_start > huge(0_id_k) - (nv - 1_id_k)) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_validate_tree: int64 overflow in global_index_start + number_of_vertices', status)
            return
         end if

      case (HTG_INDEX_EXPLICIT)
         if (tree%global_index_start /= HTG_INVALID_ID) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_validate_tree: global_index_start must be unset in EXPLICIT mode', status)
            return
         end if
         if (.not. allocated(tree%global_index_table)) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_validate_tree: global_index_table unallocated in EXPLICIT mode', status)
            return
         end if
         if (lbound(tree%global_index_table, 1) /= 0 .or. &
             ubound(tree%global_index_table, 1, kind=id_k) < nv - 1_id_k) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_validate_tree: global_index_table bounds invalid for number_of_vertices', status)
            return
         end if

         allocate(sorted_gids(0_id_k : nv - 1_id_k))
         do v = 0_id_k, nv - 1_id_k
            if (tree%global_index_table(v) < 0_id_k) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_validate_tree: unassigned/negative global index in EXPLICIT global_index_table', status)
               return
            end if
            sorted_gids(v) = tree%global_index_table(v)
         end do
         call heapsort_i64(sorted_gids, nv)
         do v = 1_id_k, nv - 1_id_k
            if (sorted_gids(v) == sorted_gids(v - 1_id_k)) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_validate_tree: duplicate global index in EXPLICIT global_index_table', status)
               return
            end if
         end do

      case default
         call htg_fail(HTG_ERR_INDEXING, &
            'htg_validate_tree: unrecognized indexing_mode', status)
         return
      end select

      call htg_set_ok(status)
   end subroutine htg_validate_tree

   ! --- Comprehensive Whole-Grid Structural Invariant Audit ---

   subroutine htg_validate_grid(grid, status, allow_unset_global_index)
      type(htg_grid_t), intent(in) :: grid
      type(htg_status_t), intent(out), optional :: status
      logical, intent(in), optional :: allow_unset_global_index

      logical :: allow_unset, any_gid_set, all_gid_unset
      integer(i32_k) :: exp_cdims(0:2), exp_ndim, exp_orient, exp_axes(0:2), exp_nchild
      integer(i32_k) :: arr_idx, nc
      integer(id_k) :: max_trees, pos, prev_tid, tid, nv, v
      integer(id_k) :: num_imp, num_exp, imp_pos, exp_pos, max_gid, t_max
      integer(id_k) :: pm_bytes
      integer(id_k), allocatable :: imp_lo(:), imp_hi(:), exp_gids(:)
      type(htg_bit_array_t) :: expected_pm
      character(len=:), allocatable :: norm_name, inter_name, act_scalars
      logical :: found_norm, found_inter
      integer(i32_k) :: norm_nc, inter_nc
      type(htg_status_t) :: st

      allow_unset = .false.
      if (present(allow_unset_global_index)) allow_unset = allow_unset_global_index

      ! 1. Dimensionality, BranchFactor & Squeezed Geometry Consistency
      if (.not. grid%dimensions_set) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_validate_grid: grid dimensions are not set', status)
         return
      end if

      if (grid%branch_factor /= 2_i32_k .and. grid%branch_factor /= 3_i32_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_validate_grid: branch_factor must be 2 or 3', status)
         return
      end if

      if (any(grid%dimensions < 1_i32_k)) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_validate_grid: point dimensions must all be >= 1', status)
         return
      end if

      if (grid%depth_limiter < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_validate_grid: negative depth_limiter', status)
         return
      end if

      call htg_compute_squeezed_geometry(grid%dimensions, grid%branch_factor, &
         exp_cdims, exp_ndim, exp_orient, exp_axes, exp_nchild)

      if (any(grid%cell_dims /= exp_cdims) .or. &
          grid%dimension /= exp_ndim .or. &
          grid%orientation /= exp_orient .or. &
          any(grid%axes /= exp_axes) .or. &
          grid%number_of_children /= exp_nchild) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_validate_grid: squeezed geometry metadata inconsistent with dimensions', &
            status)
         return
      end if

      ! 2. Rectilinear Coordinates Consistency & Strict Monotonicity
      call validate_axis_coordinates(grid%x_coordinates, grid%has_x_coordinates, &
         grid%dimensions(0), 'X', st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      call validate_axis_coordinates(grid%y_coordinates, grid%has_y_coordinates, &
         grid%dimensions(1), 'Y', st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      call validate_axis_coordinates(grid%z_coordinates, grid%has_z_coordinates, &
         grid%dimensions(2), 'Z', st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      ! 3. Sparse Forest Container & Per-Tree Topological Audit
      max_trees = grid%get_max_number_of_trees()
      if (grid%num_trees < 0_id_k .or. grid%num_trees > max_trees) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_grid: num_trees out of range [0, max_number_of_trees]', status)
         return
      end if

      if (grid%dimension == 0_i32_k .and. grid%num_trees > 0_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_validate_grid: 0D (1,1,1) grid cannot contain active HyperTrees', status)
         return
      end if

      if (grid%num_trees > 0_id_k) then
         if (.not. allocated(grid%tree_slots) .or. &
             lbound(grid%tree_slots, 1) /= 0 .or. &
             size(grid%tree_slots, kind=id_k) < grid%num_trees) then
            call htg_fail(HTG_ERR_INVALID_STATE, &
               'htg_validate_grid: tree_slots unallocated or mis-sized for num_trees', status)
            return
         end if
      end if

      any_gid_set = .false.
      all_gid_unset = .true.
      num_imp = 0_id_k
      num_exp = 0_id_k
      max_gid = HTG_INVALID_ID
      prev_tid = HTG_INVALID_ID

      do pos = 0_id_k, grid%num_trees - 1_id_k
         tid = grid%tree_slots(pos)%tree_id
         if (tid < 0_id_k .or. tid >= max_trees) then
            call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
               'htg_validate_grid: tree_slots contains out-of-bounds tree_id', status)
            return
         end if
         if (pos > 0_id_k .and. tid <= prev_tid) then
            call htg_fail(HTG_ERR_INVALID_STATE, &
               'htg_validate_grid: tree_slots are not strictly sorted by ascending tree_id', status)
            return
         end if
         prev_tid = tid

         if (.not. allocated(grid%tree_slots(pos)%tree)) then
            call htg_fail(HTG_ERR_INVALID_STATE, &
               'htg_validate_grid: unallocated tree inside active tree_slot', status)
            return
         end if

         if (grid%tree_slots(pos)%tree%get_tree_index() /= tid) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_validate_grid: tree%tree_index disagrees with slot tree_id', status)
            return
         end if

         if (grid%tree_slots(pos)%tree%get_branch_factor() /= grid%branch_factor .or. &
             grid%tree_slots(pos)%tree%get_dimension() /= grid%dimension) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_validate_grid: tree (branch_factor, dimension) disagrees with grid', status)
            return
         end if

         call htg_validate_tree(grid%tree_slots(pos)%tree, st, &
            require_global_index = (.not. allow_unset))
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if

         select case (grid%tree_slots(pos)%tree%get_indexing_mode())
         case (HTG_INDEX_UNSET)
            continue
         case (HTG_INDEX_IMPLICIT)
            any_gid_set = .true.
            all_gid_unset = .false.
            num_imp = num_imp + 1_id_k
            t_max = grid%tree_slots(pos)%tree%get_max_global_index()
            if (t_max > max_gid) max_gid = t_max
         case (HTG_INDEX_EXPLICIT)
            any_gid_set = .true.
            all_gid_unset = .false.
            num_exp = num_exp + grid%tree_slots(pos)%tree%get_number_of_vertices()
            t_max = grid%tree_slots(pos)%tree%get_max_global_index()
            if (t_max > max_gid) max_gid = t_max
         end select
      end do

      if (any_gid_set .and. .not. all_gid_unset) then
         do pos = 0_id_k, grid%num_trees - 1_id_k
            if (grid%tree_slots(pos)%tree%get_indexing_mode() == HTG_INDEX_UNSET) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_validate_grid: mixed unset and set global indexing across active trees', status)
               return
            end if
         end do
      end if

      ! 4. Forest-Wide Global Index Completeness & Zero Duplicate/Overlapping IDs
      if (num_imp > 0_id_k) then
         allocate(imp_lo(0_id_k : num_imp - 1_id_k))
         allocate(imp_hi(0_id_k : num_imp - 1_id_k))
         imp_pos = 0_id_k
         do pos = 0_id_k, grid%num_trees - 1_id_k
            if (grid%tree_slots(pos)%tree%get_indexing_mode() == HTG_INDEX_IMPLICIT) then
               nv = grid%tree_slots(pos)%tree%get_number_of_vertices()
               imp_lo(imp_pos) = grid%tree_slots(pos)%tree%get_global_index_start()
               imp_hi(imp_pos) = imp_lo(imp_pos) + nv - 1_id_k
               imp_pos = imp_pos + 1_id_k
            end if
         end do

         call heapsort_intervals(imp_lo, imp_hi, num_imp)
         do imp_pos = 1_id_k, num_imp - 1_id_k
            if (imp_lo(imp_pos) <= imp_hi(imp_pos - 1_id_k)) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_validate_grid: overlapping implicit global index intervals across trees', status)
               return
            end if
         end do
      end if

      if (num_exp > 0_id_k) then
         allocate(exp_gids(0_id_k : num_exp - 1_id_k))
         exp_pos = 0_id_k
         do pos = 0_id_k, grid%num_trees - 1_id_k
            if (grid%tree_slots(pos)%tree%get_indexing_mode() == HTG_INDEX_EXPLICIT) then
               nv = grid%tree_slots(pos)%tree%get_number_of_vertices()
               do v = 0_id_k, nv - 1_id_k
                  exp_gids(exp_pos) = grid%tree_slots(pos)%tree%global_index_table(v)
                  exp_pos = exp_pos + 1_id_k
               end do
            end if
         end do

         call heapsort_i64(exp_gids, num_exp)
         do exp_pos = 1_id_k, num_exp - 1_id_k
            if (exp_gids(exp_pos) == exp_gids(exp_pos - 1_id_k)) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_validate_grid: duplicate explicit global ID across trees', status)
               return
            end if
         end do

         ! Cross-check sorted explicit global IDs against sorted disjoint implicit intervals
         if (num_imp > 0_id_k) then
            imp_pos = 0_id_k
            exp_pos = 0_id_k
            do while (imp_pos < num_imp .and. exp_pos < num_exp)
               if (exp_gids(exp_pos) < imp_lo(imp_pos)) then
                  exp_pos = exp_pos + 1_id_k
               else if (exp_gids(exp_pos) > imp_hi(imp_pos)) then
                  imp_pos = imp_pos + 1_id_k
               else
                  call htg_fail(HTG_ERR_INDEXING, &
                     'htg_validate_grid: explicit global ID collides with implicit tree interval', status)
                  return
               end if
            end do
         end if
      end if

      ! 5. CellData Tuple Sufficiency & Structural Integrity
      if (grid%cell_data%num_arrays < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_validate_grid: negative cell_data%num_arrays', status)
         return
      end if

      if (grid%cell_data%num_arrays > 0_i32_k) then
         if (.not. allocated(grid%cell_data%arrays) .or. &
             lbound(grid%cell_data%arrays, 1) /= 0 .or. &
             size(grid%cell_data%arrays, kind=i32_k) < grid%cell_data%num_arrays) then
            call htg_fail(HTG_ERR_INVALID_STATE, &
               'htg_validate_grid: cell_data%arrays container inconsistent with num_arrays', &
               status)
            return
         end if

         if (grid%num_trees > 0_id_k .and. max_gid < 0_id_k) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_validate_grid: cell_data arrays present while trees lack global indices', &
               status)
            return
         end if

         do arr_idx = 0_i32_k, grid%cell_data%num_arrays - 1_i32_k
            if (.not. allocated(grid%cell_data%arrays(arr_idx)%item)) then
               call htg_fail(HTG_ERR_INVALID_STATE, &
                  'htg_validate_grid: unallocated item in cell_data%arrays', status)
               return
            end if
            associate (arr => grid%cell_data%arrays(arr_idx)%item)
               if (.not. arr%is_initialized() .or. &
                   arr%get_number_of_components() < 1_i32_k) then
                  call htg_fail(HTG_ERR_INVALID_STATE, &
                     'htg_validate_grid: invalid cell_data field array metadata', status)
                  return
               end if
               if (max_gid >= 0_id_k .and. &
                   arr%get_number_of_tuples() - 1_id_k < max_gid) then
                  call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
                     'htg_validate_grid: cell_data array "' // trim(arr%get_name()) // &
                     '" has fewer tuples than max_global_index + 1', status)
                  return
               end if
               if (arr%get_data_type() == HTG_TYPE_BIT) then
                  call validate_bit_array_internal(arr%bit_values, &
                     HTG_ERR_INVALID_STATE, 'htg_validate_grid(cell_data bit array)', st)
                  if (.not. st%ok()) then
                     call htg_fail(st%code, trim(st%message), status)
                     return
                  end if
               end if
            end associate
         end do
      end if

      if (grid%cell_data%has_active_scalars()) then
         act_scalars = grid%cell_data%get_active_scalars()
         if (.not. grid%cell_data%has_array(act_scalars)) then
            call htg_fail(HTG_ERR_INVALID_STATE, &
               'htg_validate_grid: active_scalars array not found in cell_data', status)
            return
         end if
      end if

      ! 6. Material Interface Consistency
      if (grid%get_has_interface()) then
         norm_name = grid%get_interface_normals_name()
         inter_name = grid%get_interface_intercepts_name()

         if (len_trim(norm_name) == 0) then
            call htg_fail(HTG_ERR_INTERFACE, &
               'htg_validate_grid: has_interface true with empty interface_normals_name', &
               status)
            return
         end if

         if (len_trim(inter_name) == 0) then
            call htg_fail(HTG_ERR_INTERFACE, &
               'htg_validate_grid: has_interface true with empty interface_intercepts_name', &
               status)
            return
         end if

         found_norm = .false.
         found_inter = .false.
         norm_nc = 0_i32_k
         inter_nc = 0_i32_k

         if (allocated(grid%cell_data%arrays)) then
            do arr_idx = 0_i32_k, grid%cell_data%num_arrays - 1_i32_k
               if (allocated(grid%cell_data%arrays(arr_idx)%item)) then
                  nc = grid%cell_data%arrays(arr_idx)%item%get_number_of_components()
                  if (grid%cell_data%arrays(arr_idx)%item%get_name() == &
                      trim(norm_name)) then
                     found_norm = .true.
                     norm_nc = nc
                  end if
                  if (grid%cell_data%arrays(arr_idx)%item%get_name() == &
                      trim(inter_name)) then
                     found_inter = .true.
                     inter_nc = nc
                  end if
               end if
            end do
         end if

         if (.not. found_norm) then
            call htg_fail(HTG_ERR_INTERFACE, &
               'htg_validate_grid: interface normals array "' // trim(norm_name) // &
               '" missing from cell_data', status)
            return
         end if

         if (norm_nc /= 3_i32_k) then
            call htg_fail(HTG_ERR_INTERFACE, &
               'htg_validate_grid: interface normals array must have 3 components', &
               status)
            return
         end if

         if (.not. found_inter) then
            call htg_fail(HTG_ERR_INTERFACE, &
               'htg_validate_grid: interface intercepts array "' // trim(inter_name) // &
               '" missing from cell_data', status)
            return
         end if

         if (inter_nc /= 3_i32_k) then
            call htg_fail(HTG_ERR_INTERFACE, &
               'htg_validate_grid: interface intercepts array must have 3 components', &
               status)
            return
         end if
      end if

      ! 7. Mask & PureMask Consistency
      if (grid%mask_active) then
         if (grid%has_mask() .and. grid%num_trees > 0_id_k .and. max_gid < 0_id_k) then
            call htg_fail(HTG_ERR_MASK, &
               'htg_validate_grid: mask is active while trees lack global indices', &
               status)
            return
         end if
         call validate_bit_array_internal(grid%mask, HTG_ERR_MASK, &
            'htg_validate_grid(mask)', st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
      end if

      if (grid%pure_mask_valid) then
         call validate_bit_array_internal(grid%pure_mask, HTG_ERR_MASK, &
            'htg_validate_grid(pure_mask)', st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if

         call grid%compute_pure_mask_into(expected_pm, st)
         if (.not. st%ok()) then
            call htg_fail(HTG_ERR_MASK, &
               'htg_validate_grid: failed recomputing PureMask: ' // trim(st%message), &
               status)
            return
         end if

         if (expected_pm%get_number_of_values() /= &
             grid%pure_mask%get_number_of_values()) then
            call htg_fail(HTG_ERR_MASK, &
               'htg_validate_grid: cached pure_mask bit length disagrees with PureMask', &
               status)
            return
         end if

         pm_bytes = expected_pm%get_number_of_bytes()
         if (pm_bytes > 0_id_k) then
            if (any(expected_pm%bytes(0_id_k : pm_bytes - 1_id_k) /= &
                    grid%pure_mask%bytes(0_id_k : pm_bytes - 1_id_k))) then
               call htg_fail(HTG_ERR_MASK, &
                  'htg_validate_grid: cached pure_mask bits disagree with PureMask', &
                  status)
               return
            end if
         end if
      end if

      call htg_set_ok(status)
   end subroutine htg_validate_grid

end module fvtkhtg_validate_m

submodule (fvtkhtg_grid_m) fvtkhtg_grid_validate_s
   use fvtkhtg_validate_m, only: htg_validate_grid
   implicit none
contains
   module procedure grid_validate
      select type (self)
      type is (htg_grid_t)
         call htg_validate_grid(self, status)
      class default
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_grid_t%validate: unsupported polymorphic dynamic type', status)
      end select
   end procedure grid_validate
end submodule fvtkhtg_grid_validate_s
