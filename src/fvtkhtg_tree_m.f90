module fvtkhtg_tree_m
   use fvtkhtg_kinds_m, only: i32_k, id_k, htg_status_t, &
      HTG_ERR_INDEXING, HTG_ERR_INVALID_ARG, HTG_ERR_INVALID_STATE, &
      HTG_ERR_MASK, HTG_ERR_OUT_OF_BOUNDS, HTG_ERR_TOPOLOGY, &
      HTG_INDEX_EXPLICIT, HTG_INDEX_IMPLICIT, HTG_INDEX_UNSET, &
      HTG_INVALID_ID, HTG_UINT32_MAX, htg_fail, htg_set_ok
   use fvtkhtg_bit_array_m, only: htg_bit_array_t
   implicit none
   private

   integer(id_k), parameter, public :: HTG_TREE_LEAF_SENTINEL = HTG_UINT32_MAX

   integer(id_k), parameter :: INITIAL_DEPTH_CAPACITY = 16_id_k
   integer(id_k), parameter :: INITIAL_ELDER_CAPACITY = 16_id_k
   integer(id_k), parameter :: INITIAL_GLOBAL_CAPACITY = 16_id_k
   integer(id_k), parameter :: INITIAL_BLOCK_CAPACITY = 16_id_k

   !> What htg_tree_t%coarsen_node did to the tree's local numbering.
   !>
   !> Local vertex ids stay contiguous: removing the children of a node leaves
   !> a gap of number_of_children ids, which is closed by relocating the
   !> sibling block with the highest ids into it. No other vertex is renumbered.
   type, public :: htg_coarsen_info_t
      integer(i32_k) :: number_of_children = 0_i32_k
      !> Local id of the coarsened node after the operation (it changes only
      !> if the node itself was in the relocated block).
      integer(id_k) :: node_id = HTG_INVALID_ID
      !> First local id of the removed child block, as numbered before.
      integer(id_k) :: removed_first = HTG_INVALID_ID
      !> First local id, before and after, of the relocated sibling block;
      !> both HTG_INVALID_ID when the removed block was already the last one.
      integer(id_k) :: moved_from = HTG_INVALID_ID
      integer(id_k) :: moved_to = HTG_INVALID_ID
   end type htg_coarsen_info_t

   type, public :: htg_tree_t
      logical :: initialized = .false.
      integer(id_k) :: tree_index = HTG_INVALID_ID
      integer(i32_k) :: branch_factor = 2_i32_k
      integer(i32_k) :: dimension = 3_i32_k
      integer(i32_k) :: number_of_children = 8_i32_k

      integer(id_k) :: number_of_levels = 0_id_k
      integer(id_k) :: number_of_vertices = 0_id_k
      integer(id_k) :: number_of_nodes = 0_id_k
      integer(id_k) :: number_of_leaves = 0_id_k

      integer(id_k) :: depth_capacity = 0_id_k
      integer(id_k), allocatable :: number_of_vertices_per_depth(:)

      integer(id_k) :: elder_active_len = 0_id_k
      integer(id_k) :: elder_capacity = 0_id_k
      integer(id_k), allocatable :: parent_to_elder_child(:)

      ! Parent of each sibling block. Children are stored in blocks of
      ! number_of_children consecutive local ids: block b holds vertices
      ! 1 + b*k .. (b + 1)*k and block_parent(b) is the local id of their
      ! parent. Entries 0 .. number_of_nodes - 1 are in use. VTK keeps no such
      ! table; it is what makes parent lookup and coarsening O(1).
      integer(id_k) :: block_parent_capacity = 0_id_k
      integer(id_k), allocatable :: block_parent(:)

      integer(i32_k) :: indexing_mode = HTG_INDEX_UNSET
      integer(id_k) :: global_index_start = HTG_INVALID_ID
      integer(id_k) :: global_index_capacity = 0_id_k
      integer(id_k), allocatable :: global_index_table(:)
   contains
      procedure :: initialize => tree_initialize
      procedure :: clear => tree_clear
      procedure :: is_initialized => tree_is_initialized

      procedure :: set_tree_index => tree_set_tree_index
      procedure :: get_tree_index => tree_get_tree_index
      procedure :: get_branch_factor => tree_get_branch_factor
      procedure :: get_dimension => tree_get_dimension
      procedure :: get_number_of_children => tree_get_number_of_children
      procedure :: get_number_of_levels => tree_get_number_of_levels
      procedure :: get_number_of_vertices => tree_get_number_of_vertices
      procedure :: get_number_of_nodes => tree_get_number_of_nodes
      procedure :: get_number_of_leaves => tree_get_number_of_leaves
      procedure :: get_elder_child_array_size => tree_get_elder_child_array_size

      procedure, private :: tree_get_vertices_at_depth_i64
      procedure, private :: tree_get_vertices_at_depth_i32
      generic :: get_number_of_vertices_at_depth => &
         tree_get_vertices_at_depth_i64, tree_get_vertices_at_depth_i32

      procedure :: is_leaf => tree_is_leaf
      procedure :: is_terminal_node => tree_is_terminal_node
      procedure :: get_elder_child_index => tree_get_elder_child_index
      procedure :: get_parent_index => tree_get_parent_index
      procedure :: get_child_index => tree_get_child_index
      procedure :: get_level => tree_get_level

      procedure, private :: tree_subdivide_leaf
      procedure, private :: tree_subdivide_leaf_i64
      generic :: subdivide_leaf => tree_subdivide_leaf, tree_subdivide_leaf_i64
      procedure :: coarsen_node => tree_coarsen_node

      procedure :: set_global_index_start => tree_set_global_index_start
      procedure :: get_global_index_start => tree_get_global_index_start
      procedure :: set_global_index_from_local => tree_set_global_index_from_local
      procedure :: get_global_index_from_local => tree_get_global_index_from_local
      procedure :: get_indexing_mode => tree_get_indexing_mode
      procedure :: get_max_global_index => tree_get_max_global_index

      procedure :: build_from_breadth_first_order_descriptor => &
         tree_build_from_breadth_first_order_descriptor
      procedure :: compute_breadth_first_order_descriptor => &
         tree_compute_breadth_first_order_descriptor

      procedure :: get_parent_to_elder_child_view => tree_get_parent_to_elder_child_view
      procedure :: get_vertices_per_depth_view => tree_get_vertices_per_depth_view
      procedure :: get_global_index_table_view => tree_get_global_index_table_view
   end type htg_tree_t

contains

   pure subroutine ensure_depth_capacity(self, min_levels)
      type(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: min_levels

      integer(id_k) :: new_cap
      integer(id_k), allocatable :: tmp(:)

      if (min_levels <= self%depth_capacity .and. allocated(self%number_of_vertices_per_depth)) then
         return
      end if

      new_cap = max(min_levels, max(INITIAL_DEPTH_CAPACITY, 2_id_k * self%depth_capacity))
      allocate(tmp(0_id_k : new_cap - 1_id_k))
      tmp = 0_id_k
      if (self%number_of_levels > 0_id_k .and. allocated(self%number_of_vertices_per_depth)) then
         tmp(0_id_k : self%number_of_levels - 1_id_k) = &
            self%number_of_vertices_per_depth(0_id_k : self%number_of_levels - 1_id_k)
      end if
      call move_alloc(tmp, self%number_of_vertices_per_depth)
      self%depth_capacity = new_cap
   end subroutine ensure_depth_capacity

   pure subroutine ensure_elder_capacity(self, min_len)
      type(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: min_len

      integer(id_k) :: new_cap
      integer(id_k), allocatable :: tmp(:)

      if (min_len <= self%elder_capacity .and. allocated(self%parent_to_elder_child)) then
         return
      end if

      new_cap = max(min_len, max(INITIAL_ELDER_CAPACITY, 2_id_k * self%elder_capacity))
      allocate(tmp(0_id_k : new_cap - 1_id_k))
      tmp = HTG_TREE_LEAF_SENTINEL
      if (self%elder_active_len > 0_id_k .and. allocated(self%parent_to_elder_child)) then
         tmp(0_id_k : self%elder_active_len - 1_id_k) = &
            self%parent_to_elder_child(0_id_k : self%elder_active_len - 1_id_k)
      end if
      call move_alloc(tmp, self%parent_to_elder_child)
      self%elder_capacity = new_cap
   end subroutine ensure_elder_capacity

   pure subroutine ensure_block_parent_capacity(self, min_blocks)
      type(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: min_blocks

      integer(id_k) :: new_cap
      integer(id_k), allocatable :: tmp(:)

      if (min_blocks <= self%block_parent_capacity .and. allocated(self%block_parent)) return

      new_cap = max(min_blocks, max(INITIAL_BLOCK_CAPACITY, 2_id_k * self%block_parent_capacity))
      allocate(tmp(0_id_k : new_cap - 1_id_k))
      tmp = HTG_INVALID_ID
      if (self%number_of_nodes > 0_id_k .and. allocated(self%block_parent)) then
         tmp(0_id_k : self%number_of_nodes - 1_id_k) = &
            self%block_parent(0_id_k : self%number_of_nodes - 1_id_k)
      end if
      call move_alloc(tmp, self%block_parent)
      self%block_parent_capacity = new_cap
   end subroutine ensure_block_parent_capacity

   pure subroutine ensure_global_table_capacity(self, min_vertices, old_valid_len)
      type(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: min_vertices
      integer(id_k), intent(in) :: old_valid_len

      integer(id_k) :: new_cap
      integer(id_k), allocatable :: tmp(:)

      if (min_vertices <= self%global_index_capacity .and. allocated(self%global_index_table)) then
         if (min_vertices > old_valid_len) then
            self%global_index_table(old_valid_len : min_vertices - 1_id_k) = HTG_INVALID_ID
         end if
         return
      end if

      new_cap = max(min_vertices, max(INITIAL_GLOBAL_CAPACITY, 2_id_k * self%global_index_capacity))
      allocate(tmp(0_id_k : new_cap - 1_id_k))
      tmp = HTG_INVALID_ID
      if (old_valid_len > 0_id_k .and. allocated(self%global_index_table)) then
         tmp(0_id_k : old_valid_len - 1_id_k) = &
            self%global_index_table(0_id_k : old_valid_len - 1_id_k)
      end if
      call move_alloc(tmp, self%global_index_table)
      self%global_index_capacity = new_cap
   end subroutine ensure_global_table_capacity

   pure logical function is_leaf_unchecked(self, local_id) result(flag)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: local_id

      if (self%number_of_vertices == 1_id_k) then
         flag = .true.
      else if (local_id >= self%elder_active_len) then
         flag = .true.
      else
         flag = (self%parent_to_elder_child(local_id) == HTG_TREE_LEAF_SENTINEL .or. &
                 self%parent_to_elder_child(local_id) < 0_id_k)
      end if
   end function is_leaf_unchecked

   pure subroutine tree_clear(self)
      class(htg_tree_t), intent(inout) :: self

      if (allocated(self%number_of_vertices_per_depth)) then
         deallocate(self%number_of_vertices_per_depth)
      end if
      if (allocated(self%parent_to_elder_child)) then
         deallocate(self%parent_to_elder_child)
      end if
      if (allocated(self%global_index_table)) then
         deallocate(self%global_index_table)
      end if
      if (allocated(self%block_parent)) then
         deallocate(self%block_parent)
      end if
      self%block_parent_capacity = 0_id_k

      self%initialized = .false.
      self%tree_index = HTG_INVALID_ID
      self%branch_factor = 2_i32_k
      self%dimension = 3_i32_k
      self%number_of_children = 8_i32_k

      self%number_of_levels = 0_id_k
      self%number_of_vertices = 0_id_k
      self%number_of_nodes = 0_id_k
      self%number_of_leaves = 0_id_k

      self%depth_capacity = 0_id_k
      self%elder_active_len = 0_id_k
      self%elder_capacity = 0_id_k

      self%indexing_mode = HTG_INDEX_UNSET
      self%global_index_start = HTG_INVALID_ID
      self%global_index_capacity = 0_id_k
   end subroutine tree_clear

   pure logical function tree_is_initialized(self) result(flag)
      class(htg_tree_t), intent(in) :: self
      flag = self%initialized
   end function tree_is_initialized

   pure subroutine tree_initialize(self, branch_factor, dimension, tree_index, status)
      class(htg_tree_t), intent(inout) :: self
      integer(i32_k), intent(in) :: branch_factor
      integer(i32_k), intent(in) :: dimension
      integer(id_k), intent(in), optional :: tree_index
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: tid

      if (branch_factor < 2_i32_k .or. branch_factor > 3_i32_k) then
         call self%clear()
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_tree_t%initialize: branch_factor must be 2 or 3', status)
         return
      end if

      if (dimension < 1_i32_k .or. dimension > 3_i32_k) then
         call self%clear()
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_tree_t%initialize: dimension must be 1, 2, or 3', status)
         return
      end if

      tid = HTG_INVALID_ID
      if (present(tree_index)) then
         if (tree_index < HTG_INVALID_ID) then
            call self%clear()
            call htg_fail(HTG_ERR_INVALID_ARG, &
               'htg_tree_t%initialize: tree_index must be >= -1', status)
            return
         end if
         tid = tree_index
      end if

      call self%clear()

      self%initialized = .true.
      self%tree_index = tid
      self%branch_factor = branch_factor
      self%dimension = dimension
      self%number_of_children = branch_factor ** dimension

      self%number_of_levels = 1_id_k
      self%number_of_vertices = 1_id_k
      self%number_of_nodes = 0_id_k
      self%number_of_leaves = 1_id_k

      call ensure_depth_capacity(self, 1_id_k)
      self%number_of_vertices_per_depth(0_id_k) = 1_id_k

      call ensure_elder_capacity(self, 1_id_k)
      self%parent_to_elder_child(0_id_k) = HTG_TREE_LEAF_SENTINEL
      self%elder_active_len = 1_id_k

      self%indexing_mode = HTG_INDEX_UNSET
      self%global_index_start = HTG_INVALID_ID

      call htg_set_ok(status)
   end subroutine tree_initialize

   pure subroutine tree_set_tree_index(self, tree_index, status)
      class(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: tree_index
      type(htg_status_t), intent(out), optional :: status

      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%set_tree_index: tree is not initialized', status)
         return
      end if

      if (tree_index < HTG_INVALID_ID) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_tree_t%set_tree_index: tree_index must be >= -1', status)
         return
      end if

      self%tree_index = tree_index
      call htg_set_ok(status)
   end subroutine tree_set_tree_index

   pure function tree_get_tree_index(self) result(tid)
      class(htg_tree_t), intent(in) :: self
      integer(id_k) :: tid
      tid = self%tree_index
   end function tree_get_tree_index

   pure function tree_get_branch_factor(self) result(bf)
      class(htg_tree_t), intent(in) :: self
      integer(i32_k) :: bf
      bf = self%branch_factor
   end function tree_get_branch_factor

   pure function tree_get_dimension(self) result(dim)
      class(htg_tree_t), intent(in) :: self
      integer(i32_k) :: dim
      dim = self%dimension
   end function tree_get_dimension

   pure function tree_get_number_of_children(self) result(nchildren)
      class(htg_tree_t), intent(in) :: self
      integer(i32_k) :: nchildren
      nchildren = self%number_of_children
   end function tree_get_number_of_children

   pure function tree_get_number_of_levels(self) result(nlevels)
      class(htg_tree_t), intent(in) :: self
      integer(id_k) :: nlevels
      nlevels = self%number_of_levels
   end function tree_get_number_of_levels

   pure function tree_get_number_of_vertices(self) result(nverts)
      class(htg_tree_t), intent(in) :: self
      integer(id_k) :: nverts
      nverts = self%number_of_vertices
   end function tree_get_number_of_vertices

   pure function tree_get_number_of_nodes(self) result(nnodes)
      class(htg_tree_t), intent(in) :: self
      integer(id_k) :: nnodes
      nnodes = self%number_of_nodes
   end function tree_get_number_of_nodes

   pure function tree_get_number_of_leaves(self) result(nleaves)
      class(htg_tree_t), intent(in) :: self
      integer(id_k) :: nleaves
      nleaves = self%number_of_leaves
   end function tree_get_number_of_leaves

   pure function tree_get_elder_child_array_size(self) result(n)
      class(htg_tree_t), intent(in) :: self
      integer(id_k) :: n
      n = self%elder_active_len
   end function tree_get_elder_child_array_size

   function tree_get_vertices_at_depth_i64(self, depth, status) result(nverts)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: depth
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: nverts

      nverts = 0_id_k
      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%get_number_of_vertices_at_depth: tree is not initialized', status)
         return
      end if

      if (depth < 0_id_k .or. depth >= self%number_of_levels) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%get_number_of_vertices_at_depth: depth out of bounds', status)
         return
      end if

      nverts = self%number_of_vertices_per_depth(depth)
      call htg_set_ok(status)
   end function tree_get_vertices_at_depth_i64

   function tree_get_vertices_at_depth_i32(self, depth, status) result(nverts)
      class(htg_tree_t), intent(in) :: self
      integer(i32_k), intent(in) :: depth
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: nverts

      nverts = self%tree_get_vertices_at_depth_i64(int(depth, id_k), status)
   end function tree_get_vertices_at_depth_i32

   logical function tree_is_leaf(self, local_id, status) result(flag)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: local_id
      type(htg_status_t), intent(out), optional :: status

      flag = .false.
      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%is_leaf: tree is not initialized', status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%is_leaf: local_id out of bounds', status)
         return
      end if

      flag = is_leaf_unchecked(self, local_id)
      call htg_set_ok(status)
   end function tree_is_leaf

   logical function tree_is_terminal_node(self, local_id, status) result(flag)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: local_id
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: elder, child_id
      integer(i32_k) :: ichild

      flag = .false.
      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%is_terminal_node: tree is not initialized', status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%is_terminal_node: local_id out of bounds', status)
         return
      end if

      if (is_leaf_unchecked(self, local_id)) then
         flag = .false.
         call htg_set_ok(status)
         return
      end if

      elder = self%parent_to_elder_child(local_id)
      do ichild = 0_i32_k, self%number_of_children - 1_i32_k
         child_id = elder + int(ichild, id_k)
         if (.not. is_leaf_unchecked(self, child_id)) then
            flag = .false.
            call htg_set_ok(status)
            return
         end if
      end do

      flag = .true.
      call htg_set_ok(status)
   end function tree_is_terminal_node

   function tree_get_elder_child_index(self, local_id, status) result(elder)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: local_id
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: elder

      elder = HTG_TREE_LEAF_SENTINEL
      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%get_elder_child_index: tree is not initialized', status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%get_elder_child_index: local_id out of bounds', status)
         return
      end if

      if (is_leaf_unchecked(self, local_id)) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_tree_t%get_elder_child_index: vertex is a leaf', status)
         return
      end if

      elder = self%parent_to_elder_child(local_id)
      call htg_set_ok(status)
   end function tree_get_elder_child_index

   !> Local id of the parent of a vertex; HTG_INVALID_ID (with an OK status)
   !> for the root. O(1).
   function tree_get_parent_index(self, local_id, status) result(parent)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: local_id
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: parent

      parent = HTG_INVALID_ID
      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%get_parent_index: tree is not initialized', status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%get_parent_index: local_id out of bounds', status)
         return
      end if

      if (local_id > 0_id_k) then
         parent = self%block_parent((local_id - 1_id_k) / int(self%number_of_children, id_k))
      end if
      call htg_set_ok(status)
   end function tree_get_parent_index

   !> Which child of its parent a vertex is (0 .. number_of_children - 1);
   !> -1 (with an OK status) for the root. O(1).
   function tree_get_child_index(self, local_id, status) result(ichild)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: local_id
      type(htg_status_t), intent(out), optional :: status
      integer(i32_k) :: ichild

      ichild = -1_i32_k
      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%get_child_index: tree is not initialized', status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%get_child_index: local_id out of bounds', status)
         return
      end if

      if (local_id > 0_id_k) then
         ichild = int(mod(local_id - 1_id_k, int(self%number_of_children, id_k)), i32_k)
      end if
      call htg_set_ok(status)
   end function tree_get_child_index

   !> Depth of a vertex (0 for the root). O(depth).
   pure subroutine level_unchecked(self, local_id, level, ok)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: local_id
      integer(id_k), intent(out) :: level
      logical, intent(out) :: ok

      integer(id_k) :: v, k

      k = int(self%number_of_children, id_k)
      v = local_id
      level = 0_id_k
      ok = .true.
      do while (v /= 0_id_k)
         ! A vertex is never deeper than number_of_levels - 1; anything else
         ! means the parent table is corrupt (and would loop forever).
         if (level >= self%number_of_levels .or. v < 0_id_k .or. &
             v >= self%number_of_vertices) then
            ok = .false.
            return
         end if
         v = self%block_parent((v - 1_id_k) / k)
         level = level + 1_id_k
      end do
   end subroutine level_unchecked

   function tree_get_level(self, local_id, status) result(level)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: local_id
      type(htg_status_t), intent(out), optional :: status
      integer(i32_k) :: level

      integer(id_k) :: lvl
      logical :: ok

      level = -1_i32_k
      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%get_level: tree is not initialized', status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%get_level: local_id out of bounds', status)
         return
      end if

      call level_unchecked(self, local_id, lvl, ok)
      if (.not. ok) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_tree_t%get_level: parent table is inconsistent', status)
         return
      end if
      level = int(lvl, i32_k)
      call htg_set_ok(status)
   end function tree_get_level

   !> Inverse of subdivide_leaf: remove the children of a refined vertex whose
   !> children are all leaves, turning it back into a leaf. VTK has no such
   !> operation.
   !>
   !> Local ids stay contiguous. The freed block of number_of_children ids is
   !> filled by relocating the sibling block with the highest ids (unless the
   !> freed block is that block); `info` reports the relocation. In EXPLICIT
   !> indexing mode global ids travel with their vertices and the removed
   !> children's entries are dropped. In IMPLICIT mode the relocated vertices
   !> get new global ids (start + local id), so the caller must move their
   !> cell data; htg_cursor_t%coarsen_node does that.
   pure subroutine tree_coarsen_node(self, local_id, info, status)
      class(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: local_id
      type(htg_coarsen_info_t), intent(out) :: info
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: k, elder, nv_old, last_first, src, dst, moved_elder
      integer(id_k) :: block_parent_of_last, child_depth, i
      logical :: level_ok, has_table

      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%coarsen_node: tree is not initialized', status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%coarsen_node: local_id out of bounds', status)
         return
      end if

      if (is_leaf_unchecked(self, local_id)) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_tree_t%coarsen_node: vertex is a leaf (nothing to coarsen)', status)
         return
      end if

      k = int(self%number_of_children, id_k)
      elder = self%parent_to_elder_child(local_id)
      do i = 0_id_k, k - 1_id_k
         if (.not. is_leaf_unchecked(self, elder + i)) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_tree_t%coarsen_node: children are not all leaves (coarsen bottom-up)', &
               status)
            return
         end if
      end do

      call level_unchecked(self, local_id, child_depth, level_ok)
      if (.not. level_ok) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_tree_t%coarsen_node: parent table is inconsistent', status)
         return
      end if
      child_depth = child_depth + 1_id_k

      nv_old = self%number_of_vertices
      last_first = nv_old - k
      has_table = (self%indexing_mode == HTG_INDEX_EXPLICIT) .and. &
         allocated(self%global_index_table)

      info%number_of_children = self%number_of_children
      info%node_id = local_id
      info%removed_first = elder

      self%parent_to_elder_child(local_id) = HTG_TREE_LEAF_SENTINEL

      if (elder /= last_first) then
         ! Relocate the last sibling block into the gap. Its parent is neither
         ! local_id (whose children are the removed block) nor one of the
         ! removed leaves, so re-pointing it afterwards is safe.
         block_parent_of_last = self%block_parent((last_first - 1_id_k) / k)
         do i = 0_id_k, k - 1_id_k
            src = last_first + i
            dst = elder + i

            moved_elder = HTG_TREE_LEAF_SENTINEL
            if (src < self%elder_active_len) then
               moved_elder = self%parent_to_elder_child(src)
               self%parent_to_elder_child(src) = HTG_TREE_LEAF_SENTINEL
            end if
            ! dst < src, so dst is inside the active elder range whenever the
            ! moved vertex is refined; a leaf beyond the range needs no entry.
            if (dst < self%elder_active_len) self%parent_to_elder_child(dst) = moved_elder
            if (moved_elder /= HTG_TREE_LEAF_SENTINEL .and. moved_elder >= 0_id_k) then
               self%block_parent((moved_elder - 1_id_k) / k) = dst
            end if

            if (has_table) self%global_index_table(dst) = self%global_index_table(src)
            if (src == local_id) info%node_id = dst
         end do
         self%parent_to_elder_child(block_parent_of_last) = elder
         self%block_parent((elder - 1_id_k) / k) = block_parent_of_last
         info%moved_from = last_first
         info%moved_to = elder
      end if

      if (has_table) self%global_index_table(last_first : nv_old - 1_id_k) = HTG_INVALID_ID

      self%number_of_vertices = nv_old - k
      self%number_of_nodes = self%number_of_nodes - 1_id_k
      self%number_of_leaves = self%number_of_leaves - (k - 1_id_k)
      self%number_of_vertices_per_depth(child_depth) = &
         self%number_of_vertices_per_depth(child_depth) - k
      ! With no vertex left at child_depth there can be none deeper either
      if (self%number_of_vertices_per_depth(child_depth) == 0_id_k) then
         self%number_of_levels = child_depth
      end if

      ! The lazy elder tail must not reach past the remaining vertices
      self%elder_active_len = min(self%elder_active_len, self%number_of_vertices)

      call htg_set_ok(status)
   end subroutine tree_coarsen_node

   pure subroutine tree_subdivide_leaf(self, local_id, level, status)
      class(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: local_id
      integer(i32_k), intent(in) :: level
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: elder, old_verts, k_children, child_depth, target_levels

      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%subdivide_leaf: tree is not initialized', status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%subdivide_leaf: local_id out of bounds', status)
         return
      end if

      if (level < 0_i32_k .or. int(level, id_k) >= self%number_of_levels) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_tree_t%subdivide_leaf: level out of valid range', status)
         return
      end if

      if ((local_id == 0_id_k .and. level /= 0_i32_k) .or. &
          (local_id > 0_id_k .and. level == 0_i32_k)) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_tree_t%subdivide_leaf: level inconsistent with root/non-root vertex', status)
         return
      end if

      if (.not. is_leaf_unchecked(self, local_id)) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_tree_t%subdivide_leaf: vertex is already refined (not a leaf)', status)
         return
      end if

      old_verts = self%number_of_vertices
      elder = old_verts
      k_children = int(self%number_of_children, id_k)
      child_depth = int(level, id_k) + 1_id_k
      target_levels = max(self%number_of_levels, child_depth + 1_id_k)

      call ensure_depth_capacity(self, target_levels)
      if (target_levels > self%number_of_levels) then
         self%number_of_vertices_per_depth(self%number_of_levels : target_levels - 1_id_k) = 0_id_k
      end if

      if (local_id >= self%elder_active_len) then
         call ensure_elder_capacity(self, local_id + 1_id_k)
         self%parent_to_elder_child(self%elder_active_len : local_id - 1_id_k) = &
            HTG_TREE_LEAF_SENTINEL
         self%elder_active_len = local_id + 1_id_k
      end if
      self%parent_to_elder_child(local_id) = elder

      ! The new children form sibling block number_of_nodes (elder == 1 + k*number_of_nodes)
      call ensure_block_parent_capacity(self, self%number_of_nodes + 1_id_k)
      self%block_parent(self%number_of_nodes) = local_id

      self%number_of_nodes = self%number_of_nodes + 1_id_k
      self%number_of_leaves = self%number_of_leaves + (k_children - 1_id_k)
      self%number_of_vertices = old_verts + k_children
      self%number_of_vertices_per_depth(child_depth) = &
         self%number_of_vertices_per_depth(child_depth) + k_children
      self%number_of_levels = target_levels

      if (self%indexing_mode == HTG_INDEX_EXPLICIT) then
         call ensure_global_table_capacity(self, self%number_of_vertices, old_verts)
      end if

      call htg_set_ok(status)
   end subroutine tree_subdivide_leaf

   pure subroutine tree_subdivide_leaf_i64(self, local_id, level, status)
      class(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: local_id
      integer(id_k), intent(in) :: level
      type(htg_status_t), intent(out), optional :: status

      if (level < 0_id_k .or. level > int(huge(0_i32_k), id_k)) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_tree_t%subdivide_leaf: level out of 32-bit range', status)
         return
      end if

      call self%tree_subdivide_leaf(local_id, int(level, i32_k), status)
   end subroutine tree_subdivide_leaf_i64

   pure subroutine tree_set_global_index_start(self, start_id, status)
      class(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: start_id
      type(htg_status_t), intent(out), optional :: status

      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%set_global_index_start: tree is not initialized', status)
         return
      end if

      if (self%indexing_mode == HTG_INDEX_EXPLICIT) then
         call htg_fail(HTG_ERR_INDEXING, &
            'htg_tree_t%set_global_index_start: cannot mix implicit and explicit global indexing', &
            status)
         return
      end if

      if (start_id < 0_id_k) then
         call htg_fail(HTG_ERR_INDEXING, &
            'htg_tree_t%set_global_index_start: start_id must be non-negative', status)
         return
      end if

      self%indexing_mode = HTG_INDEX_IMPLICIT
      self%global_index_start = start_id
      call htg_set_ok(status)
   end subroutine tree_set_global_index_start

   pure function tree_get_global_index_start(self) result(start_id)
      class(htg_tree_t), intent(in) :: self
      integer(id_k) :: start_id
      start_id = self%global_index_start
   end function tree_get_global_index_start

   pure subroutine tree_set_global_index_from_local(self, local_id, global_id, status)
      class(htg_tree_t), intent(inout) :: self
      integer(id_k), intent(in) :: local_id
      integer(id_k), intent(in) :: global_id
      type(htg_status_t), intent(out), optional :: status

      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%set_global_index_from_local: tree is not initialized', status)
         return
      end if

      if (self%indexing_mode == HTG_INDEX_IMPLICIT) then
         call htg_fail(HTG_ERR_INDEXING, &
            'htg_tree_t%set_global_index_from_local: cannot mix explicit and implicit indexing', &
            status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%set_global_index_from_local: local_id out of bounds', status)
         return
      end if

      if (global_id < 0_id_k) then
         call htg_fail(HTG_ERR_INDEXING, &
            'htg_tree_t%set_global_index_from_local: global_id must be non-negative', status)
         return
      end if

      if (self%indexing_mode == HTG_INDEX_UNSET) then
         call ensure_global_table_capacity(self, self%number_of_vertices, 0_id_k)
         self%indexing_mode = HTG_INDEX_EXPLICIT
      else
         call ensure_global_table_capacity(self, self%number_of_vertices, self%number_of_vertices)
      end if

      self%global_index_table(local_id) = global_id
      call htg_set_ok(status)
   end subroutine tree_set_global_index_from_local

   function tree_get_global_index_from_local(self, local_id, status) result(global_id)
      class(htg_tree_t), intent(in) :: self
      integer(id_k), intent(in) :: local_id
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: global_id

      global_id = HTG_INVALID_ID
      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%get_global_index_from_local: tree is not initialized', status)
         return
      end if

      if (local_id < 0_id_k .or. local_id >= self%number_of_vertices) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%get_global_index_from_local: local_id out of bounds', status)
         return
      end if

      select case (self%indexing_mode)
      case (HTG_INDEX_IMPLICIT)
         global_id = self%global_index_start + local_id
         call htg_set_ok(status)
      case (HTG_INDEX_EXPLICIT)
         if (.not. allocated(self%global_index_table)) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_tree_t%get_global_index_from_local: explicit global table unallocated', status)
            return
         end if
         global_id = self%global_index_table(local_id)
         if (global_id < 0_id_k) then
            global_id = HTG_INVALID_ID
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_tree_t%get_global_index_from_local: explicit global id unset for vertex', &
               status)
            return
         end if
         call htg_set_ok(status)
      case default
         call htg_fail(HTG_ERR_INDEXING, &
            'htg_tree_t%get_global_index_from_local: global indexing mode is unset', status)
      end select
   end function tree_get_global_index_from_local

   pure function tree_get_indexing_mode(self) result(mode)
      class(htg_tree_t), intent(in) :: self
      integer(i32_k) :: mode
      mode = self%indexing_mode
   end function tree_get_indexing_mode

   function tree_get_max_global_index(self, status) result(max_gid)
      class(htg_tree_t), intent(in) :: self
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: max_gid

      max_gid = HTG_INVALID_ID
      if (.not. self%initialized .or. self%number_of_vertices <= 0_id_k) then
         if (present(status)) then
            call htg_fail(HTG_ERR_INVALID_STATE, &
               'htg_tree_t%get_max_global_index: tree is not initialized', status)
         end if
         return
      end if

      select case (self%indexing_mode)
      case (HTG_INDEX_IMPLICIT)
         max_gid = self%global_index_start + self%number_of_vertices - 1_id_k
         call htg_set_ok(status)
      case (HTG_INDEX_EXPLICIT)
         if (.not. allocated(self%global_index_table)) then
            if (present(status)) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_tree_t%get_max_global_index: explicit global table unallocated', status)
            end if
            return
         end if
         max_gid = maxval(self%global_index_table(0_id_k : self%number_of_vertices - 1_id_k))
         if (max_gid < 0_id_k) then
            max_gid = HTG_INVALID_ID
            if (present(status)) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_tree_t%get_max_global_index: all explicit global indices are unset', status)
            end if
            return
         end if
         call htg_set_ok(status)
      case default
         if (present(status)) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_tree_t%get_max_global_index: global indexing mode is unset', status)
         end if
      end select
   end function tree_get_max_global_index

   subroutine tree_build_from_breadth_first_order_descriptor( &
      self, descriptor, number_of_bits, start_index, status)
      class(htg_tree_t), intent(inout) :: self
      type(htg_bit_array_t), intent(in) :: descriptor
      integer(id_k), intent(in) :: number_of_bits
      integer(id_k), intent(in), optional :: start_index
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: bit_st
      integer(id_k) :: s_idx, desc_len, k_children, bit_idx, offset, depth_idx
      integer(id_k) :: num_levels, num_coarse, num_verts
      integer(id_k) :: cur_depth_size, next_depth_size, pos_at_depth, old_verts
      logical :: is_coarse

      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%build_from_breadth_first_order_descriptor: tree is not initialized', &
            status)
         return
      end if

      s_idx = 0_id_k
      if (present(start_index)) s_idx = start_index

      if (s_idx < 0_id_k .or. number_of_bits < 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_tree_t%build_from_breadth_first_order_descriptor: negative start or bit count', &
            status)
         return
      end if

      desc_len = descriptor%get_number_of_values()
      if (s_idx + number_of_bits > desc_len) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_tree_t%build_from_breadth_first_order_descriptor: slice exceeds descriptor size', &
            status)
         return
      end if

      old_verts = self%number_of_vertices
      k_children = int(self%number_of_children, id_k)

      if (number_of_bits == 0_id_k) then
         call ensure_depth_capacity(self, 1_id_k)
         self%number_of_vertices_per_depth(0_id_k) = 1_id_k

         call ensure_elder_capacity(self, 1_id_k)
         self%parent_to_elder_child(0_id_k) = HTG_TREE_LEAF_SENTINEL
         self%elder_active_len = 1_id_k

         self%number_of_levels = 1_id_k
         self%number_of_nodes = 0_id_k
         self%number_of_vertices = 1_id_k
         self%number_of_leaves = 1_id_k

         if (self%indexing_mode == HTG_INDEX_EXPLICIT) then
            call ensure_global_table_capacity(self, 1_id_k, min(old_verts, 1_id_k))
         end if

         call htg_set_ok(status)
         return
      end if

      ! Pass 1: Validate descriptor completeness without mutating self
      num_levels = 1_id_k
      num_coarse = 0_id_k
      num_verts = 1_id_k
      cur_depth_size = 1_id_k
      next_depth_size = 0_id_k
      pos_at_depth = 0_id_k

      do offset = 0_id_k, number_of_bits - 1_id_k
         if (cur_depth_size == 0_id_k) then
            call htg_fail(HTG_ERR_TOPOLOGY, &
               'htg_tree_t%build_from_breadth_first_order_descriptor: trailing orphan bits', &
               status)
            return
         end if

         bit_idx = s_idx + offset
         is_coarse = descriptor%get_value(bit_idx, bit_st)
         if (.not. bit_st%ok()) then
            call htg_fail(bit_st%code, trim(bit_st%message), status)
            return
         end if

         if (is_coarse) then
            num_verts = num_verts + k_children
            num_coarse = num_coarse + 1_id_k
            next_depth_size = next_depth_size + k_children
         end if

         pos_at_depth = pos_at_depth + 1_id_k
         if (pos_at_depth == cur_depth_size) then
            num_levels = num_levels + 1_id_k
            cur_depth_size = next_depth_size
            next_depth_size = 0_id_k
            pos_at_depth = 0_id_k
         end if
      end do

      if (pos_at_depth /= 0_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_tree_t%build_from_breadth_first_order_descriptor: truncated breadth-first level', &
            status)
         return
      end if

      if (cur_depth_size <= 0_id_k) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_tree_t%build_from_breadth_first_order_descriptor: terminal leaf depth encoded', &
            status)
         return
      end if

      ! Pass 2: Materialize topology and per-depth counts
      call ensure_depth_capacity(self, num_levels)
      self%number_of_vertices_per_depth(0_id_k : num_levels - 1_id_k) = 0_id_k
      self%number_of_vertices_per_depth(0_id_k) = 1_id_k

      call ensure_elder_capacity(self, number_of_bits)
      self%elder_active_len = number_of_bits

      ! Rebuilding replaces every sibling block; drop the old parents first so
      ! the capacity helper does not copy them.
      self%number_of_nodes = 0_id_k
      call ensure_block_parent_capacity(self, num_coarse)

      num_verts = 1_id_k
      cur_depth_size = 1_id_k
      next_depth_size = 0_id_k
      pos_at_depth = 0_id_k
      depth_idx = 0_id_k

      do offset = 0_id_k, number_of_bits - 1_id_k
         is_coarse = descriptor%get_value(s_idx + offset)
         if (is_coarse) then
            self%parent_to_elder_child(offset) = num_verts
            self%block_parent((num_verts - 1_id_k) / k_children) = offset
            num_verts = num_verts + k_children
            next_depth_size = next_depth_size + k_children
         else
            self%parent_to_elder_child(offset) = HTG_TREE_LEAF_SENTINEL
         end if

         pos_at_depth = pos_at_depth + 1_id_k
         if (pos_at_depth == cur_depth_size) then
            depth_idx = depth_idx + 1_id_k
            self%number_of_vertices_per_depth(depth_idx) = next_depth_size
            cur_depth_size = next_depth_size
            next_depth_size = 0_id_k
            pos_at_depth = 0_id_k
         end if
      end do

      self%number_of_levels = num_levels
      self%number_of_nodes = num_coarse
      self%number_of_vertices = num_verts
      self%number_of_leaves = num_verts - num_coarse

      if (self%indexing_mode == HTG_INDEX_EXPLICIT) then
         call ensure_global_table_capacity( &
            self, self%number_of_vertices, min(old_verts, self%number_of_vertices))
      end if

      call htg_set_ok(status)
   end subroutine tree_build_from_breadth_first_order_descriptor

   subroutine tree_compute_breadth_first_order_descriptor( &
      self, descriptor, number_of_vertices_per_depth, breadth_first_id_map, &
      input_mask, depth_limiter, status)
      class(htg_tree_t), intent(in) :: self
      type(htg_bit_array_t), intent(inout) :: descriptor
      integer(id_k), allocatable, intent(inout) :: number_of_vertices_per_depth(:)
      integer(id_k), allocatable, intent(inout) :: breadth_first_id_map(:)
      type(htg_bit_array_t), intent(in), optional :: input_mask
      integer(i32_k), intent(in), optional :: depth_limiter
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: step_st
      integer(id_k) :: max_levels, total_cap, max_allowed_depth
      integer(id_k) :: sp, vid, d, idg, slot, elder, eff_levels
      integer(id_k) :: total_visited, total_desc_bits
      integer(id_k) :: old_nv_len, old_nv_lb, old_map_len, old_map_lb, old_desc_bits
      integer(id_k) :: write_pos, idx
      integer(i32_k) :: ichild
      logical :: masked
      integer(id_k), allocatable :: depth_offsets(:)
      integer(id_k), allocatable :: depth_counts(:)
      integer(id_k), allocatable :: flat_ids(:)
      logical, allocatable :: flat_desc(:)
      integer(id_k), allocatable :: stack_vid(:)
      integer(id_k), allocatable :: stack_depth(:)
      integer(id_k), allocatable :: tmp_i64(:)

      if (.not. self%initialized) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_tree_t%compute_breadth_first_order_descriptor: tree is not initialized', status)
         return
      end if

      if (self%indexing_mode == HTG_INDEX_UNSET) then
         call htg_fail(HTG_ERR_INDEXING, &
            'htg_tree_t%compute_breadth_first_order_descriptor: global indexing is unset', status)
         return
      end if

      max_allowed_depth = huge(0_id_k)
      if (present(depth_limiter)) then
         if (depth_limiter < 0_i32_k) then
            call htg_fail(HTG_ERR_INVALID_ARG, &
               'htg_tree_t%compute_breadth_first_order_descriptor: depth_limiter must be >= 0', &
               status)
            return
         end if
         max_allowed_depth = int(depth_limiter, id_k)
      end if

      max_levels = self%number_of_levels
      total_cap = self%number_of_vertices

      allocate(depth_offsets(0_id_k : max_levels - 1_id_k))
      allocate(depth_counts(0_id_k : max_levels - 1_id_k))
      allocate(flat_ids(0_id_k : total_cap - 1_id_k))
      allocate(flat_desc(0_id_k : total_cap - 1_id_k))
      allocate(stack_vid(0_id_k : total_cap - 1_id_k))
      allocate(stack_depth(0_id_k : total_cap - 1_id_k))

      depth_offsets(0_id_k) = 0_id_k
      do d = 1_id_k, max_levels - 1_id_k
         depth_offsets(d) = depth_offsets(d - 1_id_k) + &
            self%number_of_vertices_per_depth(d - 1_id_k)
      end do
      depth_counts = 0_id_k

      sp = 0_id_k
      stack_vid(0_id_k) = 0_id_k
      stack_depth(0_id_k) = 0_id_k

      do while (sp >= 0_id_k)
         vid = stack_vid(sp)
         d = stack_depth(sp)
         sp = sp - 1_id_k

         idg = self%get_global_index_from_local(vid, step_st)
         if (.not. step_st%ok()) then
            call htg_fail(step_st%code, trim(step_st%message), status)
            return
         end if

         masked = .false.
         if (present(input_mask)) then
            masked = input_mask%get_value(idg, step_st)
            if (.not. step_st%ok()) then
               call htg_fail(HTG_ERR_MASK, &
                  'htg_tree_t%compute_breadth_first_order_descriptor: input_mask out of bounds', &
                  status)
               return
            end if
         end if

         slot = depth_offsets(d) + depth_counts(d)
         depth_counts(d) = depth_counts(d) + 1_id_k
         flat_ids(slot) = idg

         if (.not. is_leaf_unchecked(self, vid) .and. (.not. masked) .and. &
             d < max_allowed_depth) then
            flat_desc(slot) = .true.
            elder = self%parent_to_elder_child(vid)
            do ichild = self%number_of_children - 1_i32_k, 0_i32_k, -1_i32_k
               sp = sp + 1_id_k
               stack_vid(sp) = elder + int(ichild, id_k)
               stack_depth(sp) = d + 1_id_k
            end do
         else
            flat_desc(slot) = .false.
         end if
      end do

      eff_levels = max_levels
      do while (eff_levels > 0_id_k)
         if (depth_counts(eff_levels - 1_id_k) > 0_id_k) exit
         eff_levels = eff_levels - 1_id_k
      end do

      total_visited = sum(depth_counts(0_id_k : eff_levels - 1_id_k))
      if (eff_levels > 1_id_k) then
         total_desc_bits = sum(depth_counts(0_id_k : eff_levels - 2_id_k))
      else
         total_desc_bits = 0_id_k
      end if

      ! Append per-depth counts with explicit 0: lower bound
      if (allocated(number_of_vertices_per_depth)) then
         old_nv_len = size(number_of_vertices_per_depth, kind=id_k)
         old_nv_lb = lbound(number_of_vertices_per_depth, 1, kind=id_k)
      else
         old_nv_len = 0_id_k
         old_nv_lb = 0_id_k
      end if
      allocate(tmp_i64(0_id_k : old_nv_len + eff_levels - 1_id_k))
      if (old_nv_len > 0_id_k) then
         tmp_i64(0_id_k : old_nv_len - 1_id_k) = &
            number_of_vertices_per_depth(old_nv_lb : old_nv_lb + old_nv_len - 1_id_k)
      end if
      tmp_i64(old_nv_len : old_nv_len + eff_levels - 1_id_k) = &
         depth_counts(0_id_k : eff_levels - 1_id_k)
      call move_alloc(tmp_i64, number_of_vertices_per_depth)

      ! Append breadth-first global ID map with explicit 0: lower bound
      if (allocated(breadth_first_id_map)) then
         old_map_len = size(breadth_first_id_map, kind=id_k)
         old_map_lb = lbound(breadth_first_id_map, 1, kind=id_k)
      else
         old_map_len = 0_id_k
         old_map_lb = 0_id_k
      end if
      allocate(tmp_i64(0_id_k : old_map_len + total_visited - 1_id_k))
      if (old_map_len > 0_id_k) then
         tmp_i64(0_id_k : old_map_len - 1_id_k) = &
            breadth_first_id_map(old_map_lb : old_map_lb + old_map_len - 1_id_k)
      end if
      write_pos = old_map_len
      do d = 0_id_k, eff_levels - 1_id_k
         do idx = 0_id_k, depth_counts(d) - 1_id_k
            tmp_i64(write_pos) = flat_ids(depth_offsets(d) + idx)
            write_pos = write_pos + 1_id_k
         end do
      end do
      call move_alloc(tmp_i64, breadth_first_id_map)

      ! Append descriptor bits for depths 0 .. eff_levels - 2
      if (total_desc_bits > 0_id_k) then
         old_desc_bits = descriptor%get_number_of_values()
         call descriptor%resize(old_desc_bits + total_desc_bits, step_st)
         if (.not. step_st%ok()) then
            call htg_fail(step_st%code, trim(step_st%message), status)
            return
         end if
         write_pos = old_desc_bits
         do d = 0_id_k, eff_levels - 2_id_k
            do idx = 0_id_k, depth_counts(d) - 1_id_k
               call descriptor%set_value(write_pos, flat_desc(depth_offsets(d) + idx), step_st)
               if (.not. step_st%ok()) then
                  call htg_fail(step_st%code, trim(step_st%message), status)
                  return
               end if
               write_pos = write_pos + 1_id_k
            end do
         end do
      end if

      call htg_set_ok(status)
   end subroutine tree_compute_breadth_first_order_descriptor

   subroutine tree_get_parent_to_elder_child_view(self, ptr)
      class(htg_tree_t), target, intent(inout) :: self
      integer(id_k), pointer, intent(out) :: ptr(:)

      if (.not. self%initialized .or. self%elder_active_len <= 0_id_k .or. &
          .not. allocated(self%parent_to_elder_child)) then
         ptr => null()
      else
         ptr(0_id_k : self%elder_active_len - 1_id_k) => &
            self%parent_to_elder_child(0_id_k : self%elder_active_len - 1_id_k)
      end if
   end subroutine tree_get_parent_to_elder_child_view

   subroutine tree_get_vertices_per_depth_view(self, ptr)
      class(htg_tree_t), target, intent(inout) :: self
      integer(id_k), pointer, intent(out) :: ptr(:)

      if (.not. self%initialized .or. self%number_of_levels <= 0_id_k .or. &
          .not. allocated(self%number_of_vertices_per_depth)) then
         ptr => null()
      else
         ptr(0_id_k : self%number_of_levels - 1_id_k) => &
            self%number_of_vertices_per_depth(0_id_k : self%number_of_levels - 1_id_k)
      end if
   end subroutine tree_get_vertices_per_depth_view

   subroutine tree_get_global_index_table_view(self, ptr)
      class(htg_tree_t), target, intent(inout) :: self
      integer(id_k), pointer, intent(out) :: ptr(:)

      if (.not. self%initialized .or. self%indexing_mode /= HTG_INDEX_EXPLICIT .or. &
          self%number_of_vertices <= 0_id_k .or. .not. allocated(self%global_index_table)) then
         ptr => null()
      else
         ptr(0_id_k : self%number_of_vertices - 1_id_k) => &
            self%global_index_table(0_id_k : self%number_of_vertices - 1_id_k)
      end if
   end subroutine tree_get_global_index_table_view

end module fvtkhtg_tree_m
