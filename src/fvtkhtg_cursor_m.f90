module fvtkhtg_cursor_m
   use fvtkhtg_kinds_m, only: i32_k, id_k, r64_k, htg_status_t, &
      HTG_ERR_INVALID_ARG, HTG_ERR_OUT_OF_BOUNDS, HTG_ERR_INVALID_STATE, &
      HTG_ERR_TOPOLOGY, HTG_ERR_INDEXING, HTG_ERR_MASK, HTG_INVALID_ID, &
      HTG_INDEX_IMPLICIT, HTG_INDEX_EXPLICIT, htg_set_ok, htg_fail
   use fvtkhtg_tree_m, only: htg_tree_t, htg_coarsen_info_t
   use fvtkhtg_grid_m, only: htg_grid_t
   implicit none
   private

   integer(i32_k), parameter :: INITIAL_ENTRY_CAPACITY = 8_i32_k

   ! Lightweight stack entry matching vtkHyperTreeGridEntry
   type, public :: htg_entry_t
      integer(id_k) :: vertex_id = 0_id_k
   end type htg_entry_t

   ! Build-time stack cursor matching vtkHyperTreeGridNonOrientedCursor.
   ! NOTE ON TARGET LIFETIME SAFETY:
   ! The caller MUST declare the htg_grid_t actual argument with the `target`
   ! (or `pointer`) attribute so `cursor%grid` and `cursor%tree` remain valid
   ! for the lifetime of the cursor.
   type, public :: htg_cursor_t
      type(htg_grid_t), pointer :: grid => null()
      type(htg_tree_t), pointer :: tree => null()
      integer(id_k) :: tree_id = HTG_INVALID_ID
      integer(i32_k) :: last_valid_entry = -1_i32_k
      type(htg_entry_t), allocatable :: entries(:)
   contains
      procedure :: initialize => cursor_initialize
      procedure :: reset => cursor_reset
      procedure :: has_tree => cursor_has_tree
      procedure :: get_tree => cursor_get_tree
      procedure :: get_tree_id => cursor_get_tree_id
      procedure :: get_vertex_id => cursor_get_vertex_id
      procedure :: get_global_node_index => cursor_get_global_node_index
      procedure :: get_dimension => cursor_get_dimension
      procedure :: get_number_of_children => cursor_get_number_of_children
      procedure :: get_level => cursor_get_level
      procedure :: is_leaf => cursor_is_leaf
      procedure :: is_terminal_node => cursor_is_terminal_node
      procedure :: is_root => cursor_is_root
      procedure :: to_root => cursor_to_root
      procedure :: to_parent => cursor_to_parent
      procedure, private :: cursor_to_child_i32
      procedure, private :: cursor_to_child_i64
      generic :: to_child => cursor_to_child_i32, cursor_to_child_i64
      procedure :: subdivide_leaf => cursor_subdivide_leaf
      procedure :: coarsen_node => cursor_coarsen_node
      procedure :: get_child_index => cursor_get_child_index

      ! Cell geometry (fvtkhtg_cursor_geometry_s)
      procedure :: get_origin => cursor_get_origin
      procedure :: get_size => cursor_get_size
      procedure :: get_bounds => cursor_get_bounds
      procedure :: get_point => cursor_get_point
      ! Face-neighbour lookup (fvtkhtg_cursor_geometry_s)
      procedure :: get_face_neighbor => cursor_get_face_neighbor
      procedure :: set_global_index_start => cursor_set_global_index_start
      procedure :: set_global_index_from_local => cursor_set_global_index_from_local
      procedure :: set_mask => cursor_set_mask
      procedure :: is_masked => cursor_is_masked
   end type htg_cursor_t

   public :: htg_new_non_oriented_cursor
   public :: new_non_oriented_cursor

   interface new_non_oriented_cursor
      module procedure htg_new_non_oriented_cursor
   end interface new_non_oriented_cursor

   ! Implemented in submodule fvtkhtg_cursor_geometry_s
   interface
      !> Low corner of the current cell, as vtkHyperTreeGridNonOrientedGeometryCursor::GetOrigin.
      module subroutine cursor_get_origin(self, origin, status)
         class(htg_cursor_t), intent(in) :: self
         real(r64_k), intent(out) :: origin(0:2)
         type(htg_status_t), intent(out), optional :: status
      end subroutine cursor_get_origin

      !> Edge lengths of the current cell (0 along a squeezed axis), as GetSize.
      module subroutine cursor_get_size(self, cell_size, status)
         class(htg_cursor_t), intent(in) :: self
         real(r64_k), intent(out) :: cell_size(0:2)
         type(htg_status_t), intent(out), optional :: status
      end subroutine cursor_get_size

      !> (xmin, xmax, ymin, ymax, zmin, zmax) of the current cell, as GetBounds.
      module subroutine cursor_get_bounds(self, bounds, status)
         class(htg_cursor_t), intent(in) :: self
         real(r64_k), intent(out) :: bounds(0:5)
         type(htg_status_t), intent(out), optional :: status
      end subroutine cursor_get_bounds

      !> Centre of the current cell, as GetPoint.
      module subroutine cursor_get_point(self, point, status)
         class(htg_cursor_t), intent(in) :: self
         real(r64_k), intent(out) :: point(0:2)
         type(htg_status_t), intent(out), optional :: status
      end subroutine cursor_get_point

      !> Position `neighbor` on the cell across the face of the current cell
      !> that looks along physical `axis` (HTG_AXIS_X/Y/Z) in `direction`
      !> (-1 or +1). See fvtkhtg_cursor_geometry_s for the rules.
      module subroutine cursor_get_face_neighbor(self, axis, direction, neighbor, status)
         class(htg_cursor_t), intent(in) :: self
         integer(i32_k), intent(in) :: axis
         integer(i32_k), intent(in) :: direction
         type(htg_cursor_t), intent(inout) :: neighbor
         type(htg_status_t), intent(out), optional :: status
      end subroutine cursor_get_face_neighbor
   end interface

contains

   subroutine cursor_reset(self)
      class(htg_cursor_t), intent(inout) :: self
      self%grid => null()
      self%tree => null()
      self%tree_id = HTG_INVALID_ID
      self%last_valid_entry = -1_i32_k
   end subroutine cursor_reset

   subroutine ensure_entry_capacity(self, min_idx)
      class(htg_cursor_t), intent(inout) :: self
      integer(i32_k), intent(in) :: min_idx

      integer(i32_k) :: old_cap, new_cap
      type(htg_entry_t), allocatable :: tmp(:)

      if (.not. allocated(self%entries)) then
         new_cap = max(INITIAL_ENTRY_CAPACITY, min_idx + 1_i32_k)
         allocate(self%entries(0_i32_k : new_cap - 1_i32_k))
         return
      end if

      old_cap = size(self%entries, kind=i32_k)
      if (min_idx < old_cap) return

      new_cap = max(2_i32_k * old_cap, min_idx + 1_i32_k)
      allocate(tmp(0_i32_k : new_cap - 1_i32_k))
      if (self%last_valid_entry >= 0_i32_k) then
         tmp(0_i32_k : self%last_valid_entry) = self%entries(0_i32_k : self%last_valid_entry)
      end if
      call move_alloc(tmp, self%entries)
   end subroutine ensure_entry_capacity

   ! Using `pointer, intent(in)` on `grid` allows passing either a `target`
   ! `type(htg_grid_t)` variable OR a `pointer` `type(htg_grid_t)` variable,
   ! while enforcing AT COMPILE TIME in Fortran 2018 that the caller declared
   ! `target` on non-pointer actual arguments so `self%grid` can never dangle!
   subroutine cursor_initialize(self, grid, tree_id, create, status)
      class(htg_cursor_t), intent(inout) :: self
      type(htg_grid_t), pointer, intent(in) :: grid
      integer(id_k), intent(in) :: tree_id
      logical, intent(in), optional :: create
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st

      call self%reset()
      if (.not. associated(grid)) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_cursor_t%initialize: grid pointer is not associated', status)
         return
      end if

      self%grid => grid
      self%tree_id = tree_id
      call self%grid%get_tree(tree_id, self%tree, create=create, status=st)
      if (.not. st%ok()) then
         call self%reset()
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      call ensure_entry_capacity(self, 0_i32_k)
      if (associated(self%tree)) then
         self%last_valid_entry = 0_i32_k
         self%entries(0)%vertex_id = 0_id_k
      else
         self%last_valid_entry = -1_i32_k
         self%entries(0)%vertex_id = 0_id_k
      end if

      call htg_set_ok(status)
   end subroutine cursor_initialize

   subroutine htg_new_non_oriented_cursor(grid, tree_id, cursor, create, status)
      type(htg_grid_t), pointer, intent(in) :: grid
      integer(id_k), intent(in) :: tree_id
      type(htg_cursor_t), intent(out) :: cursor
      logical, intent(in), optional :: create
      type(htg_status_t), intent(out), optional :: status

      call cursor%initialize(grid, tree_id, create=create, status=status)
   end subroutine htg_new_non_oriented_cursor

   pure logical function cursor_has_tree(self) result(flag)
      class(htg_cursor_t), intent(in) :: self
      flag = associated(self%tree)
   end function cursor_has_tree

   subroutine cursor_get_tree(self, tree_ptr)
      class(htg_cursor_t), intent(in) :: self
      type(htg_tree_t), pointer, intent(out) :: tree_ptr
      tree_ptr => self%tree
   end subroutine cursor_get_tree

   pure function cursor_get_tree_id(self) result(tid)
      class(htg_cursor_t), intent(in) :: self
      integer(id_k) :: tid
      tid = self%tree_id
   end function cursor_get_tree_id

   pure function cursor_get_vertex_id(self) result(vid)
      class(htg_cursor_t), intent(in) :: self
      integer(id_k) :: vid
      if (associated(self%tree) .and. self%last_valid_entry >= 0_i32_k .and. &
          allocated(self%entries)) then
         vid = self%entries(self%last_valid_entry)%vertex_id
      else
         vid = HTG_INVALID_ID
      end if
   end function cursor_get_vertex_id

   function cursor_get_global_node_index(self, status) result(gid)
      class(htg_cursor_t), intent(in) :: self
      type(htg_status_t), intent(out), optional :: status
      integer(id_k) :: gid

      gid = HTG_INVALID_ID
      if (.not. associated(self%tree) .or. self%last_valid_entry < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%get_global_node_index: cursor has no active tree', status)
         return
      end if

      gid = self%tree%get_global_index_from_local( &
         self%entries(self%last_valid_entry)%vertex_id, status)
   end function cursor_get_global_node_index

   pure function cursor_get_dimension(self) result(ndim)
      class(htg_cursor_t), intent(in) :: self
      integer(i32_k) :: ndim
      if (associated(self%grid)) then
         ndim = self%grid%get_dimension()
      else
         ndim = 0_i32_k
      end if
   end function cursor_get_dimension

   pure function cursor_get_number_of_children(self) result(nchild)
      class(htg_cursor_t), intent(in) :: self
      integer(i32_k) :: nchild
      if (associated(self%tree)) then
         nchild = self%tree%get_number_of_children()
      else if (associated(self%grid)) then
         nchild = self%grid%get_number_of_children()
      else
         nchild = 0_i32_k
      end if
   end function cursor_get_number_of_children

   pure function cursor_get_level(self) result(lvl)
      class(htg_cursor_t), intent(in) :: self
      integer(i32_k) :: lvl
      if (self%last_valid_entry >= 0_i32_k) then
         lvl = self%last_valid_entry
      else
         lvl = 0_i32_k
      end if
   end function cursor_get_level

   logical function cursor_is_leaf(self, status) result(flag)
      class(htg_cursor_t), intent(in) :: self
      type(htg_status_t), intent(out), optional :: status

      if (.not. associated(self%grid)) then
         flag = .false.
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%is_leaf: cursor is not initialized', status)
         return
      end if

      if (.not. associated(self%tree)) then
         flag = .true.
         call htg_set_ok(status)
         return
      end if

      if (self%last_valid_entry >= self%grid%get_depth_limiter()) then
         flag = .true.
         call htg_set_ok(status)
         return
      end if

      flag = self%tree%is_leaf(self%entries(self%last_valid_entry)%vertex_id, status)
   end function cursor_is_leaf

   logical function cursor_is_terminal_node(self, status) result(flag)
      class(htg_cursor_t), intent(in) :: self
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st
      logical :: leaf_flag

      flag = .false.
      if (.not. associated(self%tree) .or. self%last_valid_entry < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%is_terminal_node: cursor has no active tree', status)
         return
      end if

      leaf_flag = self%is_leaf(st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      if (leaf_flag) then
         flag = .false.
         call htg_set_ok(status)
         return
      end if

      flag = self%tree%is_terminal_node(self%entries(self%last_valid_entry)%vertex_id, status)
   end function cursor_is_terminal_node

   pure logical function cursor_is_root(self) result(flag)
      class(htg_cursor_t), intent(in) :: self
      if (.not. associated(self%tree)) then
         flag = (self%last_valid_entry <= 0_i32_k)
      else
         flag = (self%last_valid_entry == 0_i32_k .and. &
                 self%entries(0)%vertex_id == 0_id_k)
      end if
   end function cursor_is_root

   subroutine cursor_to_root(self, status)
      class(htg_cursor_t), intent(inout) :: self
      type(htg_status_t), intent(out), optional :: status

      if (.not. associated(self%grid)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%to_root: cursor is not initialized', status)
         return
      end if

      if (associated(self%tree)) then
         self%last_valid_entry = 0_i32_k
         self%entries(0)%vertex_id = 0_id_k
      end if
      call htg_set_ok(status)
   end subroutine cursor_to_root

   subroutine cursor_to_parent(self, status)
      class(htg_cursor_t), intent(inout) :: self
      type(htg_status_t), intent(out), optional :: status

      if (.not. associated(self%tree) .or. self%last_valid_entry < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%to_parent: cursor has no active tree', status)
         return
      end if

      if (self%is_root()) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%to_parent: cursor is already at root', status)
         return
      end if

      self%last_valid_entry = self%last_valid_entry - 1_i32_k
      call htg_set_ok(status)
   end subroutine cursor_to_parent

   subroutine cursor_to_child_i32(self, ichild, status)
      class(htg_cursor_t), intent(inout) :: self
      integer(i32_k), intent(in) :: ichild
      type(htg_status_t), intent(out), optional :: status

      integer(i32_k) :: nchild
      integer(id_k) :: cur_vid, elder
      logical :: leaf_flag, masked_flag
      type(htg_status_t) :: st

      if (.not. associated(self%grid) .or. .not. associated(self%tree) .or. &
          self%last_valid_entry < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%to_child: cursor has no active tree', status)
         return
      end if

      if (self%last_valid_entry > self%grid%get_depth_limiter()) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%to_child: cursor level exceeds depth_limiter', status)
         return
      end if

      nchild = self%tree%get_number_of_children()
      if (ichild < 0_i32_k .or. ichild >= nchild) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_cursor_t%to_child: ichild out of range [0, number_of_children-1]', status)
         return
      end if

      leaf_flag = self%is_leaf(st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      if (leaf_flag) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_cursor_t%to_child: cannot descend into child of a leaf node', status)
         return
      end if

      if (self%grid%has_mask()) then
         masked_flag = self%is_masked(st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
         if (masked_flag) then
            call htg_fail(HTG_ERR_MASK, &
               'htg_cursor_t%to_child: cannot descend into child of a masked node', status)
            return
         end if
      end if

      cur_vid = self%entries(self%last_valid_entry)%vertex_id
      elder = self%tree%get_elder_child_index(cur_vid, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      call ensure_entry_capacity(self, self%last_valid_entry + 1_i32_k)
      self%last_valid_entry = self%last_valid_entry + 1_i32_k
      self%entries(self%last_valid_entry)%vertex_id = elder + int(ichild, id_k)
      call htg_set_ok(status)
   end subroutine cursor_to_child_i32

   subroutine cursor_to_child_i64(self, ichild, status)
      class(htg_cursor_t), intent(inout) :: self
      integer(id_k), intent(in) :: ichild
      type(htg_status_t), intent(out), optional :: status

      if (ichild < 0_id_k .or. ichild > int(huge(0_i32_k), id_k)) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_cursor_t%to_child: ichild out of range', status)
         return
      end if
      call self%cursor_to_child_i32(int(ichild, i32_k), status)
   end subroutine cursor_to_child_i64

   subroutine cursor_subdivide_leaf(self, status)
      class(htg_cursor_t), intent(inout) :: self
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: cur_vid, max_gid
      logical :: leaf_flag, masked_flag
      type(htg_status_t) :: st

      if (.not. associated(self%grid) .or. .not. associated(self%tree) .or. &
          self%last_valid_entry < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%subdivide_leaf: cursor has no active tree', status)
         return
      end if

      if (self%last_valid_entry > self%grid%get_depth_limiter()) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%subdivide_leaf: cursor level exceeds depth_limiter', status)
         return
      end if

      if (self%grid%has_mask()) then
         masked_flag = self%is_masked(st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
         if (masked_flag) then
            call htg_fail(HTG_ERR_MASK, &
               'htg_cursor_t%subdivide_leaf: cannot subdivide a masked cell', status)
            return
         end if
      end if

      leaf_flag = self%is_leaf(st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      if (.not. leaf_flag) then
         call htg_set_ok(status)
         return
      end if

      cur_vid = self%entries(self%last_valid_entry)%vertex_id
      leaf_flag = self%tree%is_leaf(cur_vid, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      ! Matching VTK v9.7.2 Entry::SubdivideLeaf: coarse node is a silent no-op
      if (.not. leaf_flag) then
         call htg_set_ok(status)
         return
      end if

      call self%tree%subdivide_leaf(cur_vid, self%last_valid_entry, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      call self%grid%invalidate_pure_mask()

      if (self%grid%cell_data%get_number_of_arrays() > 0_i32_k) then
         max_gid = self%tree%get_max_global_index()
         if (max_gid >= 0_id_k) then
            if (max_gid == huge(0_id_k)) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_cursor_t%subdivide_leaf: max_global_index == huge(id_k) overflows', &
                  status)
               return
            end if
            call self%grid%cell_data%ensure_tuples(max_gid + 1_id_k, st)
            if (.not. st%ok()) then
               call htg_fail(st%code, trim(st%message), status)
               return
            end if
         end if
      end if

      call htg_set_ok(status)
   end subroutine cursor_subdivide_leaf

   !> Which child of its parent the current cell is (0 .. number_of_children - 1),
   !> or -1 at a root or when the cursor has no tree.
   function cursor_get_child_index(self) result(ichild)
      class(htg_cursor_t), intent(in) :: self
      integer(i32_k) :: ichild

      ichild = -1_i32_k
      if (.not. associated(self%tree) .or. self%last_valid_entry < 1_i32_k) return
      ichild = self%tree%get_child_index(self%entries(self%last_valid_entry)%vertex_id)
   end function cursor_get_child_index

   !> Inverse of subdivide_leaf: remove the children of the current cell,
   !> which must all be leaves, so that the current cell becomes a leaf.
   !> VTK has no such operation.
   !>
   !> The cursor stays on the same cell. Compute the cell's own value from its
   !> children BEFORE calling: the children's cell attributes are not kept.
   !>
   !> Numbering. Local ids stay contiguous, so the sibling block with the
   !> highest local ids in this tree is renumbered into the gap (unless the
   !> removed children were that block). No other cell changes.
   !>   - EXPLICIT global indexing: every surviving cell keeps its global id
   !>     and cell attributes. `freed_global_ids` lists the removed
   !>     children's ids, which are now free for the caller to reuse.
   !>   - IMPLICIT global indexing (global = start + local): the renumbered
   !>     cells get new global ids, and their CellData tuples and mask bits are
   !>     moved with them. `moved_from_global_ids(n)` -> `moved_to_global_ids(n)`
   !>     lists those moves, for arrays the caller keeps outside the grid;
   !>     `freed_global_ids` lists the ids left unused at the end of the tree's range.
   !> Freed ids have their CellData tuples zeroed and their mask bits cleared.
   !>
   !> Any other cursor on this tree may be left pointing at a renumbered or
   !> removed cell and must be re-initialised.
   subroutine cursor_coarsen_node(self, status, freed_global_ids, &
         moved_from_global_ids, moved_to_global_ids)
      class(htg_cursor_t), intent(inout) :: self
      type(htg_status_t), intent(out), optional :: status
      integer(id_k), allocatable, intent(out), optional :: freed_global_ids(:)
      integer(id_k), allocatable, intent(out), optional :: moved_from_global_ids(:)
      integer(id_k), allocatable, intent(out), optional :: moved_to_global_ids(:)

      type(htg_coarsen_info_t) :: info
      type(htg_status_t) :: st
      integer(id_k) :: cur_vid, elder, gid, start, k, i, nfreed, nmoved, vid
      integer(id_k), allocatable :: freed(:), moved_from(:), moved_to(:)
      integer(i32_k) :: mode, lvl
      logical :: leaf_flag, terminal

      if (.not. associated(self%grid) .or. .not. associated(self%tree) .or. &
          self%last_valid_entry < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%coarsen_node: cursor has no active tree', status)
         return
      end if

      if (self%last_valid_entry > self%grid%get_depth_limiter()) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%coarsen_node: cursor level exceeds depth_limiter', status)
         return
      end if

      cur_vid = self%entries(self%last_valid_entry)%vertex_id
      leaf_flag = self%tree%is_leaf(cur_vid, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if
      if (leaf_flag) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_cursor_t%coarsen_node: current cell is a leaf (nothing to coarsen)', status)
         return
      end if
      terminal = self%tree%is_terminal_node(cur_vid, st)
      if (.not. terminal) then
         call htg_fail(HTG_ERR_TOPOLOGY, &
            'htg_cursor_t%coarsen_node: children are not all leaves (coarsen bottom-up)', &
            status)
         return
      end if

      mode = self%tree%get_indexing_mode()
      k = int(self%tree%get_number_of_children(), id_k)
      elder = self%tree%get_elder_child_index(cur_vid, st)
      allocate(freed(k), moved_from(k), moved_to(k))
      nfreed = 0_id_k
      nmoved = 0_id_k

      ! EXPLICIT: the children's global ids are released (unset entries are skipped)
      if (mode == HTG_INDEX_EXPLICIT) then
         do i = 0_id_k, k - 1_id_k
            gid = self%tree%get_global_index_from_local(elder + i, st)
            if (st%ok()) then
               nfreed = nfreed + 1_id_k
               freed(nfreed) = gid
            end if
         end do
      end if

      call self%tree%coarsen_node(cur_vid, info, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      ! IMPLICIT: renumbered cells take their attributes to their new global
      ! ids, and the ids at the end of the tree's range fall out of use
      if (mode == HTG_INDEX_IMPLICIT) then
         start = self%tree%get_global_index_start()
         if (info%moved_from /= HTG_INVALID_ID) then
            do i = 0_id_k, k - 1_id_k
               nmoved = nmoved + 1_id_k
               moved_from(nmoved) = start + info%moved_from + i
               moved_to(nmoved) = start + info%moved_to + i
               call self%grid%copy_cell_attributes(moved_from(nmoved), moved_to(nmoved))
            end do
         end if
         do i = 0_id_k, k - 1_id_k
            nfreed = nfreed + 1_id_k
            freed(nfreed) = start + self%tree%get_number_of_vertices() + i
         end do
      end if

      do i = 1_id_k, nfreed
         call self%grid%clear_cell_attributes(freed(i))
      end do

      ! Keep this cursor on the same cell: the coarsened cell or any of its
      ! ancestors on the stack may be among the renumbered vertices
      if (info%moved_from /= HTG_INVALID_ID) then
         do lvl = 0_i32_k, self%last_valid_entry
            vid = self%entries(lvl)%vertex_id
            if (vid >= info%moved_from .and. vid < info%moved_from + k) then
               self%entries(lvl)%vertex_id = info%moved_to + (vid - info%moved_from)
            end if
         end do
      end if

      call self%grid%invalidate_pure_mask()

      if (present(freed_global_ids)) freed_global_ids = freed(1:nfreed)
      if (present(moved_from_global_ids)) moved_from_global_ids = moved_from(1:nmoved)
      if (present(moved_to_global_ids)) moved_to_global_ids = moved_to(1:nmoved)
      call htg_set_ok(status)
   end subroutine cursor_coarsen_node

   subroutine cursor_set_global_index_start(self, offset, status)
      class(htg_cursor_t), intent(inout) :: self
      integer(id_k), intent(in) :: offset
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: max_gid
      type(htg_status_t) :: st

      if (.not. associated(self%grid) .or. .not. associated(self%tree)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%set_global_index_start: cursor has no active tree', status)
         return
      end if

      call self%tree%set_global_index_start(offset, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      call self%grid%invalidate_pure_mask()

      if (self%grid%cell_data%get_number_of_arrays() > 0_i32_k) then
         max_gid = self%tree%get_max_global_index()
         if (max_gid >= 0_id_k) then
            if (max_gid == huge(0_id_k)) then
               call htg_fail(HTG_ERR_INDEXING, &
                  'htg_cursor_t%set_global_index_start: max_global_index overflows', &
                  status)
               return
            end if
            call self%grid%cell_data%ensure_tuples(max_gid + 1_id_k, st)
            if (.not. st%ok()) then
               call htg_fail(st%code, trim(st%message), status)
               return
            end if
         end if
      end if

      call htg_set_ok(status)
   end subroutine cursor_set_global_index_start

   subroutine cursor_set_global_index_from_local(self, global_id, status)
      class(htg_cursor_t), intent(inout) :: self
      integer(id_k), intent(in) :: global_id
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st

      if (.not. associated(self%grid) .or. .not. associated(self%tree) .or. &
          self%last_valid_entry < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%set_global_index_from_local: cursor has no active tree', status)
         return
      end if

      call self%tree%set_global_index_from_local( &
         self%entries(self%last_valid_entry)%vertex_id, global_id, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      call self%grid%invalidate_pure_mask()

      if (self%grid%cell_data%get_number_of_arrays() > 0_i32_k .and. global_id >= 0_id_k) then
         if (global_id == huge(0_id_k)) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_cursor_t%set_global_index_from_local: global_id == huge(id_k) overflows', &
               status)
            return
         end if
         call self%grid%cell_data%ensure_tuples(global_id + 1_id_k, st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
      end if

      call htg_set_ok(status)
   end subroutine cursor_set_global_index_from_local

   subroutine cursor_set_mask(self, state, status)
      class(htg_cursor_t), intent(inout) :: self
      logical, intent(in) :: state
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: gid
      type(htg_status_t) :: st

      if (.not. associated(self%grid) .or. .not. associated(self%tree) .or. &
          self%last_valid_entry < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%set_mask: cursor has no active tree', status)
         return
      end if

      gid = self%get_global_node_index(st)
      if (.not. st%ok()) then
         call htg_fail(st%code, trim(st%message), status)
         return
      end if

      call self%grid%set_cell_mask(gid, state, status)
   end subroutine cursor_set_mask

   logical function cursor_is_masked(self, status) result(masked)
      class(htg_cursor_t), intent(in) :: self
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: gid
      type(htg_status_t) :: st

      masked = .false.
      if (.not. associated(self%grid)) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%is_masked: cursor is not initialized', status)
         return
      end if

      if (associated(self%tree) .and. self%grid%has_mask()) then
         gid = self%get_global_node_index(st)
         if (.not. st%ok()) then
            call htg_fail(st%code, trim(st%message), status)
            return
         end if
         masked = self%grid%is_cell_masked(gid, status)
         return
      end if

      call htg_set_ok(status)
   end function cursor_is_masked

end module fvtkhtg_cursor_m
