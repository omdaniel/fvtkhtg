!> Tests for coarsening (htg_tree_t%coarsen_node, htg_cursor_t%coarsen_node)
!> and for the parent table behind it (get_parent_index, get_child_index,
!> get_level).
!>
!> VTK has no coarsening operation to compare with, so correctness is checked
!> structurally: the grid validates after every operation, every cell keeps
!> its own data through hundreds of random refine/coarsen steps, and a
!> coarsened tree is indistinguishable from one built directly in that shape.
program test_coarsen
   use fvtkhtg
   implicit none

   integer(i64_k) :: rng_state = 987654321_i64_k

   ! State of the random stress test (see test_random_stress)
   integer(i32_k), parameter :: ST_MAX_LEVEL = 4_i32_k
   type(htg_grid_t), target :: st_grid
   type(htg_field_array_t), pointer :: st_tag => null(), st_flag => null(), st_w => null()
   integer(i64_k), allocatable :: st_ext(:)        ! caller-side mirror of 'tag', by global id
   integer(id_k), allocatable :: st_free_ids(:)    ! explicit-id allocator: stack of freed ids
   integer :: st_nfree = 0
   integer(id_k) :: st_next_id = 0_id_k, st_total = 0_id_k
   integer(i32_k) :: st_path(ST_MAX_LEVEL)
   logical :: st_explicit = .false.
   character(len=8) :: st_mode = ''

   call test_coarsen_single_node()
   call test_relocation_numbering()
   call test_relocated_node_and_cursor_stack()
   call test_relocated_ancestor()
   call test_parent_table()
   call test_mask_follows_cells()
   call test_matches_direct_build()
   call test_errors()
   call test_random_stress(use_explicit_ids=.false.)
   call test_random_stress(use_explicit_ids=.true.)

   print '(a)', 'ALL test_coarsen TESTS PASSED'

contains

   subroutine check(cond, msg)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) error stop 'FAIL: ' // msg
   end subroutine check

   pure logical function same(a, b)
      real(r64_k), intent(in) :: a, b
      same = .not. (a < b .or. a > b)
   end function same

   function rand_int(n) result(r)
      integer, intent(in) :: n
      integer :: r
      rng_state = mod(rng_state * 1103515245_i64_k + 12345_i64_k, 2147483648_i64_k)
      r = int(mod(rng_state / 65536_i64_k, int(n, i64_k)))
   end function rand_int

   subroutine make_quad_grid(g, nx, ny)
      type(htg_grid_t), target, intent(inout) :: g
      integer(i32_k), intent(in) :: nx, ny
      type(htg_status_t) :: s
      real(r64_k), allocatable :: x(:), y(:)
      integer(i32_k) :: i

      allocate(x(0:nx), y(0:ny))
      x = [(real(i, r64_k), i = 0_i32_k, nx)]
      y = [(2.0_r64_k * real(i, r64_k), i = 0_i32_k, ny)]
      call g%initialize()
      call g%set_dimensions(nx + 1_i32_k, ny + 1_i32_k, 1_i32_k, s)
      call g%set_x_coordinates(x, s)
      call g%set_y_coordinates(y, s)
      call check(s%ok(), 'make_quad_grid')
   end subroutine make_quad_grid

   ! ---------------------------------------------------------------------------

   subroutine test_coarsen_single_node()
      type(htg_tree_t), target :: tree
      type(htg_coarsen_info_t) :: info
      type(htg_status_t) :: s
      integer(id_k), pointer :: elder(:)
      logical :: leaf

      call tree%initialize(2_i32_k, 2_i32_k, 0_id_k, s)
      call tree%subdivide_leaf(0_id_k, 0_i32_k, s)
      call check(tree%get_number_of_vertices() == 5_id_k, 'refined root has 5 vertices')

      call tree%coarsen_node(0_id_k, info, s)
      call check(s%ok(), 'coarsen the root: ' // trim(s%message))
      call check(info%number_of_children == 4_i32_k .and. info%node_id == 0_id_k .and. &
         info%removed_first == 1_id_k, 'coarsen info for the root')
      call check(info%moved_from == HTG_INVALID_ID .and. info%moved_to == HTG_INVALID_ID, &
         'removing the last block relocates nothing')

      ! Back to exactly the state of a fresh single-leaf tree
      call check(tree%get_number_of_vertices() == 1_id_k .and. tree%get_number_of_nodes() == 0_id_k &
         .and. tree%get_number_of_leaves() == 1_id_k .and. tree%get_number_of_levels() == 1_id_k, &
         'counters after coarsening the root')
      leaf = tree%is_leaf(0_id_k, s)
      call check(leaf, 'root is a leaf again')
      call tree%get_parent_to_elder_child_view(elder)
      call check(size(elder) == 1 .and. elder(0) == HTG_TREE_LEAF_SENTINEL, &
         'elder table is back to one leaf entry')
      call htg_validate_tree(tree, s)
      call check(s%ok(), 'coarsened tree validates: ' // trim(s%message))

      ! It can be refined again
      call tree%subdivide_leaf(0_id_k, 0_i32_k, s)
      call check(s%ok() .and. tree%get_number_of_vertices() == 5_id_k, 'refine again after coarsening')
      call htg_validate_tree(tree, s)
      call check(s%ok(), 'tree validates after refine-coarsen-refine')
   end subroutine test_coarsen_single_node

   !> Removing a block that is not the last one relocates the last block.
   subroutine test_relocation_numbering()
      type(htg_tree_t) :: tree
      type(htg_coarsen_info_t) :: info
      type(htg_status_t) :: s
      integer(id_k) :: gid, i
      logical :: leaf

      ! Root -> 1..4;  vertex 1 -> 5..8;  vertex 2 -> 9..12
      call tree%initialize(2_i32_k, 2_i32_k, 0_id_k, s)
      call tree%subdivide_leaf(0_id_k, 0_i32_k, s)
      call tree%subdivide_leaf(1_id_k, 1_i32_k, s)
      call tree%subdivide_leaf(2_id_k, 1_i32_k, s)
      do i = 0_id_k, 12_id_k
         call tree%set_global_index_from_local(i, 100_id_k + i, s)   ! explicit ids
      end do

      call tree%coarsen_node(1_id_k, info, s)
      call check(s%ok(), 'coarsen vertex 1: ' // trim(s%message))
      call check(info%node_id == 1_id_k .and. info%removed_first == 5_id_k .and. &
         info%moved_from == 9_id_k .and. info%moved_to == 5_id_k, 'relocation info')
      call check(tree%get_number_of_vertices() == 9_id_k .and. tree%get_number_of_nodes() == 2_id_k, &
         'counters after relocation')
      call check(tree%get_number_of_levels() == 3_id_k, 'two levels of children remain')
      leaf = tree%is_leaf(1_id_k, s)
      call check(leaf, 'vertex 1 is a leaf')
      call check(tree%get_elder_child_index(2_id_k) == 5_id_k, 'vertex 2 now owns block 5..8')
      call check(tree%get_parent_index(5_id_k) == 2_id_k, 'relocated block has parent 2 (first)')
      call check(tree%get_parent_index(8_id_k) == 2_id_k, 'relocated block has parent 2 (last)')

      ! Explicit global ids travel with their vertices; the removed ones are gone
      do i = 0_id_k, 3_id_k
         gid = tree%get_global_index_from_local(5_id_k + i, s)
         call check(gid == 109_id_k + i, 'relocated vertices keep their explicit global ids')
      end do
      gid = tree%get_global_index_from_local(1_id_k, s)
      call check(gid == 101_id_k, 'untouched vertex keeps its global id')
      call check(tree%get_max_global_index() == 112_id_k, 'max global id after coarsening')
      call htg_validate_tree(tree, s, require_global_index=.true.)
      call check(s%ok(), 'relocated tree validates: ' // trim(s%message))

      ! One level left once vertex 2 is coarsened too
      call tree%coarsen_node(2_id_k, info, s)
      call check(s%ok() .and. info%moved_from == HTG_INVALID_ID, 'coarsen the last block')
      call check(tree%get_number_of_levels() == 2_id_k .and. tree%get_number_of_vertices() == 5_id_k, &
         'level count drops with the deepest level')
      call htg_validate_tree(tree, s, require_global_index=.true.)
      call check(s%ok(), 'tree validates after second coarsening')
   end subroutine test_relocation_numbering

   !> The coarsened node can itself be in the relocated block, and a cursor's
   !> ancestors can be too; the cursor must stay on the same cell throughout.
   subroutine test_relocated_node_and_cursor_stack()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_tree_t), pointer :: tree
      type(htg_status_t) :: s
      real(r64_k) :: before(0:5), after(0:5)
      logical :: leaf

      call make_quad_grid(g, 1_i32_k, 1_i32_k)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%get_tree(tree)

      ! Build, by local id:  0 -> 1..4;  1 -> 5..8;  5 -> 9..12;  2 -> 13..16
      call c%subdivide_leaf(s)
      call c%to_child(0_i32_k, s)          ! vertex 1
      call c%subdivide_leaf(s)
      call c%to_child(0_i32_k, s)          ! vertex 5
      call c%subdivide_leaf(s)             ! 9..12
      call c%to_root(s)
      call c%to_child(1_i32_k, s)          ! vertex 2
      call c%subdivide_leaf(s)             ! 13..16

      ! Coarsen vertex 5: its block 9..12 is replaced by 13..16 (children of 2)
      call c%to_root(s)
      call c%to_child(0_i32_k, s)
      call c%to_child(0_i32_k, s)
      call check(c%get_vertex_id() == 5_id_k, 'cursor on vertex 5')
      call c%coarsen_node(s)
      call check(s%ok(), 'coarsen vertex 5: ' // trim(s%message))
      call check(tree%get_elder_child_index(2_id_k) == 9_id_k, 'children of 2 moved to 9..12')

      ! Refine vertex 9 (first child of 2) -> 13..16, then coarsen vertex 1:
      ! block 5..8 goes, block 13..16 (children of 9) moves into it
      call c%to_root(s)
      call c%to_child(1_i32_k, s)
      call c%to_child(0_i32_k, s)
      call check(c%get_vertex_id() == 9_id_k, 'cursor on vertex 9')
      call c%subdivide_leaf(s)
      call c%to_root(s)
      call c%to_child(0_i32_k, s)
      call c%coarsen_node(s)
      call check(s%ok(), 'coarsen vertex 1: ' // trim(s%message))
      call check(tree%get_elder_child_index(9_id_k) == 5_id_k, &
         'children of 9 now sit below their parent, at 5..8')
      call check(tree%get_number_of_vertices() == 13_id_k, '13 vertices remain')
      call htg_validate_grid(g, s, allow_unset_global_index=.true.)
      call check(s%ok(), 'grid with a child block below its parent validates: ' // trim(s%message))

      ! Now vertex 9 is in the LAST block (9..12) while its children are at
      ! 5..8. Coarsening it relocates vertex 9 itself to 5.
      call c%to_root(s)
      call c%to_child(1_i32_k, s)
      call c%to_child(0_i32_k, s)
      call check(c%get_vertex_id() == 9_id_k .and. c%get_level() == 2_i32_k, 'cursor on vertex 9 again')
      call c%get_bounds(before, s)
      call c%coarsen_node(s)
      call check(s%ok(), 'coarsen vertex 9: ' // trim(s%message))
      call check(c%get_vertex_id() == 5_id_k, 'the coarsened cell was renumbered 9 -> 5')
      call check(c%get_level() == 2_i32_k, 'cursor level is unchanged')
      leaf = c%is_leaf(s)
      call check(leaf, 'cursor cell is now a leaf')
      call c%get_bounds(after, s)
      call check(all_same(before, after), 'cursor is still on the same cell (same bounds)')
      call check(c%get_child_index() == 0_i32_k, 'cursor cell is still child 0 of its parent')

      ! The cursor still navigates correctly from where it is
      call c%to_parent(s)
      call check(c%get_vertex_id() == 2_id_k, 'parent of the renumbered cell is vertex 2')
      call check(tree%get_number_of_vertices() == 9_id_k, '9 vertices remain')
      call htg_validate_tree(tree, s)
      call check(s%ok(), 'final tree validates: ' // trim(s%message))
   end subroutine test_relocated_node_and_cursor_stack

   !> An ancestor of the coarsened cell is renumbered while the cell is not:
   !> every level of the cursor's stack has to follow.
   subroutine test_relocated_ancestor()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_tree_t), pointer :: tree
      type(htg_status_t) :: s
      real(r64_k) :: before(0:5), after(0:5), parent_before(0:5), parent_after(0:5)

      call make_quad_grid(g, 1_i32_k, 1_i32_k)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%get_tree(tree)

      ! 0 -> 1..4;  1 -> 5..8;  5 -> 9..12;  9 -> 13..16;  2 -> 17..20
      call c%subdivide_leaf(s)
      call c%to_child(0_i32_k, s)
      call c%subdivide_leaf(s)
      call c%to_child(0_i32_k, s)
      call c%subdivide_leaf(s)
      call c%to_child(0_i32_k, s)
      call check(c%get_vertex_id() == 9_id_k, 'ancestor case: on vertex 9')
      call c%subdivide_leaf(s)
      call c%to_root(s)
      call c%to_child(1_i32_k, s)
      call c%subdivide_leaf(s)

      ! Coarsen 9: children of 2 move from 17..20 to 13..16
      call c%to_root(s)
      call c%to_child(0_i32_k, s)
      call c%to_child(0_i32_k, s)
      call c%to_child(0_i32_k, s)
      call c%coarsen_node(s)
      call check(s%ok(), 'ancestor case: coarsen 1: ' // trim(s%message))
      call check(tree%get_elder_child_index(2_id_k) == 13_id_k, 'ancestor case: step 1')

      ! Refine 13 (child 0 of 2) -> 17..20, coarsen 5: children of 13 move to 9..12
      call c%to_root(s)
      call c%to_child(1_i32_k, s)
      call c%to_child(0_i32_k, s)
      call check(c%get_vertex_id() == 13_id_k, 'ancestor case: on vertex 13')
      call c%subdivide_leaf(s)
      call c%to_root(s)
      call c%to_child(0_i32_k, s)
      call c%to_child(0_i32_k, s)
      call c%coarsen_node(s)
      call check(s%ok(), 'ancestor case: coarsen 2: ' // trim(s%message))
      call check(tree%get_elder_child_index(13_id_k) == 9_id_k, 'ancestor case: step 2')

      ! Refine 9 (child 0 of 13) -> 17..20, coarsen 1: children of 9 move to 5..8
      call c%to_root(s)
      call c%to_child(1_i32_k, s)
      call c%to_child(0_i32_k, s)
      call c%to_child(0_i32_k, s)
      call check(c%get_vertex_id() == 9_id_k, 'ancestor case: on vertex 9 below 13')
      call c%subdivide_leaf(s)
      call c%to_root(s)
      call c%to_child(0_i32_k, s)
      call c%coarsen_node(s)
      call check(s%ok(), 'ancestor case: coarsen 3: ' // trim(s%message))
      call check(tree%get_elder_child_index(9_id_k) == 5_id_k, 'ancestor case: step 3')
      call check(tree%get_number_of_vertices() == 17_id_k, 'ancestor case: 17 vertices')

      ! Cursor path 0 -> 2 -> 13 -> 9. Vertex 13 sits in the last block
      ! (13..16) and vertex 9's children are at 5..8: coarsening 9 renumbers
      ! its PARENT 13 -> 5 and leaves 9 itself alone.
      call c%to_root(s)
      call c%to_child(1_i32_k, s)
      call c%to_child(0_i32_k, s)
      call c%get_bounds(parent_before, s)
      call c%to_child(0_i32_k, s)
      call check(c%get_vertex_id() == 9_id_k .and. c%get_level() == 3_i32_k, 'ancestor case: start')
      call c%get_bounds(before, s)
      call c%coarsen_node(s)
      call check(s%ok(), 'ancestor case: coarsen: ' // trim(s%message))
      call check(c%get_vertex_id() == 9_id_k, 'the coarsened cell keeps its id')
      call c%get_bounds(after, s)
      call check(all_same(before, after), 'cursor is on the same cell after its parent moved')
      call c%to_parent(s)
      call check(c%get_vertex_id() == 5_id_k, 'the parent on the stack was renumbered 13 -> 5')
      call c%get_bounds(parent_after, s)
      call check(all_same(parent_before, parent_after), 'the parent is the same cell')
      call c%to_parent(s)
      call check(c%get_vertex_id() == 2_id_k, 'grandparent is vertex 2')
      call check(tree%get_number_of_vertices() == 13_id_k, 'ancestor case: 13 vertices remain')
      call htg_validate_tree(tree, s)
      call check(s%ok(), 'ancestor case: tree validates: ' // trim(s%message))
   end subroutine test_relocated_ancestor

   pure logical function all_same(a, b)
      real(r64_k), intent(in) :: a(:), b(:)
      integer :: i
      all_same = (size(a) == size(b))
      if (.not. all_same) return
      do i = 1, size(a)
         if (.not. same(a(i), b(i))) all_same = .false.
      end do
   end function all_same

   !> get_parent_index / get_child_index / get_level agree with a cursor walk.
   subroutine test_parent_table()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_tree_t), pointer :: tree
      type(htg_status_t) :: s
      integer :: iter
      integer(id_k) :: parent
      integer(i32_k) :: lvl

      call make_quad_grid(g, 1_i32_k, 1_i32_k)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%get_tree(tree)
      parent = tree%get_parent_index(0_id_k, s)
      call check(parent == HTG_INVALID_ID .and. s%ok(), 'root has no parent')
      call check(tree%get_child_index(0_id_k) == -1_i32_k, 'root has no child index')
      call check(tree%get_level(0_id_k) == 0_i32_k, 'root is at level 0')

      do iter = 1, 40
         call random_refine_or_coarsen(c, 4_i32_k)
         call c%to_root(s)
         call walk_and_check(c, tree, HTG_INVALID_ID)
      end do

      parent = tree%get_parent_index(tree%get_number_of_vertices(), s)
      call check(s%code == HTG_ERR_OUT_OF_BOUNDS, 'get_parent_index rejects an id out of range')
      call check(parent == HTG_INVALID_ID, 'out-of-range parent is invalid')
      lvl = tree%get_level(-1_id_k, s)
      call check(lvl == -1_i32_k .and. s%code == HTG_ERR_OUT_OF_BOUNDS, &
         'get_level rejects a negative id')
   end subroutine test_parent_table

   recursive subroutine walk_and_check(c, tree, expected_parent)
      type(htg_cursor_t), intent(inout) :: c
      type(htg_tree_t), intent(in) :: tree
      integer(id_k), intent(in) :: expected_parent
      type(htg_status_t) :: s
      integer(id_k) :: vid
      integer(i32_k) :: ic
      logical :: leaf

      vid = c%get_vertex_id()
      call check(tree%get_parent_index(vid) == expected_parent, 'parent table matches the cursor')
      call check(tree%get_level(vid) == c%get_level(), 'get_level matches the cursor')
      call check(tree%get_child_index(vid) == c%get_child_index(), 'child index matches the cursor')
      leaf = c%is_leaf(s)
      if (leaf) return
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         call check(tree%get_child_index(c%get_vertex_id()) == ic, 'child index along the walk')
         call walk_and_check(c, tree, vid)
         call c%to_parent(s)
      end do
   end subroutine walk_and_check

   !> One random step on the tree under `c`: descend at random, then refine
   !> the leaf reached or coarsen the first terminal node met on the way.
   subroutine random_refine_or_coarsen(c, max_level)
      type(htg_cursor_t), intent(inout) :: c
      integer(i32_k), intent(in) :: max_level
      type(htg_status_t) :: s
      logical :: leaf, terminal, want_coarsen

      want_coarsen = rand_int(100) < 40
      call c%to_root(s)
      do
         leaf = c%is_leaf(s)
         if (leaf) then
            if (c%get_level() < max_level) call c%subdivide_leaf(s)
            return
         end if
         terminal = c%is_terminal_node(s)
         if (terminal .and. want_coarsen) then
            call c%coarsen_node(s)
            call check(s%ok(), 'random coarsen: ' // trim(s%message))
            return
         end if
         call c%to_child(int(rand_int(int(c%get_number_of_children())), i32_k), s)
      end do
   end subroutine random_refine_or_coarsen

   !> In implicit mode the mask bit of a renumbered cell moves with it and the
   !> bits of ids that fall out of use are cleared.
   subroutine test_mask_follows_cells()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_bit_array_t), pointer :: pm
      type(htg_status_t) :: s
      integer(id_k), allocatable :: freed(:), from(:), to(:)
      integer(id_k) :: i
      logical :: m

      ! Root -> 1..4;  vertex 1 -> 5..8;  vertex 2 -> 9..12; implicit ids = local ids
      call make_quad_grid(g, 1_i32_k, 1_i32_k)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(0_i32_k, s)
      call c%subdivide_leaf(s)
      call c%to_root(s)
      call c%to_child(1_i32_k, s)
      call c%subdivide_leaf(s)
      call g%initialize_local_index_node(s)
      call g%set_cell_mask(5_id_k, .true., s)     ! a child of vertex 1: removed, bit must go
      call g%set_cell_mask(10_id_k, .true., s)    ! a child of vertex 2: moves to id 6
      call g%set_cell_mask(12_id_k, .false., s)
      call g%get_pure_mask(pm, s)
      call check(g%is_pure_mask_valid(), 'PureMask cached before coarsening')

      call c%to_root(s)
      call c%to_child(0_i32_k, s)
      call c%coarsen_node(s, freed_global_ids=freed, moved_from_global_ids=from, moved_to_global_ids=to)
      call check(s%ok(), 'coarsen with a mask: ' // trim(s%message))
      call check(.not. g%is_pure_mask_valid(), 'coarsening invalidates the PureMask')
      call check(size(from) == 4 .and. size(to) == 4 .and. size(freed) == 4, 'reported moves and frees')
      do i = 1_id_k, 4_id_k
         call check(from(i) == 8_id_k + i .and. to(i) == 4_id_k + i, 'cells 9..12 moved to 5..8')
         call check(freed(i) == 8_id_k + i, 'ids 9..12 fell out of use')
      end do
      m = g%is_cell_masked(6_id_k, s)
      call check(m, 'the moved cell is still masked, now at id 6')
      m = g%is_cell_masked(5_id_k, s)
      call check(.not. m, 'id 5 now holds an unmasked cell; the removed cell''s bit is gone')
      m = g%is_cell_masked(7_id_k, s)
      call check(.not. m, 'unmasked moved cell stays unmasked')
      do i = 9_id_k, 12_id_k
         m = g%is_cell_masked(i, s)
         call check(.not. m, 'mask bits of unused ids are cleared')
      end do
      call g%validate(s)
      call check(s%ok(), 'masked grid validates after coarsening: ' // trim(s%message))
   end subroutine test_mask_follows_cells

   !> After any refine/coarsen history, the tree's breadth-first description
   !> and counters equal those of a tree built directly from that description.
   subroutine test_matches_direct_build()
      type(htg_grid_t), target :: g, g2
      type(htg_cursor_t) :: c
      type(htg_tree_t), pointer :: tree, tree2
      type(htg_breadth_first_t) :: bf, bf2
      type(htg_bit_array_t) :: desc, desc2
      integer(id_k), allocatable :: vpd(:), idmap(:), vpd2(:), idmap2(:)
      type(htg_status_t) :: s
      integer :: iter
      integer(id_k) :: i
      logical :: b1, b2

      call make_quad_grid(g, 1_i32_k, 1_i32_k)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%get_tree(tree)
      do iter = 1, 120
         call random_refine_or_coarsen(c, 4_i32_k)
      end do
      call g%initialize_local_index_node(s)
      call g%validate(s)
      call check(s%ok(), 'history grid validates: ' // trim(s%message))
      call tree%compute_breadth_first_order_descriptor(desc, vpd, idmap, status=s)

      call make_quad_grid(g2, 1_i32_k, 1_i32_k)
      call g2%get_tree(0_id_k, tree2, create=.true., status=s)
      call tree2%build_from_breadth_first_order_descriptor(desc, desc%get_number_of_values(), 0_id_k, s)
      call check(s%ok(), 'direct build from the descriptor: ' // trim(s%message))
      call g2%initialize_local_index_node(s)
      call g2%validate(s)
      call check(s%ok(), 'tree built from a descriptor validates: ' // trim(s%message))
      call new_non_oriented_cursor(g2, 0_id_k, c, create=.false., status=s)
      call walk_and_check(c, tree2, HTG_INVALID_ID)   ! parent table of a descriptor-built tree
      call tree2%compute_breadth_first_order_descriptor(desc2, vpd2, idmap2, status=s)

      call check(tree%get_number_of_vertices() == tree2%get_number_of_vertices() .and. &
                 tree%get_number_of_nodes() == tree2%get_number_of_nodes() .and. &
                 tree%get_number_of_leaves() == tree2%get_number_of_leaves() .and. &
                 tree%get_number_of_levels() == tree2%get_number_of_levels(), &
         'counters match a direct build')
      call check(size(vpd) == size(vpd2), 'depth count matches a direct build')
      call check(all(vpd == vpd2), 'vertices per depth match a direct build')
      call check(desc%get_number_of_values() == desc2%get_number_of_values(), 'descriptor length')
      do i = 0_id_k, desc%get_number_of_values() - 1_id_k
         b1 = desc%get_value(i)
         b2 = desc2%get_value(i)
         call check(b1 .eqv. b2, 'descriptor bits match a direct build')
      end do

      ! and the exports agree cell for cell
      call g%export_breadth_first(bf, s)
      call g2%export_breadth_first(bf2, s)
      call check(all(bf%refined .eqv. bf2%refined) .and. all(bf%xtree == bf2%xtree), &
         'breadth-first exports match')
      call check(bf2%identity_order, 'the direct build is in export order')
   end subroutine test_matches_direct_build

   subroutine test_errors()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c, unset
      type(htg_tree_t) :: raw
      type(htg_coarsen_info_t) :: info
      type(htg_status_t) :: s
      integer(id_k), allocatable :: freed(:)

      call unset%coarsen_node(s)
      call check(s%code == HTG_ERR_INVALID_STATE, 'coarsen with an unset cursor')
      call raw%coarsen_node(0_id_k, info, s)
      call check(s%code == HTG_ERR_INVALID_STATE, 'coarsen an uninitialised tree')

      call make_quad_grid(g, 1_i32_k, 1_i32_k)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%coarsen_node(s, freed_global_ids=freed)
      call check(s%code == HTG_ERR_TOPOLOGY, 'coarsening a leaf is an error')

      call c%subdivide_leaf(s)
      call c%to_child(3_i32_k, s)
      call c%subdivide_leaf(s)
      call c%to_root(s)
      call c%coarsen_node(s)
      call check(s%code == HTG_ERR_TOPOLOGY, 'coarsening needs all children to be leaves')
      call check(g%get_number_of_cells() == 9_id_k, 'a rejected coarsening changes nothing')

      ! Unset indexing: nothing to free or move, but the tree is coarsened
      call c%to_child(3_i32_k, s)
      call c%coarsen_node(s, freed_global_ids=freed)
      call check(s%ok() .and. size(freed) == 0, 'coarsen without global indices')
      call check(g%get_number_of_cells() == 5_id_k, 'cells after coarsening')
   end subroutine test_errors

   ! ---------------------------------------------------------------------------
   ! Random stress: every cell carries data derived from its own path, so any
   ! mix-up of cells, ids or tuples is caught. The stress state lives at
   ! program level (st_*), since internal procedures cannot nest.
   ! ---------------------------------------------------------------------------

   pure function path_hash(tree_id, path, n) result(h)
      integer(id_k), intent(in) :: tree_id
      integer(i32_k), intent(in) :: path(:)
      integer, intent(in) :: n
      integer(i64_k) :: h
      integer :: l
      h = tree_id + 1_i64_k
      do l = 1, n
         h = h * 7_i64_k + int(path(l), i64_k) + 1_i64_k
      end do
   end function path_hash

   function take_id() result(id)
      integer(id_k) :: id
      if (st_nfree > 0) then
         id = st_free_ids(st_nfree)          ! reuse a freed id first
         st_nfree = st_nfree - 1
      else
         id = st_next_id
         st_next_id = st_next_id + 1_id_k
      end if
   end function take_id

   !> Write the data that identifies the cell under the cursor.
   subroutine store(cur, tree_id, n)
      type(htg_cursor_t), intent(in) :: cur
      integer(id_k), intent(in) :: tree_id
      integer, intent(in) :: n
      integer(i64_k) :: h
      integer(id_k) :: id
      type(htg_status_t) :: s

      h = path_hash(tree_id, st_path, n)
      id = cur%get_global_node_index(s)
      call check(s%ok(), st_mode // ': cell has a global id')
      call st_tag%set_scalar_i64(id, h, s)
      call st_flag%set_tuple_bit(id, [mod(h, 2_i64_k) == 0_i64_k, mod(h, 3_i64_k) == 0_i64_k], s)
      call st_w%set_tuple_r64(id, [0.5_r64_k * real(h, r64_k), real(n, r64_k), -1.0_r64_k], s)
      st_ext(id) = h
   end subroutine store

   subroutine verify_all()
      type(htg_cursor_t) :: v
      logical, allocatable :: used(:)
      integer(id_k) :: tt, ncells
      type(htg_status_t) :: s

      allocate(used(0 : st_total - 1_id_k))
      used = .false.
      ncells = 0_id_k
      do tt = 0_id_k, st_grid%get_max_number_of_trees() - 1_id_k
         call new_non_oriented_cursor(st_grid, tt, v, create=.false., status=s)
         call verify_cell(v, tt, 0, used, ncells)
      end do
      call check(ncells == st_grid%get_number_of_cells(), st_mode // ': every cell was verified')
      call check(count(used, kind=id_k) == ncells, st_mode // ': global ids are unique')
      ! Ids not used by any cell hold no caller data
      call check(all((st_ext < 0_i64_k) .eqv. (.not. used)), &
         st_mode // ': mirror marks exactly the used ids')
      if (st_explicit) then
         call check(st_next_id - int(st_nfree, id_k) == ncells, st_mode // ': no explicit id leaked')
      end if
   end subroutine verify_all

   recursive subroutine verify_cell(v, tree_id, n, used, ncells)
      type(htg_cursor_t), intent(inout) :: v
      integer(id_k), intent(in) :: tree_id
      integer, intent(in) :: n
      logical, intent(inout) :: used(0:)
      integer(id_k), intent(inout) :: ncells
      integer(i64_k) :: h, got
      integer(id_k) :: id
      integer(i32_k) :: jc
      logical :: bits(0:1), is_leaf
      real(r64_k) :: wv(0:2)
      type(htg_status_t) :: s

      h = path_hash(tree_id, st_path, n)
      id = v%get_global_node_index(s)
      call check(s%ok(), st_mode // ': verified cell has a global id')
      call check(.not. used(id), st_mode // ': global id used twice')
      used(id) = .true.
      ncells = ncells + 1_id_k
      got = st_tag%get_scalar_i64(id)
      call check(got == h, st_mode // ': tag follows its cell')
      call check(st_ext(id) == h, st_mode // ': caller-side mirror follows its cell')
      call st_flag%get_tuple_bit(id, bits, s)
      call check((bits(0) .eqv. mod(h, 2_i64_k) == 0_i64_k) .and. &
                 (bits(1) .eqv. mod(h, 3_i64_k) == 0_i64_k), st_mode // ': bit tuple follows its cell')
      call st_w%get_tuple_r64(id, wv, s)
      call check(same(wv(0), 0.5_r64_k * real(h, r64_k)) .and. same(wv(1), real(n, r64_k)), &
         st_mode // ': real tuple follows its cell')
      is_leaf = v%is_leaf(s)
      if (is_leaf) return
      do jc = 0_i32_k, 3_i32_k
         call v%to_child(jc, s)
         st_path(n + 1) = jc
         call verify_cell(v, tree_id, n + 1, used, ncells)
         call v%to_parent(s)
      end do
   end subroutine verify_cell

   subroutine test_random_stress(use_explicit_ids)
      logical, intent(in) :: use_explicit_ids

      integer(id_k), parameter :: CAP = 400_id_k      ! > 1 + 4 + 16 + 64 + 256
      integer, parameter :: NSTEPS = 700

      type(htg_cursor_t) :: c
      type(htg_cell_data_t), pointer :: cd
      type(htg_status_t) :: s
      integer(id_k), allocatable :: freed(:), from(:), to(:), child_gid(:)
      integer(id_k) :: ntrees, t, n
      integer(i64_k) :: got
      integer :: step, ncoarsened, nrefined, depth, i
      integer(i32_k) :: ic, k
      real(r64_k) :: before(0:5), after(0:5)
      logical :: leaf, terminal, want_coarsen

      st_explicit = use_explicit_ids
      st_mode = merge('explicit', 'implicit', use_explicit_ids)
      call make_quad_grid(st_grid, 2_i32_k, 2_i32_k)
      ntrees = st_grid%get_max_number_of_trees()
      st_total = ntrees * CAP
      call st_grid%get_cell_data(cd)
      call cd%add_array('tag', HTG_TYPE_INT64, 1_i32_k, st_total, s)
      call cd%add_array('flag', HTG_TYPE_BIT, 2_i32_k, st_total, s)
      call cd%add_array('w', HTG_TYPE_REAL64, 3_i32_k, st_total, s)
      call cd%get_array('tag', st_tag, s)
      call cd%get_array('flag', st_flag, s)
      call cd%get_array('w', st_w, s)
      if (allocated(st_ext)) deallocate(st_ext)
      if (allocated(st_free_ids)) deallocate(st_free_ids)
      allocate(st_ext(0 : st_total - 1_id_k), st_free_ids(st_total), child_gid(4))
      st_ext = -1_i64_k
      st_nfree = 0
      st_next_id = 0_id_k

      do t = 0_id_k, ntrees - 1_id_k
         call new_non_oriented_cursor(st_grid, t, c, create=.true., status=s)
         if (use_explicit_ids) then
            call c%set_global_index_from_local(take_id(), s)
         else
            call c%set_global_index_start(t * CAP, s)     ! widely spaced implicit ranges
         end if
         call store(c, t, 0)
      end do

      ncoarsened = 0
      nrefined = 0
      do step = 1, NSTEPS
         t = int(rand_int(int(ntrees)), id_k)
         want_coarsen = rand_int(100) < 45
         call new_non_oriented_cursor(st_grid, t, c, create=.false., status=s)
         depth = 0
         do
            leaf = c%is_leaf(s)
            if (leaf) then
               if (c%get_level() < ST_MAX_LEVEL) then
                  call c%subdivide_leaf(s)
                  call check(s%ok(), st_mode // ': subdivide: ' // trim(s%message))
                  do ic = 0_i32_k, 3_i32_k
                     call c%to_child(ic, s)
                     if (use_explicit_ids) call c%set_global_index_from_local(take_id(), s)
                     st_path(depth + 1) = ic
                     call store(c, t, depth + 1)
                     call c%to_parent(s)
                  end do
                  nrefined = nrefined + 1
               end if
               exit
            end if
            terminal = c%is_terminal_node(s)
            if (terminal .and. want_coarsen) then
               k = c%get_number_of_children()
               do ic = 0_i32_k, k - 1_i32_k
                  call c%to_child(ic, s)
                  child_gid(ic + 1) = c%get_global_node_index(s)
                  call c%to_parent(s)
               end do
               call c%get_bounds(before, s)
               call c%coarsen_node(s, freed_global_ids=freed, moved_from_global_ids=from, &
                  moved_to_global_ids=to)
               call check(s%ok(), st_mode // ': coarsen: ' // trim(s%message))
               call c%get_bounds(after, s)
               call check(all_same(before, after), st_mode // ': cursor stays on the coarsened cell')
               leaf = c%is_leaf(s)
               call check(leaf .and. c%get_level() == int(depth, i32_k), &
                  st_mode // ': coarsened cell is a leaf')

               call check(size(freed) == 4, st_mode // ': four ids are freed')
               if (use_explicit_ids) then
                  call check(size(from) == 0 .and. size(to) == 0, st_mode // ': explicit ids never move')
                  call check(all(freed == child_gid), st_mode // ': the children''s ids are freed')
                  do i = 1, 4
                     st_nfree = st_nfree + 1
                     st_free_ids(st_nfree) = freed(i)
                  end do
               else
                  call check(size(from) == size(to), st_mode // ': moves come in pairs')
                  do n = 1_id_k, size(from, kind=id_k)
                     st_ext(to(n)) = st_ext(from(n))          ! mirror the reported moves
                  end do
               end if
               do i = 1, 4
                  st_ext(freed(i)) = -1_i64_k
                  got = st_tag%get_scalar_i64(freed(i))
                  call check(got == 0_i64_k, st_mode // ': freed tuples are zeroed')
               end do
               ncoarsened = ncoarsened + 1
               exit
            end if
            ic = int(rand_int(4), i32_k)
            call c%to_child(ic, s)
            depth = depth + 1
            st_path(depth) = ic
         end do

         call st_grid%validate(s)
         call check(s%ok(), st_mode // ': grid validates after every step: ' // trim(s%message))
         if (mod(step, 7) == 0 .or. step == NSTEPS) call verify_all()
      end do
      call check(ncoarsened > 50 .and. nrefined > 50, &
         st_mode // ': the stress run exercised both operations')
   end subroutine test_random_stress

end program test_coarsen
