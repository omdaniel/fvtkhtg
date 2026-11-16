!> Tests for cell geometry (cursor%get_origin / get_size / get_bounds /
!> get_point) and face-neighbour lookup (cursor%get_face_neighbor).
!>
!> Neighbours are checked against a brute-force search that knows nothing
!> about tree structure: among all cells no deeper than the query cell that
!> touch its face and cover it in the other directions, the deepest one is the
!> neighbour. Exact agreement with VTK's geometry cursor and supercursor is
!> covered by the oracle suite (test_oracle_geometry_neighbors).
program test_geometry_neighbors
   use fvtkhtg
   implicit none

   type :: cell_t
      integer(id_k) :: tree_id = -1_id_k
      integer(id_k) :: vertex_id = -1_id_k
      integer(i32_k) :: level = -1_i32_k
      logical :: leaf = .false.
      real(r64_k) :: bounds(0:5) = 0.0_r64_k
   end type cell_t

   integer(i64_k) :: rng_state = 20261116_i64_k

   call test_geometry_values()
   call test_geometry_planes()
   call test_geometry_ternary_3d()
   call test_geometry_descending_and_missing_tree()
   call test_child_position_and_axes()
   call test_neighbors_all_configurations()
   call test_neighbor_mask_and_refinement()
   call test_errors()

   print '(a)', 'ALL test_geometry_neighbors TESTS PASSED'

contains

   subroutine check(cond, msg)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) error stop 'FAIL: ' // msg
   end subroutine check

   !> Exact equality without -Wcompare-reals
   pure logical function same(a, b)
      real(r64_k), intent(in) :: a, b
      same = .not. (a < b .or. a > b)
   end function same

   pure logical function all_same(a, b)
      real(r64_k), intent(in) :: a(:), b(:)
      integer :: i
      all_same = (size(a) == size(b))
      if (.not. all_same) return
      do i = 1, size(a)
         if (.not. same(a(i), b(i))) all_same = .false.
      end do
   end function all_same

   pure logical function near(a, b, tol)
      real(r64_k), intent(in) :: a, b, tol
      near = abs(a - b) <= tol
   end function near

   !> Deterministic pseudo-random number in [0, 1)
   function rand01() result(r)
      real(r64_k) :: r
      rng_state = mod(rng_state * 1103515245_i64_k + 12345_i64_k, 2147483648_i64_k)
      r = real(rng_state, r64_k) / 2147483648.0_r64_k
   end function rand01

   function coords(n, x0, dx, stretch) result(c)
      integer(i32_k), intent(in) :: n
      real(r64_k), intent(in) :: x0, dx, stretch
      real(r64_k), allocatable :: c(:)
      integer(i32_k) :: i
      allocate(c(0 : n - 1_i32_k))
      do i = 0_i32_k, n - 1_i32_k
         c(i) = x0 + dx * real(i, r64_k) + stretch * real(i * i, r64_k)
      end do
   end function coords

   subroutine make_grid(g, dims, bf, transposed)
      type(htg_grid_t), target, intent(inout) :: g
      integer(i32_k), intent(in) :: dims(3), bf
      logical, intent(in) :: transposed
      type(htg_status_t) :: s

      call g%initialize()
      call g%set_branch_factor(bf, s)
      call g%set_transposed_root_indexing(transposed, s)
      call g%set_dimensions(dims(1), dims(2), dims(3), s)
      call check(s%ok(), 'make_grid dimensions')
      call g%set_x_coordinates(coords(dims(1), -1.0_r64_k, 1.0_r64_k, 0.125_r64_k), s)
      call g%set_y_coordinates(coords(dims(2), 2.0_r64_k, 0.5_r64_k, 0.0625_r64_k), s)
      call g%set_z_coordinates(coords(dims(3), 10.0_r64_k, 2.0_r64_k, 0.25_r64_k), s)
      call check(s%ok(), 'make_grid coordinates')
   end subroutine make_grid

   !> Refine below the cursor at random, to at most max_level.
   recursive subroutine refine_random(c, max_level, prob)
      type(htg_cursor_t), intent(inout) :: c
      integer(i32_k), intent(in) :: max_level
      real(r64_k), intent(in) :: prob
      type(htg_status_t) :: s
      integer(i32_k) :: ic

      if (c%get_level() >= max_level) return
      if (rand01() >= prob) return
      call c%subdivide_leaf(s)
      call check(s%ok(), 'refine_random subdivide')
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         call refine_random(c, max_level, prob)
         call c%to_parent(s)
      end do
   end subroutine refine_random

   !> Append every cell below the cursor to `cells`.
   recursive subroutine collect_cells(c, cells, n)
      type(htg_cursor_t), intent(inout) :: c
      type(cell_t), allocatable, intent(inout) :: cells(:)
      integer, intent(inout) :: n
      type(cell_t), allocatable :: grown(:)
      type(htg_status_t) :: s
      integer(i32_k) :: ic
      logical :: leaf

      if (n == size(cells)) then
         allocate(grown(2 * size(cells)))
         grown(1:n) = cells(1:n)
         call move_alloc(grown, cells)
      end if
      leaf = c%is_leaf(s)
      n = n + 1
      cells(n)%tree_id = c%get_tree_id()
      cells(n)%vertex_id = c%get_vertex_id()
      cells(n)%level = c%get_level()
      cells(n)%leaf = leaf
      call c%get_bounds(cells(n)%bounds, s)
      call check(s%ok(), 'collect_cells bounds')
      if (.not. leaf) then
         do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
            call c%to_child(ic, s)
            call collect_cells(c, cells, n)
            call c%to_parent(s)
         end do
      end if
   end subroutine collect_cells

   ! ---------------------------------------------------------------------------
   ! Geometry
   ! ---------------------------------------------------------------------------

   subroutine test_geometry_values()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      real(r64_k) :: origin(0:2), sz(0:2), bounds(0:5), point(0:2)

      ! 2 x 1 root cells in the xy plane: x in {0, 2, 6}, y in {1, 5}, z = 3
      call g%initialize()
      call g%set_dimensions(3_i32_k, 2_i32_k, 1_i32_k, s)
      call g%set_x_coordinates([0.0_r64_k, 2.0_r64_k, 6.0_r64_k], s)
      call g%set_y_coordinates([1.0_r64_k, 5.0_r64_k], s)
      call g%set_z_coordinates([3.0_r64_k], s)

      call g%get_level_zero_origin_and_size_from_index(1_id_k, origin, sz, s)
      call check(s%ok(), 'level-zero origin and size')
      call check(all_same(origin, [2.0_r64_k, 1.0_r64_k, 3.0_r64_k]), 'root 1 origin')
      call check(all_same(sz, [4.0_r64_k, 4.0_r64_k, 0.0_r64_k]), 'root 1 size (0 on squeezed z)')
      call g%get_level_zero_origin_from_index(0_id_k, origin, s)
      call check(all_same(origin, [0.0_r64_k, 1.0_r64_k, 3.0_r64_k]), 'root 0 origin')

      call new_non_oriented_cursor(g, 1_id_k, c, create=.true., status=s)
      call c%get_bounds(bounds, s)
      call check(s%ok(), 'root bounds')
      call check(all_same(bounds, [2.0_r64_k, 6.0_r64_k, 1.0_r64_k, 5.0_r64_k, 3.0_r64_k, 3.0_r64_k]), &
         'root 1 bounds')
      call c%get_point(point, s)
      call check(all_same(point, [4.0_r64_k, 3.0_r64_k, 3.0_r64_k]), 'root 1 centre')

      ! Child 3 (+x, +y) of the root, then child 1 (+x, -y) of that
      call c%subdivide_leaf(s)
      call c%to_child(3_i32_k, s)
      call c%get_origin(origin, s)
      call c%get_size(sz, s)
      call check(all_same(origin, [4.0_r64_k, 3.0_r64_k, 3.0_r64_k]), 'child 3 origin')
      call check(all_same(sz, [2.0_r64_k, 2.0_r64_k, 0.0_r64_k]), 'child 3 size')
      call c%subdivide_leaf(s)
      call c%to_child(1_i32_k, s)
      call c%get_bounds(bounds, s)
      call check(all_same(bounds, [5.0_r64_k, 6.0_r64_k, 3.0_r64_k, 4.0_r64_k, 3.0_r64_k, 3.0_r64_k]), &
         'grandchild bounds')
      call c%get_point(point, s)
      call check(all_same(point, [5.5_r64_k, 3.5_r64_k, 3.0_r64_k]), 'grandchild centre')
      call check(c%get_child_index() == 1_i32_k, 'cursor child index')

      ! Moving back up gives the parent's geometry again
      call c%to_parent(s)
      call c%get_origin(origin, s)
      call check(all_same(origin, [4.0_r64_k, 3.0_r64_k, 3.0_r64_k]), 'origin after to_parent')
   end subroutine test_geometry_values

   !> In every plane child bit 0 follows the lower in-plane axis and bit 1 the
   !> higher one, even in the xz plane where get_axes() reports (z, x).
   subroutine test_geometry_planes()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      real(r64_k) :: o1(0:2), o2(0:2)
      integer(i32_k) :: axes(0:2)

      ! xz plane
      call make_grid(g, [3_i32_k, 1_i32_k, 3_i32_k], 2_i32_k, .false.)
      axes = g%get_axes()
      call check(axes(0) == 2_i32_k .and. axes(1) == 0_i32_k, 'xz plane get_axes is (z, x)')
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(1_i32_k, s)
      call c%get_origin(o1, s)
      call c%to_parent(s)
      call c%to_child(2_i32_k, s)
      call c%get_origin(o2, s)
      call check(o1(0) > -1.0_r64_k .and. same(o1(2), 10.0_r64_k), 'xz plane: child 1 is offset in x')
      call check(same(o2(0), -1.0_r64_k) .and. o2(2) > 10.0_r64_k, 'xz plane: child 2 is offset in z')
      call check(same(o1(1), 2.0_r64_k) .and. same(o2(1), 2.0_r64_k), 'xz plane: y is untouched')

      ! yz plane
      call make_grid(g, [1_i32_k, 3_i32_k, 3_i32_k], 2_i32_k, .false.)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(1_i32_k, s)
      call c%get_origin(o1, s)
      call c%to_parent(s)
      call c%to_child(2_i32_k, s)
      call c%get_origin(o2, s)
      call check(o1(1) > 2.0_r64_k .and. same(o1(2), 10.0_r64_k), 'yz plane: child 1 is offset in y')
      call check(same(o2(1), 2.0_r64_k) .and. o2(2) > 10.0_r64_k, 'yz plane: child 2 is offset in z')

      ! 1D along y
      call make_grid(g, [1_i32_k, 4_i32_k, 1_i32_k], 3_i32_k, .false.)
      call new_non_oriented_cursor(g, 1_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(2_i32_k, s)
      call c%get_origin(o1, s)
      call c%get_size(o2, s)
      call check(same(o1(0), -1.0_r64_k) .and. same(o1(2), 10.0_r64_k), '1D y: x and z untouched')
      call check(o1(1) > 2.5625_r64_k, '1D y: child 2 is offset in y')
      call check(same(o2(0), 0.0_r64_k) .and. same(o2(2), 0.0_r64_k) .and. o2(1) > 0.0_r64_k, &
         '1D y: size is zero on the squeezed axes')
   end subroutine test_geometry_planes

   subroutine test_geometry_ternary_3d()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      real(r64_k) :: origin(0:2), sz(0:2), parent(0:5), bounds(0:5), lo(0:2), hi(0:2)
      integer(i32_k) :: ic, a

      call g%initialize()
      call g%set_branch_factor(3_i32_k, s)
      call g%set_dimensions(2_i32_k, 2_i32_k, 2_i32_k, s)
      call g%set_x_coordinates([0.0_r64_k, 9.0_r64_k], s)
      call g%set_y_coordinates([0.0_r64_k, 18.0_r64_k], s)
      call g%set_z_coordinates([0.0_r64_k, 27.0_r64_k], s)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)

      ! Child 23 = 2 + 3*1 + 9*2: third in x, second in y, third in z
      call c%to_child(23_i32_k, s)
      call c%get_origin(origin, s)
      call c%get_size(sz, s)
      call check(all_same(origin, [6.0_r64_k, 6.0_r64_k, 18.0_r64_k]), 'ternary child 23 origin')
      call check(all_same(sz, [3.0_r64_k, 6.0_r64_k, 9.0_r64_k]), 'ternary child 23 size')
      call check(g%get_child_position(23_i32_k, HTG_AXIS_X) == 2_i32_k .and. &
                 g%get_child_position(23_i32_k, HTG_AXIS_Y) == 1_i32_k .and. &
                 g%get_child_position(23_i32_k, HTG_AXIS_Z) == 2_i32_k, 'ternary child positions')
      call c%to_parent(s)

      ! The 27 children tile their parent
      call c%get_bounds(parent, s)
      lo = huge(1.0_r64_k)
      hi = -huge(1.0_r64_k)
      do ic = 0_i32_k, 26_i32_k
         call c%to_child(ic, s)
         call c%get_bounds(bounds, s)
         do a = 0_i32_k, 2_i32_k
            lo(a) = min(lo(a), bounds(2*a))
            hi(a) = max(hi(a), bounds(2*a + 1))
         end do
         call c%to_parent(s)
      end do
      do a = 0_i32_k, 2_i32_k
         call check(same(lo(a), parent(2*a)) .and. near(hi(a), parent(2*a + 1), 1.0e-12_r64_k), &
            'ternary children tile the parent')
      end do
   end subroutine test_geometry_ternary_3d

   subroutine test_geometry_descending_and_missing_tree()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      real(r64_k) :: bounds(0:5), sz(0:2)

      ! Descending x: sizes are negative, as in VTK
      call g%initialize()
      call g%set_dimensions(3_i32_k, 2_i32_k, 1_i32_k, s)
      call g%set_x_coordinates([4.0_r64_k, 2.0_r64_k, 0.0_r64_k], s)
      call g%set_y_coordinates([0.0_r64_k, 1.0_r64_k], s)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(1_i32_k, s)
      call c%get_bounds(bounds, s)
      call c%get_size(sz, s)
      call check(same(sz(0), -1.0_r64_k), 'descending x gives a negative size')
      call check(same(bounds(0), 3.0_r64_k) .and. same(bounds(1), 2.0_r64_k), &
         'descending x bounds run from origin to origin + size')

      ! A root cell without a tree still has geometry
      call new_non_oriented_cursor(g, 1_id_k, c, create=.false., status=s)
      call check(s%ok() .and. .not. c%has_tree(), 'cursor on a root cell without a tree')
      call c%get_bounds(bounds, s)
      call check(s%ok(), 'bounds of a root cell without a tree')
      call check(same(bounds(0), 2.0_r64_k) .and. same(bounds(1), 0.0_r64_k), &
         'root cell geometry without a tree')
   end subroutine test_geometry_descending_and_missing_tree

   subroutine test_child_position_and_axes()
      type(htg_grid_t), target :: g
      integer(i32_k) :: active(0:2), nactive

      call g%initialize()
      call g%get_active_axes(active, nactive)
      call check(nactive == 0_i32_k, 'no active axes before dimensions are set')

      call make_grid(g, [3_i32_k, 1_i32_k, 3_i32_k], 2_i32_k, .false.)
      call g%get_active_axes(active, nactive)
      call check(nactive == 2_i32_k .and. active(0) == HTG_AXIS_X .and. active(1) == HTG_AXIS_Z &
         .and. active(2) == -1_i32_k, 'xz plane active axes are (x, z)')
      call check(g%get_child_position(1_i32_k, HTG_AXIS_X) == 1_i32_k, 'xz child 1 is +x')
      call check(g%get_child_position(1_i32_k, HTG_AXIS_Z) == 0_i32_k, 'xz child 1 is -z')
      call check(g%get_child_position(2_i32_k, HTG_AXIS_Z) == 1_i32_k, 'xz child 2 is +z')
      call check(g%get_child_position(1_i32_k, HTG_AXIS_Y) == -1_i32_k, 'squeezed axis has no position')
      call check(g%get_child_position(4_i32_k, HTG_AXIS_X) == -1_i32_k, 'child index out of range')
   end subroutine test_child_position_and_axes

   ! ---------------------------------------------------------------------------
   ! Neighbours
   ! ---------------------------------------------------------------------------

   subroutine test_neighbors_all_configurations()
      integer(i32_k), parameter :: dims_list(3, 7) = reshape([ &
         5_i32_k, 1_i32_k, 1_i32_k,   1_i32_k, 4_i32_k, 1_i32_k,   1_i32_k, 1_i32_k, 4_i32_k, &
         4_i32_k, 3_i32_k, 1_i32_k,   4_i32_k, 1_i32_k, 3_i32_k,   1_i32_k, 4_i32_k, 3_i32_k, &
         3_i32_k, 3_i32_k, 2_i32_k], [3, 7])
      integer :: idims, itrans, isparse, ncases
      integer(i32_k) :: bf

      ncases = 0
      do idims = 1, 7
         do bf = 2_i32_k, 3_i32_k
            do itrans = 0, 1
               do isparse = 0, 1
                  call check_neighbors(dims_list(:, idims), bf, itrans == 1, isparse == 1)
                  ncases = ncases + 1
               end do
            end do
         end do
      end do
      call check(ncases == 56, 'all neighbour configurations ran')
   end subroutine test_neighbors_all_configurations

   subroutine check_neighbors(dims, bf, transposed, sparse)
      integer(i32_k), intent(in) :: dims(3), bf
      logical, intent(in) :: transposed, sparse

      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c, nb, back
      type(htg_status_t) :: s
      type(cell_t), allocatable :: cells(:)
      integer :: n, i, j, best, ndim
      integer(i32_k) :: active(0:2), nactive, p, t, a, dir, max_level, lvl
      integer(id_k) :: tid, ntrees
      real(r64_k) :: tol, face, far, domain_lo(0:2), domain_hi(0:2), prob
      real(r64_k), pointer :: xc(:), yc(:), zc(:)
      logical :: ok, at_edge
      character(len=64) :: tag

      write (tag, '(a,3(i0,1x),a,i0,a,l1,a,l1)') 'dims ', dims, 'bf ', bf, ' T', transposed, ' S', sparse
      call make_grid(g, dims, bf, transposed)
      call g%get_active_axes(active, nactive)
      ndim = int(nactive)
      max_level = merge(3_i32_k, 2_i32_k, bf == 2_i32_k)
      if (ndim == 3) max_level = max_level - 1_i32_k
      prob = merge(0.55_r64_k, 0.4_r64_k, bf == 2_i32_k)

      ntrees = g%get_max_number_of_trees()
      do tid = 0_id_k, ntrees - 1_id_k
         if (sparse .and. mod(tid, 3_id_k) == 1_id_k) cycle
         call new_non_oriented_cursor(g, tid, c, create=.true., status=s)
         call check(s%ok(), trim(tag) // ': create tree')
         call refine_random(c, max_level, prob)
      end do
      call g%initialize_local_index_node(s)
      call g%validate(s)
      call check(s%ok(), trim(tag) // ': grid validates: ' // trim(s%message))

      allocate(cells(64))
      n = 0
      do tid = 0_id_k, ntrees - 1_id_k
         if (.not. g%has_tree(tid)) cycle
         call new_non_oriented_cursor(g, tid, c, create=.false., status=s)
         call collect_cells(c, cells, n)
      end do
      call check(int(n, id_k) == g%get_number_of_cells(), trim(tag) // ': every cell was visited')

      call g%get_x_coordinates_view(xc)
      call g%get_y_coordinates_view(yc)
      call g%get_z_coordinates_view(zc)
      domain_lo = [xc(0), yc(0), zc(0)]
      domain_hi = [xc(ubound(xc, 1)), yc(ubound(yc, 1)), zc(ubound(zc, 1))]
      tol = 1.0e-9_r64_k

      do i = 1, n
         ! Put a cursor on cell i by replaying its path from stored ancestry:
         ! navigate by geometry instead, which keeps this test independent
         call locate(g, cells(i), c)
         do p = 0_i32_k, nactive - 1_i32_k
            a = active(p)
            do dir = -1_i32_k, 1_i32_k, 2_i32_k
               face = cells(i)%bounds(2*a + merge(1, 0, dir > 0))

               ! Brute force: deepest cell, no deeper than cell i, across the face
               best = 0
               do j = 1, n
                  if (cells(j)%level > cells(i)%level) cycle
                  far = cells(j)%bounds(2*a + merge(0, 1, dir > 0))
                  if (.not. near(far, face, tol)) cycle
                  ! and on the far side of the face, not the near side
                  if (dir > 0 .and. cells(j)%bounds(2*a + 1) <= face + tol) cycle
                  if (dir < 0 .and. cells(j)%bounds(2*a) >= face - tol) cycle
                  ok = .true.
                  do t = 0_i32_k, nactive - 1_i32_k
                     if (active(t) == a) cycle
                     if (cells(j)%bounds(2*active(t)) > cells(i)%bounds(2*active(t)) + tol) ok = .false.
                     if (cells(j)%bounds(2*active(t) + 1) < cells(i)%bounds(2*active(t) + 1) - tol) ok = .false.
                  end do
                  if (.not. ok) cycle
                  if (best == 0) then
                     best = j
                  else if (cells(j)%level > cells(best)%level) then
                     best = j
                  end if
               end do

               call c%get_face_neighbor(a, dir, nb, s)
               call check(s%ok(), trim(tag) // ': get_face_neighbor status: ' // trim(s%message))

               if (best == 0) then
                  call check(.not. nb%has_tree(), trim(tag) // ': no neighbour expected')
                  at_edge = near(face, merge(domain_hi(a), domain_lo(a), dir > 0), tol)
                  if (at_edge) then
                     call check(nb%get_tree_id() == HTG_INVALID_ID, &
                        trim(tag) // ': edge of the grid reports no root cell')
                  else
                     call check(nb%get_tree_id() >= 0_id_k, &
                        trim(tag) // ': missing tree reports its root cell')
                     call check(.not. g%has_tree(nb%get_tree_id()), &
                        trim(tag) // ': reported root cell really has no tree')
                  end if
               else
                  call check(nb%has_tree(), trim(tag) // ': neighbour expected')
                  call check(nb%get_tree_id() == cells(best)%tree_id .and. &
                             nb%get_vertex_id() == cells(best)%vertex_id, &
                     trim(tag) // ': neighbour is the brute-force cell')
                  lvl = nb%get_level()
                  call check(lvl == cells(best)%level .and. lvl <= cells(i)%level, &
                     trim(tag) // ': neighbour level')
                  if (lvl < cells(i)%level) then
                     call check(cells(best)%leaf, trim(tag) // ': a coarser neighbour is a leaf')
                  else
                     ! Same level: looking back finds the starting cell
                     call nb%get_face_neighbor(a, -dir, back, s)
                     call check(back%has_tree(), trim(tag) // ': reciprocal neighbour exists')
                     call check(back%get_tree_id() == cells(i)%tree_id .and. &
                                back%get_vertex_id() == cells(i)%vertex_id, &
                        trim(tag) // ': reciprocal neighbour is the starting cell')
                  end if
               end if
            end do
         end do
      end do
   end subroutine check_neighbors

   !> Position cursor `c` on a recorded cell by descending towards its centre.
   subroutine locate(g, cell, c)
      type(htg_grid_t), target, intent(inout) :: g
      type(cell_t), intent(in) :: cell
      type(htg_cursor_t), intent(inout) :: c
      type(htg_status_t) :: s
      real(r64_k) :: centre(0:2), b(0:5)
      integer(i32_k) :: ic, lvl, a, active(0:2), nactive, p
      logical :: inside

      call g%get_active_axes(active, nactive)
      do a = 0_i32_k, 2_i32_k
         centre(a) = 0.5_r64_k * (cell%bounds(2*a) + cell%bounds(2*a + 1))
      end do
      call new_non_oriented_cursor(g, cell%tree_id, c, create=.false., status=s)
      do lvl = 1_i32_k, cell%level
         do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
            call c%to_child(ic, s)
            call c%get_bounds(b, s)
            inside = .true.
            do p = 0_i32_k, nactive - 1_i32_k
               a = active(p)
               if (centre(a) < b(2*a) .or. centre(a) > b(2*a + 1)) inside = .false.
            end do
            if (inside) exit
            call c%to_parent(s)
         end do
      end do
      call check(c%get_vertex_id() == cell%vertex_id, 'locate found the recorded cell')
   end subroutine locate

   subroutine test_neighbor_mask_and_refinement()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c, nb
      type(htg_status_t) :: s
      integer(id_k) :: gid
      logical :: masked, leaf

      ! Two root cells side by side in x. Tree 0: root refined, child 1 (+x, -y)
      ! refined. Tree 1: root refined, child 0 (-x, -y) refined.
      call make_grid(g, [3_i32_k, 2_i32_k, 1_i32_k], 2_i32_k, .false.)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(1_i32_k, s)
      call c%subdivide_leaf(s)
      call new_non_oriented_cursor(g, 1_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(0_i32_k, s)
      call c%subdivide_leaf(s)
      call g%initialize_local_index_node(s)

      ! From tree 0 / child 1 / child 1 (level 2), +x crosses into tree 1 and
      ! lands on the level-2 cell child 0 / child 0
      call new_non_oriented_cursor(g, 0_id_k, c, create=.false., status=s)
      call c%to_child(1_i32_k, s)
      call c%to_child(1_i32_k, s)
      call c%get_face_neighbor(HTG_AXIS_X, 1_i32_k, nb, s)
      call check(s%ok() .and. nb%has_tree(), 'cross-tree neighbour exists')
      call check(nb%get_tree_id() == 1_id_k .and. nb%get_level() == 2_i32_k, &
         'cross-tree neighbour is at the same level in tree 1')

      ! +y from the same cell stays inside the tree: sibling child 3
      call c%get_face_neighbor(HTG_AXIS_Y, 1_i32_k, nb, s)
      call check(nb%get_tree_id() == 0_id_k .and. nb%get_level() == 2_i32_k, 'sibling neighbour')
      call check(nb%get_child_index() == 3_i32_k, 'sibling neighbour is child 3')

      ! Mask tree 1's child 0 (a refined cell): the search stops there
      call new_non_oriented_cursor(g, 1_id_k, nb, create=.false., status=s)
      call nb%to_child(0_i32_k, s)
      gid = nb%get_global_node_index(s)
      call g%set_cell_mask(gid, .true., s)
      call c%get_face_neighbor(HTG_AXIS_X, 1_i32_k, nb, s)
      call check(s%ok() .and. nb%has_tree(), 'neighbour behind a mask still exists')
      masked = nb%is_masked(s)
      call check(nb%get_level() == 1_i32_k .and. masked, &
         'the search stops at the masked cell and returns it')
      call g%clear_mask()

      ! The neighbour cursor is a normal cursor: refine through it (2:1 balance)
      call new_non_oriented_cursor(g, 1_id_k, c, create=.false., status=s)
      call c%to_child(0_i32_k, s)
      call c%to_child(1_i32_k, s)                      ! level 2, tree 1
      call c%get_face_neighbor(HTG_AXIS_X, 1_i32_k, nb, s)
      leaf = nb%is_leaf(s)
      call check(nb%get_level() == 1_i32_k .and. leaf, 'coarser leaf neighbour')
      call nb%subdivide_leaf(s)
      call check(s%ok(), 'refine through the neighbour cursor')
      call c%get_face_neighbor(HTG_AXIS_X, 1_i32_k, nb, s)
      call check(nb%get_level() == 2_i32_k, 'neighbour is at the same level after refining')
      call g%initialize_local_index_node(s)
      call g%validate(s)
      call check(s%ok(), 'grid validates after refining through a neighbour')
   end subroutine test_neighbor_mask_and_refinement

   subroutine test_errors()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c, nb, unset
      type(htg_status_t) :: s
      real(r64_k) :: bounds(0:5)

      call unset%get_bounds(bounds, s)
      call check(s%code == HTG_ERR_INVALID_STATE, 'geometry of an unset cursor')
      call unset%get_face_neighbor(HTG_AXIS_X, 1_i32_k, nb, s)
      call check(s%code == HTG_ERR_INVALID_STATE, 'neighbour of an unset cursor')

      ! Dimensions without coordinates on an active axis
      call g%initialize()
      call g%set_dimensions(3_i32_k, 3_i32_k, 1_i32_k, s)
      call g%set_x_coordinates([0.0_r64_k, 1.0_r64_k, 2.0_r64_k], s)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%get_bounds(bounds, s)
      call check(s%code == HTG_ERR_INVALID_STATE, 'geometry needs coordinates on every active axis')
      call g%set_y_coordinates([0.0_r64_k, 1.0_r64_k, 2.0_r64_k], s)
      call c%get_bounds(bounds, s)
      call check(s%ok(), 'geometry once coordinates are complete')

      call c%get_face_neighbor(HTG_AXIS_Z, 1_i32_k, nb, s)
      call check(s%code == HTG_ERR_INVALID_ARG, 'neighbour along a squeezed axis')
      call c%get_face_neighbor(HTG_AXIS_X, 0_i32_k, nb, s)
      call check(s%code == HTG_ERR_INVALID_ARG, 'neighbour direction must be -1 or +1')
      call c%get_face_neighbor(HTG_AXIS_X, 2_i32_k, nb, s)
      call check(s%code == HTG_ERR_INVALID_ARG, 'neighbour direction +2 is rejected')

      ! A cursor on a root cell without a tree has no path to step from
      call new_non_oriented_cursor(g, 1_id_k, c, create=.false., status=s)
      call c%get_face_neighbor(HTG_AXIS_X, 1_i32_k, nb, s)
      call check(s%code == HTG_ERR_INVALID_STATE, 'neighbour of a root cell without a tree')
   end subroutine test_errors

end program test_geometry_neighbors
