!> Contract tests for the fvtkhtg -> serializer seam (grid%export_breadth_first
!> and htg_breadth_first_t%gather).
!>
!> The export must satisfy the input contract of the fvtkhdf HyperTreeGrid
!> writer, `vtkhdf_htg_file%write_mesh_topology(tree_ids, xtree, refined [, mask])`:
!>   - size(xtree) == ntree + 1, xtree(1) == 1, xtree(ntree+1) == ncell + 1;
!>   - 1 <= tree_ids <= number of root cells;
!>   - refined describes complete breadth-first trees;
!>   - mask, if present, has one entry per cell and no refined cell is masked.
!> and what the writer then stores (DepthPerTree, NumberOfCellsPerTreeDepth,
!> Descriptors, cell order) must equal what VTK's own vtkHDFWriter stores.
!> VTK's writer uses vtkHyperTree::ComputeBreadthFirstOrderDescriptor, whose
!> fvtkhtg counterpart is verified against VTK by the oracle suite, so it is
!> the reference here.
program test_export
   use fvtkhtg
   implicit none

   call test_contract_2d()
   call test_contract_3d_ternary()
   call test_empty_and_unrefined()
   call test_explicit_indexing_order()
   call test_gather()
   call test_errors()

   print '(a)', 'ALL test_export TESTS PASSED'

contains

   !> Exact equality of copied values, without -Wcompare-reals
   pure logical function same(a, b)
      real(r64_k), intent(in) :: a, b
      same = .not. (a < b .or. a > b)
   end function same

   subroutine check(cond, msg)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) error stop 'FAIL: ' // msg
   end subroutine check

   ! ---------------------------------------------------------------------------
   ! Fixtures
   ! ---------------------------------------------------------------------------

   !> 2D quadtree forest, 3 of 6 root cells populated. Tree 0 is refined out of
   !> breadth-first order (child 2 before child 0), so storage order differs
   !> from export order; tree 2 is a single leaf; tree 5 has three levels.
   subroutine build_forest_2d(g)
      type(htg_grid_t), target, intent(inout) :: g
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s

      call g%initialize()
      call g%set_branch_factor(2_i32_k, s)
      call g%set_dimensions(4_i32_k, 3_i32_k, 1_i32_k, s)
      call g%set_x_coordinates([0.0_r64_k, 1.0_r64_k, 2.0_r64_k, 3.0_r64_k], s)
      call g%set_y_coordinates([0.0_r64_k, 1.0_r64_k, 2.0_r64_k], s)

      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(2_i32_k, s)
      call c%subdivide_leaf(s)
      call c%to_root(s)
      call c%to_child(0_i32_k, s)
      call c%subdivide_leaf(s)

      call new_non_oriented_cursor(g, 2_id_k, c, create=.true., status=s)

      call new_non_oriented_cursor(g, 5_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(1_i32_k, s)
      call c%subdivide_leaf(s)
      call c%to_child(0_i32_k, s)
      call c%subdivide_leaf(s)

      call g%initialize_local_index_node(s)
      call g%validate(s)
      call check(s%ok(), 'fixture 2D validates: ' // trim(s%message))
   end subroutine build_forest_2d

   ! ---------------------------------------------------------------------------
   ! Independent derivation, from the writer's input alone, of what a VTKHDF
   ! writer stores: per-tree depth count, cells per depth, descriptor bits
   ! (refinement flags of every level but each tree's deepest). Returns
   ! valid = .false. if the flags do not describe complete trees.
   ! ---------------------------------------------------------------------------
   subroutine derive_layout(xtree, refined, nchild, depth_per_tree, cells_per_depth, bits, valid)
      integer(id_k), intent(in) :: xtree(:)
      logical, intent(in) :: refined(:)
      integer(id_k), intent(in) :: nchild
      integer(id_k), allocatable, intent(out) :: depth_per_tree(:), cells_per_depth(:)
      logical, allocatable, intent(out) :: bits(:)
      logical, intent(out) :: valid

      integer(id_k) :: t, first, level_size, n_refined, ntree

      ntree = size(xtree, kind=id_k) - 1_id_k
      allocate(depth_per_tree(ntree), cells_per_depth(0), bits(0))
      depth_per_tree = 0_id_k
      valid = .false.
      do t = 1_id_k, ntree
         first = xtree(t)
         level_size = 1_id_k
         do
            if (first + level_size - 1_id_k > xtree(t + 1_id_k) - 1_id_k) return
            depth_per_tree(t) = depth_per_tree(t) + 1_id_k
            cells_per_depth = [cells_per_depth, level_size]
            n_refined = count(refined(first : first + level_size - 1_id_k), kind=id_k)
            if (n_refined > 0_id_k) bits = [bits, refined(first : first + level_size - 1_id_k)]
            first = first + level_size
            if (n_refined == 0_id_k) exit
            level_size = n_refined * nchild
         end do
         if (first /= xtree(t + 1_id_k)) return
      end do
      valid = .true.
   end subroutine derive_layout

   !> Full contract check of one export against the writer's rules and VTK.
   subroutine check_contract(g, label)
      type(htg_grid_t), target, intent(inout) :: g
      character(len=*), intent(in) :: label

      type(htg_breadth_first_t) :: bf, bf2
      type(htg_status_t) :: s
      type(htg_bit_array_t) :: desc, bits_in
      type(htg_bit_array_t), pointer :: mask_ptr
      type(htg_tree_t), pointer :: tree
      type(htg_grid_t), target :: g2
      integer(id_k), allocatable :: depth_per_tree(:), cells_per_depth(:), vpd(:), idmap(:)
      logical, allocatable :: bits(:)
      logical :: valid, ref_bit, masked
      integer(id_k) :: ntree, ncell, pos, tid, i, prev_levels, d0, b0, nb

      call g%export_breadth_first(bf, s)
      call check(s%ok(), label // ': export ok: ' // trim(s%message))
      ntree = bf%number_of_trees()
      ncell = bf%number_of_cells()

      ! 1. Writer input contract
      call check(ntree == g%get_number_of_non_empty_trees(), label // ': one entry per tree')
      call check(size(bf%xtree, kind=id_k) == ntree + 1_id_k, label // ': size(xtree) == ntree+1')
      call check(bf%xtree(1) == 1_id_k, label // ': xtree(1) == 1')
      call check(bf%xtree(ntree + 1_id_k) == ncell + 1_id_k, label // ': xtree(ntree+1) == ncell+1')
      call check(size(bf%global_ids, kind=id_k) == ncell, label // ': one global id per cell')
      if (ntree > 0_id_k) then
         call check(minval(bf%tree_ids) >= 1_id_k, label // ': tree ids are 1-based')
         call check(maxval(bf%tree_ids) <= g%get_max_number_of_trees(), label // ': tree ids <= nroot')
      end if
      if (allocated(bf%mask)) then
         call check(size(bf%mask, kind=id_k) == ncell, label // ': one mask entry per cell')
         call check(.not. any(bf%mask .and. bf%refined), label // ': no refined cell is masked')
      end if
      call derive_layout(bf%xtree, bf%refined, int(g%get_number_of_children(), id_k), &
         depth_per_tree, cells_per_depth, bits, valid)
      call check(valid, label // ': refined describes complete breadth-first trees')

      ! 2. Stored layout equals VTK's
      call g%get_mask(mask_ptr)
      prev_levels = 0_id_k
      do pos = 0_id_k, ntree - 1_id_k
         call g%get_tree_by_pos(pos, tid, tree, s)
         call check(bf%tree_ids(pos + 1_id_k) == tid + 1_id_k, label // ': tree_ids == tree_id + 1')
         if (associated(mask_ptr)) then
            call tree%compute_breadth_first_order_descriptor(desc, vpd, idmap, &
               input_mask=mask_ptr, depth_limiter=g%get_depth_limiter(), status=s)
         else
            call tree%compute_breadth_first_order_descriptor(desc, vpd, idmap, &
               depth_limiter=g%get_depth_limiter(), status=s)
         end if
         call check(s%ok(), label // ': reference descriptor: ' // trim(s%message))
         call check(depth_per_tree(pos + 1_id_k) == size(vpd, kind=id_k) - prev_levels, &
            label // ': DepthPerTree matches VTK')
         prev_levels = size(vpd, kind=id_k)
      end do
      if (ntree > 0_id_k) then
         call check(size(cells_per_depth) == size(vpd), label // ': NumberOfDepths matches VTK')
         call check(all(cells_per_depth == vpd), label // ': NumberOfCellsPerTreeDepth matches VTK')
         call check(size(bf%global_ids) == size(idmap), label // ': NumberOfCells matches VTK')
         call check(all(bf%global_ids == idmap), label // ': cell order matches VTK id map')
         call check(size(bits, kind=id_k) == desc%get_number_of_values(), &
            label // ': DescriptorsSize matches VTK')
         do i = 1_id_k, size(bits, kind=id_k)
            ref_bit = desc%get_value(i - 1_id_k)
            call check(bits(i) .eqv. ref_bit, label // ': descriptor bits match VTK')
         end do
      end if
      if (allocated(bf%mask)) then
         do i = 1_id_k, ncell
            masked = g%is_cell_masked(bf%global_ids(i), s)
            call check(bf%mask(i) .eqv. masked, label // ': mask follows export order')
         end do
      end if

      ! 3. Round trip: rebuild the trees the way VTK's reader does
      call g2%initialize()
      call g2%set_branch_factor(g%get_branch_factor(), s)
      call g2%set_dimensions(g%get_dimensions(), s)
      call bits_in%initialize(size(bits, kind=id_k), .false., s)
      do i = 1_id_k, size(bits, kind=id_k)
         call bits_in%set_value(i - 1_id_k, bits(i), s)
      end do
      d0 = 0_id_k
      b0 = 0_id_k
      do pos = 1_id_k, ntree
         nb = sum(cells_per_depth(d0 + 1_id_k : d0 + depth_per_tree(pos) - 1_id_k))
         call g2%get_tree(bf%tree_ids(pos) - 1_id_k, tree, create=.true., status=s)
         call check(s%ok(), label // ': round trip creates tree')
         call tree%build_from_breadth_first_order_descriptor(bits_in, nb, b0, s)
         call check(s%ok(), label // ': round trip builds tree: ' // trim(s%message))
         d0 = d0 + depth_per_tree(pos)
         b0 = b0 + nb
      end do
      call g2%initialize_local_index_node(s)
      call g2%export_breadth_first(bf2, s)
      call check(s%ok(), label // ': round trip export')
      call check(size(bf2%refined) == size(bf%refined), label // ': round trip cell count')
      call check(all(bf2%refined .eqv. bf%refined), label // ': round trip refinement')
      call check(all(bf2%tree_ids == bf%tree_ids), label // ': round trip tree ids')
      call check(all(bf2%xtree == bf%xtree), label // ': round trip offsets')
      call check(bf2%identity_order, label // ': a rebuilt grid is in export order')
   end subroutine check_contract

   ! ---------------------------------------------------------------------------
   ! Tests
   ! ---------------------------------------------------------------------------

   subroutine test_contract_2d()
      type(htg_grid_t), target :: g
      type(htg_breadth_first_t) :: bf
      type(htg_status_t) :: s

      call build_forest_2d(g)
      call check_contract(g, '2D unmasked')

      ! Exact expected export for the fixture (tree 0 has 13 cells, tree 2 one,
      ! tree 5 thirteen); child 2 of tree 0 was refined first, so its
      ! grandchildren have smaller global ids than child 0's.
      call g%export_breadth_first(bf, s)
      call check(all(bf%tree_ids == [1_id_k, 3_id_k, 6_id_k]), '2D tree_ids')
      call check(all(bf%xtree == [1_id_k, 14_id_k, 15_id_k, 28_id_k]), '2D xtree')
      call check(all(bf%global_ids(1:13) == &
         [0_id_k, 1_id_k, 2_id_k, 3_id_k, 4_id_k, 9_id_k, 10_id_k, 11_id_k, 12_id_k, &
          5_id_k, 6_id_k, 7_id_k, 8_id_k]), '2D tree 0 cell order')
      call check(.not. bf%identity_order, '2D out-of-order refinement is not identity order')
      call check(.not. allocated(bf%mask), '2D unmasked grid exports no mask')

      ! Masked coarse cell (gid 16, tree 5) and masked leaf (gid 12, tree 0)
      call g%set_cell_mask(16_id_k, .true., s)
      call g%set_cell_mask(12_id_k, .true., s)
      call g%set_cell_mask(26_id_k, .false., s)
      call check_contract(g, '2D masked')
      call g%export_breadth_first(bf, s)
      call check(bf%number_of_cells() == 19_id_k, '2D masked coarse cell drops its 8 descendants')

      call g%set_depth_limiter(1_i32_k, s)
      call check_contract(g, '2D masked, depth_limiter=1')
      call g%set_depth_limiter(0_i32_k, s)
      call check_contract(g, '2D masked, depth_limiter=0')
      call g%export_breadth_first(bf, s)
      call check(bf%number_of_cells() == 3_id_k .and. .not. any(bf%refined), &
         '2D depth_limiter=0 exports root leaves only')
   end subroutine test_contract_2d

   subroutine test_contract_3d_ternary()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s

      call g%initialize()
      call g%set_branch_factor(3_i32_k, s)
      call g%set_dimensions(2_i32_k, 2_i32_k, 2_i32_k, s)
      call g%set_x_coordinates([0.0_r64_k, 1.0_r64_k], s)
      call g%set_y_coordinates([0.0_r64_k, 1.0_r64_k], s)
      call g%set_z_coordinates([0.0_r64_k, 1.0_r64_k], s)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%to_child(26_i32_k, s)
      call c%subdivide_leaf(s)
      call c%to_root(s)
      call c%to_child(13_i32_k, s)
      call c%subdivide_leaf(s)
      call c%to_child(5_i32_k, s)
      call c%subdivide_leaf(s)
      call g%initialize_local_index_node(s)
      call check_contract(g, '3D ternary')
   end subroutine test_contract_3d_ternary

   subroutine test_empty_and_unrefined()
      type(htg_grid_t), target :: g
      type(htg_breadth_first_t) :: bf
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s

      call g%initialize()
      call g%set_dimensions(3_i32_k, 3_i32_k, 1_i32_k, s)
      call g%set_x_coordinates([0.0_r64_k, 1.0_r64_k, 2.0_r64_k], s)
      call g%set_y_coordinates([0.0_r64_k, 1.0_r64_k, 2.0_r64_k], s)

      ! A rank that owns no trees still takes part in a collective write
      call g%export_breadth_first(bf, s)
      call check(s%ok(), 'empty grid exports')
      call check(bf%number_of_trees() == 0_id_k .and. bf%number_of_cells() == 0_id_k, &
         'empty grid has no trees or cells')
      call check(size(bf%xtree) == 1 .and. bf%xtree(1) == 1_id_k, 'empty grid xtree == [1]')
      call check(bf%identity_order, 'empty export is trivially in identity order')

      call g%export_breadth_first(bf, s, force_mask=.true.)
      call check(allocated(bf%mask), 'force_mask allocates a mask on an empty grid')
      call check(size(bf%mask) == 0, 'empty forced mask has no entries')

      call new_non_oriented_cursor(g, 3_id_k, c, create=.true., status=s)
      call g%initialize_local_index_node(s)
      call check_contract(g, 'single unrefined tree')
      call g%export_breadth_first(bf, s, force_mask=.true.)
      call check(all(bf%tree_ids == [4_id_k]) .and. all(bf%xtree == [1_id_k, 2_id_k]), &
         'single unrefined tree layout')
      call check(allocated(bf%mask), 'force_mask allocates a mask on an unmasked grid')
      call check(.not. any(bf%mask), 'forced mask is all false')
   end subroutine test_empty_and_unrefined

   subroutine test_explicit_indexing_order()
      type(htg_grid_t), target :: g
      type(htg_breadth_first_t) :: bf
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      integer(i32_k) :: ch

      call g%initialize()
      call g%set_dimensions(2_i32_k, 2_i32_k, 1_i32_k, s)
      call g%set_x_coordinates([0.0_r64_k, 1.0_r64_k], s)
      call g%set_y_coordinates([0.0_r64_k, 1.0_r64_k], s)
      call new_non_oriented_cursor(g, 0_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)
      call c%set_global_index_from_local(4_id_k, s)
      do ch = 0_i32_k, 3_i32_k
         call c%to_root(s)
         call c%to_child(ch, s)
         call c%set_global_index_from_local(int(3_i32_k - ch, id_k), s)
      end do
      call check_contract(g, 'explicit reversed global ids')
      call g%export_breadth_first(bf, s)
      call check(all(bf%global_ids == [4_id_k, 3_id_k, 2_id_k, 1_id_k, 0_id_k]), &
         'explicit ids are exported in breadth-first order')
      call check(.not. bf%identity_order, 'explicit reversed ids are not identity order')
   end subroutine test_explicit_indexing_order

   subroutine test_gather()
      type(htg_grid_t), target :: g
      type(htg_breadth_first_t) :: bf
      type(htg_status_t) :: s
      type(htg_cell_data_t), pointer :: cd
      type(htg_field_array_t), pointer :: f64, f32, fi64, fi32, fbit, fvec, fbits3
      real(r64_k), allocatable :: v64(:), vvec(:, :)
      real(r32_k), allocatable :: v32(:)
      integer(i64_k), allocatable :: vi64(:)
      integer(i32_k), allocatable :: vi32(:), vi32_2(:, :)
      integer(byte_k), allocatable :: vbit(:), vbit3(:, :)
      integer(id_k) :: gid, n, ncells
      integer(i32_k) :: comp
      real(r64_k) :: expected_vec(3)

      call build_forest_2d(g)
      ncells = g%get_number_of_cells()
      call g%get_cell_data(cd)
      call cd%add_array('r64', HTG_TYPE_REAL64, 1_i32_k, ncells, s)
      call cd%add_array('r32', HTG_TYPE_REAL32, 1_i32_k, ncells, s)
      call cd%add_array('i64', HTG_TYPE_INT64, 1_i32_k, ncells, s)
      call cd%add_array('i32', HTG_TYPE_INT32, 1_i32_k, ncells, s)
      call cd%add_array('bit', HTG_TYPE_BIT, 1_i32_k, ncells, s)
      call cd%add_array('vec', HTG_TYPE_REAL64, 3_i32_k, ncells, s)
      call cd%add_array('bits3', HTG_TYPE_BIT, 3_i32_k, ncells, s)
      call cd%get_array('r64', f64, s)
      call cd%get_array('r32', f32, s)
      call cd%get_array('i64', fi64, s)
      call cd%get_array('i32', fi32, s)
      call cd%get_array('bit', fbit, s)
      call cd%get_array('vec', fvec, s)
      call cd%get_array('bits3', fbits3, s)
      do gid = 0_id_k, ncells - 1_id_k
         call f64%set_scalar_r64(gid, 1.5_r64_k * real(gid, r64_k), s)
         call f32%set_scalar_r32(gid, -real(gid, r32_k), s)
         call fi64%set_scalar_i64(gid, 1000_i64_k + gid, s)
         call fi32%set_scalar_i32(gid, int(7_id_k * gid, i32_k), s)
         call fbit%set_scalar_bit(gid, mod(gid, 3_id_k) == 0_id_k, s)
         call fvec%set_tuple_r64(gid, [real(gid, r64_k), 0.5_r64_k, -real(gid, r64_k)], s)
         call fbits3%set_tuple_bit(gid, [.true., mod(gid, 2_id_k) == 1_id_k, .false.], s)
      end do

      call g%export_breadth_first(bf, s)
      call bf%gather(f64, v64, s)
      call check(s%ok() .and. size(v64, kind=id_k) == bf%number_of_cells(), 'gather r64')
      call bf%gather(f32, v32, s)
      call check(s%ok(), 'gather r32')
      call bf%gather(fi64, vi64, s)
      call check(s%ok(), 'gather i64')
      call bf%gather(fi32, vi32, s)
      call check(s%ok(), 'gather i32')
      call bf%gather(fbit, vbit, s)
      call check(s%ok(), 'gather bit -> int8')
      call bf%gather(fvec, vvec, s)
      call check(s%ok() .and. all(shape(vvec) == [3, int(bf%number_of_cells())]), &
         'gather 3-component r64 as (ncomp, ncell)')
      call bf%gather(fi32, vi32_2, s)
      call check(s%ok() .and. all(shape(vi32_2) == [1, int(bf%number_of_cells())]), &
         'gather scalar as rank-2 (1, ncell)')
      call bf%gather(fbits3, vbit3, s)
      call check(s%ok() .and. size(vbit3, 1) == 3, 'gather 3-component bit -> int8')

      do n = 1_id_k, bf%number_of_cells()
         gid = bf%global_ids(n)
         call check(same(v64(n), 1.5_r64_k * real(gid, r64_k)), 'r64 values permuted')
         call check(same(real(v32(n), r64_k), -real(gid, r64_k)), 'r32 values permuted')
         call check(vi64(n) == 1000_i64_k + gid, 'i64 values permuted')
         call check(vi32(n) == int(7_id_k * gid, i32_k), 'i32 values permuted')
         call check(vi32_2(1, n) == vi32(n), 'rank-2 scalar matches rank-1')
         call check(vbit(n) == merge(1_byte_k, 0_byte_k, mod(gid, 3_id_k) == 0_id_k), &
            'bit values permuted as 0/1')
         expected_vec = [real(gid, r64_k), 0.5_r64_k, -real(gid, r64_k)]
         do comp = 1_i32_k, 3_i32_k
            call check(same(vvec(comp, n), expected_vec(comp)), 'vector components permuted')
         end do
         call check(all(vbit3(:, n) == [1_byte_k, merge(1_byte_k, 0_byte_k, &
            mod(gid, 2_id_k) == 1_id_k), 0_byte_k]), 'bit tuples permuted')
      end do

      ! Masking prunes cells; gather follows the shorter export
      call g%set_cell_mask(16_id_k, .true., s)
      call g%export_breadth_first(bf, s)
      call bf%gather(f64, v64, s)
      call check(s%ok() .and. size(v64) == 19, 'gather follows a pruned export')
   end subroutine test_gather

   subroutine test_errors()
      type(htg_grid_t), target :: g
      type(htg_breadth_first_t) :: bf
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      type(htg_field_array_t) :: short, vec, uninit
      real(r64_k), allocatable :: v64(:)
      real(r32_k), allocatable :: v32(:)

      call g%initialize()
      call g%set_dimensions(3_i32_k, 2_i32_k, 1_i32_k, s)
      call g%set_x_coordinates([0.0_r64_k, 1.0_r64_k, 2.0_r64_k], s)
      call g%set_y_coordinates([0.0_r64_k, 1.0_r64_k], s)
      call new_non_oriented_cursor(g, 1_id_k, c, create=.true., status=s)
      call c%subdivide_leaf(s)

      call g%export_breadth_first(bf, s)
      call check(s%code == HTG_ERR_INDEXING, 'export without global indices is rejected')

      call g%initialize_local_index_node(s)
      call g%export_breadth_first(bf, s)
      call check(s%ok(), 'export after indexing')

      call short%initialize('short', HTG_TYPE_REAL64, 1_i32_k, 3_id_k, s)
      call bf%gather(short, v64, s)
      call check(s%code == HTG_ERR_OUT_OF_BOUNDS, 'gather rejects too few tuples')
      call check(.not. allocated(v64), 'failed gather leaves output unallocated')

      call vec%initialize('vec', HTG_TYPE_REAL64, 3_i32_k, 5_id_k, s)
      call bf%gather(vec, v64, s)
      call check(s%code == HTG_ERR_INVALID_ARG, 'rank-1 gather rejects a multi-component field')

      call bf%gather(vec, v32, s)
      call check(s%code == HTG_ERR_INVALID_ARG, 'gather rejects a kind mismatch')

      call bf%gather(uninit, v64, s)
      call check(s%code == HTG_ERR_INVALID_STATE, 'gather rejects an uninitialized field')
   end subroutine test_errors

end program test_export
