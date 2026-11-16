!> Breadth-first export of an htg_grid_t (see fvtkhtg_breadth_first_m).
!>
!> Matches what VTK's own vtkHDFWriter stores for a vtkHyperTreeGrid, which
!> calls vtkHyperTree::ComputeBreadthFirstOrderDescriptor with the grid's
!> depth limiter and mask:
!>   - a masked cell is exported as a leaf and its descendants are dropped
!>     (so no refined cell is ever masked, as the fvtkhdf writer requires);
!>   - a cell at the depth limiter is exported as a leaf and its descendants
!>     are dropped;
!>   - the mask, when exported, holds each exported cell's mask bit.
submodule (fvtkhtg_grid_m) fvtkhtg_grid_export_s
   implicit none
contains

   module procedure grid_export_breadth_first
      type(htg_status_t) :: st
      integer(id_k) :: ntree, ncap, n, pos, v, gid, elder, i, ncur, nnxt, c, k
      integer(i32_k) :: depth
      logical :: has_mask, masked, leaf, is_refined
      integer(id_k), allocatable :: cur(:), nxt(:), ids(:)
      logical, allocatable :: ref(:), msk(:)

      ntree = self%num_trees
      has_mask = self%has_mask()

      do pos = 0_id_k, ntree - 1_id_k
         if (self%tree_slots(pos)%tree%get_indexing_mode() == HTG_INDEX_UNSET) then
            call htg_fail(HTG_ERR_INDEXING, &
               'htg_grid_t%export_breadth_first: every tree needs global indices ' // &
               '(call initialize_local_index_node or set them explicitly)', status)
            return
         end if
      end do

      ncap = max(self%get_number_of_vertices(), 1_id_k)
      allocate(bf%tree_ids(ntree), bf%xtree(ntree + 1_id_k))
      allocate(ids(ncap), ref(ncap), msk(ncap), cur(ncap), nxt(ncap))

      n = 0_id_k
      do pos = 0_id_k, ntree - 1_id_k
         associate (tree => self%tree_slots(pos)%tree)
            bf%tree_ids(pos + 1_id_k) = self%tree_slots(pos)%tree_id + 1_id_k
            bf%xtree(pos + 1_id_k) = n + 1_id_k
            k = int(tree%get_number_of_children(), id_k)

            ncur = 1_id_k
            cur(1) = 0_id_k
            depth = 0_i32_k
            do while (ncur > 0_id_k)
               nnxt = 0_id_k
               do i = 1_id_k, ncur
                  v = cur(i)
                  gid = tree%get_global_index_from_local(v, st)
                  if (.not. st%ok()) then
                     call htg_fail(st%code, 'htg_grid_t%export_breadth_first: ' // &
                        trim(st%message), status)
                     return
                  end if

                  masked = .false.
                  if (has_mask) then
                     masked = self%is_cell_masked(gid, st)
                     if (.not. st%ok()) then
                        call htg_fail(st%code, 'htg_grid_t%export_breadth_first: ' // &
                           trim(st%message), status)
                        return
                     end if
                  end if

                  leaf = tree%is_leaf(v, st)
                  is_refined = .not. leaf
                  if (masked) is_refined = .false.
                  if (depth >= self%depth_limiter) is_refined = .false.

                  n = n + 1_id_k
                  ids(n) = gid
                  ref(n) = is_refined
                  msk(n) = masked

                  if (is_refined) then
                     elder = tree%get_elder_child_index(v, st)
                     do c = 0_id_k, k - 1_id_k
                        nnxt = nnxt + 1_id_k
                        nxt(nnxt) = elder + c
                     end do
                  end if
               end do
               cur(1:nnxt) = nxt(1:nnxt)
               ncur = nnxt
               depth = depth + 1_i32_k
            end do
         end associate
      end do
      bf%xtree(ntree + 1_id_k) = n + 1_id_k

      bf%global_ids = ids(1:n)
      bf%refined = ref(1:n)

      if (has_mask) then
         bf%mask = msk(1:n)
      else if (present(force_mask)) then
         ! fvtkhdf requires all MPI ranks to pass a mask or none
         if (force_mask) then
            allocate(bf%mask(n))
            bf%mask = .false.
         end if
      end if

      bf%identity_order = .true.
      do i = 1_id_k, n
         if (bf%global_ids(i) /= i - 1_id_k) then
            bf%identity_order = .false.
            exit
         end if
      end do

      call htg_set_ok(status)
   end procedure grid_export_breadth_first

end submodule fvtkhtg_grid_export_s
