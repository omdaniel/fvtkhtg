!> Cell geometry and face-neighbour lookup for htg_cursor_t.
!>
!> GEOMETRY
!> Mirrors vtkHyperTreeGridNonOrientedGeometryCursor (GetOrigin, GetSize,
!> GetBounds, GetPoint) and reproduces VTK's floating-point results exactly,
!> because it performs the same operations in the same order:
!>   - the root cell's origin and size come from the grid's coordinates;
!>   - each level divides the size by the branch factor (vtkHyperTreeGridScales);
!>   - each step to child `ichild` adds digit * size along every active axis,
!>     where digit p (base branch_factor) of `ichild` belongs to the p-th
!>     active axis in ASCENDING order (x before y before z).
!> Nothing is stored in the cursor: the geometry is recomputed from the
!> cursor's path on each query, in O(level).
!>
!> Note on xz-plane grids. vtkHyperTreeGrid::GetAxes (grid%get_axes) lists the
!> axes of an xz-plane grid as (z, x), but VTK's geometry cursors number the
!> children x first in every plane. Child numbering therefore follows
!> grid%get_active_axes, never grid%get_axes.
!>
!> FACE NEIGHBOURS
!> cursor%get_face_neighbor(axis, direction, neighbor) finds the cell on the
!> other side of one face of the current cell:
!>   - if the region across the face is covered by a cell of the same level,
!>     `neighbor` is that cell (it may itself be refined: descend from it to
!>     reach the smaller cells that touch the face);
!>   - if it is covered by a coarser leaf, `neighbor` is that leaf, so
!>     neighbor%get_level() < cursor%get_level();
!>   - the search does not descend through a masked cell: a masked cell on the
!>     way is returned as the neighbour (test it with neighbor%is_masked());
!>   - at the edge of the grid there is no neighbour: neighbor%has_tree() is
!>     .false. and neighbor%get_tree_id() is HTG_INVALID_ID;
!>   - where the neighbouring root cell has no tree, neighbor%has_tree() is
!>     .false. and neighbor%get_tree_id() is that root cell's index.
!> The result is a full cursor, so its geometry, global index and children are
!> available, and it can refine the neighbour (2:1 balancing).
!>
!> This is what VTK's Von Neumann supercursor provides for each face, except
!> that fvtkhtg defines the direction geometrically. VTK 9.7.2's supercursor
!> takes its root-level neighbours along vtkHyperTreeGrid::GetAxes, which for
!> an xz-plane grid disagrees with its own child numbering, so across root
!> boundaries of such grids VTK reports the neighbour along the other axis.
submodule (fvtkhtg_cursor_m) fvtkhtg_cursor_geometry_s
   implicit none
contains

   !> Origin and size of the cell under the cursor. A cursor on a root cell
   !> that has no tree still has the geometry of that root cell.
   subroutine cell_geometry(self, origin, cell_size, caller, status)
      class(htg_cursor_t), intent(in) :: self
      real(r64_k), intent(out) :: origin(0:2)
      real(r64_k), intent(out) :: cell_size(0:2)
      character(len=*), intent(in) :: caller
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: st
      integer(i32_k) :: active(0:2), nactive, lvl, p, bf, ichild, stride, digit
      integer(id_k) :: elder
      real(r64_k) :: factor

      origin = 0.0_r64_k
      cell_size = 0.0_r64_k
      if (.not. associated(self%grid) .or. self%tree_id < 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%' // caller // ': cursor is not on a root cell of a grid', status)
         return
      end if

      call self%grid%get_level_zero_origin_and_size_from_index(self%tree_id, origin, cell_size, st)
      if (.not. st%ok()) then
         call htg_fail(st%code, 'htg_cursor_t%' // caller // ': ' // trim(st%message), status)
         return
      end if

      if (associated(self%tree)) then
         bf = self%grid%get_branch_factor()
         factor = real(bf, r64_k)
         call self%grid%get_active_axes(active, nactive)
         do lvl = 1_i32_k, self%last_valid_entry
            elder = self%tree%get_elder_child_index(self%entries(lvl - 1_i32_k)%vertex_id, st)
            if (.not. st%ok()) then
               call htg_fail(st%code, 'htg_cursor_t%' // caller // ': ' // trim(st%message), status)
               return
            end if
            ichild = int(self%entries(lvl)%vertex_id - elder, i32_k)

            cell_size = cell_size / factor
            stride = 1_i32_k
            do p = 0_i32_k, nactive - 1_i32_k
               digit = mod(ichild / stride, bf)
               origin(active(p)) = origin(active(p)) + real(digit, r64_k) * cell_size(active(p))
               stride = stride * bf
            end do
         end do
      end if

      call htg_set_ok(status)
   end subroutine cell_geometry

   module procedure cursor_get_origin
      real(r64_k) :: cell_size(0:2)

      call cell_geometry(self, origin, cell_size, 'get_origin', status)
   end procedure cursor_get_origin

   module procedure cursor_get_size
      real(r64_k) :: origin(0:2)

      call cell_geometry(self, origin, cell_size, 'get_size', status)
   end procedure cursor_get_size

   module procedure cursor_get_bounds
      real(r64_k) :: origin(0:2), cell_size(0:2)
      integer :: a

      call cell_geometry(self, origin, cell_size, 'get_bounds', status)
      do a = 0, 2
         bounds(2*a) = origin(a)
         bounds(2*a + 1) = origin(a) + cell_size(a)
      end do
   end procedure cursor_get_bounds

   module procedure cursor_get_point
      real(r64_k) :: origin(0:2), cell_size(0:2)

      call cell_geometry(self, origin, cell_size, 'get_point', status)
      point = origin + cell_size / 2.0_r64_k
   end procedure cursor_get_point

   module procedure cursor_get_face_neighbor
      type(htg_status_t) :: st
      type(htg_tree_t), pointer :: nb_tree
      integer(i32_k) :: active(0:2), nactive, p, pos, bf, stride, lvl, level, digit
      integer(i32_k) :: carry_level, start_level
      integer(i32_k), allocatable :: child(:)
      integer(id_k) :: ijk(0:2), nb_tree_id, elder, cur, gid
      integer(i32_k) :: cell_dims(0:2)
      logical :: leaf, masked, has_mask

      ! `neighbor` starts out as "no neighbour"; its entry storage is reused
      neighbor%grid => self%grid
      neighbor%tree => null()
      neighbor%tree_id = HTG_INVALID_ID
      neighbor%last_valid_entry = -1_i32_k

      if (.not. associated(self%grid) .or. .not. associated(self%tree) .or. &
          self%last_valid_entry < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_STATE, &
            'htg_cursor_t%get_face_neighbor: cursor has no active tree', status)
         return
      end if

      if (direction /= -1_i32_k .and. direction /= 1_i32_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_cursor_t%get_face_neighbor: direction must be -1 or +1', status)
         return
      end if

      call self%grid%get_active_axes(active, nactive)
      pos = -1_i32_k
      do p = 0_i32_k, nactive - 1_i32_k
         if (active(p) == axis) pos = p
      end do
      if (pos < 0_i32_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_cursor_t%get_face_neighbor: axis must be one of the grid''s non-squeezed axes', &
            status)
         return
      end if

      bf = self%grid%get_branch_factor()
      stride = bf ** pos
      level = self%last_valid_entry

      ! The cursor's path as child indices, one per level below the root
      allocate(child(1:level))
      do lvl = 1_i32_k, level
         elder = self%tree%get_elder_child_index(self%entries(lvl - 1_i32_k)%vertex_id, st)
         if (.not. st%ok()) then
            call htg_fail(st%code, 'htg_cursor_t%get_face_neighbor: ' // trim(st%message), status)
            return
         end if
         child(lvl) = int(self%entries(lvl)%vertex_id - elder, i32_k)
      end do

      ! The digits of the child indices along `axis` form the cell's integer
      ! coordinate in base branch_factor, most significant digit at level 1.
      ! Step it by `direction`, carrying towards the root. carry_level is the
      ! level that absorbs the carry; 0 means the step leaves the tree.
      carry_level = 0_i32_k
      do lvl = level, 1_i32_k, -1_i32_k
         digit = mod(child(lvl) / stride, bf) + direction
         if (digit >= 0_i32_k .and. digit < bf) then
            child(lvl) = child(lvl) + direction * stride
            carry_level = lvl
            exit
         end if
         ! Wrap this digit (bf-1 -> 0 going up, 0 -> bf-1 going down) and carry on
         child(lvl) = child(lvl) - direction * (bf - 1_i32_k) * stride
      end do

      if (carry_level == 0_i32_k) then
         ! Across a root boundary: the neighbouring root cell, if there is one
         call self%grid%get_level_zero_coordinates_from_index(self%tree_id, &
            ijk(0), ijk(1), ijk(2), st)
         if (.not. st%ok()) then
            call htg_fail(st%code, 'htg_cursor_t%get_face_neighbor: ' // trim(st%message), status)
            return
         end if
         cell_dims = self%grid%get_cell_dims()
         ijk(axis) = ijk(axis) + int(direction, id_k)
         if (ijk(axis) < 0_id_k .or. ijk(axis) >= int(cell_dims(axis), id_k)) then
            call htg_set_ok(status)   ! edge of the grid: no neighbour
            return
         end if
         call self%grid%get_index_from_level_zero_coordinates(nb_tree_id, &
            ijk(0), ijk(1), ijk(2), st)
         if (.not. st%ok()) then
            call htg_fail(st%code, 'htg_cursor_t%get_face_neighbor: ' // trim(st%message), status)
            return
         end if
         neighbor%tree_id = nb_tree_id
         call self%grid%get_tree(nb_tree_id, nb_tree, create=.false., status=st)
         if (.not. st%ok()) then
            call htg_fail(st%code, 'htg_cursor_t%get_face_neighbor: ' // trim(st%message), status)
            return
         end if
         if (.not. associated(nb_tree)) then
            call htg_set_ok(status)   ! that root cell has no tree
            return
         end if
         neighbor%tree => nb_tree
         call ensure_entry_capacity(neighbor, level)
         neighbor%entries(0)%vertex_id = 0_id_k
         start_level = 0_i32_k
      else
         ! Same tree: the path is shared down to the level above the carry
         neighbor%tree => self%tree
         neighbor%tree_id = self%tree_id
         call ensure_entry_capacity(neighbor, level)
         start_level = carry_level - 1_i32_k
         neighbor%entries(0_i32_k : start_level) = self%entries(0_i32_k : start_level)
      end if
      neighbor%last_valid_entry = start_level

      ! Follow the stepped path down, stopping at a leaf (a coarser
      ! neighbour) or at a masked cell, whose subtree is hidden
      has_mask = self%grid%has_mask()
      do lvl = start_level + 1_i32_k, level
         cur = neighbor%entries(lvl - 1_i32_k)%vertex_id
         leaf = neighbor%tree%is_leaf(cur, st)
         if (leaf) exit
         if (has_mask) then
            gid = neighbor%tree%get_global_index_from_local(cur, st)
            if (.not. st%ok()) then
               call htg_fail(st%code, 'htg_cursor_t%get_face_neighbor: ' // trim(st%message), &
                  status)
               return
            end if
            masked = self%grid%is_cell_masked(gid, st)
            if (masked) exit
         end if
         elder = neighbor%tree%get_elder_child_index(cur, st)
         neighbor%entries(lvl)%vertex_id = elder + int(child(lvl), id_k)
         neighbor%last_valid_entry = lvl
      end do

      call htg_set_ok(status)
   end procedure cursor_get_face_neighbor

end submodule fvtkhtg_cursor_geometry_s
