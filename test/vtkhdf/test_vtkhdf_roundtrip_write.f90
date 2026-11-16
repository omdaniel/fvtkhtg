!> Round trip, part 1: build grids with fvtkhtg, write them as VTKHDF
!> HyperTreeGrid files through the fvtkhtg_vtkhdf bridge, and record next to
!> each file what every cell should be.
!>
!> Part 2 (check_vtkhdf_roundtrip.py) reads the files with VTK's own
!> vtkHDFReader and compares cell by cell: structure, geometry, mask, and
!> every cell-data value.
!>
!> Output, in the working directory:
!>   manifest.txt          one line per check: <file> <expectation> <step or -1>
!>   <case>.vtkhdf         the file under test
!>   <case>[.stepN].expect what fvtkhtg says the file holds:
!>                           header lines (dims, coordinates, arrays, ...)
!>                           cell <tree> <path> <level> <refined> <masked> <6 bounds> <values...>
!>                         <path> is the child indices from the root, joined by
!>                         '.', or '-' for a root cell. Cells are listed as
!>                         exported: a masked cell or a cell at the depth
!>                         limiter is a leaf and its descendants are absent.
!>
!> The "standard" arrays hold values derived from each cell's own position
!> in its tree, so the checker can recompute them independently:
!>   h        = fold over the path of  h*29 + child + 1,  starting at tree_id + 1
!>   Pressure = 0.25*h + 0.5                  (REAL64)
!>   Velocity = (h, -h, 0.5*level)            (REAL32, 3 components)
!>   Id       = h                             (INT64)
!>   Count    = (level, last child or -1)     (INT32, 2 components)
!>   Flag     = mod(h, 3) == 0                (BIT, written as 0/1 bytes)
program test_vtkhdf_roundtrip_write
   use fvtkhtg
   use fvtkhtg_vtkhdf
   use vtkhdf_htg_file_type, only: vtkhdf_htg_file, vtkhdf_cell_data_handle, &
      VTKHDF_STATIC, VTKHDF_DYNAMIC_MESH
   implicit none

   integer, parameter :: MAX_PATH = 16
   integer(i64_k) :: rng_state = 4711_i64_k
   integer :: manifest_unit
   integer :: ncases = 0

   open (newunit=manifest_unit, file='manifest.txt', status='replace', action='write')

   call case_static_2d()
   call case_static_3d_ternary_explicit_ids()
   call case_static_xz_interface()
   call case_static_1d_descending()
   call case_depth_limited()
   call case_coarsened()
   call case_dynamic()

   close (manifest_unit)
   print '(a,i0,a)', 'ALL ', ncases, ' ROUND-TRIP FILES WRITTEN'

contains

   subroutine check(cond, msg)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) error stop 'FAIL: ' // msg
   end subroutine check

   subroutine check_status(s, what)
      type(htg_status_t), intent(in) :: s
      character(len=*), intent(in) :: what
      if (.not. s%ok()) error stop 'FAIL: ' // what // ': ' // trim(s%message)
   end subroutine check_status

   function rand01() result(r)
      real(r64_k) :: r
      rng_state = mod(rng_state * 1103515245_i64_k + 12345_i64_k, 2147483648_i64_k)
      r = real(rng_state, r64_k) / 2147483648.0_r64_k
   end function rand01

   function stretched(n, x0, dx, stretch) result(c)
      integer(i32_k), intent(in) :: n
      real(r64_k), intent(in) :: x0, dx, stretch
      real(r64_k), allocatable :: c(:)
      integer(i32_k) :: i
      allocate(c(0 : n - 1_i32_k))
      do i = 0_i32_k, n - 1_i32_k
         c(i) = x0 + dx * real(i, r64_k) + stretch * real(i * i, r64_k)
      end do
   end function stretched

   pure function path_hash(tree_id, path, n) result(h)
      integer(id_k), intent(in) :: tree_id
      integer(i32_k), intent(in) :: path(:)
      integer, intent(in) :: n
      integer(i64_k) :: h
      integer :: l
      h = tree_id + 1_i64_k
      do l = 1, n
         h = h * 29_i64_k + int(path(l), i64_k) + 1_i64_k
      end do
   end function path_hash

   function path_string(path, n) result(str)
      integer(i32_k), intent(in) :: path(:)
      integer, intent(in) :: n
      character(len=:), allocatable :: str
      character(len=12) :: digits
      integer :: l
      if (n == 0) then
         str = '-'
         return
      end if
      str = ''
      do l = 1, n
         write (digits, '(i0)') path(l)
         if (l > 1) str = str // '.'
         str = str // trim(digits)
      end do
   end function path_string

   subroutine make_grid(g, dims, bf, transposed)
      type(htg_grid_t), target, intent(inout) :: g
      integer(i32_k), intent(in) :: dims(3), bf
      logical, intent(in) :: transposed
      type(htg_status_t) :: s

      call g%initialize()
      call g%set_branch_factor(bf, s)
      call g%set_transposed_root_indexing(transposed, s)
      call g%set_dimensions(dims(1), dims(2), dims(3), s)
      call check_status(s, 'set_dimensions')
      call g%set_x_coordinates(stretched(dims(1), -1.0_r64_k, 1.0_r64_k, 0.125_r64_k), s)
      call g%set_y_coordinates(stretched(dims(2), 2.0_r64_k, 0.5_r64_k, 0.0625_r64_k), s)
      call g%set_z_coordinates(stretched(dims(3), 10.0_r64_k, 2.0_r64_k, 0.25_r64_k), s)
      call check_status(s, 'coordinates')
   end subroutine make_grid

   recursive subroutine refine_random(c, max_level, prob)
      type(htg_cursor_t), intent(inout) :: c
      integer(i32_k), intent(in) :: max_level
      real(r64_k), intent(in) :: prob
      type(htg_status_t) :: s
      integer(i32_k) :: ic

      if (c%get_level() >= max_level) return
      if (rand01() >= prob) return
      call c%subdivide_leaf(s)
      call check_status(s, 'subdivide_leaf')
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         call refine_random(c, max_level, prob)
         call c%to_parent(s)
      end do
   end subroutine refine_random

   !> Create every tree except `skip_tree` and refine it at random.
   subroutine grow_forest(g, max_level, prob, skip_tree)
      type(htg_grid_t), target, intent(inout) :: g
      integer(i32_k), intent(in) :: max_level
      real(r64_k), intent(in) :: prob
      integer(id_k), intent(in) :: skip_tree
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      integer(id_k) :: tid

      do tid = 0_id_k, g%get_max_number_of_trees() - 1_id_k
         if (tid == skip_tree) cycle
         call new_non_oriented_cursor(g, tid, c, create=.true., status=s)
         call check_status(s, 'create tree')
         call refine_random(c, max_level, prob)
      end do
   end subroutine grow_forest

   ! ---------------------------------------------------------------------------
   ! Standard arrays
   ! ---------------------------------------------------------------------------

   subroutine add_standard_fields(g, ntuples)
      type(htg_grid_t), target, intent(inout) :: g
      integer(id_k), intent(in) :: ntuples
      type(htg_cell_data_t), pointer :: cd
      type(htg_status_t) :: s

      call g%get_cell_data(cd)
      call cd%add_array('Pressure', HTG_TYPE_REAL64, 1_i32_k, ntuples, s)
      call cd%add_array('Velocity', HTG_TYPE_REAL32, 3_i32_k, ntuples, s)
      call cd%add_array('Id', HTG_TYPE_INT64, 1_i32_k, ntuples, s)
      call cd%add_array('Count', HTG_TYPE_INT32, 2_i32_k, ntuples, s)
      call cd%add_array('Flag', HTG_TYPE_BIT, 1_i32_k, ntuples, s)
      call check_status(s, 'add standard fields')
   end subroutine add_standard_fields

   subroutine fill_standard_fields(g)
      type(htg_grid_t), target, intent(inout) :: g
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      integer(id_k) :: tid
      integer(i32_k) :: path(MAX_PATH)

      do tid = 0_id_k, g%get_max_number_of_trees() - 1_id_k
         if (.not. g%has_tree(tid)) cycle
         call new_non_oriented_cursor(g, tid, c, create=.false., status=s)
         call fill_cell(g, c, tid, path, 0)
      end do
   end subroutine fill_standard_fields

   !> Set the standard values of the cell under the cursor.
   subroutine store_standard(g, c, tree_id, path, n)
      type(htg_grid_t), target, intent(inout) :: g
      type(htg_cursor_t), intent(in) :: c
      integer(id_k), intent(in) :: tree_id
      integer(i32_k), intent(in) :: path(:)
      integer, intent(in) :: n
      type(htg_cell_data_t), pointer :: cd
      type(htg_field_array_t), pointer :: f
      type(htg_status_t) :: s
      integer(i64_k) :: h
      integer(id_k) :: gid
      integer(i32_k) :: last_child

      h = path_hash(tree_id, path, n)
      gid = c%get_global_node_index(s)
      call check_status(s, 'global id of a cell')
      last_child = -1_i32_k
      if (n > 0) last_child = path(n)

      call g%get_cell_data(cd)
      call cd%get_array('Pressure', f, s)
      call f%set_scalar_r64(gid, 0.25_r64_k * real(h, r64_k) + 0.5_r64_k, s)
      call cd%get_array('Velocity', f, s)
      call f%set_tuple_r32(gid, [real(h, r32_k), -real(h, r32_k), 0.5_r32_k * real(n, r32_k)], s)
      call cd%get_array('Id', f, s)
      call f%set_scalar_i64(gid, h, s)
      call cd%get_array('Count', f, s)
      call f%set_tuple_i32(gid, [int(n, i32_k), last_child], s)
      call cd%get_array('Flag', f, s)
      call f%set_scalar_bit(gid, mod(h, 3_i64_k) == 0_i64_k, s)
      call check_status(s, 'store standard values')
   end subroutine store_standard

   recursive subroutine fill_cell(g, c, tree_id, path, n)
      type(htg_grid_t), target, intent(inout) :: g
      type(htg_cursor_t), intent(inout) :: c
      integer(id_k), intent(in) :: tree_id
      integer(i32_k), intent(inout) :: path(:)
      integer, intent(in) :: n
      type(htg_status_t) :: s
      integer(i32_k) :: ic
      logical :: leaf

      call store_standard(g, c, tree_id, path, n)
      leaf = c%is_leaf(s)
      if (leaf) return
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         path(n + 1) = ic
         call fill_cell(g, c, tree_id, path, n + 1)
         call c%to_parent(s)
      end do
   end subroutine fill_cell

   !> Mask roughly one cell in eleven, leaves and refined cells alike.
   subroutine mask_some(g, nmasked)
      type(htg_grid_t), target, intent(inout) :: g
      integer, intent(out) :: nmasked
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      integer(id_k) :: tid
      integer(i32_k) :: path(MAX_PATH)

      nmasked = 0
      do tid = 0_id_k, g%get_max_number_of_trees() - 1_id_k
         if (.not. g%has_tree(tid)) cycle
         call new_non_oriented_cursor(g, tid, c, create=.false., status=s)
         call mask_cell(c, tid, path, 0, nmasked)
      end do
   end subroutine mask_some

   recursive subroutine mask_cell(c, tree_id, path, n, nmasked)
      type(htg_cursor_t), intent(inout) :: c
      integer(id_k), intent(in) :: tree_id
      integer(i32_k), intent(inout) :: path(:)
      integer, intent(in) :: n
      integer, intent(inout) :: nmasked
      type(htg_status_t) :: s
      integer(i32_k) :: ic
      logical :: leaf

      if (n > 0 .and. mod(path_hash(tree_id, path, n), 11_i64_k) == 5_i64_k) then
         call c%set_mask(.true., s)
         call check_status(s, 'set_mask')
         nmasked = nmasked + 1
         return                       ! nothing below a masked cell is visible
      end if
      leaf = c%is_leaf(s)
      if (leaf) return
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         path(n + 1) = ic
         call mask_cell(c, tree_id, path, n + 1, nmasked)
         call c%to_parent(s)
      end do
   end subroutine mask_cell

   ! ---------------------------------------------------------------------------
   ! Expectation files
   ! ---------------------------------------------------------------------------

   subroutine write_coordinate_line(u, label, coords)
      integer, intent(in) :: u
      character(len=*), intent(in) :: label
      real(r64_k), intent(in) :: coords(:)
      integer :: i
      write (u, '(a,1x,i0)', advance='no') label, size(coords)
      do i = 1, size(coords)
         write (u, '(1x,es25.17e3)', advance='no') coords(i)
      end do
      write (u, '()')
   end subroutine write_coordinate_line

   !> Describe the grid as it will be exported. Returns the number of cells.
   subroutine write_expectation(g, filename, case_name, ncells)
      type(htg_grid_t), target, intent(inout) :: g
      character(len=*), intent(in) :: filename, case_name
      integer(id_k), intent(out) :: ncells
      type(htg_cursor_t) :: c
      type(htg_cell_data_t), pointer :: cd
      type(htg_field_array_t), pointer :: f
      type(htg_breadth_first_t) :: bf
      type(htg_status_t) :: s
      real(r64_k), pointer :: xc(:), yc(:), zc(:)
      integer(i32_k) :: path(MAX_PATH), i, dims(0:2)
      integer(id_k) :: tid
      integer :: u
      character(len=:), allocatable :: active
      character(len=6) :: type_name

      call g%get_cell_data(cd)
      call g%get_x_coordinates_view(xc)
      call g%get_y_coordinates_view(yc)
      call g%get_z_coordinates_view(zc)
      dims = g%get_dimensions()

      open (newunit=u, file=filename, status='replace', action='write')
      write (u, '(a,1x,a)') 'case', case_name
      write (u, '(a,3(1x,i0))') 'dims', dims
      write (u, '(a,1x,i0)') 'branch_factor', g%get_branch_factor()
      write (u, '(a,1x,i0)') 'transposed', merge(1, 0, g%get_transposed_root_indexing())
      if (g%get_has_interface()) then
         write (u, '(a,1x,a,1x,a)') 'interface 1', g%get_interface_normals_name(), &
            g%get_interface_intercepts_name()
      else
         write (u, '(a)') 'interface 0'
      end if
      call write_coordinate_line(u, 'xcoords', xc)
      call write_coordinate_line(u, 'ycoords', yc)
      call write_coordinate_line(u, 'zcoords', zc)
      active = cd%get_active_scalars()
      if (len(active) == 0) active = '-'
      write (u, '(a,1x,a)') 'active_scalars', active
      write (u, '(a,1x,i0)') 'has_mask', merge(1, 0, g%has_mask())
      do i = 0_i32_k, cd%get_number_of_arrays() - 1_i32_k
         call cd%get_array_by_index(i, f, s)
         select case (f%get_data_type())
         case (HTG_TYPE_REAL64)
            type_name = 'REAL64'
         case (HTG_TYPE_REAL32)
            type_name = 'REAL32'
         case (HTG_TYPE_INT64)
            type_name = 'INT64'
         case (HTG_TYPE_INT32)
            type_name = 'INT32'
         case default
            type_name = 'BIT'
         end select
         write (u, '(a,1x,a,1x,a,1x,i0)') 'array', f%get_name(), trim(type_name), &
            f%get_number_of_components()
      end do

      ncells = 0_id_k
      do tid = 0_id_k, g%get_max_number_of_trees() - 1_id_k
         if (.not. g%has_tree(tid)) cycle
         call new_non_oriented_cursor(g, tid, c, create=.false., status=s)
         call write_cell(g, c, u, tid, path, 0, ncells)
      end do
      write (u, '(a,1x,i0)') 'ncells', ncells
      close (u)

      ! The description must cover exactly the cells the export hands over
      call g%export_breadth_first(bf, s)
      call check_status(s, 'export_breadth_first')
      call check(bf%number_of_cells() == ncells, case_name // ': expectation covers the exported cells')
   end subroutine write_expectation

   recursive subroutine write_cell(g, c, u, tree_id, path, n, ncells)
      type(htg_grid_t), target, intent(inout) :: g
      type(htg_cursor_t), intent(inout) :: c
      integer, intent(in) :: u
      integer(id_k), intent(in) :: tree_id
      integer(i32_k), intent(inout) :: path(:)
      integer, intent(in) :: n
      integer(id_k), intent(inout) :: ncells
      type(htg_cell_data_t), pointer :: cd
      type(htg_field_array_t), pointer :: f
      type(htg_status_t) :: s
      real(r64_k) :: bounds(0:5), v64
      real(r32_k) :: v32
      integer(i64_k) :: vi64
      integer(i32_k) :: vi32, ic, i, comp
      integer(id_k) :: gid
      logical :: leaf, masked, refined, bit

      leaf = c%is_leaf(s)
      masked = c%is_masked(s)
      refined = (.not. leaf) .and. (.not. masked)
      call c%get_bounds(bounds, s)
      call check_status(s, 'cell bounds')
      gid = c%get_global_node_index(s)
      ncells = ncells + 1_id_k

      write (u, '(a,1x,i0,1x,a,1x,i0,1x,i0,1x,i0,6(1x,es25.17e3))', advance='no') &
         'cell', tree_id, path_string(path, n), n, merge(1, 0, refined), merge(1, 0, masked), bounds

      call g%get_cell_data(cd)
      do i = 0_i32_k, cd%get_number_of_arrays() - 1_i32_k
         call cd%get_array_by_index(i, f, s)
         do comp = 0_i32_k, f%get_number_of_components() - 1_i32_k
            select case (f%get_data_type())
            case (HTG_TYPE_REAL64)
               v64 = f%get_component_r64(gid, comp, s)
               write (u, '(1x,es25.17e3)', advance='no') v64
            case (HTG_TYPE_REAL32)
               v32 = f%get_component_r32(gid, comp, s)
               write (u, '(1x,es25.17e3)', advance='no') real(v32, r64_k)
            case (HTG_TYPE_INT64)
               vi64 = f%get_component_i64(gid, comp, s)
               write (u, '(1x,i0)', advance='no') vi64
            case (HTG_TYPE_INT32)
               vi32 = f%get_component_i32(gid, comp, s)
               write (u, '(1x,i0)', advance='no') vi32
            case default
               bit = f%get_component_bit(gid, comp, s)
               write (u, '(1x,i0)', advance='no') merge(1, 0, bit)
            end select
            call check_status(s, 'read back a cell value')
         end do
      end do
      write (u, '()')

      if (.not. refined) return
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         call check_status(s, 'to_child while describing the grid')
         path(n + 1) = ic
         call write_cell(g, c, u, tree_id, path, n + 1, ncells)
         call c%to_parent(s)
      end do
   end subroutine write_cell

   subroutine open_file(file, name, mode)
      type(vtkhdf_htg_file), intent(out) :: file
      character(len=*), intent(in) :: name
      integer, intent(in) :: mode
      integer :: stat
      character(len=:), allocatable :: errmsg

      call file%create(name, stat, errmsg, mode=mode)
      if (stat /= 0) error stop 'FAIL: cannot create ' // name // ': ' // errmsg
   end subroutine open_file

   !> Validate, write a static file through the bridge, and describe it.
   subroutine write_static_case(g, case_name)
      type(htg_grid_t), target, intent(inout) :: g
      character(len=*), intent(in) :: case_name
      type(vtkhdf_htg_file) :: file
      type(htg_status_t) :: s
      integer(id_k) :: ncells

      call g%validate(s)
      call check_status(s, case_name // ': validate')
      call open_file(file, case_name // '.vtkhdf', VTKHDF_STATIC)
      call htg_vtkhdf_write_static(file, g, s)
      call check_status(s, case_name // ': htg_vtkhdf_write_static')
      call file%close()
      call write_expectation(g, case_name // '.expect', case_name, ncells)
      write (manifest_unit, '(a,1x,a,1x,i0)') case_name // '.vtkhdf', case_name // '.expect', -1
      ncases = ncases + 1
      print '(a,a,a,i0,a,i0,a)', '  wrote ', case_name, ': ', ncells, ' cells exported of ', &
         g%get_number_of_cells(), ' in the grid'
   end subroutine write_static_case

   ! ---------------------------------------------------------------------------
   ! Cases
   ! ---------------------------------------------------------------------------

   !> 2D, sparse forest, implicit ids, all five array types, active scalars,
   !> and a mask on leaves and refined cells.
   subroutine case_static_2d()
      type(htg_grid_t), target :: g
      type(htg_cell_data_t), pointer :: cd
      type(htg_status_t) :: s
      integer :: nmasked

      call make_grid(g, [4_i32_k, 3_i32_k, 1_i32_k], 2_i32_k, .false.)
      call grow_forest(g, 3_i32_k, 0.6_r64_k, skip_tree=4_id_k)
      call g%initialize_local_index_node(s)
      call add_standard_fields(g, g%get_number_of_cells())
      call fill_standard_fields(g)
      call g%get_cell_data(cd)
      call cd%set_active_scalars('Pressure', s)
      call mask_some(g, nmasked)
      call check(nmasked > 0, 'static_2d: some cells are masked')
      call write_static_case(g, 'rt_static_2d')
   end subroutine case_static_2d

   !> 3D, branch factor 3, transposed root indexing, and explicit global ids
   !> assigned in reverse order, so cell data needs a real permutation.
   subroutine case_static_3d_ternary_explicit_ids()
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_breadth_first_t) :: bf
      type(htg_status_t) :: s
      integer(id_k) :: tid, total, counter
      integer(i32_k) :: ic

      ! Every root is refined into 27 children, some of which are refined again
      call make_grid(g, [3_i32_k, 2_i32_k, 2_i32_k], 3_i32_k, .true.)
      do tid = 0_id_k, g%get_max_number_of_trees() - 1_id_k
         call new_non_oriented_cursor(g, tid, c, create=.true., status=s)
         call c%subdivide_leaf(s)
         do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
            call c%to_child(ic, s)
            call refine_random(c, 2_i32_k, 0.2_r64_k)
            call c%to_parent(s)
         end do
      end do
      call check(g%get_number_of_levels() == 3_id_k, 'ternary: three levels were built')
      total = g%get_number_of_cells()
      counter = 0_id_k
      do tid = 0_id_k, g%get_max_number_of_trees() - 1_id_k
         call new_non_oriented_cursor(g, tid, c, create=.false., status=s)
         call number_reversed(c, total, counter)
      end do
      call add_standard_fields(g, total)
      call fill_standard_fields(g)
      call g%export_breadth_first(bf, s)
      call check(.not. bf%identity_order, 'ternary: explicit ids are not in export order')
      call write_static_case(g, 'rt_static_3d_ternary_explicit')
   end subroutine case_static_3d_ternary_explicit_ids

   recursive subroutine number_reversed(c, total, counter)
      type(htg_cursor_t), intent(inout) :: c
      integer(id_k), intent(in) :: total
      integer(id_k), intent(inout) :: counter
      type(htg_status_t) :: s
      integer(i32_k) :: ic
      logical :: leaf

      call c%set_global_index_from_local(total - 1_id_k - counter, s)
      call check_status(s, 'explicit global id')
      counter = counter + 1_id_k
      leaf = c%is_leaf(s)
      if (leaf) return
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         call number_reversed(c, total, counter)
         call c%to_parent(s)
      end do
   end subroutine number_reversed

   !> xz plane (where VTK's axis metadata and child numbering differ), with a
   !> material interface: normals and intercepts arrays named in the file.
   subroutine case_static_xz_interface()
      type(htg_grid_t), target :: g
      type(htg_cell_data_t), pointer :: cd
      type(htg_field_array_t), pointer :: normals, intercepts
      type(htg_status_t) :: s
      integer(id_k) :: gid, ncells

      call make_grid(g, [3_i32_k, 1_i32_k, 3_i32_k], 2_i32_k, .false.)
      call grow_forest(g, 3_i32_k, 0.6_r64_k, skip_tree=-1_id_k)
      call g%initialize_local_index_node(s)
      ncells = g%get_number_of_cells()
      call add_standard_fields(g, ncells)
      call fill_standard_fields(g)
      call g%get_cell_data(cd)
      call cd%add_array('Normals', HTG_TYPE_REAL64, 3_i32_k, ncells, s)
      call cd%add_array('Intercepts', HTG_TYPE_REAL64, 3_i32_k, ncells, s)
      call cd%get_array('Normals', normals, s)
      call cd%get_array('Intercepts', intercepts, s)
      do gid = 0_id_k, ncells - 1_id_k
         call normals%set_tuple_r64(gid, [1.0_r64_k, 0.0_r64_k, real(mod(gid, 2_id_k), r64_k)], s)
         call intercepts%set_tuple_r64(gid, &
            [0.25_r64_k * real(gid, r64_k), 0.0_r64_k, real(mod(gid, 3_id_k), r64_k)], s)
      end do
      call g%set_interface_normals_name('Normals')
      call g%set_interface_intercepts_name('Intercepts')
      call g%set_has_interface(.true.)
      call write_static_case(g, 'rt_static_xz_interface')
   end subroutine case_static_xz_interface

   !> 1D along y, branch factor 3, descending coordinates (negative sizes).
   subroutine case_static_1d_descending()
      type(htg_grid_t), target :: g
      type(htg_status_t) :: s

      call g%initialize()
      call g%set_branch_factor(3_i32_k, s)
      call g%set_dimensions(1_i32_k, 5_i32_k, 1_i32_k, s)
      call g%set_y_coordinates([8.0_r64_k, 6.0_r64_k, 3.0_r64_k, 1.0_r64_k, 0.0_r64_k], s)
      call check_status(s, '1D descending setup')
      call grow_forest(g, 3_i32_k, 0.6_r64_k, skip_tree=-1_id_k)
      call g%initialize_local_index_node(s)
      call add_standard_fields(g, g%get_number_of_cells())
      call fill_standard_fields(g)
      call write_static_case(g, 'rt_static_1d_descending')
   end subroutine case_static_1d_descending

   !> A depth limiter of 1: the file holds only two levels of a deeper grid.
   subroutine case_depth_limited()
      type(htg_grid_t), target :: g
      type(htg_status_t) :: s

      call make_grid(g, [3_i32_k, 3_i32_k, 1_i32_k], 2_i32_k, .false.)
      call grow_forest(g, 3_i32_k, 0.7_r64_k, skip_tree=-1_id_k)
      call check(g%get_number_of_levels() >= 3_id_k, 'depth_limited: grid is deeper than the limit')
      call g%initialize_local_index_node(s)
      call add_standard_fields(g, g%get_number_of_cells())
      call fill_standard_fields(g)
      call g%set_depth_limiter(1_i32_k, s)
      call write_static_case(g, 'rt_depth_limited')
   end subroutine case_depth_limited

   !> Refined, filled, then coarsened without touching the data again: the
   !> surviving cells must still carry their own values, in the file too.
   subroutine case_coarsened()
      integer(id_k), parameter :: CAPACITY = 400_id_k      ! > 1 + 4 + 16 + 64 + 256
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_status_t) :: s
      integer(id_k) :: tid, ntrees
      integer :: attempt, ncoarsened
      logical :: leaf, terminal

      call make_grid(g, [3_i32_k, 3_i32_k, 1_i32_k], 2_i32_k, .false.)
      ntrees = g%get_max_number_of_trees()
      do tid = 0_id_k, ntrees - 1_id_k
         call new_non_oriented_cursor(g, tid, c, create=.true., status=s)
         call c%set_global_index_start(tid * CAPACITY, s)    ! a fixed id range per tree
         call refine_random(c, 4_i32_k, 0.7_r64_k)
      end do
      call add_standard_fields(g, ntrees * CAPACITY)
      call fill_standard_fields(g)

      ncoarsened = 0
      do attempt = 1, 400
         tid = int(rand01() * real(ntrees, r64_k), id_k)
         call new_non_oriented_cursor(g, tid, c, create=.false., status=s)
         do
            leaf = c%is_leaf(s)
            if (leaf) exit
            terminal = c%is_terminal_node(s)
            if (terminal) then
               if (rand01() < 0.5_r64_k) then
                  call c%coarsen_node(s)
                  call check_status(s, 'coarsen_node')
                  ncoarsened = ncoarsened + 1
               end if
               exit
            end if
            call c%to_child(int(rand01() * 4.0_r64_k, i32_k), s)
         end do
      end do
      call check(ncoarsened > 20, 'coarsened: enough cells were coarsened')
      call write_static_case(g, 'rt_coarsened')
   end subroutine case_coarsened

   !> A changing mesh: three time steps of refining around a moving circle and
   !> coarsening behind it, with temporal cell data.
   subroutine case_dynamic()
      integer(id_k), parameter :: CAPACITY = 512_id_k
      character(len=*), parameter :: case_name = 'rt_dynamic'
      type(htg_grid_t), target :: g
      type(htg_cursor_t) :: c
      type(htg_cell_data_t), pointer :: cd
      type(htg_field_array_t), pointer :: distance, level
      type(htg_breadth_first_t) :: bf
      type(vtkhdf_htg_file) :: file
      type(vtkhdf_cell_data_handle) :: h_distance, h_level
      type(htg_status_t) :: s
      real(r64_k) :: centre(0:1)
      integer(id_k) :: tid, ntrees, ncells
      integer :: step
      character(len=64) :: expect_name

      call g%initialize()
      call g%set_dimensions(5_i32_k, 3_i32_k, 1_i32_k, s)
      call g%set_x_coordinates(stretched(5_i32_k, 0.0_r64_k, 1.0_r64_k, 0.0_r64_k), s)
      call g%set_y_coordinates(stretched(3_i32_k, 0.0_r64_k, 1.0_r64_k, 0.0_r64_k), s)
      ntrees = g%get_max_number_of_trees()
      do tid = 0_id_k, ntrees - 1_id_k
         call new_non_oriented_cursor(g, tid, c, create=.true., status=s)
         call c%set_global_index_start(tid * CAPACITY, s)
      end do
      call g%get_cell_data(cd)
      call cd%add_array('Distance', HTG_TYPE_REAL64, 1_i32_k, ntrees * CAPACITY, s)
      call cd%add_array('Level', HTG_TYPE_INT32, 1_i32_k, ntrees * CAPACITY, s)
      call cd%get_array('Distance', distance, s)
      call cd%get_array('Level', level, s)
      call cd%set_active_scalars('Distance', s)      ! registered below with the 'Scalars' role
      call check_status(s, 'dynamic setup')

      call open_file(file, case_name // '.vtkhdf', VTKHDF_DYNAMIC_MESH)
      call htg_vtkhdf_write_mesh(file, g, s)
      call htg_vtkhdf_write_geometry(file, g, s)
      call htg_vtkhdf_register_temporal_cell_data(file, distance, h_distance, s, attribute='Scalars')
      call htg_vtkhdf_register_temporal_cell_data(file, level, h_level, s)
      call check_status(s, 'dynamic file setup')

      do step = 0, 2
         centre = [0.9_r64_k + 1.1_r64_k * real(step, r64_k), 1.0_r64_k]
         do tid = 0_id_k, ntrees - 1_id_k
            call new_non_oriented_cursor(g, tid, c, create=.false., status=s)
            call coarsen_away_from_circle(c, centre)
         end do
         do tid = 0_id_k, ntrees - 1_id_k
            call new_non_oriented_cursor(g, tid, c, create=.false., status=s)
            call refine_at_circle(c, centre)
         end do
         do tid = 0_id_k, ntrees - 1_id_k
            call new_non_oriented_cursor(g, tid, c, create=.false., status=s)
            call sample(c, centre, distance, level)
         end do
         call g%validate(s)
         call check_status(s, 'dynamic validate')

         call file%start_time_step(real(step, r64_k))
         call htg_vtkhdf_write_topology(file, g, bf, s)
         call check_status(s, 'dynamic topology')
         call htg_vtkhdf_write_temporal_cell_data(file, h_distance, bf, distance, s)
         call htg_vtkhdf_write_temporal_cell_data(file, h_level, bf, level, s)
         call check_status(s, 'dynamic cell data')
         call file%finalize_time_step()

         write (expect_name, '(a,a,i0,a)') case_name, '.step', step, '.expect'
         call write_expectation(g, trim(expect_name), case_name, ncells)
         write (manifest_unit, '(a,1x,a,1x,i0)') case_name // '.vtkhdf', trim(expect_name), step
         ncases = ncases + 1
         print '(a,a,a,i0,a,i0,a)', '  wrote ', case_name, ' step ', step, ': ', ncells, ' cells'
      end do
      call file%close()
   end subroutine case_dynamic

   logical function circle_crosses(c, centre) result(crosses)
      type(htg_cursor_t), intent(in) :: c
      real(r64_k), intent(in) :: centre(0:1)
      real(r64_k), parameter :: RADIUS = 0.6_r64_k
      real(r64_k) :: b(0:5), nearest, farthest, lo, hi, d
      type(htg_status_t) :: s
      integer :: a

      call c%get_bounds(b, s)
      nearest = 0.0_r64_k
      farthest = 0.0_r64_k
      do a = 0, 1
         lo = b(2*a) - centre(a)
         hi = b(2*a + 1) - centre(a)
         d = max(lo, -hi, 0.0_r64_k)
         nearest = nearest + d*d
         farthest = farthest + max(lo*lo, hi*hi)
      end do
      crosses = (nearest < RADIUS*RADIUS) .and. (farthest > RADIUS*RADIUS)
   end function circle_crosses

   recursive subroutine refine_at_circle(c, centre)
      type(htg_cursor_t), intent(inout) :: c
      real(r64_k), intent(in) :: centre(0:1)
      type(htg_status_t) :: s
      integer(i32_k) :: ic
      logical :: leaf, needed

      leaf = c%is_leaf(s)
      if (leaf) then
         if (c%get_level() >= 4_i32_k) return
         needed = circle_crosses(c, centre)
         if (.not. needed) return
         call c%subdivide_leaf(s)
      end if
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         call refine_at_circle(c, centre)
         call c%to_parent(s)
      end do
   end subroutine refine_at_circle

   recursive subroutine coarsen_away_from_circle(c, centre)
      type(htg_cursor_t), intent(inout) :: c
      real(r64_k), intent(in) :: centre(0:1)
      type(htg_status_t) :: s
      integer(i32_k) :: ic
      logical :: leaf, terminal, needed

      leaf = c%is_leaf(s)
      if (leaf) return
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         call coarsen_away_from_circle(c, centre)
         call c%to_parent(s)
      end do
      terminal = c%is_terminal_node(s)
      needed = circle_crosses(c, centre)
      if (terminal .and. .not. needed) then
         call c%coarsen_node(s)
         call check_status(s, 'dynamic coarsen_node')
      end if
   end subroutine coarsen_away_from_circle

   recursive subroutine sample(c, centre, distance, level)
      type(htg_cursor_t), intent(inout) :: c
      real(r64_k), intent(in) :: centre(0:1)
      type(htg_field_array_t), intent(inout) :: distance, level
      type(htg_status_t) :: s
      real(r64_k) :: p(0:2)
      integer(id_k) :: gid
      integer(i32_k) :: ic
      logical :: leaf

      call c%get_point(p, s)
      gid = c%get_global_node_index(s)
      call distance%set_scalar_r64(gid, norm2(p(0:1) - centre) - 0.6_r64_k, s)
      call level%set_scalar_i32(gid, c%get_level(), s)
      leaf = c%is_leaf(s)
      if (leaf) return
      do ic = 0_i32_k, c%get_number_of_children() - 1_i32_k
         call c%to_child(ic, s)
         call sample(c, centre, distance, level)
         call c%to_parent(s)
      end do
   end subroutine sample

end program test_vtkhdf_roundtrip_write
