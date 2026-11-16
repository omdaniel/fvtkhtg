# `fvtkhtg`

**Pure Fortran 2018 `vtkHyperTreeGrid` Construction, AMR Bookkeeping, and Structural Validation Library**

`fvtkhtg` provides a modern, self-contained **100% pure Fortran 2018** implementation of VTK's `vtkHyperTreeGrid`, `vtkHyperTree`, `vtkHyperTreeGridNonOrientedCursor`, `vtkBitArray`, and cell-centered `CellData` structures—engineered for scientific AMR solvers and downstream visualization I/O libraries with **zero C, C++, `iso_c_binding`, VTK, or external runtime dependencies**.

All topological conventions, squeezed-plane orientations, MSB-first bit layouts, breadth-first descriptors, and `PureMask` semantics are differentially verified and LLDB-lockstep-tested against **C++ VTK `v9.7.2`**.

`fvtkhtg` owns **creation**: building, refining, indexing and validating a HyperTreeGrid in-process. It does not write files. Serialization belongs to [`fvtkhdf`](https://github.com/omdaniel/fvtkhdf/tree/feature/hypertreegrid), whose `vtkhdf_htg_file` writes VTKHDF HyperTreeGrid files (ParaView 6.2+ / VTK 9.7+) with MPI and time-step support. `grid%export_breadth_first` produces exactly the arrays that writer takes, and an optional bridge library calls it for you; neither core library depends on the other ([Writing VTKHDF Files](#writing-vtkhdf-files-with-fvtkhdf)).

---

## Key Features

- **Pure Fortran 2018 (`libfvtkhtg.a`)**: Zero external libraries; compiles cleanly under `flang` and `gfortran` with `-std=f2018 -Werror`.
- **0D/1D/2D/3D Rectilinear Forests**: Binary (`branch_factor = 2`) and ternary (`branch_factor = 3`) HyperTrees ($k = f^d \in \{2, 3, 4, 8, 9, 27\}$), monotonic ascending/descending coordinates, KJI and transposed IJK root indexing, and exact VTK `v9.7+` squeezed orientation axes.
- **Compile-Time Cursor Safety**: Stack-based `htg_cursor_t` (`new_non_oriented_cursor`) enforces `target` on `htg_grid_t` at compile time (Fortran 2018 C1541) while internal `move_alloc` boxing preserves tree and field addresses across dynamic growth.
- **Implicit & Explicit Global Indexing**: Forest-wide contiguous offsets (`initialize_local_index_node`), per-tree implicit offsets (`set_global_index_start`), and arbitrary per-vertex mappings (`set_global_index_from_local`).
- **Dynamic `CellData`, Masking & Material Interfaces**: Auto-growing 0-based multi-component field arrays (`REAL64`, `REAL32`, `INT64`, `INT32`, `BIT`), cell `Mask`, 3-component `Normals` & `Intercepts` interface metadata, and lazy bottom-up `PureMask`.
- **Breadth-First Descriptors & Zero-Copy Views**: Bit-stream tree construction/export (`build_from_breadth_first_order_descriptor`, `compute_breadth_first_order_descriptor`) and zero-copy `0:`-bounded pointer views into every internal topological, coordinate, MSB-first byte, and AOS `CellData` buffer.
- **Cell Geometry, Face Neighbours & Coarsening**: On the same cursor, `get_bounds` / `get_point` give a cell's box and centre (bit-identical to VTK's geometry cursor), `get_face_neighbor` finds the same-level or coarser cell across a face inside a tree or across root cells, and `coarsen_node` merges children back into their parent while `CellData` follows the cells — enough to adapt a mesh in place with 2:1 balance.
- **Whole-Grid Structural Validation**: Non-recursive `call grid%validate(status)` auditing tree topology, coordinates, disjoint global IDs, `CellData` tuple counts, interfaces, and bit-array invariants.
- **Breadth-First Export for Serializers**: `grid%export_breadth_first(bf)` and `bf%gather(field, values)` yield the per-cell breadth-first topology, mask and permuted `CellData` that the `fvtkhdf` writer takes unchanged, matching what VTK's own `vtkHDFWriter` stores.

---

## Quick Start

```fortran
program quickstart
   use fvtkhtg, only: &
      htg_grid_t, htg_cursor_t, htg_cell_data_t, htg_field_array_t, htg_status_t, &
      new_non_oriented_cursor, HTG_TYPE_REAL64, id_k, i32_k, r64_k
   implicit none

   type(htg_grid_t), target         :: grid
   type(htg_cursor_t)               :: cursor
   type(htg_cell_data_t), pointer   :: cell_data => null()
   type(htg_field_array_t), pointer :: pressure  => null()
   type(htg_status_t)               :: status
   integer(id_k)                    :: gid

   ! 1. Configure a 2x2 2D rectilinear quadtree grid in the XY plane (dz = 1)
   call grid%initialize()
   call grid%set_branch_factor(2_i32_k, status)
   call grid%set_dimensions(3_i32_k, 3_i32_k, 1_i32_k, status)
   call grid%set_x_coordinates([0.0_r64_k, 0.5_r64_k, 1.0_r64_k], status)
   call grid%set_y_coordinates([0.0_r64_k, 0.5_r64_k, 1.0_r64_k], status)

   ! 2. Instantiate root tree 0 and refine selectively with a stack cursor
   call new_non_oriented_cursor(grid, 0_id_k, cursor, create=.true., status=status)
   call cursor%subdivide_leaf(status)
   call cursor%to_child(3_i32_k, status)
   call cursor%subdivide_leaf(status)

   ! 3. Assign contiguous 0-based implicit global indices and attach CellData
   call grid%initialize_local_index_node(status)
   call grid%get_cell_data(cell_data)
   call cell_data%add_array('Pressure', HTG_TYPE_REAL64, &
      num_components=1_i32_k, num_tuples=grid%get_number_of_cells(), status=status)
   call cell_data%get_array('Pressure', pressure, status)

   gid = cursor%get_global_node_index(status)
   call pressure%set_scalar_r64(gid, 101325.0_r64_k, status)

   ! 4. Validate all grid, tree, global-index, and CellData invariants
   call grid%validate(status)
   if (.not. status%ok()) error stop trim(status%message)
end program quickstart
```

---

## Writing VTKHDF Files with `fvtkhdf`

| Concern | Owner |
| :--- | :--- |
| Trees, refinement, cursors, global IDs, `CellData`, mask, validation | `fvtkhtg` |
| Breadth-first order, pruning, the `CellData` permutation | `fvtkhtg` (`export_breadth_first`, `gather`) |
| Bit packing, HDF5 layout, MPI parts, time steps | `fvtkhdf` |

```fortran
type(htg_breadth_first_t)  :: bf
real(r64_k), allocatable   :: p(:)

call grid%export_breadth_first(bf, status)
call bf%gather(pressure, p, status)

call file%write_mesh(int(grid%get_dimensions()), int(grid%get_branch_factor()))
call file%write_mesh_geometry(grid%x_coordinates, grid%y_coordinates, grid%z_coordinates)
call file%write_mesh_topology(bf%tree_ids, bf%xtree, bf%refined, bf%mask)
call file%write_cell_data('Pressure', p)
```

**Index convention at the seam.** Everything inside `fvtkhtg` is 0-based, matching VTK. `htg_breadth_first_t` is the one exception: `tree_ids` and `xtree` are 1-based, as `fvtkhdf` expects, while `global_ids` keeps `fvtkhtg`'s 0-based global IDs and is used only for gathering.

With `-DFVTKHTG_BUILD_VTKHDF_BRIDGE=ON`, the separate `fvtkhtg::fvtkhtg_vtkhdf` library wraps these calls (`call htg_vtkhdf_write_static(file, grid, status)` writes a whole grid). The core `fvtkhtg` library never links `fvtkhdf`. The same option builds a round-trip test (`ctest -L vtkhdf`) that writes seven kinds of grid through the bridge, including masks, a material interface, a coarsened grid and a changing mesh, and reads every cell back with VTK's own `vtkHDFReader`, the reader ParaView uses. See [§4.9 of the user guide](docs/user_guide.md#49-writing-vtkhdf-files-breadth-first-export--the-fvtkhdf-seam) and [`examples/05_breadth_first_export.f90`](examples/05_breadth_first_export.f90).

---

## Building & Running Tests

### Requirements
- **CMake** `>= 3.24` and **Ninja**
- **Fortran 2018 compiler**: LLVM `flang` (`23.1+` tested) or GNU `gfortran` (`16.2+` tested)

### Pure-Fortran Build (Default — Zero C++/VTK Dependencies)

```bash
# Release preset
cmake --preset pure-fortran-release
cmake --build --preset pure-fortran-release
ctest --preset pure-fortran-release

# Debug preset
cmake --preset pure-fortran-debug
cmake --build --preset pure-fortran-debug
ctest --preset pure-fortran-debug
```

To configure manually with `gfortran` runtime checks (`-fcheck=all`) and run the example programs:

```bash
cmake -S . -B build/gfortran-debug -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DFVTKHTG_BUILD_TESTING=ON \
  -DFVTKHTG_BUILD_EXAMPLES=ON \
  -DFVTKHTG_WARNINGS_AS_ERRORS=ON
cmake --build build/gfortran-debug
ctest --test-dir build/gfortran-debug --output-on-failure

# Run only the example programs with verbose output
ctest --test-dir build/gfortran-debug -L examples -V
```

### Optional C++ VTK `v9.7.2` Differential & LLDB Lockstep Oracle Suite

For library contributors verifying bit-exact parity against C++ VTK `v9.7.2` (`VTK::CommonCore`, `VTK::CommonDataModel` installed under `.cache/vtk-oracle/install`) and LLDB DWARF lockstep scripts managed via [`uv`](https://docs.astral.sh/uv/) (`Python >=3.14,<3.15` per `pyproject.toml`):

```bash
# Sync the pinned Python 3.14 environment for LLDB verification tools
uv sync

# Run Flang + C++ VTK v9.7.2 differential & LLDB lockstep oracle suite
cmake --preset flang-debug-oracle
cmake --build --preset flang-debug-oracle
ctest --preset flang-debug-oracle

# Run gfortran + C++ VTK v9.7.2 differential & LLDB lockstep oracle suite
cmake --preset gfortran-debug-oracle
cmake --build --preset gfortran-debug-oracle
ctest --preset gfortran-debug-oracle
```

---

## Repository Layout

```text
fvtkhtg/
├── CMakeLists.txt                       # Pure-Fortran build + optional examples/tests/oracle
├── CMakePresets.json                    # Pure-Fortran & VTK v9.7.2 oracle presets
├── pyproject.toml                       # uv Python >=3.14,<3.15 config for LLDB verification tools
├── cmake/
│   ├── CompilerWarnings.cmake           # Strict F2018 flags for flang & gfortran (per target)
│   ├── fvtkhtgConfig.cmake.in           # find_package(fvtkhtg [COMPONENTS vtkhdf]) package config
│   └── VtkOracle.cmake                  # Optional VTK v9.7.2 + uv locator
├── src/
│   ├── fvtkhtg.f90                      # Public umbrella module (`use fvtkhtg`)
│   ├── fvtkhtg_kinds_m.f90              # Kinds, status/error codes, htg_status_t
│   ├── fvtkhtg_bit_array_m.f90          # MSB-first packed bit vector (htg_bit_array_t)
│   ├── fvtkhtg_cell_data_m.f90          # Dynamic 0-based multi-component CellData arrays
│   ├── fvtkhtg_breadth_first_m.f90      # Breadth-first export type & CellData gather (htg_breadth_first_t)
│   ├── fvtkhtg_tree_m.f90               # Compact k-ary HyperTree & BF descriptors (htg_tree_t)
│   ├── fvtkhtg_grid_m.f90               # Rectilinear forest, Mask, PureMask & interfaces (htg_grid_t)
│   ├── fvtkhtg_cursor_m.f90             # Build-time stack cursor (htg_cursor_t), coarsening
│   ├── fvtkhtg_cursor_geometry_s.f90    # Cell geometry & face-neighbour lookup (submodule of fvtkhtg_cursor_m)
│   ├── fvtkhtg_validate_m.f90           # Whole-grid & per-tree invariant validation
│   └── fvtkhtg_grid_export_s.f90        # grid%export_breadth_first (submodule of fvtkhtg_grid_m)
├── bridge/
│   └── vtkhdf/fvtkhtg_vtkhdf.f90        # Optional bridge to the fvtkhdf writer (FVTKHTG_BUILD_VTKHDF_BRIDGE)
├── examples/
│   ├── 01_basic_quadtree.f90            # 2D quadtree AMR, cursor navigation & validation
│   ├── 02_octree_amr_celldata.f90       # 3D sparse octree, KJI/IJK roots, implicit/explicit IDs & CellData
│   ├── 03_mask_puremask_interface.f90   # Cell masking, material interfaces & lazy PureMask caching
│   ├── 04_bf_descriptors_and_views.f90  # Breadth-first bit descriptors & zero-copy serializer views
│   ├── 05_breadth_first_export.f90      # Breadth-first export & gather for the fvtkhdf writer
│   └── 06_adaptive_refinement.f90       # Geometry-driven refinement, 2:1 balancing via neighbours, coarsening
├── docs/
│   ├── user_guide.md                    # Comprehensive end-user & downstream integration guide
│   └── vtk_oracle_audit.md              # VTK v9.7.2 C++ source & LLDB ground-truth audit
├── test/
│   ├── unit/                            # Pure-Fortran unit test suites
│   ├── vtkhdf/                          # Round trip through fvtkhdf and VTK's VTKHDF reader (with the bridge)
│   └── oracle/                          # C++17 VTK v9.7.2 differential & lockstep tests (test-only)
└── tools/
    └── lldb/                            # Python 3.14 LLDB batch lockstep & DWARF inspection scripts
```

---

## Documentation

- **[End-User & Integration Guide (`docs/user_guide.md`)](docs/user_guide.md)**: Complete walkthrough of `fvtkhtg` conventions (`0:` indexing, `target` cursor safety, `htg_status_t`), step-by-step tutorials, zero-copy view extraction for external serializers, and C++ `vtkHyperTreeGrid` $\leftrightarrow$ Fortran `fvtkhtg` API reference tables.
- **[VTK `v9.7.2` Ground-Truth Audit (`docs/vtk_oracle_audit.md`)](docs/vtk_oracle_audit.md)**: Detailed C++ source and `/usr/bin/lldb` memory/behavior audit of VTK `v9.7.2` edge cases.
