# `fvtkhtg` End-User & Integration Guide

`fvtkhtg` is a **100% pure Fortran 2018** library for programmatic construction, adaptive mesh refinement (AMR) cursor traversal, global cell indexing, cell-centered field storage (`CellData`), masking (`Mask` & `PureMask`), material interface metadata, breadth-first bit descriptor import/export, structural validation, and zero-copy buffer extraction for **VTK HyperTreeGrids (`vtkHyperTreeGrid`)**.

Every topological convention, array layout, bit-packing rule, and edge case in `fvtkhtg` is verified in lockstep against **VTK `v9.7+`** (pinned and differentially tested against C++ VTK `v9.7.2`; see [`docs/vtk_oracle_audit.md`](vtk_oracle_audit.md)), while the library itself (`libfvtkhtg.a`) has **zero C, C++, `iso_c_binding`, VTK, or third-party dependencies**.

---

## Table of Contents

1. [Overview & Design Scope](#1-overview--design-scope)
2. [Building & Integrating `fvtkhtg`](#2-building--integrating-fvtkhtg)
3. [Core Conventions & Safety Idioms](#3-core-conventions--safety-idioms)
   - [3.1 Zero-Based Indexing Everywhere (`0:`)](#31-zero-based-indexing-everywhere-0)
   - [3.2 Compile-Time Cursor & View Lifetime Safety (`target`)](#32-compile-time-cursor--view-lifetime-safety-target)
   - [3.3 Error Handling (`htg_status_t` vs. Fail-Fast)](#33-error-handling-htg_status_t-vs-fail-fast)
4. [Step-by-Step Tutorials & Annotated Patterns](#4-step-by-step-tutorials--annotated-patterns)
   - [4.1 Grid Setup, Squeezed Orientations & Rectilinear Coordinates](#41-grid-setup-squeezed-orientations--rectilinear-coordinates)
   - [4.2 Cursor-Driven AMR Refinement](#42-cursor-driven-amr-refinement)
   - [4.3 Implicit vs. Explicit Global Node Indexing](#43-implicit-vs-explicit-global-node-indexing)
   - [4.4 Dynamic `CellData` & Material Interfaces](#44-dynamic-celldata--material-interfaces)
   - [4.5 Cell Masking (`Mask`) & Lazy Bottom-Up `PureMask`](#45-cell-masking-mask--lazy-bottom-up-puremask)
   - [4.6 Breadth-First Order Bit Descriptors](#46-breadth-first-order-bit-descriptors)
   - [4.7 Structural Invariant Validation (`grid%validate`)](#47-structural-invariant-validation-gridvalidate)
   - [4.8 Zero-Copy Buffer Extraction for External Serializers](#48-zero-copy-buffer-extraction-for-external-serializers)
   - [4.9 Writing VTKHDF Files: Breadth-First Export & the `fvtkhdf` Seam](#49-writing-vtkhdf-files-breadth-first-export--the-fvtkhdf-seam)
   - [4.10 Cell Geometry](#410-cell-geometry)
   - [4.11 Face Neighbours](#411-face-neighbours)
   - [4.12 Coarsening](#412-coarsening)
5. [C++ `vtkHyperTreeGrid` to Fortran `fvtkhtg` API Reference](#5-c-vtkhypertreegrid-to-fortran-fvtkhtg-api-reference)

---

## 1. Overview & Design Scope

### What `fvtkhtg` Provides
- **Rectilinear Forest Geometry (0D/1D/2D/3D)**: Point dimensions, automatic squeezed-axis orientation detection (including VTK's cyclic right-handed 2D XZ `[2, 0]` axes), strictly monotonic ascending or descending rectilinear coordinates, and both default KJI (`i`-fastest) and transposed IJK (`k`-fastest) level-zero root indexing.
- **Sparse Multi-Tree Forest (`htg_grid_t`) & Compact Trees (`htg_tree_t`)**: Binary (`branch_factor = 2`) and ternary (`branch_factor = 3`) trees in 1D, 2D, and 3D ($k = f^d \in \{2, 3, 4, 8, 9, 27\}$ children per coarse vertex), backed by VTK `v9.7+`'s compact `ParentToElderChild` table with lazy tail sizing and $O(\log N_{\text{trees}})$ sparse root storage.
- **Build-Time Non-Oriented Stack Cursor (`htg_cursor_t`)**: Root-to-leaf traversal (`to_root`, `to_child`, `to_parent`), incremental AMR subdivision (`subdivide_leaf`), global ID assignment, and per-cell masking.
- **Cell Geometry**: Origin, size, bounds and centre of the cell under a cursor (`get_origin`, `get_size`, `get_bounds`, `get_point`), bit-identical to `vtkHyperTreeGridNonOrientedGeometryCursor`, so a solver can decide where to refine and sample fields without tracking coordinates itself.
- **Face Neighbours**: `cursor%get_face_neighbor(axis, direction, neighbor)` returns a cursor on the same-level or coarser cell across a face, inside a tree or across root cells — what a solver needs for stencils and 2:1 balancing.
- **Coarsening**: `cursor%coarsen_node` merges the children of a cell back into it (the inverse of `subdivide_leaf`, which VTK does not offer), keeping local IDs contiguous and `CellData` and mask bits attached to their cells.
- **Global Node Indexing**: Contiguous forest-wide implicit offsets (`initialize_local_index_node`), per-tree implicit offsets (`set_global_index_start`), and arbitrary per-vertex explicit mappings (`set_global_index_from_local`).
- **Dynamic `CellData` (`htg_cell_data_t` / `htg_field_array_t`)**: Named 0-based multi-component cell arrays (`REAL64`, `REAL32`, `INT64`, `INT32`, and MSB-first packed `BIT`) that grow automatically as leaves subdivide.
- **Masking, Material Interfaces & `PureMask`**: Cell blanking (`Mask`), 3-component `Normals` and `Intercepts` interface metadata, and lazy bottom-up `PureMask` evaluation with automatic cache invalidation.
- **Breadth-First Order Descriptors**: Building trees from packed bit slices (`build_from_breadth_first_order_descriptor`) and exporting forest descriptors (`compute_breadth_first_order_descriptor`) with optional `input_mask` subtree pruning and `depth_limiter`.
- **Whole-Grid Structural Validation (`grid%validate`)**: Non-recursive $O(N \log N)$ / $O(N)$ verification of tree topology, depth histograms, coordinate monotonicity, disjoint global IDs, `CellData` tuple sufficiency, interface arrays, and bit-array slack invariants.
- **Zero-Copy `0:`-Bounded Pointer Views**: Direct pointer access to every internal array and MSB-first byte buffer for external I/O writers and MPI communication layers.

### Explicit Non-Goals
1. **Post-Processing & VTK's Other Cursor Classes**: `fvtkhtg` covers what a solver needs to build, adapt and hand off a grid, not what a viewer needs to process one. It intentionally omits post-processing visualization filters (contour, clip, slice, dual-grid generation) and does not reproduce VTK's separate cursor classes (oriented, "unlimited", Moore supercursors with edge and corner neighbours). The two things a solver needs from those classes are provided on the one cursor type instead: cell geometry (Section [4.10](#410-cell-geometry)) and face neighbours (Section [4.11](#411-face-neighbours)).
2. **Built-In File Serialization**: `fvtkhtg` does not write `.htg` XML or VTKHDF/HDF5 files, keeping the core library free of XML, zlib/LZ4, HDF5, or MPI dependencies. Serialization belongs to [`fvtkhdf`](https://github.com/omdaniel/fvtkhdf/tree/feature/hypertreegrid), which writes VTKHDF HyperTreeGrid files from plain breadth-first arrays. `fvtkhtg` produces exactly those arrays with `grid%export_breadth_first` (Section [4.9](#49-writing-vtkhdf-files-breadth-first-export--the-fvtkhdf-seam)); an optional bridge library joins the two without either core library depending on the other.

### Module Architecture
Applications only need a single `use` statement:

```fortran
use fvtkhtg
```

`fvtkhtg` re-exports all public types, constants, and procedures from 8 internal implementation modules:

| Module | Primary Responsibility |
| :--- | :--- |
| `fvtkhtg_kinds_m` | Integer/real kind parameters (`id_k`, `i32_k`, `r64_k`, ...), status/type/index constants, `htg_status_t` |
| `fvtkhtg_bit_array_m` | `htg_bit_array_t`: MSB-first packed bit vector matching `vtkBitArray` |
| `fvtkhtg_cell_data_m` | `htg_field_array_t`, `htg_cell_data_t`: dynamic multi-component cell-centered arrays |
| `fvtkhtg_breadth_first_m` | `htg_breadth_first_t`: breadth-first export for serializers (`grid%export_breadth_first`, `bf%gather`) |
| `fvtkhtg_tree_m` | `htg_tree_t`: single compact $k$-ary HyperTree, global indexing, breadth-first descriptors |
| `fvtkhtg_grid_m` | `htg_grid_t`: rectilinear geometry, sparse forest, `Mask`, `PureMask`, interface metadata |
| `fvtkhtg_cursor_m` | `htg_cursor_t`, `new_non_oriented_cursor`: stack-based AMR build/traversal cursor, coarsening; cell geometry and face neighbours (submodule `fvtkhtg_cursor_geometry_s`) |
| `fvtkhtg_validate_m` | `htg_validate_grid`, `htg_validate_tree` (also bound as `grid%validate`) |

---

## 2. Building & Integrating `fvtkhtg`

### Compiler & Tooling Requirements
- **CMake**: `>= 3.24` (Ninja `>= 1.10` recommended for Fortran 2018 module dependency scanning).
- **Fortran 2018 Compiler** (zero C/C++ compiler required for `fvtkhtg`, unit tests, or examples):
  - **LLVM Flang** (`flang` / `flang-new`, tested with LLVM Flang `23.1.5`)
  - **GNU Fortran** (`gfortran`, tested with GCC `16.2.0`)

### CMake Configuration Options

| CMake Option | Default | Description |
| :--- | :---: | :--- |
| `FVTKHTG_BUILD_TESTING` | `ON` when top-level | Build pure-Fortran unit tests (`test_foundations`, `test_tree`, `test_grid_cursor_validate`, `test_export`, `test_geometry_neighbors`, `test_coarsen`) |
| `FVTKHTG_BUILD_EXAMPLES` | `ON` when top-level | Build runnable Fortran 2018 example programs (`<build>/examples/01..06`) and register them with CTest when testing is enabled |
| `FVTKHTG_WARNINGS_AS_ERRORS` | `OFF` | Treat compiler warnings as errors (`-Werror`) |
| `FVTKHTG_BUILD_ORACLE_TESTS` | `OFF` | Build optional C++17 differential & LLDB lockstep tests against VTK `v9.7.2` |
| `FVTKHTG_BUILD_VTKHDF_BRIDGE` | `OFF` | Build the optional `fvtkhtg_vtkhdf` library that writes grids with `fvtkhdf` (needs an installed `fVTKHDF` >= 0.7 with the HyperTreeGrid writer) |
| `FVTKHTG_INSTALL` | `ON` when top-level | Generate install rules and a `fvtkhtgConfig.cmake` package |

"When top-level" means the default is `ON` when `fvtkhtg` is built on its own and `OFF` when it is pulled into another project with `add_subdirectory` or `FetchContent`. Compiler flags (`-std=f2018`, warnings, `-fcheck=all`) are applied to `fvtkhtg`'s own targets only and never leak into a parent project.

### Building & Testing with CMake Presets

Pure-Fortran Release or Debug builds (zero C++/VTK dependencies):

```bash
# Release build + run unit tests and examples
cmake --preset pure-fortran-release
cmake --build --preset pure-fortran-release
ctest --preset pure-fortran-release

# Debug build + run unit tests and examples
cmake --preset pure-fortran-debug
cmake --build --preset pure-fortran-debug
ctest --preset pure-fortran-debug
```

Manual configuration with `flang` or `gfortran` (e.g., enabling full runtime bounds checking in `gfortran`):

```bash
# GNU gfortran Debug build (-std=f2018 -Wall -Wextra -fcheck=all -g -O0)
cmake -S . -B build/gfortran-debug -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DFVTKHTG_BUILD_TESTING=ON \
  -DFVTKHTG_BUILD_EXAMPLES=ON \
  -DFVTKHTG_WARNINGS_AS_ERRORS=ON
cmake --build build/gfortran-debug
ctest --test-dir build/gfortran-debug --output-on-failure

# Run only the runnable example programs via CTest label
ctest --test-dir build/gfortran-debug -L examples -V
```

### Integrating into Your CMake Project

#### Option A: Vendored Subdirectory (`add_subdirectory`)

```cmake
add_subdirectory(external/fvtkhtg)   # tests and examples are off by default here

add_executable(my_amr_solver src/main.f90)
target_link_libraries(my_amr_solver PRIVATE fvtkhtg::fvtkhtg)
```

#### Option B: CMake `FetchContent`

```cmake
include(FetchContent)

FetchContent_Declare(
  fvtkhtg
  GIT_REPOSITORY https://example.org/fvtkhtg.git
  GIT_TAG        v0.1.0
)
FetchContent_MakeAvailable(fvtkhtg)

target_link_libraries(my_amr_solver PRIVATE fvtkhtg::fvtkhtg)
```

#### Option C: Installed Package (`find_package`)

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/fvtkhtg
cmake --build build
cmake --install build
```

```cmake
find_package(fvtkhtg 0.1 CONFIG REQUIRED)                    # fvtkhtg::fvtkhtg
find_package(fvtkhtg 0.1 CONFIG REQUIRED COMPONENTS vtkhdf)  # also fvtkhtg::fvtkhtg_vtkhdf
target_link_libraries(my_amr_solver PRIVATE fvtkhtg::fvtkhtg)
```

The `vtkhdf` component exists only if `fvtkhtg` was built with `-DFVTKHTG_BUILD_VTKHDF_BRIDGE=ON`; requesting it also finds `fVTKHDF`. Fortran `.mod` files only work with the compiler that produced them, so the package warns when a consumer uses a different compiler.

`fvtkhtg` exports its compiled `.mod` directory via `PUBLIC` `target_include_directories`, so linking `fvtkhtg::fvtkhtg` automatically makes `use fvtkhtg` available to your Fortran targets.

---

## 3. Core Conventions & Safety Idioms

### 3.1 Zero-Based Indexing Everywhere (`0:`)

To maintain **exact 1-to-1 numerical parity with VTK C++** (`vtkHyperTreeGrid`), every index, coordinate, bit offset, and array view in `fvtkhtg` is **0-based**. The one deliberate exception is the serializer seam, `htg_breadth_first_t` (Section [4.9](#49-writing-vtkhdf-files-breadth-first-export--the-fvtkhdf-seam)), whose arrays follow `fvtkhdf`'s 1-based Fortran conventions:

- **Root tree IDs (`tree_id`)**: `0_id_k .. grid%get_max_number_of_trees() - 1_id_k`
- **Level-zero root coordinates `(i, j, k)`**: `0 .. cell_dims(a) - 1`
- **Active tree slot positions (`pos`)**: `0_id_k .. grid%get_number_of_non_empty_trees() - 1_id_k`
- **Local vertex IDs (`local_id`)**: `0_id_k .. tree%get_number_of_vertices() - 1_id_k` (root is always `0_id_k`)
- **Global node IDs (`global_id`)**: `0_id_k .. grid%get_max_global_index()`
- **Child indices (`ichild`)**: `0_i32_k .. number_of_children - 1_i32_k`
- **Depth levels (`level` / `depth`)**: `0 .. number_of_levels - 1` (root is at depth `0`)
- **Bit and byte indices**: `0_id_k .. num_bits - 1_id_k` and `0_id_k .. num_bytes - 1_id_k`
- **Field components & tuples**: `comp_id` in `0_i32_k .. num_components - 1_i32_k`, `tuple_id` in `0_id_k .. num_tuples - 1_id_k`, and `get_view_*` pointers bounded `(0:nc-1, 0:nt-1)`

> [!IMPORTANT]
> **Fortran Assumed-Shape Dummy Lower-Bound Rule**:
> In Fortran, an assumed-shape dummy argument declared as `arr(:)` **implicitly resets its lower bound to `1`** inside the procedure—even when the caller passes a `0:`-indexed array or `0:` pointer view! Always write explicit `0:` lower bounds on dummy arguments (`arr(0:)` or `arr(0:, 0:)`) that receive `fvtkhtg` views or tuples:
>
> ```fortran
> ! CORRECT: preserves 0-based indexing (coords(0) .. coords(n-1))
> subroutine write_axis_coords(coords)
>    real(r64_k), intent(in) :: coords(0:)
>    ! lbound(coords, 1) == 0
> end subroutine write_axis_coords
>
> ! BUGGY IF INDEXED FROM 0: Fortran resets lbound(coords, 1) to 1!
> subroutine broken_axis_coords(coords)
>    real(r64_k), intent(in) :: coords(:)
>    ! coords(0) is OUT OF BOUNDS here!
> end subroutine broken_axis_coords
> ```
>
> Likewise, when remapping or slicing a `0:` pointer into another pointer in user code, specify explicit `0:` lower bounds (`sub_ptr(0:) => x_view(0:n-1)`). All `fvtkhtg` `get_*_view` procedures already set `0:` lower bounds directly on their `pointer, intent(out)` dummy arguments.

### 3.2 Compile-Time Cursor & View Lifetime Safety (`target`)

`htg_cursor_t` stores internal Fortran pointers (`cursor%grid` and `cursor%tree`) to traverse and refine trees without copying grid state. To prevent dangling-pointer bugs, `new_non_oriented_cursor` and `cursor%initialize` declare their `grid` dummy argument as:

```fortran
type(htg_grid_t), pointer, intent(in) :: grid
```

Under **Fortran 2018 Rule C1541**, you pass a regular `htg_grid_t` variable directly (no manual pointer temporary needed), **provided that variable carries the `target` (or `pointer`) attribute**. Forgetting `target` produces an immediate compile-time error instead of a silent runtime memory bug:

```fortran
type(htg_grid_t), target :: grid   ! <-- `target` is enforced at compile time
type(htg_cursor_t)       :: cursor
type(htg_status_t)       :: status

call grid%initialize()
call grid%set_branch_factor(2_i32_k, status)
call grid%set_dimensions(3_i32_k, 3_i32_k, 1_i32_k, status)
call new_non_oriented_cursor(grid, 0_id_k, cursor, create=.true., status=status)
```

**Pointer Stability Across Dynamic Growth**:
- Internally, `htg_grid_t` boxes each sparse `htg_tree_t` inside an allocatable `htg_tree_slot_t%tree`, and `htg_cell_data_t` boxes each `htg_field_array_t` inside an allocatable `htg_field_box_t%item`. When the forest or `CellData` container grows, slots are transferred via `move_alloc` (Fortran 2018 §16.9.137), **preserving the heap addresses** of all existing trees and field arrays so active cursors and `get_array` pointers remain valid.
- **View Target Rule**: Declare `target` on any standalone `htg_grid_t`, `htg_tree_t`, `htg_bit_array_t`, `htg_cell_data_t`, or `htg_field_array_t` variable before calling zero-copy pointer getters (`get_*_view`, `get_raw_bytes`, `get_tree`, `get_cell_data`, `get_mask`, `get_pure_mask`). Note that low-level data-buffer pointers (`get_view_r64`, `get_parent_to_elder_child_view`, `get_raw_bytes`) point directly into growable value arrays and should be re-queried if further subdivisions or resizes occur.

### 3.3 Error Handling (`htg_status_t` vs. Fail-Fast)

Every fallible `fvtkhtg` procedure accepts an optional trailing `status` argument of `type(htg_status_t)`:

1. **Checked Mode (`status` supplied)**: Never aborts execution. On success, `status%ok()` is `.true.` (`status%code == HTG_OK`). On failure, `status%code` holds a specific diagnostic code and `trim(status%message)` explains the violated invariant.
2. **Fail-Fast Mode (`status` omitted)**: Immediately terminates via standard Fortran 2018 `error stop msg` if any precondition or structural check fails.

```fortran
type(htg_status_t) :: status

! Checked mode: inspect status%ok(), status%code, and status%message
call grid%set_branch_factor(5_i32_k, status)
if (.not. status%ok()) then
   write(*, '(a,i0,a,a)') 'Error [', status%code, ']: ', trim(status%message)
end if

! Fail-fast mode: omitting `status` invokes `error stop` on invalid input
call grid%set_branch_factor(2_i32_k)
```

#### Kind Parameters & Error Code Constants

| Symbol | Kind / Value | Purpose |
| :--- | :--- | :--- |
| `id_k`, `i64_k` | `iso_fortran_env::int64` | 64-bit tree, vertex, global node, tuple, and bit IDs (`vtkIdType`) |
| `i32_k` | `iso_fortran_env::int32` | 32-bit dimensions, branch factors, levels, child & component indices |
| `byte_k` | `iso_fortran_env::int8` | Raw MSB-first packed bit-array storage byte |
| `r64_k`, `r32_k` | `real64`, `real32` | Double- and single-precision coordinates and `CellData` fields |
| `HTG_OK` | `0_i32_k` | Operation succeeded |
| `HTG_ERR_INVALID_ARG` | `1_i32_k` | Invalid dimension, branch factor, non-monotonic coordinate, or bad slice |
| `HTG_ERR_OUT_OF_BOUNDS` | `2_i32_k` | `tree_id`, `(i,j,k)`, `local_id`, `ichild`, `tuple_id`, or `bit_id` out of range |
| `HTG_ERR_INVALID_STATE` | `3_i32_k` | Uninitialized object, missing dimensions, or mutating frozen grid geometry |
| `HTG_ERR_TOPOLOGY` | `4_i32_k` | Descending into a leaf, malformed BF descriptor, or broken tree invariant |
| `HTG_ERR_INDEXING` | `5_i32_k` | Missing, mixed implicit/explicit, negative, or duplicate/overlapping global IDs |
| `HTG_ERR_MASK` | `6_i32_k` | Subdividing/descending into a masked cell, or invalid `Mask` / `PureMask` state |
| `HTG_ERR_INTERFACE` | `7_i32_k` | Missing or non-3-component `Normals` / `Intercepts` arrays when `has_interface` is true |
| `HTG_ERR_NOT_FOUND` | `8_i32_k` | Named `CellData` array or `tree_id` not present during lookup/removal |

---

## 4. Step-by-Step Tutorials & Annotated Patterns

### 4.1 Grid Setup, Squeezed Orientations & Rectilinear Coordinates

A `vtkHyperTreeGrid` is defined at level zero by **point counts** `(dx, dy, dz)` where every axis satisfies `d >= 1`. Passing `1` along an axis squeezes (collapses) that axis so `cell_dims(a) = 1`; passing `d > 1` creates `d - 1` root tree cells along that axis.

> [!WARNING]
> Always call `set_branch_factor(2_i32_k)` (or `3_i32_k`) and `set_dimensions(...)` (and optional `set_transposed_root_indexing`) **before** instantiating any trees. Once non-empty trees exist in `grid`, changing `dimensions`, `branch_factor`, or `transposed_root_indexing` fails with `HTG_ERR_INVALID_STATE`.

#### Squeezed Geometry & VTK `v9.7+` 2D XZ Cyclic Axis Convention

`fvtkhtg` computes `dimension`, `orientation`, and active plane axes identically to VTK `v9.7.2`:

| Grid Configuration | `dimensions` | `cell_dims` | `get_dimension()` | `get_orientation()` | `get_1d_axis` / `get_2d_axes` | `get_axes()` |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: |
| 0D degenerate | `[1, 1, 1]` | `[1, 1, 1]` | `0` | `-1` | *(invalid state)* | `[-1, -1, -1]` |
| 1D along X | `[nx, 1, 1]` | `[nx-1, 1, 1]` | `1` | `0` (`HTG_AXIS_X`) | `axis = 0` | `[0, -1, -1]` |
| 1D along Y | `[1, ny, 1]` | `[1, ny-1, 1]` | `1` | `1` (`HTG_AXIS_Y`) | `axis = 1` | `[1, -1, -1]` |
| 1D along Z | `[1, 1, nz]` | `[1, 1, nz-1]` | `1` | `2` (`HTG_AXIS_Z`) | `axis = 2` | `[2, -1, -1]` |
| 2D in YZ plane | `[1, ny, nz]` | `[1, ny-1, nz-1]` | `2` | `0` (normal X) | `axis1 = 1, axis2 = 2` | `[1, 2, -1]` |
| **2D in XZ plane** | **`[nx, 1, nz]`** | **`[nx-1, 1, nz-1]`** | **`2`** | **`1` (normal Y)** | **`axis1 = 2, axis2 = 0`** | **`[2, 0, -1]`** |
| 2D in XY plane | `[nx, ny, 1]` | `[nx-1, ny-1, 1]` | `2` | `2` (normal Z) | `axis1 = 0, axis2 = 1` | `[0, 1, -1]` |
| 3D volume | `[nx, ny, nz]` | `[nx-1, ny-1, nz-1]` | `3` | `-1` | *(invalid state)* | `[0, 1, 2]` |

> [!NOTE]
> For a **2D XZ** grid (`dimensions = [nx, 1, nz]`), VTK `v9.7+` orders the active plane axes cyclically right-handed around normal Y (`orientation = 1`): `(orientation + 1) mod 3 = 2` (Z) first and `(orientation + 2) mod 3 = 0` (X) second, so `get_2d_axes` returns `(2, 0)` rather than `(0, 2)`.
>
> This order is grid metadata only. **Child numbering does not follow it**: in every plane, the lowest digit of a child index belongs to the lower-numbered in-plane axis, so in an XZ grid child 1 is the `+x` child and child 2 the `+z` child. Use `grid%get_active_axes` (always ascending) and `grid%get_child_position(ichild, axis)` when relating child indices to directions (Section [4.10](#410-cell-geometry)).

```fortran
use fvtkhtg, only: htg_grid_t, htg_status_t, id_k, i32_k, r64_k

type(htg_grid_t), target :: grid
type(htg_status_t)       :: status
real(r64_k)              :: x(0:3), z(0:2)
integer(i32_k)           :: ax1, ax2
integer(id_k)            :: tree_id

call grid%initialize()
call grid%set_branch_factor(2_i32_k, status)        ! Binary quadtrees (k = 2^2 = 4)
call grid%set_dimensions(4_i32_k, 1_i32_k, 3_i32_k, status) ! 2D XZ: 3 x 1 x 2 = 6 trees

! Coordinates must be finite and strictly monotonic (ascending OR descending).
! Squeezed axis Y (dim == 1) automatically defaults to [0.0_r64_k] if omitted.
x = [0.0_r64_k, 1.0_r64_k, 2.5_r64_k, 5.0_r64_k]
z = [10.0_r64_k, 5.0_r64_k, 0.0_r64_k]              ! Descending coordinates are valid
call grid%set_x_coordinates(x, status)
call grid%set_z_coordinates(z, status)

call grid%get_2d_axes(ax1, ax2, status)             ! ax1 == 2 (Z), ax2 == 0 (X)

! Map 0-based level-zero cell coordinates (i=2, j=0, k=1) -> tree_id:
!   Default KJI (transposed_root_indexing = .false., i fastest):
!     tree_id = i + j*nx + k*nx*ny = 2 + 0*3 + 1*3*1 = 5
!   Transposed IJK (transposed_root_indexing = .true., k fastest):
!     tree_id = k + j*nz + i*nz*ny
call grid%get_index_from_level_zero_coordinates(tree_id, 2_i32_k, 0_i32_k, 1_i32_k, status)
```

---

### 4.2 Cursor-Driven AMR Refinement

Use `new_non_oriented_cursor(grid, tree_id, cursor, create=.true., status=status)` to instantiate a root HyperTree on demand and position a stack cursor at its root (`vertex_id = 0`, `level = 0`). Pass `create=.false.` (the default) to inspect a slot without allocating an empty tree (`cursor%has_tree()` will be `.false.`).

```fortran
use fvtkhtg, only: htg_grid_t, htg_cursor_t, htg_status_t, &
   new_non_oriented_cursor, id_k, i32_k, r64_k

type(htg_grid_t), target :: grid
type(htg_cursor_t)       :: cursor
type(htg_status_t)       :: status

call grid%initialize()
call grid%set_branch_factor(2_i32_k, status)
call grid%set_dimensions(3_i32_k, 3_i32_k, 1_i32_k, status)
call grid%set_x_coordinates([0.0_r64_k, 0.5_r64_k, 1.0_r64_k], status)
call grid%set_y_coordinates([0.0_r64_k, 0.5_r64_k, 1.0_r64_k], status)

! Create tree 0 and subdivide root -> spawns 4 child leaves (local IDs 1..4)
call new_non_oriented_cursor(grid, 0_id_k, cursor, create=.true., status=status)
call cursor%subdivide_leaf(status)

! Descend into quadrant 3 (0 <= ichild < 4), refine it, and return to parent
call cursor%to_child(3_i32_k, status)
if (cursor%is_leaf(status)) then
   call cursor%subdivide_leaf(status)
end if
call cursor%to_parent(status)
```

Key cursor behaviors matching VTK `v9.7.2`:
- **`cursor%subdivide_leaf` vs. `tree%subdivide_leaf`**: Calling `cursor%subdivide_leaf` on an already-refined coarse vertex is a **silent no-op** (`HTG_OK`), matching `vtkHyperTreeGridEntry::SubdivideLeaf`. Calling `tree%subdivide_leaf(local_id, level)` directly on a coarse vertex fails with `HTG_ERR_TOPOLOGY`.
- **`is_terminal_node`**: Returns `.true.` if and only if the current vertex is a **coarse** (`not is_leaf`) node whose $k$ children are all leaves.
- **Depth Limiter (`grid%set_depth_limiter(limit)`)**: Caps cursor traversal depth so `cursor%is_leaf()` reports `.true.` once `cursor%get_level() >= limit`.

See runnable program: [`examples/01_basic_quadtree.f90`](../examples/01_basic_quadtree.f90).

---

### 4.3 Implicit vs. Explicit Global Node Indexing

Global node IDs map every tree's local vertices (`0 .. nv - 1`) into a single 0-based forest-wide address space used by `CellData`, `Mask`, `PureMask`, and `breadth_first_id_map`. Each tree operates in one of three modes (`tree%get_indexing_mode()`):

1. **`HTG_INDEX_UNSET`**: Initial state before global indices are assigned.
2. **`HTG_INDEX_IMPLICIT`**: Stores a single `global_index_start` offset per tree so `global_id = start + local_id` with zero per-vertex memory overhead.
   - **Whole-forest helper**: `call grid%initialize_local_index_node(status)` (alias `shift_all_trees_implicit_global_index`) walks active trees in ascending `tree_id` order and assigns contiguous `[0 .. N_cells - 1]` implicit offsets.
   - **Per-tree offset**: `call cursor%set_global_index_start(offset, status)` (or `call tree%set_global_index_start(offset, status)`).
3. **`HTG_INDEX_EXPLICIT`**: Stores a per-vertex `global_index_table(0:nv-1)` allowing non-contiguous or permuted global IDs (e.g., matching an external AMR solver's numbering).
   - **Via cursor (1 ID argument)**: `call cursor%set_global_index_from_local(global_id, status)` sets the global ID of the cursor's current vertex.
   - **Via tree (2 ID arguments)**: `call tree%set_global_index_from_local(local_id, global_id, status)`.

> [!NOTE]
> Implicit and explicit modes cannot be mixed **within the same tree** (`HTG_ERR_INDEXING`), though different trees in the same `htg_grid_t` may use implicit and explicit modes simultaneously as long as their assigned global IDs do not overlap.

---

### 4.4 Dynamic `CellData` & Material Interfaces

`htg_grid_t` owns a `htg_cell_data_t` container (`call grid%get_cell_data(cell_data)`) of named cell-centered arrays (`htg_field_array_t`), indexed by 0-based `global_id`:

- Supported types: `HTG_TYPE_REAL64`, `HTG_TYPE_REAL32`, `HTG_TYPE_INT64`, `HTG_TYPE_INT32`, `HTG_TYPE_BIT`.
- **Automatic Tuple Growth**: Whenever `CellData` contains at least one array and a tree has active global indexing, calling `cursor%subdivide_leaf`, `cursor%set_global_index_start`, `cursor%set_global_index_from_local`, `grid%initialize_local_index_node`, or `grid%set_tree` automatically invokes `cell_data%ensure_tuples(max_gid + 1_id_k)` with geometric reallocation. (If you add a new array *after* global IDs have already been finalized, pass `num_tuples = grid%get_max_global_index() + 1_id_k` to `add_array` or call `cell_data%ensure_tuples(...)`.)

```fortran
type(htg_cell_data_t),   pointer :: cell_data => null()
type(htg_field_array_t), pointer :: density   => null()
type(htg_field_array_t), pointer :: velocity  => null()
real(r64_k)                      :: vel(0:2)
integer(id_k)                    :: gid, n_tuples

! Add arrays before `initialize_local_index_node` / `subdivide_leaf` for automatic
! growth, or size explicitly to `max_global_index + 1` on an already-indexed grid:
n_tuples = grid%get_max_global_index() + 1_id_k
call grid%get_cell_data(cell_data)
call cell_data%add_array('Density',  HTG_TYPE_REAL64, &
   num_components=1_i32_k, num_tuples=n_tuples, status=status)
call cell_data%add_array('Velocity', HTG_TYPE_REAL64, &
   num_components=3_i32_k, num_tuples=n_tuples, status=status)
call cell_data%set_active_scalars('Density', status)

call cell_data%get_array('Density',  density,  status)
call cell_data%get_array('Velocity', velocity, status)

! Write scalar and 0-based 3-component vector values at the cursor's global ID
gid = cursor%get_global_node_index(status)
call density%set_scalar_r64(gid, 1.225_r64_k, status)
vel = [10.0_r64_k, -2.5_r64_k, 0.0_r64_k]
call velocity%set_tuple_r64(gid, vel, status)
```

#### Material Interface Convention (`Normals` & `Intercepts`)

To describe sub-cell planar material interfaces matching VTK `v9.7+`:
1. Add two **3-component** arrays (e.g., `'Normals'` and `'Intercepts'`) to `CellData`.
2. Enable interface metadata on the grid:
   ```fortran
   call grid%set_has_interface(.true.)
   call grid%set_interface_normals_name('Normals')
   call grid%set_interface_intercepts_name('Intercepts')
   ```
3. Store `[nx, ny, nz]` in `Normals` and `[d1, d2, itype]` in `Intercepts` (`0:2` components). Following VTK `vtkHyperTreeGrid::RecursivelyInitializePureMask`, component `2` (`itype`) classifies unmasked leaf cells:
   - **`Intercepts(2) < 2.0`** (`-1`, `0`, or `1` in VTK convention): leaf contains an **active interface** (mixed/impure cell).
   - **`Intercepts(2) >= 2.0`** (canonically `2.0`): leaf is a **pure single-material cell**.

See runnable program: [`examples/02_octree_amr_celldata.f90`](../examples/02_octree_amr_celldata.f90).

---

### 4.5 Cell Masking (`Mask`) & Lazy Bottom-Up `PureMask`

#### Cell Masking (`Mask`)
Blanked/inactive cells are tracked in an MSB-first `htg_bit_array_t` indexed by `global_id`:
- Set per-cell via `call cursor%set_mask(.true., status)` or `call grid%set_cell_mask(gid, .true., status)`, or bulk-assign via `call grid%set_mask(mask_bits, status)`.
- Query via `grid%has_mask()` (`.true.` iff a mask is attached **and** has `> 0` bits, matching `vtkHyperTreeGrid::HasMask()`), `cursor%is_masked(status)`, or `grid%is_cell_masked(gid, status)`.
- **Preconditions**: Calling `cursor%subdivide_leaf` on a masked leaf or `cursor%to_child` on a masked coarse vertex fails with `HTG_ERR_MASK` (matching VTK `v9.7.2`'s `assert(!IsMasked)` preconditions). When `grid%has_mask()` is `.false.`, `cursor%subdivide_leaf` and `cursor%to_child` short-circuit without querying global node indices, so you can refine trees freely before assigning global IDs.

#### Lazy `PureMask` (`get_pure_mask`)
`call grid%get_pure_mask(pure_mask_ptr, status)` lazily computes and caches a `max_global_index + 1`-bit mask (`1` bit `[0]` on an empty grid) where:
- **`.false.` (`0`) = PURE cell** (unmasked leaf with no active material interface, or unmasked coarse cell whose $k$ children are all pure).
- **`.true.` (`1`) = IMPURE cell** (masked cell, unmasked leaf with `Intercepts[gid][2] < 2.0`, or unmasked coarse ancestor of any impure cell).

Exact VTK `v9.7.2` bottom-up evaluation rules:
1. **Masked Short-Circuit First**: If a coarse or leaf vertex is masked (`is_cell_masked(gid)`), `PureMask(gid)` is set to `.true.` immediately **without** descending into its children (leaving unreachable descendant bits at `.false.`).
2. **Unmasked Leaf**: `PureMask(gid)` is `.true.` iff valid 3-component interface arrays exist and `Intercepts(2, gid) < 2.0`.
3. **Unmasked Coarse Vertex**: Visits **all** $k$ children and sets `PureMask(gid)` to the logical OR of its children's `PureMask` bits.

`fvtkhtg` automatically invalidates its cached `PureMask` (`grid%is_pure_mask_valid() == .false.`) whenever tree topology, global indices, cell masks, or interface metadata change. (If you mutate `Intercepts` values directly inside `CellData` after computing `PureMask`, call `call grid%invalidate_pure_mask()` before the next `get_pure_mask` query.)

See runnable program: [`examples/03_mask_puremask_interface.f90`](../examples/03_mask_puremask_interface.f90).

---

### 4.6 Breadth-First Order Bit Descriptors

VTK's `.htg` and VTKHDF formats represent HyperTree topologies compactly as Breadth-First Order bit streams where `1` (`.true.`) denotes a coarse vertex and `0` (`.false.`) denotes a leaf vertex.

> [!IMPORTANT]
> **Deepest Level Omission & VTK `v9.7+` 0-Bit Single-Leaf Root Convention**:
> Because every vertex on the deepest level `number_of_levels - 1` is necessarily a leaf in the output tree, VTK writes descriptor bits **only for depths `0 .. number_of_levels - 2`**.
> Consequently, a **1-level single-leaf root tree** (`number_of_levels == 1`, `number_of_vertices == 1`) is canonically represented by **`number_of_bits = 0_id_k`**!

#### Building a Tree from a Descriptor Slice
```fortran
type(htg_bit_array_t), target :: desc
type(htg_tree_t),     pointer :: tree_ptr => null()

! Build a 2D quadtree from 5 bits (`1 | 1 0 0 1`) starting at bit offset 0:
!   depth 0: [1]          -> root subdivided into 4 children
!   depth 1: [1, 0, 0, 1] -> children 0 and 3 subdivided (+8 leaves at depth 2)
!   depth 2: (8 leaves omitted) -> 3 levels, 13 vertices, 10 leaves
call desc%initialize(5_id_k, .false., status)
call desc%set_value(0_id_k, .true., status)
call desc%set_value(1_id_k, .true., status)
call desc%set_value(4_id_k, .true., status)

call grid%get_tree(0_id_k, tree_ptr, create=.true., status=status)
call tree_ptr%build_from_breadth_first_order_descriptor( &
   desc, 5_id_k, start_index=0_id_k, status=status)

! Build a 1-level single-leaf root in tree 1 using `number_of_bits = 0_id_k`
call grid%get_tree(1_id_k, tree_ptr, create=.true., status=status)
call tree_ptr%build_from_breadth_first_order_descriptor( &
   desc, 0_id_k, start_index=5_id_k, status=status)
```

#### Exporting Forest Descriptors (`compute_breadth_first_order_descriptor`)
`compute_breadth_first_order_descriptor` **appends** to `descriptor`, `number_of_vertices_per_depth` (`0:` allocatable), and `breadth_first_id_map` (`0:` allocatable of **global** node IDs, requiring global indexing on the tree), allowing a single loop across `get_tree_by_pos` to assemble forest-wide serialization buffers:

```fortran
type(htg_bit_array_t), target :: forest_desc
integer(id_k), allocatable    :: forest_vpd(:), forest_id_map(:)
integer(id_k)                 :: pos, tid

do pos = 0_id_k, grid%get_number_of_non_empty_trees() - 1_id_k
   call grid%get_tree_by_pos(pos, tid, tree_ptr, status)
   call tree_ptr%compute_breadth_first_order_descriptor( &
      forest_desc, forest_vpd, forest_id_map, status=status)
end do
```

Optional `input_mask` (prunes descendants of masked coarse nodes while recording the masked node itself as a `0` leaf bit) and `depth_limiter` (caps exported depth) are also supported as trailing keyword arguments.

---

### 4.7 Structural Invariant Validation (`grid%validate`)

Before handing an in-memory `htg_grid_t` to a solver or file serializer, run a complete structural audit:

```fortran
call grid%validate(status)
if (.not. status%ok()) then
   error stop 'Invalid HyperTreeGrid: ' // trim(status%message)
end if
```

`grid%validate(status)` (and standalone `htg_validate_grid(grid, [status], [allow_unset_global_index])` / `htg_validate_tree(tree, [status], [require_global_index])`) verifies in $O(N \log N)$ / $O(N)$ time with zero recursion:

1. **Grid Geometry**: `dimensions >= 1`, `branch_factor` $\in \{2, 3\}$, and exact consistency of cached `cell_dims`, `dimension`, `orientation`, `axes`, and `number_of_children`.
2. **Coordinates**: Presence on non-squeezed axes, `0:` lower bounds, exact length match with `dimensions`, `ieee_is_finite` on every coordinate, and strict ascending or descending monotonicity.
3. **Sparse Forest & Tree Topology**: Strictly sorted `tree_id` slots in `[0, max_trees - 1]`, $k$-ary identities $N_v = 1 + k N_n$ and $N_l = 1 + (k - 1) N_n$, `ParentToElderChild` bounds, full iterative DFS reachability (no cycles, no shared children, no orphan vertices), exact match with `number_of_vertices_per_depth`, and agreement of the parent table (used by `get_parent_index`, `get_level` and `coarsen_node`) with `ParentToElderChild`.
4. **Global Indexing**: No mixed unset/set trees, no negative IDs, and zero overlapping or duplicate global IDs across all implicit intervals and explicit index tables.
5. **`CellData` Integrity**: `num_tuples >= max_global_index + 1` on every field array, valid `active_scalars` reference, and zeroed trailing slack bits in `HTG_TYPE_BIT` arrays.
6. **Material Interfaces**: Existence and `num_components == 3` on both `interface_normals_name` and `interface_intercepts_name` whenever `has_interface` is `.true.`.
7. **`Mask` & `PureMask`**: Zeroed trailing slack bits in `mask`, and bit-for-bit agreement of any valid cached `pure_mask` with a fresh bottom-up recomputation.

---

### 4.8 Zero-Copy Buffer Extraction for External Serializers

Downstream file writers (`.htg` XML, VTKHDF/HDF5) and MPI halo-exchange libraries can extract direct `0:`-bounded Fortran pointers into `fvtkhtg`'s internal buffers without copying:

- **Coordinate Views**: `grid%get_x_coordinates_view(x_ptr)`, `get_y_coordinates_view(y_ptr)`, `get_z_coordinates_view(z_ptr)` $\rightarrow$ `real(r64_k), pointer :: x_ptr(:)` bounded `0 : dims(a) - 1`.
- **Tree Topology Views**:
  - `tree_ptr%get_parent_to_elder_child_view(elder_ptr)` $\rightarrow$ `integer(id_k), pointer :: elder_ptr(:)` bounded `0 : tree_ptr%get_elder_child_array_size() - 1` (leaf entries hold `HTG_TREE_LEAF_SENTINEL = 4294967295_id_k`; vertices `>= size` are leaves by VTK's lazy tail rule).
  - `tree_ptr%get_vertices_per_depth_view(vpd_ptr)` $\rightarrow$ `integer(id_k), pointer :: vpd_ptr(:)` bounded `0 : number_of_levels - 1`.
  - `tree_ptr%get_global_index_table_view(gid_ptr)` $\rightarrow$ `integer(id_k), pointer :: gid_ptr(:)` bounded `0 : number_of_vertices - 1` when `get_indexing_mode() == HTG_INDEX_EXPLICIT` (`null()` otherwise).
- **MSB-First Packed Bit Bytes (`vtkBitArray` Layout)**:
  - `bits%get_raw_bytes(byte_ptr)` $\rightarrow$ `integer(byte_k), pointer :: byte_ptr(:)` bounded `0 : bits%get_number_of_bytes() - 1` (`null()` when `num_bits == 0`).
  - Bit `id` resides at `byte_ptr(id / 8)` in bit position `7 - mod(id, 8)` (`0x80 >> (id % 8)`), with all unused trailing bits in the final active byte guaranteed `0`.
- **2D `CellData` Field Views (AOS Memory Layout)**:
  - `field_ptr%get_view_r64(r64_ptr, status)` (and `_r32`, `_i64`, `_i32`) $\rightarrow$ `ptr(0:nc-1, 0:nt-1)`.
  - Because Fortran stores 2D arrays column-major, `ptr(0:nc-1, 0:nt-1)` places components `0 .. nc - 1` of tuple `0` contiguously in memory, followed by tuple `1`, etc.—**byte-for-byte identical** to C++ VTK's `vtkAOSDataArrayTemplate` tuple buffer!

```fortran
real(r64_k),     pointer :: x_view(:)      => null()
integer(id_k),   pointer :: elder_view(:)  => null()
integer(byte_k), pointer :: desc_bytes(:)  => null()
real(r64_k),     pointer :: p_view(:, :)   => null()

call grid%get_x_coordinates_view(x_view)             ! x_view(0 : nx - 1)
call tree_ptr%get_parent_to_elder_child_view(elder_view) ! elder_view(0 : n_active - 1)
call forest_desc%get_raw_bytes(desc_bytes)           ! desc_bytes(0 : n_bytes - 1), MSB-first
call pressure%get_view_r64(p_view, status)           ! p_view(0 : nc - 1, 0 : nt - 1)
```

See runnable program: [`examples/04_bf_descriptors_and_views.f90`](../examples/04_bf_descriptors_and_views.f90).

### 4.9 Writing VTKHDF Files: Breadth-First Export & the `fvtkhdf` Seam

`fvtkhtg` builds and validates grids; [`fvtkhdf`](https://github.com/omdaniel/fvtkhdf/tree/feature/hypertreegrid)'s `vtkhdf_htg_file` writes them as VTKHDF HyperTreeGrid files (readable by ParaView 6.2+ / VTK 9.7+), including MPI partitions and time steps. Neither core library depends on the other: the writer takes plain arrays, and `fvtkhtg` produces exactly those arrays.

| Concern | Owner |
| :--- | :--- |
| Building and refining trees, cursors, global IDs, `CellData` | `fvtkhtg` |
| `Mask`, depth limiter and interface metadata as a model; `grid%validate` | `fvtkhtg` |
| Breadth-first cell order, pruning under masked cells and at the depth limiter, the `CellData` permutation | `fvtkhtg` (`grid%export_breadth_first`, `bf%gather`) |
| Bit packing, per-tree depth counts, offsets, HDF5 datasets, MPI parts, time steps | `fvtkhdf` |
| 1-based tree numbers $\rightarrow$ the file's 0-based `TreeIds` | `fvtkhdf` |
| Calling one from the other | optional `fvtkhtg_vtkhdf` bridge |

#### The Export (`htg_breadth_first_t`)

```fortran
type(htg_breadth_first_t) :: bf
real(r64_k), allocatable  :: pressure_bf(:)

call grid%validate(status)
call grid%export_breadth_first(bf, status [, force_mask])
call bf%gather(pressure_field, pressure_bf, status)   ! CellData in export order
```

| Component | Shape | Contents |
| :--- | :--- | :--- |
| `bf%tree_ids` | `(ntree)` | **1-based** root-cell number of each non-empty tree (`tree_id + 1`), ascending |
| `bf%xtree` | `(ntree+1)` | **1-based** start of each tree's cells; `xtree(1) == 1`, `xtree(ntree+1) == ncell + 1` |
| `bf%refined` | `(ncell)` | `.true.` for a refined cell, per cell in breadth-first order |
| `bf%mask` | `(ncell)` | Mask bit per cell; allocated only if the grid has a mask or `force_mask = .true.` |
| `bf%global_ids` | `(ncell)` | The cell's **0-based** `fvtkhtg` global ID (its `CellData` tuple) |
| `bf%identity_order` | scalar | `.true.` when `global_ids(n) == n - 1`, i.e. `CellData` needs no permutation |

Within a tree, cells are listed level by level; the children of each refined cell are contiguous, in child order, and follow the order of their parents. This is the order of VTK's `ComputeBreadthFirstOrderDescriptor`, which VTK's own `vtkHDFWriter` uses, and the order `fvtkhdf` requires. As in VTK's writer, a **masked cell** or a cell **at the depth limiter** is exported as a leaf and its descendants are dropped, so no refined cell is ever masked (an `fvtkhdf` requirement). Trees must have global indices before export.

`bf%gather(field, values, status)` permutes a `CellData` array from global-ID order into export order: rank-1 `values(ncell)` for single-component fields, rank-2 `values(ncomp, ncell)` for any field. `BIT` fields gather into `integer(byte_k)` 0/1 values, since `fvtkhdf` has no bit type.

#### Writing with `fvtkhdf` Directly

```fortran
use vtkhdf_htg_file_type
type(vtkhdf_htg_file) :: file

call file%create('amr.vtkhdf', stat, errmsg)                      ! serial fvtkhdf build
call file%write_mesh(int(grid%get_dimensions()), int(grid%get_branch_factor()), &
                     grid%get_transposed_root_indexing())
call file%write_mesh_geometry(grid%x_coordinates, grid%y_coordinates, grid%z_coordinates)
call file%write_mesh_topology(bf%tree_ids, bf%xtree, bf%refined, bf%mask)
call file%write_cell_data('Pressure', pressure_bf)
call file%close()
```

An unallocated `bf%mask` is passed as an absent optional argument.

#### Writing with the `fvtkhtg_vtkhdf` Bridge

Configure with `-DFVTKHTG_BUILD_VTKHDF_BRIDGE=ON` (and `CMAKE_PREFIX_PATH` pointing at an `fVTKHDF` install), link `fvtkhtg::fvtkhtg_vtkhdf`, then:

```fortran
use fvtkhtg_vtkhdf
call htg_vtkhdf_write_static(file, grid, status)   ! validate, mesh, geometry, topology, all CellData
```

For time-dependent output the caller keeps control of `start_time_step` / `finalize_time_step` and uses the pieces:

| Routine | Wraps |
| :--- | :--- |
| `htg_vtkhdf_write_mesh(file, grid)` | `write_mesh` (once per file) |
| `htg_vtkhdf_write_geometry(file, grid)` | `write_mesh_geometry` |
| `htg_vtkhdf_write_topology(file, grid, bf [, status, force_mask])` | `export_breadth_first` + `write_mesh_topology`; returns `bf` |
| `htg_vtkhdf_write_cell_data(file, bf, field [, status, name, attribute])` | `gather` + static `write_cell_data` |
| `htg_vtkhdf_write_all_cell_data(file, grid, bf [, status])` | every array; active scalars get the `Scalars` role |
| `htg_vtkhdf_register_temporal_cell_data(file, field, handle [, ...])` | `register_temporal_cell_data` with a matching mold |
| `htg_vtkhdf_write_temporal_cell_data(file, handle, bf, field [, status])` | `gather` + temporal `write_cell_data` |

In MPI builds of `fvtkhdf` every call is collective: each rank passes the trees it owns (no tree may be split between ranks), all ranks must agree on dimensions, branch factor, root indexing and the names, order, types and component counts of their `CellData` arrays, and if any rank has a mask all ranks must pass `force_mask = .true.`.

> [!NOTE]
> **Current limitations of the seam**:
> - `fvtkhdf` does not yet write the `InterfaceNormalsName` / `InterfaceInterceptsName` attributes that VTK's reader uses for material interfaces; the `Normals` and `Intercepts` arrays are written as ordinary cell data.
> - `fvtkhtg` cannot coarsen a tree, so a solver writing a changing mesh (`VTKHDF_DYNAMIC_MESH`) rebuilds any tree that coarsens; the export handles any refinement order.

See runnable program: [`examples/05_breadth_first_export.f90`](../examples/05_breadth_first_export.f90).

### 4.10 Cell Geometry

The cursor knows where its cell is. Four queries mirror `vtkHyperTreeGridNonOrientedGeometryCursor` and return **bit-identical** values (the oracle suite compares them with `==`):

```fortran
real(r64_k) :: origin(0:2), cell_size(0:2), bounds(0:5), centre(0:2)

call cursor%get_origin(origin, status)       ! low corner (x, y, z)
call cursor%get_size(cell_size, status)      ! edge lengths; 0 along a squeezed axis
call cursor%get_bounds(bounds, status)       ! (xmin, xmax, ymin, ymax, zmin, zmax)
call cursor%get_point(centre, status)        ! cell centre
```

- **No per-cursor storage**: the geometry is recomputed from the cursor's path on each call, in O(level). The root cell comes from the grid's coordinates (`grid%get_level_zero_origin_and_size_from_index(tree_id, origin, cell_size)`), each level divides the size by the branch factor, and each step to a child adds `position * size` along every active axis — the same operations, in the same order, as VTK.
- **Coordinates are required** on every non-squeezed axis (`HTG_ERR_INVALID_STATE` otherwise). With descending coordinates the size is negative and `bounds` runs from origin to origin + size, as in VTK.
- **A root cell without a tree still has geometry**: a cursor opened with `create=.false.` on an empty root cell reports that root cell's box.

**Which child is where.** Digit `p` (base `branch_factor`) of a child index is the child's position along the `p`-th active axis in ascending order:

```fortran
integer(i32_k) :: active(0:2), nactive, pos

call grid%get_active_axes(active, nactive)            ! e.g. (0, 2, -1) for an XZ grid
pos = grid%get_child_position(ichild, HTG_AXIS_Z)     ! 0 .. branch_factor-1, or -1 if squeezed
ichild = cursor%get_child_index()                     ! which child of its parent the cell is
```

### 4.11 Face Neighbours

`get_face_neighbor` puts a second cursor on the cell across one face:

```fortran
type(htg_cursor_t) :: neighbor

call cursor%get_face_neighbor(HTG_AXIS_X, +1_i32_k, neighbor, status)   ! across the +x face
if (neighbor%has_tree()) then
   ! neighbor is an ordinary cursor: level, global id, geometry, children, ...
   if (neighbor%get_level() < cursor%get_level()) ...   ! a coarser leaf
end if
```

`axis` is a physical axis (`HTG_AXIS_X/Y/Z`) that the grid does not squeeze, `direction` is `-1` or `+1`, and `neighbor` must be a different variable from the cursor itself.

| Across the face there is ... | `neighbor` is ... |
| :--- | :--- |
| a cell of the same level | that cell. It may itself be refined: descend into the children on the shared face (`grid%get_child_position`) to reach the smaller cells |
| a coarser leaf | that leaf (`neighbor%get_level() < cursor%get_level()`) |
| a masked cell on the way down | that masked cell; the search never descends through a mask (`neighbor%is_masked()`) |
| the edge of the grid | nothing: `has_tree()` is `.false.` and `get_tree_id()` is `HTG_INVALID_ID` |
| a root cell without a tree | nothing: `has_tree()` is `.false.` and `get_tree_id()` is that root cell's index |

The lookup works inside a tree and across root cells alike, honours transposed root indexing, and costs O(level). It is what VTK's `vtkHyperTreeGridNonOrientedVonNeumannSuperCursor` holds for each face; the oracle suite checks every neighbour of every cell against it.

> [!NOTE]
> **XZ-plane grids.** `fvtkhtg` defines the direction geometrically, so the neighbour always touches the face you asked for. VTK `v9.7.2`'s supercursor looks up root-level neighbours along `GetAxes()`, which for an XZ grid is `(z, x)` while its child numbering is `(x, z)`; across root-cell boundaries of such grids it therefore holds the neighbour along the *other* axis. Inside a tree, and in every other plane, the two agree exactly.

**2:1 balancing** needs nothing more. Looked at from the finer side, a grid is unbalanced where the face neighbour of a leaf is a leaf two or more levels coarser; refine it through the neighbour cursor and repeat until nothing changes:

```fortran
call cursor%get_face_neighbor(axis, direction, neighbor, status)
if (neighbor%has_tree()) then
   neighbor_is_leaf = neighbor%is_leaf(status)
   if (neighbor_is_leaf .and. neighbor%get_level() < cursor%get_level() - 1_i32_k) then
      call neighbor%subdivide_leaf(status)
   end if
end if
```

### 4.12 Coarsening

`coarsen_node` is the inverse of `subdivide_leaf`: on a refined cell whose children are all leaves, it removes the children and leaves the cursor on the same cell, now a leaf. VTK has no equivalent.

```fortran
! Set the cell's own value from its children first: their data is not kept
call cursor%coarsen_node(status [, freed_global_ids, moved_from_global_ids, moved_to_global_ids])
```

Coarsen bottom-up: a cell with refined children is rejected with `HTG_ERR_TOPOLOGY`, as is a leaf.

**What happens to numbering.** Local vertex IDs stay contiguous (`0 .. nv-1`), as every other operation and VTK assume. Removing a block of `number_of_children` IDs leaves a gap, which is closed by renumbering the sibling block with the highest local IDs in that tree into it (unless the removed block was that block). No other cell changes. What that means for you depends on the tree's indexing mode:

| Indexing mode | Global IDs of surviving cells | `CellData` and mask | Optional outputs |
| :--- | :--- | :--- | :--- |
| **Explicit** (`set_global_index_from_local`) | Unchanged for every cell | Untouched for surviving cells | `freed_global_ids`: the removed children's IDs, free for you to reuse |
| **Implicit** (`global = start + local`) | The renumbered block gets new IDs | Moved with the renumbered cells automatically | `moved_from_global_ids(n)` $\rightarrow$ `moved_to_global_ids(n)`, for arrays you keep outside the grid; `freed_global_ids`: the IDs left unused at the end of the tree's range |

In both modes the freed IDs have their `CellData` tuples reset to zero and their mask bits cleared, so a later `subdivide_leaf` that reuses them starts clean. The cached `PureMask` is invalidated.

> [!IMPORTANT]
> - **Other cursors on the same tree** may be left pointing at a renumbered or removed cell. Re-initialise them after coarsening. The cursor that coarsened stays valid, including its whole path to the root.
> - **Implicit indexing and growth.** With `initialize_local_index_node` the trees' ID ranges are packed back to back, so a tree that grows after indexing runs into the next tree's range (`grid%validate` reports the overlap). For a mesh that refines and coarsens over time, either use explicit IDs with your own free list, or give each tree a fixed range with `cursor%set_global_index_start(tree_id * capacity)`, as [`examples/06_adaptive_refinement.f90`](../examples/06_adaptive_refinement.f90) does. The breadth-first export gathers only the IDs in use, so gaps cost nothing in the file.

Tree-level access, for code that works with local IDs rather than cursors:

| Call | Result |
| :--- | :--- |
| `tree%get_parent_index(local_id, [status])` | Parent's local ID (`HTG_INVALID_ID` for the root), O(1) |
| `tree%get_child_index(local_id, [status])` | Which child of its parent (`-1` for the root), O(1) |
| `tree%get_level(local_id, [status])` | Depth, O(level) |
| `call tree%coarsen_node(local_id, info, [status])` | Coarsens; `info` (`htg_coarsen_info_t`) reports `node_id`, `removed_first`, `moved_from`, `moved_to`. Moving cell data in implicit mode is then the caller's job |

See runnable program: [`examples/06_adaptive_refinement.f90`](../examples/06_adaptive_refinement.f90), which refines around a moving circle from cell bounds, restores 2:1 balance through face neighbours, and coarsens what the circle leaves behind.

---

## 5. C++ `vtkHyperTreeGrid` to Fortran `fvtkhtg` API Reference

### 5.1 Classes, Types & Sentinels

| C++ VTK `v9.7+` | Fortran `fvtkhtg` (`use fvtkhtg`) | Notes |
| :--- | :--- | :--- |
| `vtkHyperTreeGrid` | `type(htg_grid_t)` | Declare with `target` when creating cursors or pointer views |
| `vtkHyperTree` | `type(htg_tree_t)` | Concrete compact $k$-ary tree (replaces legacy `vtkCompactHyperTree`) |
| `vtkHyperTreeGridNonOrientedCursor` | `type(htg_cursor_t)` | Stack-based build/traversal cursor; no heap refcounting overhead |
| `vtkHyperTreeGridEntry` | `type(htg_entry_t)` | Stack frame holding `vertex_id` |
| `vtkBitArray` | `type(htg_bit_array_t)` | MSB-first bit packing (`0x80 >> (id % 8)`), trailing slack bits zeroed |
| `vtkCellData` | `type(htg_cell_data_t)` | Owned by `htg_grid_t`; accessed via `call grid%get_cell_data(ptr)` |
| `vtkDataArray` / `vtkDoubleArray` / ... | `type(htg_field_array_t)` | Tagged with `HTG_TYPE_{REAL64,REAL32,INT64,INT32,BIT}` |
| `vtkIdType` | `integer(id_k)` (`int64`) | 0-based 64-bit integer across all IDs and counts |
| `UINT_MAX` (`4294967295u`) | `HTG_TREE_LEAF_SENTINEL` (`HTG_UINT32_MAX`) | Leaf sentinel in `ParentToElderChild` |

### 5.2 `vtkHyperTreeGrid` $\leftrightarrow$ `htg_grid_t`

| C++ `vtkHyperTreeGrid` Method | Fortran `htg_grid_t` Binding | VTK `v9.7+` Parity & Usage Notes |
| :--- | :--- | :--- |
| `Initialize()` | `call grid%initialize()` / `call grid%clear()` | Resets grid; default `branch_factor` in Fortran is `2` (in C++ VTK `Initialize()` leaves `BranchFactor = 0`) |
| `SetDimensions(const int[3])` / `(i,j,k)` | `call grid%set_dimensions(ni, nj, nk, [status])` or `(dims, [status])` | Point counts `>= 1`; `1` marks a squeezed axis (`i32_k` or `i64_k`) |
| `GetDimensions(int[3])` | `dims = grid%get_dimensions()` | Returns `integer(i32_k) :: dims(0:2)` |
| `GetCellDims(int[3])` / `GetGridSizes` | `cdims = grid%get_cell_dims()` / `get_grid_sizes()` | Returns `integer(i32_k) :: cdims(0:2)` (`max(1, dims - 1)`) |
| `GetDimension()` | `ndim = grid%get_dimension()` | Active spatial dimensionality `0`, `1`, `2`, or `3` |
| `GetOrientation()` | `orient = grid%get_orientation()` | Active axis (`1D`), normal axis (`2D`), or `HTG_INVALID_ORIENTATION` (`-1`) in `0D`/`3D` |
| `Get1DAxis(unsigned int&)` | `call grid%get_1d_axis(axis, [status])` | Valid when `get_dimension() == 1` |
| `Get2DAxes(unsigned int&, unsigned int&)` | `call grid%get_2d_axes(axis1, axis2, [status])` | Cyclic right-handed: YZ `(1,2)`, **XZ `(2,0)`**, XY `(0,1)` |
| `SetBranchFactor(unsigned int)` | `call grid%set_branch_factor(factor, [status])` | `2` or `3`; must be set before trees exist |
| `GetBranchFactor()` / `GetNumberOfChildren()` | `grid%get_branch_factor()` / `grid%get_number_of_children()` | $k = f^d$ (`1` when `0D`) |
| `SetTransposedRootIndexing(bool)` | `call grid%set_transposed_root_indexing(flag, [status])` | `.false.` = KJI (`i` fastest, default); `.true.` = IJK (`k` fastest) |
| `SetDepthLimiter(unsigned int)` | `call grid%set_depth_limiter(limit, [status])` | Default `huge(0_i32_k)` |
| `Set{X,Y,Z}Coordinates(vtkDataArray*)` | `call grid%set_{x,y,z}_coordinates(coords, [status])` | `real(r64_k) :: coords(0:)`, strictly monotonic ascending or descending |
| `Get{X,Y,Z}Coordinates()` | `call grid%get_{x,y,z}_coordinates_view(ptr)` | Zero-copy `0:` pointer view (`has_{x,y,z}_coords()` checks presence) |
| `GetIndexFromLevelZeroCoordinates(idx,i,j,k)` | `call grid%get_index_from_level_zero_coordinates(tree_id, i, j, k, [status])` | Generic `i32_k` / `id_k` for `(i, j, k)` |
| `GetLevelZeroCoordinatesFromIndex(idx,i,j,k)` | `call grid%get_level_zero_coordinates_from_index(tree_id, i, j, k, [status])` | Generic `i32_k` / `id_k` for `(i, j, k)` |
| `GetLevelZeroOriginAndSizeFromIndex(idx, origin, size)` / `GetLevelZeroOriginFromIndex` | `call grid%get_level_zero_origin_and_size_from_index(tree_id, origin, cell_size, [status])` / `call grid%get_level_zero_origin_from_index(tree_id, origin, [status])` | Root cell box; size is `0` along a squeezed axis; needs coordinates on every active axis |
| *(none)* | `call grid%get_active_axes(active, nactive)` / `grid%get_child_position(ichild, axis)` | Non-squeezed axes in ascending order, and a child's position along one of them; this, not `get_axes()`, is the order child indices use |
| *(none)* | `call grid%copy_cell_attributes(from_id, to_id, [status])` / `call grid%clear_cell_attributes(global_id, [status])` | Move or reset one cell's `CellData` tuples and mask bit; used by `coarsen_node`, and useful when recycling explicit global IDs |
| `GetMaxNumberOfTrees()` | `n = grid%get_max_number_of_trees()` | Product `cell_dims(0) * cell_dims(1) * cell_dims(2)` |
| `GetNumberOfNonEmptyTrees()` | `n = grid%get_number_of_non_empty_trees()` | Count of instantiated sparse trees |
| `GetTree(vtkIdType, bool create=false)` | `call grid%get_tree(tree_id, tree_ptr, [create], [status])` | `has_tree(tree_id)` queries existence; `get_tree_by_pos(pos, tid, ptr)` iterates `0 .. num_trees - 1` |
| `SetTree(vtkIdType, vtkHyperTree*)` | `call grid%set_tree(tree_id, tree, [status])` | Copies `tree` into `tree_id` (also `remove_tree(tree_id, [status])`) |
| `GetNumberOfCells()` | `n = grid%get_number_of_cells()` (alias `get_number_of_vertices()`) | Note: C++ `vtkHyperTreeGrid` names this `GetNumberOfCells()`, not `GetNumberOfVertices()` |
| `GetNumberOfLeaves()` / `GetNumberOfLevels()` | `grid%get_number_of_leaves()` / `grid%get_number_of_levels()` | Sum of leaves / max of levels across active trees (`0` when empty) |
| `GetGlobalNodeIndexMax()` | `grid%get_global_node_index_max()` / `grid%get_max_global_index()` | `get_global_node_index_max()` returns `0` when empty/unset (matching VTK); `get_max_global_index()` returns `-1_id_k` |
| `InitializeLocalIndexNode()` | `call grid%initialize_local_index_node([status])` | Assigns contiguous implicit global offsets `[0 .. N_cells - 1]` |
| `SetMask(vtkBitArray*)` / `HasMask()` / `GetMask()` | `set_mask`, `clear_mask`, `has_mask()`, `get_mask(ptr)` | `has_mask()` is `.true.` iff mask is active **and** `num_bits > 0` |
| `GetPureMask()` | `call grid%get_pure_mask(pure_mask_ptr, [status])` | Lazy cached bottom-up mask (`.false.` = pure, `.true.` = impure); auto-invalidated on mutations |
| `GetCellData()` | `call grid%get_cell_data(cell_data_ptr)` | Returns pointer to `grid%cell_data` |
| `SetHasInterface(bool)` / `GetHasInterface()` | `call grid%set_has_interface(flag)` / `grid%get_has_interface()` | Enables 3-component material interface check in `PureMask` & `validate` |
| `SetInterfaceNormalsName(const char*)` | `call grid%set_interface_normals_name(name)` | Must name a 3-component `CellData` array when `has_interface` is true |
| `SetInterfaceInterceptsName(const char*)` | `call grid%set_interface_intercepts_name(name)` | Leaf is impure iff `Intercepts(2, gid) < 2.0` |
| *(manual assertions only in VTK)* | `call grid%validate([status])` | Comprehensive 7-stage topological, coordinate, index, mask & field audit |

### 5.3 `vtkHyperTreeGridNonOrientedCursor` $\leftrightarrow$ `htg_cursor_t`

| C++ Cursor Method | Fortran `htg_cursor_t` Binding | Notes |
| :--- | :--- | :--- |
| `grid->NewNonOrientedCursor(treeId, create)` | `call new_non_oriented_cursor(grid, tree_id, cursor, [create], [status])` | `grid` actual argument must have `target` (or `pointer`); or `call cursor%initialize(...)` |
| `HasTree()` / `GetTree()` | `cursor%has_tree()` / `call cursor%get_tree(tree_ptr)` | `.false.` / `null()` when `create=.false.` on an empty root slot |
| `GetTreeId()`*(via Tree)* / `GetVertexId()` | `cursor%get_tree_id()` / `cursor%get_vertex_id()` | 0-based `tree_id` and local `vertex_id` |
| `GetGlobalNodeIndex()` | `gid = cursor%get_global_node_index([status])` | Fails with `HTG_ERR_INDEXING` if global indexing is unset |
| `GetDimension()` / `GetNumberOfChildren()` | `cursor%get_dimension()` / `cursor%get_number_of_children()` | Matches grid/tree dimensionality and $k = f^d$ |
| `GetLevel()` / `IsRoot()` | `cursor%get_level()` / `cursor%is_root()` | Root is `level == 0`, `vertex_id == 0` |
| `IsLeaf()` / `IsTerminalNode()` | `cursor%is_leaf([status])` / `cursor%is_terminal_node([status])` | `is_leaf` honors `grid%get_depth_limiter()` |
| `ToRoot()` / `ToParent()` / `ToChild(ichild)` | `to_root([status])`, `to_parent([status])`, `to_child(ichild, [status])` | `0 <= ichild < k`; `to_child` rejects leaves (`HTG_ERR_TOPOLOGY`) & masked coarse nodes (`HTG_ERR_MASK`) |
| `SubdivideLeaf()` | `call cursor%subdivide_leaf([status])` | Rejects masked leaves (`HTG_ERR_MASK`); silent no-op on already-coarse nodes; auto-grows `CellData` |
| `SetGlobalIndexStart(offset)` | `call cursor%set_global_index_start(offset, [status])` | Sets implicit start offset for the cursor's tree |
| `SetGlobalIndexFromLocal(globalId)` | `call cursor%set_global_index_from_local(global_id, [status])` | **1 ID argument** (assigns explicit global ID to cursor's current vertex) |
| `SetMask(bool)` / `IsMasked()` | `call cursor%set_mask(state, [status])` / `cursor%is_masked([status])` | Short-circuits to `.false.` without global ID lookup when `has_mask()` is `.false.` |
| `vtkHyperTreeGridNonOrientedGeometryCursor::GetOrigin()` / `GetSize()` | `call cursor%get_origin(origin, [status])` / `call cursor%get_size(cell_size, [status])` | `(0:2)` arrays; bit-identical to VTK; recomputed from the path in O(level) |
| `...GeometryCursor::GetBounds(double[6])` / `GetPoint(double[3])` | `call cursor%get_bounds(bounds, [status])` / `call cursor%get_point(point, [status])` | `bounds(0:5)` = `(xmin, xmax, ymin, ymax, zmin, zmax)`; `point` is the cell centre |
| `vtkHyperTreeGridNonOrientedVonNeumannSuperCursor` (one face slot) | `call cursor%get_face_neighbor(axis, direction, neighbor, [status])` | Physical `axis`, `direction = -1` or `+1`; same-level or coarser cell; stops at masks; geometric direction (differs from VTK across root cells of XZ grids) |
| *(none)* | `ichild = cursor%get_child_index()` | Which child of its parent the current cell is; `-1` at a root |
| *(none: VTK cannot coarsen)* | `call cursor%coarsen_node([status], [freed_global_ids], [moved_from_global_ids], [moved_to_global_ids])` | Inverse of `subdivide_leaf`; children must all be leaves; moves `CellData`/mask in implicit mode |

### 5.4 `vtkHyperTree` $\leftrightarrow$ `htg_tree_t`

| C++ `vtkHyperTree` Method | Fortran `htg_tree_t` Binding | Notes |
| :--- | :--- | :--- |
| `Initialize(branchFactor, dimension)` | `call tree%initialize(branch_factor, dimension, [tree_index], [status])` | Validates `branch_factor` $\in \{2,3\}$, `dimension` $\in \{1,2,3\}$ |
| `GetNumberOfLevels()` / `Vertices` / `Nodes` / `Leaves` | `get_number_of_levels()`, `get_number_of_vertices()`, `get_number_of_nodes()`, `get_number_of_leaves()` | $N_v = 1 + k N_n$, $N_l = N_v - N_n$ |
| `IsLeaf(index)` / `IsTerminalNode(index)` | `tree%is_leaf(local_id, [status])` / `tree%is_terminal_node(local_id, [status])` | `is_terminal_node` checks `.not. is_leaf` first (avoids VTK `v9.7.2` C++ bug on leaf slots `< elder_size`) |
| `GetElderChildIndex(index)` | `elder = tree%get_elder_child_index(local_id, [status])` | Fails with `HTG_ERR_TOPOLOGY` if `local_id` is a leaf |
| `GetElderChildIndexArray(size_t&)` | `call tree%get_parent_to_elder_child_view(ptr)` | Active length `tree%get_elder_child_array_size()` |
| `SubdivideLeaf(index, depth)` | `call tree%subdivide_leaf(local_id, level, [status])` | Fails with `HTG_ERR_TOPOLOGY` if `local_id` is already coarse |
| *(none: VTK cannot coarsen)* | `call tree%coarsen_node(local_id, info, [status])` | Removes the (leaf) children of `local_id`; relocates the last sibling block into the gap and reports it in `info` |
| *(none: VTK keeps no parent links)* | `tree%get_parent_index(local_id, [status])`, `tree%get_child_index(local_id, [status])`, `tree%get_level(local_id, [status])` | Backed by a per-sibling-block parent table that `grid%validate` checks |
| `SetGlobalIndexStart(start)` | `call tree%set_global_index_start(start_id, [status])` | `tree%get_indexing_mode()` returns `HTG_INDEX_{UNSET,IMPLICIT,EXPLICIT}` (note: C++ `IsGlobalIndexImplicit()` is inverted in VTK `v9.7.2`) |
| `SetGlobalIndexFromLocal(local, global)` | `call tree%set_global_index_from_local(local_id, global_id, [status])` | **2 ID arguments** `(local_id, global_id)` |
| `BuildFromBreadthFirstOrderDescriptor(desc, nBits, start)` | `call tree%build_from_breadth_first_order_descriptor(descriptor, number_of_bits, [start_index], [status])` | Canonical single-leaf root uses `number_of_bits = 0_id_k` |
| `ComputeBreadthFirstOrderDescriptor(limiter, mask, vpd, desc, idMap)` | `call tree%compute_breadth_first_order_descriptor(descriptor, number_of_vertices_per_depth, breadth_first_id_map, [input_mask], [depth_limiter], [status])` | Note Fortran argument order: required outputs `(descriptor, vpd, id_map)` first, optional `input_mask`, `depth_limiter`, `status` last |

### 5.5 `vtkBitArray` & `vtkCellData` $\leftrightarrow$ `htg_bit_array_t`, `htg_cell_data_t`, `htg_field_array_t`

| C++ VTK Method | Fortran `fvtkhtg` Binding | Notes |
| :--- | :--- | :--- |
| `vtkBitArray::SetNumberOfValues(n)` | `call bits%initialize([num_bits], [init_value], [status])` / `call bits%resize(num_bits, [status])` | Alias `set_number_of_values`; clears trailing slack bits in final byte |
| `vtkBitArray::SetValue(id, v)` / `InsertValue(id, v)` | `call bits%set_value(id, val, [status])` / `call bits%insert_value(id, val, [status])` | `logical :: val`; `insert_next_value(val, [inserted_id], [status])` appends |
| `vtkBitArray::GetValue(id)` | `val = bits%get_value(id, [status])` | Returns `logical` |
| `vtkBitArray::GetPointer(0)` | `call bits%get_raw_bytes(byte_ptr)` | `integer(byte_k), pointer :: byte_ptr(:)` (`0:` bounds); `bits` needs `target` |
| `vtkCellData::AddArray(arr)` | `call cd%add_array(name, data_type, [num_components], [num_tuples], [status])` | Boxed via `move_alloc` so `get_array` pointers never dangle on growth |
| `vtkCellData::GetArray(name)` / `(idx)` | `call cd%get_array(name, ptr, [status])` / `call cd%get_array_by_index(idx, ptr, [status])` | `0 <= idx < cd%get_number_of_arrays()` |
| `vtkCellData::RemoveArray(name)` | `call cd%remove_array(name, [status])` | Clears `active_scalars` if it named the removed array |
| `vtkCellData::SetActiveScalars(name)` | `call cd%set_active_scalars(name, [status])` | Query via `has_active_scalars()` and `get_active_scalars()` |
| `vtkDataArray::SetComponent(t, c, v)` / `GetComponent` | `set_component_<t>(tuple_id, comp_id, val, [status])` / `get_component_<t>(...)` | `<t>` $\in$ `{r64, r32, i64, i32, bit}`; also `set_scalar_<t>` / `get_scalar_<t>` |
| `vtkDataArray::SetTuple(t, vals)` / `GetTuple` | `set_tuple_<t>(tuple_id, vals, [status])` / `get_tuple_<t>(tuple_id, vals, [status])` | `vals(0:)` must have `size(vals) == num_components` |
| `vtkAOSDataArrayTemplate::GetPointer(0)` | `call field%get_view_<t>(ptr, [status])` | `ptr(0:nc-1, 0:nt-1)` zero-copy AOS buffer (`get_view_bit` / `get_view_bit_bytes` for `BIT`) |
