# VTK `v9.7.2` C++ Source & LLDB Ground-Truth Audit

## 1. Pinned VTK Oracle Release & Build Configuration

- **Repository**: `https://gitlab.kitware.com/vtk/vtk.git`
- **Pinned Release Tag**: `v9.7.2` (annotated tag `695092e25b8aff988ae0f116f03284923f24bb33` $\rightarrow$ commit `0558d4d4e4bc86c4b4bf7af2f2488011a4de9e0f`)
- **Runtime Version Verification**: `vtkVersion::GetVTKVersion() == "9.7.2"` (`9.7.2`, strictly `>= 9.7.0`)
- **Install Prefix**: `.cache/vtk-oracle/install` under the repository root (where `cmake/VtkOracle.cmake` looks for it)
- **Enabled Modules**: `VTK::CommonCore`, `VTK::CommonDataModel` (all Rendering, StandAlone, Python, and Testing disabled)
- **Debug Compilation Flags**: `-DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG="-g -O0 -fno-inline"` with AppleClang `clang++` (`arm64-apple-darwin`) and Ninja `1.13.2`

---

## 2. `vtkBitArray` MSB-First Bit Packing (`Common/Core/vtkBitArray`)

- **Storage Type**: `unsigned char* Array` (`ValueType = unsigned char`).
- **Bit Indexing Formula** (`vtkBitArray::SetValue(vtkIdType id, int value)`):
  - Byte index: `bitsetId = id / 8`
  - Bit offset within byte: `bitId = id % 8`
  - Mask: `0x80 >> bitId` (MSB-first: bit `0` is `0x80` = `-128_int8`, bit `7` is `0x01` = `1_int8`).
  - Fortran 2018 `int8` equivalent: byte `bytes(id / 8_id_k)`, bit index `7 - int(mod(id, 8_id_k), i32_k)` via `ibset` / `ibclr` / `btest`.
- **Verified Probe Byte Layout**:
  - Setting bits `{0, 7, 8}` in a 10-bit `vtkBitArray` (`bits = 1000000110`) yields `bytes = [0x81, 0x80]`.
- **Unused Trailing Bits Invariant**:
  - `SetNumberOfValues(vtkIdType number)` calls `InitializeUnusedBitsInLastByte()`, zeroing all trailing bits `(number % 8) .. 7` in the final active byte `(number + 7) / 8 - 1`.
  - Verified by shrinking `10000001` (`8` bits, `0x81`) to `5` bits: byte `0` becomes `0x80` (bit `7` cleared to `0`).

---

## 3. `vtkHyperTree` & `vtkHyperTreeData` (Replaces Legacy `vtkCompactHyperTree`)

### 3.1 Concrete Class & Struct Layout
In VTK `v9.7.2`, `vtkCompactHyperTree` is consolidated directly into concrete `vtkHyperTree` (`vtkStandardNewMacro(vtkHyperTree)`), holding `std::shared_ptr<vtkHyperTreeData> Datas`:

```cpp
struct vtkHyperTreeData
{
  vtkIdType TreeIndex;                        // init: -1
  unsigned int NumberOfLevels;                // init: 1
  vtkIdType NumberOfVertices;                 // init: 1 (coarse + leaves)
  vtkIdType NumberOfNodes;                    // init: 0 (coarse internal nodes)
  vtkIdType GlobalIndexStart;                 // init: -1
  std::vector<unsigned int> ParentToElderChild; // init: size=1, [0]=UINT_MAX
  std::vector<vtkIdType> GlobalIndexTable;    // init: empty (size=0)
};
```

- **`Initialize(unsigned char branchFactor, unsigned char dimension) -> bool`**:
  - Validates `branchFactor` $\in \{2, 3\}$ and `dimension` $\in \{1, 2, 3\}$; returns `false` with warning otherwise.
  - Computes `NumberOfChildren = pow(branchFactor, dimension)` $\in \{2, 3, 4, 8, 9, 27\}$.
  - Note: Default `vtkHyperTree()` constructor calls `Initialize(2, 3)`.
- **`GetNumberOfLeaves()`**: Not stored separately in `vtkHyperTreeData`; computed dynamically as `NumberOfVertices - NumberOfNodes` (`1` on initial 1-leaf root).
- **`NumberOfVerticesPerDepth`**: Not stored inside C++ `vtkHyperTreeData`; emitted on demand by `ComputeBreadthFirstOrderDescriptor` as `vtkTypeInt64Array*`. Maintaining `number_of_vertices_per_depth(0:)` eagerly in Fortran `htg_tree_t` is fast and enables O(1) depth queries and continuous validation.

### 3.2 `SubdivideLeaf(vtkIdType index, unsigned int depth)` & Lazy `ParentToElderChild` Sizing
Verified in `vtkHyperTree.cxx:403-425` and live `/usr/bin/lldb` stepping:

```cpp
assert(index < NumberOfVertices);
assert(IsLeaf(index));
if (static_cast<vtkIdType>(ParentToElderChild.size()) <= index)
{
  ParentToElderChild.resize(index + 1, std::numeric_limits<unsigned int>::max());
}
ParentToElderChild[index] = static_cast<unsigned int>(NumberOfVertices);
if (depth + 1 == NumberOfLevels)
{
  ++NumberOfLevels;
}
NumberOfNodes += 1;
NumberOfVertices += NumberOfChildren;
```

- **Lazy Tail Truncation**: During incremental `SubdivideLeaf`, `ParentToElderChild` only grows up to `max_coarse_index + 1` (minimum size `1` from initialization). Newly created leaf children at `[elder .. elder + k - 1]` are **not** appended to `ParentToElderChild` until one of them (or a higher-indexed leaf) is later subdivided, at which point intermediate leaf slots are filled with `UINT_MAX` (`4294967295u`).
- **LLDB Step-Through Trace (`f=2, d=2`, `k=4`)**:
  1. Initial 1-leaf root: `Levels=1, Vertices=1, Nodes=0, Leaves=1`, `ParentToElderChild` (`size=1`): `[4294967295]`.
  2. `SubdivideLeaf(0, 0)`: `Levels=2, Vertices=5, Nodes=1, Leaves=4`, `ParentToElderChild` (`size=1`): `[1]`.
  3. `SubdivideLeaf(3, 1)`: `Levels=3, Vertices=9, Nodes=2, Leaves=7`, `ParentToElderChild` (`size=4`): `[1, 4294967295, 4294967295, 5]`.
- **Differential Comparison Rule for `@oracle`**:
  - In C++, `const unsigned int* arr = tree->GetElderChildIndexArray(nbElements)` may have `nbElements < NumberOfVertices` (after `SubdivideLeaf`) or `nbElements == max(1, numberOfBits)` (after `BuildFromBreadthFirstOrderDescriptor`).
  - For any local vertex `0 <= v < NumberOfVertices`, C++ treats `v` as a **leaf** iff `v >= nbElements || arr[v] == UINT_MAX || NumberOfVertices == 1`; otherwise `v` is **coarse** with elder child `arr[v]`.
- **`IsTerminalNode` Nuance**:
  - Intended semantics (`vtkHyperTreeGridEntry::IsTerminalNode` and Fortran `tree%is_terminal_node`): returns `true` iff vertex is a **coarse** node (`!IsLeaf`) whose `NumberOfChildren` children are **all leaves**.
  - *C++ `vtkHyperTree::IsTerminalNode(index)` quirk*: `vtkHyperTree::IsTerminalNode` forgets to check `!IsLeaf(index)` when `index < ParentToElderChild.size()` and `ParentToElderChild[index] == UINT_MAX` (computing `UINT_MAX + ichild >= size` $\rightarrow$ `true`!). Always check `!IsLeaf(index)` first in Fortran `is_terminal_node`, and in C++ oracle comparisons guard `tree->IsTerminalNode(v)` with `!tree->IsLeaf(v)`.

---

## 4. Global Indexing (`GlobalIndexStart` vs `GlobalIndexTable`)

- **Initial Unset State**: `GlobalIndexStart = -1`, `GlobalIndexTable.empty() == true`.
- **Implicit Mode (`SetGlobalIndexStart(vtkIdType start)`)**:
  - Sets `Datas->GlobalIndexStart = start`.
  - `GetGlobalIndexFromLocal(index)` (when `GlobalIndexTable.empty()`): asserts `GlobalIndexStart >= 0` and `index >= 0`, returns `GlobalIndexStart + index`.
  - `GetGlobalNodeIndexMax()`: returns `GlobalIndexStart + NumberOfVertices - 1`.
- **Explicit Mode (`SetGlobalIndexFromLocal(vtkIdType index, vtkIdType global)`)**:
  - Asserts `GlobalIndexStart < 0` (cannot mix implicit `GlobalIndexStart` and explicit `GlobalIndexTable`).
  - On `vtkHyperTree`, takes 2 arguments `(local_index, global_index)`; on `vtkHyperTreeGridNonOrientedCursor`, `SetGlobalIndexFromLocal(global_index)` takes **1 argument** (for the cursor's current vertex).
  - Resizes `GlobalIndexTable` to `index + 1` initialized with `-1` if `index >= GlobalIndexTable.size()`, and sets `GlobalIndexTable[index] = global` (overwriting prior explicit values is allowed).
  - `GetGlobalIndexFromLocal(index)`: asserts `0 <= index < GlobalIndexTable.size()` and `GlobalIndexTable[index] >= 0`, returns `GlobalIndexTable[index]`.
  - `GetGlobalNodeIndexMax()`: returns `*std::max_element(GlobalIndexTable.begin(), GlobalIndexTable.end())`.
- **VTK `v9.7.2` Quirk Alert (`IsGlobalIndexImplicit`)**:
  - `vtkHyperTree::IsGlobalIndexImplicit()` in `vtkHyperTree.cxx:319` is implemented as `return this->Datas->GlobalIndexStart == -1;` (inverted!). Do **not** call `tree->IsGlobalIndexImplicit()` in C++ oracle tests; check `tree->GetGlobalIndexStart() >= 0` instead.

---

## 5. Breadth-First Order Descriptors

### 5.1 Exact VTK `v9.7.2` Signatures
```cpp
virtual void BuildFromBreadthFirstOrderDescriptor(
  vtkBitArray* descriptor, vtkIdType numberOfBits, vtkIdType startIndex = 0);

virtual void ComputeBreadthFirstOrderDescriptor(
  unsigned int depthLimiter,
  vtkBitArray* inputMask,
  vtkTypeInt64Array* numberOfVerticesPerDepth,
  vtkBitArray* descriptor,
  vtkIdList* breadthFirstIdMap);
```

### 5.2 Exact Encoding & Edge-Case Semantics
1. **Last Depth Omitted & Single-Leaf Root (`numberOfBits = 0`)**:
   - `ComputeBreadthFirstOrderDescriptor` writes bits only for depths `0 .. maxDepth - 2` (`for (int idepth = 0; idepth < maxDepth - 1; ++idepth)`), because all vertices at the deepest active level `maxDepth - 1` have no children in the output topology.
   - Therefore, for a **1-level single-leaf root tree** (`NumberOfLevels = 1`):
     - `ComputeBreadthFirstOrderDescriptor` appends **`0` bits** to `descriptor`, `[1]` to `numberOfVerticesPerDepth`, and `[root_global_id]` to `breadthFirstIdMap`.
     - `BuildFromBreadthFirstOrderDescriptor(desc, 0, startIndex)` (`numberOfBits == 0`) cleanly constructs the 1-level single-leaf root (`Levels = 1, Vertices = 1, Nodes = 0, ParentToElderChild = [UINT_MAX]`).
     - *Warning*: Passing `numberOfBits = 1` with bit `[0]` to C++ `BuildFromBreadthFirstOrderDescriptor` causes VTK's depth counter loop (`if (++currentPositionAtDepth == currentDepthSize) ++numberOfDepths;`) to set `NumberOfLevels = 2` despite `NumberOfVertices = 1`! Canonical single-leaf root descriptor length in VTK `v9.7.2` is `numberOfBits = 0`.
2. **Coarse (`1`) vs Leaf (`0`) Bits & Appending**:
   - `1` (`true`) = coarse vertex (spawns `NumberOfChildren` vertices at `depth + 1`); `0` (`false`) = leaf vertex.
   - All three output containers (`numberOfVerticesPerDepth`, `descriptor`, `breadthFirstIdMap`) **append** to existing contents so forest-wide arrays can be concatenated across trees.
3. **`breadthFirstIdMap` Stores Global IDs**:
   - Maps breadth-first traversal slot $\rightarrow$ `idg = GetGlobalIndexFromLocal(local_index)` (**global** node index, requires global indexing to be set on the tree).
   - Verified on out-of-order bintree (`GlobalIndexStart = 100`, `ParentToElderChild = [1, 5, 3, 7]`): `vpd = [1, 2, 4, 2]`, `descriptor = 1110010` (`7` bits), `breadthFirstIdMap = [100, 101, 102, 105, 106, 103, 104, 107, 108]`.
4. **`inputMask` Pruning**:
   - Recursion rule (`vtkHyperTree.cxx:485-503`):
     - `idg = GetGlobalIndexFromLocal(index)` is **always** recorded in `breadthFirstOrderIdMapPerDepth[depth]` when `index` is visited (because its unmasked parent requires all `NumberOfChildren` child slots!).
     - Subdivision continues iff `!IsLeaf(index) && !(inputMask && inputMask->GetValue(idg)) && depth < depthLimiter`.
     - Thus, if a coarse vertex `v` is masked (`inputMask[idg] == 1`), `v` itself appears at its depth with descriptor bit `0`, and **its entire descendant subtree is omitted**.
     - Trailing empty depths (where all subtrees were masked out) are squeezed off `maxDepth` so `numberOfVerticesPerDepth` never ends with `0` entries.
5. **`depthLimiter`**:
   - Stops recursion once `depth == depthLimiter` (emitting `0` at `depthLimiter - 1` if `depthLimiter == 0`, or omitting depth `depthLimiter`'s all-zero bits as the new deepest level). For `depthLimiter = 1` on the bintree above: `vpd = [1, 2]`, `idMap = [100, 101, 102]`, `descriptor = [1]` (`1` bit).
6. **`BuildFromBreadthFirstOrderDescriptor` Layout**:
   - Assigns local vertex IDs strictly in breadth-first order (`0`, then `1 .. k`, etc.).
   - Populates `ParentToElderChild` to exact size `max(1, numberOfBits)` (`elder_id` for `1` bits, `UINT_MAX` for `0` bits).

---

## 6. `vtkHyperTreeGrid`, Cursor, `Mask`, `PureMask`, and Material Interfaces

### 6.1 Initialization & Deprecated Methods
- **Default Constructor & `Initialize()`**:
  - **CRITICAL**: Both `vtkNew<vtkHyperTreeGrid>` and `htg->Initialize()` leave **`BranchFactor = 0`**! Callers **must** call `SetBranchFactor(2)` or `SetBranchFactor(3)` (and `SetDimensions(...)`) after `Initialize()` before creating any trees via `NewNonOrientedCursor(treeId, true)`.
  - Defaults: `TransposedRootIndexing = false`, `DepthLimiter = UINT_MAX`, `HasInterface = false`, `InterfaceNormalsName = nullptr`, `InterfaceInterceptsName = nullptr`, `Mask = nullptr`, `PureMask = nullptr`.
- **`ModeSqueeze`**: `SetModeSqueeze` / `GetModeSqueeze` are marked `VTK_DEPRECATED_IN_9_6_0("No effect anymore, do not use.")` (`GetModeSqueeze()` returns `nullptr`).

### 6.2 `SetDimensions` Squeeze Rules
Given point counts `(dx, dy, dz)` with `dx >= 1, dy >= 1, dz >= 1` (`1` indicates a squeezed/degenerate axis; do **not** pass `0`):
- `Dimensions = [dx, dy, dz]` (`unsigned int[3]`)
- `CellDims[a] = (Dimensions[a] == 1) ? 1 : (Dimensions[a] - 1)` (`unsigned int[3]`)
- `MaxNumberOfTrees = CellDims[0] * CellDims[1] * CellDims[2]`
- Exact `Dimension`, `Orientation`, and `Axis[2]` (`unsigned int[2]`) verified via `probe` and `lldb`:

| Grid Configuration | `Dimensions` | `CellDims` | `Dimension` | `Orientation` | `Axis[0]` | `Axis[1]` |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: |
| 0D degenerate | `(1, 1, 1)` | `(1, 1, 1)` | `0` | `UINT_MAX` | `UINT_MAX` | `UINT_MAX` |
| 1D along X | `(nx, 1, 1)` | `(nx-1, 1, 1)` | `1` | `0` | `0` | `UINT_MAX` |
| 1D along Y | `(1, ny, 1)` | `(1, ny-1, 1)` | `1` | `1` | `1` | `UINT_MAX` |
| 1D along Z | `(1, 1, nz)` | `(1, 1, nz-1)` | `1` | `2` | `2` | `UINT_MAX` |
| 2D in YZ plane | `(1, ny, nz)` | `(1, ny-1, nz-1)` | `2` | `0` (normal X) | `1` | `2` |
| 2D in XZ plane | `(nx, 1, nz)` | `(nx-1, 1, nz-1)` | `2` | `1` (normal Y) | **`2`** | **`0`** |
| 2D in XY plane | `(nx, ny, 1)` | `(nx-1, ny-1, 1)` | `2` | `2` (normal Z) | `0` | `1` |
| 3D volume | `(nx, ny, nz)` | `(nx-1, ny-1, nz-1)` | `3` | `UINT_MAX` | `UINT_MAX` | `UINT_MAX` |

> [!IMPORTANT]
> Notice that for **2D XZ (`Dimensions = (nx, 1, nz)`)**, `Orientation = 1` and **`Axis = (2, 0)`** (cyclic right-handed `(Orientation + 1) % 3`, `(Orientation + 2) % 3`), **not** `(0, 2)`!

- **Coordinates (`SetXCoordinates`, `SetYCoordinates`, `SetZCoordinates`)**:
  - VTK accepts monotonic coordinate arrays (both ascending and descending coordinates were verified without error).

### 6.3 Root Level-Zero Indexing Formulas
Let `(ni, nj, nk) = (CellDims[0], CellDims[1], CellDims[2])` and `0 <= i < ni`, `0 <= j < nj`, `0 <= k < nk`:
- **`TransposedRootIndexing == false`** (`SetIndexingModeToKJI()`, default — `i` varies fastest):
  $$\text{tree\_id} = (k \cdot n_j + j) \cdot n_i + i = i + j \cdot n_i + k \cdot n_i \cdot n_j$$
  $$i = \text{tree\_id} \bmod n_i, \qquad j = (\text{tree\_id} / n_i) \bmod n_j, \qquad k = \text{tree\_id} / (n_i \cdot n_j)$$
- **`TransposedRootIndexing == true`** (`SetIndexingModeToIJK()` — `k` varies fastest):
  $$\text{tree\_id} = (i \cdot n_j + j) \cdot n_k + k = k + j \cdot n_k + i \cdot n_k \cdot n_j$$
  $$k = \text{tree\_id} \bmod n_k, \qquad j = (\text{tree\_id} / n_k) \bmod n_j, \qquad i = \text{tree\_id} / (n_k \cdot n_j)$$
- Both formulas hold uniformly across 0D, 1D, 2D, and 3D because squeezed axes have `CellDims[a] == 1` and coordinate `0`.

### 6.4 Grid-Wide Aggregations
- **Total vertices method name**: `vtkHyperTreeGrid::GetNumberOfCells()` (there is **no** `GetNumberOfVertices()` on `vtkHyperTreeGrid` in VTK `v9.7.2`!). Sums `tree->GetNumberOfVertices()` across non-empty trees (`0` when empty).
- `GetNumberOfLeaves()`: sums `tree->GetNumberOfLeaves()` across non-empty trees (`0` when empty).
- `GetNumberOfLevels()`: max `tree->GetNumberOfLevels()` across non-empty trees (`0` when empty).
- `GetGlobalNodeIndexMax()`: max `tree->GetGlobalNodeIndexMax()` across non-empty trees (`0` when empty).

### 6.5 `Mask`, `IsMasked`, `ToChild`, and `SubdivideLeaf` Preconditions
- **`IsMasked` (`vtkHyperTreeGridEntry::IsMasked(grid, tree)`)**:
  - Short-circuits to `false` immediately if `tree == nullptr` or `!grid->HasMask()` (**does not** query `GetGlobalNodeIndex()` when `HasMask() == false`, allowing `SubdivideLeaf()` before `SetGlobalIndexStart()`!).
  - When `grid->HasMask()` is `true`: computes `gid = GetGlobalNodeIndex(tree)` and safely bounds-checks: returns `(gid < Mask->GetNumberOfValues()) ? (Mask->GetValue(gid) != 0) : false`.
- **`ToChild` Precondition on Masked Coarse Nodes (`vtkHyperTreeGridEntry.cxx:126`)**:
  - `vtkHyperTreeGridEntry::ToChild` executes `assert("pre: is_masked" && !this->IsMasked(grid, tree))` (in addition to `assert(!this->IsLeaf(grid, tree, level))` and `ichild < tree->GetNumberOfChildren()`). Calling `cursor->ToChild(ichild)` on a **masked coarse node** is illegal in VTK `v9.7.2` Debug builds.
- **`SubdivideLeaf` Behavior on Masked or Already-Coarse Cells**:
  1. **Masked Leaf**: `vtkHyperTreeGridEntry::SubdivideLeaf` executes `assert("pre: is_masked" && !this->IsMasked(grid, tree))` — subdividing a masked cell is **illegal** (aborts in C++ Debug; must fail with `HTG_ERR_MASK` in Fortran).
  2. **Already-Coarse Node**:
     - On `vtkHyperTreeGridNonOrientedCursor::SubdivideLeaf()` (`Entry::SubdivideLeaf`), `if (this->IsLeaf(grid, tree, level))` wraps the call — calling `cursor->SubdivideLeaf()` on an already-coarse node is a **silent no-op** in C++ VTK `v9.7.2`.
     - Direct `vtkHyperTree::SubdivideLeaf(index, depth)` executes `assert("pre: not_leaf" && this->IsLeaf(index))` — subdividing a coarse vertex directly on `vtkHyperTree` is **illegal**.

### 6.6 `GetPureMask()` Exact Bottom-Up Logic & Material Interface Interaction
Verified across `probe`, `probe2`, `probe3`, `vtkHyperTreeGrid.cxx:1172-1285`, and `lldb`:
- **Bit Polarity**: `PureMask` is a mask where **`0` = PURE** and **`1` = IMPURE (masked or mixed/interfaced)**.
- **Allocation Size**: `numBits = GetGlobalNodeIndexMax() + 1` (initialized to all `0`, including unreferenced global ID gaps in sparse/explicit forests; `1` bit `[0]` on an empty grid).
- **Material Interface Presence & 3-Component Check (`vtkHyperTreeGrid.cxx:1270, 1277`)**:
  - Interface impurity checking is enabled iff `HasInterface == true` **and** `InterfaceNormalsName` names an existing array in `CellData` with **`GetNumberOfComponents() == 3`** **and** `InterfaceInterceptsName` names an existing array (`intercepts`) in `CellData` with **`GetNumberOfComponents() == 3`**.
  - Checked **only on unmasked leaf vertices**: leaf `gid` has an active interface iff **`intercepts->GetTuple(gid)[2] < 2.0`** (in VTK HTG convention, component `2` is `-1`, `0`, or `1` for interfaced cells and `2` for pure cells; coarse `intercepts` values are ignored).
- **Exact Recursion Order (`RecursivelyInitializePureMask`, `vtkHyperTreeGrid.cxx:1181-1208`)**:
  1. **Masked Short-Circuit First (Both Coarse & Leaf)**:
     Checks `bool mask = cursor->IsMasked();` **before** `IsLeaf()` or `ToChild()`. If `mask == true`, immediately sets `PureMask[gid] = true` (`1`) and returns `true` **without** descending into children! (This avoids hitting `ToChild`'s `!IsMasked` assertion on masked coarse nodes and leaves any unvisited descendant `PureMask` bits at their initial `0` value.)
  2. **Unmasked Leaf Vertex (`gid`)**:
     Computes `impure = hasValidInterface && (Intercepts[gid][2] < 2.0)`, sets `PureMask[gid] = impure`, and returns `impure`.
  3. **Unmasked Coarse Vertex (`gid`)**:
     Starts with `impure = false`, visits **all** `NumberOfChildren` children `c` setting `impure |= recurse(c)` (bitwise `|=` without short-circuiting, so every reachable unmasked-parent descendant's `PureMask` bit is populated), sets `PureMask[gid] = impure`, and returns `impure`.
- **Caching Caveat for `@oracle`**:
  - C++ `vtkHyperTreeGrid` caches `vtkBitArray* PureMask` until `CleanPureMask()` is called (`SetMask(differentPtr)`, `SetHasInterface`, `SetInterfaceNormalsName`, `SetInterfaceInterceptsName`, `Initialize`). Calling `cursor->SetMask(state)` mutates `this->Mask` in-place via `InsertValue` **without** calling `CleanPureMask()`. In C++ oracle tests, finish `cursor->SetMask` mutations before calling `GetPureMask()` (or reset `SetMask`). In Fortran `fvtkhtg`, invalidate `pure_mask_valid` on every topology, mask, or interface mutation.

---

## 7. Cell Geometry, Face Neighbours, and Coarsening

### 7.1 Geometry Cursor Arithmetic (`vtkHyperTreeGridNonOrientedGeometryCursor`)
Verified in `vtkHyperTreeGridGeometryEntry.{h,cxx}`, `vtkHyperTreeGridScales.h`, `vtkHyperTreeGrid.cxx:1483-1540`, and by `test_oracle_geometry_neighbors` (bit-exact over 5484 cells in 128 configurations):
- **Root cell** (`GetLevelZeroOriginAndSizeFromIndex`): `Origin = (X[i], Y[j], Z[k])`; `Size[a] = Coord[a][i+1] - Origin[a]`, and exactly `0.` when `Dimensions[a] == 1`. Descending coordinates give a negative size; nothing is reordered.
- **Scales** (`vtkHyperTreeGridScales::Update`): `scale(level) = scale(level - 1) / BranchFactor`, computed iteratively as `double` division for all three components.
- **Descent** (`vtkHyperTreeGridGeometryEntry::ToChild`): `Origin[axis] += digit * sizeChild[axis]` for every active axis, where `sizeChild = scale(level + 1)` and `digit` is the base-`BranchFactor` digit of `ichild` for that axis. The addition is always performed, including `0 * size`.
- **Queries**: `GetBounds` = `(Origin[a], Origin[a] + size[a])`; `GetPoint` = `Origin[a] + size[a] / 2.`.
- `fvtkhtg` recomputes this from the cursor's path on demand with the same operations in the same order, so results are bit-identical without storing an origin per stack entry.

### 7.2 Child Numbering Uses Ascending Axes, Not `GetAxes()`
- `vtkHyperTreeGridGeometryEntry::ToChild` (and `vtkHyperTreeGridGeometryLevelEntry::ToChild`) hard-code the in-plane axes from the orientation: YZ plane `(1, 2)`, **XZ plane `(0, 2)`**, XY plane `(0, 1)`. Digit 0 of a child index is always the lower-numbered axis.
- `vtkHyperTreeGrid::GetAxes()` for the XZ plane is `(2, 0)` (`vtkHyperTreeGrid.cxx:575-579` swaps the axes "to have a direct frame"). This is grid metadata; it does **not** describe child numbering.
- Probe on VTK `9.7.2`, `SetDimensions(3, 1, 3)`, root refined: child 1 has origin offset in **x**, child 2 in **z**.
- `fvtkhtg`: `grid%get_axes()` mirrors `GetAxes()`; `grid%get_active_axes()` and `grid%get_child_position()` give the ascending order that child indices actually use.

### 7.3 Von Neumann Supercursor Semantics (`vtkHyperTreeGridNonOrientedVonNeumannSuperCursor`)
- **What a face slot holds** (`vtkHyperTreeGridNonOrientedSuperCursor::ToChild`, lines 353-415): the cell at the neighbouring position if it exists at the central cell's level; otherwise the coarser leaf that covers it. A neighbour entry is **not** descended when it has no tree, is a leaf, or is **masked** (`Mask->GetValue(globalIndex)`), so a masked cell stops the descent and is itself reported.
- **Root level** (`...VonNeumannSuperCursor::Initialize`, lines 150-300): neighbour trees come from `GetShiftedLevelZeroIndex`, honouring `TransposedRootIndexing`; a slot has no tree at the edge of the grid or where the neighbouring root cell has none.
- `cursor%get_face_neighbor` reproduces this per face. `test_oracle_geometry_neighbors` identifies each VTK slot by geometry (which face of the cell it touches) and compares tree, local id, global id, level, leaf and mask state, and bounds: 18761 neighbours matched (3411 coarser, 2110 across root cells), and 3389 faces had no neighbour in both.

### 7.4 VTK `v9.7.2` Inconsistency on XZ-Plane Grids with Several Root Cells
- For `Dimension == 2` the supercursor takes its root-level "W/E" neighbours along `GetAxes()[0]` and "S/N" along `GetAxes()[1]` (`...VonNeumannSuperCursor.cxx:185-232`, `GetShiftedLevelZeroIndex` at `vtkHyperTreeGrid.cxx:1414-1428`). Inside a tree, the traversal tables tie "W/E" to digit 0 of the child index and "S/N" to digit 1.
- For the XZ plane `GetAxes()` is `(z, x)` but digit 0 is `x` (Section 7.2). Across a root-cell boundary the supercursor therefore holds the neighbour along the **other** axis: a cell that does not touch the central cell.
- Measured with VTK's own bounds on `SetDimensions(4, 1, 3)` grids: between 0 and 48 of the cross-tree neighbour slots per configuration do not touch the central cell (for example 39 of 58 for branch factor 2 with a mask). Single-root XZ grids, and all XY, YZ, 1D and 3D grids, show no such slot.
- `fvtkhtg` defines the direction geometrically, so `get_face_neighbor` always returns the cell that touches the requested face. On XZ grids with several root cells the oracle test checks geometric adjacency instead of comparing with VTK across root boundaries, and prints the count of non-touching VTK slots.

### 7.5 Coarsening Has No VTK Counterpart
- `vtkHyperTree`, `vtkHyperTreeGrid` and the cursors expose no operation that removes children; VTK trees only grow (or are rebuilt from a descriptor).
- `fvtkhtg` adds `tree%coarsen_node` / `cursor%coarsen_node`. To keep VTK's invariants (contiguous local ids, children in blocks of `NumberOfChildren` starting at `1 + m * NumberOfChildren`), the block with the highest ids is relocated into the freed block. Children may then have lower ids than their parent, which VTK's own construction never produces but none of its data layout forbids.
- A per-block parent table (`block_parent`, absent in VTK) makes the relocation and `get_parent_index` O(1); `htg_validate_tree` checks it against `parent_to_elder_child`.
- Verified by `test_oracle_geometry_neighbors` (729 cells coarsened back: structure, geometry, vertex and level counts and breadth-first descriptors equal a VTK grid that was never refined further) and by `test_coarsen` (random refine/coarsen with per-cell data in both indexing modes).
