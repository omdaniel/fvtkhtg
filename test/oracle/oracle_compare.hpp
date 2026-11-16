#ifndef FVTKHTG_ORACLE_COMPARE_HPP
#define FVTKHTG_ORACLE_COMPARE_HPP

#include "fvtkhtg_test_c_bridge.h"

#include <vtkBitArray.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkHyperTree.h>
#include <vtkHyperTreeGrid.h>
#include <vtkHyperTreeGridNonOrientedCursor.h>
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkTypeInt64Array.h>
#include <vtkVersion.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#if (VTK_MAJOR_VERSION < 9) || (VTK_MAJOR_VERSION == 9 && VTK_MINOR_VERSION < 7)
#error "fvtkhtg oracle tests strictly require VTK >= 9.7.0!"
#endif

namespace oracle {

inline void fail_msg(const char* file, int line, const std::string& detail)
{
  std::cerr << "\n[ORACLE MISMATCH] " << file << ":" << line << "\n  " << detail << std::endl;
  std::exit(1);
}

#define ORACLE_ASSERT(cond, msg_expr)                                                              \
  do                                                                                               \
  {                                                                                                \
    if (!(cond))                                                                                   \
    {                                                                                              \
      std::ostringstream _oss;                                                                     \
      _oss << msg_expr;                                                                            \
      ::oracle::fail_msg(__FILE__, __LINE__, _oss.str());                                          \
    }                                                                                              \
  } while (0)

#define ORACLE_EQ(lhs, rhs, ctx)                                                                   \
  do                                                                                               \
  {                                                                                                \
    auto _l = (lhs);                                                                               \
    auto _r = (rhs);                                                                               \
    if (_l != _r)                                                                                  \
    {                                                                                              \
      std::ostringstream _oss;                                                                     \
      _oss << ctx << " | expected (VTK) = " << _l << ", actual (Fortran) = " << _r;                \
      ::oracle::fail_msg(__FILE__, __LINE__, _oss.str());                                          \
    }                                                                                              \
  } while (0)

inline void verify_vtk_oracle_version()
{
  const int major = vtkVersion::GetVTKMajorVersion();
  const int minor = vtkVersion::GetVTKMinorVersion();
  const char* ver = vtkVersion::GetVTKVersion();
  ORACLE_ASSERT(major > 9 || (major == 9 && minor >= 7),
    "VTK runtime version must be >= 9.7.0, got: " << (ver ? ver : "null"));
  std::cout << "[Oracle] Verified VTK runtime version: " << ver << std::endl;
}

inline void set_rectilinear_coords_both(vtkHyperTreeGrid* vtk_grid, void* ftn_grid,
  const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z)
{
  vtkNew<vtkDoubleArray> vx, vy, vz;
  vx->SetNumberOfValues(static_cast<vtkIdType>(x.size()));
  for (size_t i = 0; i < x.size(); ++i)
  {
    vx->SetValue(static_cast<vtkIdType>(i), x[i]);
  }
  vy->SetNumberOfValues(static_cast<vtkIdType>(y.size()));
  for (size_t i = 0; i < y.size(); ++i)
  {
    vy->SetValue(static_cast<vtkIdType>(i), y[i]);
  }
  vz->SetNumberOfValues(static_cast<vtkIdType>(z.size()));
  for (size_t i = 0; i < z.size(); ++i)
  {
    vz->SetValue(static_cast<vtkIdType>(i), z[i]);
  }
  vtk_grid->SetXCoordinates(vx);
  vtk_grid->SetYCoordinates(vy);
  vtk_grid->SetZCoordinates(vz);

  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_set_x_coordinates(ftn_grid, x.data(), static_cast<int64_t>(x.size())),
    "ftn_htg_grid_set_x_coordinates");
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_set_y_coordinates(ftn_grid, y.data(), static_cast<int64_t>(y.size())),
    "ftn_htg_grid_set_y_coordinates");
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_set_z_coordinates(ftn_grid, z.data(), static_cast<int64_t>(z.size())),
    "ftn_htg_grid_set_z_coordinates");
}

inline void compute_vtk_vpd_recursive(
  const vtkHyperTree* tree, vtkIdType v, size_t depth, std::vector<int64_t>& vpd)
{
  if (depth >= vpd.size())
  {
    vpd.resize(depth + 1, 0);
  }
  vpd[depth] += 1;
  if (!tree->IsLeaf(v))
  {
    const vtkIdType elder = tree->GetElderChildIndex(static_cast<unsigned int>(v));
    const vtkIdType nch = tree->GetNumberOfChildren();
    for (vtkIdType c = 0; c < nch; ++c)
    {
      compute_vtk_vpd_recursive(tree, elder + c, depth + 1, vpd);
    }
  }
}

inline void compare_cursors_recursive(vtkHyperTreeGrid* vtk_grid,
  vtkHyperTreeGridNonOrientedCursor* vtk_cur, void* ftn_cur, bool check_global_idx,
  bool check_masks, const std::string& label)
{
  int32_t st = 0;
  ORACLE_EQ(static_cast<int64_t>(vtk_cur->GetVertexId()), ftn_htg_cursor_get_vertex_id(ftn_cur),
    label << " cursor GetVertexId");
  ORACLE_EQ(static_cast<int32_t>(vtk_cur->GetLevel()), ftn_htg_cursor_get_level(ftn_cur),
    label << " cursor GetLevel");
  ORACLE_EQ(vtk_cur->IsRoot(), ftn_htg_cursor_is_root(ftn_cur), label << " cursor IsRoot");

  const bool vtk_leaf = vtk_cur->IsLeaf();
  const bool ftn_leaf = ftn_htg_cursor_is_leaf(ftn_cur, &st);
  ORACLE_EQ(FTN_HTG_OK, st, label << " cursor is_leaf status");
  ORACLE_EQ(vtk_leaf, ftn_leaf, label << " cursor IsLeaf at v=" << vtk_cur->GetVertexId());

  const bool vtk_term = !vtk_leaf && vtk_cur->GetTree()->IsTerminalNode(vtk_cur->GetVertexId());
  const bool ftn_term = ftn_htg_cursor_is_terminal_node(ftn_cur, &st);
  ORACLE_EQ(FTN_HTG_OK, st, label << " cursor is_terminal_node status");
  ORACLE_EQ(
    vtk_term, ftn_term, label << " cursor IsTerminalNode at v=" << vtk_cur->GetVertexId());

  if (check_global_idx)
  {
    const int64_t vtk_gid = static_cast<int64_t>(vtk_cur->GetGlobalNodeIndex());
    const int64_t ftn_gid = ftn_htg_cursor_get_global_node_index(ftn_cur, &st);
    ORACLE_EQ(FTN_HTG_OK, st, label << " cursor get_global_node_index status");
    ORACLE_EQ(
      vtk_gid, ftn_gid, label << " cursor GetGlobalNodeIndex at v=" << vtk_cur->GetVertexId());
  }

  if (check_masks)
  {
    const bool vtk_masked = vtk_cur->IsMasked();
    const bool ftn_masked = ftn_htg_cursor_is_masked(ftn_cur, &st);
    ORACLE_EQ(FTN_HTG_OK, st, label << " cursor is_masked status");
    ORACLE_EQ(
      vtk_masked, ftn_masked, label << " cursor IsMasked at v=" << vtk_cur->GetVertexId());
  }

  // Note: In VTK Debug builds, `vtkHyperTreeGridEntry::ToChild` asserts `!IsMasked()`.
  // So when walking via cursor, only descend into unmasked coarse nodes (matching VTK contract);
  // direct tree-level checks below still verify ALL vertices even beneath masked coarse nodes!
  const bool can_descend = !vtk_leaf && !vtk_cur->IsMasked();
  if (can_descend)
  {
    const int32_t nch = static_cast<int32_t>(vtk_cur->GetNumberOfChildren());
    ORACLE_EQ(nch, ftn_htg_cursor_get_number_of_children(ftn_cur), label << " cursor nchildren");
    for (int32_t c = 0; c < nch; ++c)
    {
      vtk_cur->ToChild(static_cast<unsigned char>(c));
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(ftn_cur, c), label << " to_child");
      compare_cursors_recursive(
        vtk_grid, vtk_cur, ftn_cur, check_global_idx, check_masks, label);
      vtk_cur->ToParent();
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_parent(ftn_cur), label << " to_parent");
    }
  }
}

inline std::string bytes_to_hex(const uint8_t* ptr, int64_t n)
{
  std::ostringstream oss;
  for (int64_t i = 0; i < n; ++i)
  {
    if (i > 0)
      oss << ' ';
    oss << std::hex << std::setw(2) << std::setfill('0') << static_cast<unsigned int>(ptr[i]);
  }
  return oss.str();
}

struct CompareOptions
{
  bool check_coordinates = true;
  bool check_global_indexing = false;
  bool check_masks = false;
  bool check_pure_mask = false;
  bool check_validate = false;
};

inline void assert_grids_equivalent(vtkHyperTreeGrid* vtk_grid, void* ftn_grid,
  const std::string& label, const CompareOptions& opts = CompareOptions{})
{
  int vtk_dims[3] = { 0, 0, 0 };
  int ftn_dims[3] = { 0, 0, 0 };
  vtk_grid->GetDimensions(vtk_dims);
  ftn_htg_grid_get_dimensions(ftn_grid, ftn_dims);
  for (int a = 0; a < 3; ++a)
  {
    ORACLE_EQ(vtk_dims[a], ftn_dims[a], label << " Dimensions[" << a << "]");
  }

  int vtk_cdims[3] = { 0, 0, 0 };
  int ftn_cdims[3] = { 0, 0, 0 };
  vtk_grid->GetCellDims(vtk_cdims);
  ftn_htg_grid_get_cell_dims(ftn_grid, ftn_cdims);
  for (int a = 0; a < 3; ++a)
  {
    ORACLE_EQ(vtk_cdims[a], ftn_cdims[a], label << " CellDims[" << a << "]");
  }

  const int32_t dim = static_cast<int32_t>(vtk_grid->GetDimension());
  ORACLE_EQ(dim, ftn_htg_grid_get_dimension(ftn_grid), label << " Dimension");
  ORACLE_EQ(static_cast<int32_t>(vtk_grid->GetBranchFactor()),
    ftn_htg_grid_get_branch_factor(ftn_grid), label << " BranchFactor");
  ORACLE_EQ(static_cast<int32_t>(vtk_grid->GetNumberOfChildren()),
    ftn_htg_grid_get_number_of_children(ftn_grid), label << " NumberOfChildren");
  ORACLE_EQ(vtk_grid->GetTransposedRootIndexing(),
    ftn_htg_grid_get_transposed_root_indexing(ftn_grid), label << " TransposedRootIndexing");

  // Orientation & Axes in 1D / 2D
  if (dim == 1 || dim == 2)
  {
    ORACLE_EQ(static_cast<int32_t>(vtk_grid->GetOrientation()),
      ftn_htg_grid_get_orientation(ftn_grid), label << " Orientation");
    const unsigned int* vtk_axes = vtk_grid->GetAxes();
    int32_t ftn_axes[3] = { -1, -1, -1 };
    int32_t naxes = 0;
    ftn_htg_grid_get_axes(ftn_grid, ftn_axes, &naxes);
    for (int32_t a = 0; a < dim; ++a)
    {
      ORACLE_EQ(
        static_cast<int32_t>(vtk_axes[a]), ftn_axes[a], label << " Axis[" << a << "]");
    }
  }

  // Coordinates zero-copy views
  if (opts.check_coordinates)
  {
    vtkDataArray* vcoords[3] = { vtk_grid->GetXCoordinates(), vtk_grid->GetYCoordinates(),
      vtk_grid->GetZCoordinates() };
    const double* fptr[3] = { nullptr, nullptr, nullptr };
    int64_t flen[3] = { 0, 0, 0 };
    ftn_htg_grid_get_x_coordinates_view(ftn_grid, &fptr[0], &flen[0]);
    ftn_htg_grid_get_y_coordinates_view(ftn_grid, &fptr[1], &flen[1]);
    ftn_htg_grid_get_z_coordinates_view(ftn_grid, &fptr[2], &flen[2]);
    for (int a = 0; a < 3; ++a)
    {
      ORACLE_ASSERT(vcoords[a] != nullptr, label << " VTK coord array null on axis " << a);
      ORACLE_EQ(static_cast<int64_t>(vcoords[a]->GetNumberOfTuples()), flen[a],
        label << " Coord length axis " << a);
      for (int64_t idx = 0; idx < flen[a]; ++idx)
      {
        const double vval = vcoords[a]->GetTuple1(static_cast<vtkIdType>(idx));
        const double fval = fptr[a][idx];
        ORACLE_ASSERT(std::abs(vval - fval) <= 1e-15,
          label << " Coord[" << a << "][" << idx << "] VTK=" << vval << " FTN=" << fval);
      }
    }
  }

  // Forest-level counts
  ORACLE_EQ(static_cast<int64_t>(vtk_grid->GetMaxNumberOfTrees()),
    ftn_htg_grid_get_max_number_of_trees(ftn_grid), label << " MaxNumberOfTrees");
  ORACLE_EQ(static_cast<int64_t>(vtk_grid->GetNumberOfNonEmptyTrees()),
    ftn_htg_grid_get_number_of_non_empty_trees(ftn_grid), label << " NumberOfNonEmptyTrees");
  ORACLE_EQ(static_cast<int64_t>(vtk_grid->GetNumberOfCells()),
    ftn_htg_grid_get_number_of_vertices(ftn_grid), label << " NumberOfVertices/Cells");
  ORACLE_EQ(static_cast<int64_t>(vtk_grid->GetNumberOfLeaves()),
    ftn_htg_grid_get_number_of_leaves(ftn_grid), label << " NumberOfLeaves");
  ORACLE_EQ(static_cast<int64_t>(vtk_grid->GetNumberOfLevels()),
    ftn_htg_grid_get_number_of_levels(ftn_grid), label << " Grid NumberOfLevels");

  // Deterministic ascending tree iteration & per-tree deep comparison
  vtkHyperTreeGrid::vtkHyperTreeGridIterator it;
  vtk_grid->InitializeTreeIterator(it);
  vtkIdType vtk_tree_id = -1;
  int64_t pos = 0;
  while (vtkHyperTree* vtk_tree = it.GetNextTree(vtk_tree_id))
  {
    int64_t ftn_tree_id = -1;
    void* ftn_tree = nullptr;
    ORACLE_EQ(FTN_HTG_OK,
      ftn_htg_grid_get_tree_by_pos(ftn_grid, pos, &ftn_tree_id, &ftn_tree),
      label << " get_tree_by_pos(" << pos << ")");
    ORACLE_EQ(static_cast<int64_t>(vtk_tree_id), ftn_tree_id,
      label << " active tree_id at pos=" << pos);
    ORACLE_ASSERT(ftn_tree != nullptr, label << " ftn_tree null at pos=" << pos);

    const std::string tlabel = label + " [tree " + std::to_string(ftn_tree_id) + "]";
    ORACLE_EQ(static_cast<int64_t>(vtk_tree->GetTreeIndex()),
      ftn_htg_tree_get_tree_index(ftn_tree), tlabel << " GetTreeIndex");
    ORACLE_EQ(static_cast<int32_t>(vtk_tree->GetBranchFactor()),
      ftn_htg_tree_get_branch_factor(ftn_tree), tlabel << " GetBranchFactor");
    ORACLE_EQ(static_cast<int32_t>(vtk_tree->GetDimension()),
      ftn_htg_tree_get_dimension(ftn_tree), tlabel << " GetDimension");
    ORACLE_EQ(static_cast<int32_t>(vtk_tree->GetNumberOfChildren()),
      ftn_htg_tree_get_number_of_children(ftn_tree), tlabel << " GetNumberOfChildren");

    const int64_t nlevels = static_cast<int64_t>(vtk_tree->GetNumberOfLevels());
    const int64_t nverts = static_cast<int64_t>(vtk_tree->GetNumberOfVertices());
    const int64_t nnodes = static_cast<int64_t>(vtk_tree->GetNumberOfNodes());
    const int64_t nleaves = static_cast<int64_t>(vtk_tree->GetNumberOfLeaves());

    ORACLE_EQ(nlevels, ftn_htg_tree_get_number_of_levels(ftn_tree), tlabel << " NumberOfLevels");
    ORACLE_EQ(
      nverts, ftn_htg_tree_get_number_of_vertices(ftn_tree), tlabel << " NumberOfVertices");
    ORACLE_EQ(nnodes, ftn_htg_tree_get_number_of_nodes(ftn_tree), tlabel << " NumberOfNodes");
    ORACLE_EQ(nleaves, ftn_htg_tree_get_number_of_leaves(ftn_tree), tlabel << " NumberOfLeaves");

    // Per-depth vertex count comparison
    std::vector<int64_t> expected_vpd(static_cast<size_t>(nlevels), 0);
    compute_vtk_vpd_recursive(vtk_tree, 0, 0, expected_vpd);
    const int64_t* ftn_vpd_ptr = nullptr;
    int64_t ftn_vpd_len = 0;
    ftn_htg_tree_get_vertices_per_depth_view(ftn_tree, &ftn_vpd_ptr, &ftn_vpd_len);
    ORACLE_EQ(nlevels, ftn_vpd_len, tlabel << " vertices_per_depth_view len");
    for (int64_t d = 0; d < nlevels; ++d)
    {
      int32_t st = 0;
      const int64_t gd =
        ftn_htg_tree_get_number_of_vertices_at_depth(ftn_tree, static_cast<int32_t>(d), &st);
      ORACLE_EQ(FTN_HTG_OK, st, tlabel << " get_number_of_vertices_at_depth status");
      ORACLE_EQ(expected_vpd[static_cast<size_t>(d)], gd, tlabel << " VPD(" << d << ")");
      ORACLE_EQ(expected_vpd[static_cast<size_t>(d)], ftn_vpd_ptr[d],
        tlabel << " VPD view[" << d << "]");
    }

    // ParentToElderChild & every vertex topology slot
    size_t vtk_ptec_len = 0;
    const unsigned int* vtk_ptec = vtk_tree->GetElderChildIndexArray(vtk_ptec_len);
    const int64_t* ftn_ptec = nullptr;
    int64_t ftn_ptec_len = 0;
    ftn_htg_tree_get_parent_to_elder_child_view(ftn_tree, &ftn_ptec, &ftn_ptec_len);
    ORACLE_EQ(static_cast<int64_t>(vtk_ptec_len), ftn_ptec_len,
      tlabel << " ParentToElderChild active size");

    for (int64_t v = 0; v < nverts; ++v)
    {
      int32_t st = 0;
      const bool v_leaf = vtk_tree->IsLeaf(static_cast<vtkIdType>(v));
      const bool f_leaf = ftn_htg_tree_is_leaf(ftn_tree, v, &st);
      ORACLE_EQ(FTN_HTG_OK, st, tlabel << " is_leaf status v=" << v);
      ORACLE_EQ(v_leaf, f_leaf, tlabel << " IsLeaf(" << v << ")");

      const bool v_term = !v_leaf && vtk_tree->IsTerminalNode(static_cast<vtkIdType>(v));
      const bool f_term = ftn_htg_tree_is_terminal_node(ftn_tree, v, &st);
      ORACLE_EQ(FTN_HTG_OK, st, tlabel << " is_terminal_node status v=" << v);
      ORACLE_EQ(v_term, f_term, tlabel << " IsTerminalNode(" << v << ")");

      if (!v_leaf)
      {
        const int64_t v_elder =
          static_cast<int64_t>(vtk_tree->GetElderChildIndex(static_cast<unsigned int>(v)));
        const int64_t f_elder = ftn_htg_tree_get_elder_child_index(ftn_tree, v, &st);
        ORACLE_EQ(FTN_HTG_OK, st, tlabel << " get_elder_child_index status v=" << v);
        ORACLE_EQ(v_elder, f_elder, tlabel << " GetElderChildIndex(" << v << ")");
        ORACLE_ASSERT(v < ftn_ptec_len, tlabel << " coarse v=" << v << " >= ftn_ptec_len");
        ORACLE_EQ(v_elder, ftn_ptec[v], tlabel << " ParentToElderChild[" << v << "]");
      }
      else if (v < ftn_ptec_len)
      {
        // Leaf slot inside ParentToElderChild must be leaf sentinel (UINT32_MAX or -1)
        const unsigned int vp = vtk_ptec[v];
        ORACLE_EQ(std::numeric_limits<unsigned int>::max(), vp, tlabel << " VTK leaf sentinel");
        ORACLE_ASSERT(
          ftn_ptec[v] == static_cast<int64_t>(std::numeric_limits<uint32_t>::max()) ||
            ftn_ptec[v] == -1,
          tlabel << " Fortran leaf sentinel at v=" << v << " got " << ftn_ptec[v]);
      }

      if (opts.check_global_indexing)
      {
        const int64_t v_gid =
          static_cast<int64_t>(vtk_tree->GetGlobalIndexFromLocal(static_cast<vtkIdType>(v)));
        const int64_t f_gid = ftn_htg_tree_get_global_index_from_local(ftn_tree, v, &st);
        ORACLE_EQ(FTN_HTG_OK, st, tlabel << " get_global_index_from_local status v=" << v);
        ORACLE_EQ(v_gid, f_gid, tlabel << " GetGlobalIndexFromLocal(" << v << ")");
      }
    }

    // Simultaneous NonOrientedCursor walk over active tree
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vtk_cur;
    vtk_grid->InitializeNonOrientedCursor(vtk_cur, vtk_tree_id, false);
    int32_t cst = 0;
    void* ftn_cur = ftn_htg_cursor_create(ftn_grid, ftn_tree_id, false, &cst);
    ORACLE_EQ(FTN_HTG_OK, cst, tlabel << " cursor_create");
    compare_cursors_recursive(
      vtk_grid, vtk_cur, ftn_cur, opts.check_global_indexing, opts.check_masks, tlabel);
    ftn_htg_cursor_destroy(ftn_cur);

    ++pos;
  }
  ORACLE_EQ(pos, ftn_htg_grid_get_number_of_non_empty_trees(ftn_grid),
    label << " iterated tree count");

  if (opts.check_global_indexing && pos > 0)
  {
    ORACLE_EQ(static_cast<int64_t>(vtk_grid->GetGlobalNodeIndexMax()),
      ftn_htg_grid_get_global_node_index_max(ftn_grid), label << " GetGlobalNodeIndexMax");
  }

  // Raw MSB-first Mask byte buffer comparison
  if (opts.check_masks)
  {
    const bool v_has = vtk_grid->HasMask();
    const bool f_has = ftn_htg_grid_has_mask(ftn_grid);
    ORACLE_EQ(v_has, f_has, label << " HasMask");
    if (v_has)
    {
      vtkBitArray* vmask = vtk_grid->GetMask();
      const int64_t v_nbits = static_cast<int64_t>(vmask->GetNumberOfValues());
      const int64_t v_nbytes = (v_nbits + 7) / 8;
      const uint8_t* v_bytes = reinterpret_cast<const uint8_t*>(vmask->GetPointer(0));

      const uint8_t* f_bytes = nullptr;
      int64_t f_nbits = 0, f_nbytes = 0;
      ftn_htg_grid_get_mask_raw_bytes(ftn_grid, &f_bytes, &f_nbits, &f_nbytes);
      ORACLE_EQ(v_nbits, f_nbits, label << " Mask num_bits");
      ORACLE_EQ(v_nbytes, f_nbytes, label << " Mask active_bytes");
      ORACLE_ASSERT(std::memcmp(v_bytes, f_bytes, static_cast<size_t>(v_nbytes)) == 0,
        label << " Mask raw MSB-first bytes mismatch!\n  VTK: [" << bytes_to_hex(v_bytes, v_nbytes)
              << "]\n  FTN: [" << bytes_to_hex(f_bytes, f_nbytes) << "]");
    }
  }

  // PureMask bit-for-bit & MSB-first byte-for-byte comparison
  if (opts.check_pure_mask && pos > 0)
  {
    vtkBitArray* vpure = vtk_grid->GetPureMask();
    ORACLE_ASSERT(vpure != nullptr, label << " VTK GetPureMask returned null");
    const int64_t v_nbits = static_cast<int64_t>(vpure->GetNumberOfValues());
    const int64_t v_nbytes = (v_nbits + 7) / 8;
    const uint8_t* v_bytes = reinterpret_cast<const uint8_t*>(vpure->GetPointer(0));

    const uint8_t* f_bytes = nullptr;
    int64_t f_nbits = 0, f_nbytes = 0;
    ORACLE_EQ(FTN_HTG_OK,
      ftn_htg_grid_get_pure_mask_raw_bytes(ftn_grid, &f_bytes, &f_nbits, &f_nbytes),
      label << " ftn_htg_grid_get_pure_mask_raw_bytes");
    ORACLE_EQ(v_nbits, f_nbits, label << " PureMask num_bits");
    ORACLE_EQ(v_nbytes, f_nbytes, label << " PureMask active_bytes");

    // Check every active tree vertex's PureMask bit first for clear per-vertex diagnostics
    vtkHyperTreeGrid::vtkHyperTreeGridIterator pit;
    vtk_grid->InitializeTreeIterator(pit);
    vtkIdType ptid = -1;
    while (vtkHyperTree* ptree = pit.GetNextTree(ptid))
    {
      const vtkIdType nv = ptree->GetNumberOfVertices();
      for (vtkIdType lv = 0; lv < nv; ++lv)
      {
        const vtkIdType gid = ptree->GetGlobalIndexFromLocal(lv);
        const bool vb = (vpure->GetValue(gid) != 0);
        int32_t pst = 0;
        const bool fb = ftn_htg_grid_get_pure_mask_bit(ftn_grid, static_cast<int64_t>(gid), &pst);
        ORACLE_EQ(FTN_HTG_OK, pst, label << " get_pure_mask_bit status gid=" << gid);
        ORACLE_EQ(vb, fb,
          label << " PureMask bit mismatch at tree=" << ptid << " local=" << lv
                << " global=" << gid);
      }
    }

    // Also compare raw MSB-first packed bytes directly
    if (v_nbytes > 0)
    {
      ORACLE_ASSERT(std::memcmp(v_bytes, f_bytes, static_cast<size_t>(v_nbytes)) == 0,
        label << " PureMask raw MSB-first bytes mismatch!\n  VTK: ["
              << bytes_to_hex(v_bytes, v_nbytes) << "]\n  FTN: ["
              << bytes_to_hex(f_bytes, f_nbytes) << "]");
    }
  }

  // Whole-grid Fortran structural invariant validation
  if (opts.check_validate)
  {
    char err_buf[512] = { 0 };
    const int32_t vst = ftn_htg_grid_validate(ftn_grid, err_buf, sizeof(err_buf));
    ORACLE_EQ(FTN_HTG_OK, vst, label << " htg%validate failed: " << err_buf);
  }
}

} // namespace oracle

#endif /* FVTKHTG_ORACLE_COMPARE_HPP */
