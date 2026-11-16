#include "oracle_compare.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

std::vector<double> make_stretched_coords(int npts, double x0, double scale)
{
  std::vector<double> c(static_cast<size_t>(npts));
  if (npts == 1)
  {
    c[0] = x0;
    return c;
  }
  for (int i = 0; i < npts; ++i)
  {
    const double t = static_cast<double>(i) / static_cast<double>(npts - 1);
    c[static_cast<size_t>(i)] = x0 + scale * (0.65 * t + 0.35 * t * t);
  }
  return c;
}

struct GridCase
{
  const char* name;
  int dims[3];
};

void run_case(const GridCase& gc, int branch_factor, bool transposed)
{
  const std::string tag = std::string(gc.name) + "_f" + std::to_string(branch_factor) + "_" +
    (transposed ? "IJK" : "KJI");

  vtkNew<vtkHyperTreeGrid> vgrid;
  void* fgrid = ftn_htg_grid_create();

  vgrid->SetBranchFactor(static_cast<unsigned int>(branch_factor));
  vgrid->SetTransposedRootIndexing(transposed);
  vgrid->SetDimensions(gc.dims[0], gc.dims[1], gc.dims[2]);

  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_branch_factor(fgrid, branch_factor), tag << " set_bf");
  ftn_htg_grid_set_transposed_root_indexing(fgrid, transposed);
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_set_dimensions(fgrid, gc.dims[0], gc.dims[1], gc.dims[2]), tag << " set_dims");

  const auto xc = make_stretched_coords(gc.dims[0], -2.5, 5.0);
  const auto yc = make_stretched_coords(gc.dims[1], 1.0, 3.75);
  const auto zc = make_stretched_coords(gc.dims[2], -0.5, 8.25);
  oracle::set_rectilinear_coords_both(vgrid, fgrid, xc, yc, zc);

  // Verify initial empty forest metadata
  oracle::CompareOptions opts;
  opts.check_coordinates = true;
  opts.check_global_indexing = false;
  opts.check_validate = true;
  oracle::assert_grids_equivalent(vgrid, fgrid, tag + "/empty", opts);

  // Exhaustively verify GetIndexFromLevelZeroCoordinates & GetLevelZeroCoordinatesFromIndex
  int cdims[3] = { 0, 0, 0 };
  vgrid->GetCellDims(cdims);
  for (int k = 0; k < cdims[2]; ++k)
  {
    for (int j = 0; j < cdims[1]; ++j)
    {
      for (int i = 0; i < cdims[0]; ++i)
      {
        vtkIdType vid = -1;
        vgrid->GetIndexFromLevelZeroCoordinates(
          vid, static_cast<unsigned int>(i), static_cast<unsigned int>(j),
          static_cast<unsigned int>(k));

        int64_t fid = -1;
        ORACLE_EQ(FTN_HTG_OK,
          ftn_htg_grid_get_index_from_level_zero_coordinates(fgrid, i, j, k, &fid),
          tag << " get_index_from_level_zero_coordinates(" << i << "," << j << "," << k << ")");
        ORACLE_EQ(static_cast<int64_t>(vid), fid,
          tag << " root index (" << i << "," << j << "," << k << ")");

        unsigned int vi = 999, vj = 999, vk = 999;
        vgrid->GetLevelZeroCoordinatesFromIndex(vid, vi, vj, vk);
        int32_t fi = 999, fj = 999, fk = 999;
        ORACLE_EQ(FTN_HTG_OK,
          ftn_htg_grid_get_level_zero_coordinates_from_index(fgrid, fid, &fi, &fj, &fk),
          tag << " get_level_zero_coordinates_from_index(" << fid << ")");
        ORACLE_EQ(static_cast<int32_t>(vi), fi, tag << " inverse i for tid=" << fid);
        ORACLE_EQ(static_cast<int32_t>(vj), fj, tag << " inverse j for tid=" << fid);
        ORACLE_EQ(static_cast<int32_t>(vk), fk, tag << " inverse k for tid=" << fid);
        ORACLE_EQ(i, fi, tag << " roundtrip i");
        ORACLE_EQ(j, fj, tag << " roundtrip j");
        ORACLE_EQ(k, fk, tag << " roundtrip k");
      }
    }
  }

  // Instantiate sparse subset of trees & refine with NonOrientedCursor
  const int64_t max_trees = static_cast<int64_t>(vgrid->GetMaxNumberOfTrees());
  const std::vector<int64_t> active_roots = { 0, max_trees / 2, max_trees - 1 };
  for (size_t idx = 0; idx < active_roots.size(); ++idx)
  {
    const int64_t tid = active_roots[idx];
    if (vgrid->GetTree(static_cast<vtkIdType>(tid), false) != nullptr)
    {
      continue;
    }
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vcur;
    vgrid->InitializeNonOrientedCursor(vcur, static_cast<vtkIdType>(tid), true);
    int32_t cst = 0;
    void* fcur = ftn_htg_cursor_create(fgrid, tid, true, &cst);
    ORACLE_EQ(FTN_HTG_OK, cst, tag << " cursor_create tid=" << tid);

    if (idx % 2 == 0)
    {
      vcur->SubdivideLeaf();
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), tag << " subdivide root");
      const int32_t nch = static_cast<int32_t>(vcur->GetNumberOfChildren());
      const int32_t pick = (nch - 1) / 2;
      vcur->ToChild(static_cast<unsigned char>(pick));
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fcur, pick), tag << " to_child");
      vcur->SubdivideLeaf();
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), tag << " subdivide lvl1");
      vcur->ToRoot();
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_root(fcur), tag << " to_root");
    }
    ftn_htg_cursor_destroy(fcur);
  }

  vgrid->InitializeLocalIndexNode();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_initialize_local_index_node(fgrid),
    tag << " initialize_local_index_node");

  opts.check_global_indexing = true;
  opts.check_validate = true;
  oracle::assert_grids_equivalent(vgrid, fgrid, tag + "/refined", opts);

  ftn_htg_grid_destroy(fgrid);
}

} // namespace

int main()
{
  oracle::verify_vtk_oracle_version();

  const GridCase cases[] = {
    { "1D_X", { 6, 1, 1 } },
    { "1D_Y", { 1, 7, 1 } },
    { "1D_Z", { 1, 1, 5 } },
    { "2D_XY", { 5, 4, 1 } },
    { "2D_YZ", { 1, 5, 6 } },
    { "2D_XZ", { 6, 1, 4 } },
    { "3D_XYZ", { 4, 3, 5 } },
    { "3D_RECT", { 7, 2, 4 } },
  };

  int count = 0;
  for (const auto& gc : cases)
  {
    for (int bf : { 2, 3 })
    {
      for (bool transposed : { false, true })
      {
        run_case(gc, bf, transposed);
        ++count;
      }
    }
  }

  std::cout << "[PASS] test_oracle_dimensions_coords verified all " << count
            << " (dimension/orientation x branch_factor x root_indexing) cases!" << std::endl;
  return 0;
}
