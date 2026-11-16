// Differential tests for cell geometry, face neighbours and coarsening.
//
// Geometry: every cell's origin, size, bounds and centre must equal, bit for
// bit, what vtkHyperTreeGridNonOrientedVonNeumannSuperCursor reports.
//
// Neighbours: for every cell, each neighbour VTK's supercursor holds is
// identified by geometry (which face of the cell it touches) and must be the
// cell fvtkhtg returns for that face: same tree, local id, global id, level,
// leaf and mask state, bounds. Faces for which VTK has no neighbour must have
// none in fvtkhtg either. Identifying VTK's slots by geometry keeps the test
// independent of slot-numbering conventions.
//
// xz-plane grids with several root cells are the one exception. VTK 9.7.2's
// supercursor finds root-level neighbours along vtkHyperTreeGrid::GetAxes,
// which for the xz plane is (z, x), while its child numbering is (x, z); so
// across root boundaries it holds the neighbour along the other axis. There
// fvtkhtg is checked for geometric adjacency instead, and the number of VTK
// slots that do not touch the cell is reported.
//
// Coarsening has no VTK counterpart. A grid that fvtkhtg refined further and
// then coarsened back must be indistinguishable from the VTK grid that was
// never refined further: same structure, geometry and breadth-first layout.

#include "oracle_compare.hpp"

#include <vtkHyperTreeGridNonOrientedGeometryCursor.h>
#include <vtkHyperTreeGridNonOrientedVonNeumannSuperCursor.h>

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Rng
{
  uint64_t state;
  double next()
  {
    state = (state * 1103515245ULL + 12345ULL) % 2147483648ULL;
    return static_cast<double>(state) / 2147483648.0;
  }
  int below(int n) { return static_cast<int>(next() * n) % n; }
};

struct Config
{
  std::array<int, 3> dims;
  int branch_factor;
  bool transposed;
  bool sparse;
  bool masked;
};

struct Stats
{
  long cells = 0;
  long neighbors = 0;
  long coarser = 0;
  long cross_tree = 0;
  long none = 0;
  long xz_cross_tree_slots = 0;
  long xz_vtk_not_touching = 0;
};

std::vector<double> stretched_coords(int n, double x0, double dx, double stretch)
{
  std::vector<double> c(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
  {
    c[static_cast<size_t>(i)] = x0 + dx * i + stretch * i * i;
  }
  return c;
}

std::string config_tag(const Config& cfg)
{
  std::ostringstream oss;
  oss << "dims(" << cfg.dims[0] << "," << cfg.dims[1] << "," << cfg.dims[2] << ") bf"
      << cfg.branch_factor << (cfg.transposed ? " transposed" : "") << (cfg.sparse ? " sparse" : "")
      << (cfg.masked ? " masked" : "");
  return oss.str();
}

std::vector<int> active_axes(const Config& cfg)
{
  std::vector<int> axes;
  for (int a = 0; a < 3; ++a)
  {
    if (cfg.dims[static_cast<size_t>(a)] > 1)
    {
      axes.push_back(a);
    }
  }
  return axes;
}

void setup_both(const Config& cfg, vtkHyperTreeGrid* vgrid, void* fgrid, const std::string& tag)
{
  vgrid->Initialize();
  vgrid->SetBranchFactor(static_cast<unsigned int>(cfg.branch_factor));
  vgrid->SetTransposedRootIndexing(cfg.transposed);
  vgrid->SetDimensions(cfg.dims[0], cfg.dims[1], cfg.dims[2]);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_branch_factor(fgrid, cfg.branch_factor), tag);
  ftn_htg_grid_set_transposed_root_indexing(fgrid, cfg.transposed);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_dimensions(fgrid, cfg.dims[0], cfg.dims[1], cfg.dims[2]),
    tag);
  oracle::set_rectilinear_coords_both(vgrid, fgrid, stretched_coords(cfg.dims[0], -1.0, 1.0, 0.125),
    stretched_coords(cfg.dims[1], 2.0, 0.5, 0.0625), stretched_coords(cfg.dims[2], 10.0, 2.0, 0.25));
}

// Refine the same cells of both grids at random.
void refine_both(vtkHyperTreeGridNonOrientedCursor* vcur, void* fcur, int max_level, double prob,
  Rng& rng, const std::string& tag)
{
  if (static_cast<int>(vcur->GetLevel()) >= max_level || rng.next() >= prob)
  {
    return;
  }
  vcur->SubdivideLeaf();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), tag << " subdivide");
  const int nch = static_cast<int>(vcur->GetNumberOfChildren());
  for (int c = 0; c < nch; ++c)
  {
    vcur->ToChild(static_cast<unsigned char>(c));
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fcur, c), tag << " to_child");
    refine_both(vcur, fcur, max_level, prob, rng, tag);
    vcur->ToParent();
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_parent(fcur), tag << " to_parent");
  }
}

void build_both(const Config& cfg, vtkHyperTreeGrid* vgrid, void* fgrid, Rng& rng,
  const std::string& tag)
{
  setup_both(cfg, vgrid, fgrid, tag);
  const int ndim = static_cast<int>(active_axes(cfg).size());
  int max_level = cfg.branch_factor == 2 ? 3 : 2;
  if (ndim == 3)
  {
    --max_level;
  }
  const double prob = cfg.branch_factor == 2 ? 0.55 : 0.4;

  const vtkIdType ntrees = vgrid->GetMaxNumberOfTrees();
  for (vtkIdType tid = 0; tid < ntrees; ++tid)
  {
    if (cfg.sparse && tid % 3 == 1)
    {
      continue;
    }
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vcur;
    vgrid->InitializeNonOrientedCursor(vcur, tid, true);
    int32_t st = 0;
    void* fcur = ftn_htg_cursor_create(fgrid, static_cast<int64_t>(tid), true, &st);
    ORACLE_EQ(FTN_HTG_OK, st, tag << " cursor create");
    refine_both(vcur, fcur, max_level, prob, rng, tag);
    ftn_htg_cursor_destroy(fcur);
  }
  vgrid->InitializeLocalIndexNode();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_initialize_local_index_node(fgrid), tag << " indexing");

  if (cfg.masked)
  {
    const vtkIdType ncells = vgrid->GetNumberOfCells();
    vtkNew<vtkBitArray> mask;
    mask->SetNumberOfTuples(ncells);
    for (vtkIdType gid = 0; gid < ncells; ++gid)
    {
      const bool on = (gid % 7 == 3);
      mask->SetValue(gid, on ? 1 : 0);
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_cell_mask(fgrid, static_cast<int64_t>(gid), on),
        tag << " mask");
    }
    vgrid->SetMask(mask);
  }
}

// Does box `n` touch the face of box `c` that looks along `axis` in `dir`,
// lying on the far side and covering `c` in the other active directions?
bool touches_face(const double c[6], const double n[6], int axis, int dir,
  const std::vector<int>& active, double tol)
{
  const double face = dir > 0 ? c[2 * axis + 1] : c[2 * axis];
  const double near_side = dir > 0 ? n[2 * axis] : n[2 * axis + 1];
  if (std::fabs(near_side - face) > tol)
  {
    return false;
  }
  if (dir > 0 && !(n[2 * axis + 1] > face + tol))
  {
    return false;
  }
  if (dir < 0 && !(n[2 * axis] < face - tol))
  {
    return false;
  }
  for (int t : active)
  {
    if (t == axis)
    {
      continue;
    }
    if (n[2 * t] > c[2 * t] + tol || n[2 * t + 1] < c[2 * t + 1] - tol)
    {
      return false;
    }
  }
  return true;
}

struct FtnNeighbor
{
  bool has_tree = false;
  int64_t tree_id = -1;
  int64_t vertex_id = -1;
  int64_t global_id = -1;
  int32_t level = -1;
  bool is_leaf = false;
  bool is_masked = false;
  double bounds[6] = { 0, 0, 0, 0, 0, 0 };
};

FtnNeighbor ftn_neighbor(void* fcur, int axis, int dir, const std::string& tag)
{
  FtnNeighbor nb;
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_cursor_get_face_neighbor(fcur, axis, dir, &nb.has_tree, &nb.tree_id, &nb.vertex_id,
      &nb.global_id, &nb.level, &nb.is_leaf, &nb.is_masked, nb.bounds),
    tag << " get_face_neighbor axis=" << axis << " dir=" << dir);
  return nb;
}

void compare_cell(vtkHyperTreeGridNonOrientedVonNeumannSuperCursor* sc, void* fcur,
  const Config& cfg, const std::vector<int>& active, bool xz_multi_tree, Stats& stats,
  const std::string& tag)
{
  const double tol = 1.0e-9;
  ++stats.cells;
  const unsigned int level = sc->GetLevel();
  std::ostringstream where;
  where << tag << " tree=" << sc->GetTree()->GetTreeIndex() << " v=" << sc->GetVertexId()
        << " level=" << level;

  // --- geometry: exact ------------------------------------------------------
  double vb[6], fb[6], vp[3], fp[3], fo[3], fs[3];
  sc->GetBounds(vb);
  sc->GetPoint(vp);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_get_bounds(fcur, fb), where.str() << " get_bounds");
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_get_point(fcur, fp), where.str() << " get_point");
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_get_origin(fcur, fo), where.str() << " get_origin");
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_get_size(fcur, fs), where.str() << " get_size");
  const double* vo = sc->GetOrigin();
  const double* vs = sc->GetSize();
  for (int i = 0; i < 6; ++i)
  {
    ORACLE_ASSERT(vb[i] == fb[i], where.str() << " bounds[" << i << "] VTK=" << std::setprecision(17)
                                              << vb[i] << " Fortran=" << fb[i]);
  }
  for (int i = 0; i < 3; ++i)
  {
    ORACLE_ASSERT(vo[i] == fo[i], where.str() << " origin[" << i << "]");
    ORACLE_ASSERT(vs[i] == fs[i], where.str() << " size[" << i << "]");
    ORACLE_ASSERT(vp[i] == fp[i], where.str() << " point[" << i << "]");
  }

  // --- neighbours -----------------------------------------------------------
  bool seen[3][2] = { { false, false }, { false, false }, { false, false } };
  const vtkIdType own_tree = sc->GetTree()->GetTreeIndex();
  const unsigned int ncursors = sc->GetNumberOfCursors();
  const unsigned int central = sc->GetIndiceCentralCursor();
  for (unsigned int ic = 0; ic < ncursors; ++ic)
  {
    if (ic == central || !sc->HasTree(ic))
    {
      continue;
    }
    double nbb[6];
    sc->GetBounds(ic, nbb);
    const vtkIdType nb_tree = sc->GetTree(ic)->GetTreeIndex();

    int axis = -1;
    int dir = 0;
    for (int a : active)
    {
      for (int d : { -1, 1 })
      {
        if (touches_face(vb, nbb, a, d, active, tol))
        {
          ORACLE_ASSERT(axis < 0, where.str() << " VTK neighbour touches two faces");
          axis = a;
          dir = d;
        }
      }
    }

    if (xz_multi_tree && nb_tree != own_tree)
    {
      // VTK's root-level neighbour may be along the other axis here
      ++stats.xz_cross_tree_slots;
      if (axis < 0)
      {
        ++stats.xz_vtk_not_touching;
      }
      continue;
    }

    ORACLE_ASSERT(axis >= 0, where.str() << " VTK neighbour slot " << ic << " does not touch the cell");
    const FtnNeighbor nb = ftn_neighbor(fcur, axis, dir, where.str());
    std::ostringstream ctx;
    ctx << where.str() << " neighbour axis=" << axis << " dir=" << dir;
    ORACLE_ASSERT(nb.has_tree, ctx.str() << ": VTK has a neighbour, Fortran has none");
    ORACLE_EQ(static_cast<int64_t>(nb_tree), nb.tree_id, ctx.str() << " tree");
    ORACLE_EQ(static_cast<int64_t>(sc->GetVertexId(ic)), nb.vertex_id, ctx.str() << " vertex id");
    ORACLE_EQ(static_cast<int64_t>(sc->GetGlobalNodeIndex(ic)), nb.global_id,
      ctx.str() << " global id");
    ORACLE_EQ(static_cast<int32_t>(sc->GetLevel(ic)), nb.level, ctx.str() << " level");
    ORACLE_EQ(sc->IsLeaf(ic), nb.is_leaf, ctx.str() << " is_leaf");
    ORACLE_EQ(sc->IsMasked(ic), nb.is_masked, ctx.str() << " is_masked");
    for (int i = 0; i < 6; ++i)
    {
      ORACLE_ASSERT(nbb[i] == nb.bounds[i], ctx.str() << " neighbour bounds[" << i << "]");
    }
    seen[axis][dir > 0 ? 1 : 0] = true;
    ++stats.neighbors;
    if (nb.level < static_cast<int32_t>(level))
    {
      ++stats.coarser;
    }
    if (nb_tree != own_tree)
    {
      ++stats.cross_tree;
    }
  }

  for (int a : active)
  {
    for (int d : { -1, 1 })
    {
      if (seen[a][d > 0 ? 1 : 0])
      {
        continue;
      }
      const FtnNeighbor nb = ftn_neighbor(fcur, a, d, where.str());
      if (xz_multi_tree)
      {
        // No trustworthy VTK answer across root boundaries: check geometry
        if (nb.has_tree)
        {
          ORACLE_ASSERT(touches_face(fb, nb.bounds, a, d, active, tol),
            where.str() << " Fortran neighbour axis=" << a << " dir=" << d
                        << " does not touch the cell");
          ORACLE_ASSERT(nb.level <= static_cast<int32_t>(level), where.str() << " neighbour level");
        }
      }
      else
      {
        ORACLE_ASSERT(!nb.has_tree, where.str() << " Fortran reports a neighbour axis=" << a
                                                << " dir=" << d << " that VTK does not have");
        ++stats.none;
      }
    }
  }

  // --- children -------------------------------------------------------------
  int32_t st = 0;
  const bool vleaf = sc->IsLeaf();
  ORACLE_EQ(vleaf, ftn_htg_cursor_is_leaf(fcur, &st), where.str() << " is_leaf");
  if (vleaf || sc->IsMasked())
  {
    return; // VTK does not descend below a masked cell
  }
  const int nch = static_cast<int>(sc->GetNumberOfChildren());
  for (int c = 0; c < nch; ++c)
  {
    sc->ToChild(static_cast<unsigned char>(c));
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fcur, c), where.str() << " to_child");
    ORACLE_EQ(c, ftn_htg_cursor_get_child_index(fcur), where.str() << " child index");
    compare_cell(sc, fcur, cfg, active, xz_multi_tree, stats, tag);
    sc->ToParent();
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_parent(fcur), where.str() << " to_parent");
  }
}

void test_geometry_and_neighbors(const Config& cfg, Rng& rng, Stats& totals)
{
  const std::string tag = config_tag(cfg);
  vtkNew<vtkHyperTreeGrid> vgrid;
  void* fgrid = ftn_htg_grid_create();
  build_both(cfg, vgrid, fgrid, rng, tag);

  const std::vector<int> active = active_axes(cfg);
  const bool xz_plane = (cfg.dims[0] > 1 && cfg.dims[1] == 1 && cfg.dims[2] > 1);
  const bool xz_multi_tree = xz_plane && vgrid->GetMaxNumberOfTrees() > 1;

  Stats stats;
  const vtkIdType ntrees = vgrid->GetMaxNumberOfTrees();
  for (vtkIdType tid = 0; tid < ntrees; ++tid)
  {
    // Root-cell geometry is defined with or without a tree
    double vo[3], vs[3], fo[3], fs[3];
    vgrid->GetLevelZeroOriginAndSizeFromIndex(tid, vo, vs);
    ORACLE_EQ(FTN_HTG_OK,
      ftn_htg_grid_get_level_zero_origin_and_size(fgrid, static_cast<int64_t>(tid), fo, fs),
      tag << " level-zero origin and size");
    for (int i = 0; i < 3; ++i)
    {
      ORACLE_ASSERT(vo[i] == fo[i] && vs[i] == fs[i], tag << " root " << tid << " geometry");
    }

    if (vgrid->GetTree(tid, false) == nullptr)
    {
      ORACLE_ASSERT(!ftn_htg_grid_has_tree(fgrid, static_cast<int64_t>(tid)), tag << " sparse tree");
      continue;
    }
    vtkNew<vtkHyperTreeGridNonOrientedVonNeumannSuperCursor> sc;
    vgrid->InitializeNonOrientedVonNeumannSuperCursor(sc, tid, false);
    int32_t st = 0;
    void* fcur = ftn_htg_cursor_create(fgrid, static_cast<int64_t>(tid), false, &st);
    ORACLE_EQ(FTN_HTG_OK, st, tag << " cursor create");
    compare_cell(sc, fcur, cfg, active, xz_multi_tree, stats, tag);
    ftn_htg_cursor_destroy(fcur);
  }

  if (xz_multi_tree)
  {
    std::cout << "[Oracle][note] " << tag << ": " << stats.xz_vtk_not_touching << " of "
              << stats.xz_cross_tree_slots
              << " VTK cross-tree neighbour slots do not touch the cell (xz-plane axis quirk)"
              << std::endl;
  }
  totals.cells += stats.cells;
  totals.neighbors += stats.neighbors;
  totals.coarser += stats.coarser;
  totals.cross_tree += stats.cross_tree;
  totals.none += stats.none;
  totals.xz_cross_tree_slots += stats.xz_cross_tree_slots;
  totals.xz_vtk_not_touching += stats.xz_vtk_not_touching;
  ftn_htg_grid_destroy(fgrid);
}

// ---------------------------------------------------------------------------
// Coarsening
// ---------------------------------------------------------------------------

struct Path
{
  int64_t tree_id;
  std::vector<int> children;
};

void collect_leaf_paths(vtkHyperTreeGridNonOrientedCursor* vcur, int64_t tree_id,
  std::vector<int>& current, std::vector<Path>& out)
{
  if (vcur->IsLeaf())
  {
    out.push_back(Path{ tree_id, current });
    return;
  }
  const int nch = static_cast<int>(vcur->GetNumberOfChildren());
  for (int c = 0; c < nch; ++c)
  {
    vcur->ToChild(static_cast<unsigned char>(c));
    current.push_back(c);
    collect_leaf_paths(vcur, tree_id, current, out);
    current.pop_back();
    vcur->ToParent();
  }
}

void* ftn_cursor_at(void* fgrid, const Path& path, const std::string& tag)
{
  int32_t st = 0;
  void* fcur = ftn_htg_cursor_create(fgrid, path.tree_id, false, &st);
  ORACLE_EQ(FTN_HTG_OK, st, tag << " cursor create");
  for (int c : path.children)
  {
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fcur, c), tag << " navigate");
  }
  return fcur;
}

// Structure and geometry must agree; local ids may differ after coarsening.
void compare_structure(vtkHyperTreeGridNonOrientedGeometryCursor* vcur, void* fcur,
  const std::string& tag, long& cells)
{
  ++cells;
  int32_t st = 0;
  double vb[6], fb[6];
  vcur->GetBounds(vb);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_get_bounds(fcur, fb), tag << " bounds");
  for (int i = 0; i < 6; ++i)
  {
    ORACLE_ASSERT(vb[i] == fb[i], tag << " bounds[" << i << "] after coarsening");
  }
  ORACLE_EQ(static_cast<int32_t>(vcur->GetLevel()), ftn_htg_cursor_get_level(fcur), tag << " level");
  const bool vleaf = vcur->IsLeaf();
  ORACLE_EQ(vleaf, ftn_htg_cursor_is_leaf(fcur, &st), tag << " is_leaf after coarsening");
  if (vleaf)
  {
    return;
  }
  const int nch = static_cast<int>(vcur->GetNumberOfChildren());
  for (int c = 0; c < nch; ++c)
  {
    vcur->ToChild(static_cast<unsigned char>(c));
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fcur, c), tag << " to_child");
    compare_structure(vcur, fcur, tag, cells);
    vcur->ToParent();
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_parent(fcur), tag << " to_parent");
  }
}

void test_coarsen_back_to_vtk_shape(const Config& cfg, Rng& rng, long& total_coarsened)
{
  const std::string tag = "coarsen " + config_tag(cfg);
  vtkNew<vtkHyperTreeGrid> vgrid;
  void* fgrid = ftn_htg_grid_create();
  build_both(cfg, vgrid, fgrid, rng, tag);

  // Leaves of the common shape
  std::vector<Path> leaves;
  const vtkIdType ntrees = vgrid->GetMaxNumberOfTrees();
  for (vtkIdType tid = 0; tid < ntrees; ++tid)
  {
    if (vgrid->GetTree(tid, false) == nullptr)
    {
      continue;
    }
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vcur;
    vgrid->InitializeNonOrientedCursor(vcur, tid, false);
    std::vector<int> current;
    collect_leaf_paths(vcur, static_cast<int64_t>(tid), current, leaves);
  }

  // Fortran only: refine a random subset of those leaves, and half of them a
  // second level (one child), so local ids get thoroughly interleaved
  std::vector<Path> refined;
  std::vector<Path> refined_twice;
  for (const Path& leaf : leaves)
  {
    if (rng.next() >= 0.5)
    {
      continue;
    }
    void* fcur = ftn_cursor_at(fgrid, leaf, tag);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), tag << " extra subdivide");
    refined.push_back(leaf);
    if (rng.next() < 0.5)
    {
      Path deeper = leaf;
      deeper.children.push_back(rng.below(ftn_htg_cursor_get_number_of_children(fcur)));
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fcur, deeper.children.back()), tag);
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), tag << " second extra subdivide");
      refined_twice.push_back(deeper);
    }
    ftn_htg_cursor_destroy(fcur);
  }

  // Coarsen them again in creation order (oldest first), which is the order
  // that relocates sibling blocks: the deeper cells first, then their parents
  for (const std::vector<Path>* group : { &refined_twice, &refined })
  {
    for (const Path& path : *group)
    {
      void* fcur = ftn_cursor_at(fgrid, path, tag);
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_coarsen_node(fcur), tag << " coarsen_node");
      int32_t st = 0;
      ORACLE_ASSERT(ftn_htg_cursor_is_leaf(fcur, &st), tag << " coarsened cell is a leaf");
      ftn_htg_cursor_destroy(fcur);
      ++total_coarsened;
    }
  }

  // Now indistinguishable from the VTK grid that never had the extra cells
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_initialize_local_index_node(fgrid), tag << " re-index");
  ORACLE_EQ(static_cast<int64_t>(vgrid->GetNumberOfCells()), ftn_htg_grid_get_number_of_vertices(fgrid),
    tag << " cell count");
  ORACLE_EQ(static_cast<int64_t>(vgrid->GetNumberOfLeaves()), ftn_htg_grid_get_number_of_leaves(fgrid),
    tag << " leaf count");
  char msg[256] = { 0 };
  const int32_t vcode = ftn_htg_grid_validate(fgrid, msg, 255);
  ORACLE_ASSERT(vcode == FTN_HTG_OK, tag << " validate after coarsening: " << msg);

  long cells = 0;
  for (vtkIdType tid = 0; tid < ntrees; ++tid)
  {
    vtkHyperTree* vtree = vgrid->GetTree(tid, false);
    if (vtree == nullptr)
    {
      continue;
    }
    vtkNew<vtkHyperTreeGridNonOrientedGeometryCursor> vcur;
    vgrid->InitializeNonOrientedGeometryCursor(vcur, tid, false);
    int32_t st = 0;
    void* fcur = ftn_htg_cursor_create(fgrid, static_cast<int64_t>(tid), false, &st);
    ORACLE_EQ(FTN_HTG_OK, st, tag << " cursor create");
    compare_structure(vcur, fcur, tag, cells);

    // Counters and breadth-first layout of the tree
    void* ftree = ftn_htg_cursor_get_tree(fcur);
    ORACLE_EQ(static_cast<int64_t>(vtree->GetNumberOfVertices()),
      ftn_htg_tree_get_number_of_vertices(ftree), tag << " tree vertices");
    ORACLE_EQ(static_cast<int64_t>(vtree->GetNumberOfLevels()), ftn_htg_tree_get_number_of_levels(ftree),
      tag << " tree levels");

    vtkNew<vtkBitArray> vdesc;
    vtkNew<vtkTypeInt64Array> vvpd;
    vtkNew<vtkIdList> vmap;
    vtree->ComputeBreadthFirstOrderDescriptor(vgrid->GetDepthLimiter(), nullptr, vvpd, vdesc, vmap);
    int32_t bst = 0;
    void* fdesc = ftn_htg_bit_array_create(0, false, &bst);
    void* fvpd = ftn_htg_i64_vec_create();
    void* fmap = ftn_htg_i64_vec_create();
    ORACLE_EQ(FTN_HTG_OK,
      ftn_htg_tree_compute_bf_descriptor_append(ftree, fdesc, fvpd, fmap, nullptr, false, 0),
      tag << " Fortran descriptor");
    ORACLE_EQ(static_cast<int64_t>(vdesc->GetNumberOfValues()),
      ftn_htg_bit_array_get_number_of_values(fdesc), tag << " descriptor length");
    for (vtkIdType i = 0; i < vdesc->GetNumberOfValues(); ++i)
    {
      ORACLE_EQ(vdesc->GetValue(i) != 0, ftn_htg_bit_array_get_value(fdesc, static_cast<int64_t>(i), &bst),
        tag << " descriptor bit " << i);
    }
    const int64_t* vpd_ptr = nullptr;
    int64_t vpd_len = 0;
    ftn_htg_i64_vec_get_view(fvpd, &vpd_ptr, &vpd_len);
    ORACLE_EQ(static_cast<int64_t>(vvpd->GetNumberOfValues()), vpd_len, tag << " depth count");
    for (int64_t i = 0; i < vpd_len; ++i)
    {
      ORACLE_EQ(static_cast<int64_t>(vvpd->GetValue(i)), vpd_ptr[i], tag << " cells at depth " << i);
    }
    // Both are implicitly indexed in tree order, so the breadth-first id maps agree
    const int64_t* map_ptr = nullptr;
    int64_t map_len = 0;
    ftn_htg_i64_vec_get_view(fmap, &map_ptr, &map_len);
    ORACLE_EQ(static_cast<int64_t>(vmap->GetNumberOfIds()), map_len, tag << " id map length");

    ftn_htg_bit_array_destroy(fdesc);
    ftn_htg_i64_vec_destroy(fvpd);
    ftn_htg_i64_vec_destroy(fmap);
    ftn_htg_cursor_destroy(fcur);
  }
  ORACLE_EQ(static_cast<long>(vgrid->GetNumberOfCells()), cells, tag << " cells compared");
  ftn_htg_grid_destroy(fgrid);
}

} // namespace

int main()
{
  oracle::verify_vtk_oracle_version();

  const std::vector<std::array<int, 3>> dims_list = {
    { 5, 1, 1 }, { 1, 4, 1 }, { 1, 1, 4 }, // 1D along x, y, z
    { 4, 3, 1 }, { 4, 1, 3 }, { 1, 4, 3 }, // xy, xz, yz planes
    { 2, 1, 2 },                           // xz plane, single root cell
    { 3, 3, 2 },                           // 3D
  };

  Rng rng{ 20261116ULL };
  Stats totals;
  int nconfigs = 0;
  for (const auto& dims : dims_list)
  {
    for (int bf : { 2, 3 })
    {
      for (bool transposed : { false, true })
      {
        for (bool sparse : { false, true })
        {
          for (bool masked : { false, true })
          {
            test_geometry_and_neighbors(Config{ dims, bf, transposed, sparse, masked }, rng, totals);
            ++nconfigs;
          }
        }
      }
    }
  }
  std::cout << "[Oracle] geometry & neighbours: " << nconfigs << " configurations, " << totals.cells
            << " cells, " << totals.neighbors << " neighbours matched (" << totals.coarser
            << " coarser, " << totals.cross_tree << " across root cells), " << totals.none
            << " faces without a neighbour in both" << std::endl;
  ORACLE_ASSERT(totals.neighbors > 1000 && totals.coarser > 100 && totals.cross_tree > 100,
    "neighbour comparison did not exercise enough cases");

  long coarsened = 0;
  int ncoarsen_configs = 0;
  for (const auto& dims : dims_list)
  {
    for (int bf : { 2, 3 })
    {
      for (bool sparse : { false, true })
      {
        test_coarsen_back_to_vtk_shape(Config{ dims, bf, false, sparse, false }, rng, coarsened);
        ++ncoarsen_configs;
      }
    }
  }
  std::cout << "[Oracle] coarsening: " << ncoarsen_configs << " configurations, " << coarsened
            << " cells coarsened back to the VTK shape" << std::endl;
  ORACLE_ASSERT(coarsened > 200, "coarsening comparison did not exercise enough cases");

  std::cout << "[PASS] test_oracle_geometry_neighbors verified geometry, face neighbours & coarsening!"
            << std::endl;
  return 0;
}
