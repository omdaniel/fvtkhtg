#include "oracle_compare.hpp"

#include <iostream>
#include <vector>

namespace {

std::vector<double> uniform_coords(int n, double start, double step)
{
  std::vector<double> c(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
  {
    c[static_cast<size_t>(i)] = start + step * static_cast<double>(i);
  }
  return c;
}

void test_deep_asymmetric_amr_step_by_step(int branch_factor, int dim_mode)
{
  const std::string tag =
    "AMR_f" + std::to_string(branch_factor) + "_d" + std::to_string(dim_mode);

  vtkNew<vtkHyperTreeGrid> vgrid;
  void* fgrid = ftn_htg_grid_create();

  int dims[3] = { 12, 10, 8 };
  if (dim_mode == 2)
  {
    dims[2] = 1;
  }
  else if (dim_mode == 1)
  {
    dims[1] = 1;
    dims[2] = 1;
    dims[0] = 25;
  }

  vgrid->SetBranchFactor(static_cast<unsigned int>(branch_factor));
  vgrid->SetDimensions(dims[0], dims[1], dims[2]);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_branch_factor(fgrid, branch_factor), tag);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_dimensions(fgrid, dims[0], dims[1], dims[2]), tag);

  oracle::set_rectilinear_coords_both(vgrid, fgrid, uniform_coords(dims[0], 0.0, 1.0),
    uniform_coords(dims[1], -1.0, 0.5), uniform_coords(dims[2], 2.0, 0.25));

  const int64_t max_trees = static_cast<int64_t>(vgrid->GetMaxNumberOfTrees());
  // Insert trees in non-monotonic order to exercise sparse map sorting + cursor pointer stability
  std::vector<int64_t> insert_order = {
    max_trees - 1,
    max_trees / 3,
    0,
    (2 * max_trees) / 3,
    max_trees / 2,
  };

  oracle::CompareOptions step_opts;
  step_opts.check_coordinates = true;
  step_opts.check_global_indexing = false;
  step_opts.check_validate = false;

  // Keep a live cursor on tree `max_trees - 1` while inserting other trees to verify stability!
  vtkNew<vtkHyperTreeGridNonOrientedCursor> pinned_vcur;
  vgrid->InitializeNonOrientedCursor(
    pinned_vcur, static_cast<vtkIdType>(insert_order[0]), true);
  int32_t pst = 0;
  void* pinned_fcur = ftn_htg_cursor_create(fgrid, insert_order[0], true, &pst);
  ORACLE_EQ(FTN_HTG_OK, pst, tag << " pinned cursor create");

  // Refine tree `insert_order[0]` deeply to depth 6+ (7 levels total!)
  // For f=3 d=3 (27 children), depth 6 is 1 + 6*27 = 163 vertices on a single spire — super fast!
  const int target_depth = 6;
  const int32_t nch = static_cast<int32_t>(pinned_vcur->GetNumberOfChildren());
  for (int lvl = 0; lvl < target_depth; ++lvl)
  {
    pinned_vcur->SubdivideLeaf();
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(pinned_fcur),
      tag << " spire subdivide lvl=" << lvl);
    oracle::assert_grids_equivalent(
      vgrid, fgrid, tag + "/spire_step_" + std::to_string(lvl), step_opts);

    // Also split a sibling at levels 1 and 3 for rich asymmetric branching
    if (lvl == 1 || lvl == 3)
    {
      const int32_t sib = (lvl * 5 + 1) % nch;
      pinned_vcur->ToChild(static_cast<unsigned char>(sib));
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(pinned_fcur, sib), tag << " to_sib");
      pinned_vcur->SubdivideLeaf();
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(pinned_fcur), tag << " split_sib");
      oracle::assert_grids_equivalent(
        vgrid, fgrid, tag + "/sib_step_" + std::to_string(lvl), step_opts);
      pinned_vcur->ToParent();
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_parent(pinned_fcur), tag << " back_from_sib");
    }

    const int32_t next_ch = (lvl * 3 + (nch - 1)) % nch;
    pinned_vcur->ToChild(static_cast<unsigned char>(next_ch));
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(pinned_fcur, next_ch), tag << " descend");
  }

  // Instantiate remaining sparse trees while `pinned_fcur` stays alive at depth 6!
  for (size_t i = 1; i < insert_order.size(); ++i)
  {
    const int64_t tid = insert_order[i];
    if (vgrid->GetTree(static_cast<vtkIdType>(tid), false) != nullptr)
    {
      continue;
    }
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vcur;
    vgrid->InitializeNonOrientedCursor(vcur, static_cast<vtkIdType>(tid), true);
    int32_t cst = 0;
    void* fcur = ftn_htg_cursor_create(fgrid, tid, true, &cst);
    ORACLE_EQ(FTN_HTG_OK, cst, tag << " cursor create i=" << i);

    if (i == 1)
    {
      // Single-leaf root (leave unrefined)
    }
    else if (i == 2)
    {
      // Balanced 2-level refinement (all root children refined once)
      vcur->SubdivideLeaf();
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), tag);
      for (int32_t c = 0; c < nch; ++c)
      {
        vcur->ToChild(static_cast<unsigned char>(c));
        ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fcur, c), tag);
        vcur->SubdivideLeaf();
        ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), tag);
        oracle::assert_grids_equivalent(
          vgrid, fgrid, tag + "/tree2_child_" + std::to_string(c), step_opts);
        vcur->ToParent();
        ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_parent(fcur), tag);
      }
    }
    else
    {
      // 1-split terminal root
      vcur->SubdivideLeaf();
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), tag);
    }
    ftn_htg_cursor_destroy(fcur);
  }

  // Verify `pinned_fcur` is STILL valid after sparse forest growth and subdivide at depth 6 -> 7!
  ORACLE_EQ(static_cast<int64_t>(pinned_vcur->GetVertexId()),
    ftn_htg_cursor_get_vertex_id(pinned_fcur), tag << " pinned cursor vertex after insertions");
  pinned_vcur->SubdivideLeaf();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(pinned_fcur),
    tag << " pinned cursor subdivide depth 6->7");
  ftn_htg_cursor_destroy(pinned_fcur);

  // Remove one intermediate tree (`insert_order[3]`) from both grids and verify parity
  const int64_t rem_tid = insert_order[3];
  vgrid->RemoveTree(static_cast<vtkIdType>(rem_tid));
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_remove_tree(fgrid, rem_tid), tag << " remove_tree");
  oracle::assert_grids_equivalent(vgrid, fgrid, tag + "/after_remove", step_opts);

  // Assign Implicit Global Indexing via InitializeLocalIndexNode & validate full forest
  vgrid->InitializeLocalIndexNode();
  ORACLE_EQ(
    FTN_HTG_OK, ftn_htg_grid_initialize_local_index_node(fgrid), tag << " init_local_index");

  step_opts.check_global_indexing = true;
  step_opts.check_validate = true;
  oracle::assert_grids_equivalent(vgrid, fgrid, tag + "/implicit_indexed", step_opts);

  ftn_htg_grid_destroy(fgrid);
}

void test_explicit_and_custom_implicit_global_indexing()
{
  // Part A: Custom non-contiguous `SetGlobalIndexStart` across multi-tree forest
  {
    const std::string tag = "CustomImplicitOffsets";
    vtkNew<vtkHyperTreeGrid> vgrid;
    void* fgrid = ftn_htg_grid_create();
    vgrid->SetBranchFactor(2);
    vgrid->SetDimensions(4, 4, 1);
    ftn_htg_grid_set_branch_factor(fgrid, 2);
    ftn_htg_grid_set_dimensions(fgrid, 4, 4, 1);
    oracle::set_rectilinear_coords_both(vgrid, fgrid, uniform_coords(4, 0.0, 1.0),
      uniform_coords(4, 0.0, 1.0), { 0.0 });

    const int64_t tids[3] = { 1, 5, 8 };
    const int64_t starts[3] = { 100, 250, 500 };
    for (int t = 0; t < 3; ++t)
    {
      vtkNew<vtkHyperTreeGridNonOrientedCursor> vcur;
      vgrid->InitializeNonOrientedCursor(vcur, static_cast<vtkIdType>(tids[t]), true);
      int32_t cst = 0;
      void* fcur = ftn_htg_cursor_create(fgrid, tids[t], true, &cst);
      ORACLE_EQ(FTN_HTG_OK, cst, tag);

      vcur->SubdivideLeaf();
      ftn_htg_cursor_subdivide_leaf(fcur);
      if (t == 1)
      {
        vcur->ToChild(2);
        ftn_htg_cursor_to_child(fcur, 2);
        vcur->SubdivideLeaf();
        ftn_htg_cursor_subdivide_leaf(fcur);
        vcur->ToRoot();
        ftn_htg_cursor_to_root(fcur);
      }
      vcur->SetGlobalIndexStart(static_cast<vtkIdType>(starts[t]));
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_set_global_index_start(fcur, starts[t]), tag);
      ftn_htg_cursor_destroy(fcur);
    }

    oracle::CompareOptions opts;
    opts.check_coordinates = true;
    opts.check_global_indexing = true;
    opts.check_validate = true;
    oracle::assert_grids_equivalent(vgrid, fgrid, tag, opts);
    ftn_htg_grid_destroy(fgrid);
  }

  // Part B: Explicit `SetGlobalIndexFromLocal` with custom permuted non-contiguous IDs
  {
    const std::string tag = "ExplicitPermutedGlobalIndexing";
    vtkNew<vtkHyperTreeGrid> vgrid;
    void* fgrid = ftn_htg_grid_create();
    vgrid->SetBranchFactor(2);
    vgrid->SetDimensions(3, 3, 3);
    ftn_htg_grid_set_branch_factor(fgrid, 2);
    ftn_htg_grid_set_dimensions(fgrid, 3, 3, 3);
    oracle::set_rectilinear_coords_both(vgrid, fgrid, uniform_coords(3, 0.0, 1.0),
      uniform_coords(3, 0.0, 1.0), uniform_coords(3, 0.0, 1.0));

    // Build two trees (tree 0 and tree 6) via cursor, then assign explicit non-contiguous global IDs
    int64_t next_gid = 17;
    auto assign_explicit_dfs = [&](auto& self, vtkHyperTreeGridNonOrientedCursor* vc, void* fc) -> void
    {
      const int64_t gid = next_gid;
      next_gid += 5; // non-contiguous stride-5 unique global IDs
      vc->SetGlobalIndexFromLocal(static_cast<vtkIdType>(gid));
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_set_global_index_from_local(fc, gid),
        tag << " cursor set_global_index_from_local");
      if (!vc->IsLeaf())
      {
        const int32_t nch = static_cast<int32_t>(vc->GetNumberOfChildren());
        for (int32_t c = 0; c < nch; ++c)
        {
          vc->ToChild(static_cast<unsigned char>(c));
          ftn_htg_cursor_to_child(fc, c);
          self(self, vc, fc);
          vc->ToParent();
          ftn_htg_cursor_to_parent(fc);
        }
      }
    };

    for (int64_t tid : { int64_t(0), int64_t(6) })
    {
      vtkNew<vtkHyperTreeGridNonOrientedCursor> vc;
      vgrid->InitializeNonOrientedCursor(vc, static_cast<vtkIdType>(tid), true);
      int32_t cst = 0;
      void* fc = ftn_htg_cursor_create(fgrid, tid, true, &cst);
      ORACLE_EQ(FTN_HTG_OK, cst, tag);

      vc->SubdivideLeaf();
      ftn_htg_cursor_subdivide_leaf(fc);
      vc->ToChild(static_cast<unsigned char>(tid == 0 ? 3 : 7));
      ftn_htg_cursor_to_child(fc, tid == 0 ? 3 : 7);
      vc->SubdivideLeaf();
      ftn_htg_cursor_subdivide_leaf(fc);
      vc->ToRoot();
      ftn_htg_cursor_to_root(fc);

      assign_explicit_dfs(assign_explicit_dfs, vc, fc);

      // Verify zero-copy global_index_table_view on Fortran tree
      void* ftree = ftn_htg_cursor_get_tree(fc);
      const int64_t* gtab = nullptr;
      int64_t gtab_len = 0;
      ftn_htg_tree_get_global_index_table_view(ftree, &gtab, &gtab_len);
      ORACLE_EQ(ftn_htg_tree_get_number_of_vertices(ftree), gtab_len, tag << " gtab_len");
      ORACLE_EQ(FTN_HTG_INDEX_EXPLICIT, ftn_htg_tree_get_indexing_mode(ftree), tag << " mode");

      ftn_htg_cursor_destroy(fc);
    }

    oracle::CompareOptions opts;
    opts.check_coordinates = true;
    opts.check_global_indexing = true;
    opts.check_validate = true;
    oracle::assert_grids_equivalent(vgrid, fgrid, tag, opts);
    ftn_htg_grid_destroy(fgrid);
  }
}

} // namespace

int main()
{
  oracle::verify_vtk_oracle_version();

  for (int bf : { 2, 3 })
  {
    for (int d : { 1, 2, 3 })
    {
      test_deep_asymmetric_amr_step_by_step(bf, d);
    }
  }
  test_explicit_and_custom_implicit_global_indexing();

  std::cout << "[PASS] test_oracle_amr_indexing verified deep AMR & implicit/explicit indexing!"
            << std::endl;
  return 0;
}
