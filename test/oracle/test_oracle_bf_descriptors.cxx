#include "oracle_compare.hpp"

#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

namespace {

std::vector<double> lin_coords(int n)
{
  std::vector<double> c(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
  {
    c[static_cast<size_t>(i)] = static_cast<double>(i);
  }
  return c;
}

// Helper to append bits to both a `vtkBitArray` and a Fortran `htg_bit_array_t`
void append_bits_both(vtkBitArray* vb, void* fb, const std::vector<bool>& bits)
{
  for (bool b : bits)
  {
    vb->InsertNextValue(b ? 1 : 0);
    int64_t fid = -1;
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_bit_array_insert_next_value(fb, b, &fid), "append_bits_both");
  }
}

void test_build_from_bf_descriptor_all_combinations(int branch_factor, int dim_mode)
{
  const std::string tag =
    "BuildBF_f" + std::to_string(branch_factor) + "_d" + std::to_string(dim_mode);

  vtkNew<vtkHyperTreeGrid> vgrid;
  void* fgrid = ftn_htg_grid_create();

  int dims[3] = { 5, 1, 1 };
  if (dim_mode == 2)
  {
    dims[1] = 2; // 4 x 1 = 4 root cells in 2D XY
  }
  else if (dim_mode == 3)
  {
    dims[1] = 2;
    dims[2] = 2; // 4 x 1 x 1 = 4 root cells in 3D XYZ
  }

  vgrid->SetBranchFactor(static_cast<unsigned int>(branch_factor));
  vgrid->SetDimensions(dims[0], dims[1], dims[2]);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_branch_factor(fgrid, branch_factor), tag);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_dimensions(fgrid, dims[0], dims[1], dims[2]), tag);

  oracle::set_rectilinear_coords_both(
    vgrid, fgrid, lin_coords(dims[0]), lin_coords(dims[1]), lin_coords(dims[2]));

  const int k = static_cast<int>(vgrid->GetNumberOfChildren());

  // Construct a combined multi-tree forest bit array for 4 trees:
  // Tree 0: single-leaf root -> `numberOfBits = 0`!
  // Tree 1: 2-level terminal root -> `[1]` (`numberOfBits = 1`)
  // Tree 2: 3-level asymmetric -> root `1` + children `[1, 0, ..., 1]` (`numberOfBits = 1 + k`)
  // Tree 3: 4-level asymmetric -> root `1` (1 bit)
  //         + depth-1: child 0 `1`, child (k-1) `1`, others `0` (k bits, 2 coarse)
  //         + depth-2: 2*k bits with 3 coarse nodes (`numberOfBits = 1 + k + 2*k`)
  std::vector<bool> t0_bits;
  std::vector<bool> t1_bits = { true };
  std::vector<bool> t2_bits(static_cast<size_t>(1 + k), false);
  t2_bits[0] = true;
  t2_bits[1] = true;
  t2_bits[static_cast<size_t>(k)] = true;

  std::vector<bool> t3_bits(static_cast<size_t>(1 + k + 2 * k), false);
  t3_bits[0] = true;
  t3_bits[1] = true;                       // depth-1 child 0 coarse
  t3_bits[static_cast<size_t>(k)] = true;  // depth-1 child (k-1) coarse
  t3_bits[static_cast<size_t>(1 + k)] = true;
  t3_bits[static_cast<size_t>(1 + k + (k / 2))] = true;
  t3_bits[static_cast<size_t>(1 + 2 * k)] = true;

  const std::vector<std::vector<bool>> tree_descs = { t0_bits, t1_bits, t2_bits, t3_bits };

  vtkNew<vtkBitArray> vforest_desc;
  int32_t bst = 0;
  void* fforest_desc = ftn_htg_bit_array_create(0, false, &bst);
  ORACLE_EQ(FTN_HTG_OK, bst, tag);

  // Add a 3-bit dummy prefix so even Tree 1 starts at a non-zero unaligned bit offset!
  append_bits_both(vforest_desc, fforest_desc, { true, false, true });
  int64_t cur_offset = 3;

  for (int64_t tid = 0; tid < 4; ++tid)
  {
    const auto& bits = tree_descs[static_cast<size_t>(tid)];
    const int64_t nbits = static_cast<int64_t>(bits.size());
    append_bits_both(vforest_desc, fforest_desc, bits);

    vtkHyperTree* vtree = vgrid->GetTree(static_cast<vtkIdType>(tid), true);
    int32_t tst = 0;
    void* ftree = ftn_htg_grid_get_tree(fgrid, tid, true, &tst);
    ORACLE_EQ(FTN_HTG_OK, tst, tag << " get_tree tid=" << tid);

    vtree->BuildFromBreadthFirstOrderDescriptor(
      vforest_desc, static_cast<vtkIdType>(nbits), static_cast<vtkIdType>(cur_offset));
    ORACLE_EQ(FTN_HTG_OK,
      ftn_htg_tree_build_from_bf_descriptor(ftree, fforest_desc, nbits, cur_offset),
      tag << " build_from_bf tid=" << tid << " offset=" << cur_offset << " nbits=" << nbits);

    cur_offset += nbits;
  }

  vgrid->InitializeLocalIndexNode();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_initialize_local_index_node(fgrid), tag);

  oracle::CompareOptions opts;
  opts.check_coordinates = true;
  opts.check_global_indexing = true;
  opts.check_validate = true;
  oracle::assert_grids_equivalent(vgrid, fgrid, tag, opts);

  ftn_htg_bit_array_destroy(fforest_desc);
  ftn_htg_grid_destroy(fgrid);
}

void test_compute_bf_descriptor_and_roundtrip(int branch_factor, int dim_mode, bool explicit_idx)
{
  const std::string tag = "ComputeBF_f" + std::to_string(branch_factor) + "_d" +
    std::to_string(dim_mode) + (explicit_idx ? "_explicit" : "_implicit");

  vtkNew<vtkHyperTreeGrid> vgrid;
  void* fgrid = ftn_htg_grid_create();

  int dims[3] = { 4, 1, 1 };
  if (dim_mode == 2)
  {
    dims[1] = 2; // 3 root trees
  }
  else if (dim_mode == 3)
  {
    dims[1] = 2;
    dims[2] = 2; // 3 root trees
  }

  vgrid->SetBranchFactor(static_cast<unsigned int>(branch_factor));
  vgrid->SetDimensions(dims[0], dims[1], dims[2]);
  ftn_htg_grid_set_branch_factor(fgrid, branch_factor);
  ftn_htg_grid_set_dimensions(fgrid, dims[0], dims[1], dims[2]);
  oracle::set_rectilinear_coords_both(
    vgrid, fgrid, lin_coords(dims[0]), lin_coords(dims[1]), lin_coords(dims[2]));

  const int32_t k = static_cast<int32_t>(vgrid->GetNumberOfChildren());

  // Build 3 trees via DFS `NonOrientedCursor::SubdivideLeaf` so local vertex IDs are DFS-ordered
  // (meaning `breadthFirstIdMap` is a genuine non-identity permutation!)
  // Tree 0: single-leaf root
  // Tree 1: deep asymmetric 4-level DFS tree
  // Tree 2: 3-level multi-branch DFS tree
  {
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc0;
    vgrid->InitializeNonOrientedCursor(vc0, 0, true);
    int32_t st = 0;
    void* fc0 = ftn_htg_cursor_create(fgrid, 0, true, &st);
    ftn_htg_cursor_destroy(fc0);
  }
  {
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc1;
    vgrid->InitializeNonOrientedCursor(vc1, 1, true);
    int32_t st = 0;
    void* fc1 = ftn_htg_cursor_create(fgrid, 1, true, &st);
    vc1->SubdivideLeaf();
    ftn_htg_cursor_subdivide_leaf(fc1);
    // Dive down child 0 first (DFS assigns low local IDs deep down child 0!)
    vc1->ToChild(0);
    ftn_htg_cursor_to_child(fc1, 0);
    vc1->SubdivideLeaf();
    ftn_htg_cursor_subdivide_leaf(fc1);
    vc1->ToChild(static_cast<unsigned char>(k - 1));
    ftn_htg_cursor_to_child(fc1, k - 1);
    vc1->SubdivideLeaf();
    ftn_htg_cursor_subdivide_leaf(fc1);
    vc1->ToRoot();
    ftn_htg_cursor_to_root(fc1);
    // Now split child (k - 1) at depth 1 AFTER child 0's grandchildren!
    vc1->ToChild(static_cast<unsigned char>(k - 1));
    ftn_htg_cursor_to_child(fc1, k - 1);
    vc1->SubdivideLeaf();
    ftn_htg_cursor_subdivide_leaf(fc1);
    ftn_htg_cursor_destroy(fc1);
  }
  {
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc2;
    vgrid->InitializeNonOrientedCursor(vc2, 2, true);
    int32_t st = 0;
    void* fc2 = ftn_htg_cursor_create(fgrid, 2, true, &st);
    vc2->SubdivideLeaf();
    ftn_htg_cursor_subdivide_leaf(fc2);
    vc2->ToChild(static_cast<unsigned char>(k / 2));
    ftn_htg_cursor_to_child(fc2, k / 2);
    vc2->SubdivideLeaf();
    ftn_htg_cursor_subdivide_leaf(fc2);
    ftn_htg_cursor_destroy(fc2);
  }

  if (!explicit_idx)
  {
    vgrid->InitializeLocalIndexNode();
    ftn_htg_grid_initialize_local_index_node(fgrid);
  }
  else
  {
    // Non-contiguous explicit global IDs (stride 3 starting at 11)
    int64_t next_gid = 11;
    for (int64_t tid = 0; tid < 3; ++tid)
    {
      vtkHyperTree* vt = vgrid->GetTree(static_cast<vtkIdType>(tid), false);
      int32_t st = 0;
      void* ft = ftn_htg_grid_get_tree(fgrid, tid, false, &st);
      const int64_t nv = ftn_htg_tree_get_number_of_vertices(ft);
      for (int64_t lv = 0; lv < nv; ++lv)
      {
        vt->SetGlobalIndexFromLocal(static_cast<vtkIdType>(lv), static_cast<vtkIdType>(next_gid));
        ORACLE_EQ(FTN_HTG_OK, ftn_htg_tree_set_global_index_from_local(ft, lv, next_gid), tag);
        next_gid += 3;
      }
    }
  }

  const int64_t max_gid = static_cast<int64_t>(vgrid->GetGlobalNodeIndexMax());

  // Build an optional `inputMask` that masks:
  // - In Tree 1: child 0's grandchild coarse node (reducing Tree 1's effective maxDepth from 4->3!)
  // - In Tree 2: child (k/2) coarse node (reducing Tree 2's effective maxDepth from 3->2!)
  vtkNew<vtkBitArray> vmask;
  vmask->SetNumberOfTuples(static_cast<vtkIdType>(max_gid + 1));
  for (int64_t i = 0; i <= max_gid; ++i)
  {
    vmask->SetValue(static_cast<vtkIdType>(i), 0);
  }
  int32_t mst = 0;
  void* fmask = ftn_htg_bit_array_create(max_gid + 1, false, &mst);

  {
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc;
    vgrid->InitializeNonOrientedCursor(vc, 1, false);
    vc->ToChild(0);
    vc->ToChild(static_cast<unsigned char>(k - 1));
    const int64_t gid1 = static_cast<int64_t>(vc->GetGlobalNodeIndex());
    vmask->SetValue(static_cast<vtkIdType>(gid1), 1);
    ftn_htg_bit_array_set_value(fmask, gid1, true);

    vgrid->InitializeNonOrientedCursor(vc, 2, false);
    vc->ToChild(static_cast<unsigned char>(k / 2));
    const int64_t gid2 = static_cast<int64_t>(vc->GetGlobalNodeIndex());
    vmask->SetValue(static_cast<vtkIdType>(gid2), 1);
    ftn_htg_bit_array_set_value(fmask, gid2, true);
  }

  struct ExportConfig
  {
    const char* sub_name;
    bool use_mask;
    bool use_limiter;
    int32_t limiter;
  };

  const ExportConfig configs[] = {
    { "unmasked_unlimited", false, false, 0 },
    { "masked_unlimited", true, false, 0 },
    { "unmasked_limit0", false, true, 0 },
    { "unmasked_limit1", false, true, 1 },
    { "unmasked_limit2", false, true, 2 },
    { "masked_limit2", true, true, 2 },
  };

  for (const auto& cfg : configs)
  {
    const std::string ctag = tag + "/" + cfg.sub_name;

    // Concatenate across all 3 trees in forest order into shared output buffers!
    vtkNew<vtkBitArray> v_desc;
    vtkNew<vtkTypeInt64Array> v_vpd;
    vtkNew<vtkIdList> v_bfid;

    int32_t dst = 0;
    void* f_desc = ftn_htg_bit_array_create(0, false, &dst);
    void* f_vpd = ftn_htg_i64_vec_create();
    void* f_bfid = ftn_htg_i64_vec_create();

    std::vector<int64_t> per_tree_bit_counts(3, 0);
    for (int64_t tid = 0; tid < 3; ++tid)
    {
      const int64_t before_bits = static_cast<int64_t>(v_desc->GetNumberOfValues());
      vtkHyperTree* vt = vgrid->GetTree(static_cast<vtkIdType>(tid), false);
      int32_t tst = 0;
      void* ft = ftn_htg_grid_get_tree(fgrid, tid, false, &tst);

      const unsigned int v_lim =
        cfg.use_limiter ? static_cast<unsigned int>(cfg.limiter)
                        : std::numeric_limits<unsigned int>::max();
      vt->ComputeBreadthFirstOrderDescriptor(
        v_lim, cfg.use_mask ? vmask.GetPointer() : nullptr, v_vpd, v_desc, v_bfid);

      ORACLE_EQ(FTN_HTG_OK,
        ftn_htg_tree_compute_bf_descriptor_append(ft, f_desc, f_vpd, f_bfid,
          cfg.use_mask ? fmask : nullptr, cfg.use_limiter, cfg.limiter),
        ctag << " compute_bf tid=" << tid);

      per_tree_bit_counts[static_cast<size_t>(tid)] =
        static_cast<int64_t>(v_desc->GetNumberOfValues()) - before_bits;
    }

    // 1. Compare concatenated `descriptor` bit count & raw MSB-first bytes
    const int64_t v_nbits = static_cast<int64_t>(v_desc->GetNumberOfValues());
    const int64_t f_nbits = ftn_htg_bit_array_get_number_of_values(f_desc);
    ORACLE_EQ(v_nbits, f_nbits, ctag << " descriptor num_bits");
    if (v_nbits > 0)
    {
      const int64_t nbytes = (v_nbits + 7) / 8;
      const uint8_t* vraw = reinterpret_cast<const uint8_t*>(v_desc->GetPointer(0));
      const uint8_t* fraw = nullptr;
      int64_t fnbytes = 0;
      ftn_htg_bit_array_get_raw_bytes(f_desc, &fraw, &fnbytes);
      ORACLE_EQ(nbytes, fnbytes, ctag << " descriptor nbytes");
      ORACLE_ASSERT(std::memcmp(vraw, fraw, static_cast<size_t>(nbytes)) == 0,
        ctag << " descriptor raw bytes mismatch\n  VTK: [" << oracle::bytes_to_hex(vraw, nbytes)
             << "]\n  FTN: [" << oracle::bytes_to_hex(fraw, fnbytes) << "]");
    }

    // 2. Compare concatenated `numberOfVerticesPerDepth` element-for-element
    const int64_t* fvpd_ptr = nullptr;
    int64_t fvpd_len = 0;
    ftn_htg_i64_vec_get_view(f_vpd, &fvpd_ptr, &fvpd_len);
    ORACLE_EQ(static_cast<int64_t>(v_vpd->GetNumberOfValues()), fvpd_len, ctag << " vpd len");
    for (int64_t i = 0; i < fvpd_len; ++i)
    {
      ORACLE_EQ(static_cast<int64_t>(v_vpd->GetValue(static_cast<vtkIdType>(i))), fvpd_ptr[i],
        ctag << " vpd[" << i << "]");
    }

    // 3. Compare concatenated `breadthFirstIdMap` element-for-element
    const int64_t* fbfid_ptr = nullptr;
    int64_t fbfid_len = 0;
    ftn_htg_i64_vec_get_view(f_bfid, &fbfid_ptr, &fbfid_len);
    ORACLE_EQ(static_cast<int64_t>(v_bfid->GetNumberOfIds()), fbfid_len, ctag << " bfid len");
    for (int64_t i = 0; i < fbfid_len; ++i)
    {
      ORACLE_EQ(static_cast<int64_t>(v_bfid->GetId(static_cast<vtkIdType>(i))), fbfid_ptr[i],
        ctag << " bfid[" << i << "]");
    }

    // 4. Round-trip rebuild a fresh VTK & Fortran forest from the concatenated `descriptor`!
    vtkNew<vtkHyperTreeGrid> v_rebuilt;
    void* f_rebuilt = ftn_htg_grid_create();
    v_rebuilt->SetBranchFactor(static_cast<unsigned int>(branch_factor));
    v_rebuilt->SetDimensions(dims[0], dims[1], dims[2]);
    ftn_htg_grid_set_branch_factor(f_rebuilt, branch_factor);
    ftn_htg_grid_set_dimensions(f_rebuilt, dims[0], dims[1], dims[2]);
    oracle::set_rectilinear_coords_both(
      v_rebuilt, f_rebuilt, lin_coords(dims[0]), lin_coords(dims[1]), lin_coords(dims[2]));

    int64_t start_off = 0;
    for (int64_t tid = 0; tid < 3; ++tid)
    {
      const int64_t nb = per_tree_bit_counts[static_cast<size_t>(tid)];
      vtkHyperTree* vrt = v_rebuilt->GetTree(static_cast<vtkIdType>(tid), true);
      int32_t rst = 0;
      void* frt = ftn_htg_grid_get_tree(f_rebuilt, tid, true, &rst);
      vrt->BuildFromBreadthFirstOrderDescriptor(
        v_desc, static_cast<vtkIdType>(nb), static_cast<vtkIdType>(start_off));
      ORACLE_EQ(FTN_HTG_OK,
        ftn_htg_tree_build_from_bf_descriptor(frt, f_desc, nb, start_off),
        ctag << " roundtrip rebuild tid=" << tid);
      start_off += nb;
    }
    v_rebuilt->InitializeLocalIndexNode();
    ftn_htg_grid_initialize_local_index_node(f_rebuilt);

    oracle::CompareOptions ropts;
    ropts.check_coordinates = true;
    ropts.check_global_indexing = true;
    ropts.check_validate = true;
    oracle::assert_grids_equivalent(v_rebuilt, f_rebuilt, ctag + "/roundtrip", ropts);

    ftn_htg_grid_destroy(f_rebuilt);
    ftn_htg_i64_vec_destroy(f_bfid);
    ftn_htg_i64_vec_destroy(f_vpd);
    ftn_htg_bit_array_destroy(f_desc);
  }

  ftn_htg_bit_array_destroy(fmask);
  ftn_htg_grid_destroy(fgrid);
}

} // namespace

int main()
{
  oracle::verify_vtk_oracle_version();

  for (int bf : { 2, 3 })
  {
    for (int d : { 1, 2, 3 })
    {
      test_build_from_bf_descriptor_all_combinations(bf, d);
      test_compute_bf_descriptor_and_roundtrip(bf, d, false);
      test_compute_bf_descriptor_and_roundtrip(bf, d, true);
    }
  }

  std::cout << "[PASS] test_oracle_bf_descriptors verified all BF descriptor build/compute/mask/"
               "limiter/roundtrip combinations!"
            << std::endl;
  return 0;
}
