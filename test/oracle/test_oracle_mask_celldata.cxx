#include "oracle_compare.hpp"

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkTypeInt32Array.h>
#include <vtkTypeInt64Array.h>

#include <cstring>
#include <iostream>
#include <vector>

namespace {

void test_bit_array_msb_parity()
{
  const std::string tag = "BitArrayMSBParity";
  vtkNew<vtkBitArray> vbits;
  int32_t st = 0;
  void* fbits = ftn_htg_bit_array_create(0, false, &st);
  ORACLE_EQ(FTN_HTG_OK, st, tag << " create");

  const int64_t test_ids[] = { 0, 7, 8, 13, 15, 16, 31, 32, 63, 64, 99, 127, 128, 130 };
  for (size_t i = 0; i < sizeof(test_ids) / sizeof(test_ids[0]); ++i)
  {
    const bool bit_val = (i % 2 == 0);
    vbits->InsertValue(static_cast<vtkIdType>(test_ids[i]), bit_val ? 1 : 0);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_bit_array_insert_value(fbits, test_ids[i], bit_val), tag);

    const int64_t v_n = static_cast<int64_t>(vbits->GetNumberOfValues());
    const int64_t f_n = ftn_htg_bit_array_get_number_of_values(fbits);
    ORACLE_EQ(v_n, f_n, tag << " num_bits after insert " << test_ids[i]);

    const int64_t nbytes = (v_n + 7) / 8;
    const uint8_t* vraw = reinterpret_cast<const uint8_t*>(vbits->GetPointer(0));
    const uint8_t* fraw = nullptr;
    int64_t fnbytes = 0;
    ftn_htg_bit_array_get_raw_bytes(fbits, &fraw, &fnbytes);
    ORACLE_EQ(nbytes, fnbytes, tag << " active_bytes");
    ORACLE_ASSERT(std::memcmp(vraw, fraw, static_cast<size_t>(nbytes)) == 0,
      tag << " raw MSB bytes mismatch after insert id=" << test_ids[i] << "\n  VTK: ["
          << oracle::bytes_to_hex(vraw, nbytes) << "]\n  FTN: ["
          << oracle::bytes_to_hex(fraw, fnbytes) << "]");
  }

  // Test InsertNextValue
  for (int k = 0; k < 19; ++k)
  {
    const bool bv = ((k * 3 + 1) % 4 == 0);
    const vtkIdType vid = vbits->InsertNextValue(bv ? 1 : 0);
    int64_t fid = -1;
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_bit_array_insert_next_value(fbits, bv, &fid), tag);
    ORACLE_EQ(static_cast<int64_t>(vid), fid, tag << " inserted_id");
  }

  // Test shrinking via SetNumberOfValues (verifies trailing slack bits in last byte are cleared!)
  vbits->SetNumberOfValues(67);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_bit_array_resize(fbits, 67), tag << " shrink to 67");
  {
    const int64_t nbytes = (67 + 7) / 8;
    const uint8_t* vraw = reinterpret_cast<const uint8_t*>(vbits->GetPointer(0));
    const uint8_t* fraw = nullptr;
    int64_t fnbytes = 0;
    ftn_htg_bit_array_get_raw_bytes(fbits, &fraw, &fnbytes);
    ORACLE_EQ(nbytes, fnbytes, tag << " shrunk nbytes");
    ORACLE_ASSERT(std::memcmp(vraw, fraw, static_cast<size_t>(nbytes)) == 0,
      tag << " shrunk raw bytes mismatch\n  VTK: [" << oracle::bytes_to_hex(vraw, nbytes)
          << "]\n  FTN: [" << oracle::bytes_to_hex(fraw, fnbytes) << "]");
  }

  ftn_htg_bit_array_destroy(fbits);
}

void test_mask_pure_mask_and_material_interfaces(int branch_factor)
{
  const std::string tag = "MaskPureMaskInterfaces_f" + std::to_string(branch_factor);

  vtkNew<vtkHyperTreeGrid> vgrid;
  void* fgrid = ftn_htg_grid_create();

  vgrid->SetBranchFactor(static_cast<unsigned int>(branch_factor));
  vgrid->SetDimensions(5, 2, 2); // 4 x 1 x 1 = 4 root trees in 3D (or 5x3x2 -> 4x2x1 = 8 trees)
  vgrid->SetDimensions(5, 3, 2); // 4 * 2 * 1 = 8 root trees in 3D

  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_branch_factor(fgrid, branch_factor), tag);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_dimensions(fgrid, 5, 3, 2), tag);

  oracle::set_rectilinear_coords_both(
    vgrid, fgrid, { 0.0, 1.0, 2.0, 3.0, 4.0 }, { 0.0, 1.0, 2.0 }, { 0.0, 1.0 });

  // Build 4 trees:
  // Tree 0: refined 2 levels, purely unmasked & pure material
  // Tree 1: refined 2 levels, partial leaf masked
  // Tree 2: refined 2 levels, one coarse internal child masked
  // Tree 3: single-leaf root, masked
  // Tree 5: refined 1 level, unmasked in `Mask` but mixed material interface on child 1!
  const int64_t tree_ids[5] = { 0, 1, 2, 3, 5 };
  for (int idx = 0; idx < 5; ++idx)
  {
    const int64_t tid = tree_ids[idx];
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc;
    vgrid->InitializeNonOrientedCursor(vc, static_cast<vtkIdType>(tid), true);
    int32_t cst = 0;
    void* fc = ftn_htg_cursor_create(fgrid, tid, true, &cst);
    ORACLE_EQ(FTN_HTG_OK, cst, tag);

    if (tid != 3)
    {
      vc->SubdivideLeaf();
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fc), tag);
      if (tid != 5)
      {
        const int32_t sub_ch = (tid == 2) ? 0 : 1;
        vc->ToChild(static_cast<unsigned char>(sub_ch));
        ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fc, sub_ch), tag);
        vc->SubdivideLeaf();
        ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fc), tag);
        vc->ToRoot();
        ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_root(fc), tag);
      }
    }
    ftn_htg_cursor_destroy(fc);
  }

  // Assign implicit global indices across forest
  vgrid->InitializeLocalIndexNode();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_initialize_local_index_node(fgrid), tag);

  const int64_t total_cells = static_cast<int64_t>(vgrid->GetGlobalNodeIndexMax() + 1);

  // 0-tuple SetMask parity check: both VTK and Fortran report HasMask() == false on empty bit array
  {
    vtkNew<vtkBitArray> vempty_mask;
    vgrid->SetMask(vempty_mask);
    int32_t est = 0;
    void* fempty_mask = ftn_htg_bit_array_create(0, false, &est);
    ORACLE_EQ(FTN_HTG_OK, est, tag << " create empty mask");
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_mask(fgrid, fempty_mask), tag << " set empty mask");
    ftn_htg_bit_array_destroy(fempty_mask);

    ORACLE_EQ(vgrid->HasMask(), ftn_htg_grid_has_mask(fgrid), tag << " 0-tuple HasMask()");
    ORACLE_EQ(false, ftn_htg_grid_has_mask(fgrid), tag << " 0-tuple HasMask() is false");
    char verr[256] = { 0 };
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_validate(fgrid, verr, sizeof(verr)),
      tag << " validate 0-tuple mask: " << verr);
  }

  // Attach initial all-false Mask sized to `total_cells`
  vtkNew<vtkBitArray> vmask;
  vmask->SetNumberOfTuples(static_cast<vtkIdType>(total_cells));
  for (int64_t i = 0; i < total_cells; ++i)
  {
    vmask->SetValue(static_cast<vtkIdType>(i), 0);
  }
  vgrid->SetMask(vmask);

  int32_t bst = 0;
  void* fmask = ftn_htg_bit_array_create(total_cells, false, &bst);
  ORACLE_EQ(FTN_HTG_OK, bst, tag);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_mask(fgrid, fmask), tag);
  ftn_htg_bit_array_destroy(fmask);

  // Set masks via NonOrientedCursor!
  // 1) Tree 1 -> child 1 -> grandchild (nch - 1) leaf masked
  {
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc;
    vgrid->InitializeNonOrientedCursor(vc, 1, false);
    int32_t cst = 0;
    void* fc = ftn_htg_cursor_create(fgrid, 1, false, &cst);
    const int32_t nch = static_cast<int32_t>(vc->GetNumberOfChildren());
    vc->ToChild(1);
    ftn_htg_cursor_to_child(fc, 1);
    vc->ToChild(static_cast<unsigned char>(nch - 1));
    ftn_htg_cursor_to_child(fc, nch - 1);
    vc->SetMask(true);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_set_mask(fc, true), tag << " mask leaf");
    ftn_htg_cursor_destroy(fc);
  }

  // 2) Tree 2 -> child 0 (coarse internal node!) masked
  {
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc;
    vgrid->InitializeNonOrientedCursor(vc, 2, false);
    int32_t cst = 0;
    void* fc = ftn_htg_cursor_create(fgrid, 2, false, &cst);
    vc->ToChild(0);
    ftn_htg_cursor_to_child(fc, 0);
    vc->SetMask(true);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_set_mask(fc, true), tag << " mask coarse");
    ftn_htg_cursor_destroy(fc);
  }

  // 3) Tree 3 -> root leaf masked
  {
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc;
    vgrid->InitializeNonOrientedCursor(vc, 3, false);
    int32_t cst = 0;
    void* fc = ftn_htg_cursor_create(fgrid, 3, false, &cst);
    vc->SetMask(true);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_set_mask(fc, true), tag << " mask root");
    ftn_htg_cursor_destroy(fc);
  }

  // Verify Mask & PureMask BEFORE enabling Material Interfaces!
  oracle::CompareOptions opts;
  opts.check_coordinates = true;
  opts.check_global_indexing = true;
  opts.check_masks = true;
  opts.check_pure_mask = true;
  opts.check_validate = true;
  oracle::assert_grids_equivalent(vgrid, fgrid, tag + "/mask_only", opts);

  // Now attach CellData arrays (Pressure r64, MaterialId i64, Zone i32, Temp r32, Flag bit)
  // AND 3-component Material Interface arrays (`InterfaceNormals` & `InterfaceIntercepts`)!
  const char* n_press = "Pressure";
  const char* n_mat = "MaterialId";
  const char* n_zone = "Zone";
  const char* n_temp = "Temp";
  const char* n_flag = "FlagBit";
  const char* n_norm = "InterfaceNormals";
  const char* n_intr = "InterfaceIntercepts";

  vtkNew<vtkDoubleArray> v_press;
  v_press->SetName(n_press);
  v_press->SetNumberOfComponents(1);
  v_press->SetNumberOfTuples(static_cast<vtkIdType>(total_cells));

  vtkNew<vtkTypeInt64Array> v_mat;
  v_mat->SetName(n_mat);
  v_mat->SetNumberOfComponents(1);
  v_mat->SetNumberOfTuples(static_cast<vtkIdType>(total_cells));

  vtkNew<vtkDoubleArray> v_norm;
  v_norm->SetName(n_norm);
  v_norm->SetNumberOfComponents(3);
  v_norm->SetNumberOfTuples(static_cast<vtkIdType>(total_cells));

  vtkNew<vtkDoubleArray> v_intr;
  v_intr->SetName(n_intr);
  v_intr->SetNumberOfComponents(3);
  v_intr->SetNumberOfTuples(static_cast<vtkIdType>(total_cells));

  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_cell_data_add_array(
      fgrid, n_press, static_cast<int32_t>(std::strlen(n_press)), FTN_HTG_TYPE_REAL64, 1,
      total_cells),
    tag);
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_cell_data_add_array(
      fgrid, n_mat, static_cast<int32_t>(std::strlen(n_mat)), FTN_HTG_TYPE_INT64, 1, total_cells),
    tag);
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_cell_data_add_array(
      fgrid, n_zone, static_cast<int32_t>(std::strlen(n_zone)), FTN_HTG_TYPE_INT32, 1, total_cells),
    tag);
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_cell_data_add_array(
      fgrid, n_temp, static_cast<int32_t>(std::strlen(n_temp)), FTN_HTG_TYPE_REAL32, 1,
      total_cells),
    tag);
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_cell_data_add_array(
      fgrid, n_flag, static_cast<int32_t>(std::strlen(n_flag)), FTN_HTG_TYPE_BIT, 2, total_cells),
    tag);
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_cell_data_add_array(
      fgrid, n_norm, static_cast<int32_t>(std::strlen(n_norm)), FTN_HTG_TYPE_REAL64, 3,
      total_cells),
    tag);
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_cell_data_add_array(
      fgrid, n_intr, static_cast<int32_t>(std::strlen(n_intr)), FTN_HTG_TYPE_REAL64, 3,
      total_cells),
    tag);

  // Initialize all cells as pure material (`intercepts[2] = 2.0`), then mark Tree 5 child 1 mixed!
  for (int64_t g = 0; g < total_cells; ++g)
  {
    const double pval = 101.325 + 0.5 * static_cast<double>(g);
    const int64_t mval = 1000 + g;
    const int32_t zval = static_cast<int32_t>(g % 7);
    const float tval = static_cast<float>(300.0 + g);
    const double nrm[3] = { 1.0, 0.0, 0.0 };
    const double itc[3] = { 0.25, 0.75, 2.0 }; // type 2.0 => pure cell

    v_press->SetTuple1(static_cast<vtkIdType>(g), pval);
    v_mat->SetValue(static_cast<vtkIdType>(g), mval);
    v_norm->SetTuple(static_cast<vtkIdType>(g), nrm);
    v_intr->SetTuple(static_cast<vtkIdType>(g), itc);

    ftn_htg_grid_cell_data_set_tuple_r64(
      fgrid, n_press, static_cast<int32_t>(std::strlen(n_press)), g, &pval, 1);
    ftn_htg_grid_cell_data_set_tuple_i64(
      fgrid, n_mat, static_cast<int32_t>(std::strlen(n_mat)), g, &mval, 1);
    ftn_htg_grid_cell_data_set_tuple_i32(
      fgrid, n_zone, static_cast<int32_t>(std::strlen(n_zone)), g, &zval, 1);
    ftn_htg_grid_cell_data_set_tuple_r32(
      fgrid, n_temp, static_cast<int32_t>(std::strlen(n_temp)), g, &tval, 1);
    ftn_htg_grid_cell_data_set_component_bit(
      fgrid, n_flag, static_cast<int32_t>(std::strlen(n_flag)), g, 0, (g % 2 == 0));
    ftn_htg_grid_cell_data_set_component_bit(
      fgrid, n_flag, static_cast<int32_t>(std::strlen(n_flag)), g, 1, (g % 3 == 0));
    ftn_htg_grid_cell_data_set_tuple_r64(
      fgrid, n_norm, static_cast<int32_t>(std::strlen(n_norm)), g, nrm, 3);
    ftn_htg_grid_cell_data_set_tuple_r64(
      fgrid, n_intr, static_cast<int32_t>(std::strlen(n_intr)), g, itc, 3);
  }

  // Mark Tree 5 -> child 1 as mixed material (`intercepts[2] = 1.0 < 2.0`)!
  {
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc;
    vgrid->InitializeNonOrientedCursor(vc, 5, false);
    int32_t cst = 0;
    void* fc = ftn_htg_cursor_create(fgrid, 5, false, &cst);
    vc->ToChild(1);
    ftn_htg_cursor_to_child(fc, 1);

    const int64_t mixed_gid = static_cast<int64_t>(vc->GetGlobalNodeIndex());
    ORACLE_EQ(mixed_gid, ftn_htg_cursor_get_global_node_index(fc, &cst), tag);

    const double mixed_nrm[3] = { 0.6, 0.8, 0.0 };
    const double mixed_itc[3] = { 0.125, 0.0, 1.0 }; // type 1.0 (< 2) => mixed cell!
    v_norm->SetTuple(static_cast<vtkIdType>(mixed_gid), mixed_nrm);
    v_intr->SetTuple(static_cast<vtkIdType>(mixed_gid), mixed_itc);
    ftn_htg_grid_cell_data_set_tuple_r64(
      fgrid, n_norm, static_cast<int32_t>(std::strlen(n_norm)), mixed_gid, mixed_nrm, 3);
    ftn_htg_grid_cell_data_set_tuple_r64(
      fgrid, n_intr, static_cast<int32_t>(std::strlen(n_intr)), mixed_gid, mixed_itc, 3);

    ftn_htg_cursor_destroy(fc);
  }

  vgrid->GetCellData()->AddArray(v_press);
  vgrid->GetCellData()->AddArray(v_mat);
  vgrid->GetCellData()->AddArray(v_norm);
  vgrid->GetCellData()->AddArray(v_intr);
  vgrid->SetHasInterface(true);
  vgrid->SetInterfaceNormalsName(n_norm);
  vgrid->SetInterfaceInterceptsName(n_intr);
  // Note: In VTK, SetMask calls CleanPureMask() so GetPureMask() recomputes with HasInterface!
  vgrid->SetMask(vmask);

  ftn_htg_grid_set_has_interface(fgrid, true);
  ftn_htg_grid_set_interface_normals_name(
    fgrid, n_norm, static_cast<int32_t>(std::strlen(n_norm)));
  ftn_htg_grid_set_interface_intercepts_name(
    fgrid, n_intr, static_cast<int32_t>(std::strlen(n_intr)));

  // Verify PureMask invalidation & exact parity with Material Interfaces active!
  oracle::assert_grids_equivalent(vgrid, fgrid, tag + "/with_interface", opts);

  // Verify Fortran PureMask automatic invalidation upon cursor SetMask mutations (without calling
  // ftn_htg_grid_set_mask!):
  // (a) Masking a single grandchild in pure Tree 0 -> coarse parent flips 0 (pure) -> 1 (impure);
  // (b) Unmasking the sole masked grandchild in Tree 1 -> coarse parent flips 1 (impure) -> 0 (pure).
  {
    int32_t cst = 0;
    // (a) Tree 0 -> child 1 starts pure (0)
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc0;
    vgrid->InitializeNonOrientedCursor(vc0, 0, false);
    void* fc0 = ftn_htg_cursor_create(fgrid, 0, false, &cst);
    vc0->ToChild(1);
    ftn_htg_cursor_to_child(fc0, 1);
    const int64_t t0_coarse_gid = static_cast<int64_t>(vc0->GetGlobalNodeIndex());
    ORACLE_EQ(false, vgrid->GetPureMask()->GetValue(static_cast<vtkIdType>(t0_coarse_gid)) != 0,
      tag << " VTK Tree 0 child 1 initially pure");
    ORACLE_EQ(false, ftn_htg_grid_get_pure_mask_bit(fgrid, t0_coarse_gid, &cst),
      tag << " FTN Tree 0 child 1 initially pure");
    vc0->ToChild(0);
    ftn_htg_cursor_to_child(fc0, 0);
    vc0->SetMask(true);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_set_mask(fc0, true), tag << " mask Tree 0 grandchild");
    ftn_htg_cursor_destroy(fc0);

    // (b) Tree 1 -> child 1 starts impure (1); unmasking grandchild (nch - 1) restores pure (0)
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc1;
    vgrid->InitializeNonOrientedCursor(vc1, 1, false);
    void* fc1 = ftn_htg_cursor_create(fgrid, 1, false, &cst);
    const int32_t nch = static_cast<int32_t>(vc1->GetNumberOfChildren());
    vc1->ToChild(1);
    ftn_htg_cursor_to_child(fc1, 1);
    const int64_t t1_coarse_gid = static_cast<int64_t>(vc1->GetGlobalNodeIndex());
    ORACLE_EQ(true, vgrid->GetPureMask()->GetValue(static_cast<vtkIdType>(t1_coarse_gid)) != 0,
      tag << " VTK Tree 1 child 1 initially impure");
    ORACLE_EQ(true, ftn_htg_grid_get_pure_mask_bit(fgrid, t1_coarse_gid, &cst),
      tag << " FTN Tree 1 child 1 initially impure");
    vc1->ToChild(static_cast<unsigned char>(nch - 1));
    ftn_htg_cursor_to_child(fc1, nch - 1);
    vc1->SetMask(false);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_set_mask(fc1, false), tag << " unmask Tree 1 grandchild");
    ftn_htg_cursor_destroy(fc1);

    // Reset VTK's Mask pointer so VTK calls CleanPureMask() after in-place cursor mutations;
    // Fortran fgrid gets NO grid_set_mask call, proving automatic pure_mask_valid invalidation!
    vgrid->SetMask(nullptr);
    vgrid->SetMask(vmask);

    ORACLE_EQ(true, ftn_htg_grid_get_pure_mask_bit(fgrid, t0_coarse_gid, &cst),
      tag << " Tree 0 child 1 flips to impure (1) in PureMask");
    ORACLE_EQ(false, ftn_htg_grid_get_pure_mask_bit(fgrid, t1_coarse_gid, &cst),
      tag << " Tree 1 child 1 flips to pure (0) in PureMask");
    oracle::assert_grids_equivalent(vgrid, fgrid, tag + "/pure_mask_invalidation", opts);
  }

  // Verify CellData zero-copy view & tuple values match VTK CellData at every global index
  const double* f_intr_view = nullptr;
  int32_t f_intr_ncomp = 0;
  int64_t f_intr_ntup = 0;
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_cell_data_get_view_r64(fgrid, n_intr, static_cast<int32_t>(std::strlen(n_intr)),
      &f_intr_view, &f_intr_ncomp, &f_intr_ntup),
    tag << " get_view_r64");
  ORACLE_EQ(3, f_intr_ncomp, tag << " intercepts ncomp");
  ORACLE_ASSERT(f_intr_ntup >= total_cells, tag << " intercepts ntuples");
  for (int64_t g = 0; g < total_cells; ++g)
  {
    double vtup[3] = { 0.0, 0.0, 0.0 };
    v_intr->GetTuple(static_cast<vtkIdType>(g), vtup);
    for (int c = 0; c < 3; ++c)
    {
      ORACLE_ASSERT(std::abs(vtup[c] - f_intr_view[g * 3 + c]) <= 1e-15,
        tag << " intercepts mismatch at g=" << g << " comp=" << c);
    }
  }

  // Test dynamic CellData capacity growth during further subdivision!
  {
    vtkNew<vtkHyperTreeGridNonOrientedCursor> vc;
    vgrid->InitializeNonOrientedCursor(vc, 0, false);
    int32_t cst = 0;
    void* fc = ftn_htg_cursor_create(fgrid, 0, false, &cst);
    vc->ToChild(0);
    ftn_htg_cursor_to_child(fc, 0);
    // Clear mask during structural expansion, subdivide leaf, re-index & grow mask/celldata
    vgrid->SetMask(nullptr);
    ftn_htg_grid_clear_mask(fgrid);
    vc->SubdivideLeaf();
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fc), tag << " dynamic subdivide");
    ftn_htg_cursor_destroy(fc);

    vgrid->InitializeLocalIndexNode();
    ftn_htg_grid_initialize_local_index_node(fgrid);
    const int64_t grown_cells = static_cast<int64_t>(vgrid->GetGlobalNodeIndexMax() + 1);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_cell_data_ensure_tuples(fgrid, grown_cells), tag);

    opts.check_masks = false;
    opts.check_pure_mask = false;
    opts.check_validate = true;
    oracle::assert_grids_equivalent(vgrid, fgrid, tag + "/after_dynamic_growth", opts);

    // Verify DepthLimiter cursor IsLeaf() parity on VTK vs Fortran
    const unsigned int orig_vlim = vgrid->GetDepthLimiter();
    const int32_t orig_flim = ftn_htg_grid_get_depth_limiter(fgrid);
    vgrid->SetDepthLimiter(1);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_depth_limiter(fgrid, 1), tag << " SetDepthLimiter(1)");
    {
      vtkNew<vtkHyperTreeGridNonOrientedCursor> vc_lim;
      vgrid->InitializeNonOrientedCursor(vc_lim, 0, false);
      int32_t lst = 0;
      void* fc_lim = ftn_htg_cursor_create(fgrid, 0, false, &lst);
      ORACLE_EQ(vc_lim->IsLeaf(), ftn_htg_cursor_is_leaf(fc_lim, &lst), tag << " root IsLeaf dlim=1");
      ORACLE_EQ(false, vc_lim->IsLeaf(), tag << " root IsLeaf false");
      vc_lim->ToChild(0);
      ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fc_lim, 0), tag << " ToChild(0) at dlim=1");
      ORACLE_EQ(
        vc_lim->IsLeaf(), ftn_htg_cursor_is_leaf(fc_lim, &lst), tag << " level-1 IsLeaf dlim=1");
      ORACLE_EQ(true, vc_lim->IsLeaf(), tag << " level-1 treated as leaf at DepthLimiter=1");
      ftn_htg_cursor_destroy(fc_lim);
    }
    vgrid->SetDepthLimiter(orig_vlim);
    ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_depth_limiter(fgrid, orig_flim), tag);
  }

  ftn_htg_grid_destroy(fgrid);
}

} // namespace

int main()
{
  oracle::verify_vtk_oracle_version();
  test_bit_array_msb_parity();
  for (int bf : { 2, 3 })
  {
    test_mask_pure_mask_and_material_interfaces(bf);
  }
  std::cout << "[PASS] test_oracle_mask_celldata verified MSB bit arrays, Mask, PureMask & "
               "CellData/Interfaces!"
            << std::endl;
  return 0;
}
