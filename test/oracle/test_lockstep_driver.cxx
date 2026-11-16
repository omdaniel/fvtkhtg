#include "oracle_compare.hpp"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>

#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE_ATTR __attribute__((noinline))
#else
#define NOINLINE_ATTR
#endif

struct LockstepSnapshotView
{
  int32_t step_id;
  const char* phase_name;
  int64_t tree_id;
  int64_t vtk_num_levels;
  int64_t ftn_num_levels;
  int64_t vtk_num_vertices;
  int64_t ftn_num_vertices;
  int64_t vtk_num_nodes;
  int64_t ftn_num_nodes;
  const unsigned int* vtk_ptec_ptr;
  size_t vtk_ptec_len;
  const int64_t* ftn_ptec_ptr;
  int64_t ftn_ptec_len;
  const int64_t* ftn_vpd_ptr;
  int64_t ftn_vpd_len;
  const uint8_t* vtk_bits_ptr;
  int64_t vtk_bits_nbytes;
  const uint8_t* ftn_bits_ptr;
  int64_t ftn_bits_nbytes;
};

extern "C" NOINLINE_ATTR void lockstep_on_step(const LockstepSnapshotView* snap)
{
  // Prevent compiler from optimizing away `snap` in `-O0` debug builds
  volatile int32_t sid = snap ? snap->step_id : -1;
  (void)sid;
}

static void emit_tree_step_snapshot(int32_t step_id, const char* phase, vtkHyperTree* vt,
  void* ft, const uint8_t* vbits = nullptr, int64_t vnbytes = 0, const uint8_t* fbits = nullptr,
  int64_t fnbytes = 0)
{
  LockstepSnapshotView snap{};
  snap.step_id = step_id;
  snap.phase_name = phase;
  if (vt && ft)
  {
    snap.tree_id = static_cast<int64_t>(vt->GetTreeIndex());
    snap.vtk_num_levels = static_cast<int64_t>(vt->GetNumberOfLevels());
    snap.ftn_num_levels = ftn_htg_tree_get_number_of_levels(ft);
    snap.vtk_num_vertices = static_cast<int64_t>(vt->GetNumberOfVertices());
    snap.ftn_num_vertices = ftn_htg_tree_get_number_of_vertices(ft);
    snap.vtk_num_nodes = static_cast<int64_t>(vt->GetNumberOfNodes());
    snap.ftn_num_nodes = ftn_htg_tree_get_number_of_nodes(ft);
    snap.vtk_ptec_ptr = vt->GetElderChildIndexArray(snap.vtk_ptec_len);
    ftn_htg_tree_get_parent_to_elder_child_view(ft, &snap.ftn_ptec_ptr, &snap.ftn_ptec_len);
    ftn_htg_tree_get_vertices_per_depth_view(ft, &snap.ftn_vpd_ptr, &snap.ftn_vpd_len);
  }
  snap.vtk_bits_ptr = vbits;
  snap.vtk_bits_nbytes = vnbytes;
  snap.ftn_bits_ptr = fbits;
  snap.ftn_bits_nbytes = fnbytes;
  lockstep_on_step(&snap);
}

int main()
{
  oracle::verify_vtk_oracle_version();

  vtkNew<vtkHyperTreeGrid> vgrid;
  void* fgrid = ftn_htg_grid_create();

  // -------------------------------------------------------------------------
  // Step 1: SetBranchFactor(2) & SetDimensions(4, 1, 3) -> 2D XZ squeezed!
  //   Orientation = 1, Axis = {2, 0}, CellDims = {3, 1, 2} (6 root slots)
  // -------------------------------------------------------------------------
  vgrid->SetBranchFactor(2);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_branch_factor(fgrid, 2), "lockstep set_branch_factor");

  vgrid->SetDimensions(4, 1, 3);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_dimensions(fgrid, 4, 1, 3), "lockstep set_dimensions");

  oracle::set_rectilinear_coords_both(
    vgrid, fgrid, { 0.0, 1.0, 2.5, 4.0 }, { 0.0 }, { -1.0, 0.0, 2.0 });
  emit_tree_step_snapshot(1, "SetDimensions_2D_XZ", nullptr, nullptr);

  // -------------------------------------------------------------------------
  // Step 2..4: NonOrientedCursor SubdivideLeaf & ToChild on Tree 1 (quadtree, k=4)
  // -------------------------------------------------------------------------
  vtkNew<vtkHyperTreeGridNonOrientedCursor> vcur;
  vgrid->InitializeNonOrientedCursor(vcur, 1, true);
  int32_t cst = 0;
  void* fcur = ftn_htg_cursor_create(fgrid, 1, true, &cst);
  ORACLE_EQ(FTN_HTG_OK, cst, "lockstep cursor_create");

  vtkHyperTree* vt1 = vgrid->GetTree(1, false);
  void* ft1 = ftn_htg_cursor_get_tree(fcur);

  // Step 2: Subdivide root (vertex 0, depth 0) -> ParentToElderChild = [1]
  vcur->SubdivideLeaf();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), "lockstep subdivide 0");
  emit_tree_step_snapshot(2, "SubdivideLeaf_root", vt1, ft1);

  // Step 3: ToChild(2) -> SubdivideLeaf (vertex 3, depth 1)
  //   -> ParentToElderChild = [1, UINT_MAX, UINT_MAX, 5], VPD = [1, 4, 4]
  vcur->ToChild(2);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fcur, 2), "lockstep to_child(2)");
  vcur->SubdivideLeaf();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), "lockstep subdivide child 2");
  emit_tree_step_snapshot(3, "SubdivideLeaf_child2", vt1, ft1);

  // Step 4: ToChild(1) -> SubdivideLeaf (vertex 6, depth 2)
  //   -> ParentToElderChild = [1, UINT_MAX, UINT_MAX, 5, UINT_MAX, UINT_MAX, 9]
  vcur->ToChild(1);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_child(fcur, 1), "lockstep to_child(1)");
  vcur->SubdivideLeaf();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_subdivide_leaf(fcur), "lockstep subdivide grandchild");
  vcur->ToParent();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_parent(fcur), "lockstep to_parent");
  vcur->ToRoot();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_cursor_to_root(fcur), "lockstep to_root");
  emit_tree_step_snapshot(4, "SubdivideLeaf_grandchild", vt1, ft1);
  ftn_htg_cursor_destroy(fcur);

  // -------------------------------------------------------------------------
  // Step 5: BuildFromBreadthFirstOrderDescriptor on Tree 4 (`startIndex = 2`)
  //   Prefix bits [1, 0], then 5-bit 3-level quadtree descriptor [1, 0, 1, 0, 0]
  //   -> Levels=3, Vertices=9, Nodes=2, Leaves=7, ParentToElderChild=[1, MAX, 5, MAX, MAX]
  // -------------------------------------------------------------------------
  vtkNew<vtkBitArray> vdesc;
  int32_t bst = 0;
  void* fdesc = ftn_htg_bit_array_create(0, false, &bst);
  const bool raw_desc_bits[] = { true, false, /* start=2: */ true, false, true, false, false };
  for (bool b : raw_desc_bits)
  {
    vdesc->InsertNextValue(b ? 1 : 0);
    int64_t ins = -1;
    ftn_htg_bit_array_insert_next_value(fdesc, b, &ins);
  }
  vtkHyperTree* vt4 = vgrid->GetTree(4, true);
  int32_t tst4 = 0;
  void* ft4 = ftn_htg_grid_get_tree(fgrid, 4, true, &tst4);
  ORACLE_EQ(FTN_HTG_OK, tst4, "lockstep get_tree(4)");

  vt4->BuildFromBreadthFirstOrderDescriptor(vdesc, 5, 2);
  ORACLE_EQ(
    FTN_HTG_OK, ftn_htg_tree_build_from_bf_descriptor(ft4, fdesc, 5, 2), "lockstep build_bf");

  {
    const uint8_t* vb_ptr = reinterpret_cast<const uint8_t*>(vdesc->GetPointer(0));
    const uint8_t* fb_ptr = nullptr;
    int64_t fb_nbytes = 0;
    ftn_htg_bit_array_get_raw_bytes(fdesc, &fb_ptr, &fb_nbytes);
    emit_tree_step_snapshot(5, "BuildFromBreadthFirstOrderDescriptor", vt4, ft4, vb_ptr, 1, fb_ptr,
      fb_nbytes);
  }
  ftn_htg_bit_array_destroy(fdesc);

  // -------------------------------------------------------------------------
  // Step 6: Global Indexing + Mask MSB-first bits + Material Interface + GetPureMask()
  // -------------------------------------------------------------------------
  vgrid->InitializeLocalIndexNode();
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_initialize_local_index_node(fgrid), "lockstep init_gid");

  const int64_t ncells = static_cast<int64_t>(vgrid->GetGlobalNodeIndexMax() + 1);
  vtkNew<vtkBitArray> vmask;
  vmask->SetNumberOfTuples(static_cast<vtkIdType>(ncells));
  for (int64_t i = 0; i < ncells; ++i)
  {
    vmask->SetValue(static_cast<vtkIdType>(i), 0);
  }
  void* fmask = ftn_htg_bit_array_create(ncells, false, &bst);

  // Mask leaf vertex gid=2 (in Tree 1) and leaf vertex gid=19 (in Tree 4, across byte 2!)
  vmask->SetValue(2, 1);
  vmask->SetValue(19, 1);
  ftn_htg_bit_array_set_value(fmask, 2, true);
  ftn_htg_bit_array_set_value(fmask, 19, true);

  vgrid->SetMask(vmask);
  ORACLE_EQ(FTN_HTG_OK, ftn_htg_grid_set_mask(fgrid, fmask), "lockstep set_mask");
  ftn_htg_bit_array_destroy(fmask);

  // Attach Material Interface arrays marking leaf gid=10 (Tree 1 grandchild) as mixed (< 2.0)!
  const char* n_norm = "Normals";
  const char* n_intr = "Intercepts";
  vtkNew<vtkDoubleArray> vnorm, vintr;
  vnorm->SetName(n_norm);
  vnorm->SetNumberOfComponents(3);
  vnorm->SetNumberOfTuples(static_cast<vtkIdType>(ncells));
  vintr->SetName(n_intr);
  vintr->SetNumberOfComponents(3);
  vintr->SetNumberOfTuples(static_cast<vtkIdType>(ncells));

  ftn_htg_grid_cell_data_add_array(fgrid, n_norm, 7, FTN_HTG_TYPE_REAL64, 3, ncells);
  ftn_htg_grid_cell_data_add_array(fgrid, n_intr, 10, FTN_HTG_TYPE_REAL64, 3, ncells);
  for (int64_t g = 0; g < ncells; ++g)
  {
    const double nrm[3] = { 0.0, 0.0, 1.0 };
    const double itc[3] = { 0.5, 0.0, (g == 10) ? 1.0 : 2.0 };
    vnorm->SetTuple(static_cast<vtkIdType>(g), nrm);
    vintr->SetTuple(static_cast<vtkIdType>(g), itc);
    ftn_htg_grid_cell_data_set_tuple_r64(fgrid, n_norm, 7, g, nrm, 3);
    ftn_htg_grid_cell_data_set_tuple_r64(fgrid, n_intr, 10, g, itc, 3);
  }
  vgrid->GetCellData()->AddArray(vnorm);
  vgrid->GetCellData()->AddArray(vintr);
  vgrid->SetHasInterface(true);
  vgrid->SetInterfaceNormalsName(n_norm);
  vgrid->SetInterfaceInterceptsName(n_intr);
  vgrid->SetMask(vmask); // ensure PureMask cache clean in VTK

  ftn_htg_grid_set_has_interface(fgrid, true);
  ftn_htg_grid_set_interface_normals_name(fgrid, n_norm, 7);
  ftn_htg_grid_set_interface_intercepts_name(fgrid, n_intr, 10);

  // Trigger GetPureMask() in lockstep: VTK first, then Fortran!
  vtkBitArray* vpure = vgrid->GetPureMask();
  const uint8_t* fpure_bytes = nullptr;
  int64_t fpure_nbits = 0, fpure_nbytes = 0;
  ORACLE_EQ(FTN_HTG_OK,
    ftn_htg_grid_get_pure_mask_raw_bytes(fgrid, &fpure_bytes, &fpure_nbits, &fpure_nbytes),
    "lockstep get_pure_mask");

  emit_tree_step_snapshot(6, "GetPureMask_with_Interface", vt1, ft1,
    reinterpret_cast<const uint8_t*>(vpure->GetPointer(0)), (vpure->GetNumberOfValues() + 7) / 8,
    fpure_bytes, fpure_nbytes);

  // Full structural & byte-for-byte verification
  oracle::CompareOptions opts;
  opts.check_coordinates = true;
  opts.check_global_indexing = true;
  opts.check_masks = true;
  opts.check_pure_mask = true;
  opts.check_validate = true;
  oracle::assert_grids_equivalent(vgrid, fgrid, "lockstep_final", opts);

  ftn_htg_grid_destroy(fgrid);
  std::cout << "[PASS] test_lockstep_driver completed all 6 lockstep phases!" << std::endl;
  return 0;
}
