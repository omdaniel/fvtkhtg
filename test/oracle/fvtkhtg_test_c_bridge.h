#ifndef FVTKHTG_TEST_C_BRIDGE_H
#define FVTKHTG_TEST_C_BRIDGE_H

/*
 * Test-only C/C++ header for the Fortran `iso_c_binding` shim in
 * `test/oracle/fvtkhtg_test_c_bridge.f90`.
 *
 * NOTE: `libfvtkhtg` in `src/` is 100% pure Fortran 2018 with ZERO C/C++
 * bindings. This bridge exists exclusively inside `test/oracle/` so a single
 * Debug executable can drive C++ VTK 9.7+ (`vtkHyperTreeGrid`) and pure
 * Fortran (`htg_grid_t`) in lockstep and inspect both via zero-copy views
 * and `/usr/bin/lldb`.
 */

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Status & data-type constants matching `fvtkhtg_kinds_m` */
enum {
  FTN_HTG_OK = 0,
  FTN_HTG_ERR_INVALID_ARG = 1,
  FTN_HTG_ERR_OUT_OF_BOUNDS = 2,
  FTN_HTG_ERR_INVALID_STATE = 3,
  FTN_HTG_ERR_TOPOLOGY = 4,
  FTN_HTG_ERR_INDEXING = 5,
  FTN_HTG_ERR_MASK = 6,
  FTN_HTG_ERR_INTERFACE = 7,
  FTN_HTG_ERR_NOT_FOUND = 8
};

enum {
  FTN_HTG_INDEX_UNSET = 0,
  FTN_HTG_INDEX_IMPLICIT = 1,
  FTN_HTG_INDEX_EXPLICIT = 2
};

enum {
  FTN_HTG_TYPE_REAL64 = 1,
  FTN_HTG_TYPE_REAL32 = 2,
  FTN_HTG_TYPE_INT64 = 3,
  FTN_HTG_TYPE_INT32 = 4,
  FTN_HTG_TYPE_BIT = 5
};

/* -------------------------------------------------------------------------
 * BitArray (`htg_bit_array_t`) handle operations & zero-copy raw byte view
 * ------------------------------------------------------------------------- */
void* ftn_htg_bit_array_create(int64_t num_bits, bool init_val, int32_t* status_out);
void ftn_htg_bit_array_destroy(void* bit_arr);
int32_t ftn_htg_bit_array_resize(void* bit_arr, int64_t num_bits);
int64_t ftn_htg_bit_array_get_number_of_values(void* bit_arr);
int64_t ftn_htg_bit_array_get_number_of_bytes(void* bit_arr);
int32_t ftn_htg_bit_array_set_value(void* bit_arr, int64_t id, bool val);
int32_t ftn_htg_bit_array_insert_value(void* bit_arr, int64_t id, bool val);
int32_t ftn_htg_bit_array_insert_next_value(void* bit_arr, bool val, int64_t* inserted_id_out);
bool ftn_htg_bit_array_get_value(void* bit_arr, int64_t id, int32_t* status_out);
void ftn_htg_bit_array_fill(void* bit_arr, bool val);
void ftn_htg_bit_array_get_raw_bytes(
  void* bit_arr, const uint8_t** ptr_out, int64_t* nbytes_out);

/* -------------------------------------------------------------------------
 * Grid (`htg_grid_t`) lifecycle, dimensions, coordinates & root indexing
 * ------------------------------------------------------------------------- */
void* ftn_htg_grid_create(void);
void ftn_htg_grid_destroy(void* grid);
void ftn_htg_grid_initialize(void* grid);

int32_t ftn_htg_grid_set_dimensions(void* grid, int32_t ni, int32_t nj, int32_t nk);
void ftn_htg_grid_get_dimensions(void* grid, int32_t dims_out[3]);
void ftn_htg_grid_get_cell_dims(void* grid, int32_t cell_dims_out[3]);
int32_t ftn_htg_grid_get_dimension(void* grid);
int32_t ftn_htg_grid_get_orientation(void* grid);
void ftn_htg_grid_get_axes(void* grid, int32_t axes_out[3], int32_t* naxes_out);

int32_t ftn_htg_grid_set_branch_factor(void* grid, int32_t factor);
int32_t ftn_htg_grid_get_branch_factor(void* grid);
int32_t ftn_htg_grid_get_number_of_children(void* grid);
int32_t ftn_htg_grid_set_depth_limiter(void* grid, int32_t limit);
int32_t ftn_htg_grid_get_depth_limiter(void* grid);

void ftn_htg_grid_set_transposed_root_indexing(void* grid, bool flag);
bool ftn_htg_grid_get_transposed_root_indexing(void* grid);

int32_t ftn_htg_grid_set_x_coordinates(void* grid, const double* coords, int64_t n);
int32_t ftn_htg_grid_set_y_coordinates(void* grid, const double* coords, int64_t n);
int32_t ftn_htg_grid_set_z_coordinates(void* grid, const double* coords, int64_t n);
void ftn_htg_grid_get_x_coordinates_view(void* grid, const double** ptr_out, int64_t* len_out);
void ftn_htg_grid_get_y_coordinates_view(void* grid, const double** ptr_out, int64_t* len_out);
void ftn_htg_grid_get_z_coordinates_view(void* grid, const double** ptr_out, int64_t* len_out);

int32_t ftn_htg_grid_get_index_from_level_zero_coordinates(
  void* grid, int32_t i, int32_t j, int32_t k, int64_t* tree_id_out);
int32_t ftn_htg_grid_get_level_zero_coordinates_from_index(
  void* grid, int64_t tree_id, int32_t* i_out, int32_t* j_out, int32_t* k_out);

/* -------------------------------------------------------------------------
 * Sparse Forest & Aggregations
 * ------------------------------------------------------------------------- */
int64_t ftn_htg_grid_get_max_number_of_trees(void* grid);
int64_t ftn_htg_grid_get_number_of_non_empty_trees(void* grid);
int64_t ftn_htg_grid_get_number_of_vertices(void* grid);
int64_t ftn_htg_grid_get_number_of_leaves(void* grid);
int64_t ftn_htg_grid_get_number_of_levels(void* grid);
int64_t ftn_htg_grid_get_global_node_index_max(void* grid);
int32_t ftn_htg_grid_initialize_local_index_node(void* grid);

bool ftn_htg_grid_has_tree(void* grid, int64_t tree_id);
void* ftn_htg_grid_get_tree(void* grid, int64_t tree_id, bool create, int32_t* status_out);
int32_t ftn_htg_grid_remove_tree(void* grid, int64_t tree_id);
int32_t ftn_htg_grid_get_tree_by_pos(
  void* grid, int64_t pos, int64_t* tree_id_out, void** tree_ptr_out);

/* -------------------------------------------------------------------------
 * Tree (`htg_tree_t`) topology, indexing, BF descriptors & zero-copy views
 * ------------------------------------------------------------------------- */
int64_t ftn_htg_tree_get_tree_index(void* tree);
int32_t ftn_htg_tree_get_branch_factor(void* tree);
int32_t ftn_htg_tree_get_dimension(void* tree);
int32_t ftn_htg_tree_get_number_of_children(void* tree);
int64_t ftn_htg_tree_get_number_of_levels(void* tree);
int64_t ftn_htg_tree_get_number_of_vertices(void* tree);
int64_t ftn_htg_tree_get_number_of_nodes(void* tree);
int64_t ftn_htg_tree_get_number_of_leaves(void* tree);
int64_t ftn_htg_tree_get_number_of_vertices_at_depth(
  void* tree, int32_t depth, int32_t* status_out);

bool ftn_htg_tree_is_leaf(void* tree, int64_t local_id, int32_t* status_out);
bool ftn_htg_tree_is_terminal_node(void* tree, int64_t local_id, int32_t* status_out);
int64_t ftn_htg_tree_get_elder_child_index(void* tree, int64_t local_id, int32_t* status_out);
int32_t ftn_htg_tree_subdivide_leaf(void* tree, int64_t local_id, int32_t level);

int32_t ftn_htg_tree_set_global_index_start(void* tree, int64_t start_id);
int64_t ftn_htg_tree_get_global_index_start(void* tree);
int32_t ftn_htg_tree_set_global_index_from_local(void* tree, int64_t local_id, int64_t global_id);
int64_t ftn_htg_tree_get_global_index_from_local(
  void* tree, int64_t local_id, int32_t* status_out);
int32_t ftn_htg_tree_get_indexing_mode(void* tree);
int64_t ftn_htg_tree_get_max_global_index(void* tree);

void ftn_htg_tree_get_parent_to_elder_child_view(
  void* tree, const int64_t** ptr_out, int64_t* len_out);
void ftn_htg_tree_get_vertices_per_depth_view(
  void* tree, const int64_t** ptr_out, int64_t* len_out);
void ftn_htg_tree_get_global_index_table_view(
  void* tree, const int64_t** ptr_out, int64_t* len_out);

int32_t ftn_htg_tree_build_from_bf_descriptor(
  void* tree, void* descriptor_bit_arr, int64_t num_bits, int64_t start_index);
void* ftn_htg_i64_vec_create(void);
void ftn_htg_i64_vec_destroy(void* vec);
void ftn_htg_i64_vec_get_view(void* vec, const int64_t** ptr_out, int64_t* len_out);
int32_t ftn_htg_tree_compute_bf_descriptor_append(void* tree, void* inout_descriptor_bit_arr,
  void* inout_vpd_vec, void* inout_bfid_vec, void* opt_input_mask_bit_arr, bool has_depth_limiter,
  int32_t depth_limiter);

/* -------------------------------------------------------------------------
 * Build-time NonOrientedCursor (`htg_cursor_t`)
 * ------------------------------------------------------------------------- */
void* ftn_htg_cursor_create(void* grid, int64_t tree_id, bool create, int32_t* status_out);
void ftn_htg_cursor_destroy(void* cursor);
int32_t ftn_htg_cursor_initialize(void* cursor, void* grid, int64_t tree_id, bool create);

bool ftn_htg_cursor_has_tree(void* cursor);
void* ftn_htg_cursor_get_tree(void* cursor);
int64_t ftn_htg_cursor_get_vertex_id(void* cursor);
int64_t ftn_htg_cursor_get_global_node_index(void* cursor, int32_t* status_out);
int32_t ftn_htg_cursor_get_dimension(void* cursor);
int32_t ftn_htg_cursor_get_number_of_children(void* cursor);
int32_t ftn_htg_cursor_get_level(void* cursor);

bool ftn_htg_cursor_is_leaf(void* cursor, int32_t* status_out);
bool ftn_htg_cursor_is_terminal_node(void* cursor, int32_t* status_out);
bool ftn_htg_cursor_is_root(void* cursor);

int32_t ftn_htg_cursor_to_root(void* cursor);
int32_t ftn_htg_cursor_to_parent(void* cursor);
int32_t ftn_htg_cursor_to_child(void* cursor, int32_t ichild);
int32_t ftn_htg_cursor_subdivide_leaf(void* cursor);

int32_t ftn_htg_cursor_set_global_index_start(void* cursor, int64_t offset);
int32_t ftn_htg_cursor_set_global_index_from_local(void* cursor, int64_t global_id);
int32_t ftn_htg_cursor_set_mask(void* cursor, bool state);
bool ftn_htg_cursor_is_masked(void* cursor, int32_t* status_out);

/* -------------------------------------------------------------------------
 * Grid Mask, PureMask, CellData, Material Interfaces & Structural Validate
 * ------------------------------------------------------------------------- */
int32_t ftn_htg_grid_set_mask(void* grid, void* bit_arr);
int32_t ftn_htg_grid_clear_mask(void* grid);
bool ftn_htg_grid_has_mask(void* grid);
void ftn_htg_grid_get_mask_raw_bytes(
  void* grid, const uint8_t** ptr_out, int64_t* nbits_out, int64_t* nbytes_out);
int32_t ftn_htg_grid_set_cell_mask(void* grid, int64_t global_id, bool state);
bool ftn_htg_grid_is_cell_masked(void* grid, int64_t global_id, int32_t* status_out);

int32_t ftn_htg_grid_get_pure_mask_raw_bytes(
  void* grid, const uint8_t** ptr_out, int64_t* nbits_out, int64_t* nbytes_out);
bool ftn_htg_grid_get_pure_mask_bit(void* grid, int64_t global_id, int32_t* status_out);

void ftn_htg_grid_set_has_interface(void* grid, bool flag);
bool ftn_htg_grid_get_has_interface(void* grid);
void ftn_htg_grid_set_interface_normals_name(void* grid, const char* name, int32_t name_len);
void ftn_htg_grid_set_interface_intercepts_name(void* grid, const char* name, int32_t name_len);

int32_t ftn_htg_grid_cell_data_add_array(void* grid, const char* name, int32_t name_len,
  int32_t data_type, int32_t num_components, int64_t num_tuples);
int32_t ftn_htg_grid_cell_data_ensure_tuples(void* grid, int64_t min_tuples);

int32_t ftn_htg_grid_cell_data_set_tuple_r64(void* grid, const char* name, int32_t name_len,
  int64_t tuple_id, const double* vals, int32_t ncomp);
int32_t ftn_htg_grid_cell_data_get_tuple_r64(void* grid, const char* name, int32_t name_len,
  int64_t tuple_id, double* vals_out, int32_t ncomp);
int32_t ftn_htg_grid_cell_data_set_tuple_r32(void* grid, const char* name, int32_t name_len,
  int64_t tuple_id, const float* vals, int32_t ncomp);
int32_t ftn_htg_grid_cell_data_get_tuple_r32(void* grid, const char* name, int32_t name_len,
  int64_t tuple_id, float* vals_out, int32_t ncomp);
int32_t ftn_htg_grid_cell_data_set_tuple_i64(void* grid, const char* name, int32_t name_len,
  int64_t tuple_id, const int64_t* vals, int32_t ncomp);
int32_t ftn_htg_grid_cell_data_get_tuple_i64(void* grid, const char* name, int32_t name_len,
  int64_t tuple_id, int64_t* vals_out, int32_t ncomp);
int32_t ftn_htg_grid_cell_data_set_tuple_i32(void* grid, const char* name, int32_t name_len,
  int64_t tuple_id, const int32_t* vals, int32_t ncomp);
int32_t ftn_htg_grid_cell_data_get_tuple_i32(void* grid, const char* name, int32_t name_len,
  int64_t tuple_id, int32_t* vals_out, int32_t ncomp);
int32_t ftn_htg_grid_cell_data_set_component_bit(
  void* grid, const char* name, int32_t name_len, int64_t tuple_id, int32_t comp_id, bool val);
bool ftn_htg_grid_cell_data_get_component_bit(void* grid, const char* name, int32_t name_len,
  int64_t tuple_id, int32_t comp_id, int32_t* status_out);

int32_t ftn_htg_grid_cell_data_get_view_r64(void* grid, const char* name, int32_t name_len,
  const double** ptr_out, int32_t* ncomp_out, int64_t* ntuples_out);

int32_t ftn_htg_grid_validate(void* grid, char* err_msg_out, int32_t err_msg_cap);

/* -------------------------------------------------------------------------
 * Cell geometry, face neighbours & coarsening
 * ------------------------------------------------------------------------- */
int32_t ftn_htg_grid_get_level_zero_origin_and_size(
  void* grid, int64_t tree_id, double origin_out[3], double size_out[3]);
int32_t ftn_htg_cursor_get_origin(void* cursor, double origin_out[3]);
int32_t ftn_htg_cursor_get_size(void* cursor, double size_out[3]);
int32_t ftn_htg_cursor_get_bounds(void* cursor, double bounds_out[6]);
int32_t ftn_htg_cursor_get_point(void* cursor, double point_out[3]);
int32_t ftn_htg_cursor_get_child_index(void* cursor);

/* Face neighbour along physical `axis` (0 = x, 1 = y, 2 = z) in `direction`
 * (-1 or +1). Without a neighbour *has_tree_out is false and *tree_id_out is
 * -1 at the edge of the grid, or the index of the root cell that has no tree. */
int32_t ftn_htg_cursor_get_face_neighbor(void* cursor, int32_t axis, int32_t direction,
  bool* has_tree_out, int64_t* tree_id_out, int64_t* vertex_id_out, int64_t* global_id_out,
  int32_t* level_out, bool* is_leaf_out, bool* is_masked_out, double bounds_out[6]);

int32_t ftn_htg_cursor_coarsen_node(void* cursor);

#ifdef __cplusplus
}
#endif

#endif /* FVTKHTG_TEST_C_BRIDGE_H */
