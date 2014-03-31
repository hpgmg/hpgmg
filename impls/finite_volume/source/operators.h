//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
  void                  apply_op(level_type * level, int Ax_id,  int x_id, double a, double b);
  void                  residual(level_type * level, int res_id, int x_id, int rhs_id, double a, double b);
  void                    smooth(level_type * level, int phi_id, int rhs_id, double a, double b);
  void          rebuild_operator(level_type * level, level_type *fromLevel, double a, double b);
//------------------------------------------------------------------------------------------------------------------------------
  void               restriction(level_type * level_c, int id_c, level_type *level_f, int id_f, int restrictionType);
  void      interpolation_vcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c);
  void      interpolation_fcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c);
//------------------------------------------------------------------------------------------------------------------------------
  void         exchange_boundary(level_type * level, int id_a, int justFaces);
  void        apply_BCs_linear(  level_type * level, int x_id);
  void        apply_BCs_4thOrder(level_type * level, int x_id);
  void        apply_BCs_2ndOrder(level_type * level, int x_id);
//------------------------------------------------------------------------------------------------------------------------------
double                       dot(level_type * level, int id_a, int id_b);
double                      norm(level_type * level, int id_a);
double                      mean(level_type * level, int id_a);
double                     error(level_type * level, int id_a, int id_b);
  void                 add_grids(level_type * level, int id_c, double scale_a, int id_a, double scale_b, int id_b);
  void                scale_grid(level_type * level, int id_c, double scale_a, int id_a);
  void                 zero_grid(level_type * level, int id_a);
  void                shift_grid(level_type * level, int id_c, int id_a, double shift_a);
  void                 mul_grids(level_type * level, int id_c, double scale, int id_a, int id_b);
  void               invert_grid(level_type * level, int id_c, double scale_a, int id_a);
  void initialize_grid_to_scalar(level_type * level, int id_a, double scalar);
//------------------------------------------------------------------------------------------------------------------------------
  void      project_cell_to_face(level_type * level, int id_cell, int id_face, int dir);
//------------------------------------------------------------------------------------------------------------------------------
  void              matmul_grids(level_type * level, double *C, int *id_A, int *id_B, int rows, int cols, int A_equals_B_transpose);
//------------------------------------------------------------------------------------------------------------------------------
  void        initialize_problem(level_type * level, double hLevel, double a, double b);
  void   initialize_valid_region(level_type * level);
//------------------------------------------------------------------------------------------------------------------------------
