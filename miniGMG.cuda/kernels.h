//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void ConfigureGPU();
//------------------------------------------------------------------------------------------------------------------------------
void                  apply_op(domain_type *domain, int level, int     Ax_id, int      x_id, double a, double b, double hLevel);
void                    smooth(domain_type *domain, int level, int    phi_id, int    rhs_id, double a, double b, double hLevel, int s);
void                  residual(domain_type *domain, int level, int    res_id, int    phi_id, int rhs_id, double a, double b, double hLevel);
void            rebuild_lambda(domain_type *domain, int level, double a, double b, double hLevel);
//------------------------------------------------------------------------------------------------------------------------------
void                 zero_grid(domain_type *domain, int level, int grid_id);
void initialize_grid_to_scalar(domain_type *domain, int level, int grid_id, double h, double scalar);
void                scale_grid(domain_type *domain, int level, int id_c, double scale_a, int id_a);
void                 add_grids(domain_type *domain, int level, int id_c, double scale_a, int id_a, double scale_b, int id_b);
void                shift_grid(domain_type *domain, int level, int id_c, int id_a, double scale_a);
void                 mul_grids(domain_type *domain, int level, int id_c, double scale, int id_a, int id_b);
void               restriction(domain_type *domain, int level, int coarse_id, int   fine_id);
void         restriction_betas(domain_type *domain, int level);
void             interpolation(domain_type *domain, int level, int   fine_id, int coarse_id);
void               norm_on_gpu(domain_type *domain, int level, int grid_id, double *gpu_norm);
void                dot_on_gpu(domain_type *domain, int level, int id_a, int id_b, double *gpu_dot);
void               mean_on_gpu(domain_type *domain, int level, int id_a, double *gpu_mean);
void      project_cell_to_face(domain_type *domain, int level, int id_cell, int id_face, int dir);
//------------------------------------------------------------------------------------------------------------------------------
void   initialize_exact_on_gpu(domain_type * domain, int level, double hLevel, double a, double b);
//------------------------------------------------------------------------------------------------------------------------------
void surface_buffers_to_send_buffer(domain_type *domain, int level, int grid_id);
void   recv_buffer_to_ghost_buffers(domain_type *domain, int level, int grid_id);
//------------------------------------------------------------------------------------------------------------------------------
void grid_to_surface_buffers(domain_type *domain, int level, int grid_id);
void   ghost_buffers_to_grid(domain_type *domain, int level, int grid_id);
void surface_buffers_to_ghost_buffers(domain_type *domain, int level, int grid_id);
//------------------------------------------------------------------------------------------------------------------------------



