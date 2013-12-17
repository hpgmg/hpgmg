//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include "cuda.h"
#include "cuda_runtime.h"
#define max_cudaEvents      100000
#define cudaEvent_smooth         0
#define cudaEvent_apply_op       1
#define cudaEvent_residual       2
#define cudaEvent_restriction    3
#define cudaEvent_interpolation  4
#define cudaEvent_blas1          5
#define cudaEvent_blas3          6
#define cudaEvent_PCIe           7
#define cudaEvent_pack           8
#define cudaEvent_bufcopy        9
#define cudaEvent_unpack        10
#define cudaEvent_s2buf         11
#define cudaEvent_buf2g         12
#define cudaEvent_communication 13
//------------------------------------------------------------------------------------------------------------------------------
uint64_t CycleTime();
//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  int level;
  int type;
  cudaEvent_t event;
} cudaMGEvent_type;

//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  int rank;                                                     // MPI rank of remote process
  int local_index;                                              // index in subdomains[] on remote process
  #ifdef __MPI
  struct{int buf;struct{int faces,edges,corners;}offset;}send;  // i.e. calculate offset as faceSize*faces + edgeSize*edges + cornerSize*corners
  struct{int buf;struct{int faces,edges,corners;}offset;}recv;  // i.e. calculate offset as faceSize*faces + edgeSize*edges + cornerSize*corners
  #endif
} neighbor_type;

//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  struct {int i, j, k;}low;                     // global coordinates of the first (non-ghost) element of subdomain at the finest resolution
  struct {int i, j, k;}dim;                     // subdomain dimensions at finest resolution
  int numLevels;                                // number of levels in MG v-cycle.  1=no restrictions
  int ghosts;                                   // ghost zone depth
  neighbor_type neighbors[27];                  // MPI rank and local index (on remote process) of each subdomain neighboring this subdomain
  box_type * levels;                            // pointer to an array of all coarsenings of this box
} subdomain_type;
//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  // timing information...
  struct {
    uint64_t build;   // total time spent building the coefficients...
    uint64_t MGSolve; // total time spent in MGSolve
    uint64_t vcycles; // total time spent in all vcycles (all CycleMG)
    uint64_t          send[10];
    uint64_t          recv[10];
    uint64_t          wait[10];
    uint64_t   collectives[10];
    float    communication[10];
    float            Total[10]; // in seconds ?
    float           smooth[10];
    float         apply_op[10];
    float         residual[10];
    float      restriction[10];
    float    interpolation[10];
    float            blas1[10];
    float            blas3[10];
    float            s2buf[10];
    float             pack[10];
    float          bufcopy[10];
    float           unpack[10];
    float            buf2g[10];
    float             PCIe[10];
    cudaEvent_t communicationStart;
    cudaEvent_t communicationEnd;
  }cycles;
  int vcycles_performed, BiCGStab_iterations;


  int                        rank_of_neighbor[27]; // = MPI rank of the neighbors of this process's subdomains (presumes rectahedral packing)
  #ifdef __MPI
  double *  __restrict__          send_buffer[27]; // = MPI send buffers (one per neighbor)    (buffers on the host)
  double *  __restrict__          recv_buffer[27]; // = MPI recieve buffer (one per neighbor)  (buffers on the host)
  double *                    gpu_send_buffer[27]; // =   host array of pointers to gpu send buffers
  double *                    gpu_recv_buffer[27]; // =   host array of pointers to gpu recv buffers
  double **   gpu_pointers_to_gpu_send_buffer;     // = device array of pointers to gpu send buffers
  double **   gpu_pointers_to_gpu_recv_buffer;     // = device array of pointers to gpu send buffers
  struct{int faces,edges,corners;}buffer_size[27]; // = MPI buffer size (one per neighbor) in the units of faces/edges/corners
  #endif

// n.b. i=unit stride
  struct {int i, j, k;}dim;                     // global dimensions at finest resolution
  struct {int i, j, k;}subdomains_in;           // total number of subdomains in i,j,k
  struct {int i, j, k;}subdomains_per_rank_in;  // number of subdomains in i,j,k
  int rank;                                     // MPI rank of this process
  int numsubdomains;                            // number of subdomains owned by this process
  int numLevels;                                // number of levels in MG v-cycle.  1=no restrictions
  int numGrids;                                 // number of grids (variables)
  int ghosts;                                   // ghost zone depth
//double h[MaxLevels];                          // h at each level

  subdomain_type *     subdomains;              // pointer to a list of all subdomains owned by this process
  subdomain_type * gpu_subdomains;              // GPU shadow copy of subdomains

  double  * gpu_reduction_buffer;		// CPU pointer to buffer on the GPU for reductions

 int num_cudaEvents;
 cudaMGEvent_type * cudaEvents;
 //struct{int level,type; cudaEvent_t event;}cudaEvents[max_cudaEvents];
} domain_type;


//------------------------------------------------------------------------------------------------------------------------------
 int create_subdomain(subdomain_type * box, 
                      int subdomain_low_i, int subdomain_low_j, int subdomain_low_k,
                      int subdomain_dim_i, int subdomain_dim_j, int subdomain_dim_k,
                      int numGrids, int ghosts, int numLevels);
void destroy_domain(domain_type * domain);
 int create_domain(domain_type * domain,
                   int subdomain_dim_i, int subdomain_dim_j, int subdomain_dim_k,
                   int subdomains_per_rank_in_i, int subdomains_per_rank_in_j, int subdomains_per_rank_in_k,
                   int ranks_in_i, int ranks_in_j, int ranks_in_k,
                   int rank, int numGrids, int ghosts, int numLevels);
void MGSolve(domain_type * domain, int u_id, int F_id, double a, double b, double h, double desired_mg_norm);
void exchange_boundary(domain_type *domain, int level, int grid_id, int exchange_faces, int exchange_edges, int exchange_vertices);
double norm(domain_type * domain, int level, int grid_id);
double mean(domain_type * domain, int level, int id_a);
double  dot(domain_type * domain, int level, int id_a, int id_b);
void print_timing(domain_type *domain);


//------------------------------------------------------------------------------------------------------------------------------
