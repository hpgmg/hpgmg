 HPGMG: High-performance Geometric Multigrid
===========================================

HPGMG implements full multigrid (FMG) algorithms using finite-volume and
finite-element methods.  Different algorithmic variants adjust the
arithmetic intensity and architectural properties that are tested.
These FMG methods converge up to discretization error in one F-cycle,
thus may be considered direct solvers.  An F-cycle visits the finest
level a total of two times, the first coarsening (8x smaller) 4 times,
the second coarsening 6 times, etc.

#General installation

Run configure, and then make as instructed:

    $ ./configure --CC=/path/to/mpicc

This will create the finite volume solver.
The finite volume solver can be disabled by configuring with `--no-fv`.
The finite element solver can be enabled by configuring with `--fe`.

# HPGMG-FE: Finite Element FAS solver

The finite-element full approximation scheme (FAS) solver has higher
arithmetic intensity due to quadrature, metric terms, and higher
approximation order.  HPGMG-FE configuration is via run-time options
described below.

## Installation

The finite-element solver currently requires PETSc, obtained via

```
$ git clone https://bitbucket.org/petsc/petsc
$ cd petsc
$ export PETSC_DIR=`pwd` PETSC_ARCH=name-of-your-choice
$ ./configure --with-debugging=0  ...
$ make
```

After the above, you can build HPGMG-FE with `make`.  If you are using a
non-batch system, the test suite can be run with `make test`.

## Running

Typical runs sample across a range of problem sizes in understand the
range of problem sizes that can be solved efficiently.

* `-op_type <poisson1,poisson2,poisson2affine>` specifies the operator
  type and basis function degree (2=Q_2, biquadratic elements).

* `-local min,max` specifies the minimum and maximum number of elements
  per MPI rank.

* `-maxsamples <6>` is the maximum number of samples across the range of
  problem sizes.

* `-repeat <5>` number of samples for each problem size.  Some machines
  have significant performance variability, requiring several solves to
  extract meaningful timing information.

* `-mintime <1>` minimum interval (in seconds) for repeatedly solving
  each problem size. A negative or zero value will disable this constraint.
  If both `-repeat` and `-mintime` are specified, sampling will continue
  until both conditions are satisfied.

The run solves the smallest problem size first to provide instant
feedback about incompatible configuration.  Then it solves the largest
problem size to ensure that the whole run will not exceed machine memory
and to "warm up" the machine (e.g., faulting memory pages).  Finally,
timing data is collected across the range of problem sizes.  Timing is
reported in seconds per solve, manually-counted gigaflops (GF), and
millions of equations solved per second (MEq/s).

```
$ mpiexec -n 4 mpich-opt/bin/hpgmg-fe sample -local 50,10000 -repeat 3 -op_type poisson2
Finite Element FAS Performance Sampler on process grid [1 2 2] = 4
Max memory per MPI rank: 0.006230 GB
Small Test G[]
Q2 G[   4   6   8] P[  1  2  2]  3.825e-03 s    4.686462 GF    0.510524 MEq/s
Large Test
Q2  0 e_max 9.98e-01(0.0) e_L2 1.00e+00(0.0) r_2 0.00e+00(0.0) G[   1   1   1] L[  1  1  1] P[  1  1  1]
Q2  1 e_max 3.11e-01(1.7) e_L2 3.59e-01(1.5) r_2 1.15e-02(0.0) G[   2   2   2] L[  2  2  2] P[  1  1  1]
Q2  2 e_max 1.10e-01(1.5) e_L2 5.77e-02(2.6) r_2 6.25e-03(0.9) G[   4   4   4] L[  4  4  4] P[  1  1  1]
Q2  3 e_max 1.53e-02(2.8) e_L2 7.67e-03(2.9) r_2 3.75e-03(0.7) G[   8   8   8] L[  8  4  4] P[  1  2  2]
Q2  4 e_max 1.90e-03(3.0) e_L2 9.68e-04(3.0) r_2 1.47e-03(1.4) G[  16  16  16] L[ 16  8  8] P[  1  2  2]
Q2  5 e_max 2.39e-04(3.0) e_L2 1.21e-04(3.0) r_2 6.30e-04(1.2) G[  32  32  32] L[ 32 16 16] P[  1  2  2]
Q2 G[  32  32  32] P[  1  2  2]  2.298e-01 s   11.780769 GF    1.195290 MEq/s
Max memory per MPI rank: 0.020877 GB
Starting performance sampling
Q2 G[   4   6   8] P[  1  2  2]  3.628e-03 s    5.021407 GF    0.547011 MEq/s
Q2 G[   4   6   8] P[  1  2  2]  3.591e-03 s    5.071860 GF    0.552628 MEq/s
Q2 G[   4   6   8] P[  1  2  2]  3.603e-03 s    5.053782 GF    0.550658 MEq/s
Q2 G[   8   8  12] P[  1  2  2]  9.157e-03 s    5.596428 GF    0.786602 MEq/s
Q2 G[   8   8  12] P[  1  2  2]  7.466e-03 s    6.857400 GF    0.964108 MEq/s
Q2 G[   8   8  12] P[  1  2  2]  7.420e-03 s    6.900650 GF    0.970189 MEq/s
Q2 G[  12  12  16] P[  1  2  2]  1.638e-02 s    9.708319 GF    1.254641 MEq/s
Q2 G[  12  12  16] P[  1  2  2]  1.620e-02 s    9.825944 GF    1.270171 MEq/s
Q2 G[  12  12  16] P[  1  2  2]  2.060e-02 s    7.698593 GF    0.995174 MEq/s
Q2 G[  16  16  24] P[  1  2  2]  3.397e-02 s   11.023745 GF    1.499027 MEq/s
Q2 G[  16  16  24] P[  1  2  2]  3.521e-02 s   11.132441 GF    1.514220 MEq/s
Q2 G[  16  16  24] P[  1  2  2]  3.494e-02 s   11.175234 GF    1.520040 MEq/s
Q2 G[  16  24  32] P[  1  2  2]  6.892e-02 s   11.438116 GF    1.525116 MEq/s
Q2 G[  16  24  32] P[  1  2  2]  6.944e-02 s   11.348496 GF    1.513570 MEq/s
Q2 G[  16  24  32] P[  1  2  2]  6.867e-02 s   11.467359 GF    1.529423 MEq/s
Q2 G[  32  32  32] P[  1  2  2]  1.769e-01 s   11.647631 GF    1.541970 MEq/s
Q2 G[  32  32  32] P[  1  2  2]  1.585e-01 s   13.053139 GF    1.728495 MEq/s
Q2 G[  32  32  32] P[  1  2  2]  1.594e-01 s   12.884064 GF    1.706106 MEq/s
```

### Debugging and testing convergence
HPGMG-FE has several modes for debugging and testing convergence.

```
$ mpiexec -n 4 build/bin/hpgmg-fe fmg -mg_monitor -p 1,2,2 -M 16,24,32 -op_type poisson2
Q2  0 e_max 2.85e-01(0.0) e_L2 2.37e-01(0.0) r_2 1.46e-13(0.0) G[   2   3   4] L[  2  3  4] P[  1  1  1]
Q2  1 e_max 5.13e-02(2.5) e_L2 3.53e-02(2.7) r_2 1.82e-03(0.0) G[   4   6   8] L[  4  6  8] P[  1  1  1]
Q2  2 e_max 7.79e-03(2.7) e_L2 4.62e-03(2.9) r_2 9.77e-04(0.9) G[   8  12  16] L[  8  6  8] P[  1  2  2]
Q2  3 e_max 9.80e-04(3.0) e_L2 5.84e-04(3.0) r_2 4.37e-04(1.2) G[  16  24  32] L[ 16 12 16] P[  1  2  2]
F(2,3)  0: |e|_2/|u|_2 2.97e-05  |r|_2/|f|_2 4.37e-04
V(2,3)  1: |e|_2/|u|_2 1.93e-05  |r|_2/|f|_2 4.06e-06
V(2,3)  2: |e|_2/|u|_2 1.95e-05  |r|_2/|f|_2 1.88e-07
```

# HPGMG-FV: Finite Volume solver

The finite-volume solver uses cell-centered methods with constant or
variable coefficients.  This implementation requires OpenMP and cannot
be configured at run-time.  See `./configure --help` for configuration
options.  Be sure to pass suitable OpenMP flags, e.g.,

    $ ./configure --CC=/path/to/mpicc --CFLAGS=-fopenmp

## Running

Using the Cray XC-30 at NERSC, Edison, with 96 core in an interactive shell, HPGMG-FV
generates the following output:

```
$ export OMP_NUM_THREADS=8                                                                             
$ aprun -n 8 -d 12  -N  2  -S 1  -ss  -cc numa_node ./arch-xc30-opt64/bin/hpgmg-fv      7  8           
Requested MPI_THREAD_FUNNELED, got MPI_THREAD_FUNNELED                                                                              
8 MPI Tasks of 8 threads (OMP_NESTED=FALSE)                                                                                         

attempting to create a   512^3 level using a   4^3 grid of 128^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=8.000, max=8                   
  rebuilding operator for level...  h=1.953125e-03  eigenvalue_max<2.000000e+00

attempting to create a   256^3 level using a   4^3 grid of  64^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=8.000, max=8                   

attempting to create a   128^3 level using a   4^3 grid of  32^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=8.000, max=8                   

attempting to create a    64^3 level using a   4^3 grid of  16^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=8.000, max=8                   

attempting to create a    32^3 level using a   4^3 grid of   8^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=8.000, max=8                   

attempting to create a    16^3 level using a   4^3 grid of   4^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=8.000, max=8                   

attempting to create a     8^3 level using a   2^3 grid of   4^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=1.000, max=1                   

attempting to create a     4^3 level using a   1^3 grid of   4^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=0.125, max=1                   

attempting to create a     2^3 level using a   1^3 grid of   2^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=0.125, max=1                   

attempting to create a     1^3 level using a   1^3 grid of   1^3 boxes...
  OMP_NESTED=FALSE OMP_NUM_THREADS=8 ... 1 teams of 8 threads            
  calculating boxes per process... target=0.125, max=1                   
 Building MPI subcommunicator for level 1...done                        
  Building MPI subcommunicator for level 2...done                        
  Building MPI subcommunicator for level 3...done                        
  Building MPI subcommunicator for level 4...done                        
  Building MPI subcommunicator for level 5...done                        
  Building MPI subcommunicator for level 6...done                        
  Building MPI subcommunicator for level 7...done                        
  Building MPI subcommunicator for level 8...done                        
  Building MPI subcommunicator for level 9...done                        
  rebuilding operator for level...  h=3.906250e-03  eigenvalue_max<2.000000e+00
  rebuilding operator for level...  h=7.812500e-03  eigenvalue_max<2.000000e+00
  rebuilding operator for level...  h=1.562500e-02  eigenvalue_max<2.000000e+00
  rebuilding operator for level...  h=3.125000e-02  eigenvalue_max<2.000000e+00
  rebuilding operator for level...  h=6.250000e-02  eigenvalue_max<2.000000e+00
  rebuilding operator for level...  h=1.250000e-01  eigenvalue_max<2.000000e+00
  rebuilding operator for level...  h=2.500000e-01  eigenvalue_max<2.000000e+00
  rebuilding operator for level...  h=5.000000e-01  eigenvalue_max<1.286089e+00
  rebuilding operator for level...  h=1.000000e+00  eigenvalue_max<1.000000e+00
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...                                                                    
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)                
done                                                                           
FMGSolve...
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)
done
FMGSolve...
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)
done
FMGSolve...
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)
done
FMGSolve...
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)
done
FMGSolve...
f-cycle,    norm=0.00000000000457238665 (4.572386651768147e-12)
done
h = 1.953125e-03, error = 3.696307577501308e-09
                                     0            1            2            3            4            5            6            7            8            9
box dimension                    128^3         64^3         32^3         16^3          8^3          4^3          4^3          4^3          2^3          1^3        total
------------------        ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------
smooth                        0.233412     0.064273     0.005881     0.001720     0.001034     0.000923     0.000162     0.000178     0.000162     0.000000     0.307743
residual                      0.052603     0.006971     0.000547     0.000169     0.000111     0.000105     0.000019     0.000022     0.000019     0.000021     0.060587
applyOp                       0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000023     0.000023
BLAS1                         0.030815     0.000512     0.000169     0.000084     0.000065     0.000075     0.000012     0.000015     0.000016     0.000275     0.032039
BLAS3                         0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000
Boundary Conditions           0.006946     0.004180     0.001653     0.000625     0.000238     0.000111     0.000036     0.000068     0.000048     0.000034     0.013938
Restriction                   0.008624     0.001549     0.000156     0.000108     0.000109     0.000117     0.000537     0.000027     0.000029     0.000000     0.011256
  local restriction           0.008594     0.001539     0.000151     0.000103     0.000104     0.000110     0.000018     0.000020     0.000022     0.000000     0.010662
  pack MPI buffers            0.000016     0.000005     0.000001     0.000002     0.000002     0.000002     0.000002     0.000002     0.000002     0.000000     0.000034
  unpack MPI buffers          0.000011     0.000001     0.000001     0.000001     0.000001     0.000001     0.000113     0.000001     0.000001     0.000000     0.000132
  MPI_Isend                   0.000001     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000003
  MPI_Irecv                   0.000001     0.000001     0.000001     0.000001     0.000000     0.000001     0.000013     0.000001     0.000001     0.000000     0.000019
  MPI_Waitall                 0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000389     0.000000     0.000000     0.000000     0.000390
Interpolation                 0.016691     0.002704     0.000414     0.000158     0.000124     0.000124     0.000247     0.000027     0.000033     0.000000     0.020523
  local interpolation         0.016674     0.002696     0.000409     0.000154     0.000120     0.000120     0.000020     0.000021     0.000025     0.000000     0.020238
  pack MPI buffers            0.000005     0.000002     0.000002     0.000001     0.000001     0.000001     0.000121     0.000001     0.000002     0.000000     0.000137
  unpack MPI buffers          0.000010     0.000003     0.000001     0.000001     0.000001     0.000001     0.000002     0.000002     0.000002     0.000000     0.000023
  MPI_Isend                   0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000056     0.000000     0.000000     0.000000     0.000057
  MPI_Irecv                   0.000000     0.000000     0.000000     0.000000     0.000001     0.000000     0.000000     0.000001     0.000001     0.000000     0.000004
  MPI_Waitall                 0.000000     0.000001     0.000000     0.000000     0.000000     0.000000     0.000045     0.000000     0.000000     0.000000     0.000048
Ghost Zone Exchange           0.024454     0.009990     0.003605     0.001738     0.001092     0.000874     0.000682     0.000200     0.000222     0.000057     0.042915
  local exchange              0.004445     0.001670     0.000481     0.000287     0.000172     0.000138     0.000016     0.000146     0.000159     0.000039     0.007552
  pack MPI buffers            0.006175     0.002598     0.000722     0.000252     0.000150     0.000125     0.000122     0.000011     0.000013     0.000004     0.010172
  unpack MPI buffers          0.002692     0.001118     0.000496     0.000276     0.000166     0.000129     0.000130     0.000013     0.000018     0.000005     0.005044
  MPI_Isend                   0.000259     0.000360     0.000176     0.000185     0.000294     0.000176     0.000146     0.000003     0.000003     0.000001     0.001603
  MPI_Irecv                   0.000072     0.000070     0.000028     0.000027     0.000028     0.000025     0.000035     0.000002     0.000004     0.000001     0.000293
  MPI_Waitall                 0.010801     0.004161     0.001691     0.000698     0.000265     0.000264     0.000213     0.000004     0.000003     0.000001     0.018100
MPI_collectives               0.000657     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000037     0.000694
------------------        ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------
Total by level                0.368311     0.089306     0.012201     0.004555     0.002779     0.002270     0.001697     0.000551     0.000548     0.000436     0.482656

   Total time in MGBuild      0.171733 seconds
   Total time in MGSolve      0.482693 seconds
      number of v-cycles             1
Bottom solver iterations            10

            Performance      2.781e+08 DOF/s
```

