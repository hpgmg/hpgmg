HPGMG: High-performance Geometric Multigrid
===========================================

HPGMG implements full multigrid (FMG) algorithms using finite-volume and
finite-element methods.  Different algorithmic variants adjust the
arithmetic intensity and architectural properties that are tested.
These FMG methods converge up to discretization error in one F-cycle,
thus may be considered direct solvers.  An F-cycle visits the finest
level a total of two times, the first coarsening (8x smaller) 4 times,
the second coarsening 6 times, etc.

# HPGMG-FV: Finite Volume solver

The finite-volume solver uses cell-centered methods with constant or
variable coefficients.  This implementation requires OpenMP and cannot
be configured at run-time.  See `./configure --help` for configuration
options.  Be sure to pass suitable OpenMP flags, e.g.,

    $ ./configure CC=/path/to/mpicc --CFLAGS=-fopenmp

The finite volume solver can be disabled by configuring with `--no-fv`.

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
Q2  0 e_∞ 9.98e-01(0.0) e_L2 1.00e+00(0.0) r_2 0.00e+00(0.0) G[   1   1   1] L[  1  1  1] P[  1  1  1]
Q2  1 e_∞ 3.11e-01(1.7) e_L2 3.59e-01(1.5) r_2 1.15e-02(0.0) G[   2   2   2] L[  2  2  2] P[  1  1  1]
Q2  2 e_∞ 1.10e-01(1.5) e_L2 5.77e-02(2.6) r_2 6.25e-03(0.9) G[   4   4   4] L[  4  4  4] P[  1  1  1]
Q2  3 e_∞ 1.53e-02(2.8) e_L2 7.67e-03(2.9) r_2 3.75e-03(0.7) G[   8   8   8] L[  8  4  4] P[  1  2  2]
Q2  4 e_∞ 1.90e-03(3.0) e_L2 9.68e-04(3.0) r_2 1.47e-03(1.4) G[  16  16  16] L[ 16  8  8] P[  1  2  2]
Q2  5 e_∞ 2.39e-04(3.0) e_L2 1.21e-04(3.0) r_2 6.30e-04(1.2) G[  32  32  32] L[ 32 16 16] P[  1  2  2]
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
Q2  0 e_∞ 2.85e-01(0.0) e_L2 2.37e-01(0.0) r_2 1.46e-13(0.0) G[   2   3   4] L[  2  3  4] P[  1  1  1]
Q2  1 e_∞ 5.13e-02(2.5) e_L2 3.53e-02(2.7) r_2 1.82e-03(0.0) G[   4   6   8] L[  4  6  8] P[  1  1  1]
Q2  2 e_∞ 7.79e-03(2.7) e_L2 4.62e-03(2.9) r_2 9.77e-04(0.9) G[   8  12  16] L[  8  6  8] P[  1  2  2]
Q2  3 e_∞ 9.80e-04(3.0) e_L2 5.84e-04(3.0) r_2 4.37e-04(1.2) G[  16  24  32] L[ 16 12 16] P[  1  2  2]
F(2,3)  0: |e|_2/|u|_2 2.97e-05  |r|_2/|f|_2 4.37e-04
V(2,3)  1: |e|_2/|u|_2 1.93e-05  |r|_2/|f|_2 4.06e-06
V(2,3)  2: |e|_2/|u|_2 1.95e-05  |r|_2/|f|_2 1.88e-07
```
