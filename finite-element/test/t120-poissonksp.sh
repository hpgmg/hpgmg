#!/bin/sh

test_description='Test Poisson solve using KSP'

. ./hpgmg-sharness.sh

# Error norm is converging at second order
test_expect_stdout 'FE Poisson KSP solve fedegree=1 serial' 1 'hpgmg-fe test-kspsolve -op_type poisson1 -M 8,12,16 -L 1,1,1 -ksp_converged_reason -ksp_view -ksp_type chebyshev -ksp_chebyshev_eigenvalues 0.2,2 -pc_type jacobi -poisson_solution sine' '
[0] Level 2: [   0:   8,   0:  12,   0:  16] of [   8,  12,  16] on [  1,  1,  1]
[0] Level 1: [   0:   4,   0:   6,   0:   8] of [   4,   6,   8] on [  1,  1,  1]
[0] Level 0: [   0:   2,   0:   3,   0:   4] of [   2,   3,   4] on [  1,  1,  1]
Linear solve converged due to CONVERGED_RTOL iterations 20
KSP Object: 1 MPI processes
  type: chebyshev
    eigenvalue estimates used:  min = 0.2, max = 2.
  maximum iterations=10000, initial guess is zero
  tolerances:  relative=1e-05, absolute=1e-50, divergence=10000.
  left preconditioning
  using PRECONDITIONED norm type for convergence test
PC Object: 1 MPI processes
  type: jacobi
  linear system matrix = precond matrix:
  Mat Object: 1 MPI processes
    type: shell
    rows=1989, cols=1989
|v-u|_2/|u|_2 = 0.0393899
'

test_expect_stdout 'FE Poisson KSP solve fedegree=1 parallel' 4 'hpgmg-fe test-kspsolve -op_type poisson1 -M 8,12,16 -L 1,1,1 -ksp_converged_reason -ksp_view -ksp_type chebyshev -ksp_chebyshev_eigenvalues 0.2,2 -pc_type jacobi -p 1,2,2 -poisson_solution sine' '
[0] Level 2: [   0:   8,   0:   6,   0:   8] of [   8,  12,  16] on [  1,  2,  2]
[0] Level 1: [   0:   4,   0:   3,   0:   4] of [   4,   6,   8] on [  1,  2,  2]
[0] Level 0: [   0:   2,   0:   3,   0:   4] of [   2,   3,   4] on [  1,  1,  1]
[1] Level 2: [   0:   8,   0:   6,   8:  16] of [   8,  12,  16] on [  1,  2,  2]
[1] Level 1: [   0:   4,   0:   3,   4:   8] of [   4,   6,   8] on [  1,  2,  2]
[2] Level 2: [   0:   8,   6:  12,   0:   8] of [   8,  12,  16] on [  1,  2,  2]
[2] Level 1: [   0:   4,   3:   6,   0:   4] of [   4,   6,   8] on [  1,  2,  2]
[3] Level 2: [   0:   8,   6:  12,   8:  16] of [   8,  12,  16] on [  1,  2,  2]
[3] Level 1: [   0:   4,   3:   6,   4:   8] of [   4,   6,   8] on [  1,  2,  2]
Linear solve converged due to CONVERGED_RTOL iterations 20
KSP Object: 4 MPI processes
  type: chebyshev
    eigenvalue estimates used:  min = 0.2, max = 2.
  maximum iterations=10000, initial guess is zero
  tolerances:  relative=1e-05, absolute=1e-50, divergence=10000.
  left preconditioning
  using PRECONDITIONED norm type for convergence test
PC Object: 4 MPI processes
  type: jacobi
  linear system matrix = precond matrix:
  Mat Object: 4 MPI processes
    type: shell
    rows=1989, cols=1989
|v-u|_2/|u|_2 = 0.0393899
'

test_done
