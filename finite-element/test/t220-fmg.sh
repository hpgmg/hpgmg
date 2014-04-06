#!/bin/sh

test_description='Test Poisson solve using FMG'

. ./hpgmg-sharness.sh

# Error norm is converging at second order
test_expect_stdout 'FE Poisson FMG solve fedegree=1 serial' 1 'hpgmg-fe fmg -op_type poisson1 -M 8,16,24 -smooth 3,3 -mg_eig_target 2,0.2 -poisson_solution sine' '
F(3,3)  0: |e|_2/|u|_2 2.26e-02  |r|_2/|f|_2 3.37e-02
V(3,3)  1: |e|_2/|u|_2 2.58e-02  |r|_2/|f|_2 2.05e-03
V(3,3)  2: |e|_2/|u|_2 2.60e-02  |r|_2/|f|_2 1.25e-04
'

test_expect_stdout 'FE Poisson FMG solve fedegree=1 parallel' 4 'hpgmg-fe fmg -op_type poisson1 -M 8,16,24 -p 1,2,2 -smooth 3,3 -mg_eig_target 2,0.2 -poisson_solution sine' '
F(3,3)  0: |e|_2/|u|_2 2.26e-02  |r|_2/|f|_2 3.37e-02
V(3,3)  1: |e|_2/|u|_2 2.58e-02  |r|_2/|f|_2 2.05e-03
V(3,3)  2: |e|_2/|u|_2 2.60e-02  |r|_2/|f|_2 1.25e-04
'

test_done
