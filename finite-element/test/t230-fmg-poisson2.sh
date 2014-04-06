#!/bin/sh

test_description='Test Q2 Poisson solve using FMG'

. ./hpgmg-sharness.sh

# Error norm is converging at fourth order (superconvergent at Lagrange nodes)
test_expect_stdout 'FE Poisson FMG solve fedegree=2 serial' 1 'hpgmg-fe fmg -op_type poisson2 -M 4,4,6 -smooth 4,3' '
F(4,3)  0: |e|_2/|u|_2 9.08e-03  |r|_2/|f|_2 3.35e-04
V(4,3)  1: |e|_2/|u|_2 9.17e-03  |r|_2/|f|_2 8.27e-07
V(4,3)  2: |e|_2/|u|_2 9.17e-03  |r|_2/|f|_2 5.54e-09
'

test_expect_stdout 'FE Poisson FMG solve fedegree=2 parallel' 4 'hpgmg-fe fmg -op_type poisson2 -M 4,4,6 -smooth 4,3 -p 1,2,2' '
F(4,3)  0: |e|_2/|u|_2 9.08e-03  |r|_2/|f|_2 3.35e-04
V(4,3)  1: |e|_2/|u|_2 9.17e-03  |r|_2/|f|_2 8.27e-07
V(4,3)  2: |e|_2/|u|_2 9.17e-03  |r|_2/|f|_2 5.54e-09
'

test_done
