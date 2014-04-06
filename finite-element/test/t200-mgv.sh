#!/bin/sh

test_description='Test Poisson solve using MG V-cycles'

. ./hpgmg-sharness.sh

# Error norm is converging at second order
test_expect_stdout 'FE Poisson MG V-cycle solve fedegree=1 serial' 1 'hpgmg-fe mgv -op_type poisson1 -M 16,20,24 -L 1,1,1 -smooth 2,2 -mg_eig_target 2,0.2 -poisson_solution sine' '
[0] Level 2: [   0:  16,   0:  20,   0:  24] of [  16,  20,  24] on [  1,  1,  1]
[0] Level 1: [   0:   8,   0:  10,   0:  12] of [   8,  10,  12] on [  1,  1,  1]
[0] Level 0: [   0:   4,   0:   5,   0:   6] of [   4,   5,   6] on [  1,  1,  1]
V(2,2)  0: |e|_2/|u|_2 1.50e-02  |r|_2/|f|_2 2.25e-01
V(2,2)  1: |e|_2/|u|_2 1.06e-02  |r|_2/|f|_2 5.40e-02
V(2,2)  2: |e|_2/|u|_2 1.27e-02  |r|_2/|f|_2 1.31e-02
V(2,2)  3: |e|_2/|u|_2 1.34e-02  |r|_2/|f|_2 3.22e-03
V(2,2)  4: |e|_2/|u|_2 1.35e-02  |r|_2/|f|_2 7.91e-04
'

test_expect_stdout 'FE Poisson MG V-cycle solve fedegree=1 parallel' 4 'hpgmg-fe mgv -op_type poisson1 -M 16,20,24 -L 1,1,1 -p 1,2,2 -cmax 240 -smooth 2,2 -mg_eig_target 2,0.2 -poisson_solution sine' '
[0] Level 2: [   0:  16,   0:  10,   0:  12] of [  16,  20,  24] on [  1,  2,  2]
[0] Level 1: [   0:   8,   0:   5,   0:   6] of [   8,  10,  12] on [  1,  2,  2]
[0] Level 0: [   0:   4,   0:   5,   0:   6] of [   4,   5,   6] on [  1,  1,  1]
[1] Level 2: [   0:  16,   0:  10,  12:  24] of [  16,  20,  24] on [  1,  2,  2]
[1] Level 1: [   0:   8,   0:   5,   6:  12] of [   8,  10,  12] on [  1,  2,  2]
[2] Level 2: [   0:  16,  10:  20,   0:  12] of [  16,  20,  24] on [  1,  2,  2]
[2] Level 1: [   0:   8,   5:  10,   0:   6] of [   8,  10,  12] on [  1,  2,  2]
[3] Level 2: [   0:  16,  10:  20,  12:  24] of [  16,  20,  24] on [  1,  2,  2]
[3] Level 1: [   0:   8,   5:  10,   6:  12] of [   8,  10,  12] on [  1,  2,  2]
V(2,2)  0: |e|_2/|u|_2 1.50e-02  |r|_2/|f|_2 2.25e-01
V(2,2)  1: |e|_2/|u|_2 1.06e-02  |r|_2/|f|_2 5.40e-02
V(2,2)  2: |e|_2/|u|_2 1.27e-02  |r|_2/|f|_2 1.31e-02
V(2,2)  3: |e|_2/|u|_2 1.34e-02  |r|_2/|f|_2 3.22e-03
V(2,2)  4: |e|_2/|u|_2 1.35e-02  |r|_2/|f|_2 7.91e-04
'

test_done
