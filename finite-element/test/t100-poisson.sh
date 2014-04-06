#!/bin/sh

test_description='Test Poisson solver'

. ./hpgmg-sharness.sh

# We expect second-order convergence on the residual for fedegree=1.
test_expect_stdout 'FE Poisson fedegree=1 serial' 1 'hpgmg-fe test-opapply -op_type poisson1 -M 4,8,12 -L 1,1,1 -poisson_solution sine' '
[0] Level 2: [   0:   4,   0:   8,   0:  12] of [   4,   8,  12] on [  1,  1,  1]
[0] Level 1: [   0:   2,   0:   4,   0:   6] of [   2,   4,   6] on [  1,  1,  1]
[0] Level 0: [   0:   1,   0:   2,   0:   3] of [   1,   2,   3] on [  1,  1,  1]
|A u - F|_max/|F|_max = 0.0978195
'

test_expect_stdout 'FE Poisson fedegree=1 serial refined' 1 'hpgmg-fe test-opapply -op_type poisson1 -M 8,16,24 -L 1,1,1 -poisson_solution sine' '
[0] Level 3: [   0:   8,   0:  16,   0:  24] of [   8,  16,  24] on [  1,  1,  1]
[0] Level 2: [   0:   4,   0:   8,   0:  12] of [   4,   8,  12] on [  1,  1,  1]
[0] Level 1: [   0:   2,   0:   4,   0:   6] of [   2,   4,   6] on [  1,  1,  1]
[0] Level 0: [   0:   1,   0:   2,   0:   3] of [   1,   2,   3] on [  1,  1,  1]
|A u - F|_max/|F|_max = 0.0253888
'

test_expect_stdout 'FE Poisson fedegree=1 parallel refined' 4 'hpgmg-fe test-opapply -op_type poisson1 -M 8,16,24 -L 1,1,1 -p 1,2,2 -cmax 48 -poisson_solution sine' '
[0] Level 3: [   0:   8,   0:   8,   0:  12] of [   8,  16,  24] on [  1,  2,  2]
[0] Level 2: [   0:   4,   0:   4,   0:   6] of [   4,   8,  12] on [  1,  2,  2]
[0] Level 1: [   0:   2,   0:   4,   0:   6] of [   2,   4,   6] on [  1,  1,  1]
[0] Level 0: [   0:   1,   0:   2,   0:   3] of [   1,   2,   3] on [  1,  1,  1]
[1] Level 3: [   0:   8,   0:   8,  12:  24] of [   8,  16,  24] on [  1,  2,  2]
[1] Level 2: [   0:   4,   0:   4,   6:  12] of [   4,   8,  12] on [  1,  2,  2]
[2] Level 3: [   0:   8,   8:  16,   0:  12] of [   8,  16,  24] on [  1,  2,  2]
[2] Level 2: [   0:   4,   4:   8,   0:   6] of [   4,   8,  12] on [  1,  2,  2]
[3] Level 3: [   0:   8,   8:  16,  12:  24] of [   8,  16,  24] on [  1,  2,  2]
[3] Level 2: [   0:   4,   4:   8,   6:  12] of [   4,   8,  12] on [  1,  2,  2]
|A u - F|_max/|F|_max = 0.0253888
'

test_done
