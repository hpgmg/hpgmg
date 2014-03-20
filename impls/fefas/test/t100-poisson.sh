#!/bin/sh

test_description='Test Poisson solver'

. ./fefas-sharness.sh

test_expect_stdout 'FE Poisson fedegree=1 serial' 1 'fefas fmg -op_type poisson1 -M 4,6,4 -L 2,2,2' '
Create poisson1
[0] Level 1: [   0:   4,   0:   6,   0:   4] of [   4,   6,   4] on [  1,  1,  1]
[0] Level 0: [   0:   2,   0:   3,   0:   2] of [   2,   3,   2] on [  1,  1,  1]
'

test_done
