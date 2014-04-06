#!/bin/sh

test_description='Test Poisson diagonal'

. ./hpgmg-sharness.sh

test_expect_stdout 'FE Poisson diagonal fedegree=1 serial' 1 'hpgmg-fe test-opdiagonal -op_type poisson1 -M 8,12,16 -L 1,1,1' '
[0] Level 2: [   0:   8,   0:  12,   0:  16] of [   8,  12,  16] on [  1,  1,  1]
[0] Level 1: [   0:   4,   0:   6,   0:   8] of [   4,   6,   8] on [  1,  1,  1]
[0] Level 0: [   0:   2,   0:   3,   0:   4] of [   2,   3,   4] on [  1,  1,  1]
|D|_1 = 310.139  |D|_2 = 9.12568  |D|_max = 0.268519
'

test_expect_stdout 'FE Poisson diagonal fedegree=1 parallel' 4 'hpgmg-fe test-opdiagonal -op_type poisson1 -M 8,12,16 -L 1,1,1 -p 2,1,2' '
[0] Level 2: [   0:   4,   0:  12,   0:   8] of [   8,  12,  16] on [  2,  1,  2]
[0] Level 1: [   0:   2,   0:   6,   0:   4] of [   4,   6,   8] on [  2,  1,  2]
[0] Level 0: [   0:   2,   0:   3,   0:   4] of [   2,   3,   4] on [  1,  1,  1]
[1] Level 2: [   0:   4,   0:  12,   8:  16] of [   8,  12,  16] on [  2,  1,  2]
[1] Level 1: [   0:   2,   0:   6,   4:   8] of [   4,   6,   8] on [  2,  1,  2]
[2] Level 2: [   4:   8,   0:  12,   0:   8] of [   8,  12,  16] on [  2,  1,  2]
[2] Level 1: [   2:   4,   0:   6,   0:   4] of [   4,   6,   8] on [  2,  1,  2]
[3] Level 2: [   4:   8,   0:  12,   8:  16] of [   8,  12,  16] on [  2,  1,  2]
[3] Level 1: [   2:   4,   0:   6,   4:   8] of [   4,   6,   8] on [  2,  1,  2]
|D|_1 = 310.139  |D|_2 = 9.12568  |D|_max = 0.268519
'

test_done
