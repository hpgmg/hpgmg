#!/bin/sh

test_description='Test FE coarsening and injection'

. ./hpgmg-sharness.sh

test_expect_stdout 'FE Inject fedegree=1' 4 'hpgmg-fe test-feinject -M 4,2,6 -p 2,1,2 -L 4,2,6' '
coarse u[ 0] =        0.0 at  0.0  0.0  0.0
coarse u[ 1] =        2.0 at  0.0  0.0  2.0
coarse u[ 2] =        4.0 at  0.0  0.0  4.0
coarse u[ 3] =        6.0 at  0.0  0.0  6.0
coarse u[ 4] =     2000.0 at  0.0  2.0  0.0
coarse u[ 5] =     2002.0 at  0.0  2.0  2.0
coarse u[ 6] =     2004.0 at  0.0  2.0  4.0
coarse u[ 7] =     2006.0 at  0.0  2.0  6.0
coarse u[ 8] =  2000000.0 at  2.0  0.0  0.0
coarse u[ 9] =  2000002.0 at  2.0  0.0  2.0
coarse u[10] =  2000004.0 at  2.0  0.0  4.0
coarse u[11] =  2000006.0 at  2.0  0.0  6.0
coarse u[12] =  2002000.0 at  2.0  2.0  0.0
coarse u[13] =  2002002.0 at  2.0  2.0  2.0
coarse u[14] =  2002004.0 at  2.0  2.0  4.0
coarse u[15] =  2002006.0 at  2.0  2.0  6.0
coarse u[16] =  4000000.0 at  4.0  0.0  0.0
coarse u[17] =  4000002.0 at  4.0  0.0  2.0
coarse u[18] =  4000004.0 at  4.0  0.0  4.0
coarse u[19] =  4000006.0 at  4.0  0.0  6.0
coarse u[20] =  4002000.0 at  4.0  2.0  0.0
coarse u[21] =  4002002.0 at  4.0  2.0  2.0
coarse u[22] =  4002004.0 at  4.0  2.0  4.0
coarse u[23] =  4002006.0 at  4.0  2.0  6.0
'

test_done
