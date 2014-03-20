#!/bin/sh

test_description='Test FE restriction'

. ./fefas-sharness.sh

test_expect_stdout 'FE Restriction fedegree=1 serial' 1 'fefas test-ferestrict -M 4,4,6 -L 4,4,6' '
|u_c - I_h^H u_f|_max =     0
'

test_expect_stdout 'FE Interpolation fedegree=1 parallel' 4 'fefas test-feinterp -M 6,4,10 -L 6,4,10 -p 2,1,2' '
|u - I Ihat u|_max =     0
'

test_done
