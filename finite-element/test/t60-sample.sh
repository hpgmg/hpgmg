#!/bin/sh

test_description='Test grid creation'

. ./hpgmg-sharness.sh

test_expect_stdout 'list of samples' 1 'hpgmg-fe test-sampler -local 100,1e9 -maxsamples 10 -nranks 192' '
Processors: [   4    6    8] = 192
Filtered Grid: L11 [4096 6144 6144] = 154618822656
Filtered Grid: L12 [4096 4096 4096] =  68719476736
Filtered Grid: L10 [2048 2048 3072] =  12884901888
Filtered Grid: L 9 [1024 1024 1536] =   1610612736
Filtered Grid: L 8 [ 512  512  768] =    201326592
Filtered Grid: L 7 [ 256  256  384] =     25165824
Filtered Grid: L 6 [ 128  128  192] =      3145728
Filtered Grid: L 5 [  64   64   96] =       393216
Filtered Grid: L 4 [  32   32   48] =        49152
Filtered Grid: L 3 [  24   24   32] =        18432
'

test_done
