#!/bin/sh

test_description='Test FE creation and scatters'

. ./hpgmg-sharness.sh

test_expect_stdout 'FE GlobalToLocal fedegree=1' 4 'hpgmg-fe test-fespace -p 2,2,1 -M 4,4,2' '
Vec Object: 1 MPI processes
  type: seq
0.
1.
2.
3.
4.
5.
12.
13.
14.
6.
7.
8.
9.
10.
11.
21.
22.
23.
30.
31.
32.
33.
34.
35.
48.
49.
50.
'

test_expect_stdout 'FE Gradient/coordinates fedegree=1' 4 'hpgmg-fe test-fegrad -M 6,2,10 -p 2,1,2 -L 7,11,13' '
'

test_done
