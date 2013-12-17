//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

//------------------------------------------------------------------------------------------------------------------------------
#include "timer.h"
#include "defines.h"
#include "box.h"
#include "mg.h"
#include "operators.h"
//------------------------------------------------------------------------------------------------------------------------------
#include "operators.ompif/exchange_boundary.c"
#include "operators.27pt/lambda.c"
#include "operators.27pt/jacobi.c"
#include "operators.27pt/apply_op.c"
#include "operators.27pt/residual.c"
#include "operators.ompif/restriction.c"
#include "operators.ompif/interpolation.c"
#include "operators.ompif/misc.c"
#include "operators.ompif/matmul.c"
#include "operators.ompif/problem0.c"
//------------------------------------------------------------------------------------------------------------------------------
