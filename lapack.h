/*
 *  lapack.h
 *  RVM-Rob
 *
 *  Created by Robert Lowe on 23/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <Accelerate/Accelerate.h>
#include "matrix.h"

int matrixprod(const matrix &A,const matrix &B, matrix &C,int trans,double alpha);

int matrixvprod(const matrix &A,const matrix &B, matrix &C, int trans, double alpha);
int linalg(const matrix &A,const matrix &B, matrix &C, int trans);
int chol(const matrix &A, matrix &B);
int inverse(const matrix &A, matrix&B);