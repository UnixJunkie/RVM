/*
 *  lapack.cpp
 *  RVM-Rob
 *
 *  Created by Robert Lowe on 23/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "lapack.h"
#include <iostream>

using namespace std;

int matrixprod(const matrix &A,const matrix &B,matrix &C, int trans,double alpha){
	
	int nrows=A.rows;
	int ncols=A.cols;
	int bcols=B.cols;
	int brows=B.rows;
	
	int lda=A.cols;
	int ldb=B.cols;
	

	if (trans==0) {
		C.reset(nrows,bcols);
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nrows, bcols, ncols, alpha, A.data, lda, B.data, ldb, 0.0, C.data, bcols);
		
	}
	
	else if(trans==1){

		C.reset(ncols,bcols);

		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, ncols, bcols, brows, alpha, A.data, lda, B.data, ldb, 0.0, C.data, bcols);
		

	}
	
	else if(trans==2){
		C.reset(nrows,brows);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nrows, brows, ncols, alpha, A.data, lda, B.data, ldb, 0.0, C.data, brows);
		
	
	}
	else if(trans==3){
		C.reset(ncols,brows);
		
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, ncols, brows, nrows, alpha, A.data, lda, B.data, ldb, 0.0, C.data, brows);
		
	}
	
	return 0;
}


int matrixvprod(const matrix &A,const  matrix &B,matrix &C, int trans, double alpha){
	
	int nrows=A.rows;
	int ncols=A.cols;
	int lda=A.cols;
	
	if (trans==0){
		C.reset(nrows,1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, nrows, ncols, alpha, A.data, lda, B.data, 1, 0.0, C.data, 1);
	}
	else if (trans==1){
		C.reset(ncols,1);
		cblas_dgemv(CblasRowMajor, CblasTrans, nrows, ncols, alpha, A.data, lda, B.data, 1, 0.0, C.data, 1);
	}
	
	return(0);
	
}


int linalg(const matrix &A,const matrix &B,matrix &C, int trans){

	matrix Atemp(A.rows,A.cols);
	matrix Btemp(B.rows,B.cols);
	
	for (int i=0; i<A.rows; i++) {
			for (int k=0; k<A.cols; k++) {
				Atemp.data[(i*A.cols)+k]=A.data[(i*A.cols)+k];
			}
		}

	
	for (int i=0; i<B.rows; i++) {
		for (int k=0; k<B.cols; k++) {
			Btemp.data[(i*B.cols)+k]=B.data[(i*B.cols)+k];
		}
	}
	__CLPK_integer nrowsA=Atemp.rows;
	__CLPK_integer ncolsA=Atemp.cols;
	
	__CLPK_integer nrhs=Btemp.cols;
	
	__CLPK_integer lda = Atemp.rows;
	__CLPK_integer ldb = Btemp.rows;
	__CLPK_integer info;
	//Not sure on what work is so set size to A
	double work[A.rows*A.cols];
	__CLPK_integer lwork=A.rows*A.cols;
	

	if (trans==0){
	char T='N';
	dgels_(&T, &nrowsA, &ncolsA, &nrhs, Atemp.data, &lda, Btemp.data, &ldb, work, &lwork, &info);
	
	C.reset(A.cols,1);
	
	for(int i=0; i<A.cols; i++){
			C.data[i]=Btemp.data[i];
	}
	
	}
	else {
		char T='T';
		dgels_(&T, &ncolsA, &nrowsA, &nrhs, Atemp.data, &lda, Btemp.data, &ldb, work, &lwork, &info);
		C.reset(A.rows,1);
		
		for(int i=0; i<A.rows; i++){
			C.data[i]=Btemp.data[i];
		}
		
	}
	return(0);

}


int chol(const matrix &A,matrix &B){
	
	__CLPK_integer N=A.rows;
	char T='L';
	__CLPK_integer info;
	
	B.reset(A.rows,A.cols);
	
	for (int i=0; i<A.rows; i++) {
		for (int k=0; k<A.cols; k++) {
			B.data[(i*A.cols)+k]=A.data[(i*A.cols)+k];
		}
	}
	
	
	dpotrf_(&T, &N, B.data, &N, &info);
	
	for (int i=1; i<N; i++) {
		for (int k=0; k<i; k++) {
			B.data[i*B.cols+k]=0.0;
		}
	}
	
	return(info);
}


int inverse(const matrix &A, matrix&B){
	
	__CLPK_integer M=A.rows;
	__CLPK_integer N=A.cols;
	
	__CLPK_integer  *IPIV=new __CLPK_integer[A.rows];
	
	B.reset(A.rows,A.cols);
	
	for (int i=0; i<A.rows; i++) {
		for (int k=0; k<A.cols; k++) {
			B.data[(i*A.cols)+k]=A.data[(i*A.cols)+k];
		}
	}
	
	__CLPK_integer info;
	double *work = new double[A.rows*A.cols];
	__CLPK_integer lwork=A.rows*A.cols;

	
	dgetrf_(&M, &N, B.data, &M, IPIV, &info);
	dgetri_(&N, B.data, &M, IPIV, work, &lwork, &info);
	
	delete[] IPIV;
	delete[] work;
	
	return(info);
}
	
	
	
