/*
 *  SparseBayes.h
 *  RVM-Speed
 *
 *  Created by Robert Lowe on 27/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include "matrix.h"
#include <string>
#include <vector>
#include <iostream>
#include<fstream>
#include<math.h>
#include "lapack.h"
#include "fullstatistics.h"
#include <mach/mach_time.h>


using namespace std;

typedef std::vector<double> DOUBLE;
typedef std::vector<string> LINE;

template<class T> struct index_cmp {
	index_cmp(const T arr) : arr(arr) {}
	bool operator()(const size_t a, const size_t b) const
	{ return arr[a] < arr[b]; }
	const T arr;
};


void kernelfunction_gauss(matrix &BASIS,const std::vector<DOUBLE> &data,double basisWidth,const std::vector<int> &dataclass,matrix &Targets);
void kernelfunction_cauch(matrix &BASIS,const std::vector<DOUBLE> &data,double basisWidth,const std::vector<int> &dataclass,matrix &Targets);
void kernelfunction_binary(matrix &BASIS,const std::vector<DOUBLE> &data,double basisWidth,const std::vector<int> &dataclass,matrix &Targets);
double SparseBayes(int LIKELIHOOD,int ItNum,int monitor_its,double MinDeltaLogAlpha,double AlignmentMax,std::vector<DOUBLE> data,std::vector<int> dataclass,int PriorityAddition,int PriorityDeletion
				,int BasisAlignmentTest,std::vector<int> &PARAMATERrev,matrix &PARAMATERval,int kernel,double basisWidth);
void Initialisation(int LIKELIHOOD, matrix &BASIS, std::vector<int> Targets,matrix &Alpha,matrix &beta,matrix &Mu, matrix &Phi, std::vector<int> &Used, matrix &Scales);


