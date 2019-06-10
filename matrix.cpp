/*
 *  matrix.cpp
 *  RVM-Speed
 *
 *  Created by Robert Lowe on 27/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "matrix.h"

matrix::matrix(int a, int b){
	rows=a;
	cols=b;
	data = new double[a*b];
}



matrix::matrix(){
	rows=0;
	cols=0;
	data = new double[1];
};

void matrix::reset(int a ,int b){
	if(data)delete[] data;
	rows=a;
	cols=b;
	data = new double[a*b];
}	


void matrix::resize(int a ,int b){
	double *tmpcpy=new double[a*b];
	
	for (int i=0; i<rows; i++) {
		for (int k=0; k<cols; k++) {
			tmpcpy[i*cols+k]=data[i*cols+k];
		}
	}
	rows=a;
	cols=b;
	delete[] data;
	data = tmpcpy;
}	
void  matrix::AddColumn(const matrix &A,int column){
	

	if (rows!=0 && cols!=0){
		double *tmpcpy=new double[rows*(cols+1)];
		for (int i=0; i<rows; i++) {
			for (int k=0; k<cols; k++) {
				tmpcpy[i*(cols+1)+k]=data[i*cols+k];
			}
			tmpcpy[i*(cols+1)+(cols)]=A.data[i*A.cols+column];
		}
		cols+=1;
		delete[] data;
		data=tmpcpy;

	}
	else{
		delete[] data;
		rows=A.rows;
		cols=1;
		data = new double[rows*cols];
		for	(int i=0; i<A.rows; i++){
			data[i]=A.data[(i*A.cols)+column];
		}
	}
}

void matrix::RemoveColumn(int a){
	
	double *tmpcpy=new double[rows*(cols-1)];
	
	for (int i=0; i<rows; i++) {
		for (int k=0; k<cols-1; k++){
			if(k<a)
				tmpcpy[i*(cols-1)+k]=data[i*cols+k];
			else
 				tmpcpy[i*(cols-1)+k]=data[i*cols+(k+1)];
		}

	}
	
	cols-=1;
	delete[] data;
	data=tmpcpy;
}

void matrix::RemoveRow(int a){
	
	double *tmpcpy=new double[(rows-1)*cols];

	for (int i=0; i<(rows-1); i++) {
		for (int k=0; k<cols; k++){
			if(i<a)
				tmpcpy[i*cols+k]=data[i*cols+k];
			else
 				tmpcpy[i*cols+k]=data[(i+1)*cols+k];
		}
		
	}
	
	rows-=1;
	delete[] data;
	data=tmpcpy;
}

	

/*matrix matrix::diag(){
	
	matrix out(rows,rows);
	
    for (unsigned i = 0; i < rows; ++ i){
		for (unsigned j = 0; j <rows; ++ j){
			if(i==j)
				out.data[(i*rows)+j]=data[i];
			else
				out.data[(i*rows)+j]=0.0;
		}
	}
	
	
	return out;
}
 */
/*
matrix matrix::operator=(const matrix &a){
	matrix temp(a.rows,a.cols);
	
	for(int i=0; i<a.rows; i++){
		for(int k=0; k<a.cols; k++){
			temp.data[i*a.cols+k]=a.data[i*a.cols+k];
			std:: cout << i << " " <<k << std::endl;
		}
	}
	std::cout << "USING OVERLOAD" << std::endl;
	return(temp);
}
*/