/*
 *  matrix.h
 *  RVM-Speed
 *
 *  Created by Robert Lowe on 27/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <iostream>

class matrix {
    int x, y;
public:
	int cols;
	int rows;
	double* data;
	matrix();
	matrix(int,int);
	void reset(int,int);
	void resize(int,int);
	void AddColumn(const matrix &,int);
	void nme(std::string);
	void RemoveRow(int);
	void RemoveColumn(int);
	~matrix(){if(data)delete[] data;};
};


#endif