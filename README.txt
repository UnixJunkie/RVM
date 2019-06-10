1) RVM SOURCE
2) DESCRIPTOR CALCULATION
3) MDDR DATA

***************************************************************************
***************************************************************************
1) RVM SOURCE

Currently there is only support for Mac OS X as we are using the Accelerate Framework for the BLAS libraries. 

To Compile go to RVM Source Directory and type
g++ main.cpp matrix.cpp SparseBayes.cpp lapack.cpp fullstatistics.cpp -o RVM -framework accelerate -O3

This will compile a binary for your system named RVM

You can then run the program as 
./RVM

Arguments:
-k STRING Must be [Gaus,Cauch,Binary] This defines which kernel to use
-b DOUBLE This defines the hyperparameter for that kernel
-t STRING Directory of training file. The file must be tab separated and have the last column as the class label 1 or 0
-v STRING Directory of Validation file
-s STRING Directory of Test file
-i INT	  Iteration limit for algorithm
-o STRING Directory of output file
-r STRING Not used
-y	  Output probabilities

******************************************************************************************************************************************************
2) DESCRIPTOR CALCULATION

The descriptors used were MOLPRINT 2D descriptors. These can be calculated using the free available software from http://molprint.com/

******************************************************************************************************************************************************
3) MDDR DATA

The MDDR data is licensed so can not be included in the supporting information. We can supply the MDDR Codes used and these are attached in a separate file MDDR_IDS

******************************************************************************************************************************************************
