# RVM
Ranking Vector Machine

supplementary materials of:
"Lowe, R., Mussa, H. Y., Mitchell, J. B., & Glen, R. C. (2011). Classifying molecules using a sparse probabilistic kernel binary classifier. Journal of chemical information and modeling, 51(7), 1539-1544."

SI: https://pubs.acs.org/doi/full/10.1021/ci200128w

zip: https://pubs.acs.org/doi/suppl/10.1021/ci200128w/suppl_file/ci200128w_si_001.zip

***************************************************************************

1) RVM SOURCE
2) DESCRIPTOR CALCULATION
3) MDDR DATA

***************************************************************************
1) RVM SOURCE

Currently there is only support for Mac OS X as we are using the Accelerate Framework for the BLAS libraries.

To Compile, type make.
This will compile a binary for your system named RVM

You can then run the program as

./RVM

Arguments:

-k STRING Must be [Gaus,Cauch,Binary] This defines which kernel to use

-b DOUBLE This defines the hyperparameter for that kernel

-t STRING Directory of training file. The file must be tab separated and have the last column as the class label 1 or 0

-v STRING Directory of Validation file

-s STRING Directory of Test file

-i INT    Iteration limit for algorithm

-o STRING Directory of output file

-r STRING Not used

-y        Output probabilities

***************************************************************************
2) DESCRIPTOR CALCULATION

The descriptors used were MOLPRINT 2D descriptors. These can be calculated using the free available software from http://molprint.com/

***************************************************************************
3) MDDR DATA

The MDDR data is licensed so can not be included in the supporting information. We can supply the MDDR Codes used and these are attached in a separate file MDDR_IDS
