RVM:
	g++ main.cpp matrix.cpp SparseBayes.cpp lapack.cpp fullstatistics.cpp -o RVM -framework accelerate -O3

clean:
	rm -f RVM
