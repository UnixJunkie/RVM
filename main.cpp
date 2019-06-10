#include "SparseBayes.h"
#include "matrix.h"
#include <iostream>
#include <fstream>

typedef std::vector<double> DOUBLE;
typedef std::vector<string> LINE;

static mach_timebase_info_data_t	sTimebaseInfo;

double 
timeInMilliseconds_main(uint64_t time)
{    
	if ( sTimebaseInfo.denom == 0 ) {
        (void) mach_timebase_info(&sTimebaseInfo);
    }
	
    return (time * sTimebaseInfo.numer / sTimebaseInfo.denom)/1000000;
}

void readfile(string filename, std::vector<DOUBLE> &data, std::vector<int> &dataclass){
	
	string deliminator="\t";
	std::string line;
	int pos;
	
	ifstream myfile (filename.c_str());
	if (myfile.is_open()){
		while(getline(myfile,line)) /* read a record */
		{
			LINE ln;
			std::vector<double> ln2;
			while( (pos = line.find(deliminator)) > 0)
			{
				string field = line.substr(0,pos);
				line = line.substr(pos+1);
				double test= strtod(field.c_str(),NULL);
				ln2.push_back(test);
				
			}
			data.push_back(ln2);
			string field = line.substr(0,pos);
			line = line.substr(pos+1);
			dataclass.push_back(atoi(field.c_str()));
		}
		myfile.close();
	}
	else {
		std::cout << "Unable to open file";
		std::cout << " " << filename;
		exit(1);
	}
}

void Sigmoid(matrix &A,bool yout){
	
	ofstream myfile_y;
	if(yout){
		string output2="output_y.txt";
		myfile_y.open(output2.c_str(),ios::out);
	}

	for (int i=0; i<A.rows; i++){
		if(yout)
			myfile_y << (1.0/(1.0+exp(-A.data[i]))) << endl;
		if( (1.0/(1.0+exp(-A.data[i])))>0.5)
			A.data[i]=1.0;
		else
			A.data[i]=0.0;
	}	
	
	if (yout)
		myfile_y.close();
	
}



void kernelfunction_cauch_test(matrix &BASIS,const std::vector<DOUBLE> &data,std::vector<DOUBLE> &datatest,double basisWidth,std::vector<int> PARAMETERrev){
	
	//Calculate Basis
	matrix X2(datatest.size(),PARAMETERrev.size());
	matrix Y2(datatest.size(),PARAMETERrev.size());
	
	matrix X(datatest.size(),data[0].size());
	matrix Y(PARAMETERrev.size(),data[0].size());
	
	
	for (int i=0; i<datatest.size(); i++) {
		double sumofsquares=0.0;
		for(int k=0; k<datatest[i].size(); k++){
			sumofsquares+=pow(datatest[i][k],2);
			X.data[(i*datatest[i].size())+k]=datatest[i][k];
		}
		for(int p=0; p<X2.cols; p++){
			X2.data[(i*X2.cols)+p]=sumofsquares;
		}
	}
	
	
	for (int i=0; i<PARAMETERrev.size(); i++) {
		double sumofsquares=0.0;
		for (int k=0; k<data[0].size(); k++) {
			sumofsquares+=(data[PARAMETERrev[i]][k]*data[PARAMETERrev[i]][k]);
			Y.data[i*Y.cols+k]=data[PARAMETERrev[i]][k];

		}
		for(int p=0; p<X2.rows; p++){
			Y2.data[p*X2.cols+i]=sumofsquares;
		}
	}

	matrixprod(X, Y,BASIS, 2,-2.0);
	
	
	for (int i=0; i<BASIS.rows; i++) {
		for (int k=0; k<BASIS.cols; k++) {
			BASIS.data[(i*BASIS.cols)+k]=1.0/(1.0+((X2.data[(i*BASIS.cols)+k]+Y2.data[(i*BASIS.cols)+k]+BASIS.data[(i*BASIS.cols)+k])*basisWidth));
		}
	}

	
}

void kernelfunction_gauss_test(matrix &BASIS,const std::vector<DOUBLE> &data,std::vector<DOUBLE> &datatest,double basisWidth,std::vector<int> PARAMETERrev){
	
	//Calculate Basis
	matrix X2(datatest.size(),PARAMETERrev.size());
	matrix Y2(datatest.size(),PARAMETERrev.size());
	
	matrix X(datatest.size(),data[0].size());
	matrix Y(PARAMETERrev.size(),data[0].size());
	
	
	for (int i=0; i<datatest.size(); i++) {
		double sumofsquares=0.0;
		for(int k=0; k<datatest[i].size(); k++){
			sumofsquares+=pow(datatest[i][k],2);
			X.data[(i*datatest[i].size())+k]=datatest[i][k];
		}
		for(int p=0; p<X2.cols; p++){
			X2.data[(i*X2.cols)+p]=sumofsquares;
		}
	}
	
	
	for (int i=0; i<PARAMETERrev.size(); i++) {
		double sumofsquares=0.0;
		for (int k=0; k<data[0].size(); k++) {
			sumofsquares+=(data[PARAMETERrev[i]][k]*data[PARAMETERrev[i]][k]);
			Y.data[i*Y.cols+k]=data[PARAMETERrev[i]][k];
			
		}
		for(int p=0; p<X2.rows; p++){
			Y2.data[p*X2.cols+i]=sumofsquares;
		}
	}
	
	matrixprod(X, Y,BASIS, 2,-2.0);
	
	
	for (int i=0; i<BASIS.rows; i++) {
		for (int k=0; k<BASIS.cols; k++) {
			BASIS.data[(i*BASIS.cols)+k]=exp((X2.data[(i*BASIS.cols)+k]+Y2.data[(i*BASIS.cols)+k]+BASIS.data[(i*BASIS.cols)+k])*-basisWidth);
		}
	}
	
	
}

void kernelfunction_bin_test(matrix &BASIS,const std::vector<DOUBLE> &data,std::vector<DOUBLE> &datatest,double basisWidth,std::vector<int> PARAMETERrev){
	
	BASIS.reset(datatest.size(), PARAMETERrev.size());

	for (int i=0; i<datatest.size(); i++) {
		for (int j=0; j<PARAMETERrev.size(); j++) {
			double sum=0.0;
			for (int k=0; k<data[0].size(); k++) {
				sum+=abs(datatest[i][k]-data[PARAMETERrev[j]][k]);
			}
			BASIS.data[i*BASIS.cols+j]=pow(basisWidth,(data[0].size()-sum))*pow(1-basisWidth,sum);
		}
	}
}


int main (int argc, char * const argv[]) {
	
	
	int ItNum=1000;
	double MinDeltaLogAlpha=1e-3,MinDeltaLogBeta = 1e-6,AlignmentMax=1-1e-3;
	//Reporting on the iterations
	int monitor_its=10;
	
	bool PriorityAddition=0,PriorityDeletion=1,BasisAlignmentTest=1;
	
	double basisWidth=0.015625;
	//kernel 1- Gaus 2- Cauch 3- Binary
	string kernel="Gaus";
	string train="";
	string test="";
	string vals="";
	int kern=0;
	string runno="";
	
	string output="output_rvm.txt";
	
	bool yout=0;
	
	char ch;

	while ((ch = getopt(argc, argv, "k:b:t:v:s:i:r:o:y")) != -1) {
		switch (ch) {
			case 'k':
				kernel=optarg;
				cout << "Setting kernel to " << kernel << endl; 
				break;
			case 'b':
				basisWidth=atof(optarg);
				break;
			case 't':
				train=optarg;
				break;
			case 'v':
				vals=optarg;
				break;				
			case 's':
				test=optarg;
				break;
			case 'i':
				ItNum=atoi(optarg);
				break;
			case 'o':
				output=optarg;
				break;
			case 'r':
				runno.append("c");
				runno.append(optarg);
				break;
			case 'y':
				yout=1;
				break;
		}
	}
				
	
	std::vector<int> PARAMETERrev;
	matrix PARAMETERval;
	

	if(kernel=="Gaus")
		kern=1;
	else if(kernel=="Cauch")
		kern=2;
	else if (kernel=="Binary")
		kern=3;
	else {
		cout << "Kernel not known. Choose Gaus, Cauch or Binary" << endl;
		exit(0);
	}

	//READ in data
	std::vector<DOUBLE> data;
	std::vector<int> dataclass;
	
	string s=train;
	cout << s << endl;

	readfile(s,data,dataclass);
	double timeran;
	
	//First number is likelihood - 1) Bernoulli in this case, 0) Gaussian.
	timeran=SparseBayes(1,ItNum,monitor_its,MinDeltaLogAlpha,AlignmentMax,data,dataclass, PriorityAddition, PriorityDeletion, BasisAlignmentTest,PARAMETERrev,PARAMETERval,kern,basisWidth);
	
	//READ in data
	std::vector<DOUBLE> datatest;
	std::vector<int> datatestclass;
	
	s=vals;
	readfile(s,datatest,datatestclass);
	matrix BASIS;
	if (kern==1)
		kernelfunction_gauss_test(BASIS, data, datatest, basisWidth, PARAMETERrev);	
	else if(kern==2)
		kernelfunction_cauch_test(BASIS, data, datatest, basisWidth, PARAMETERrev);	
	else if(kern==3)
		kernelfunction_bin_test(BASIS, data, datatest, basisWidth, PARAMETERrev);	

	matrix y;
	matrixprod(BASIS, PARAMETERval, y, 0, 1.0);
	
	Sigmoid(y,yout);
	

	
	ofstream myfile (output.c_str(),ios::app);

	if (myfile.is_open())
	{
		double TP=0,FP=0,TN=0,FN=0;
	
		for (int i=0; i<datatestclass.size(); i++) {
			if (datatestclass[i]==1 and y.data[i]==1)
				TP+=1;
			else if (datatestclass[i]==0 and y.data[i]==1)
				FP+=1;
			else if (datatestclass[i]==0 and y.data[i]==0)
				TN+=1;
			else if (datatestclass[i]==1 and y.data[i]==0)
				FN+=1;
		}
	
		cout << TP << "\t" << FP << endl;
		cout << FN << "\t" << TN << endl;
		
		myfile << "Basiswidth: " << basisWidth << "\tNo Relevance Vectors: " << PARAMETERrev.size() << "\tTime: " << timeran <<endl;
		myfile << "Vals:\t" << TP << "\t" << FP << "\t" << FN << "\t" << TN << endl;



			datatest.clear();
			datatestclass.clear();
			s=test;
			readfile(s,datatest,datatestclass);

			if (kern==1)
				kernelfunction_gauss_test(BASIS, data, datatest, basisWidth, PARAMETERrev);	
			else if(kern==2)
				kernelfunction_cauch_test(BASIS, data, datatest, basisWidth, PARAMETERrev);	
			else if(kern==3)
				kernelfunction_bin_test(BASIS, data, datatest, basisWidth, PARAMETERrev);	
			
			matrixprod(BASIS, PARAMETERval, y, 0, 1.0);
	
			Sigmoid(y,yout);
		
	
			TP=0,FP=0,TN=0,FN=0;
	
			for (int i=0; i<datatestclass.size(); i++) {
				if (datatestclass[i]==1 and y.data[i]==1)
					TP+=1;
				else if (datatestclass[i]==0 and y.data[i]==1)
					FP+=1;
				else if (datatestclass[i]==0 and y.data[i]==0)
					TN+=1;
				else if (datatestclass[i]==1 and y.data[i]==0)
						FN+=1;
			}
			

	
		cout << TP << "\t" << FP << endl;
		cout << FN << "\t" << TN << endl;
		myfile << "Test:\t" << TP << "\t" << FP << "\t" << FN << "\t" << TN << endl;
		myfile.close();
	}
	else cout << "Unable to open file";
	
	string output3="model";
	ofstream myfile_model(output3.c_str(),ios::out);
	for(int i=0; i<PARAMETERrev.size(); i++){
		myfile_model << PARAMETERrev[i] << "\t";
	}
	myfile_model << endl;
	for(int i=0; i<PARAMETERval.rows; i++){
		myfile_model << PARAMETERval.data[i] << "\t";
	}
	myfile_model <<endl;
	
	myfile_model.close();

	return 0;
}
