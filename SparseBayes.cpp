/*
 *  SparseBayes.cpp
 *  RVM-Speed
 *
 *  Created by Robert Lowe on 27/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "SparseBayes.h"

static mach_timebase_info_data_t	sTimebaseInfo;



double 
timeInMilliseconds(uint64_t time)
{    
	if ( sTimebaseInfo.denom == 0 ) {
        (void) mach_timebase_info(&sTimebaseInfo);
    }
	
    return (time * sTimebaseInfo.numer / sTimebaseInfo.denom)/1000000;
}


double SparseBayes(int LIKELIHOOD,int ItNum,int monitor_its,double MinDeltaLogAlpha,double AlignmentMax,std::vector<DOUBLE> data
				,std::vector<int> dataclass,int PriorityAddition,int PriorityDeletion,int BasisAlignmentTest,std::vector<int> &PARAMATERrev,matrix &PARAMATERval,int kernel,double basisWidth){
	

	//Calculate Basis
	double startTime = mach_absolute_time();

	matrix BASIS(data.size(),data.size());	
	matrix Targets(dataclass.size(),1);
	
	matrix BASISSCALES;

	cout << "Calculating Kernel...." <<endl;
	if (kernel==1)
		kernelfunction_gauss(BASIS, data,basisWidth,dataclass,Targets);
	else if(kernel==2)
		kernelfunction_cauch(BASIS, data,basisWidth,dataclass,Targets);
	else if(kernel==3)
		kernelfunction_binary(BASIS, data,basisWidth,dataclass,Targets);
	cout << "Done" << endl;
	cout << "Kernel time : " << timeInMilliseconds(mach_absolute_time() - startTime) <<endl;
	startTime = mach_absolute_time();
	
	matrix beta(dataclass.size(),1),Alpha,Mu,PHI,SIGMA,Factor;
	double logML;
	std::vector<int> Used;
	Initialisation(LIKELIHOOD,BASIS,dataclass,Alpha,beta,Mu,PHI,Used,BASISSCALES);
	matrix BASIS2(BASIS.rows,BASIS.cols);
	for (int i=0; i<BASIS.rows; i++) {
		for (int k=0; k<BASIS.cols; k++) {
			BASIS2.data[i*BASIS.cols+k]=(BASIS.data[i*BASIS.cols+k]*BASIS.data[i*BASIS.cols+k]);
		}
	}
	matrix BASIS_PHI,BASIS_B_PHI;
	

	matrix BASIS_Targets;
	matrixprod(BASIS,Targets,BASIS_Targets,1,1.0);

	matrix S_out,Q_in,S_in,Q_out,Gamma;
	
	fullstatistics(LIKELIHOOD, PHI, BASIS,BASIS2,beta,SIGMA,Mu,Alpha,logML,Targets,Used,Factor,S_out,Q_in,S_in,Q_out,BASIS_B_PHI,Gamma);


	
	int N=BASIS.rows;
	int M_full=BASIS.cols;
	int M=PHI.rows;
	
	int addCount=0;
	int deleteCount=0;
	int updateCount=0;
	
	//Control not present for betaupdatestart and BetaUpdateFrequency
	int maxLogSize=ItNum+10+(int)(ItNum/5);
	
	matrix logMarginalLog(maxLogSize,1);
	int count=0;
	
	std::vector<double> Aligned_out,Aligned_in;
	int alignDeferCount=0;
	
	//Action Codes
	const int ACTION_REESTIMATE=0;
	const int ACTION_ADD=1;
	const int ACTION_DELETE=-1;
	
	const int ACTION_TERMINATE=10;
	const int ACTION_NOISE_ONLY=11;
	const int ACTION_ALIGNMENT_SKIP=12;
	
	int selectedAction;
	/*
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 %%
	 %% MAIN LOOP
	 %%
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 */
	int iteration_counter=0;
	bool LAST_ITERATION=1;
	
	while (LAST_ITERATION) {
		
		iteration_counter+=1;
		//NEVER Update  iteration??
		bool UpdateIteration=0;
		/*
		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		 
		 %% DECISION PHASE
		 %%
		 %% Assess all potential actions
		 %%
		 */
		std::vector<double> GoodFactor(M_full);
		matrix DeltaML(M_full,1);
		matrix Action(M_full,1);

		for(int i=0; i<M_full; i++){
			DeltaML.data[i]=0;
			if (Factor.data[i]>1e-12) {
				GoodFactor[i]=1;
			}
			Action.data[i]*=ACTION_REESTIMATE;
		}
		

		matrix UsedFactor(Used.size(),1);
		std::vector<int> index_c,index_c_neg;
		std::vector<int> iu,iu_neg;

		for (int i=0; i<UsedFactor.rows; i++) {
			UsedFactor.data[i]=Factor.data[Used[i]];
			GoodFactor[Used[i]]=0;
			if (Factor.data[Used[i]]>1e-12){
				index_c.push_back(Used[i]);
				iu.push_back(i);
			}
			else {
				index_c_neg.push_back(Used[i]);
				iu_neg.push_back(i);
			}
			
		}
		if (BasisAlignmentTest){
			for (int i=0; i<Aligned_out.size(); i++) {
				GoodFactor[Aligned_out[i]]=0;
			}
		}
		
		matrix NewAlpha(index_c.size(),1);
		matrix Delta(index_c.size(),1);


		for(int i=0; i<index_c.size(); i++){
			NewAlpha.data[i]=(S_out.data[index_c[i]]*S_out.data[index_c[i]])/Factor.data[index_c[i]];
			Delta.data[i]=(1.0/NewAlpha.data[i])-(1.0/Alpha.data[iu[i]]);
			DeltaML.data[index_c[i]]= (Delta.data[i]*(Q_in.data[index_c[i]]*Q_in.data[index_c[i]])/(Delta.data[i]*S_in.data[index_c[i]]+1)-log(1+S_in.data[index_c[i]]*Delta.data[i]))/2.0;
		}
		

		//FREE BASIS OPTION NOT AVAILABLE
		bool anytoDelete=0;
		if(index_c_neg.size()!=0 and M>1){
			for(int i=0; i<index_c_neg.size(); i++){
				DeltaML.data[index_c_neg[i]]= -(Q_out.data[index_c_neg[i]]*Q_out.data[index_c_neg[i]]/(S_out.data[index_c_neg[i]]+Alpha.data[iu_neg[i]])
											 -log(1+S_out.data[index_c_neg[i]]/Alpha.data[iu_neg[i]]))/2.0;
				Action.data[index_c_neg[i]]=ACTION_DELETE;
				anytoDelete=1;
			}
		}
		

		//ANYTHING TO ADD
		index_c.clear();
		for(int i=0; i<GoodFactor.size(); i++){
			if(GoodFactor[i]==1){
				index_c.push_back(i);
			}
		}
		bool anytoADD=0;
		if(index_c.size()!=0){
			matrix quot(index_c.size(),1);
			for(int i=0; i<index_c.size(); i++){
				quot.data[i]=Q_in.data[index_c[i]]*Q_in.data[index_c[i]]/S_in.data[index_c[i]];
				DeltaML.data[index_c[i]]=(quot.data[i]-1-log(quot.data[i]))/2.0;
				Action.data[index_c[i]]=ACTION_ADD;
				anytoADD=1;
			}
		}
		

		//NOT TESTED? 
		if ((anytoADD && PriorityAddition) || (anytoDelete && PriorityDeletion)){
			//We won't perform re-estimation this iteration, which we achieve by
			//zero-ing out the delta
			for(int i=0; i<Action.rows; i++){
				if (Action.data[i]==ACTION_REESTIMATE)
					DeltaML.data[i]	= 0;
				//Furthermore, we should enforce ADD if preferred and DELETE is not
				// - and vice-versa
				if (anytoADD && PriorityAddition && PriorityDeletion){
					if (Action.data[i]==ACTION_DELETE)
						DeltaML.data[i]	= 0;
					
				}
				if (anytoDelete && PriorityDeletion && !PriorityAddition){
					if (Action.data[i]==ACTION_ADD)
						DeltaML.data[i]	= 0;	
					
				}
			}
		}
		
		//Choose the one with largest likelihood
		double deltaLogMallrginal=0.0;
		int nu=0;
		selectedAction=0;
		bool anyWorthwhileAction;
		for (int i=0; i<DeltaML.rows; i++) {
			if (DeltaML.data[i]>deltaLogMallrginal){
				//Updated 23/03/10 so that we do not delete when only one relevance vector
				if(Action.data[i]==-1 && Mu.rows==1){
					cout << "TRYING TO DELETE WHEN MU=1" << endl;
				}
				else{
					deltaLogMallrginal=DeltaML.data[i];
					nu=i;
					selectedAction=Action.data[i];
				}
			}
		}
		
		int j=0;
		anyWorthwhileAction=deltaLogMallrginal>0;
		if (selectedAction==ACTION_REESTIMATE || selectedAction==ACTION_DELETE){
			for	(int i=0; i<Used.size(); i++){
				if (Used[i]==nu)
					j=i;
			}
		}
		
		//DIFFERENCE TO MATLAB WITH NU being selected./RVM-Speed -b 0.99 -k Binary -d kbd420
		//MATLAB 72	0.7476216971706305	0.7476216971706305
		//C++ 
		/*
		if (iteration_counter>220) {
			printf("%d\t%.16f\t%.16f\n",nu,deltaLogMallrginal ,DeltaML.data[71]);
 		}
		*/
		
		matrix Phi;
		string act;
		Phi.AddColumn(BASIS, nu);
		
		double newAlpha=S_out.data[nu]*S_out.data[nu]/Factor.data[nu];

		if (!anyWorthwhileAction || (selectedAction==ACTION_REESTIMATE && abs(log(newAlpha)-log(Alpha.data[j]))<MinDeltaLogAlpha && !anytoDelete)){
			selectedAction=ACTION_TERMINATE;
			act="potential termination";
		}

		if (BasisAlignmentTest){
			if (selectedAction==ACTION_ADD){
				matrix p;
				matrixprod(Phi,PHI,p,1,1.0);
				if (p.rows>1){
					cout << "BELIEVED ERROR: check p, should be a single row" << endl;
					exit(1);
				}
				std::vector<int> findAligned;
				for (int i=0; i<p.cols; i++){
					if(p.data[i]>AlignmentMax){
						findAligned.push_back(i);
					}
				}
				int numAligned=findAligned.size();
				if (numAligned>0){
					selectedAction=ACTION_ALIGNMENT_SKIP;
					act="alignment-deferred addition";
					alignDeferCount+=1;
					for(int i=0; i<numAligned; i++){
						Aligned_out.push_back(nu);
						Aligned_in.push_back(Used[findAligned[i]]);
					}
				}
			}
			if (selectedAction==ACTION_DELETE){
				std::vector<int> findAligned;
				for(int i=0; i<Aligned_in.size(); i++){
					if(Aligned_in[i]==nu){
						findAligned.push_back(i);
					}
				}
				int numAligned=findAligned.size();
				//COuld include DIAGNOSTICS here
				if(numAligned>0){
					for (int i=(numAligned-1); i>-1; i--) {
						Aligned_in.erase(Aligned_in.begin()+findAligned[i]);
						Aligned_out.erase(Aligned_out.begin()+findAligned[i]);
					}
				}
			}
			
		}


		/*		
		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		 
		 %% ACTION PHASE
		 %%
		 %% Implement above decision
		 %%
		 */
		bool UPDATE_REQUIRED=0;
		matrix SIGMANEW;
		switch (selectedAction) {
			case ACTION_REESTIMATE:{
				
				double oldAlpha=Alpha.data[j];
				Alpha.data[j]=newAlpha;
				matrix s_j;
				
				s_j.AddColumn(SIGMA, j);
				double deltaInv=1.0/(newAlpha-oldAlpha);
				double kappa=1.0/(SIGMA.data[j*SIGMA.cols+j]+deltaInv);

				matrix tmp(s_j.rows,s_j.cols);
				matrix deltaMu(tmp.rows,tmp.cols);

				for (int i=0; i<s_j.rows; i++) {
					for (int k=0; k<s_j.cols; k++) {
						tmp.data[i*s_j.cols+k]=s_j.data[i*s_j.cols+k]*kappa;
						deltaMu.data[i*s_j.cols+k]=tmp.data[i*s_j.cols+k]*-Mu.data[j];
						Mu.data[i*s_j.cols+k]+=deltaMu.data[i*s_j.cols+k];
					}
				}
				matrix SIGMANEW;
				matrixprod(tmp, s_j, SIGMANEW, 2, 1.0);
				
				for (int i=0; i<SIGMA.rows; i++) {
					for (int k=0; k<SIGMA.cols; k++) {
						SIGMANEW.data[i*SIGMANEW.cols+k]=SIGMA.data[i*SIGMA.cols+k]-SIGMANEW.data[i*SIGMANEW.cols+k];
					}
				}
				
				if (UpdateIteration){
					matrix tempbbpsj;
					matrixprod(BASIS_B_PHI,s_j,tempbbpsj,0,1.0);
					for (int i=0; i<tempbbpsj.rows; i++){
						for(int k=0; k<tempbbpsj.cols; k++){
							tempbbpsj.data[i*tempbbpsj.cols+k]=(tempbbpsj.data[i*tempbbpsj.cols+k]*tempbbpsj.data[i*tempbbpsj.cols+k]);
						}
					}
					for (int i=0; i<S_in.rows; i++) {
						S_in.data[i]=S_in.data[i]+kappa*(tempbbpsj.data[i]);

					}
					//printoutMatrix(S_in);
					matrix tmpq;
					matrixprod(BASIS_B_PHI, deltaMu, tmpq, 0, 1.0);
					for (int i=0; i<Q_in.rows; i++) {
						Q_in.data[i]=Q_in.data[i]-tmpq.data[i];
						
					}
				}
				updateCount+=1;
				act="re-estimation";
				UPDATE_REQUIRED=1;
			}
				break;
			case ACTION_ADD:{
				matrix B_Phi(beta.rows,beta.cols);
				for(int i=0; i<B_Phi.rows; i++){
					for (int k=0; k<B_Phi.cols; k++) {
						B_Phi.data[i*B_Phi.cols+k]=Phi.data[i*B_Phi.cols+k]*beta.data[i*B_Phi.cols+k];
					}
				}
				
				matrix BASIS_B_phi;

				matrixprod(BASIS,B_Phi,BASIS_B_phi,1,1.0);
			
				matrix tmp0;
				matrix tmp;
				matrixprod(B_Phi, PHI, tmp0, 1, 1.0);
				matrixprod(tmp0, SIGMA, tmp, 0, 1.0);
				tmp.rows=tmp.cols;
				tmp.cols=1;

				//cout << "tmp "<<endl;
				double *tmpalphadata=new double[Alpha.rows*Alpha.cols];
				
				
				Alpha.resize(Alpha.rows+1, Alpha.cols);
				Alpha.data[Alpha.rows-1]=newAlpha;
				PHI.AddColumn(Phi, 0);

				double s_ii=1/(newAlpha+S_in.data[nu]);
				matrix s_i(tmp.rows,1);
				for (int i=0; i<tmp.rows; i++) {
					s_i.data[i]=-s_ii*tmp.data[i];
				}
				matrix TAU;
				matrixprod(s_i,tmp,TAU,2,-1.0);
				SIGMANEW.resize(s_i.rows+1, s_i.rows+1);

				for(int i=0; i<SIGMANEW.rows; i++){
					for(int k=0; k<SIGMANEW.cols; k++){
						if(i<SIGMA.rows and k<SIGMA.cols){
							SIGMANEW.data[i*SIGMANEW.cols+k]=SIGMA.data[i*SIGMA.cols+k]+TAU.data[i*TAU.cols+k];
						}
						else if (i==SIGMA.rows and k<s_i.rows){
							SIGMANEW.data[i*SIGMANEW.cols+k]=s_i.data[k];
						}
						else if(k==SIGMA.rows and i<s_i.rows){
							SIGMANEW.data[i*SIGMANEW.cols+k]=s_i.data[i];
						}
						else if(i==SIGMA.rows and k==SIGMA.rows){
							SIGMANEW.data[i*SIGMANEW.cols+k]=s_ii;
						}
						else {
							cout << "ERROR IN ACTION ADD STATEMENT" << endl;
							exit(1);
						}
					}
				}
				double mu_i=s_ii*Q_in.data[nu];
				matrix deltaMu(tmp.rows+1,1);
				for (int i=0; i<deltaMu.rows-1; i++) {
					deltaMu.data[i]=-mu_i*tmp.data[i];
				}
				deltaMu.data[deltaMu.rows-1]=mu_i;
				//cout << "MU is being updated" <<endl;
				//printoutMatrix(deltaMu);
				Mu.resize(Mu.rows+1, Mu.cols);
				Mu.data[Mu.rows-1]=0;
				for (int i=0; i<Mu.rows; i++) {
					Mu.data[i]+=deltaMu.data[i];

				}
				//printoutMatrix(Mu);
				if(UpdateIteration){
					matrix mctmp;
					matrixprod(BASIS_B_PHI,tmp,mctmp,0,1.0);
					matrix mCi(mctmp.rows,mctmp.cols);
					matrix mCi2=mCi;

					for (int i=0; i<mCi.rows; i++) {
						for (int k=0; k<mCi.cols; k++) {
							mCi.data[i*mCi.rows+k]=BASIS_B_phi.data[i*mCi.rows+k]-mCi.data[i*mCi.rows+k];
							mCi2.data[i*mCi.rows+k]=mCi.data[i*mCi.rows+k]*mCi.data[i*mCi.rows+k];
						}
					}
					for (int i=0; i<S_in.rows; i++) {
						S_in.data[i]=S_in.data[i]-s_ii*mCi2.data[i];
						
					}
					//printoutMatrix(S_in);
					for (int i=0; i<Q_in.rows; i++) {
						Q_in.data[i]=Q_in.data[i]-mu_i*mCi.data[i];
						
					}

				}
				Used.push_back(nu);
				addCount+=1;
				act="addition";
				UPDATE_REQUIRED=true;
			}
				break;
			case ACTION_DELETE:{
				PHI.RemoveColumn(j);
				Alpha.RemoveRow(j);
				double s_jj=SIGMA.data[j*SIGMA.cols+j];
				matrix s_j;
				s_j.AddColumn(SIGMA, j);

				matrix tmp(s_j.rows,s_j.cols);
				for(int i=0; i<s_j.rows; i++){
						tmp.data[i]=s_j.data[i]/s_jj;
				}
				matrix SIGMANEW;

				matrixprod(tmp, s_j, SIGMANEW, 2, 1.0);
				
				for(int i=0; i<SIGMANEW.rows; i++){
					for (int k=0; k<SIGMANEW.cols; k++) {
						SIGMANEW.data[i*SIGMANEW.cols+k]=SIGMA.data[i*SIGMA.cols+k]-SIGMANEW.data[i*SIGMANEW.cols+k];
					}
				}
				SIGMANEW.RemoveRow(j);
				SIGMANEW.RemoveColumn(j);
				
				matrix deltaMu(tmp.rows,1);
				for (int i=0; i<tmp.rows; i++) {
					deltaMu.data[i]=-Mu.data[j]*tmp.data[i];
					Mu.data[i]+=deltaMu.data[i];

				}
				double mu_j=Mu.data[j];
				Mu.RemoveRow(j);
				if (UpdateIteration){
					matrix jPm;
					matrixprod(BASIS_B_PHI,s_j,jPm,0,1.0);
					matrix jPm2(jPm.rows,jPm.cols);
					for (int i=0; i<jPm.rows; i++) {
						for (int k=0; k<jPm.cols; k++) {
							jPm2.data[i*jPm.cols+k]=jPm.data[i*jPm.cols+k]*jPm.data[i*jPm.cols+k];
						}
					}
					for (int i=0; i<S_in.rows; i++) {
						S_in.data[i]=S_in.data[i]+jPm2.data[i]/s_jj;
						
					}
					//printoutMatrix(S_in);
					for (int i=0; i<Q_in.rows; i++) {
						Q_in.data[i]=Q_in.data[i]+jPm.data[i]*(mu_j/s_jj);
						
					}

				}
				Used.erase(Used.begin()+(j));
				deleteCount+=1;
				act="deletion";
				UPDATE_REQUIRED=true;
			}
				break;
				
			default:
				break;
		}
		M=Used.size();
		
		//cout << "ACTION: " << act << " of " << nu << " ("<<deltaLogMallrginal<<")"<<endl;
		

		/*
		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		 
		 %% UPDATE STATISTICS
		 
		 % If we've performed a meaningful action,
		 % update the relevant variables
		 % 
		 */
		if (UPDATE_REQUIRED){
			if(UpdateIteration){
				cout << "UPDATE ITERATION NOT IMPLEMENTED " << endl;
				/*
				S_out=S_in;
				Q_out=Q_in;
				matrix tmp(Used.size(),1);
				for (int i=0; i<Used.size(); i++) {
					tmp.data[i]=Alpha.data[i]/(Alpha.data[i]-S_in.data[Used[i]]);
					S_out.data[Used[i]]=tmp.data[i]*S_in.data[Used[i]];
					Q_out.data[Used[i]]=tmp.data[i]*Q_in.data[Used[i]];
				}
				for (int i=0; i<Q_out.size1(); i++) {
					Factor.data[i]=(Q_out.data[i]*Q_out.data[i])-S_out.data[i];
				}
				SIGMA=SIGMANEW;
				for (int i=0; i<Alpha.size1(); i++) {
					Gamma.data[i]=1-Alpha.data[i]*SIGMA.data[i*SIGMA.cols+i];
				}
				matrix temp3(beta.size1(),M);		
				for(int i=0; i<beta.size1(); i++){
					for (int k=0; k<M; k++){
						temp3.data[i*temp3.cols+k]=PHI.data[i*PHI+k]*beta.data[i];
					}
				}
				BASIS_B_PHI=trans(matrixprod(temp3,BASIS,1));
				*/
			}
			else{
				double newLogML=0.0;

				fullstatistics(LIKELIHOOD, PHI, BASIS,BASIS2,beta,SIGMA,Mu,Alpha,newLogML,Targets,Used,Factor,S_out,Q_in,S_in,Q_out,BASIS_B_PHI,Gamma);

				deltaLogMallrginal=newLogML-logML;
			}
			if(UpdateIteration && deltaLogMallrginal<0){
				cout << "** Alert **  DECREASE IN LIKELIHOOD !!" << endl;
			}
			logML=logML+deltaLogMallrginal;
			count+=1;
			logMarginalLog.data[count]=logML;
		}
		
		/*
		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		 
		 % END OF CYCLE PROCESSING
		 % 
		 % Check if termination still specified, and output diagnostics
		 % 
		 */
		double sumgamma=0.0;
		for (int i=0; i<Gamma.rows; i++) {
			sumgamma+=Gamma.data[i];
		}
		
		double sqrtbeta=0.0;
		for (int i=0; i<beta.rows; i++) {
			sqrtbeta+=sqrt(1.0/beta.data[i]);
		}		
		if(selectedAction==ACTION_TERMINATE){
			cout << "** Stopping at iteration "<< iteration_counter << " (Max_delta_ml="<<deltaLogMallrginal<<") **"<<endl;
			printf("'%4d> L = %.6f\t Gamma = %.2f (M = %d)\n",iteration_counter, logML/N,sumgamma, M);
			break;
		}
		
		//NO TIME LIMIT AS OF YET!!!!!!
		

		if ((iteration_counter%monitor_its==0 || iteration_counter==1)){
			printf("%5d> L = %.6f\t Gamma = %.2f (M = %d)\n",iteration_counter, logML/N, sumgamma, M);
		}
		
		if(iteration_counter==ItNum){
			break;
		}
		
	}
	/*
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%
	%% POST-PROCESSING
	%%
	
	%
	% Warn if we exited the main loop without terminating automatically
		% 
		*/
	if (selectedAction!=ACTION_TERMINATE){
		printf("** Iteration limit: algorithm did not converge\n");
	}
	int total	= addCount + deleteCount + updateCount;
	if(BasisAlignmentTest)
		total	= total+alignDeferCount;

	printf("Action Summary\n");
	printf("==============\n");
	printf("Added\t\t\t%6d (%.0f%%)\n",addCount, 100*addCount/(total*1.0));
	printf("Deleted\t\t%6d (%.0f%%)\n",deleteCount, 100*deleteCount/(total*1.0));
	printf("Reestimated\t%6d (%.0f%%)\n",updateCount, 100*updateCount/(total*1.0));
	if (BasisAlignmentTest && alignDeferCount){
		printf("--------------\n");
		printf("Deferred\t%6d (%.0f%%)\n",alignDeferCount, 100*alignDeferCount/(total*1.0));
	}
	printf("==============\n");
	printf("Total of %d likelihood updates\n", count);
	cout << "Time to run: " <<  timeInMilliseconds(mach_absolute_time() - startTime)/1000.0 << "s" << endl;

	PARAMATERrev=Used;

	vector<size_t> index;
	for (unsigned i = 0; i < PARAMATERrev.size(); ++i)
		index.push_back(i);
	sort(index.begin(), index.end(), index_cmp<vector<int>&>(PARAMATERrev));
	sort (PARAMATERrev.begin(), PARAMATERrev.end());
	
	PARAMATERval.reset(index.size(),1);
	for (int i=0; i<index.size(); i++) {
		PARAMATERval.data[i]=(Mu.data[index[i]]/(BASISSCALES.data[Used[index[i]]]));
	}

	return(timeInMilliseconds(mach_absolute_time() - startTime)/1000.0);

}





void kernelfunction_gauss(matrix &BASIS,const std::vector<DOUBLE> &data,double basisWidth,const std::vector<int> &dataclass,matrix &Targets){
	
	//Calculate Basis
	matrix X2(data.size(),data.size());
	matrix Y2(data.size(),data.size());

	matrix X(data.size(),data[0].size());

	
	for (int i=0; i<data.size(); i++) {
		double sumofsquares=0.0;
		for(int k=0; k<data[i].size(); k++){
			sumofsquares+=pow(data[i][k],2);
			X.data[(i*data[i].size())+k]=data[i][k];
		}
		Targets.data[i]=dataclass[i];
		for(int p=0; p<X2.cols; p++){
			X2.data[(i*X2.cols)+p]=sumofsquares;
			Y2.data[(p*X2.cols)+i]=sumofsquares;
		}
	}
	
	matrixprod(X, X,BASIS, 2,-2.0);
	
	
	for (int i=0; i<BASIS.rows; i++) {
		for (int k=0; k<BASIS.cols; k++) {
			BASIS.data[(i*BASIS.cols)+k]=exp((X2.data[(i*BASIS.cols)+k]+Y2.data[(i*BASIS.cols)+k]+BASIS.data[(i*BASIS.cols)+k])*-basisWidth);
		}
	}
	

}

void kernelfunction_cauch(matrix &BASIS,const std::vector<DOUBLE> &data,double basisWidth,const std::vector<int> &dataclass,matrix &Targets){
	//Calculate Basis
	matrix X2(data.size(),data.size());
	matrix Y2(data.size(),data.size());
	
	matrix X(data.size(),data[0].size());
	
	
	for (int i=0; i<data.size(); i++) {
		double sumofsquares=0.0;
		for(int k=0; k<data[i].size(); k++){
			sumofsquares+=pow(data[i][k],2);
			X.data[(i*data[i].size())+k]=data[i][k];
		}
		Targets.data[i]=dataclass[i];
		for(int p=0; p<X2.cols; p++){
			X2.data[(i*X2.cols)+p]=sumofsquares;
			Y2.data[(p*X2.cols)+i]=sumofsquares;
		}
	}
	
	matrixprod(X, X,BASIS, 2,-2.0);
	
	
	for (int i=0; i<BASIS.rows; i++) {
		for (int k=0; k<BASIS.cols; k++) {
			BASIS.data[(i*BASIS.cols)+k]=1.0/(1.0+((X2.data[(i*BASIS.cols)+k]+Y2.data[(i*BASIS.cols)+k]+BASIS.data[(i*BASIS.cols)+k])*basisWidth));
		}
	}
	
	
}

void kernelfunction_binary(matrix &BASIS,const std::vector<DOUBLE> &data,double basisWidth,const std::vector<int> &dataclass,matrix &Targets){

	BASIS.reset(data.size(), data.size());
	for (int i=0; i<data.size(); i++) {
		Targets.data[i]=dataclass[i];
		for (int j=i; j<data.size(); j++) {
			double sum=0.0;
			for (int k=0; k<data[i].size(); k++) {
				sum+=abs(data[i][k]-data[j][k]);
			}
			BASIS.data[i*BASIS.cols+j]=pow(basisWidth,(data[0].size()-sum))*pow(1-basisWidth,sum);
			if (j!=i)
				BASIS.data[j*BASIS.cols+i]=pow(basisWidth,(data[0].size()-sum))*pow(1-basisWidth,sum);
		}
	}
}

		 

void Initialisation(int LIKELIHOOD, matrix &BASIS, std::vector<int> Targets,matrix &Alpha,matrix &beta,matrix &Mu, matrix &Phi, std::vector<int> &Used, matrix &Scales){
	//Done for initialisation using Bernoulli 23/08/2010
	double GAUSSIAN_SNR_INIT	= 0.1;
	double	INIT_ALPHA_MAX	= 1e3;
	double INIT_ALPHA_MIN	= 1e-3;
	
	//Normalise each Basis vector
	Scales.reset(1,BASIS.cols);
	
	for (int k=0; k<BASIS.rows; k++){
		Scales.data[k]=0;
		for (int i=0; i<BASIS.cols; i++){
			Scales.data[k]+=pow(BASIS.data[(i*BASIS.cols)+k],2);
		}
		Scales.data[k]=sqrt(Scales.data[k]);
		if (Scales.data[k]==0){
			Scales.data[k]=1;
		}
	}
	
	for (int i=0; i<BASIS.rows; i++){
		for (int k=0; k<BASIS.cols; k++){
			BASIS.data[(i*BASIS.cols)+k]=BASIS.data[(i*BASIS.cols)+k]/Scales.data[k];
		}
	}
	
	matrix LogOut(Targets.size(),1);
	matrix TargetsPseudoLinear(Targets.size(),1);

	for (int i=0; i<Targets.size(); i++){
		TargetsPseudoLinear.data[i]=(2*Targets[i]-1);
		LogOut.data[i]	= (TargetsPseudoLinear.data[i]*0.9+1)/2.0;
		LogOut.data[i]=log(LogOut.data[i]/(1-LogOut.data[i]));
	}
	
	matrix proj;
	matrixvprod(BASIS,TargetsPseudoLinear,proj,1,1.0);

	double max=0.0;
	int maxindex=0;
	for(int i=0; i<proj.rows; i++){
		if (fabs(proj.data[i])>max) {
			max=fabs(proj.data[i]);
			maxindex=i;
		}
	}
	
	cout << "Initialising with maximally aligned basis vector (" << maxindex+1 << ")" << endl;
	Used.push_back(maxindex);
	Phi.AddColumn(BASIS, maxindex);
	
	linalg(Phi, LogOut,Mu,0);
	
	Alpha.data=new double[1];
	Alpha.rows=1;
	Alpha.cols=1;
	
	cout << "Initial Alpha = ";
	if (Mu.data[0]==0)
		Alpha.data[0]=1;
	else{
		double value=1/(Mu.data[0]*Mu.data[0]);
		if (value<INIT_ALPHA_MIN) 
			Alpha.data[0]=INIT_ALPHA_MIN;
		else if(value>INIT_ALPHA_MAX)
			Alpha.data[0]=INIT_ALPHA_MAX;
		else
			Alpha.data[0]=value;
	}
	cout << Alpha.data[0] << endl;

}

										  
		
