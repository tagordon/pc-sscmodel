#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "blazar.h"

int main(int argc, char **argv){
	
//--------------------------------------------------------------------------------------------------------------------------------------------

	printf("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	printf("                BLAZARS!                 \n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
	printf("INPUT PARAMETERS:\n\n");
	
	// get input and set parameters for the simulation
	getArgs(argv[1]);
	getParams();
	allocateArrays();
	
	printf("\nBEGINNING COMPUTATION...\n\n");

	for(i=0;i<Nebins;i++){
		lE = lEbinmin+(i+0.5)*dlE;
		lNe_array[i] = lNe0(lE);
		lNe0_array[i] = lNe0(lE);
		lE_array[i] = lE;
		lA_array[i] = lA;
	}
	for(i=0;i<N;i++){
		x_array[i] = 1+i*dx;
	}
	for(i=0;i<Nv;i++){
		lv_sync_array[i] = lvmin_sync+i*dlv_sync;
	}
	for(i=0;i<Nv;i++){
		lP_sync_array[i] = (double *)malloc(N*sizeof(double));
		k_array[i] = (double *)malloc(N*sizeof(double));
	}
	
//--------------------------------------------------------------------------------------------------------------------------------------------
// walk down jet slice by slice computing synchrtron emissions and electron losses
	
	for(i=0;i<N;i++){
		x = x_array[i];
		lx = log(x);
		printf("\r                                 ");
		printf("\rJet Slice %d of %d",i+1,N);
		fflush(stdout);
		
		for(j=0;j<Nv;j++){
			lv = lv_sync_array[j];
			v = exp(lv);
			Ei = getIndex_Ee_log(lEe(log(x),lv));
			lAx = lA;
			if(Ei<Nebins){
				lAx = lA_array[Ei];
			}
			lP_sync_array[j][i] = lP(log(x),lv,lAx); 
			k_array[j][i] = exp(lk(log(x),lv,lAx));
		}
		for(k=0;k<Nebins;k++){
			double lE = lE_array[k];
			double lvp = get_lv(lx,lE);
			double losses_nf = losses(log(x),lvp,lA_array[k]);
			if(losses_nf < lNe_array[k]){
				lNe_array[k] = log(exp(lNe_array[k]) - exp(losses_nf));
			}
			else{
				lNe_array[k] = -100000;
			}
			lA_array[k] = lA + lNe_array[k] - lNe0_array[k];
		}
	}
//---------------------------------------------------------------------------------------------------------------------------------------------
// integrate synchrotron emissions along jet axis
	
	printf("\n");
	for(i=0;i<Nv;i++){
		printf("\r                                 ");
		printf("\rFrequency %d of %d",i+1,Nv);
		fflush(stdout);
		lv = lv_sync_array[i];
		P_obs = 0;
		for(j=0;j<N;j++){
			//printf("\nslice %d of frequency %d out of %d slices per frequency and %d frequencies.",j,i,N,Nv);
			x = x_array[j];
			tau_x = exp(ltau(lv,log(x),j,k_array[i]));
			P_obs = P_obs + exp(lP_sync_array[i][j]-tau_x);
		}
		lv = lvboost(lv);	
		// doppler boost emissions
		double lP_obs_boost = lboost(log(P_obs));
		double Ephot = lplanck10+lconvert*lv-le_charge10;
		double vF = lconvert*lv+(lconvert*lP_obs_boost)-lflux_factor_sync-3;
		fprintf(synchrotron,"%f\t%f\n",Ephot,vF);
	}

//--------------------------------------------------------------------------------------------------------------------------------------------
// print initial and final electron populations
		
	for(i=0;i<Nebins;i++){
		double noloss = lNe0_array[i];
		double loss = lNe_array[i];
		double E = lE_array[i];
		noloss = (lconvert*noloss);
		loss = (lconvert*loss);
		E = lconvert*(E -le_charge);
		fprintf(initPop,"%f\t%f\n",E,noloss);
		fprintf(population,"%f\t%f\n",E,loss);
	}
	
//--------------------------------------------------------------------------------------------------------------------------------------------
	
	printf("\n");
	fclose(synchrotron);
	fclose(population);
	fclose(initPop);
	return 0;
	
}