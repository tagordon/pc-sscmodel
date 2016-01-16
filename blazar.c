#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "blazar.h"

//--------------------------------------------------------------------------------------------------------------------------------------------
// thread function
	
	void *thread_ic(void *vargp){
		int *args = (int *)vargp;
		int start = args[0];
		int finish = args[1];
		int jt,kt,lt,mt;
		double phi2t,thetat,E_gammat,lEt;
		for(jt=start;jt<finish;jt++){
			ic_index++;
			phi2t = j*dphi2;
			for(kt=0;kt<Ntheta;kt++){
				thetat = kt*dtheta;
				for(lt=0;lt<NE_gamma;lt++){
					E_gamma = E_gamma_min + lt*dE_gamma;
					for(mt=0;mt<Nebins;mt++){
						lEt = lE_array[mt];
					}
				}
			}
		}
		return NULL;
	}
	
	void *thread_sync(void *vargp){
		int *args = (int *)vargp;
		int start = args[0];
		int finish = args[1];
		int it,jt;
		double lxt;
		double lvt;
		double tau_xt;
		for(it=start;it<finish;it++){
			sync_index++;
			printf("\nIntegrating Frequency %d of %d",sync_index,Nv);
			lvt = lv_sync_array[it];
			double P_obs = 0;
			for(jt=0;jt<N;jt++){
				//printf("\nslice %d of frequency %d out of %d slices per frequency and %d frequencies.",j,i,N,Nv);
				lxt = lx_array[jt];
				tau_xt = exp(ltau(lvt,lxt,jt,k_array[it]));
				P_obs = P_obs + exp(lP_sync_array[it][jt]-tau_xt);
			}
			lvt = lvboost(lvt);	
			// doppler boost emissions
			double lP_obs_boost = lboost(log(P_obs));
			double Ephot = lplanck10+lconvert*lvt-le_charge10;
			double vF = lconvert*lvt+(lconvert*lP_obs_boost)-lflux_factor_sync-3;
			fprintf(synchrotron,"%f\t%f\n",Ephot,vF);
		}
		
		return NULL;
	}

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
		lx_array[i] = log(1+i*dx);
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
		lx = lx_array[i];
		printf("\r                                 ");
		printf("\rJet Slice %d of %d",i+1,N);
		fflush(stdout);
		
		// synchrotron computation
		if(strcmp(do_sync,"yes")==0){
			for(j=0;j<Nv;j++){
				lv = lv_sync_array[j];
				Ei = getIndex_Ee_log(lEe(lx,lv));
				lAx = lA;
				if(Ei<Nebins){
					lAx = lA_array[Ei];
				}	
				lP_sync_array[j][i] = lP(lx,lv,lAx); 
				k_array[j][i] = exp(lk(lx,lv,lAx));
			}	
			for(k=0;k<Nebins;k++){
				double lE = lE_array[k];
				double lvp = get_lv(lx,lE);
				double losses_nf = losses(lx,lvp,lA_array[k]);
				if(losses_nf < lNe_array[k]){
					lNe_array[k] = log(exp(lNe_array[k]) - exp(losses_nf));
				}	
				else{
					lNe_array[k] = -100000;
				}	
				lA_array[k] = lA + lNe_array[k] - lNe0_array[k];
			}
		}
		
		// split ic computation among nthreads threads.
		if(strcmp(do_ic,"yes")==0){ 
			int thread_index;
			pthread_t *threads = (pthread_t *)malloc(nthreads*sizeof(pthread_t));
			int j_per_thread = Nphi2/nthreads;
			for(thread_index = 0;thread_index<nthreads;thread_index++){
				int *args = (int *)malloc(2*sizeof(int));
				args[0] = j_per_thread*thread_index;
				args[1] = args[0]+j_per_thread;
				pthread_create(&threads[thread_index],NULL,thread_ic,(void *)args);
			}
		
			for(thread_index = 0;thread_index<nthreads;thread_index++){
				pthread_join(threads[thread_index],NULL);
			}
			printf("\n");
		}	
	}

//---------------------------------------------------------------------------------------------------------------------------------------------
// integrate synchrotron emissions along jet axis
	if(strcmp(do_sync,"yes")==0){
		int thread_index_sync;
		pthread_t *threads_sync = (pthread_t *)malloc(nthreads*sizeof(pthread_t));
		int i_per_thread = Nv/nthreads;
		for(thread_index_sync = 0;thread_index_sync<nthreads;thread_index_sync++){
			int *args = (int *)malloc(2*sizeof(int));
			args[0] = i_per_thread*thread_index_sync;
			args[1] = args[0]+i_per_thread;
			pthread_create(&threads_sync[thread_index_sync],NULL,thread_sync,(void *)args);
		}
	
		for(thread_index_sync = 0;thread_index_sync<nthreads;thread_index_sync++){
			pthread_join(threads_sync[thread_index_sync],NULL);
		}
		printf("\n");
	}	
	
	// if(strcmp(do_sync,"yes")==0){
// 		printf("\n");
// 		for(i=0;i<Nv;i++){
// 			printf("\r                                 ");
// 			printf("\rFrequency %d of %d",i+1,Nv);
// 			fflush(stdout);
// 			lv = lv_sync_array[i];
// 			P_obs = 0;
// 			for(j=0;j<N;j++){
// 				//printf("\nslice %d of frequency %d out of %d slices per frequency and %d frequencies.",j,i,N,Nv);
// 				lx = lx_array[j];
// 				tau_x = exp(ltau(lv,lx,j,k_array[i]));
// 				P_obs = P_obs + exp(lP_sync_array[i][j]-tau_x);
// 			}
// 			lv = lvboost(lv);
// 			// doppler boost emissions
// 			double lP_obs_boost = lboost(log(P_obs));
// 			double Ephot = lplanck10+lconvert*lv-le_charge10;
// 			double vF = lconvert*lv+(lconvert*lP_obs_boost)-lflux_factor_sync-3;
// 			fprintf(synchrotron,"%f\t%f\n",Ephot,vF);
// 		}
// 	}

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