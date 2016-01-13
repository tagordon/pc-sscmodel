#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "blazar.h"

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// synchrotron functions

// returns the log of the radius of the jet at position x
double lR(double lx){
	double x = exp(lx);
	return log(exp(lR0)+x*tan(theta_opening));
}           

// returns the magnetic field in the jet at position x                
double lB(double lx){
	return lB0+lR0-lR(lx);
}              

// the quantity epsilon for finding the electron energy that radiates at v in slice x                      
double leps(double lx){
	return log(4.0)+lpi+(3.0*lme)+(4.0*lc)-(log(3.0)+le_charge+lB(lx));
}   

// returns the log of the frequency of synchrotron emission associated with a given energy and x-slice
double getlv(double lx, double lE){
	return 2*lE - leps(lx);
}

// electron energy that emits synchrotron radiation at the critical frequency                       
double lEe(double lx, double lv){
	return 0.5*(lv+leps(lx));
}

double get_lv(double lx, double lEe){
	double lR_nf = log(exp(lR0)+x*tan(theta_opening));
	double lB_nf = lB0+lR0-lR_nf;
	double leps_nf = log(4.0)+lpi+(3.0*lme)+(4.0*lc)-(log(3.0)+le_charge+lB_nf);
	return 2*lEe-leps_nf;
}

// for calculating opacity - remember that you might change the last term of this...
double lj(double lx, double lv, double lAx){
	double lE_electron = 0.5*(leps(lx)+lv);
	gamma2 = exp(2*lE_electron-2*lme-4*lc);
	lbeta = 0.5*log(1-(1/gamma2));
	constants = log(2.0)+lsigma_t-log(3.0)-lpi-lmu0-(2.0*lme)-(3.0*lc);
	return constants+(2.0*lbeta)+(2.0*lB(lx))+leps(lx)+lAx+((1.0-alpha)*lEe(lx,lv))-(2.0*lR(lx))-(exp(lEe(lx,lv)-lEmax)); 
}

// opacity for photons of frequency v at position x in the jet                                            
double lk(double lx, double lv, double lAx){
	double x = exp(lx);
	double lR_nf = log(exp(lR0)+x*tan(theta_opening));
	double lB_nf = lB0+lR0-lR_nf;
	double leps_nf = log(4.0)+lpi+(3.0*lme)+(4.0*lc)-(log(3.0)+le_charge+lB_nf);
	double lE_electron = 0.5*(leps_nf+lv);
	gamma2 = exp(2*lE_electron-2*lme-4*lc);
	lbeta = 0.5*log(1-(1/gamma2));
	constants = log(2.0)+lsigma_t-log(3.0)-lpi-lmu0-(2.0*lme)-(3.0*lc);
	double lj0_nf =  constants+lAx+leps_nf+((1-alpha)/2.0)*leps_nf+2*lB_nf+2*lbeta-2*lR_nf-exp(lE_electron-lEmax);
	double lk_nf = lj0_nf+2*lc-((alpha+4)/2.0)*lv-log(2)-(0.5*leps_nf);
	return lk_nf;  
}

// optical depth of slice x for frequency v - i is the ith division in x and k_array is an array holding all the opacities for the frequency v.                                    
double ltau(double lv, double lx, int i, double *k_array){
	double beta = sqrt(1-(1/(gamma_bulk*gamma_bulk)));
	return (2*lgamma_bulk)+log((1.0/cos(theta_observe))-beta)+log(simple_integral_simpson(k_array, i, N-1, dx));  
}

// power lost at slice x by frequency v                             
double lP(double lx, double lv, double lAx){
	double x = exp(lx);
	double lR_nf = log(exp(lR0)+x*tan(theta_opening));
	double lB_nf = lB0+lR0-lR_nf;
	double leps_nf = log(4.0)+lpi+(3.0*lme)+(4.0*lc)-(log(3.0)+le_charge+lB_nf);
	double lE_electron = 0.5*(leps_nf+lv);
	gamma2 = exp(2*lE_electron-2*lme-4*lc);
	lbeta = 0.5*log(1-(1/gamma2));
	constants = log(2.0)+lsigma_t-log(3.0)-lpi-lmu0-(2.0*lme)-(3.0*lc);
	double lj0_nf =  constants+lAx+leps_nf+((1-alpha)/2.0)*leps_nf+2*lB_nf+2*lbeta-2*lR_nf-exp(lE_electron-lEmax);
	double lk_nf = lj0_nf+2*lc-((alpha+4)/2.0)*lv-log(2)-(0.5*leps_nf);
	double factor2 = log(1-exp(-1*exp(lk_nf+lR_nf)));
	if(lk_nf+lR_nf<-10){
		return lpi+(lR_nf)+ldx-(2.0*lc)+log(2.0)+(0.5*leps_nf)+(2.5*lv)+lk_nf+lR_nf;
	}
	return lpi+(lR_nf)+ldx-(2.0*lc)+log(2.0)+(0.5*leps_nf)+(2.5*lv)+log(1-exp(-1*exp(lk_nf+lR_nf)));
}

// log of initial (x=0) electron population
double lNe0(double lE){
	return lA-(alpha*lE)-exp(lE-lEmax);
}

// dopplar boost emissions from fluid frame to observer frame
double lboost(double lP){
	double beta = sqrt(1-(1/(gamma_bulk*gamma_bulk)));
	return lP-4*log(gamma_bulk*(1-beta*cos(theta_observe)));
}

double lvboost(double lv){
	double beta = sqrt(1-(1/pow(gamma_bulk,2)));
	return lv-log(gamma_bulk*(1-beta*cos(theta_observe)));
}

// integrates an array using the midpoint rule                                           
double simple_integral(double *f,int i, int n, double dx){
	double sum=0;
	        while(i<n){
	                i++;
	                sum = sum+dx*(f[i-1]+f[i])/2.0;
	        }
	return sum;
}  

double simple_integral_simpson(double *f,int i, int n, double dx){
	double sum=0;
	        while(i<n-2){
	                i=i+2;
	                sum = sum+f[i-2]+4*f[i-1]+f[i];
	        }
	return (dx/3.0)*sum;
}  

double simple_integral_comp_simpson(double *f,int i, int n, double dx){
	        float x;
	        float r;
	        char m = 0;
	        float s = 0.0;

	        for (i = 0; i <= n; i++) {
	                r = f[i];
	                if (i == 0 || i == N) {
	                        s += r;
	                } else {
	                        m = !m;
	                        s += r * (m+1) * 2.0;
	                }
	        }
	        return s * (dx/3.0);
}  

// returns the index of the electron energy array for a given energy
int getIndex_Ee(double lE){
	return floor((exp(lE)-exp(lEbinmin))/dE);
}

int getIndex_Ee_log(double lE){
	return floor((lE-lEbinmin)/dlE);
}

// returns the log of the energy for a given energy index
double getlE_Ee(int Ei){
	return lE_array[Ei];
}

double losses(double lx,double lv, double lAx){
	double x = exp(lx);
	double lR_nf = log(exp(lR0)+x*tan(theta_opening));
	double lB_nf = lB0+lR0-lR_nf;
	double leps_nf = log(4.0)+lpi+(3.0*lme)+(4.0*lc)-(log(3.0)+le_charge+lB_nf);
	double lE = 0.5*(leps_nf+lv);
	
	return lP(lx,lv,lAx)+log(2)-lc-leps_nf;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// inverse compton functions

// returns the index of the scattered photon energy for a given energy
int getIndex_scatt(double lE){
	return floor((lE-log(E_scatt_min))/dE_scatt);
}

double logdE(int i){
	return log(exp(lEbinmin+(Ei+1)*dlE)-exp(lEbinmin+(Ei)*dlE));
}

// returns the scattered photon energy for a given index (not log)
double getlE_scatt(int Ei){
	return E_gamma_scatt_array_E[Ei];
}

double get_ldE_scatt(int Ei){
	return log(exp(log(E_scatt_min)+(i+1)*dE_scatt)-exp(log(E_scatt_min)+i*dE_scatt));
}

// returns the energy of a ic scattered photon 
double lE_scatt(double lE, double lEe, double theta, double phi2){	
	double gamma2 = exp(2*lEe-2*lme-4*lc);
	double beta = sqrt(1-1/(gamma2));
	
	double phi = acos((cos(theta)+beta)/(1+beta*cos(theta)));
	double cos_psi_1 = (beta-cos(phi2+phi))/(1-(beta*cos(phi2+phi)));
	
	double factor = (1+beta*cos(theta))/(1-(beta*cos_psi_1));
	return lE+log(factor);
}

// used to find the energy of an ic scattered photon (the angle through which the photon is scattered)
double phi_scatt(double lE,double lEe, double theta, double phi2){
	double gamma2 = exp(2*lEe-2*lme-4*lc);
	double beta = sqrt(1-1/(gamma2));
	double phi = acos((cos(theta)+beta)/(1+beta*cos(theta)));
	return acos((cos(phi2+phi)-beta)/(1-beta*cos(phi2+phi)));
}

// returns the photon population in the optically thin case for ic calculation
double ln_gamma_thin(double lx, double lv){
	double x = exp(lx);
	double lR_nf = log(exp(lR0)+x*tan(theta_opening));
	
	double lP_for_n_gamma = lP(lx,lv,lA);
	return lP_for_n_gamma-log(2*pi)-lplanck-lv-lR_nf-ldx;
}

// returns the photon population in the optically thick case for ic calculation
double ln_gamma_thick(double lx, double lv){
	double factor = exp(lplanck + lv - lEe(lx,lv));
	return log(8*pi) + 2*lv - 3*lc - log(exp(factor)-1);
}

// this is actually the function P in potter & cotter, but P is power other places, used for ic calculation
double Q(double lE_gamma, double phi){
	return 1.0/(1.0+(exp(lE_gamma-lme-(2.0*lc))*(1+cos(phi))));
}

// returns the differential cross section for ic calculation 
double ldsig_dom2(double lE_gamma, double phi){
	return log(0.5)+2.0*lalphafs+2.0*lrc+2.0*log(Q(lE_gamma,phi))+log(Q(lE_gamma,phi)+(1.0/Q(lE_gamma,phi))-1+(pow(cos(phi),2)));
}

// weight function for ic scattering
double lweight(double lNe, double lE_gamma, double phi, double theta, double lx, double lE){	
	double gamma2 = exp(2*lE-2*lme-4*lc);
	double beta = sqrt(1-1/(gamma2));
	double lE_gamma_prime = log(exp(lE_gamma)*sqrt(gamma2)*(1-(beta*cos(theta))));
	double lvprime = lE_gamma - lplanck;
	
	double Q_nf = 1.0/(1.0+(exp(lE_gamma_prime-lme-(2.0*lc))*(1+cos(phi))));
	double ldsig_dom2_nf = log(0.5)+2.0*lalphafs+2.0*lrc+2.0*log(Q_nf)+log(Q_nf+(1.0/Q_nf)-1+(pow(cos(phi),2)));
		
	return lNe+ldsig_dom2_nf+ln_gamma_thin(lx,lvprime)-log(4*pi)+lc+log(1-beta*cos(theta))+2*log(pi)+log(fabs(sin(phi)))+log(fabs(sin(theta)));
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// computes initial values of parameters
int getParams(){

	lEj = lWj - lc;                       														// log of jet energy per meter in first meter of jet (in fluid frame)
	lUb = lA_equi+lEj-((2*lgamma_bulk)+log(A_equi+1));                       					// log of magnetic energy in jet (in fluid frame)
	lUe = lUb-lA_equi; 																			// log of electron energy in jet
	if(alpha==2){                      															// Different logs of A for different cases to avoid NaN
		lA = log(exp(lUe)/(lEmax-lEmin));                        								
	}																						
	else if(alpha > 2){
		lA = log(alpha-2.0)+lUe-log(pow(Emin,2.0-alpha)-pow(Emax,2.0-alpha)); 					
	}
	else{
		lA = log(2.0-alpha)+lUe-log(+pow(Emax,2.0-alpha)-pow(Emin,2.0-alpha));
	}
	lR0 = 0.5*(log(2.0)+lEj+lA_equi+lmu0-((2.0*lgamma_bulk)+lpi+(2.0*lB0)+log(1.0+A_equi)));    // log of radius of base of jet
	ldx = lL-lN;																				// log thickness of x-slice
	dx = L/N; 																					// thickness of x-slice
	lAx = lA;																					// log of A
	lEbinmin = -30;
	lEbinmax = -17;
	E_scatt_min = exp(lvmin_ic+lplanck);
	E_scatt_max = exp(lvmax_ic+lplanck);
	lE_scatt_min = lvmin_ic+lplanck;
	lE_scatt_max = lvmax_ic+lplanck;
	dE_scatt = (log(E_scatt_max)-log(E_scatt_min))/Nv;											// spacing between energy bins for scattered photons
	dE = (exp(lEbinmax)-exp(lEbinmin))/Nebins;
	ldE = log(dE);
	dtheta = (2*pi)/(Ntheta);																	// spacing in theta for ic calculation (angle between photon and electron before collision, not transformed)
	ldtheta = log(dtheta);																		// log of dtheta
	E_gamma0 = exp(lvmin_seed+lplanck);															// lowest photon energy for ic seed photons
	E_gamma_max = exp(lvmax_seed+lplanck);														// highest photon energy for ic seed photons
	dE_gamma = (E_gamma_max - E_gamma0)/NE_gamma;												// photon energy spacing
	ldE_gamma = log(dE_gamma);																	// log photon energy spacing
	dphi2 =2*pi/Nphi2;																			// spacing of phi2
	ldphi2 = log(dphi2);																		// log of dphi2
	lflux_factor_sync = log10(pi)+2*log10(3.08567758e22)+2*log10(d)+2*log10(tan(theta_opening))-2*log10(gamma_bulk);						
	lflux_factor_ic = log10(2*pi)+2*log10(d)+2*log10(3.08567758e22)+log10(sin(0.7922)*dphi2);
	dlE = (lEbinmax-lEbinmin)/Nebins;
	

	return 0;
} 

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// get input from input.dat
int getArgs(char *file){

        FILE *params;
        char *line = (char *)malloc(50*sizeof(char));

        params = fopen (file, "rt");
        while(fgets(line, 50, params) != NULL ){
                char *name = strtok(line, " \t\n");
                char *value = strtok(NULL, " \t\n");
                if(strcmp(name,"#") == 0){
                		// do nothing, its a comment
                }
                else if(strcmp(name,"gamma_bulk") == 0){
                        gamma_bulk = atof(value);
                }
                else if(strcmp(name,"L") == 0){
                        L = atof(value);
                }
                else if(strcmp(name,"theta_opening") == 0){
                        theta_opening = atof(value);
                }
                else if(strcmp(name,"Wj") == 0){
                        Wj = atof(value);
                }
                else if(strcmp(name,"A_equi") == 0){
                        A_equi = atof(value);
                }
                else if(strcmp(name,"alpha") == 0){
                        alpha = atof(value);
                }
                else if(strcmp(name,"B0") == 0){
                        B0 = atof(value);
                }
                else if(strcmp(name,"N") == 0){
                        N = atoi(value);
                }
                else if(strcmp(name,"Nv") == 0){
                        Nv = atoi(value);
                }
                else if(strcmp(name,"theta_observe") == 0){
                        theta_observe = atof(value);
                }
                else if(strcmp(name,"vmin_sync") == 0){
                        vmin_sync = atof(value);
                }
                else if(strcmp(name,"vmin_ic") == 0){
                        vmin_ic = atof(value);
                }
                else if(strcmp(name,"vmax_sync") == 0){
                        vmax_sync = atof(value);
                }
                else if(strcmp(name,"vmax_ic") == 0){
                        vmax_ic = atof(value);
                }
                else if(strcmp(name,"vmin_seed") == 0){
                        vmin_seed = atof(value);
                }
                else if(strcmp(name,"vmax_seed") == 0){
                        vmax_seed = atof(value);
                }
                else if(strcmp(name,"Emin") == 0){
                		Emin = atof(value);
                		Emin = Emin*1.6e-19;
				}
                else if(strcmp(name,"Emax") == 0){
                        Emax = atof(value);
                        Emax = Emax*1.6e-19;
                }
                else if(strcmp(name,"d") == 0){
                		d = atof(value);
                }
                else if(strcmp(name,"Nebins") == 0){
                		Nebins = atoi(value);
                }
                else if(strcmp(name,"Ntheta") == 0){
                        Ntheta = atof(value);
                }
                else if(strcmp(name,"Nphi2") == 0){
                        Nphi2 = atof(value);
                }
                else if(strcmp(name,"NE_gamma") == 0){
                        NE_gamma = atof(value);
                }
                else{
                        printf("Unrecognized input variable: %s.\n",name);
                }
                if(strcmp(name,"#") == 0){
                	// again, do nothing
                }
                else{
                	printf("%s = %s\n",name,value);
                }
        }

    // logarithmify everything

	lWj = log(Wj);
	lL = log(L);
	lB0 = log(B0);
	lA_equi = log(A_equi);
	lEmin = log(Emin);
	lEmax = log(Emax);
	lalpha = log(alpha);
	theta_opening = (pi/180.0)*theta_opening;		// in radians
	ltheta_opening = log(theta_opening);
	theta_observe = (pi/180.0)*theta_observe;		// in radians
	phi2 = (pi/180.0)*phi2;							// in radians
	ltheta_observe = log(theta_observe);
	lgamma_bulk = log(gamma_bulk);
	lN = log(N);
	lNv = log(Nv);
	lvmin_ic = log(vmin_ic);
	lvmax_ic = log(vmax_ic);
	lvmin_sync = log(vmin_sync);
	lvmax_sync = log(vmax_sync);
	lvmin_seed = log(vmin_seed);
	lvmax_seed = log(vmax_seed);
	dlv_sync = (lvmax_sync-lvmin_sync)/Nv;
	
	return 0;
}

// allocate arrays
int allocateArrays(){
	k_array = (double **)malloc(Nv*sizeof(double *));					// opacity by frequency
	lP_sync_array = (double **)malloc(Nv*sizeof(double *));				// synchrotron power by frequency (log)
	lNe_array = (double *)malloc(Nebins*sizeof(double));				// evolving electron population by energy (log)
	lNe0_array = (double *)malloc(Nebins*sizeof(double));				// initial electron population by energy (log)
	lE_array = (double *)malloc(Nebins*sizeof(double));					// electrion energies in an array (log)
	E_gamma_scatt_array_P = (double *)malloc(Nv*sizeof(double));		// ic photon powers by frequency 
	E_gamma_scatt_array_E = (double *)malloc(Nv*sizeof(double));		// energies associated with powers in ic photon power array
	x_array = (double *)malloc(N*sizeof(double));						// x slices along the jet
	v_sync_array = (double *)malloc(Nv*sizeof(double));					// frequencies at which to compute synchrotron powers
	lA_array = (double *)malloc(Nebins*sizeof(double));
	
	synchrotron = fopen("sync.dat","w");
	inverse_compton = fopen("ic.dat","w");
	population = fopen("pop-evol.dat","w");
	initPop = fopen("pop-init.dat","w");
	
	
	return 0;
}
