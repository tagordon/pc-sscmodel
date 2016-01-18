// constants

#define euler 2.7182818284										// e
#define pi 3.14159265358979										// pi
#define mu0 1.2566370614e-6 									// magnetic constant
#define c 2.99e8												// speed of light
#define me 9.11e-31												// mass of electron
#define e_charge 1.602e-19										// fundamental charge
#define sigma_t 6.6524e-29										// the thomson cross section

#define lplanck10 -33.1787439778								// log10 of planck's constant
#define le_charge10 -18.7952896101								// log10 of fundamental charge

#define leuler  1 												// log e
#define lpi  1.144729886										// log pi
#define lmu0  -13.5870714										// log magnetic constant
#define lc  19.518601											// log c
#define lme  -69.1708328401										// log electron mass
#define le_charge  -43.2777536746								// log fundamental charge
#define lsigma_t  -64.87999001									// log thompson cross section
#define lplanck -76.3969										// log planck constant
#define lrc -28.58329											// log compton radius
#define lalphafs -4.92024366									// log fine structure constant
#define lconvert 0.4342944819									// multiply to change from natural log to base 10 log
#define l2 0.6931471806
#define l3 1.098612289
#define l4 1.386294361

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// file pointers

FILE *synchrotron;
FILE *inverse_compton;
FILE *population;
FILE *initPop;

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// array declarations and allocations

double **k_array;					// opacity by frequency
double **lP_sync_array;				// synchrotron power by frequency (log)
double *lNe_array;					// evolving electron population by energy (log)
double *lNe0_array;					// initial electron population by energy (log)
double *lE_array;					// electrion energies in an array (log)
double *lx_array;					// x slices along the jet
double *lv_sync_array;				// frequencies at which to compute synchrotron powers
double *lA_array;
double *lE_scatt_array;

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// variables set from input.dat

int nthreads;													// number of threads to use
double gamma_bulk;              								// lorenz factor of jet fluid
double L;                        								// length of jet
double theta_opening;            								// opening angle of the jet
double Wj;                       								// total power of jet in lab frame in 10^36 W
double alpha;                    								// other constant in electron population
double B0;                       								// magnetic field at the base of the jet in 10^-5 Tesla
int  N;                         								// the number of slices of width dx in the jet
int Nv;                         								// the number of frequencies to compute for both synchrotron and inverse compton spectra
int Nebins;														// number of bins for electron energies
double theta_observe;            								// the observation angle for the jet in degrees
double Emin;                     								// minimum electron energy
double Emax;                     								// maximum electron energy
double A_equi;													// ratio of electron to magnetic field energy
double vmin_sync;												// minimum frequency for synchrtron spectrum
double vmax_sync;												// maximum frequency for synchrotron spectrum
double E_gamma_max;												// maximum seed photon energy
double E_gamma_min;												// minimum seed photon energy
double lE_scatt_min;
double lE_scatt_max;

int Ntheta;														// number of theta divisions
int Nphi2;														// number of phi2 divisions
int NE_gamma;													// number of photon energy divisions

double lflux_factor_sync;										// conversion factor from power to flux 
double lflux_factor_ic;

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// logs of parameters set by input.dat

double lgamma_bulk;              								
double lL;                        								
double ltheta_opening;            								
double lWj;                       								
double lalpha;                    								
double lB0;                       								
double lN;   	                      						
double lNv;                         							
double lNebins;													
double ltheta_observe;            								
double lEmin;                     								
double lEmax;                     								
double lA_equi;													
double lvmin_sync;
double lvmin_ic;
double lvmax_sync;
double lvmax_ic;
char *do_sync;
char *do_ic;

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// other declarations

int i,j,k,l,m,n;		// just some integers...
double theta;			// initial photon-electron angle for ic scattering in fluid frame
double theta;			// angle between electron velocity and seed photon velocity
double phi2;			// angle that photon is scattered through for ic scattering in electron rest frame
double E_gamma;			// seed photon energy
double lE;				// initial electron energy for ic scattering in fluid frame
double ldE;				// log spacing in electron energies
double dlE;				// spacing of log in electron energies
double vmin_seed;		// smallest ic seed photon energy
double vmax_seed;		// largest ic seed photon energy
double lvmin_seed;		// log smallest ic seed photon energy
double lvmax_seed;		// log largest ic seed photon energy

int sync_index;
int ic_index;

double x;
double lx;
double v;
double lv;
int Ei;
double tau_x;
double factor;

double lE_scattered;	// scattered photon energy for ic scattering
int scattIndex;			// index of scattered photon energy in E_gamma_scatt_arrays
double lweightd;		// log of weight function (with infintesimals)
double dtheta;
double dphi2;
double dE_gamma;
double dE_scatt;

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculated parameters

double Ej;                       								// jet energy per meter in first meter of jet (in fluid frame)
double Ub;                       								// magnetic energy in jet (in fluid frame)
double Ue;                       								// electron energy in jet (in fluid frame)
double A;                        								// contant in electron population equation
double R0;                       								// radius of base of jet
double dx;                       								// width of slabs the jet is divided into
double dv;                       								// width of frequency divisions
double dE;														// width of electron energy bins
double lEbinmin;												// minimum electron energy bin
double lEbinmax;												// maximum electron energy bin
double beta;                     								// v/c
double gamma2;													// square of relativistic factor gamma
double constants;												// just so I don't have to declare it in a function
double Ax;														// constant A which changes with x as electron population evolves
double d;														// distance to source in megaparsecs

double lEj;                       								// log of jet energy per meter in first meter of jet (in fluid frame)
double lUb;                       								// log of magnetic energy in jet (in fluid frame)
double lUe;                       								// log of electron energy in jet (in fluid frame)
double lA;                        								// log of contant in electron population equation
double lR0;                       								// log of radius of base of jet
double ldx;                       								// log of width of slabs the jet is divided into
double ldv;                       								// log of width of frequency divisions
double dlv_sync;												// width of log frequency divisions. Haha. 
double lbeta;                     								// log of v/c
double lAx;														// log of Ax

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Mathematical function prototypes

double lR(double x);                                                								// returns the radius of the jet at position x
double lB(double x);                                                								// returns the magnetic field in the jet at position x
double leps(double x);                                    											// the quantity epsilon that shows up in several equations
double lEe(double x,double v);                                     									// electron energy that emits synchrotron radiation at the critical frequency
double lj(double x, double v, double Ax);                          									// j_0 for calculating k_v(v,x)
double lk(double x, double v, double Ax);                           								// opacity for photons of frequency v at position x in the jet
double ltau(double v, double x, int i, double *k_array);            								// optical depth of slice x for frequency v
double lP(double x, double v, double Ax);                           								// power lost at slice x by frequency v
double lNe0(double v);                                              								// electron energy at x=0 as a function of frequency
double lboost(double lP);																			// dopplar boosts synchrotron emissions
double simple_integral(double *f,int i, int n, double dx);         									// integrates an array using the midpoint rule
double simple_integral_simpson(double *f,int i, int n, double dx);         							// integrates an array using simpson's rule
double simple_integral_comp_simpson(double *f,int i, int n, double dx);								// integrates an array using composite simpsons
int getIndex_Ee(double lE);																			// returns the index of the electron energy bin for the given energy
int getIndex_scatt(double lE);																		// returns the index of the scattered photon energy bin for the given energy
double getlE_Ee(int Ei);																			// returns energy of a given electron energy bin index
double getlE_scatt(int Ei);																			// returns the power scattered by a photon in the given scattered photon energy bin
double getlv(double lx, double lE);																	// returns frequency at which electrons emit synchrotron radiation for E,x
double ln_gamma_thin(double lx, double lv);															// photon population for optically thin case
double ln_gamma_thick(double lx, double lv);														// photon population for optically thick case
double Q(double lE_gamma, double phi);																// function for cross section
double ldsig_dom2(double lE_gamma, double phi);														// differential cross section for electron, photon scattering
double lweight(double lNe, double lE_gamma, double phi, double theta, double lx, double lE);		// weight function for ic scattering
double lE_scatt(double lE, double lEe, double theta, double phi2);									// returns scattered photon energy for ic scattering taking place at angles theta and phi2 with energies E (photon) and Ee (electron)
double phi_scatt(double lE,double lEe, double theta, double phi2);									// returns the scattered angle phi for ic scattering with the given parameters
int getIndex_Ee_log(double lE);
double logdE(int Ei);
double get_ldE_scatt(int i);
double lvboost(double lv);
double lwidth(int i);
double losses(double lx,double lv,double lAx);
double get_lv(double lx, double lEe);
void *thread_ic(void *vargp);
void *thread_sync(void *vargp);



//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Other functions

int allocateArrays();											// allocates arrays
int getArgs();                  								// gets input from file
int getParams();                								// calculates parameters
int computeSynchrotron();       								// calculates synchrotron spectrum