/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: The header file of the class used to handle the cleartext (or raw) 
						  simulation data with the initializations (i.e., rawplain.h) 
	Note				: Commenting and cleaning up has been on July 2022
*/

#ifndef RAWPLAIN_H
#define RAWPLAIN_H

/* Main struct for holding the simulation data */
struct simulationMatrixData {

	// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =====
	// ==== State-Input Dimensions, Vectors and Matrices  ==== 
	// ==== ==== ==== ==== ==== ==== ====  ==== ==== ==== ====
	// States, inputs, number of iterations 
	int m; 
	int n;
	int tMax;
	// State matrices 
	double ** AA;
 	double ** BB; 
	double ** CC; 
	double ** FF; 
	double ** WW;
	double ** VV;
 	double ** KK; 
	double ** LL; 
	double ** xrxr; 
	double ** urur; 
	double ** GAMMA; 
	double ** xGxG; 
	double ** uGuG; 
	double ** vv; 
	double ** TAU; 
	double ** ALARMSYS; 
	double ** ss; 
	double ** xx; 
	double ** xexe; 
	double ** xpxp;
	double ** yy;
	// Pre computed Matrices and Vectors
	double ** KGKG; 	// Late change based on Andreea's homomorphic optimization in 3rd equation (proposed quite close to the initial submission)
	double ** KLKL; 	// Late change based on Andreea's homomorphic optimization in 3rd equation (proposed quite close to the initial submission)	
	double ** KxuGKxuG; // Late change based on Andreea's homomorphic optimization in 3rd equation (proposed quite close to the initial submission)
	double ** ACL;		
	double ** KxKx;       
	// Noise vectors 
	double ** xNoise;
	double ** yNoise;

	// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =======
	// ==== Chebyshev  Approximation Constants and Computation Arrays ==== 
	// ==== ==== ==== ==== ==== ==== ====  ==== ==== ==== ==== ==== ======
	
	// Degrees of Chebyshev Polynomials   
	int chebDegEq8;
	int chebDegEq9;
	// Upper and lower bounds of two function's observed data	
	double alpbetLowBouEq8; // The maximum (or RELU) function's observed data's min value (namely a) 
	double alpbetUpBouEq8;  // The maximum (or RELU) function's observed data's max value (namely b) 
	double alpbetLowBouEq9; // The Indicator function's observed data's min value (namely a)
	double alpbetUpBouEq9;  // The Indicator function's observed data's min value (namely b)
	// 'Alpha' and 'beta + 1' values computed based on the observed data's max (a) and min (b) values 
	double alpEq8; // Max Equation 	  
	double betEq8; // Max Equation 
	double alpEq9; // Indicator subtraction method	  
	double betEq9; // Indicator subtraction method

	// ==========================================
	// Chebyshev Approximation Computation Arrays
	// ==========================================   	
	// Vector One Array
	double ** One; // Last Addition 06.03.2022	// Note: N could be added 	 
	// Equation-8-Max Arrays 
	double ** eq8maxAppx_PS_FT_D12_y_10_u_2;    // First Term  
	double ** eq8maxAppx_PS_Coeff_D12_y_10_u_2; // Coefficient Array	
	// Equation-9-Ind-Subtraction Arrays
	double ** eq9ISubAppx_PS_FT_D12_y_10_u_2;    // First Term
	double ** eq9ISubAppx_PS_Coeff_D12_y_10_u_2; // Coefficient Array	

	// ==== ==== ==== ==== ==== ==== ==== ==
	// ==== Experimental Result Vectors ==== 
	// ==== ==== ==== ==== ==== ==== ======= 
	// Experimental Result Vectors
	double ** xe_Res;
	double ** u_Res;
	double ** xp_Res;
	double ** residue_Res;
	double ** sBar_Res;
	double ** indInp_Res; 
	double ** alarm_Res;
	double ** s_Res;
	double ** x_Res;
	double ** y_Res;
};

// ==== ==== ==== ==== ==== ==== ==== ==== ==== ======= 
// ==== Functions for the Simulation Data Creation ==== 
// ==== ==== ==== ==== ==== ==== ==== ==== ==== =======
/* The function for creating the simulation matrix data with the empty slots */
void create_SimulationMatrixData(struct simulationMatrixData *smd);

/* The function for assigning the column and row dimensions for the simulation (i.e., m and n respectively) */
void assDimToSmd(struct simulationMatrixData *smd, char * DimFolderDir, int numIter);

/* The function for reading and assigning the values of the initial raw (cleartext) vectors and matrices */
void assignValMatrixDatabyFileRead(struct simulationMatrixData *smd, char * folderPath);

/* The function for initializing the remaining (not-read from the records received) vectors including the sensor measurement vector y 
   Status: Possibly not-used */
void initRemainVec(struct simulationMatrixData *smd);

/* The function for assigning the Chebyshev Approximation constants and arrays for the 8th Equation (Max Function) and 9th Equation (Subtraction-based Indicator Function) */
void assignCUSUMChebyshevAppxParams(struct simulationMatrixData *smd);

// ==== ==== ==== ==== ==== ==== ==== ==== ==== ====== 
// ==== ==== Functions for Reading Data ==== ==== ==== 
// ==== ==== ==== ==== ==== ==== ==== ==== ==== ======
/* The function for reading a file  (Probably not-used) */
void readFile(double **data, std::string input, int *size);

/* The function for reading and filling the simulation vectors and matrices  */
void readMatrix(struct simulationMatrixData *smd, char * dataPath, char * matname, int dim1, int dim2); 

/* The function for reading and filling the recorded y vector data */
void readYData(struct simulationMatrixData *smd, char * directory);

/* The function for checking the content of the y vector */
void readYMCheck(struct simulationMatrixData *smd, char * matname, int dim1, int dim2);

/* The function for printing the experimental result (both the sensor and CUSUM statistics) */
void printExperimentalResult(struct simulationMatrixData *smd);

// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ===== 		 
// ==== Functions Used for Sensor Measurement and Noise additions ==  
// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =====
/* The function for assigning the resulting x vector */
void assignXResult(struct simulationMatrixData *smd, double ** result, int rowDm, int colDm);

/* The function for poerforming the matrix multiplication */
void matrixMult(double ** result, double ** vect1, double ** vect2, int rowDm1, int colDm1, int rowDm2, int colDm2);

/* The function for adding three matrices in a row  */
void matrixAdditionThree(double ** result,  double ** vect1, double ** vect2, double ** vect3, int rowDm, int colDm);

/* The function for copying a result vector to a target vector */
void assignResult(double ** result, double ** target, int rowDm, int colDm);

/* The function for printing a matrix content */
void printMatrix(double **matrix, int rowDm, int colDm);

#endif

/* 
// Exemplary Column and Row Dimensions of the simulation matrices and vectors when m = 2, n = 10
A[10][10] 		= n * n 
B[10][2]  		= n * m 
C[10][10] 		= n * n
F[10][4]  		= n * 4
W[4][4]   		= 4 * 4
V[10][10] 		= n * n
K[2][10]  		= m * n 
L[10][10] 		= n * n 
xr[10][1] 		= n * 1
ur[2][1]      	= m * 1
Gamma[10][10] 	= n * n 
xG[10][1] 		= n * 1 
uG[2][1]  		= m * 1 
v[10][1]  		= n * 1 
tau[10][1] 		= n * 1 
alarmSys[10][1] = n * 1 
s[10][1] 		= n * 1 
x[10][1] 		= n * 1 
xe[10][1] 		= n * 1 
xp[10][1] 		= n * 1 
int m 			= 2;   // It might change to 4, 10
int tMax 		= 120; // A selected number of iterations 
int n 			= 10;  // It might change to 20, 50 
*/
