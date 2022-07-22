/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: Header file of the class used to convert the "Raw Plain Vectors" into the matrix-row-packing format (i.e., prepareVecMatMRP.h)   
	Abbreviation 		: MRP - Matrix Row Packing Format
	Note				: Commenting and cleaning up has been on July 2022  
*/

#ifndef PREPAREVECMATMRP_H
#define PREPAREVECMATMRP_H

/* Import the other classes' header files*/
#include "examples.h" 			  // The class for the console menu and guiding the user to the preferred application     	
#include "printCont.h" 			  // The class containing the printing functions for control purposes  
#include "rawplain.h"  			  // The class which reads and stores the plain matrices                    
#include "generateplaintextMRP.h" // The functions used to prepare plain matrices in MRP format
#include "prepareVecMatMRP.h"     // The class which prepares the read matrices in MRP format  		
#include "initializationPLCP.h"   // The class containing the functions which convert the prepared matrices to Plaintext and Ciphertext objects.            
#include "plcpOperations.h"       // The class containing the functions which does Ciphertext-Ciphertext and Plaintext-Ciphertext arithmetic and algebraic operations  	    
#include "encryptedAppx.h"        // The class containing the functions which does the Chebyshev approximation with different assumptions     
#include "secretShare.h" 		  // The class containing the functions which does the secret sharing  

/* Import the important selected C libraries*/
#include <array>
#include <cmath>
#include <vector>
#include <iostream>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <unistd.h>

/* Call main namespaces */
using namespace std;
using namespace seal;

/* Main struct for holding the simulation data in MRP format*/
struct simulationMatrixMRP{

	// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =====
	// ==== State-Input Dimensions, Vectors and Matrices  ==== 
	// ==== ==== ==== ==== ==== ==== ====  ==== ==== ==== ====
	// States, inputs, Encrypted single long row size 
	int m; // Inputs
	int n; // States
	int N; // Encypted Row Size	
	// System matrices and vectors
	std::vector<double> * AA_MRP;
 	std::vector<double> * BB_MRP; 
	std::vector<double> * CC_MRP; 
	std::vector<double> * FF_MRP; 
	std::vector<double> * WW_MRP;
	std::vector<double> * VV_MRP;
 	std::vector<double> * KK_MRP;
	std::vector<double> * LL_MRP; 
	std::vector<double> * xrxr_MRP; 
	std::vector<double> * urur_MRP; 
	std::vector<double> * GAMMA_MRP; 
	std::vector<double> * xGxG_MRP; 
	std::vector<double> * uGuG_MRP; 
	std::vector<double> * uGuG_AS_MRP; 
	std::vector<double> * vv_MRP; 
	std::vector<double> * TAU_MRP; 
	std::vector<double> * ALARMSYS_MRP; 
	std::vector<double> * ss_MRP; 
	std::vector<double> * xx_MRP; 
	std::vector<double> * xexe_MRP; 
	std::vector<double> * xpxp_MRP;
	std::vector<double> * yy_MRP;
	std::vector<double> * yyAS_MRP;
	// Pre-computed matrices
	std::vector<double> * KGKG_MRP;     // [2][10]  when m = 2, n = 10 
	std::vector<double> * KLKL_MRP;     // [2][10]  when m = 2, n = 10  	
	std::vector<double> * KxuGKxuG_MRP; // [2][1]   when m = 2, n = 10
	std::vector<double> * ACL_MRP;		// [10][10] when m = 2, n = 10 	
	std::vector<double> * KK_Minus_MRP;	 	
	std::vector<double> * KxKx_MRP;	

	// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =======
	// ==== Chebyshev  Approximation Constants and Computation Arrays ==== 
	// ==== ==== ==== ==== ==== ==== ====  ==== ==== ==== ==== ==== ======

	// Degrees of Chebyshev Polynomials   
	int chebDegEq8;
	int chebDegEq9;
	// Chebyshev Approximation Upper and Lower bounds of data
	double alpbetLowBouEq8;
	double alpbetUpBouEq8; 
	double alpbetLowBouEq9;
	double alpbetUpBouEq9;	
	// Chebyshev Approximation Alpha and Beta values
	double alpEq8; // Max Equation 	  
	double betEq8; // Max Equation 
	double alpEq9; // Indicator subtraction method	  
	double betEq9; // Indicator subtraction method

	// Chebyshev aplha and beta vectors  
	std::vector<double> * alpEq8_MRP; // Max alpha array 	
	std::vector<double> * betEq8_MRP; // Max beta array 
	std::vector<double> * alpEq9_MRP; // Indicator Subtraction alpha array  	
	std::vector<double> * betEq9_MRP; // Indicator Subtraction beta array 
	// Chebyshev One Vector
	std::vector<double> * One_MRP; // One vector (T0 in Chebyshev Appx. computation)  	
	// Chebyshev First (1st) Term Vectors (to-be-used in the final addition)   	
	std::vector<double> * chebPowSerFT_Eq8_MRP; // Chebyshev Equation-8 Vector  	
	std::vector<double> * chebPowSerFT_Eq9_MRP; // Chebyshev Equation-9 Vector  	
	// Chebyshev Power Series Coefficients    	
	std::vector<double> * chebPowSer_Coeff_Eq8_MRP; // Chebyshev Equation-8 Vector
	std::vector<double> * chebPowSer_Coeff_Eq9_MRP; // Chebyshev Equation-9 Vector

};

/* Constructor (function) of the class */
void prepareVecMatMRP();

/* The function for calculating the each row vector size including the zero trials */
void calculateVectorSize(size_t x_vector_row_size, size_t * calculated_x_row_size);

/* The function for creating the simulation matrices to convert  the simulation matrix data into MRP  */
void create_SimulationMatrixDataMRP(struct simulationMatrixData *smd, struct simulationMatrixMRP *smrp);

/* The function for generating the the simulation matrices in MRP based on the simulation data */
void assignValMatrixDataMRP(struct simulationMatrixData *smd, struct simulationMatrixMRP *smrp);

#endif

/* 
// Exemplary Column and Row Dimensions of the simulation matrices and vectors when m = 2, n = 10, N = 16
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
