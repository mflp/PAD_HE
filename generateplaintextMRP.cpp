/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which is used for converting the "Raw Plain (or cleartext) Vectors" in the matrix-row-packing format again in cleartext   	
	Abbreviation 		: MRP - Matrix Row Packing Format
	Note				: Commenting and cleaning up has been on July 2022
*/

/* Import the other classes' header files */
#include "examples.h" 			  // The class for the console menu and guiding the user to the preferred application     	
#include "rawplain.h"  			  // The class which reads and stores the plain matrices                    
#include "generateplaintextMRP.h" // The functions used to prepare plain matrices in MRP format
#include "prepareVecMatMRP.h"     // The class which prepares the read matrices in MRP format  		
#include "initializationPLCP.h"   // The class containing the functions which convert the prepared matrices to Plaintext and Ciphertext objects.            
#include "plcpOperations.h"       // The class containing the functions which does Ciphertext-Ciphertext and Plaintext-Ciphertext arithmetic and algebraic operations  	    
#include "encryptedAppx.h"        // The class containing the functions which does the Chebyshev approximation with different assumptions     
#include "secretShare.h" 		  // The class containing the functions which does the secret sharing  
#include "printCont.h" 			  // The class containing the printing functions for control purposes  

/* Import the important selected C libraries */
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
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

/*
	Function for generating the s vector in MRP for the homomorphic addition (i.e., s vector)  	
*/
void genVecAddOperMRP_RPL(std::vector<double> *x_add_vector, struct simulationMatrixData *smd, size_t x_vector_row_size, size_t x_vector_repeat, size_t calculated_x_row_sizeVal, char *matname){
	
	// Read and fill the entries (ending with trailing zeros)
	for (size_t i = 0; i < calculated_x_row_sizeVal * x_vector_repeat; i++){
		if(i % calculated_x_row_sizeVal == 0){ // Fill the vector based on the read simulation matrices in a way that each entry in the vector would be in the first indices of the row segments (e.g., 0, 16, 32 so and so forth)	
			int quot = i/calculated_x_row_sizeVal;			
			if(strcmp(matname, "s") == 0){
				smd->ss[quot][0];	
			}
		}else{ // Fill the trailing zeros for the remaining column/row segment entries 
			x_add_vector->at(i) = 0;	
		} 		
	}
}

/*
	The function for generating the x vector of generic Ax multiplication in MRP (i.e., x^e, y, uGamma vectors)
*/
void genRepXVecMRP_RPL(std::vector<double> *x_vector, struct simulationMatrixData *smd, size_t x_vector_col_size, size_t x_vector_repeat, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted){
		
	// Read the x vector (of the generic Ax multiplication) content into the first n indices based on the selected vector name
	for (size_t i = 0; i < x_vector_col_size; i++){
		if(strcmp(matname, "xe") == 0)
			x_vector->at(i) = smd->xexe[i][0]; 	
		if(strcmp(matname, "y") == 0)
			x_vector->at(i) = smd->yy[0][i];
		if(strcmp(matname, "uG") == 0)
			x_vector->at(i) = smd->uGuG[i][0];

	}

	// Fill the trailing zeros for the remaining column/row segment entries
	for (size_t j = x_vector_col_size; j < calculated_x_row_sizeVal; j++)		
		x_vector->at(j) = 0; 
 		
	// Repeat the same content as many as number of rows times
	for (size_t i = calculated_x_row_sizeVal; i < calculated_x_row_sizeVal * x_vector_repeat; i++)	
		x_vector->at(i) = x_vector->at(i % calculated_x_row_sizeVal);
	
	// If there is a request for vector content checking, print the vector content 	
	if(isPrinted)
		printVector(x_vector, calculated_x_row_sizeVal, x_vector_repeat, matname);
}

/*
	The function for generating the y vector (representing the sensor measurement) as the x vector of generic Ax multiplication in MRP (i.e., y vector)
*/
void genRepXVecMRP_RPL_v0(std::vector<double> *x_vector, double ** yy, size_t x_vector_col_size, size_t x_vector_repeat, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted){
		
	// Read the entries of y vector
	for (size_t i = 0; i < x_vector_col_size; i++)
		x_vector->at(i) = yy[i][0];

	// Fill the trailing zeros for the remaining column/row segment entries
	for (size_t j = x_vector_col_size; j < calculated_x_row_sizeVal; j++)		
		x_vector->at(j) = 0;// x_vector->push_back(0); // Vector class is used as undefined size??? 		

	// Repeat the same content as many as number of rows times
	for (size_t i = calculated_x_row_sizeVal; i < calculated_x_row_sizeVal * x_vector_repeat; i++)	
		x_vector->at(i) = x_vector->at(i % calculated_x_row_sizeVal);

	// If there is a request for vector content checking, print the vector content 	
	if(isPrinted)
		printVector(x_vector, calculated_x_row_sizeVal, x_vector_repeat, matname);
}

/*
	The function for generating the matrices vector as the A matrix of generic Ax multiplication in MRP (i.e., Gamma, L, ACL, B matrices)
	Important note: This packing function is valid for n * n, n * m  matrices where m < n (e.g., m = 2, n = 10) 
*/
void genMatMRP_RPL(std::vector<double> *matPtr, struct simulationMatrixData *smd, size_t mat_row_size, size_t mat_col_size, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted) {

	// Read and fill the entries (ending with trailing zeros)
	for (size_t i = 0; i < mat_row_size; i++){
		for (size_t j = 0; j < calculated_x_row_sizeVal; j++){
			// Read the entries of the selected matrix
			if( j < mat_col_size){
				if(strcmp(matname, "Gamma") == 0)
					matPtr->at( i * calculated_x_row_sizeVal + j) = smd->GAMMA[i][j];
				if(strcmp(matname, "L") == 0)
					matPtr->at( i * calculated_x_row_sizeVal + j) = smd->LL[i][j];
				if(strcmp(matname, "ACL") == 0)
					matPtr->at( i * calculated_x_row_sizeVal + j) = smd->ACL[i][j];
				if(strcmp(matname, "B") == 0)	
					matPtr->at( i * calculated_x_row_sizeVal + j) = smd->BB[i][j];
			}else{ // Fill the trailing zeros for the remaining column/row segment entries
				matPtr->at(i * calculated_x_row_sizeVal + j) = 0;
			}
		}	
	}

	// If there is a request for vector content checking, print the vector content 	
	if(isPrinted)
		printMatrix(matPtr, mat_row_size, calculated_x_row_sizeVal, matname);	
}

/*
	The function for generating the matrices vector as the A matrix of generic Ax multiplication in MRP (i.e., KGamma, KL, KMinus)
	Important note: This packing function is valid for m * n matrices where m < n (e.g., m = 2, n = 10)  
*/
void genMatforMtimesNMRP_RPL(std::vector<double> *matPtr, struct simulationMatrixData *smd, size_t numCol, size_t numRow, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted) {
	
	// Read and fill the entries (ending with trailing zeros)  
	for (size_t i = 0; i < numCol; i++){
		for (size_t j = 0; j < calculated_x_row_sizeVal; j++){
			if(i < numRow && j < numCol){ // Read the entries of the selected matrices for the number of rows times  
				if(strcmp(matname, "KG") == 0)
					matPtr->at(i * calculated_x_row_sizeVal + j) = smd->KGKG[i][j];
				if(strcmp(matname, "KL") == 0)
					matPtr->at( i * calculated_x_row_sizeVal + j) = smd->KLKL[i][j];
				if(strcmp(matname, "KMinus") == 0)
					matPtr->at( i * calculated_x_row_sizeVal + j) = -1 * smd->KK[i][j];
			}else{ // Fill the trailing zeros 
				matPtr->at(i * calculated_x_row_sizeVal + j) = 0;
			}
		}	
	}

	// If there is a request for vector content checking, print the vector content 	
	if(isPrinted)
		printMatrix(matPtr, numCol, calculated_x_row_sizeVal, matname);	
}

/*
	Function for generating the xGamma vector in MRP for homomorphic addition (i.e., xGamma vector)  	
*/
void genXVecAddOperMRP_RPL(std::vector<double> *x_add_vector, struct simulationMatrixData *smd,  size_t x_vector_row_size, size_t x_vector_repeat, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted){

	// Read and fill the entries (ending with trailing zeros)		
	for(size_t i = 0; i < calculated_x_row_sizeVal * x_vector_repeat; i++){
		if(i % calculated_x_row_sizeVal == 0){ // Fill the first indices of the row/column segments with the entries of the selected vector 
			int quot = i/calculated_x_row_sizeVal;
			if(strcmp(matname, "xG") == 0){
				x_add_vector->at(i) = smd->xGxG[quot][0]; 
			}
		}else{ // Fill the trailing zeros 
			x_add_vector->at(i) = 0;	
		} 		
	}

	// If there is a request for vector content checking, print the vector content 	
	if(isPrinted)
		printVector(x_add_vector, calculated_x_row_sizeVal, x_vector_repeat, matname);
}

/*
	Function for generating the vectors in MRP for homomorphic addition (i.e., KxuG, uG, Kx vectors)  	
*/
void genUVecAddOperMRP_RPL(std::vector<double> *u_add_vector, struct simulationMatrixData *smd, size_t numRow, size_t u_vector_repeat, size_t calculated_u_row_sizeVal, char *matname, bool isPrinted){

	// Read and fill the entries (ending with trailing zeros)
	for (size_t i = 0; i < calculated_u_row_sizeVal * u_vector_repeat; i++){
		if(i < calculated_u_row_sizeVal * numRow){ // Fill the first indices of the row/column segments with the entries of the selected vector 
			if(i % calculated_u_row_sizeVal == 0){  
				int quot = i/calculated_u_row_sizeVal;			
				if(strcmp(matname, "Kxug") == 0){
					u_add_vector->at(i) = smd->KxuGKxuG[quot][0];	
				}
				if(strcmp(matname, "uG") == 0){
					u_add_vector->at(i) = smd->uGuG[quot][0];	
				}
				if(strcmp(matname, "Kx") == 0){
					u_add_vector->at(i) = smd->KxKx[quot][0];	
				}
			}else{ // Fill the trailing zeros 
				u_add_vector->at(i) = 0;	
			} 		
		}else{// Fill the trailing zeros 
			u_add_vector->at(i) = 0;
		}
	}

	// If there is a request for vector content checking, print the vector content	
	if(isPrinted)
		printVector(u_add_vector, calculated_u_row_sizeVal, u_vector_repeat, matname);

}

/*
	Function for generating the vectors in MRP for homomorphic addition (i.e., y, xp, s, v, tau, one, Chebyshev vectors)  	
*/
void genYVecAddOperMRP_RPL(std::vector<double> *y_add_vector, struct simulationMatrixData *smd, size_t numIter, size_t y_vector_row_size, size_t y_vector_repeat, size_t calculated_y_row_sizeVal, char *matname, bool isPrinted){
	
	// Read and fill the entries (ending with trailing zeros)	
	for(size_t i = 0; i < calculated_y_row_sizeVal * y_vector_repeat; i++){
		// Fill the first indices of the row/column segments with the entries of the selected vector 
		if(i % calculated_y_row_sizeVal == 0){ 
			int quot = i/calculated_y_row_sizeVal;
			if(strcmp(matname, "y") == 0){
				y_add_vector->at(i) = smd->yy[numIter][quot];					
			}
			if(strcmp(matname, "xp") == 0){
				y_add_vector->at(i) = smd->xpxp[quot][0];					
			}
			if(strcmp(matname, "s") == 0){
				y_add_vector->at(i) = smd->ss[quot][0];				
			}
			if(strcmp(matname, "v") == 0){
				y_add_vector->at(i) = smd->vv[quot][0];				
			}
			if(strcmp(matname, "tau") == 0){
				y_add_vector->at(i) = smd->TAU[quot][0];				
			}
			if(strcmp(matname, "one") == 0){		
				y_add_vector->at(i) = smd->One[quot][0];		
			}
			if(strcmp(matname, "chb_1st_T_D_12_m_2_n_10") == 0){		
				y_add_vector->at(i) = smd->eq8maxAppx_PS_FT_D12_y_10_u_2[quot][0];
			}
			if(strcmp(matname, "chb_Coeff_D_12_m_2_n_10") == 0){		
				y_add_vector->at(i) = smd->eq8maxAppx_PS_Coeff_D12_y_10_u_2[quot][0];
			}
			if(strcmp(matname, "chb_2nd_T_D_12_m_2_n_10") == 0){		
				y_add_vector->at(i) = smd->eq9ISubAppx_PS_FT_D12_y_10_u_2[quot][0];
			}
			
		}else{ // Fill the trailing zeros 
			y_add_vector->at(i) = 0;	
		} 		
	}

	// If there is a request for vector content checking, print the vector content	
	if(isPrinted)
		printVector(y_add_vector, calculated_y_row_sizeVal, y_vector_repeat, matname);
}

/*
	Function for generating the y vector in MRP for homomorphic addition (i.e., y vector)  	
*/
void genYVecAddOperMRP_RPL_v0(std::vector<double> *y_add_vector, double ** yy, size_t y_vector_row_size, size_t y_vector_repeat, size_t calculated_y_row_sizeVal, char *matname, bool isPrinted){

	// Read and fill the entries (ending with trailing zeros)	
	for(size_t i = 0; i < calculated_y_row_sizeVal * y_vector_repeat; i++){
		// Fill the first indices of the row/column segments with the entries of the y vector	
		if(i % calculated_y_row_sizeVal == 0){
			int quot = i/calculated_y_row_sizeVal;
			y_add_vector->at(i) = yy[quot][0];						
		}else{ // Fill the trailing zeros 
			y_add_vector->at(i) = 0;	
		} 		
	}

	// If there is a request for vector content checking, print the vector content	
	if(isPrinted)
		printVector(y_add_vector, calculated_y_row_sizeVal, y_vector_repeat, matname);
}

/*
	Function for generating the range transformation vector of a Chebyshev Approximation in MRP for homomorphic addition (i.e., Chebyshev Approximation vector)  	
*/
void genRangTransfVecChebApprx(std::vector<double> * Alpha_vector, std::vector<double> *Beta_vector, double lowerbound, double upperbound, size_t calculated_mask_vector_sizeVal, size_t cheb_vector_repeat){

	// Chebyshev range transformation values from [a,b] to [-1, 1]
	// double a 	 = -3.369997; // Exemplary values 
	// double b 	 = 32.798565; // Exemplary values
    double a 		= lowerbound;
	double b 		= upperbound;
	double alpha 	= 2 / (b - a);
	double beta  	= 2 * a / (b - a);
		
	 	
	for(int i = 0; i < calculated_mask_vector_sizeVal * cheb_vector_repeat; i++){
		if(i % calculated_mask_vector_sizeVal == 0) // Fill the first indices of the row/column segments with the alpha value of the Chebyshev approximation
			Alpha_vector->at(i) = alpha;
		else  // Fill the trailing zeros 
			Alpha_vector->at(i) = 0; 		
	}

	// Fill the first indices of the row/column segments with the alpha value of the Chebyshev approximation
	for(int i = 0; i < calculated_mask_vector_sizeVal * cheb_vector_repeat; i++){
		if(i % calculated_mask_vector_sizeVal == 0)
			Beta_vector->at(i) = beta + 1;
		else // Fill the trailing zeros 
			Beta_vector->at(i) = 0; 		
	}

}

/* 
// Exemplary Column and Row Dimensions of the simulation matrices and vectors when m = 2, n = 10

int m 			= 2;   // It might change to 4, 10
int tMax 		= 120; // A selected number of iterations 
int n 			= 10;  // It might change to 20, 50 

A[10][10] 		= n * n 
V[10][10] 		= n * n
L[10][10] 		= n * n 
Gamma[10][10] 	= n * n 
C[10][10] 		= n * n

B[10][2]  		= n * m 
K[2][10]  		= m * n 
W[4][4]   		= 4 * 4
F[10][4]  		= n * 4

ur[2][1]      	= m * 1
uG[2][1]  		= m * 1 

xr[10][1] 		= n * 1
xG[10][1] 		= n * 1 
v[10][1]  		= n * 1 
s[10][1] 		= n * 1 
x[10][1] 		= n * 1 
xe[10][1] 		= n * 1 
xp[10][1] 		= n * 1 

tau[10][1] 		= n * 1 
alarmSys[10][1] = n * 1 
*/
