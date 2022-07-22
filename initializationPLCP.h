/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: The header file of the class to initialize the simulation data (i.e., initializationPLCP.cpp) 
	Note				: Commenting and cleaning up has been on July 2022
*/

#ifndef INITIALIZATIONPLCP_H
#define INITIALIZATIONPLCP_H

/* Import the other classes' header files */
#include "examples.h"
#include "rawplain.h"             // 1   - The class which reads and stores the plain matrices        
#include "generateplaintextMRP.h" // 2.0 - The functions used to prepare plain matrices in MRP format
#include "prepareVecMatMRP.h"     // 2.1 - The class which prepares the read matrices in MRP format  		
#include "initializationPLCP.h"   // 3   - The class containing the functions which convert the prepared matrices to Plaintext and Ciphertext objects.            
#include "plcpOperations.h"       // 4.0 - The class containing the functions which does Ciphertext-Ciphertext and Plaintext-Ciphertext arithmetic and algebraic operations  	    
#include "encryptedAppx.h"        // 4.1 - The class containing the functions which does the Chebyshev approximation with different assumptions     
#include "secretShare.h" 		  // 4.2 - The class containing the functions which does the secret sharing  
#include "applyPLCPSimulation.h"  // 6   - The class containing the crypto application functions for each targeted equation

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

/* Call the important namespaces */
using namespace std;
using namespace seal;

/* Data struct created for the Plaintext-Ciphertext setting, which includes only the data for Equation-2   */
struct simulationMatrixPLCP{
	// Variables	
	int m;
	int tMax; 
	int n;

	// Pointers	
	seal:: Plaintext * xexe_PL;  // 2nd Equation-Part1
	seal:: Plaintext * GAMMA_PL; // 2nd Equation-Part1   		
	seal:: Plaintext * LL_PL;    // 2nd Equation-Part2
	seal:: Ciphertext * yy_CP;   // 2nd Equation-Part2
	seal:: Ciphertext * xexe_CP; // 2nd Equation-Part1	
	seal:: Plaintext * xGxG_PL;  // 2nd Equation-Part3
};

/* The function for preparing the power series coefficient vector in plaintext */
void preparePwSrCoeffVec(vector<Plaintext> * powerSeriesCoeffVec_PL, struct simulationMatrixData *smd, size_t chebyshevDegree, size_t calcRowSize, double scale, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, bool isMaxFunc);

/* The function for initializing the simulation data solely for the Equation-2 in the Plaintext-Ciphertext setting */
void initSimDatPLCP(struct simulationMatrixPLCP *splcp, struct simulationMatrixMRP *smrp, double scale, Decryptor *decryptorPtr, Encryptor *encryptorPtr, CKKSEncoder *encoderPtr);

/* The function for encoding a vector into Plaintext */
void makePlaintextMatRowPacking(std::vector<double> *x_vector, double scale, seal::Plaintext *plain_xePtr, seal::CKKSEncoder *encoderPtr);

/* The function for encrypting a vector into Ciphertext */
void encryptXVectorMatRowPacking(vector<double> *x_vector,double scale, Ciphertext *x_vector_EncPtr, Encryptor *encryptorPtr, CKKSEncoder *encoderPtr);

/* The function for encrypting a matrix into Ciphertext */
void encryptMatrixMatRowPacking(vector<double> *matPtr, double scale, Ciphertext *mat_EncPtr, Encryptor *encryptorPtr, CKKSEncoder *encoderPtr);

#endif
