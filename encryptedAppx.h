/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: The header file for the class used for the encrypted Chebyshev approximation (encryptedAppx.cpp) 
	Note				: Commenting and cleaning up has been on July 2022
*/

#ifndef ENCRYPTEDAPPX_H
#define ENCRYPTEDAPPX_H

/* Import the other classes' header files */
#include "examples.h" 			  // The class for the console menu and guiding the user to the preferred application     	
#include "rawplain.h"  			  // The class which reads and stores the plain matrices                    
#include "initializationPLCP.h"   // The class containing the functions which convert the prepared matrices to Plaintext and Ciphertext objects.                
#include "encryptedAppx.h"        // The class containing the functions which does the Chebyshev approximation with different assumptions     

/* Import the important selected C libraries */
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <time.h>
#include <stdbool.h>

/* Call main namespaces */
using namespace std;
using namespace seal;

/* The function for performing the ciphertext-plaintext-mixed Chebyshev Approximation */
void applyChebyshevPolynomials(
	Ciphertext *sumOutputPtr, 
	Ciphertext *vectorNotKnownPosNegPtr, 
	Ciphertext *vectorOnePtr, 
	Ciphertext *firsCoeffArrPtr,
	double *coeffArr, 
	int numCoeff, 
	double scale, 
	SEALContext *contextPtr, 
	Evaluator *evaluatorPtr, 
	CKKSEncoder *encoderPtr, 	
	GaloisKeys *gal_keysPtr,
	RelinKeys *relin_keysPtr, 
	Decryptor *decryptorPtr);

/* The function for performing the multiple polynomials addition in the ciphertext-plaintext-mixed setting (actually all the terms except the first term) */
void addManyVectorsInaRow(Ciphertext *sumOutputPtr, vector <Ciphertext> * chebVectorPtr, int exactNumCoeff, double scale, SEALContext *contextPtr, Evaluator *evaluatorPtr, CKKSEncoder *encoderPtr, GaloisKeys *gal_keysPtr, RelinKeys *relin_keysPtr);

/* The function for performing the ciphertext-only Chebyshev Approximation */	
void makeChebyshevPolynAppxPLCP(
	struct simulationMatrixMRP *smrp,
	Ciphertext * sumOutputPtr, 
	Ciphertext * vecTobeAppx_CP, 
	Plaintext  * vectorOnePtr_PL, 
	Plaintext  * firstPowerSeriesTermPtr_PL,
	vector<Plaintext> * powSerCoeffArr_PL, 
	int numCoeff, 
	double scale, 
	SEALContext *contextPtr, 
	Evaluator *evaluatorPtr, 
	CKKSEncoder *encoderPtr, 	
	GaloisKeys *gal_keysPtr,
	RelinKeys *relin_keysPtr, 
	Decryptor *decryptorPtr);	

/* The function for performing the multiple polynomials addition in the ciphertext-only setting (actually all the terms except the first term)  */	
void addManyVectorsInaRow_PL(Ciphertext *sumOutputPtr, vector <Ciphertext> * chebVectorPtr, int exactNumCoeff, double scale, SEALContext *contextPtr, Evaluator *evaluatorPtr, CKKSEncoder *encoderPtr, GaloisKeys *gal_keysPtr, RelinKeys *relin_keysPtr);

#endif

