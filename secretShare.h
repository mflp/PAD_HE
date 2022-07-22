/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: The header file for the class used for applying the secret share operation (i.e., secretShare.cpp)   	
	Note				: Commenting and cleaning up has been on July 2022
*/

#ifndef SECRETSHARE_H
#define SECRETSHARE_H

/* Import the other classes' header files*/
#include "examples.h" 			  // The class for the console menu and guiding the user to the preferred application     	
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

/* The function for performing the secret share in which summation for each index is applied during the secret sharing, which is not preferable (thus, omitted in the systen) */
void secretShare(
	Ciphertext *x_vector_EncPtr,
	SEALContext *contextPtr,  
	Encryptor *encryptorPtr, 
	Evaluator *evaluatorPtr, 
	Decryptor *decryptorPtr, 
	CKKSEncoder *encoderPtr, 
	double scale, 	
	int x_vector_repeat,
	int x_vector_col_size,
	int calculated_x_row_sizeVal);

/* The function for performing the secret share in which no index summation is applied during the secret sharing, which is the preferable one (thus, chosen in the systen) */
void secretSharev2(
	Ciphertext *x_vector_EncPtr,
	SEALContext *contextPtr,  
	Encryptor *encryptorPtr, 
	Evaluator *evaluatorPtr, 
	Decryptor *decryptorPtr, 
	CKKSEncoder *encoderPtr, 
	double scale, 	
	int x_vector_repeat,
	int x_vector_col_size,
	int calculated_x_row_sizeVal);

/* The function for performing the secret share specifically for the CUSUM Statistics' Computation */
void secretShareCUSUMParamSum(
	Ciphertext *x_vector_EncPtr,
	SEALContext *contextPtr,  
	Encryptor *encryptorPtr, 
	Evaluator *evaluatorPtr, 
	Decryptor *decryptorPtr, 
	CKKSEncoder *encoderPtr, 
	double scale, 	
	int x_vector_repeat,
	int x_vector_col_size,
	int calculated_x_row_sizeVal);

/* The function for performing the secret share specifically for the Estimation Computation */
void secretShareEstimation(
	Ciphertext *x_vector_EncPtr,
	SEALContext *contextPtr,  
	Encryptor *encryptorPtr, 
	Evaluator *evaluatorPtr, 
	Decryptor *decryptorPtr, 
	CKKSEncoder *encoderPtr, 
	double scale, 	
	int x_vector_repeat,
	int x_vector_col_size,
	int calculated_x_row_sizeVal);

#endif
