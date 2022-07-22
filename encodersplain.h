/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: Header file of encodersplain.cpp (the main programming file of the SEAL crypto application)   	
	Note				: Commenting and cleaning up has been on July 2022
   Abbreviation/Acronym	: # _RPL: Raw plain 
						  # _mrp: Matrix row packing
*/

#ifndef ENCODERSPLAIN_H 
#define ENCODERSPLAIN_H

/* Import the other classes' header files */
#include "examples.h" 			  // The class for the console menu and guiding the user to the preferred application     	
#include "rawplain.h"  			  // The class which reads and stores the plain matrices                    
#include "generateplaintextMRP.h" // The functions used to prepare plain matrices in MRP format
#include "prepareVecMatMRP.h"     // The class which prepares the read matrices in MRP format  		
#include "initializationPLCP.h"   // The class containing the functions which convert the prepared matrices to Plaintext and Ciphertext objects.            
#include "plcpOperations.h"       // The class containing the functions which does Ciphertext-Ciphertext and Plaintext-Ciphertext arithmetic and algebraic operations  	    
#include "encryptedAppx.h"        // The class containing the functions which does the Chebyshev approximation with different assumptions     
#include "secretShare.h" 		  // The class containing the functions which does the secret sharing  
#include "applyPLCPSimulation.h"  // The class containing the crypto application functions for each targeted equation

/* Import the important selected C libraries */
#include <iostream>
#include <array>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>
#include <cmath>
#include <vector>

/* Call main namespaces */
using namespace std;
using namespace seal;


/* The function for performing the cyberphysical system's functionalities (one after another) */
void performMultipleMatrixVectorMultiplicationsPlain(EncryptionParameters *parmsPtr, 
SEALContext *contextPtr, auto *secret_keyPtr, PublicKey *public_keyPtr, RelinKeys *relin_keysPtr, 
GaloisKeys *gal_keysPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, 
FILE *fp,  double scale,  size_t numberOfIterations);

/* The function for setting the configurations of the SEAL crypto application and calling the main system function */
void ckks_encoder_modify_matrix_row_packing_functional();


#endif
