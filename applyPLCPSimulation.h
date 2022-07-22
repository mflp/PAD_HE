/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: The header class of the the main operational functions(i.e., applyPLCPSimulation.cpp)  
	Note				: Commenting and cleaning up has been on July 2022
*/

#ifndef APPLYPLCPSIMULATION_H
#define APPLYPLCPSIMULATION_H

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
#include "applyPLCPSimulation.h"  // The class containing the crypto application functions for each targeted equation

/* Import the important selected C libraries*/
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* Call main namespaces */
using namespace std;
using namespace seal;

/* The function for performing the cyberphysical system's sensor measurement functionality   */
void sense_Encrypt_y(struct simulationMatrixData *smd, struct simulationMatrixMRP *smrp, int numiter, Ciphertext *yy_CP, Ciphertext *yyAS_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, CKKSEncoder *encoderPtr);

/* The function for performing the cyberphysical system's estimation functionality   */
void applyEquation_2_PLCP(Ciphertext * secEqRes_CP, struct simulationMatrixMRP *smrp, Plaintext *GAMMA_PL, Plaintext *LL_PL, Ciphertext * xGxG_CP,  Ciphertext * yy_CP, Ciphertext * xexe_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr);

/* The function for performing the cyberphysical system's control action functionality   */
void applyEquation_3_PLCP(Ciphertext * thirdEqRes_CP, struct simulationMatrixMRP *smrp, Plaintext *KGKG_PL, Plaintext *KLKL_PL, Ciphertext * KxugKxug_CP,  Ciphertext * yy_CP, Ciphertext * xexe_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr);

/* The function for performing the cyberphysical system's control action functionality at the very first iteration   */
void applyEquation_3_fiter_PLCP(Ciphertext * thirdEqRes_CP, struct simulationMatrixMRP *smrp, Plaintext *KxKx_PL, Ciphertext * uGuG_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr);

/* The function for performing the cyberphysical system's prediction functionality   */
void applyEquation_4_5_PLCP(Ciphertext * fourthfifthEqRes_CP, struct simulationMatrixMRP *smrp, Plaintext *ACL_PL, Plaintext *BB_PL, Ciphertext * uGuG_CP, Ciphertext * xexe_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr);

/* The function for performing the cyberphysical system's residues computation functionality */
void applyEquation_6_PLCP(Ciphertext * sixthEqRes_CP, struct simulationMatrixMRP * smrp, Ciphertext * fourthfifthEqRes_CP, Ciphertext * yyAS_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr);

/* The function for performing the cyberphysical system's residues computation functionality at the very first iteration */
void applyEquation_6_fiter_PLCP(Ciphertext * sixthEqRes_CP, struct simulationMatrixMRP * smrp, Plaintext * xpxp_PL, Ciphertext * yyAS_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr);

/* The function for performing the cyberphysical system's 1st part of CUSUM Computation and Alarm Computation */
void applyEquation_CUSUM_PLCP(Ciphertext * eigthEqRes_CP, Ciphertext * ninthEqRes_CP, struct simulationMatrixMRP * smrp, struct simulationMatrixData * smd, int numiter,
							Plaintext * ss_PL, Plaintext * vv_PL, Plaintext * TAU_PL, Plaintext * alpha_Eq8_PL, Plaintext * beta_Eq8_PL, 
							Plaintext * alpha_Eq9_PL, Plaintext * beta_Eq9_PL, 
							Plaintext *vectorOnePtr_PL,  
							Plaintext *firstPowerSeriesTermEq8Ptr_PL, vector<Plaintext> * powSerCoeffArr_Eq8_PL,
							Plaintext *firstPowerSeriesTermEq9Ptr_PL, vector <Plaintext> * powSerCoeffArr_Eq9_PL,   
							Ciphertext * sixthEqRes_CP, Ciphertext * ss_CP, 
							double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, 
							Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr, bool isFirstIter);

/* The function for performing the cyberphysical system's 2nd and last part of CUSUM Computation */
void applyEquation_10_PLCP(Ciphertext * tenthEqRes_CP, struct simulationMatrixMRP * smrp, Ciphertext * eigthEqRes_CP, Ciphertext * ninthEqRes_CP, Ciphertext *vectorOnePtr_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr);

/* The function for performing the cyberphysical system's sensor measurement simulation in a naturally noisy environment */
void applyXVecNoiseAddition(Ciphertext * controlAction_CP, int numIter, struct  simulationMatrixMRP * smrp, struct  simulationMatrixData * smd, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr);

/* The function for extracting the experimental results */
void extractExpRes(Ciphertext * res_CP, struct simulationMatrixMRP * smrp, struct simulationMatrixData * smd, int numIter, int numOfRows, double scale, SEALContext *contextPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, char * matname);

#endif
