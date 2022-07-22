/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: Header file of the class to handle the fundamental homomorphic operations (i.e., plcpOperations.cpp) 
	Note				: Commenting and cleaning up has been on July 2022
*/

#ifndef PLCPOPERATIONS_H
#define PLCPOPERATIONS_H

/* Import the other classes' header files*/
#include "examples.h" 			  // The class for the console menu and guiding the user to the preferred application     	
#include "rawplain.h"  			  // The class which reads and stores the plain matrices                    
#include "initializationPLCP.h"   // The class containing the functions which convert the prepared matrices to Plaintext and Ciphertext objects.            

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
#include <unistd.h>
#include <time.h>

/* Call the important namespaces */
using namespace std;
using namespace seal;


/* The function for finding the lowest number of slots equaling to a binary power for the matrix row packing format (one after another)  */
int findBinaryPower(int indexVectorLength);

// ##################################################################################
// ##################################################################################
// ############################PLAINTEXT-CIPHERTEXT##################################
// ##################################################################################
// ##################################################################################

/* The function for equalizing the chain indexes and scales of a plaintext and a ciphertext to make a homomorphic operation (e.g., Addition, multiplication, etc.) */
void adjustScaleandChainParametersPlaintextAndVectors(Ciphertext *CiphertextPtr, Plaintext *plaintextPtr, double scale, SEALContext *contextPtr, Evaluator *evaluatorPtr);

/* The function for performing a not-in-place homomorphic multiplication with Plaintext&Ciphertext */
void plnCprAxMult_mrp(Ciphertext *res, Ciphertext *cpr, Plaintext *pln, double scale, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, SEALContext *contextPtr, Evaluator *evaluatorPtr, RelinKeys *relin_keysPtr);

/* The function for performing an in-place homomorphic multiplication with Plaintext&Ciphertext */
void plnCprAxMult_InPlac_mrp(
	Ciphertext *res, 
	Plaintext *pln,
	double scale, 
	Decryptor *decryptorPtr,
	CKKSEncoder *encoderPtr,
	SEALContext *contextPtr,
	Evaluator *evaluatorPtr,
	RelinKeys *relin_keysPtr);
/* The function for performing an not-in-place homomorphic addition or subtraction with Plaintext&Ciphertext */
void addSubtractPLCPVector(
	Ciphertext *matrixVecAddRes, 
	Ciphertext *v_1_Ciph,
	Plaintext  *v_2_Pln,
	double scale, 	
	SEALContext *contextPtr,	
	Evaluator *evaluatorPtr, 
	bool isAddition);


// ##################################################################################
// ##################################################################################
// ############################CIPHERTEXT-CIPHERTEXT#################################
// ##################################################################################
// ##################################################################################
/* The function for equalizing the chain indexes and scales of two ciphertexts to make a homomorphic operation (e.g., Addition, multiplication, etc.) */
void adjustScaleandChainParametersTwoVectors( 
	Ciphertext *Cipher1,
	Ciphertext *Cipher2,
	Ciphertext *CipherRes,	
	double scale, 
	SEALContext *contextPtr,	
	Evaluator *evaluatorPtr);

/* The function for performing the rotation through a Ciphertext to sum the first n entries into the 1st index */
void rotateVector(
	Ciphertext *sum_Output, 
	Ciphertext *v_to_Rotate,
	double scale, 	
	int indexVectorLength, 
	SEALContext *contextPtr,		
	Evaluator *evaluatorPtr,
	Decryptor *decryptorPtr, 
	CKKSEncoder *encoderPtr,
	GaloisKeys *gal_keysPtr, 
	RelinKeys *relin_keysPtr);

/* The function for performing a not-in-place homomorphic multiplication with Ciphertext&Ciphertext */
void matrixVectorMultMatRowPacking(
	Ciphertext *matrixVecMultRes, 
	Ciphertext *matrixCipher,
	Ciphertext *vectorCipher,
	double scale, 
	SEALContext *contextPtr,
	Evaluator *evaluatorPtr,
	RelinKeys *relin_keysPtr);

/* The function for performing a in-place homomorphic squaring operation with a ciphertext */
void vectorSquaringRowPacking(
	Ciphertext *vectorCipherSquareRes,	
	Ciphertext *vectorCipher,
	double scale, 
	SEALContext *contextPtr,
	Evaluator *evaluatorPtr,
	RelinKeys *relin_keysPtr);

/* The function for performing a homomorphic addition of three ciphertexts */
void addThreeVector(Ciphertext *matrixVecAddRes, 
	Ciphertext *v_1_Ciph,
	Ciphertext *v_2_Ciph,
	Ciphertext *v_3_Ciph,
	double scale, 	
	SEALContext *contextPtr,	
	Evaluator *evaluatorPtr,
	RelinKeys *relin_keysPtr);

/* The function for performing an not-in-place homomorphic addition or subtraction with Ciphertext&Ciphertext */
void addSubtractTwoVector(Ciphertext *matrixVecAddRes, 
	Ciphertext *v_1_Ciph,
	Ciphertext *v_2_Ciph,
	double scale, 	
	SEALContext *contextPtr,	
	Evaluator *evaluatorPtr, 
	bool isAddition);

# endif
