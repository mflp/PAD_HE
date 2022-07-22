/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which is used for applying the secret share operation in the row-matrix packing format  
	Note				: Commenting and cleaning up has been on July 2022 	
*/

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

/* Define the constant */
#define numRandBits 20

/* Call main namespaces */
using namespace std;
using namespace seal;

/*
	The function for performing the secret share in which summation for each index is applied during the secret sharing, which is not preferable (thus, omitted in the systen)   
*/
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
	int calculated_x_row_sizeVal){

		/* 
		**	Create addition, subtraction vectors and do the initial random addition
		*/
		// Create cleartext random vector for the eventual subtraction 
		vector<double> rand_double_sub(x_vector_repeat * calculated_x_row_sizeVal);		
		for (size_t  i = 0; i < x_vector_repeat * calculated_x_row_sizeVal; i++){				
			if(i < calculated_x_row_sizeVal){
				if(i < x_vector_col_size) {
					int randomNumSamplingInterval = pow(2, numRandBits);   // 4194304; // 134217728; // 4194304; // 500; // 4294967296;
					double ran 	 	= rand() % randomNumSamplingInterval;		
					rand_double_sub[i] 	= ran;				 
				}else{
					rand_double_sub[i] 	= 0;	
				}						
			}else{
				rand_double_sub[i] = rand_double_sub[i % calculated_x_row_sizeVal];			
			}	
		}	
		// Create cleartext random vector for the initial addition (same with the subtraction vector)
		vector<double> rand_double_add(x_vector_repeat * calculated_x_row_sizeVal);
		for (size_t  i = 0; i < x_vector_repeat * calculated_x_row_sizeVal; i++){				
				if(i % calculated_x_row_sizeVal == 0)
					rand_double_add[i] = rand_double_sub[i/calculated_x_row_sizeVal];
				else
					rand_double_add[i] = 0; 		
		}	
		// Encrypt the subtraction vector	
		Plaintext plainRandSub;
		Ciphertext randSubEnc;
		Ciphertext *randSubEncPtr;
		randSubEncPtr = &randSubEnc; 							
		encoderPtr->encode(rand_double_sub, scale, plainRandSub);
		encryptorPtr->encrypt(plainRandSub, randSubEnc);	
		// Encrypt the addition vector	
		Plaintext plainRandAdd;
		Ciphertext randAddEnc;
		Ciphertext *randAddEncPtr;
		randAddEncPtr = &randAddEnc; 							
		encoderPtr->encode(rand_double_add, scale, plainRandAdd);
		encryptorPtr->encrypt(plainRandAdd, randAddEnc);	
		// Add encrypted noise to the x vector.
		Ciphertext x_Random_Noise_Added;  
		Ciphertext *x_Random_Noise_AddedPtr = &x_Random_Noise_Added;
		// Do initial encrypted addition
		addSubtractTwoVector(x_Random_Noise_AddedPtr, x_vector_EncPtr, randAddEncPtr, scale, contextPtr, evaluatorPtr, true);

		/* 
		**	Decrpyt, obtain index sums, and add trailing zeros to each row segment 
		*/
		// Decrypt the noise added vector 							 			
		Plaintext plain_SecretShare_Decrypt_Result;	
		vector<double> decryptedVec;
		decryptorPtr->decrypt(x_Random_Noise_Added, plain_SecretShare_Decrypt_Result);
		encoderPtr->decode(plain_SecretShare_Decrypt_Result, decryptedVec);						
		// Obtain each x^e[k-1]'s index sum		
		vector<double> xIndSumDecVectors(x_vector_repeat * calculated_x_row_sizeVal);			
		for(size_t i = 0; i < x_vector_col_size; i++){
			double sumIthEntry = 0.0;
			for(size_t j = 0; j < calculated_x_row_sizeVal; j++)
				sumIthEntry += decryptedVec[i * calculated_x_row_sizeVal + j];			
			xIndSumDecVectors[i] = sumIthEntry;  
		}
		// Assign decrypted index sum and trailing zeros to each row segment (i.e., the format of "x1x2x3...x1x2x3...x1x2x3...x1x2x3...x1x2x3... ...")  	
		for(size_t i = x_vector_col_size; i < x_vector_repeat * calculated_x_row_sizeVal;  i++){
			if (i % calculated_x_row_sizeVal >= x_vector_col_size)
				xIndSumDecVectors[i] = 0; 
			else 
				xIndSumDecVectors[i] = xIndSumDecVectors[i % calculated_x_row_sizeVal]; 	
		}

		/* 
		**	Encrpyt the obtained vectors, do the eventual subtraction, and assign the fresh (or recrypted) vector 
		*/
		// Encrypt the index summed and trailing-zeros-added vector of x^e[k-1]
		Plaintext plainText; 
		Ciphertext freshCipherX;
		Ciphertext *freshCipherXPtr = &freshCipherX;
		encoderPtr->encode(xIndSumDecVectors, scale, plainText);
		encryptorPtr->encrypt(plainText, freshCipherX);
		// Do eventual encrypted subtraction from encrypted noise
		Ciphertext x_Random_Noise_Subtracted;  
		Ciphertext *x_Random_Noise_SubtractedPtr = &x_Random_Noise_Subtracted;
		addSubtractTwoVector(x_Random_Noise_SubtractedPtr, freshCipherXPtr, randSubEncPtr, scale, contextPtr, evaluatorPtr, false);
		// Assign the fresh (or recrypted) vector fresh encrypted vector	
		*x_vector_EncPtr = x_Random_Noise_Subtracted;
}


/*
	The function for performing the secret share in which no index summation is applied during the secret sharing, which is the preferable one (thus, chosen in the systen)   
*/
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
	int calculated_x_row_sizeVal){

		/* 
		**	Create addition, subtraction vectors and make the initial random addition
		*/
		// Create cleartext random vector for the initial addition
		vector<double> rand_double_add(x_vector_repeat * calculated_x_row_sizeVal);		
		for (size_t  i = 0; i < x_vector_repeat * calculated_x_row_sizeVal; i++){	
			int randomNumSamplingInterval = pow(2, numRandBits);   // 4194304; // 134217728; // 4194304; // 500; // 4294967296;
			double ran 	 	= rand() % randomNumSamplingInterval;
			rand_double_add[i] = ran;					
		}
		// Create cleartext random vector for the eventual subtraction (same with the addition vector)
		vector<double> rand_double_sub(x_vector_repeat * calculated_x_row_sizeVal);
		for (size_t  i = 0; i < x_vector_repeat * calculated_x_row_sizeVal; i++)				
			rand_double_sub[i] = rand_double_add[i]; 		
		// Create encrypted random eventual subtraction vector	
		Plaintext plainRandSub;
		Ciphertext randSubEnc;
		Ciphertext *randSubEncPtr;
		randSubEncPtr = &randSubEnc; 							
		encoderPtr->encode(rand_double_sub, scale, plainRandSub);
		encryptorPtr->encrypt(plainRandSub, randSubEnc);
		// Create encrypted random initial addition vector	
		Plaintext plainRandAdd;
		Ciphertext randAddEnc;
		Ciphertext *randAddEncPtr;
		randAddEncPtr = &randAddEnc; 							
		encoderPtr->encode(rand_double_add, scale, plainRandAdd);
		encryptorPtr->encrypt(plainRandAdd, randAddEnc);
		// Add encrypted noise to the x vector.
		Ciphertext x_Random_Noise_Added;  
		Ciphertext *x_Random_Noise_AddedPtr = &x_Random_Noise_Added;
		// Do initial encrypted addition
		addSubtractTwoVector(x_Random_Noise_AddedPtr, x_vector_EncPtr, randAddEncPtr, scale, contextPtr, evaluatorPtr, true);

		/* 
		**	Decrpyt, and add trailing zeros to each row segment 
		*/
		// Decrypt the noise added vector 							 			
		Plaintext plain_SecretShare_Decrypt_Result;	
		vector<double> decryptedVec;
		// cout << "Decrypt and decode the secret shared vector " << endl;
		decryptorPtr->decrypt(x_Random_Noise_Added, plain_SecretShare_Decrypt_Result);
		encoderPtr->decode(plain_SecretShare_Decrypt_Result, decryptedVec);				
		// Encrypt the index summed decryption vector of x^e[k-1]
		Plaintext plainText; 
		Ciphertext freshCipherX;
		Ciphertext *freshCipherXPtr = &freshCipherX;
		encoderPtr->encode(decryptedVec, scale, plainText);
		encryptorPtr->encrypt(plainText, freshCipherX);
		
		/* 
		**	Do the eventual subtraction, and assign the fresh (or recrypted) vector 
		*/	
		// Do eventual encrypted subtraction from encrypted noise
		Ciphertext x_Random_Noise_Subtracted;  
		Ciphertext *x_Random_Noise_SubtractedPtr = &x_Random_Noise_Subtracted;
		addSubtractTwoVector(x_Random_Noise_SubtractedPtr, freshCipherXPtr, randSubEncPtr, scale, contextPtr, evaluatorPtr, false);
		// Assign fresh encrypted vector	
		*x_vector_EncPtr = x_Random_Noise_Subtracted;
}

/*
	The function for performing the secret share specifically for the CUSUM Statistics' Computation   
*/
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
	int calculated_x_row_sizeVal){

		/* 
		**	Create addition, subtraction vectors and do the initial random addition
		*/
		// Create cleartext random Noise vector for the initial addition
		vector<double> rand_double_add(x_vector_repeat * calculated_x_row_sizeVal);		
		for (size_t  i = 0; i < x_vector_repeat * calculated_x_row_sizeVal; i++){				
			int randomNumSamplingInterval = pow(2, numRandBits);   // 4194304; // 134217728; // 4194304; // 500; // 4294967296;
			double ran = rand() % randomNumSamplingInterval;		
			rand_double_add[i] = ran;					
		}
		// Create cleartext random Noise vector for the eventual subtraction (same with the addition vector)
		vector<double> rand_double_sub(x_vector_repeat * calculated_x_row_sizeVal);
		for (size_t  i = 0; i < x_vector_repeat * calculated_x_row_sizeVal; i++)				
			rand_double_sub[i] = rand_double_add[i]; 		
		// Create encrypted random eventual subtraction vector	
		Plaintext plainRandSub;
		Ciphertext randSubEnc;
		Ciphertext *randSubEncPtr;
		randSubEncPtr = &randSubEnc; 							
		encoderPtr->encode(rand_double_sub, scale, plainRandSub);
		encryptorPtr->encrypt(plainRandSub, randSubEnc);
		// Create encrypted random initial addition vector	
		Plaintext plainRandAdd;
		Ciphertext randAddEnc;
		Ciphertext *randAddEncPtr;
		randAddEncPtr = &randAddEnc; 							
		encoderPtr->encode(rand_double_add, scale, plainRandAdd);
		encryptorPtr->encrypt(plainRandAdd, randAddEnc);
		// Add encrypted noise to the x vector.
		Ciphertext x_Random_Noise_Added;  
		Ciphertext *x_Random_Noise_AddedPtr = &x_Random_Noise_Added;
		// Do initial encrypted addition
		addSubtractTwoVector(x_Random_Noise_AddedPtr, x_vector_EncPtr, randAddEncPtr, scale, contextPtr, evaluatorPtr, true);
		
		/* 
		**	Decrpyt, do the eventual subtraction, and assign the fresh (or recrypted) vector 
		*/	
		// Decrypt the noise added vector 							 			
		Plaintext plain_SecretShare_Decrypt_Result;	
		vector<double> decryptedVec;
		decryptorPtr->decrypt(x_Random_Noise_Added, plain_SecretShare_Decrypt_Result);
		encoderPtr->decode(plain_SecretShare_Decrypt_Result, decryptedVec);				
		// Encrypt the already-index-summed decryption vector
		Plaintext plainText; 
		Ciphertext freshCipherX;
		Ciphertext *freshCipherXPtr = &freshCipherX;
		encoderPtr->encode(decryptedVec, scale, plainText);
		encryptorPtr->encrypt(plainText, freshCipherX);
		// Do eventual encrypted subtraction from encrypted noise
		Ciphertext x_Random_Noise_Subtracted;  
		Ciphertext *x_Random_Noise_SubtractedPtr = &x_Random_Noise_Subtracted;
		addSubtractTwoVector(x_Random_Noise_SubtractedPtr, freshCipherXPtr, randSubEncPtr, scale, contextPtr, evaluatorPtr, false);
		// Assign fresh encrypted vector	
		*x_vector_EncPtr = x_Random_Noise_Subtracted;
}

/*
	The function for performing the secret share specifically for the Estimation Computation   
*/
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
	int calculated_x_row_sizeVal){

		/* 
		**	Create addition, subtraction vectors and do the initial random addition
		*/
		// Create cleartext random vector for the eventual subtraction (w.r.t. repeated-format (Evr) addition) 
		vector<double> rand_double_sub(x_vector_repeat * calculated_x_row_sizeVal);		
		for (size_t  i = 0; i < x_vector_repeat * calculated_x_row_sizeVal; i++){				
			if(i < calculated_x_row_sizeVal){
				if(i < x_vector_col_size) {
					int randomNumSamplingInterval = pow(2, numRandBits);   // 4194304; // 134217728; // 4194304; // 500; // 4294967296;
					double ran = rand() % randomNumSamplingInterval;		
					rand_double_sub[i] 	= ran;				 
				}else{
					rand_double_sub[i] 	= 0;	
				}						
			}else{
				rand_double_sub[i] = rand_double_sub[i % calculated_x_row_sizeVal];			
			}	
		}	
		// Create cleartext random vector for the initial addition (w.r.t. not-repeated format (Evo) addition)
		vector<double> rand_double_add(x_vector_repeat * calculated_x_row_sizeVal);
		for (size_t  i = 0; i < x_vector_repeat * calculated_x_row_sizeVal; i++){				
				if(i % calculated_x_row_sizeVal == 0)
					rand_double_add[i] = rand_double_sub[i/calculated_x_row_sizeVal];
				else
					rand_double_add[i] = 0; 		
		}	
		// Create encrypted random eventual subtraction vector	
		Plaintext plainRandSub;
		Ciphertext randSubEnc;
		Ciphertext *randSubEncPtr;
		randSubEncPtr = &randSubEnc; 							
		encoderPtr->encode(rand_double_sub, scale, plainRandSub);
		encryptorPtr->encrypt(plainRandSub, randSubEnc);
		// Create encrypted random initial addition vector	
		Plaintext plainRandAdd;
		Ciphertext randAddEnc;
		Ciphertext *randAddEncPtr;
		randAddEncPtr = &randAddEnc; 							
		encoderPtr->encode(rand_double_add, scale, plainRandAdd);
		encryptorPtr->encrypt(plainRandAdd, randAddEnc);
		// Add encrypted noise to the x vector.
		Ciphertext x_Random_Noise_Added;  
		Ciphertext *x_Random_Noise_AddedPtr = &x_Random_Noise_Added;
		// Do initial encrypted addition
		addSubtractTwoVector(x_Random_Noise_AddedPtr, x_vector_EncPtr, randAddEncPtr, scale, contextPtr, evaluatorPtr, true);

		/* 
		**	Decrpyt, re-arrange the estimation vector, add trailing zeros to each row segment, and encrypt the already-index-summed vector 
		*/	
		// Decrypt the noise added vector 							 			
		Plaintext plain_SecretShare_Decrypt_Result;	
		vector<double> decryptedVec;
		decryptorPtr->decrypt(x_Random_Noise_Added, plain_SecretShare_Decrypt_Result);
		encoderPtr->decode(plain_SecretShare_Decrypt_Result, decryptedVec);				
		// Re-arrange the estimation vector and add trailing zeros to each row segment	 			
		vector<double> xE_Rearranged(x_vector_repeat * calculated_x_row_sizeVal);		
		// Obtain each index sum		
		for(size_t i = 0; i < x_vector_repeat * calculated_x_row_sizeVal; i++){ 
			if(i < calculated_x_row_sizeVal){
				if(i < x_vector_col_size)
					xE_Rearranged[i] = decryptedVec[i * calculated_x_row_sizeVal]; 
				else
					xE_Rearranged[i] = 0;	
			}else{
				xE_Rearranged[i] = xE_Rearranged[i % calculated_x_row_sizeVal];			
			}			
		}
		// Encrypt the already-index-summed vector
		Plaintext plainText; 
		Ciphertext freshCipherX;
		Ciphertext *freshCipherXPtr = &freshCipherX;
		encoderPtr->encode(xE_Rearranged, scale, plainText);
		encryptorPtr->encrypt(plainText, freshCipherX);

		/* 
		**	Do the eventual subtraction, and assign the fresh (or recrypted) vector 
		*/	
		// Do eventual encrypted subtraction from encrypted noise
		Ciphertext x_Random_Noise_Subtracted;  
		Ciphertext *x_Random_Noise_SubtractedPtr = &x_Random_Noise_Subtracted;
		addSubtractTwoVector(x_Random_Noise_SubtractedPtr, freshCipherXPtr, randSubEncPtr, scale, contextPtr, evaluatorPtr, false);
		// Assign fresh encrypted vector	
		*x_vector_EncPtr = x_Random_Noise_Subtracted;
}
