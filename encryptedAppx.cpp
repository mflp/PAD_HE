/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which is used for applying the encrypted approximation via the use of the Chebyshev Polynomials appx.
   Note					: There are two implementations which are very close to each other, 
						  One of them is ciphertext-plaintext-mixed Chebyshev Approximation, and the other is ciphertext-only Chebyshev Approximation. 
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
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <time.h>
#include <stdbool.h>

/* Call main namespaces */
using namespace std;
using namespace seal;

/*
	The function for performing the ciphertext-plaintext-mixed Chebyshev Approximation   
*/
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
	Decryptor *decryptorPtr){

	// Create Chebyshev Polynomials Array (parallel to the cleartext implementation of the Chebyshev Approximation)		  	
	vector<Ciphertext> chebyshevPolynomials(numCoeff + 1);
	chebyshevPolynomials[1] = *vecTobeAppx_CP; // First index becomes zero
	
	// Evaluate the Chebyshev Polynomials for the given u 
	for(size_t i = 2; i < numCoeff + 1; i++){
		if(i % 2 == 0){ // Once the index is even 
			// Create the ciphertext pointer for the targeted chebyshev polynomial   
			Ciphertext *target_ChebPolyPtr 				 = &chebyshevPolynomials[i];
			Ciphertext *target_ChebPoly_MultiplicandPtr1 = &chebyshevPolynomials[i/2];
			Ciphertext *target_ChebPoly_MultiplicandPtr2 = &chebyshevPolynomials[i/2]; 	
			// Multiply T[n] * T[n]
			matrixVectorMultMatRowPacking(target_ChebPolyPtr, target_ChebPoly_MultiplicandPtr1, target_ChebPoly_MultiplicandPtr2, scale, contextPtr, evaluatorPtr, relin_keysPtr);
			// Obtain  2 * T[n] * T[n] by a single addition 						
			Ciphertext twoTimesResChebyPoly;   
			Ciphertext *twoTimesResChebyPolyPtr = &twoTimesResChebyPoly;
			addSubtractTwoVector(twoTimesResChebyPolyPtr, target_ChebPolyPtr, target_ChebPolyPtr, scale, contextPtr, evaluatorPtr, true);			
			// Make the last subtraction 
			Ciphertext resChebyPolynomial;  
			Ciphertext *resChebyPolynomialPtr = &resChebyPolynomial;					
			addSubtractPLCPVector(resChebyPolynomialPtr, twoTimesResChebyPolyPtr, vectorOnePtr_PL, scale, 	contextPtr,	evaluatorPtr, false); 
			chebyshevPolynomials[i] = resChebyPolynomial; 
		}else{// Once the index is even 	 	
			// Create the ciphertext pointer for the targeted chebyshev polynomial
			Ciphertext *target_ChebPolyPtr 					= &chebyshevPolynomials[i];
			Ciphertext *target_ChebPoly_Multiplicand1Ptr 	= &chebyshevPolynomials[i/2];				
			Ciphertext *target_ChebPoly_Multiplicand2Ptr 	= &chebyshevPolynomials[i/2 + 1];
		   	// Multiply T[n] * T[n + 1]
			matrixVectorMultMatRowPacking(target_ChebPolyPtr, target_ChebPoly_Multiplicand1Ptr, target_ChebPoly_Multiplicand2Ptr, scale, contextPtr, evaluatorPtr, relin_keysPtr);
			// Obtain 2 * T[n] * T[n + 1] by a single addition
			Ciphertext twoTimesResChebyPoly;   
			Ciphertext *twoTimesResChebyPolyPtr = &twoTimesResChebyPoly;
			addSubtractTwoVector(twoTimesResChebyPolyPtr, target_ChebPolyPtr, target_ChebPolyPtr, scale, contextPtr, evaluatorPtr, true);
			// Make the last subtraction 
			Ciphertext resChebyPolynomial;  
			Ciphertext *resChebyPolynomialPtr = &resChebyPolynomial;
			addSubtractTwoVector(resChebyPolynomialPtr, twoTimesResChebyPolyPtr, vecTobeAppx_CP, scale, contextPtr, evaluatorPtr, false); 
			chebyshevPolynomials[i] = resChebyPolynomial;
		}																																																	
	}

	// Multiply the obtained Chebyshev Polynomials with the coeffients of the polynomials 
	for(size_t i = 1; i < numCoeff + 1; i++){ 
		Ciphertext * respChebTermPtr = & chebyshevPolynomials[i]; 
		plnCprAxMult_InPlac_mrp(respChebTermPtr, &powSerCoeffArr_PL->at(i), scale, decryptorPtr, encoderPtr, contextPtr, evaluatorPtr, relin_keysPtr);
	}
		
	// Add each product on top of each other: Call the add_many function implemented in a customized way 	
	vector <Ciphertext> * chebVectorPtr = & chebyshevPolynomials;  
	Ciphertext sumWoutFirstTerms;	
	Ciphertext *sumWoutFirstTermsPtr = &sumWoutFirstTerms; 				
	addManyVectorsInaRow_PL(sumWoutFirstTermsPtr, chebVectorPtr, numCoeff, scale, contextPtr, evaluatorPtr, encoderPtr, gal_keysPtr, relin_keysPtr);		
				
	// Do last plaintext-ciphertext addition to add the very first term of the Chebyshev power series 
	addSubtractPLCPVector(sumOutputPtr, sumWoutFirstTermsPtr, firstPowerSeriesTermPtr_PL, scale, contextPtr,	evaluatorPtr, true);	
}

/*
	The function for performing the multiple polynomials addition in the ciphertext-plaintext-mixed setting (actually all the terms except the first term)   
*/
void addManyVectorsInaRow_PL(Ciphertext *sumOutputPtr, vector <Ciphertext> * chebVectorPtr, int exactNumCoeff, double scale, SEALContext *contextPtr, Evaluator *evaluatorPtr, CKKSEncoder *encoderPtr, GaloisKeys *gal_keysPtr, RelinKeys *relin_keysPtr){
		
	// Find the minimum chain index
	int min_chain_index = 999;
	parms_id_type min_parms_id;
	for(size_t i = 1; i < exactNumCoeff + 1; i++){
		int chain_index = (*contextPtr->get_context_data(chebVectorPtr->at(i).parms_id())).chain_index();
		if(chain_index < min_chain_index){
			min_chain_index = chain_index;
			min_parms_id 	= chebVectorPtr->at(i).parms_id();
		} 
	}	

	// Adjust each chain index to the minimum chain index and scale to the already-setted scale, so that they could be added on top of another 
	for(size_t i = 1; i < exactNumCoeff + 1; i++){
    	evaluatorPtr->mod_switch_to_inplace(chebVectorPtr->at(i), min_parms_id);	
		chebVectorPtr->at(i).scale() = scale;	
	}
	
	// Add each Chebyshev polynomial on top of another     
	Ciphertext *ChebyshevPtr1 = &chebVectorPtr->at(1);
	Ciphertext *ChebyshevPtr2 = &chebVectorPtr->at(2);
    addSubtractTwoVector(sumOutputPtr, ChebyshevPtr1, ChebyshevPtr2, scale, contextPtr, evaluatorPtr, true);
	for(size_t i = 3; i < exactNumCoeff; i++)
		evaluatorPtr->add_inplace(*sumOutputPtr, chebVectorPtr->at(i));
}

/*
	The function for performing the ciphertext-only Chebyshev Approximation   
*/
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
	Decryptor *decryptorPtr){

	// Create Chebyshev Polynomials Array (parallel to the cleartext implementation of the Chebyshev Approximation)		  	
	vector<Ciphertext> chebyshevPolynomials(numCoeff + 1);
	chebyshevPolynomials[0] = *vectorOnePtr;
	chebyshevPolynomials[1] = *vectorNotKnownPosNegPtr;
	
	// Evaluate the Chebyshev Polynomials for the given u 
	for(size_t i = 2; i < numCoeff + 1; i++){
		if(i % 2 == 0){// Once the index is even 
			// Create the ciphertext pointer for the targeted chebyshev polynomial
			Ciphertext *target_ChebPolyPtr 				 = &chebyshevPolynomials[i];
			Ciphertext *target_ChebPoly_MultiplicandPtr1 = &chebyshevPolynomials[i/2];
			Ciphertext *target_ChebPoly_MultiplicandPtr2 = &chebyshevPolynomials[i/2]; 	
			// Multiply T[n] * T[n]
			matrixVectorMultMatRowPacking(target_ChebPolyPtr, target_ChebPoly_MultiplicandPtr1, target_ChebPoly_MultiplicandPtr2, scale, contextPtr, evaluatorPtr, relin_keysPtr);
			// Obtain  2 * T[n] * T[n] by a single addition 						
			Ciphertext twoTimesResChebyPoly;   
			Ciphertext *twoTimesResChebyPolyPtr = &twoTimesResChebyPoly;
			addSubtractTwoVector(twoTimesResChebyPolyPtr, target_ChebPolyPtr, target_ChebPolyPtr, scale, contextPtr, evaluatorPtr, true);			
			// Make the last subtraction 
			Ciphertext resChebyPolynomial;  
			Ciphertext *resChebyPolynomialPtr = &resChebyPolynomial;
			addSubtractTwoVector(resChebyPolynomialPtr, twoTimesResChebyPolyPtr, vectorOnePtr, scale, 	contextPtr,	evaluatorPtr, false);					
			chebyshevPolynomials[i] = resChebyPolynomial; 
		}else{// Once the index is odd 
			// Create the ciphertext pointer for the targeted chebyshev polynomial	
			Ciphertext *target_ChebPolyPtr 				 = &chebyshevPolynomials[i];
			Ciphertext *target_ChebPoly_Multiplicand1Ptr = &chebyshevPolynomials[i/2];				
			Ciphertext *target_ChebPoly_Multiplicand2Ptr = &chebyshevPolynomials[i/2 + 1];
		   	// Multiply T[n] * T[n + 1]
			matrixVectorMultMatRowPacking(target_ChebPolyPtr, target_ChebPoly_Multiplicand1Ptr, target_ChebPoly_Multiplicand2Ptr, scale, contextPtr, evaluatorPtr, relin_keysPtr);
			// Obtain 2 * T[n] * T[n + 1] by a single addition
			Ciphertext twoTimesResChebyPoly;   
			Ciphertext *twoTimesResChebyPolyPtr = &twoTimesResChebyPoly;
			addSubtractTwoVector(twoTimesResChebyPolyPtr, target_ChebPolyPtr, target_ChebPolyPtr, scale, contextPtr, evaluatorPtr, true);
			// Make the last subtraction 
			Ciphertext resChebyPolynomial;  
			Ciphertext *resChebyPolynomialPtr = &resChebyPolynomial;
			addSubtractTwoVector(resChebyPolynomialPtr, twoTimesResChebyPolyPtr, vectorNotKnownPosNegPtr, scale, contextPtr, evaluatorPtr, false); // chebyshevPolynomials[1]
			chebyshevPolynomials[i] = resChebyPolynomial;
		}																																																	
	}
	chebyshevPolynomials[0] = *firsCoeffArrPtr;
		
	// Multiply the Chebyshev Polynomials with the coeffients of the polynomials 
	for(size_t i = 1; i < numCoeff + 1; i++){
		double coeff 		= coeffArr[i-1]; 
		Plaintext plain_coeff; 
		Plaintext * plaintextCoeffPtr 		= &plain_coeff; 
		encoderPtr-> encode(coeff, scale, plain_coeff);
		Ciphertext * respectiveChebyshevPtr = & chebyshevPolynomials[i]; 
		adjustScaleandChainParametersPlaintextAndVectors(respectiveChebyshevPtr, plaintextCoeffPtr, scale, contextPtr, evaluatorPtr);
		evaluatorPtr->multiply_plain_inplace(chebyshevPolynomials[i], plain_coeff); 
		evaluatorPtr->rescale_to_next_inplace(chebyshevPolynomials[i]);
	}

	// Add each product on top of each other
	// evaluatorPtr->add_many(chebyshevPolynomials, *sumOutputPtr); // BUGGY: DOES NOT WORK DUE TO PARAMETER MISMATCH	 	
	vector <Ciphertext> * chebVectorPtr = &chebyshevPolynomials;  
	addManyVectorsInaRow(sumOutputPtr, chebVectorPtr, numCoeff + 1, scale, contextPtr, evaluatorPtr, encoderPtr, gal_keysPtr, relin_keysPtr);	
}

/*
	The function for performing the multiple polynomials addition in the ciphertext-only setting (actually all the terms except the first term)   
*/
void addManyVectorsInaRow(Ciphertext *sumOutputPtr, vector <Ciphertext> * chebVectorPtr, int exactNumCoeff, double scale, SEALContext *contextPtr, Evaluator *evaluatorPtr, CKKSEncoder *encoderPtr, GaloisKeys *gal_keysPtr, RelinKeys *relin_keysPtr){
		
	// Find the minimum chain index	
	int min_chain_index = 999;
	parms_id_type min_parms_id;
	for(size_t i = 0; i < exactNumCoeff; i++){
		int chain_index = (*contextPtr->get_context_data(chebVectorPtr->at(i).parms_id())).chain_index();
		if(chain_index < min_chain_index){
			min_chain_index = chain_index;
			min_parms_id 	= chebVectorPtr->at(i).parms_id();
		} 
	}	

	// Adjust each chain index to the minimum chain index and scale to the already-setted scale, so that they could be added on top of another 
	for(size_t i = 0; i < exactNumCoeff; i++){
    	evaluatorPtr->mod_switch_to_inplace(chebVectorPtr->at(i), min_parms_id);	
		chebVectorPtr->at(i).scale() = scale;	
	}

    // Add each Chebyshev polynomial on top of another
	Ciphertext *ChebyshevPtr0 = &chebVectorPtr->at(0);
	Ciphertext *ChebyshevPtr1 = &chebVectorPtr->at(1);
    addSubtractTwoVector(sumOutputPtr, ChebyshevPtr0, ChebyshevPtr1, scale, contextPtr, evaluatorPtr, true);
	for(size_t i = 2; i < exactNumCoeff; i++){
		evaluatorPtr->add_inplace(*sumOutputPtr, chebVectorPtr->at(i));
	}
}

