/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which is used to handle the fundamental homomorphic operations such as addition, multiplication etc. 	
					      between Plaintext&Ciphertext--Ciphertext&Ciphertext or through only a ciphertext 
	Note				: Commenting and cleaning up has been on July 2022
*/

/* Import the other classes' header files*/
#include "examples.h" 			  // The class for the console menu and guiding the user to the preferred application     	
#include "rawplain.h"  			  // The class which reads and stores the plain matrices                    
#include "generateplaintextMRP.h" // The functions used to prepare plain matrices in MRP format
#include "prepareVecMatMRP.h"     // The class which prepares the read matrices in MRP format  		
#include "initializationPLCP.h"   // The class containing the functions which convert the prepared matrices to Plaintext and Ciphertext objects.            
#include "plcpOperations.h"       // The class containing the functions which does Ciphertext-Ciphertext and Plaintext-Ciphertext arithmetic and algebraic operations  	    
#include "encryptedAppx.h"        // The class containing the functions which does the Chebyshev approximation with different assumptions     
#include "secretShare.h" 		  // The class containing the functions which does the secret sharing  
#include "printCont.h" 			  // The class containing the printing functions for control purposes  

/* Import the important selected C libraries*/
#include <iostream>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

// ##################################################################################
// ##################################################################################
// ############################PLAINTEXT-CIPHERTEXT##################################
// ##################################################################################
// ##################################################################################

/*
	The function for finding the lowest number of slots equaling to a binary power for the matrix row packing format (one after another)   
*/
int findBinaryPower(int indexVectorLength){
	int power = 0;
	if(indexVectorLength != 1){
		int startLength = 1;		
		while (startLength != indexVectorLength){
			startLength = startLength * 2; 
			power++; 	
		}
		return power;			
	}else{
		return 0;	
	}
}

/*
	The function for equalizing the chain indexes and scales of a plaintext and a ciphertext to make a homomorphic operation (e.g., Addition, multiplication, etc.)      
*/
void adjustScaleandChainParametersPlaintextAndVectors( 
	Ciphertext *CiphertextPtr,	
	Plaintext *plaintextPtr,
	double scale, 
	SEALContext *contextPtr,	
	Evaluator *evaluatorPtr){

	// Obtain Chain indexes		
	int first_chain_index 	= (*contextPtr->get_context_data(CiphertextPtr->parms_id())).chain_index();
	int second_chain_index 	= (*contextPtr->get_context_data(plaintextPtr->parms_id())).chain_index();	
	
	// Adjust Chain parameters if needed 
	if(first_chain_index < second_chain_index){
		parms_id_type min_parms_id = CiphertextPtr->parms_id();
    	evaluatorPtr->mod_switch_to_inplace(*plaintextPtr, min_parms_id);
	}else if(second_chain_index < first_chain_index){
		parms_id_type min_parms_id = plaintextPtr->parms_id();
    	evaluatorPtr->mod_switch_to_inplace(*CiphertextPtr, min_parms_id);
	} 
	
	// Adjust the scales if needed 	
	if(CiphertextPtr->scale() != plaintextPtr->scale()){
		CiphertextPtr->scale() = scale;
		plaintextPtr->scale() = scale;	
	}

}

/*
	The function for performing a not-in-place homomorphic multiplication with Plaintext&Ciphertext       
*/
void plnCprAxMult_mrp(
	Ciphertext *res, 
	Ciphertext *cpr,
	Plaintext *pln,
	double scale, 
	Decryptor *decryptorPtr,
	CKKSEncoder *encoderPtr,
	SEALContext *contextPtr,
	Evaluator *evaluatorPtr,
	RelinKeys *relin_keysPtr){
		
	// Adjust the scales if needed
	adjustScaleandChainParametersPlaintextAndVectors(cpr, pln, scale, contextPtr, evaluatorPtr);

	// Do the encrypted multiplication	
	evaluatorPtr->multiply_plain(*cpr, *pln, *res);
	evaluatorPtr->relinearize_inplace(*res, *relin_keysPtr);
	evaluatorPtr->rescale_to_next_inplace(*res);

}

/*
	The function for performing an in-place homomorphic multiplication with Plaintext&Ciphertext       
*/
void plnCprAxMult_InPlac_mrp(
	Ciphertext *res, 
	Plaintext *pln,
	double scale, 
	Decryptor *decryptorPtr,
	CKKSEncoder *encoderPtr,
	SEALContext *contextPtr,
	Evaluator *evaluatorPtr,
	RelinKeys *relin_keysPtr){
		
	// Adjust the scales if needed
	adjustScaleandChainParametersPlaintextAndVectors(res, pln, scale, contextPtr, evaluatorPtr);

	// Do the encrypted multiplication	
	evaluatorPtr->multiply_plain_inplace(*res, *pln);
	evaluatorPtr->relinearize_inplace(*res, *relin_keysPtr);
	evaluatorPtr->rescale_to_next_inplace(*res);

}

/*
	The function for performing an not-in-place homomorphic addition or subtraction with Plaintext&Ciphertext       
*/
void addSubtractPLCPVector(
	Ciphertext *matrixVecAddRes, 
	Ciphertext *v_1_Ciph,
	Plaintext  *v_2_Pln,
	double scale, 	
	SEALContext *contextPtr,	
	Evaluator *evaluatorPtr, 
	bool isAddition){

	// Adjust the scales if needed
	adjustScaleandChainParametersPlaintextAndVectors(v_1_Ciph, v_2_Pln, scale, contextPtr, evaluatorPtr);
	
	// Perform the selected operation
	if(isAddition == true) 
		evaluatorPtr->add_plain(*v_1_Ciph, *v_2_Pln, *matrixVecAddRes); // Make addition operation
	else
		evaluatorPtr->sub_plain(*v_1_Ciph, *v_2_Pln, *matrixVecAddRes); // Make subtraction operation

}


// ##################################################################################
// ##################################################################################
// ############################CIPHERTEXT-CIPHERTEXT#################################
// ##################################################################################
// ##################################################################################


/*
	The function for performing the rotation through a Ciphertext to sum the first n entries into the 1st index      
*/
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
	RelinKeys *relin_keysPtr){

	// Find the binary power
	int powerBinary = findBinaryPower(indexVectorLength);

	// Create the vector of Ciphertexts for the rotation operation 
	vector<Ciphertext> rotations_output(powerBinary + 1);
	rotations_output[0] = *v_to_Rotate;
	int rotOutputInd = 0; 
	for(size_t steps = indexVectorLength/2; steps >= 1; steps = steps / 2){
		
		// Rotate the prior ciphertext   	
		Ciphertext rotated_result;
		Ciphertext *rotated_resultPtr = &rotated_result;
		evaluatorPtr->rotate_vector(rotations_output[rotOutputInd], steps, *gal_keysPtr, rotated_result);
		
		// Prepare the previous rotated vector 
		Ciphertext previousRotation; 
		Ciphertext *previousRotationPtr = &previousRotation;
		previousRotation = rotations_output[rotOutputInd];
	
		// Add the new rotated vector and the previous rotated vector: 
		// Addition:  rotations_output[rotOutputInd] + rotated_result
		Ciphertext rotated_added_result;
		Ciphertext *rotated_added_resultPtr = &rotated_added_result;					
		addSubtractTwoVector(rotated_added_resultPtr, rotated_resultPtr, previousRotationPtr, scale, contextPtr, evaluatorPtr, true);
		
		// Assign the obtained sum	
		rotations_output[++rotOutputInd] = rotated_added_result; 
	}

	// Assign the concluding summation to the output variable 
	*sum_Output = rotations_output[powerBinary]; 
}


/*
	The function for equalizing the chain indexes and scales of two ciphertexts to make a homomorphic operation (e.g., Addition, multiplication, etc.)      
*/
void adjustScaleandChainParametersTwoVectors( 
	Ciphertext *Cipher1,
	Ciphertext *Cipher2,
	Ciphertext *CipherRes,	
	double scale, 
	SEALContext *contextPtr,	
	Evaluator *evaluatorPtr){

	// Obtain Chain indexes		
	int first_chain_index 	= (*contextPtr->get_context_data(Cipher1->parms_id())).chain_index();
	int second_chain_index 	= (*contextPtr->get_context_data(Cipher2->parms_id())).chain_index();	
	
	// Adjust Chain parameters if needed 
	if(first_chain_index < second_chain_index){
		parms_id_type min_parms_id = Cipher1->parms_id();
    	evaluatorPtr->mod_switch_to_inplace(*Cipher2, min_parms_id);
	}else if(second_chain_index < first_chain_index){
		parms_id_type min_parms_id = Cipher2->parms_id();
    	evaluatorPtr->mod_switch_to_inplace(*Cipher1, min_parms_id);
	} 
	
	// Adjust the scales if still  needed 	
	if(Cipher1->scale() != Cipher2->scale()){
		Cipher1->scale() = scale;
		Cipher2->scale() = scale;	
	}

}

/*
	The function for performing a not-in-place homomorphic multiplication with Ciphertext&Ciphertext       
*/
void matrixVectorMultMatRowPacking(
	Ciphertext *matrixVecMultRes, 
	Ciphertext *matrixCipher,
	Ciphertext *vectorCipher,
	double scale, 
	SEALContext *contextPtr,
	Evaluator *evaluatorPtr,
	RelinKeys *relin_keysPtr){
		
	// Adjust the scales if needed
	adjustScaleandChainParametersTwoVectors(matrixCipher, vectorCipher, matrixVecMultRes, scale, contextPtr, evaluatorPtr);

	// Do the encrypted multiplication	
	evaluatorPtr->multiply(*matrixCipher, *vectorCipher, *matrixVecMultRes);
	evaluatorPtr->relinearize_inplace(*matrixVecMultRes, *relin_keysPtr);
	evaluatorPtr->rescale_to_next_inplace(*matrixVecMultRes);

}

/*
	The function for performing a in-place homomorphic squaring operation with a ciphertext       
*/
void vectorSquaringRowPacking(
	Ciphertext *vectorCipherSquareRes,	
	Ciphertext *vectorCipher,
	double scale, 
	SEALContext *contextPtr,
	Evaluator *evaluatorPtr,
	RelinKeys *relin_keysPtr){
		
	// Do the encrypted squaring operation	
	evaluatorPtr->square(*vectorCipher, *vectorCipherSquareRes);
	evaluatorPtr->relinearize_inplace(*vectorCipherSquareRes, *relin_keysPtr);
	evaluatorPtr->rescale_to_next_inplace(*vectorCipherSquareRes);

}

/*
	The function for performing a homomorphic addition of three ciphertexts       
*/
void addThreeVector(Ciphertext *matrixVecAddRes, 
	Ciphertext *v_1_Ciph,
	Ciphertext *v_2_Ciph,
	Ciphertext *v_3_Ciph,
	double scale, 	
	SEALContext *contextPtr,	
	Evaluator *evaluatorPtr,
	RelinKeys *relin_keysPtr){

	// Adjust the chain indexes 
	int first_chain_index 	= (*contextPtr->get_context_data(v_1_Ciph->parms_id())).chain_index();
	int second_chain_index 	= (*contextPtr->get_context_data(v_2_Ciph->parms_id())).chain_index();	
	int third_chain_index   = (*contextPtr->get_context_data(v_3_Ciph->parms_id())).chain_index();	
	
	if(first_chain_index > second_chain_index){
		if(third_chain_index > second_chain_index){
			parms_id_type second_parms_id = v_2_Ciph->parms_id();
    		evaluatorPtr->mod_switch_to_inplace(*v_1_Ciph, second_parms_id);
			evaluatorPtr->mod_switch_to_inplace(*v_3_Ciph, second_parms_id);
		}else{
			parms_id_type third_parms_id = v_3_Ciph->parms_id();
    		evaluatorPtr->mod_switch_to_inplace(*v_1_Ciph, third_parms_id);
			evaluatorPtr->mod_switch_to_inplace(*v_2_Ciph, third_parms_id);
		}		
	}else{
		if(third_chain_index > first_chain_index){
			parms_id_type first_parms_id = v_1_Ciph->parms_id();
    		evaluatorPtr->mod_switch_to_inplace(*v_2_Ciph, first_parms_id);
			evaluatorPtr->mod_switch_to_inplace(*v_3_Ciph, first_parms_id);
		}else{
			parms_id_type third_parms_id = v_3_Ciph->parms_id();
    		evaluatorPtr->mod_switch_to_inplace(*v_1_Ciph, third_parms_id);
			evaluatorPtr->mod_switch_to_inplace(*v_2_Ciph, third_parms_id);
		}
	}		
	
	// Adjust the scales 
	v_1_Ciph->scale() = scale;
	v_2_Ciph->scale() = scale;	
	v_3_Ciph->scale() = scale; 	

	// Make the addition operation
	Ciphertext matrixVecAddRes_Temp;
	evaluatorPtr->add(*v_1_Ciph, *v_2_Ciph, matrixVecAddRes_Temp);
	evaluatorPtr->add(matrixVecAddRes_Temp, *v_3_Ciph, *matrixVecAddRes);
}

/*
	The function for performing an not-in-place homomorphic addition or subtraction with Ciphertext&Ciphertext       
*/
void addSubtractTwoVector(Ciphertext *matrixVecAddRes, 
	Ciphertext *v_1_Ciph,
	Ciphertext *v_2_Ciph,
	double scale, 	
	SEALContext *contextPtr,	
	Evaluator *evaluatorPtr, 
	bool isAddition){

	// Adjust the scales if needed
	adjustScaleandChainParametersTwoVectors(v_1_Ciph, v_2_Ciph, matrixVecAddRes, scale, contextPtr, evaluatorPtr);

	// Perform the selected operation
	if(isAddition == true) 
		evaluatorPtr->add(*v_1_Ciph, *v_2_Ciph, *matrixVecAddRes); // Make addition operation
	else
		evaluatorPtr->sub(*v_1_Ciph, *v_2_Ciph, *matrixVecAddRes); // Make subtraction operation
}	
