/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which is used for implementing the main operational functions   
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
#include "applyPLCPSimulation.h"  // The class containing the crypto application functions for each targeted equation

/* Import the important selected C libraries*/
#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* Define the constant */
#define numRandBits 20

/* Call main namespaces */
using namespace std;
using namespace seal;

/*
	The function for performing the cyberphysical system's sensor measurement functionality   
*/
void sense_Encrypt_y(struct simulationMatrixData *smd, struct simulationMatrixMRP *smrp, int numiter, Ciphertext *yy_CP, Ciphertext *yyAS_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, CKKSEncoder *encoderPtr){

	// Create the respective sensor data 	
	double ** ySensorData = (double **) calloc(smd->n, sizeof(double *));
	for(int i = 0; i < smrp->n; i++)
		ySensorData[i] = (double *) calloc(1, sizeof(double));	
	
	// Obtain the sensor measurement data 
	for(int i = 0; i < smrp->n; i++)
		ySensorData[i][0] = smd->xx[i][0] + smd->yNoise[numiter][i];

	// Anomaly attack detection generation through the sensor measurement 
	if(numiter - 1 > 79 && numiter - 1 < 82){ //{
		cout << numiter << "th iteration one time attack" << endl;    
		ySensorData[4][0] = ySensorData[4][0] + 4;	
	}
		
	// Define the proper y vectors to be encoded and encrypted 
	vector<double> * ySensorDataRep = new std::vector<double>(smrp->N * smrp->n); // 2nd equation (Estimation)  
	vector<double> * ySensorDataZer = new std::vector<double>(smrp->N * smrp->n); // 2nd equation (Estimation)
	genRepXVecMRP_RPL_v0(ySensorDataRep, ySensorData, smrp->n, smrp->n, smrp->N, "Sensor-Measurement-Y-Rep", false);
	genYVecAddOperMRP_RPL_v0(ySensorDataZer, ySensorData, smrp->n, smrp->n, smrp->N, "Sensor-Measurement-Y-Add", false);
	
	// Encode and encrypt the y vectors
	encryptXVectorMatRowPacking(ySensorDataRep, scale, yy_CP, encryptorPtr, encoderPtr);
	encryptXVectorMatRowPacking(ySensorDataZer, scale, yyAS_CP, encryptorPtr, encoderPtr);

}

/*
	The function for performing the cyberphysical system's estimation functionality   
*/
void applyEquation_2_PLCP(Ciphertext * secEqRes_CP, struct simulationMatrixMRP *smrp, Plaintext *GAMMA_PL, Plaintext *LL_PL, Ciphertext * xGxG_CP,  Ciphertext * yy_CP, Ciphertext * xexe_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr){

	/* 	
		Aim		: Generated for the 2nd Equation (i.e., The Estimation Computation)  
		Equation 2 -> x̂e[k] = (A − LA − BK + LBK)x̂ e [k − 1] + Ly[k] + K(B − LB)xr + (B − LB)ur
		Equation 2 -> x̂e[k] = Γx̂e[k−1] + Ly[k] + xΓ  
		Equation 2 -> x̂e[k] = Γx̂e[k − 1] (1st part) + Ly[k] (2nd part) + xΓ (3rd part) 
	*/
	
	// 1st part of the estimation equation 
	Ciphertext _2nd_eq_1p;
	Ciphertext *_2nd_eq_1p_Pt = &_2nd_eq_1p; 	 	
	plnCprAxMult_mrp(_2nd_eq_1p_Pt, xexe_CP, GAMMA_PL, scale, decryptorPtr, encoderPtr, contextPtr, evaluatorPtr, relin_keysPtr);
	
	// 2nd part of the estimation equation
	Ciphertext _2nd_eq_2p;
	Ciphertext *_2nd_eq_2p_Pt = &_2nd_eq_2p; 
	plnCprAxMult_mrp(_2nd_eq_2p_Pt, yy_CP, LL_PL, scale, decryptorPtr, encoderPtr, contextPtr, evaluatorPtr, relin_keysPtr);

	// 3rd part		
	Ciphertext _2nd_eq_3p;
	Ciphertext *_2nd_eq_3p_Pt = &_2nd_eq_3p; 
	addThreeVector(_2nd_eq_3p_Pt, _2nd_eq_1p_Pt, _2nd_eq_2p_Pt, xGxG_CP, scale, contextPtr, evaluatorPtr, relin_keysPtr); // Do encrypted final addition
	
	// Do rotation and addition
	rotateVector(secEqRes_CP, _2nd_eq_3p_Pt, scale, smrp->N, contextPtr, evaluatorPtr, decryptorPtr, encoderPtr, gal_keysPtr, relin_keysPtr);
}

/*
	The function for performing the cyberphysical system's control action functionality   
*/
void applyEquation_3_PLCP(Ciphertext * thirdEqRes_CP, struct simulationMatrixMRP *smrp, Plaintext *KGKG_PL, Plaintext *KLKL_PL, Ciphertext * KxugKxug_CP,  Ciphertext * yy_CP, Ciphertext * xexe_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr){
	/* 	
		Aim		: Generated for the 3rd Equation (i.e., The Control Action Computation)  
		Equation 3 -> u[k]  = -KΓx̂ e [k − 1] -KLy[k] - KxΓ + uΓ   
		Equation 3 -> u[k]  = KG*x̂e[k−1] + KL * y[k] + KXUΓ  
		Equation 3 -> u[k]  = KG*x̂e[k−1] (1st part) + KL * y[k] (2nd part) + KXUΓ (3rd part) 
	*/
	
	// 1st part
	Ciphertext _3rd_eq_1p;
	Ciphertext *_3rd_eq_1p_Pt = &_3rd_eq_1p;
	plnCprAxMult_mrp(_3rd_eq_1p_Pt, xexe_CP, KGKG_PL, scale, decryptorPtr, encoderPtr, contextPtr, evaluatorPtr, relin_keysPtr);	

	// 2nd part 
	Ciphertext _3rd_eq_2p;
	Ciphertext *_3rd_eq_2p_Pt = &_3rd_eq_2p;
	plnCprAxMult_mrp(_3rd_eq_2p_Pt, yy_CP, KLKL_PL, scale, decryptorPtr, encoderPtr, contextPtr, evaluatorPtr, relin_keysPtr);

	// 3rd part
	Ciphertext _3rd_eq_3p;
	Ciphertext *_3rd_eq_3p_Pt = &_3rd_eq_3p;
	addThreeVector(_3rd_eq_3p_Pt, _3rd_eq_1p_Pt, _3rd_eq_2p_Pt, KxugKxug_CP, scale, contextPtr, evaluatorPtr, relin_keysPtr); // Do encrypted final addition
	
	// Do rotation and addition
	rotateVector(thirdEqRes_CP, _3rd_eq_3p_Pt, scale, smrp->N, contextPtr, evaluatorPtr, decryptorPtr, encoderPtr, gal_keysPtr, relin_keysPtr);
	
}

/*
	The function for performing the cyberphysical system's control action functionality at the very first iteration   
*/
void applyEquation_3_fiter_PLCP(Ciphertext * thirdEqRes_CP, struct simulationMatrixMRP *smrp, Plaintext *KxKx_PL, Ciphertext * uGuG_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr){

	/* 	
		Aim		: Generated for the 3rd Equation (i.e., The Control Action Computation)  
		Equation 3 -> u[k]  = -K*xe(:, k) + uG;      
		Equation 3 -> x̂e[k] = KG*x̂e[k−1] + KL * y[k] + KXUΓ  
		Equation 3 -> x̂e[k] = KG*x̂e[k−1] (1st part) + KL * y[k] (2nd part) + KXUΓ (3rd part) 
	*/

	addSubtractPLCPVector(thirdEqRes_CP, uGuG_CP, KxKx_PL, scale, contextPtr, evaluatorPtr, true);	
}

/*
	The function for performing the cyberphysical system's prediction functionality   
*/
void applyEquation_4_5_PLCP(Ciphertext * fourthfifthEqRes_CP, struct simulationMatrixMRP *smrp, Plaintext *ACL_PL, Plaintext *BB_PL, Ciphertext * uGuG_CP, Ciphertext * xexe_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr){	
	/* 	
		Aim		: Generated for the 4-5th Equation (i.e., The Prediction Computation)  
		Equation 4-5 -> x̂p[k] = Acl*x̂e[k − 1] + B*uΓ    
		Equation 4-5 -> x̂p[k] = Acl*x̂e[k − 1] (1st part) + B*uΓ (2nd part) 
	*/

	// 1st part 
	Ciphertext _4_5th_eq_1p;
	Ciphertext *_4_5th_eq_1p_Pt = &_4_5th_eq_1p;
	plnCprAxMult_mrp(_4_5th_eq_1p_Pt, xexe_CP, ACL_PL, scale, decryptorPtr, encoderPtr, contextPtr, evaluatorPtr, relin_keysPtr);	

	// 2nd part
	Ciphertext _4_5th_eq_2p;
	Ciphertext *_4_5th_eq_2p_Pt = &_4_5th_eq_2p;
	plnCprAxMult_mrp(_4_5th_eq_2p_Pt, uGuG_CP, BB_PL, scale, decryptorPtr, encoderPtr, contextPtr, evaluatorPtr, relin_keysPtr);
	
	// 3rd part 
	Ciphertext _4_5th_eq_3p;
	Ciphertext *_4_5th_eq_3p_Pt = &_4_5th_eq_3p;
	addSubtractTwoVector(_4_5th_eq_3p_Pt, _4_5th_eq_1p_Pt, _4_5th_eq_2p_Pt, scale, 	contextPtr,	evaluatorPtr, true);

	// Do rotation and addition
	rotateVector(fourthfifthEqRes_CP, _4_5th_eq_3p_Pt, scale, smrp->N, contextPtr, evaluatorPtr, decryptorPtr, encoderPtr, gal_keysPtr, relin_keysPtr);	
}

/*
	The function for performing the cyberphysical system's residues computation functionality   
*/
void applyEquation_6_PLCP(Ciphertext * sixthEqRes_CP, struct simulationMatrixMRP * smrp, Ciphertext * fourthfifthEqRes_CP, Ciphertext * yyAS_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr){
	/* 	
		Aim		: Generated for the 6th Equation (i.e., The Residues Computation)  
		Equation 6 ->  (y[k] - x̂p[k])^2      
	*/

	// 1st part
	Ciphertext _6th_eq_1p;
	Ciphertext *_6th_eq_1p_Pt = &_6th_eq_1p;
	addSubtractTwoVector(_6th_eq_1p_Pt, yyAS_CP, fourthfifthEqRes_CP, scale, contextPtr, evaluatorPtr, false);

	// 2nd part
	vectorSquaringRowPacking(sixthEqRes_CP, _6th_eq_1p_Pt, scale, contextPtr, evaluatorPtr, relin_keysPtr);	
}

/*
	The function for performing the cyberphysical system's residues computation functionality at the very first iteration   
*/
void applyEquation_6_fiter_PLCP(Ciphertext * sixthEqRes_CP, struct simulationMatrixMRP * smrp, Plaintext * xpxp_PL, Ciphertext * yyAS_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr){

	/* 	
		Aim		: Generated for the 6th Equation (i.e., The Residues Computation)  
		Equation 6 ->  (y[k] - x̂p[k])^2      
	*/

	// 1st part
	Ciphertext _6th_eq_1p;
	Ciphertext *_6th_eq_1p_Pt = &_6th_eq_1p;
	addSubtractPLCPVector(_6th_eq_1p_Pt, yyAS_CP, xpxp_PL, scale, contextPtr, evaluatorPtr, false);

	// 2nd part
	vectorSquaringRowPacking(sixthEqRes_CP, _6th_eq_1p_Pt, scale, contextPtr, evaluatorPtr, relin_keysPtr);	
}

/*
	The function for performing the cyberphysical system's 1st part of CUSUM Computation and Alarm Computation    
*/
void applyEquation_CUSUM_PLCP(Ciphertext * eigthEqRes_CP, Ciphertext * ninthEqRes_CP, struct simulationMatrixMRP * smrp, struct simulationMatrixData * smd, int numiter, 
							Plaintext * ss_PL, Plaintext * vv_PL, Plaintext * TAU_PL, Plaintext * alpha_Eq8_PL, Plaintext * beta_Eq8_PL, 
							Plaintext * alpha_Eq9_PL, Plaintext * beta_Eq9_PL, 
							Plaintext *vectorOnePtr_PL,  
							Plaintext *firstPowerSeriesTermEq8Ptr_PL, vector <Plaintext> * powSerCoeffArr_Eq8_PL,  
							Plaintext *firstPowerSeriesTermEq9Ptr_PL, vector <Plaintext> * powSerCoeffArr_Eq9_PL, 
							Ciphertext * sixthEqRes_CP, Ciphertext * ss_CP, 
							double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, 
							Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr, bool isFirstIter){

	/*
			Aim		: Generated for the 8th and 9th Equation  
			Equation 8 ->  s̄[k + 1]    = max(r[k](i) + s[k](i) - v(i), 0) (i.e., The 1st part of the CUSUM Computation)   
			Equation 9 ->  alarm[k](i) = Ind(r[k](i) + s[k](i) - v(i) - tau(i)) if it is greater than 0 -> 1, else -> 0 (i.e., The Alarm Computation)   
	*/

	// ======================================================================================== 
	// =============== The 1st part of the CUSUM Computation-EQUATION-8 ======================= 
	// ========================================================================================   
	Ciphertext _8th_eq_1p;
	Ciphertext *_8th_eq_1p_Pt = &_8th_eq_1p;
	if(isFirstIter)
		addSubtractPLCPVector(_8th_eq_1p_Pt, sixthEqRes_CP, ss_PL, scale, contextPtr, evaluatorPtr, true);	
	else 
		addSubtractTwoVector(_8th_eq_1p_Pt, sixthEqRes_CP, ss_CP, scale, contextPtr, evaluatorPtr, true);		
	
	// 2nd part 
	Ciphertext _8th_eq_2p;
	Ciphertext * _8th_eq_2p_Pt = &_8th_eq_2p;
	addSubtractPLCPVector(_8th_eq_2p_Pt, _8th_eq_1p_Pt, vv_PL, scale, contextPtr, evaluatorPtr, false);	
		
	// Record the RELU Input	
	// extractExpRes(_8th_eq_2p_Pt, smrp, smd, numiter, smd->n, scale, contextPtr, decryptorPtr, encoderPtr, "sBar");

	// APPLY CHEBYSHEV APPROXIMATION FOR THE MAXIMUM COMPUTATION   
	// 3rd part-1: Multiply with Alpha
	// printf("\n\n===Plaintext-Encrypted or Encrypted max(r[k](i) + s[k](i) - v(i), 0)-1st-Part-Range Transformation-8theq-3rdpart-Beginning===\n\n");
	Ciphertext _8th_eq_3p_1;	
	Ciphertext * _8th_eq_3p_1_Pt = &_8th_eq_3p_1;
	plnCprAxMult_mrp(_8th_eq_3p_1_Pt, _8th_eq_2p_Pt, alpha_Eq8_PL, scale, decryptorPtr, encoderPtr, contextPtr, evaluatorPtr, relin_keysPtr);
	// 3rd part-2: Subtract beta + 1
	// printf("\n\n===Plaintext-Encrypted or Encrypted max(r[k](i) + s[k](i) - v(i), 0)-2nd-Part-Range Transformation-8theq-3rdpart-Beginning===\n\n");
	Ciphertext _8th_eq_3p_2;	
	Ciphertext * _8th_eq_3p_2_Pt = &_8th_eq_3p_2;
	addSubtractPLCPVector(_8th_eq_3p_2_Pt, _8th_eq_3p_1_Pt, beta_Eq8_PL, scale, contextPtr, evaluatorPtr, false);
	// 3rd part-3: Apply Chebyshev Polynomial Appx.  
	makeChebyshevPolynAppxPLCP(smrp, eigthEqRes_CP, _8th_eq_3p_2_Pt,  vectorOnePtr_PL,  firstPowerSeriesTermEq8Ptr_PL,  powSerCoeffArr_Eq8_PL, smrp->chebDegEq8, scale, contextPtr, evaluatorPtr, encoderPtr, gal_keysPtr, relin_keysPtr, decryptorPtr);


	// ======================================================================================== 
	// ========================= The Alarm Computation-EQUATION-9 ============================= 
	// ======================================================================================== 

	// 1st part
	Ciphertext _9th_eq_1p;
	Ciphertext *_9th_eq_1p_Pt = &_9th_eq_1p;
	addSubtractPLCPVector(_9th_eq_1p_Pt, ss_CP, TAU_PL, scale, contextPtr, evaluatorPtr, false);	

	// APPLY CHEBYSHEV APPROXIMATION FOR THE MAXIMUM COMPUTATION   
	// 2nd part-1: Multiply with Alpha
	Ciphertext _9th_eq_2p_1;	
	Ciphertext * _9th_eq_2p_1_Pt = &_9th_eq_2p_1;
	plnCprAxMult_mrp(_9th_eq_2p_1_Pt, _9th_eq_1p_Pt, alpha_Eq9_PL, scale, decryptorPtr, encoderPtr, contextPtr, evaluatorPtr, relin_keysPtr);
	// 2nd part-2: Subtract beta + 1 
	Ciphertext _9th_eq_2p_2;	
	Ciphertext * _9th_eq_2p_2_Pt = &_9th_eq_2p_2;
	addSubtractPLCPVector(_9th_eq_2p_2_Pt, _9th_eq_2p_1_Pt, beta_Eq9_PL, scale, contextPtr, evaluatorPtr, false);
	// Prepare vector one for the 9th equation Chebyshev   
	Plaintext PL_vecOne_n;
	Plaintext *vecOne_n_PL;
	vecOne_n_PL = & PL_vecOne_n;
	makePlaintextMatRowPacking(smrp->One_MRP, scale, vecOne_n_PL, encoderPtr);
	// 2nd part-3: Apply Chebyshev Polynomial Appx.	 
	makeChebyshevPolynAppxPLCP(smrp, ninthEqRes_CP, _9th_eq_2p_2_Pt,  vecOne_n_PL,  firstPowerSeriesTermEq9Ptr_PL,  powSerCoeffArr_Eq9_PL, smrp->chebDegEq9, scale, contextPtr, evaluatorPtr, encoderPtr, gal_keysPtr, relin_keysPtr, decryptorPtr);	
}

/*
	The function for performing the cyberphysical system's 2nd and last part of CUSUM Computation    
*/
void applyEquation_10_PLCP(Ciphertext * tenthEqRes_CP, struct simulationMatrixMRP * smrp, Ciphertext *eigthEqRes_CP, Ciphertext * ninthEqRes_CP, Ciphertext *vectorOnePtr_CP, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, RelinKeys *relin_keysPtr, GaloisKeys *gal_keysPtr){
		
	/*
		Aim		: Generated for the 10th Equation (i.e., The 2nd and last part of the CUSUM Computation)  
		Equation 10 ->  s[k + 1] = s̄[k + 1] ⊙ (1 − alarm [k])
		Equation 10 ->  s[k + 1] = s̄[k + 1] (2nd part) ⊙ (1 − alarm [k]) (1st part) 
		Note	: Special secret sharing which zeroes out the indices where alarm = 1 is done 
	*/ 

	// 1st part: Decryption and Rounding of alarm[k]   
	Plaintext plain_alarm;	
	vector<double> alarmRes;
	// Decrypt the ciphertext to the plaintext 
	decryptorPtr->decrypt(*ninthEqRes_CP, plain_alarm);
	// Decode the plaintext to the not-encrpted result
	encoderPtr->decode(plain_alarm, alarmRes);
	for (size_t i = 0; i < smrp->n; i++){	
		for(size_t j = 0; j < smrp->N; j++){
			double roundedAlarmValue = round(alarmRes[i * smrp->N + j]); 	
			if(roundedAlarmValue > 0.4)
				alarmRes[i * smrp->N + j] = 1;	
			else 
				alarmRes[i * smrp->N + j] = 0;			
		}
	}

	// 2nd part: Determination of CUSUM Multiplicator 
	// Initialize the alarm multiplicator with ones  		 
	int * CUSUMParamSumMultiplicator = (int *) calloc(smrp->n, sizeof(int));  	
	for (size_t i = 0; i < smrp->n; i++)
		CUSUMParamSumMultiplicator[i] = 1;
	// Update the alarm multiplicator based on the alarm value - Assign zero when alarm = 1
	for (size_t i = 0; i < smrp->n * smrp->N; i++){
		int quot = i / smrp->N;  
		if(i % smrp->N == 0)
			if(alarmRes[i] == 1)
				CUSUMParamSumMultiplicator[quot] = 0;
		// Sanity Check
		// if(i == smrp->N * 5 || i == smrp->N * 7) // Added for sanity check, to be commented after testing 
		//	CUSUMParamSumMultiplicator[quot] = 0; // Added for sanity check, to be commented after testing 				
	} 		

	// 3rd part: Special Secret Sharing  
	// Create Noise vector for the initial addition
	vector<double> rand_double_add(smrp->n * smrp->N);		
	for (size_t  i = 0; i < smrp->n * smrp->N; i++){	
		int randomNumSamplingInterval = pow(2, numRandBits); // Check for different values == // 4194304; // 134217728; // 4194304; // 500; // 4294967296; // == 
		double ran = rand() % randomNumSamplingInterval;		
		rand_double_add[i] = ran;					
	}
	// Create and specially arrange noise vector for the eventual subtraction based on the alarm results (different than the addition vector)
	vector<double> rand_double_sub(smrp->n * smrp->N);
	for (size_t  i = 0; i < smrp->n * smrp->N; i++){
		// int quot = i / smrp->N; // Added for sanity check (to be removed after testing) 
		if(i % smrp->N == 0){
			if(alarmRes[i] == 1) 
				rand_double_sub[i] = 0.0000000000000000;
			else 			
				rand_double_sub[i] = rand_double_add[i];
		}else {
		 	rand_double_sub[i] = rand_double_add[i];
		}
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
	addSubtractTwoVector(x_Random_Noise_AddedPtr, eigthEqRes_CP, randAddEncPtr, scale, contextPtr, evaluatorPtr, true);
	
	// Decrypt the noise added vector 							 			
	Plaintext plain_SecretShare_Decrypt_Result;	
	vector<double> decryptedVec;
	decryptorPtr->decrypt(x_Random_Noise_Added, plain_SecretShare_Decrypt_Result);
	encoderPtr->decode(plain_SecretShare_Decrypt_Result, decryptedVec);	
	
	// Zeroing out the decrypted CUSUM RELU computation
	for (size_t i = 0; i < smrp->n; i++){ 
		if(CUSUMParamSumMultiplicator[i] == 0)
			decryptedVec[i * smrp->N] = 0;
	}

	// Encrypt the index summed decryption vector of x^e[k-1]
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
	*tenthEqRes_CP = x_Random_Noise_Subtracted;	
}

/*
	The function for performing the cyberphysical system's sensor measurement simulation in a naturally noisy environment    
*/
void applyXVecNoiseAddition(Ciphertext * controlAction_CP, int numIter, struct  simulationMatrixMRP * smrp, struct  simulationMatrixData * smd, double scale, SEALContext *contextPtr, Encryptor *encryptorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr){

	/*
		Aim	  : Generated for the updated sensor measurements (i.e., Apply noise addition operation at the end of an online iteration)   
		Input : Cleartext/Plaintext A, Cleartext/Plaintext x, 
				Cleartext/Plaintext B, Ciphertext u 
				Cleartext/Plaintext recorded noise from the cleartext simulation (i.e., F*mvnrnd(zeros(4,1),W)')		
		Output: Updated x vector (with noise)
		Last Equation ->  x(:, k+1) = A*x(:, k) + B*u(:, k) + F*mvnrnd(zeros(4,1),W)'
		Last Equation ->  x(:, k+1) = A*x(:, k) (1st part) + B*u(:, k) (2nd part) + F*mvnrnd(zeros(4,1),W)'(3rd part)  
	*/ 

	// Initialize the pointers (i.e., Axk, u, Bu, AxkBuNoise) 	
	// Initialize and define the outer double pointers
	double ** Axk 		 = (double **) calloc(smrp->n, sizeof(double *)); 
	double ** u 		 = (double **) calloc(smrp->m, sizeof(double *));
	double ** Bu 		 = (double **) calloc(smrp->n, sizeof(double *));
	double ** xSpecNoise = (double **) calloc(smrp->n, sizeof(double *));
	double ** AxkBuNoise = (double **) calloc(smrp->n, sizeof(double *));
	// Initialize and define the inner single pointers
	for(int i = 0; i < smrp->n; i++){
		Axk[i] 		  = (double *) calloc(1, sizeof(double));
		if(i < smd->m)			 	
			u[i] = (double *) calloc(1, sizeof(double));
		Bu[i] 		  = (double *) calloc(1, sizeof(double));
		xSpecNoise[i] = (double *) calloc(1, sizeof(double));	 		
		AxkBuNoise[i] = (double *) calloc(1, sizeof(double)); 		
	}

	// Decrypt the control action and obtain cleartext u vector 
	Plaintext pl_cont;	
	vector<double> controlRes;	
	decryptorPtr->decrypt(*controlAction_CP, pl_cont);  
	encoderPtr->decode(pl_cont, controlRes);
	for (size_t i = 0; i < smrp->m; i++){	
		for(size_t j = 0; j < smrp->N; j++){
			if(j == 0)
				u[i][j] = controlRes[i * smrp->N + j];
		}	
	}

	// 1st part: Multiply A and xk in cleartext
	matrixMult(Axk, smd->AA, smd->xx, smd->n, smd->n, smd->n, 1); 
	// 2nd part: Multiply B and u in cleartext			
	matrixMult(Bu, smd->BB, u, smd->n, smd->m, smd->m, 1);	
	// 3rd part: Add each and every part of the measurement 
	// Achieve the xNoise vector w.r.t. number of iterations	
	for (size_t i = 0; i < smrp->n; i++)
		xSpecNoise[i][0] = smd->xNoise[numIter][i];
	// Add the the all three terms and assign it simulation data struct
	matrixAdditionThree(AxkBuNoise, Axk, Bu, xSpecNoise, smd->n, 1);
	assignXResult(smd, AxkBuNoise, smd->n, 1);
} 

/*
	The function for extracting the experimental results
*/
void extractExpRes(Ciphertext * res_CP, struct simulationMatrixMRP * smrp, struct simulationMatrixData * smd, int numIter, int numOfRows, double scale, SEALContext *contextPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, char * matname){

	// Print the number of iterations
	cout << "Iteration: " << numIter << endl;	
	cout << "Experimented Result vector name: " << matname << endl;	
	
	// Decrypt the result 
	Plaintext Res_PL;	
	vector<double> Res_Vec;
	decryptorPtr->decrypt(*res_CP, Res_PL);
	encoderPtr->decode(Res_PL, Res_Vec);	
	
	// Assign the content of the respective matrix to the simulation data struct 	
	for(int i = 0; i < numOfRows; i++){
		double i_Ind_Res = Res_Vec[i * smrp->N];			
		if(strcmp(matname, "xe") == 0)
			smd->xe_Res[numIter][i] = i_Ind_Res;  	
		if(strcmp(matname, "u") == 0)
			smd->u_Res[numIter][i]  = i_Ind_Res;
		if(strcmp(matname, "xp") == 0)
			smd->xp_Res[numIter][i] = i_Ind_Res;		
		if(strcmp(matname, "residue") == 0)
			smd->residue_Res[numIter][i] = i_Ind_Res;
		if(strcmp(matname, "sBar") == 0)
			smd->sBar_Res[numIter][i] = i_Ind_Res; 	  	
		if(strcmp(matname, "indInp") == 0)
			smd->indInp_Res[numIter][i] = i_Ind_Res;
		if(strcmp(matname, "s") == 0)
			smd->s_Res[numIter][i] = i_Ind_Res; 		  
		if(strcmp(matname, "x") == 0)
			smd->x_Res[numIter][i] = i_Ind_Res;		  
		if(strcmp(matname, "y") == 0)
			smd->y_Res[numIter][i] = i_Ind_Res;
		if(strcmp(matname, "alarm") == 0) {
			double roundedAlarmValue = round(i_Ind_Res); 	
			if(roundedAlarmValue > 0.4)
				smd->alarm_Res[numIter][i] = 1;
			else 
				smd->alarm_Res[numIter][i] = 0;
			printf("%f ", smd->alarm_Res[numIter][i]);
		}
		if(strcmp(matname, "alarm") != 0) 
			printf("%.11f ", i_Ind_Res);
	}	
	printf("\n");		
}
