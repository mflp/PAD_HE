/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which is used to model the online and offline application phases of the Private Anomaly Detection Application
   						  The main assumption is that the cloud (or the controller) knows the model, therefore, we should hold each matrix in plaintext.    	 
   						  The row-matrix packing format is used. 
   						  The application is implemented with the major modifications on the encoders.cpp example of SEAL library version 3.6.2    
	Note				: Commenting and cleaning up has been on July 2022
   Abbreviation/Acronym	: # _RPL: Raw plain 
						  # _MRP, _mrp: Matrix row packing format
*/

/* Define  several constants */
#define PI 3.141592653589793
#define numRandBits 20

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
#include "encodersplain.h"        // The main application class of the crypto application

/* Import the important selected C libraries*/
#include <iostream>
#include <array>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <stdbool.h>
#include <cmath>
#include <stdio.h>

/* Call main namespaces */
using namespace std;
using namespace seal;

/*
	Define and initialize the performance/time measurement parameters for each operation in the control system  
*/
double encParamInitTotalTime  = 0; 
double encPlMatrVectTotalTime = 0;
double SensorMeasurementTime  = 0;

double decryptionUTime  	  = 0;
double encryptionYTime  	  = 0;

double Eq2Time   			 = 0;
double Eq3Time   			 = 0;
double Eq4_5Time 			 = 0;
double Eq6Time   			 = 0;
double CUSUMTime 			 = 0; 
double Eq10Time  			 = 0; 
double EqLastSecretShareTime = 0;

/*
	The function for performing the cyberphysical system's functionalities (one after another)   
*/
void performMultipleMatrixVectorMultiplicationsPlain(EncryptionParameters *parmsPtr, 
SEALContext *contextPtr, auto *secret_keyPtr, PublicKey *public_keyPtr, RelinKeys *relin_keysPtr, 
GaloisKeys *gal_keysPtr, Encryptor *encryptorPtr, Evaluator *evaluatorPtr, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, 
FILE *fp,  double scale,  size_t numberOfIterations){
  
	/* 
	** =============================================================== 
	** ===============================================================
	** ==================== SYSTEM INITIALIZATION ==================== 
	** ===============================================================
	** =============================================================== 
	*/
	// Start of the initialization
	clock_t st_Crypto_Matrix_Vec_Init_Start = clock(); // # Ciphertext Initialization-Start #	
	// Define the all the simulation file paths based on the matrix dimensions 
	char * folderPath_y10_u2 	= "/home/mfcustben/SEAL_1/native/examples/ZZZ_Final_HE_Correctness_Issue/Luis_Very_Last_25_03_2022/all_data/y10_u2";
	char * folderPath_y20_u4 	= "/home/mfcustben/SEAL_1/native/examples/ZZZ_Final_HE_Correctness_Issue/Luis_Very_Last_25_03_2022/all_data/y20_u4";
	char * folderPath_y50_u10 	= "/home/mfcustben/SEAL_1/native/examples/ZZZ_Final_HE_Correctness_Issue/Luis_Very_Last_25_03_2022/all_data/y50_u10";		
	// Create and build simulationMatrixData struct holding the  	
	struct simulationMatrixData * smd = (struct simulationMatrixData *) calloc(1, sizeof(struct simulationMatrixData)); 	
	create_SimulationMatrixData(smd);
	// Read the recorded system matrices and vectors
	// assignValMatrixDatabyFileRead(smd, folderPath_y10_u2); // The system matrices and vectors when n = 10, m = 2  	
	// assignValMatrixDatabyFileRead(smd, folderPath_y20_u4); // The system matrices and vectors when n = 20, m = 4
	assignValMatrixDatabyFileRead(smd, folderPath_y50_u10); // The system matrices and vectors when n = 50, m = 10
	// Initialize the remaining vectors including the sensor measurement vector y
	initRemainVec(smd);	
	// Initialize the CUSUM Chebyshev Appx. parameters and arrays
	assignCUSUMChebyshevAppxParams(smd);  
	// Convert simulationMatrixData to vector in the proper Matrix-Row-Packing (MRP) format
	struct simulationMatrixMRP *smrp   = (struct simulationMatrixMRP *) calloc(1, sizeof(struct simulationMatrixMRP));
	create_SimulationMatrixDataMRP(smd, smrp);
	assignValMatrixDataMRP(smd, smrp);	


	// ================ ================ ================ ================ ================
	// ========= Estimation Phase-Equation-2 Plaintext-Ciphertext Initialization ========== 
	// ================ ================ ================ ================ ================	
	// Equation-2 plaintext and ciphertext pointers initialization
	// Create plaintext and ciphertext pointers   
	Plaintext PL_GAMMA, PL_LL; 
	Plaintext * GAMMA_PL,  * LL_PL;
	Ciphertext CP_xexe, CP_xGxG; 	
	Ciphertext * xexe_CP, * xGxG_CP; 
	// Initialize the pointers 
	GAMMA_PL = &PL_GAMMA;
	LL_PL	 = &PL_LL;
	xGxG_CP	 = &CP_xGxG;
	xexe_CP  = &CP_xexe;
	// Encode and encrypt the plaintexts and ciphertexts respectively
	makePlaintextMatRowPacking(smrp->GAMMA_MRP, scale, GAMMA_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->LL_MRP, scale, LL_PL, encoderPtr);	
	encryptXVectorMatRowPacking(smrp->xGxG_MRP, scale, xGxG_CP, encryptorPtr, encoderPtr);	
	encryptXVectorMatRowPacking(smrp->xexe_MRP, scale, xexe_CP, encryptorPtr, encoderPtr);	
		
	// ================ ================ ================ ================ ================
	// ================== Equation-3 Plaintext-Ciphertext Initialization ================== 
	// ================ ================ ================ ================ ================	
	// Create plaintext and ciphertext pointers  
	Plaintext PL_KGKG, PL_KLKL, PL_KxKx; // PL_KMinus; 
	Plaintext *KGKG_PL, *KLKL_PL, *KxKx_PL; // *KMinus_PL;
	Ciphertext CP_KxuGKxuG, CP_uGuG_AS; // , CP_uGuG; 
	Ciphertext * KxuGKxuG_CP, *uGuG_AS_CP; // , *uGuG_CP; 
	// Initialize the pointers 
	KGKG_PL 	= &PL_KGKG;
	KLKL_PL 	= &PL_KLKL;
	KxuGKxuG_CP	= &CP_KxuGKxuG;
	KxKx_PL 	= &PL_KxKx;  	
	uGuG_AS_CP 	= &CP_uGuG_AS; 
	// Encode and encrypt the plaintexts and ciphertexts respectively
	makePlaintextMatRowPacking(smrp->KGKG_MRP, scale, KGKG_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->KLKL_MRP, scale, KLKL_PL, encoderPtr);	
	makePlaintextMatRowPacking(smrp->KxKx_MRP, scale, KxKx_PL, encoderPtr); // Used in the very first iteration 	 
	encryptXVectorMatRowPacking(smrp->KxuGKxuG_MRP, scale, KxuGKxuG_CP, encryptorPtr, encoderPtr);
	encryptXVectorMatRowPacking(smrp->uGuG_AS_MRP, scale, uGuG_AS_CP, encryptorPtr, encoderPtr); // Used in the very first iteration  
		
	// ================ ================ ================ ================ ================
	// ================== Equation-4-5 Plaintext-Ciphertext Initialization ================
	// ================ ================ ================ ================ ================	
	// Create plaintexts and ciphertext pointers
	Plaintext PL_ACL, PL_BB; 
	Plaintext *ACL_PL, *BB_PL;
	Ciphertext CP_uGuG; 
	Ciphertext *uGuG_CP; 
	// Initialize the pointers 
	ACL_PL = &PL_ACL;   
	BB_PL  = &PL_BB; 	
	uGuG_CP = &CP_uGuG;
	// Encode and encrypt the plaintexts and ciphertexts respectively
	makePlaintextMatRowPacking(smrp->ACL_MRP, scale, ACL_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->BB_MRP, scale, BB_PL, encoderPtr);
	encryptXVectorMatRowPacking(smrp->uGuG_MRP, scale, uGuG_CP, encryptorPtr, encoderPtr);	

	// ================ ================ ================ ================ ================
	// ================== Equation-6 Plaintext Initialization ============================= 
	// ================ ================ ================ ================ ================	
	// Define plaintext pointers
	Plaintext PL_xpxp; 
	Plaintext *xpxp_PL = &PL_xpxp;
	// Encode the plaintexts	
	makePlaintextMatRowPacking(smrp->ACL_MRP, scale, ACL_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->xpxp_MRP, scale, xpxp_PL, encoderPtr);
		
	// ================ ================ ================ ================ ================
	// ================== Equation-8 Plaintext-Ciphertext Initialization ================== 
	// ================ ================ ================ ================ ================	
	// Create plaintext and ciphertext pointers
	Plaintext PL_ss, PL_vv, PL_alpha_Eq8, PL_beta_Eq8; 
	Plaintext *ss_PL, *vv_PL, *alpha_Eq8_PL, *beta_Eq8_PL; 
	Ciphertext CP_ss; 	
	Ciphertext * ss_CP; 
	// Assign Chebyshev Approximation degrees
	int chebDegMax  = 12;
	int chebDegISub = 12;	
	// Create the remaining plaintext (vector) pointers   
	Plaintext PL_firstChebPowSerTerm_max;
	Plaintext *firstChebPowSerTerm_max_PL;
	vector <Plaintext> * chebPwSrCoefVec_Eq8_PL =  new vector<Plaintext>(smd->chebDegEq8 + 1);
	// Initialize the pointers 
	ss_PL 		 = &PL_ss;   
	vv_PL 		 = &PL_vv;
	alpha_Eq8_PL = &PL_alpha_Eq8;
	beta_Eq8_PL  = &PL_beta_Eq8; 
	ss_CP 		 = &CP_ss; 
	firstChebPowSerTerm_max_PL = & PL_firstChebPowSerTerm_max;	
	// Encode and encrypt the plaintexts and ciphertexts respectively
	makePlaintextMatRowPacking(smrp->ss_MRP, scale, ss_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->vv_MRP, scale, vv_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->alpEq8_MRP, scale, alpha_Eq8_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->betEq8_MRP, scale, beta_Eq8_PL, encoderPtr);
	encryptXVectorMatRowPacking(smrp->ss_MRP, scale, ss_CP, encryptorPtr, encoderPtr);  	 			  	
	makePlaintextMatRowPacking(smrp->chebPowSerFT_Eq8_MRP, scale, firstChebPowSerTerm_max_PL, encoderPtr); 	 			
	preparePwSrCoeffVec(chebPwSrCoefVec_Eq8_PL, smd, smd->chebDegEq8 + 1, smrp->N, scale, decryptorPtr, encoderPtr, true);

	// ================ ================ ================ ================ ================
	// ================== Equation-9 Plaintext-Ciphertext initialization ================== 
	// ================ ================ ================ ================ ================	

	// Create plaintext and ciphertext pointers
	Plaintext PL_TAU, PL_alpha_Eq9, PL_beta_Eq9;
	Plaintext * TAU_PL, *alpha_Eq9_PL, *beta_Eq9_PL;
 	// Initialize the pointers 
	TAU_PL 		 = &PL_TAU;  
	alpha_Eq9_PL = &PL_alpha_Eq9;
	beta_Eq9_PL  = &PL_beta_Eq9;
	// Create the first term   
	Plaintext PL_firstChebPowSerTerm_ISub;
	Plaintext * firstChebPowSerTerm_ISub_PL;  
	// Initialize the first term  pointer	
	firstChebPowSerTerm_ISub_PL = &PL_firstChebPowSerTerm_ISub;
	// Initialize the Chenbyshev coefficient array 
	vector <Plaintext> * chebPwSrCoefVec_Eq9_PL =  new vector<Plaintext>(smd->chebDegEq9 + 1);	
	// Encode and encrypt the plaintexts and ciphertexts respectively
 	makePlaintextMatRowPacking(smrp->TAU_MRP, scale, TAU_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->alpEq9_MRP, scale, alpha_Eq9_PL, encoderPtr); // Uncomment after rawplain.cpp and prepareVecMatMRP.cpp modifications
	makePlaintextMatRowPacking(smrp->betEq9_MRP, scale, beta_Eq9_PL, encoderPtr);  // Uncomment after rawplain.cpp and prepareVecMatMRP.cpp modifications
	makePlaintextMatRowPacking(smrp->chebPowSerFT_Eq9_MRP, scale, firstChebPowSerTerm_ISub_PL, encoderPtr); // Prepare the very first term of Cheb. Polyn. for Eq-8   		
	preparePwSrCoeffVec(chebPwSrCoefVec_Eq9_PL, smd, smd->chebDegEq9 + 1, smrp->N, scale, decryptorPtr, encoderPtr, false); // Prepare coeff vector of Cheb. Polyn. for Eq-9
	// Last multiplication of Alarm 	
	Ciphertext CP_vecOne_last;
	Ciphertext * vecOne_last_CP;
	vecOne_last_CP  = &CP_vecOne_last;
	// makePlaintextMatRowPacking(smrp->One_MRP, scale, vecOne_PL_last, encoderPtr);
	encryptXVectorMatRowPacking(smrp->One_MRP, scale, vecOne_last_CP, encryptorPtr, encoderPtr);
	
	// End of the initialization	
	clock_t st_Crypto_Matrix_Vec_Init_End 	= clock(); // # Ciphertext Initialization-End #
	// Calculate The Initialization Time  
	double cryp_Mat_Vect_Init_Time   		= (double) (st_Crypto_Matrix_Vec_Init_End - st_Crypto_Matrix_Vec_Init_Start) / CLOCKS_PER_SEC;
	encPlMatrVectTotalTime  			   += cryp_Mat_Vect_Init_Time;
	printf("Crypto Matr-Vect Init Total Time-Indiv. Measur.: %f-%f\n", encPlMatrVectTotalTime, cryp_Mat_Vect_Init_Time);	
	// exit(0); // 2nd check point
			
	/* 
	** ====================================================================== 
	** ======================================================================
	** ==================== CRYPTO ONLINE ITERATION LOOP ==================== 
	** ======================================================================
	** ====================================================================== 
	*/
	bool isFirstIter = true;		
	for(int k = 0; k < smd->tMax; k++){
 		// Print the online iteration banner for time measurement tracking  
		cout << "=====================" << endl;		
		cout << "Online Iteration: " << k + 1 << endl;   		
		cout << "=====================" << endl;		

		// Sensor measurement (System- Physical Plant) - Cleartext Part     	
		// Receive the encrypted noisy sensor measurement
		clock_t start_SensMeasurement = clock(); 
		Ciphertext CP_yy, CP_yyAS; 
		Ciphertext *yy_CP, *yyAS_CP; 	
		yy_CP 	 = &CP_yy;	
		yyAS_CP  = &CP_yyAS;		
		sense_Encrypt_y(smd, smrp, k, yy_CP, yyAS_CP, scale, contextPtr, encryptorPtr, encoderPtr); // Fulfill the sensing duty of client		  
		clock_t end_SensMeasurement = clock();	 	
	
		// Compute Sensor Measurement Phase Duration 
		double sens_Meas_Indv  = (double) (end_SensMeasurement - start_SensMeasurement) / CLOCKS_PER_SEC;
		SensorMeasurementTime += sens_Meas_Indv;	
		encryptionYTime  	  += sens_Meas_Indv;

		// Define ciphertexts for the results of the crypto-estimation, control action, prediction, residues functions   
		Ciphertext CP_secEqRes, CP_thirdEqRes, CP_fourthfifthEqRes, CP_sixthEqRes;		
		Ciphertext * secEqRes_CP, * thirdEqRes_CP, * fourthfifthEqRes_CP, * sixthEqRes_CP;  
		secEqRes_CP   = &CP_secEqRes; 		
		thirdEqRes_CP = &CP_thirdEqRes; 		
		fourthfifthEqRes_CP = &CP_fourthfifthEqRes; 
		sixthEqRes_CP = & CP_sixthEqRes; 

		// Define the time measurement variables for each separate phase  
		double Eq2_Meas_Indv 	= 0;
		double Eq3_Meas_Indv 	= 0; 
		double Eq4_5_Meas_Indv 	= 0;
		double Eq6_Meas_Indv 	= 0;
	
		// ================ ================ ================ ================ ================
		// ================== The iterations when k > 0 (after the very first iteration) ====== 
		// ================ ================ ================ ================ ================		
		if(k > 0){ // or if (k >= 1)
			isFirstIter = false;
			// Estimation Phase (Equation-2) 
			clock_t start_2ndEquation = clock();
			applyEquation_2_PLCP(secEqRes_CP, smrp, GAMMA_PL,  LL_PL, xGxG_CP, yy_CP, xexe_CP, scale, contextPtr, encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, relin_keysPtr, gal_keysPtr);
			clock_t end_2ndEquation = clock();	
			// Compute Estimate Phase Duration 
			Eq2_Meas_Indv  = (double) (end_2ndEquation - start_2ndEquation) / CLOCKS_PER_SEC;
			Eq2Time    += Eq2_Meas_Indv;
			printf("Estimation x^e- Eq2 Measur. Time-Indiv. Measur.: %f-%f\n", Eq2Time, Eq2_Meas_Indv);	
			
			// Control Action Phase (Equation-3)  
			clock_t start_3rdEquation = clock();			
			applyEquation_3_PLCP(thirdEqRes_CP, smrp, KGKG_PL, KLKL_PL, KxuGKxuG_CP, yy_CP, xexe_CP, scale, contextPtr, encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, relin_keysPtr, gal_keysPtr);
			clock_t end_3rdEquation = clock();
			// Compute Control Action Phase Duration  
			Eq3_Meas_Indv  = (double) (end_3rdEquation - start_3rdEquation) / CLOCKS_PER_SEC;
			Eq3Time    += Eq3_Meas_Indv;
			printf("Control u- Eq3 Measur. Time-Indiv. Measur.: %f-%f\n", Eq3Time, Eq3_Meas_Indv);	

			// Prediction Phase (Equation-4-5) 
			clock_t start_4_5_th_Equation = clock();
			applyEquation_4_5_PLCP(fourthfifthEqRes_CP, smrp, ACL_PL, BB_PL, uGuG_CP, xexe_CP, scale, contextPtr, encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, relin_keysPtr, gal_keysPtr);
			clock_t end_4_5_th_Equation = clock();
			// Compute Prediction Phase Duration
			Eq4_5_Meas_Indv  = (double) (end_4_5_th_Equation - start_4_5_th_Equation) / CLOCKS_PER_SEC;
			Eq4_5Time    += Eq4_5_Meas_Indv;
			printf("Prediction x^p- Eq4_5 Measur. Time-Indiv. Measur.: %f-%f\n", Eq4_5Time, Eq4_5_Meas_Indv);	
			// exit(0);

			// Residues Phase (Equation-6) 
			clock_t start_6_th_Equation = clock();		
			applyEquation_6_PLCP(sixthEqRes_CP, smrp, fourthfifthEqRes_CP, yyAS_CP, scale, contextPtr, encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, relin_keysPtr, gal_keysPtr); 
			clock_t end_6_th_Equation = clock();	
			// Compute Residues Phase Duration 
			Eq6_Meas_Indv  = (double) (end_6_th_Equation - start_6_th_Equation) / CLOCKS_PER_SEC;
			Eq6Time    += Eq6_Meas_Indv;
			printf("Residues r- Eq6 Measur. Time-Indiv. Measur.: %f-%f\n", Eq6Time, Eq6_Meas_Indv);	
 		} 

		// ================ ================ ================ ================ ================
		// ================== The iterations when k = 0 (the very first iteration) ============ 
		// ================ ================ ================ ================ ================		
		else {
			isFirstIter = true;

			// Control Action Phase (Equation-3)
			clock_t start_3rdEquation = clock();		
			applyEquation_3_fiter_PLCP(thirdEqRes_CP, smrp, KxKx_PL, uGuG_AS_CP, scale, contextPtr, encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, relin_keysPtr, gal_keysPtr);
			clock_t end_3rdEquation = clock();
			// Compute Control Action Phase Duration
			Eq3_Meas_Indv  = (double) (end_3rdEquation - start_3rdEquation) / CLOCKS_PER_SEC;
			Eq3Time    += Eq3_Meas_Indv;
			printf("Control Action u - Eq3 Measur. Time-Indiv. Measur.: %f-%f\n", Eq3Time, Eq3_Meas_Indv);	

			// Residues Phase (Equation-6)
			clock_t start_6_th_Equation = clock();		
			applyEquation_6_fiter_PLCP(sixthEqRes_CP, smrp, xpxp_PL, yyAS_CP, scale, contextPtr, encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, relin_keysPtr, gal_keysPtr); 
			clock_t end_6_th_Equation = clock();	
			// Compute Residues Phase Duration
			Eq6_Meas_Indv  = (double) (end_6_th_Equation - start_6_th_Equation) / CLOCKS_PER_SEC;
			Eq6Time    += Eq6_Meas_Indv;
			printf("Residues r- Eq6 Measur. Time-Indiv. Measur.: %f-%f\n", Eq6Time, Eq6_Meas_Indv);	
		}

		/* 
		** ====================================================================== 
		** ======================================================================
		** ======================= CUSUM-ALARM COMPUTATION ====================== 
		** ======================================================================
		** ====================================================================== 
		*/
		// Initial-CUSUM-Alarm Phase (Equation-8-9)	
		clock_t start_8_9_th_Equation = clock();				
		Plaintext PL_vecOne; 
		Plaintext *vecOne_PL = & PL_vecOne;
		makePlaintextMatRowPacking(smrp->One_MRP, scale, vecOne_PL, encoderPtr);
		// Apply Chebyshev Appx. for both Equation-8-9      
		Ciphertext CP_eigthEqRes, CP_ninthEqRes;
		Ciphertext * eigthEqRes_CP = &CP_eigthEqRes;		
		Ciphertext * ninthEqRes_CP = &CP_ninthEqRes;		 	
		applyEquation_CUSUM_PLCP(eigthEqRes_CP, ninthEqRes_CP, smrp, smd, k, ss_PL, vv_PL, 
								 TAU_PL, alpha_Eq8_PL, beta_Eq8_PL, alpha_Eq9_PL, beta_Eq9_PL,
								 vecOne_PL, 
								 firstChebPowSerTerm_max_PL, chebPwSrCoefVec_Eq8_PL,
								 firstChebPowSerTerm_ISub_PL, chebPwSrCoefVec_Eq9_PL,
								 sixthEqRes_CP, ss_CP, scale, 
								 contextPtr, encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, relin_keysPtr, gal_keysPtr, isFirstIter);	
		clock_t end_8_9_th_Equation = clock();		
		// Compute CUSUM-Alarm Phase Duration 
		double CUSUM_Meas_Indv  = (double) (end_8_9_th_Equation - start_8_9_th_Equation) / CLOCKS_PER_SEC;
		CUSUMTime   += CUSUM_Meas_Indv;
		printf("CUSUM Time-Indiv. Measur.: %f-%f\n", CUSUMTime, CUSUM_Meas_Indv);
	
		// CUSUM-Parametric Sum Phase 
		clock_t start_10_th_Equation = clock();	
		applyEquation_10_PLCP(ss_CP, smrp, eigthEqRes_CP, ninthEqRes_CP, vecOne_last_CP, scale, contextPtr, encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, relin_keysPtr, gal_keysPtr);
		clock_t end_10_th_Equation = clock();
		// Compute CUSUM-Parametric Sum Phase Duration 
		double Eq10_Meas_Indv  = (double) (end_10_th_Equation - start_10_th_Equation) / CLOCKS_PER_SEC;
		Eq10Time   += Eq10_Meas_Indv;
		printf("Alarm Computation Eq10 Measur. Time-Indiv. Measur.: %f-%f\n", Eq10Time, Eq10_Meas_Indv);	

		/* 
		** ====================================================================== 
		** ======================================================================
		** =========== END-OF-LOOP-SECRET SHARE COMPUTATION ===================== 
		** ======================================================================
		** ====================================================================== 
		*/
		// Secret-share Phase (for next iteration) 
		clock_t start_Last_Secret_Share = clock(); 				
		if(k > 0){
			secretShareEstimation(secEqRes_CP, contextPtr,  encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, scale, smrp->n, smrp->n,  smrp->N);
			*xexe_CP = *secEqRes_CP;
		}
		clock_t end_Last_Secret_Share = clock();	
		// Compute Secret Share Phase Duration 
		double LastSecretShare_Meas_Indv  = (double) (end_Last_Secret_Share - start_Last_Secret_Share) / CLOCKS_PER_SEC;
		EqLastSecretShareTime    		 += LastSecretShare_Meas_Indv;
		printf("Last Secret Share of x^e Measur. Time-Indiv. Measur.: %f-%f\n", EqLastSecretShareTime, LastSecretShare_Meas_Indv);	

		// Apply process noise phase 	
		clock_t start_Process_Noise_Addition = clock();
		applyXVecNoiseAddition(thirdEqRes_CP, k, smrp, smd, scale, contextPtr, encryptorPtr, decryptorPtr, encoderPtr);	
		clock_t end_Process_Noise_Addition  = clock();
		// Compute process noise phase duration 
		double Process_Noise_Addition_Meas_Indv  = (double) (end_Process_Noise_Addition - start_Process_Noise_Addition)/CLOCKS_PER_SEC;
		SensorMeasurementTime    		 += Process_Noise_Addition_Meas_Indv;
		decryptionUTime += Process_Noise_Addition_Meas_Indv;
		

		// =========================== =========================== ===============
		// =========== PRINT TIME MEASUREMENT COMPUTATIONS =======================  
	 	// =========================== =========================== =============== 		 	
		printf("Sensor Measurement Time-Indiv. Measur.: %f-Enc-%f-Proc-Noise-%f\n", SensorMeasurementTime, sens_Meas_Indv, Process_Noise_Addition_Meas_Indv);
		cout << "Iteration" << k << " Individual Times " << endl;
		if (k == 0)
			cout  <<  sens_Meas_Indv << "," << 0 << "," << Eq3_Meas_Indv << "," << 0 << "," << Eq6_Meas_Indv << "," << CUSUM_Meas_Indv << "," << Eq10_Meas_Indv << "," << LastSecretShare_Meas_Indv << "," <<	Process_Noise_Addition_Meas_Indv << endl;
		else 
			cout <<  sens_Meas_Indv << "," << Eq2_Meas_Indv << "," << Eq3_Meas_Indv << "," << Eq4_5_Meas_Indv << "," << Eq6_Meas_Indv << "," << CUSUM_Meas_Indv << "," << Eq10_Meas_Indv << "," << LastSecretShare_Meas_Indv << "," <<	Process_Noise_Addition_Meas_Indv << endl;
		cout << "Iteration" << k << " Total Times" << endl;
		cout << SensorMeasurementTime << "," << encryptionYTime << "," << decryptionUTime << "," << Eq2Time << "," << Eq3Time <<  "," <<  Eq4_5Time << "," << Eq6Time <<  "," <<  CUSUMTime <<  "," << Eq10Time << "," << EqLastSecretShareTime << endl;   

		
		// =========================== =========================== ==========================
		// ============= DATA EXTRACTION AND EXPERIMENTAL RECORDING-REPORTING =============== 
	 	// =========================== =========================== ========================== 
		// ## Four targets: y, u, alarm, s 	
		// /*						
		extractExpRes(yyAS_CP, smrp, smd, k, smd->n, scale, contextPtr, decryptorPtr, encoderPtr, "y"); // Sensor Measurement
		// Note: Changes based on the first iteration	
		// extractExpRes(secEqRes_CP, smrp, smd, k, smd->n, scale, contextPtr, decryptorPtr, encoderPtr, "xe"); // Estimation		
		extractExpRes(thirdEqRes_CP, smrp, smd, k, smd->n, scale, contextPtr, decryptorPtr, encoderPtr, "u"); // Control 		
		// Note: Changes based on the first iteration			
		// extractExpRes(fourthfifthEqRes_CP, smrp, smd, k, smd->n, scale, contextPtr, decryptorPtr, encoderPtr, "xp");		
		// extractExpRes(sixthEqRes_CP, smrp, smd, k, smd->n, scale, contextPtr, decryptorPtr, encoderPtr, "residue"); // Residue		
		// extractExpRes(eigthEqRes_CP, smrp, smd, k + 1, smd->n, scale, contextPtr, decryptorPtr, encoderPtr, "alarm"); // RELU Appx. Function  
		extractExpRes(ninthEqRes_CP, smrp, smd, k, smd->n, scale, contextPtr, decryptorPtr, encoderPtr, "alarm"); // Alarm Appx. Func. (Indicator funct)		
		extractExpRes(ss_CP, smrp, smd, k + 1, smd->n, scale, contextPtr, decryptorPtr, encoderPtr, "s"); // CUSUM parametric SUM
		// */
	}
}

/*
	The function for setting the configurations of the SEAL crypto application and calling the main system function   
*/
void ckks_encoder_modify_matrix_row_packing_functional()
{
	// Print the introduction banner 
	print_example_banner("Example: Encoders / CKKS Encoder");
	
	/* 
	** ====================================================================== 
	** ======================================================================
	** ======================= CRYPTO-PARAMETER CONFIGURATON ================ 
	** ======================================================================
	** ====================================================================== 
	*/

	// Produce the timestamp for the start of the encryption parameters initialization  			 			
	clock_t st_EncPar_Init_Start = clock();

	// Create encryption parameters for the CKKS scheme			
	EncryptionParameters *parmsPtr;
	EncryptionParameters parms(scheme_type::ckks);
	parmsPtr = &parms;
	size_t poly_modulus_degree = 32768; // Determine the polynomial modulus degree		
	cout << "Maximum number of bits: "  << CoeffModulus::MaxBitCount(poly_modulus_degree) << "\n"; // Print maximum bit count for the given poly_modulus degree
	size_t bitsizesparam = 50; // Determine bitsize parameter w.r.t. the polynomial modulus degree     
	parms.set_poly_modulus_degree(poly_modulus_degree);	
	// parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, bitsizesparam, 60})); // # 1
	// parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, bitsizesparam, bitsizesparam, 60})); // # 2
	// parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, bitsizesparam, bitsizesparam, bitsizesparam, 60})); // # 3 
	// parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, 60})); // # 4 
	// parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, 60})); // # 5 
	// parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, 60})); // # 6 
	// parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, 60})); // # 8
	parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, bitsizesparam, 60})); // # 9
	
	// Create the SEAL Context and print the context parameters.	
	SEALContext *contextPtr;	
	SEALContext context(parms);
	contextPtr = &context; 

	print_parameters(context);
	cout << endl;
	print_line(__LINE__);
	
	// 	Create SEAL library homomorphic encryption keys
	KeyGenerator keygen(context);
	auto secret_key = keygen.secret_key();
	auto *secret_keyPtr = &secret_key; 
	
	PublicKey public_key;
	keygen.create_public_key(public_key);
	PublicKey *public_keyPtr;
	public_keyPtr = &public_key;

	RelinKeys relin_keys;
	keygen.create_relin_keys(relin_keys);
	RelinKeys *relin_keysPtr;
	relin_keysPtr = &relin_keys;

	GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
	GaloisKeys *gal_keysPtr;
    gal_keysPtr = &gal_keys; 

	// Set up an Encryptor, Evaluator, and Decryptor
	Encryptor *encryptorPtr;
	Encryptor encryptor(context, public_key);
	encryptorPtr = &encryptor;

	Evaluator *evaluatorPtr;
	Evaluator evaluator(context);
	evaluatorPtr = &evaluator; 

	Decryptor *decryptorPtr;
	Decryptor decryptor(context, secret_key);
	decryptorPtr = &decryptor;

	// Create a special encoder to create CKKS plaintexts (i.e., the CKKSEncoder encodes 
    // vectors of real or complex numbers into Plaintext objects, which can subsequently be 
    // encrypted.) 
	CKKSEncoder *encoderPtr; 
	CKKSEncoder encoder(context);
	encoderPtr = &encoder;

	// Determine the parameter scale by which the floating-point coefficients of `input will be scaled up 	
	double scale = pow(2.0, bitsizesparam);

	// In CKKS, the number of slots is poly_modulus_degree/2 and each slot encodes one real or complex number.
	size_t slot_count 	= encoder.slot_count();
	size_t row_size 	= slot_count / 2;
	cout << "Number of slots: " << slot_count << endl; // Print the number of slots
	
	// 	Open the file for recording the matrix computation errors   
	FILE *fp;
   	fp = fopen("./DifferenceResults_Trial.txt", "w");

	// 	Print the beginner app for tracking the application 
	printAppBeginner();
	
	// Produce the timestamp for the end of the encryption parameters initialization  			
	clock_t st_EncPar_Init_End = clock();
	
	// Compute the crypto-parameter configuration phase duration
	double cryp_Param_Init_Time  = (double) (st_EncPar_Init_End - st_EncPar_Init_Start) / CLOCKS_PER_SEC;
	encParamInitTotalTime += cryp_Param_Init_Time; 
	printf("Crypto Param Init Total Time-Indiv. Measur.: %f-%f\n", encParamInitTotalTime, cryp_Param_Init_Time);
	
	// Perform the multiple matrix-vector multiplications (i.e., Call the function in which the main functioning of the target cyberphysical system is implemented) 	
	size_t numberOfIterations = 2; 
	performMultipleMatrixVectorMultiplicationsPlain(parmsPtr, contextPtr, secret_keyPtr, public_keyPtr, relin_keysPtr, gal_keysPtr, encryptorPtr, evaluatorPtr, decryptorPtr, encoderPtr, fp, scale, numberOfIterations);
	
	// Close the file for recording the matrix computation errors   	
	fclose(fp);		
}	

/* Main function of the SEAL crypto application */
void example_encoders()
{
    print_example_banner("Example: Encoders");
	for (int i = 0; i < 1; i++)
		ckks_encoder_modify_matrix_row_packing_functional();	
}

