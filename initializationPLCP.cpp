/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which is used to initialize the simulation data for the Plaintext-Ciphertext setting 
	Note				: Commenting and cleaning up has been on July 2022 	 
*/

/* Import the other classes' header files */
#include "examples.h" 			  // 0	- The class for the console menu and guiding the user to the preferred application     	
#include "rawplain.h"             // 1   - The class which reads and stores the plain matrices        
#include "generateplaintextMRP.h" // 2.0 - The functions used to prepare plain matrices in MRP format
#include "prepareVecMatMRP.h"     // 2.1 - The class which prepares the read matrices in MRP format  		
#include "initializationPLCP.h"   // 3   - The class containing the functions which convert the prepared matrices to Plaintext and Ciphertext objects.            
#include "plcpOperations.h"       // 4.0 - The class containing the functions which does Ciphertext-Ciphertext and Plaintext-Ciphertext arithmetic and algebraic operations  	    
#include "encryptedAppx.h"        // 4.1 - The class containing the functions which does the Chebyshev approximation with different assumptions     
#include "secretShare.h" 		  // 4.2 - The class containing the functions which does the secret sharing  
#include "applyPLCPSimulation.h"  // 5   - The class containing the crypto application functions for each targeted equation       

/* Import the important selected C libraries */
#include <iostream>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <unistd.h>

/* Call the important namespaces */
using namespace std;
using namespace seal;

/*
 	The function for preparing the power series coefficient vector in plaintext 
*/
void preparePwSrCoeffVec(vector<Plaintext> * powerSeriesCoeffVec_PL, struct simulationMatrixData *smd, size_t chebyshevDegree, size_t calcRowSize, double scale, Decryptor *decryptorPtr, CKKSEncoder *encoderPtr, bool isMaxFunc){

	for (size_t i = 1; i < chebyshevDegree; i++){
		double powSerCoeff;		
		if (isMaxFunc) 
		 powSerCoeff= smd->eq8maxAppx_PS_Coeff_D12_y_10_u_2[i][0];
		else 	
		 powSerCoeff= smd->eq9ISubAppx_PS_Coeff_D12_y_10_u_2[i][0];
	
		Plaintext powSerCoeff_PL;
		encoderPtr->encode(powSerCoeff, scale, powSerCoeff_PL);
		powerSeriesCoeffVec_PL->at(i) = powSerCoeff_PL;						
	}
}

/*
 	The function for initializing the simulation data solely for the Equation-2 in the Plaintext-Ciphertext setting
*/
void initSimDatPLCP(struct simulationMatrixPLCP *splcp, struct simulationMatrixMRP *smrp, double scale, Decryptor *decryptorPtr, Encryptor *encryptorPtr, CKKSEncoder *encoderPtr){
	cout << "initSimDatPLCP-In" << "\n"; 
	
   // Define equation-2 pointers 
	Plaintext PL_xexe;
	Plaintext PL_GAMMA;
	Plaintext PL_LL;
	Plaintext PL_xGxG;
	Ciphertext CP_xexe; 
	Ciphertext CP_yy; 
	splcp->xexe_CP   = (Ciphertext *) malloc (1 * sizeof(Ciphertext)); 
	splcp->xexe_CP 	 =  new seal::Ciphertext();
	splcp->xexe_PL   = &PL_xexe;
	splcp->GAMMA_PL  = &PL_GAMMA;
	splcp->LL_PL	 = &PL_LL;
	splcp->xGxG_PL	 = &PL_xGxG;
	splcp->xexe_CP 	 = &CP_xexe;	
	splcp->yy_CP 	 = &CP_yy;		

	// Encode the cleartext vectors into plaintext and ciphertext whichever is preferred
	makePlaintextMatRowPacking(smrp->xexe_MRP, scale, splcp->xexe_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->GAMMA_MRP, scale, splcp->GAMMA_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->LL_MRP, scale, splcp->GAMMA_PL, encoderPtr);
	makePlaintextMatRowPacking(smrp->xGxG_MRP, scale, splcp->xGxG_PL, encoderPtr);	
	encryptXVectorMatRowPacking(smrp->xexe_MRP, scale, splcp->xexe_CP, encryptorPtr, encoderPtr);	
	encryptXVectorMatRowPacking(smrp->yy_MRP, scale, splcp->yy_CP, encryptorPtr, encoderPtr);
}

/*
 	The function for encoding a vector into Plaintext 
*/
void makePlaintextMatRowPacking(std::vector<double> *x_vector, 
								double scale, 
								seal::Plaintext *plain_xePtr, 
								seal::CKKSEncoder *encoderPtr){
		encoderPtr->encode(*x_vector, scale, *plain_xePtr);
}

/*
 	The function for encrypting a vector into Ciphertext 
*/
void encryptXVectorMatRowPacking(vector<double> *x_vector,
								 double scale, 
								 Ciphertext *x_vector_EncPtr, 					   
                    			 Encryptor *encryptorPtr, 
								 CKKSEncoder *encoderPtr){
		Plaintext plain_x;
		encoderPtr->encode(*x_vector, scale, plain_x);
		encryptorPtr->encrypt(plain_x, *x_vector_EncPtr);
}

/*
 	The function for encrypting a matrix into Ciphertext 
*/
void encryptMatrixMatRowPacking(vector<double> *matPtr,
								 double scale, 
								 Ciphertext *mat_EncPtr, 					   
                    			 Encryptor *encryptorPtr, 
								 CKKSEncoder *encoderPtr){
	Plaintext plain_x;
	encoderPtr->encode(*matPtr, scale, plain_x);
	encryptorPtr->encrypt(plain_x, *mat_EncPtr);
}
