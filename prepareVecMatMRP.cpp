/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which is used for converting the "Raw Plain Vectors" into the matrix-row-packing format   
	Abbreviation 		: MRP - Matrix Row Packing Format  		
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

/* Import the important selected C libraries */
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

/* 
	Constructor (function) of the class 
*/
void prepareVecMatMRP(){

}

/*
	The function for calculating the each row vector size including the zero trials   
*/
void calculateVectorSize(size_t x_vector_row_size, size_t * calculated_x_row_size){

 	// Initialize the variables   
	*calculated_x_row_size 	= 0;
	int powTwo 				= 1;
	int power				= 0; 	 	
	int difference 			= x_vector_row_size - powTwo; 
	
	// Find and assign the column dimension closest to a number equaling the power of two  		
	while(difference > powTwo){ 	
		powTwo 		*= 2; 	
		power  		+= 1;
		difference 	= x_vector_row_size - powTwo; 
	}
	if(difference <=  powTwo)
		*calculated_x_row_size = powTwo * 2; 
	if(*calculated_x_row_size == 0)
		cout << "There is a row dimension calculation error. Please fix it." << endl;	  
	size_t calculated_x_row_sizeVal = *calculated_x_row_size;	
	
}

/* 
	The function for creating the simulation matrices to convert the simulation matrix data into MRP  
*/
void create_SimulationMatrixDataMRP(struct simulationMatrixData *smd, struct simulationMatrixMRP *smrp){

	/* 
	** =============================================================== 
	** ===============================================================
	** ============ MATRIX DIMENSION INITIALIZATION ================== 
	** ===============================================================
	** =============================================================== 
	*/
	// Initialize the fundamental matrix dimension parameters
	smrp->n = smd->n;	
	smrp->m	= smd->m;	
	// Calculate encrypted row size parameter  
	size_t calculated_EncRowSize = 0; 
	size_t *calculated_EncRowSizePtr = &calculated_EncRowSize;
	calculateVectorSize(smrp->n, calculated_EncRowSizePtr);
	smrp->N  = calculated_EncRowSize;

	/* 
	** =============================================================== 
	** ===============================================================
	** ============ CHEBSYHEV PARAMETER INITIALIZATION =============== 
	** ===============================================================
	** =============================================================== 
	*/		
	// Assign Chebsyhev Appx. lower and upper bounds of the observed data for 8th equation
	smrp->alpbetLowBouEq8 = smd->alpbetLowBouEq8;
	smrp->alpbetUpBouEq8  = smd->alpbetUpBouEq8; 
	// Assign Chebsyhev Appx. lower and upper bounds of the observed data for 9th equation
	smrp->alpbetLowBouEq9 = smd->alpbetLowBouEq9;
	smrp->alpbetUpBouEq9  = smd->alpbetUpBouEq9;  
	// Assign Chebyshev Appx. Alpha and Beta values for 8th equation
	smrp->alpEq8 = smd->alpEq8;  	  
	smrp->betEq8 = smd->betEq8;  
	// Assign Chebyshev Appx. Alpha and Beta values for 9th equation
	smrp->alpEq9 = smd->alpEq9; 	  
	smrp->betEq9 = smd->betEq9;
	// Assign Chebyshev Appx. Degress for 8th and 9th equation
	smrp->chebDegEq8 = smd->chebDegEq8;
	smrp->chebDegEq9 = smd->chebDegEq9;

	/* 
	** =============================================================== 
	** ===============================================================
	** ============ SIMULATION VECTOR CREATION  ====================== 
	** ===============================================================
	** =============================================================== 
	*/		
	// Create the cleartext vectors corresponding to cyberphsical system's simulation matrices and vectors in matrix row packing format  		
	// 2nd Equation 	
	smrp->xexe_MRP  = new std::vector<double>(smrp->N * smrp->n); // 2nd equation (Estimation) 
	smrp->GAMMA_MRP = new std::vector<double>(smrp->N * smrp->n); // 2nd equation (Estimation)	
	smrp->LL_MRP  	= new std::vector<double>(smrp->N * smrp->n); // 2nd equation (Estimation)
	smrp->xGxG_MRP  = new std::vector<double>(smrp->N * smrp->n); // 2nd equation (Estimation)
	smrp->yy_MRP 	= new std::vector<double>(smrp->N * smrp->n); // 2nd equation (Estimation)
	
	// 3rd Equation			
	smrp->KGKG_MRP 		= new std::vector<double>(smrp->N * smrp->n); // 3rd equation (Control Action) // Precomputed: -1 * K * Gamma, Dim = [m][n]
	smrp->KLKL_MRP 		= new std::vector<double>(smrp->N * smrp->n); // 3rd equation (Control Action) // Precomputed: -1 * K * L, Dim = [m][n]   	
	smrp->KxuGKxuG_MRP  = new std::vector<double>(smrp->N * smrp->n); // 3rd equation (Control Action) // Precomputed: -1 * (K * xG) + uG, Dim = [m][1] 
	smrp->uGuG_MRP 	= new std::vector<double>(smrp->N * smrp->n); // 3rd Equation (Prediction) // Precompute: ...,  Dim = [m][1] (for the very first iteration)
	smrp->KxKx_MRP    = new std::vector<double>(smrp->N * smrp->n); // 3rd equation (Control Action) // Precomputed: K * x^e (for the very first iteration)
	smrp->uGuG_AS_MRP = new std::vector<double>(smrp->N * smrp->n); // 3rd Equation (Prediction) // Precompute: ...,  Dim = [m][1] (for the very first iteration)

	// 4-5th Equation
	smrp->BB_MRP 	= new std::vector<double>(smrp->N * smrp->n); // 4-5th Equation (Prediction) // Plain-Known: B,  Dim = [n][m]   	
	smrp->ACL_MRP 	= new std::vector<double>(smrp->N * smrp->n); // 4-5th Equation (Prediction) // Precomputed: A-B*K,  Dim = [n][n]
	
	// 6th Equation
	smrp->yyAS_MRP  = new std::vector<double>(smrp->N * smrp->n); // 6th Equation (Residues) // Sensor Measurement (updated at each iteration): y, Dim: [n][1]  
	smrp->xpxp_MRP  = new std::vector<double>(smrp->N * smrp->n); // 6th Equation (Residues) // The very first predicition: x^p[1], Dim: [n][1]  

	// 8th Equation
	smrp->ss_MRP	= new std::vector<double>(smrp->N * smrp->n); // 8th equation (CUSUM-RELU Approximation)  
	smrp->vv_MRP    = new std::vector<double>(smrp->N * smrp->n); // 8th equation (CUSUM-RELU Approximation)

	smrp->alpEq8_MRP = new std::vector<double>(smrp->N * smrp->n); // 8th equation (CUSUM-RELU Approximation) == alpha array  		
	smrp->betEq8_MRP = new std::vector<double>(smrp->N * smrp->n); // 8th equation (CUSUM-RELU Approximation) == beta array  	
	
	smrp->One_MRP 	 		   = new std::vector<double>(smrp->N * smrp->n); // 8th equation (CUSUM-RELU Approximation)	
	smrp->chebPowSerFT_Eq8_MRP = new std::vector<double>(smrp->N * smrp->n); // 8th equation (CUSUM-RELU Approximation)	
	
 
	// 9th Equation
	smrp-> TAU_MRP  = new std::vector<double>(smrp->N * smrp->n); // 9th equation (CUSUM-Indicator Approximation)

	smrp->alpEq9_MRP = new std::vector<double>(smrp->N * smrp->n);; // 9th equation (CUSUM-Subtraction based Indicator Function Approximation) == alpha array  	
	smrp->betEq9_MRP = new std::vector<double>(smrp->N * smrp->n);; // 9th equation (CUSUM-Subtraction based Indicator Function Approximation) == beta array 

	smrp->chebPowSerFT_Eq9_MRP = new std::vector<double>(smrp->N * smrp->n);	

}

/*
	The function for generating the the simulation matrices in MRP based on the simulation data 
*/
void assignValMatrixDataMRP(struct simulationMatrixData *smd, struct simulationMatrixMRP *smrp){

	// 2nd Equation
	genRepXVecMRP_RPL(smrp->xexe_MRP, smd, smrp->n, smrp->n, smrp->N, "xe", false);
	genMatMRP_RPL(smrp->GAMMA_MRP, smd, smrp->n, smrp->n, smrp->N, "Gamma", false);
	genMatMRP_RPL(smrp->LL_MRP, smd, smrp->n, smrp->n, smrp->N, "L", false);			
	genXVecAddOperMRP_RPL(smrp->xGxG_MRP, smd, smrp->n, smrp->n, smrp->N, "xG", false); // Precomputed Vector xG
	
	// 3rd Equation
	genMatforMtimesNMRP_RPL(smrp->KGKG_MRP, smd, smrp->n, smrp->m, smrp->N, "KG", false); // Precomputed Matrix KG 
	genMatforMtimesNMRP_RPL(smrp->KLKL_MRP, smd, smrp->n, smrp->m, smrp->N, "KL", false); // Precomputed Matrix KL 
 	genUVecAddOperMRP_RPL(smrp->KxuGKxuG_MRP, smd, smrp->m, smrp->n, smrp->N, "Kxug", false); // Precomputed Vector KxuG  
	genUVecAddOperMRP_RPL(smrp->KxKx_MRP, smd, smrp->m, smrp->n, smrp->N, "Kx", false); // Precomputed Vector KxuG	
	genUVecAddOperMRP_RPL(smrp->uGuG_AS_MRP, smd, smrp->m, smrp->n, smrp->N, "uG", false); // Precomputed Vector KxuG
	
	// 4-5th Equation
	genMatMRP_RPL(smrp->ACL_MRP, smd, smrp->n, smrp->n, smrp->N, "ACL", false);
	genMatMRP_RPL(smrp->BB_MRP, smd, smrp->n, smrp->m, smrp->N, "B", false);
	genRepXVecMRP_RPL(smrp->uGuG_MRP, smd, smrp->m, smrp->n, smrp->N, "uG", false);
	 	
	// 6th equation
	// genYVecAddOperMRP_RPL(smrp->yyAS_MRP, smd, 0, smrp->n, smrp->n, smrp->N, "y", false);
	genYVecAddOperMRP_RPL(smrp->xpxp_MRP, smd, 0, smrp->n, smrp->n, smrp->N, "xp", false);

	// 8th equation
	genYVecAddOperMRP_RPL(smrp->ss_MRP, smd, 0, smrp->n, smrp->n, smrp->N, "s", false);
	genYVecAddOperMRP_RPL(smrp->vv_MRP, smd, 0, smrp->n, smrp->n, smrp->N, "v", false);
	// Chebyshev Appx.  
	genYVecAddOperMRP_RPL(smrp->One_MRP, smd, 0, smrp->n, smrp->n, smrp->N, "one", false); // One Vector 

	// Chebyshev Appx. - 8th Equation		
	genRangTransfVecChebApprx(smrp->alpEq8_MRP, smrp->betEq8_MRP, smrp->alpbetLowBouEq8, smrp->alpbetUpBouEq8, smrp->N, smrp->n); // alpha, beta
	genYVecAddOperMRP_RPL(smrp->chebPowSerFT_Eq8_MRP, smd, 0, smrp->n, smrp->n, smrp->N, "chb_1st_T_D_12_m_2_n_10", false);

	// 9th Equation
	genYVecAddOperMRP_RPL(smrp->TAU_MRP, smd, 0, smrp->n, smrp->n, smrp->N, "tau", false);
	// Chebyshev Appx. - 9th Equation
	genRangTransfVecChebApprx(smrp->alpEq9_MRP, smrp->betEq9_MRP, smrp->alpbetLowBouEq9, smrp->alpbetUpBouEq9, smrp->N, smrp->n); // alpha, beta
	genYVecAddOperMRP_RPL(smrp->chebPowSerFT_Eq9_MRP, smd, 0, smrp->n, smrp->n, smrp->N, "chb_2nd_T_D_12_m_2_n_10", false);
}
