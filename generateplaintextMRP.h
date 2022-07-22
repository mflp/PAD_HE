/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: Header file of the class used for generating "Raw Plain (or cleartext) Vectors" in the MRP (i.e., generateplaintext.cpp)   	
   Abbreviation/Acronym	: # _MRP: Matrix-Row Packing 
	Note				: Commenting and cleaning up has been on July 2022
*/

#ifndef GENERATEPLAINTEXTMRP_H
#define GENERATEPLAINTEXTMRP_H
/* Function for generating the s vector in MRP for the homomorphic addition (i.e., s vector) */
void genVecAddOperMRP_RPL(std::vector<double> *x_add_vector, struct simulationMatrixData *smd, size_t x_vector_row_size, size_t x_vector_repeat, size_t calculated_x_row_sizeVal, char *matname);

/* The function for generating the x vector of generic Ax multiplication in MRP (i.e., x^e, y, uGamma vectors */
void genRepXVecMRP_RPL(std::vector<double> *x_vector, struct simulationMatrixData *smd, size_t x_vector_row_size, size_t x_vector_repeat, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted);

/* The function for generating the y vector (representing the sensor measurement) as the x vector of generic Ax multiplication in MRP (i.e., y vector) */
void genRepXVecMRP_RPL_v0(std::vector<double> *x_vector, double ** yy, size_t x_vector_col_size, size_t x_vector_repeat, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted);

/*
	The function for generating the matrices vector as the A matrix of generic Ax multiplication in MRP (i.e., Gamma, L, ACL, B matrices)
	Important note: This packing function is valid for n * n, n * m  matrices where m < n (e.g., m = 2, n = 10) 
*/
void genMatMRP_RPL(std::vector<double> *matPtr, struct simulationMatrixData *smd, size_t mat_row_size, size_t mat_col_size, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted);

/*
	The function for generating the matrices vector as the A matrix of generic Ax multiplication in MRP (i.e., KGamma, KL, KMinus)
	Important note: This packing function is valid for m * n matrices where m < n (e.g., m = 2, n = 10)  
*/
void genMatforMtimesNMRP_RPL(std::vector<double> *matPtr, struct simulationMatrixData *smd, size_t numCol, size_t numRow, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted);

/* Function for generating the xGamma vector in MRP for homomorphic addition (i.e., xGamma vector) */
void genXVecAddOperMRP_RPL(std::vector<double> *x_add_vector, struct simulationMatrixData *smd,  size_t x_vector_row_size, size_t x_vector_repeat, size_t calculated_x_row_sizeVal, char *matname, bool isPrinted);

/* Function for generating the vectors in MRP for homomorphic addition (i.e., KxuG, uG, Kx vectors) */
void genUVecAddOperMRP_RPL(std::vector<double> *u_add_vector, struct simulationMatrixData *smd, size_t numRow, size_t u_vector_repeat, size_t calculated_u_row_sizeVal, char *matname, bool isPrinted);

/* Function for generating the vectors in MRP for homomorphic addition (i.e., y, xp, s, v, tau, one, Chebyshev vectors) */
void genYVecAddOperMRP_RPL(std::vector<double> *y_add_vector, struct simulationMatrixData *smd, size_t numIter, size_t y_vector_row_size, size_t y_vector_repeat, size_t calculated_y_row_sizeVal, char *matname, bool isPrinted);

/* Function for generating the y vector in MRP for homomorphic addition (i.e., y vector) */
void genYVecAddOperMRP_RPL_v0 ( std::vector<double> *y_add_vector, double ** yy, size_t y_vector_row_size, size_t y_vector_repeat, size_t calculated_y_row_sizeVal, char *matname, bool isPrinted);

/* Function for generating the range transformation vector of a Chebyshev Approximation in MRP for homomorphic addition (i.e., Chebyshev Approximation vector) */
void genRangTransfVecChebApprx(std::vector<double> * Alpha_vector, std::vector<double> *Beta_vector, double lowerbound, double upperbound, size_t calculated_mask_vector_sizeVal, size_t cheb_vector_repeat);

#endif
