/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: The header file of the containing the printing functions (i.e., printCont.cpp)
	Note				: Commenting and cleaning up has been on July 2022
*/

#ifndef PRINTCONT_H
#define PRINTCONT_H

/* The function for printing the modulus chains and exact scales of the targeted ciphertexts for the tracking purposes */
void printScalesAndModulus( seal::Ciphertext *x3_encrypted, seal::Ciphertext *x1_encrypted, seal::Ciphertext *x2_encrypted, seal::SEALContext *contextPtr);

/* The function for printing a raw (cleartext) vector content for the tracking purposes */
void printVector(std::vector <double> *x_vector, size_t x_vector_row_size, size_t x_vector_repeat, char* vectorname);

/* The function for printing a raw (cleartext) matrix content for the tracking purposes */
void printMatrix(std::vector<double> *matPtr, size_t mat_row_size, size_t mat_col_size, char *matname);

/* The function for printing the application beginner banner for the tracking purposes */
void printAppBeginner();

/* The function for printing the equation beginner banner for the tracking purposes */
void printEqBeginner(char* eqNum);

/* The function for printing the equation end banner for the tracking purposes */
void printEqEnd(char* eqNum);		

#endif
