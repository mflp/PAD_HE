/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which contains the printing functions for tracking the Matrix-Vector Multiplications and Vector Additions
	Note				: Commenting and cleaning up has been on July 2022
*/

/* Import the other classes' header files*/
#include "examples.h" 			  // The class for the console menu and guiding the user to the preferred application   

/* Import the important selected C libraries*/
#include <float.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>

/* Call main namespaces */
using namespace std;
using namespace seal;

/* 
	The function for printing the modulus chains and exact scales of the targeted ciphertexts for the tracking purposes  
*/
void printScalesAndModulus(
	seal::Ciphertext *x3_encrypted, 
	seal::Ciphertext *x1_encrypted,
	seal::Ciphertext *x2_encrypted,
	seal::SEALContext *contextPtr){

	// Print modulus chain indexes  
	// Important note: Encrypted addition, multiplication, and subtraction requires that the modulus chains of the inputs are the same 
    cout << endl;
    print_line(__LINE__);
    cout << "Parameters used by all three terms are different." << endl;
    cout << "    + Modulus chain index for 			vector-1: " << (*contextPtr->get_context_data(x1_encrypted->parms_id())).chain_index() << endl;
    cout << "    + Modulus chain index for 			vector-2: " << (*contextPtr->get_context_data(x2_encrypted->parms_id())).chain_index() << endl;

	// Print scales of the ciphertexts  
	// Important note: Encrypted addition, multiplication, and subtraction requires that the exact scales of the inputs are the same
    print_line(__LINE__);
    cout << "The exact scales of all three terms are different:" << endl;
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);
    cout << fixed << setprecision(10);
    cout << "    + Exact scale in  		vector-1: " << x1_encrypted->scale() << endl;
    cout << "    + Exact scale in  		vector-2: " << x2_encrypted->scale() << endl;
    cout << endl;
    cout.copyfmt(old_fmt);
}

/* 
	The function for printing a raw (cleartext) vector content for the tracking purposes   
*/
void printVector(std::vector <double> *x_vector, size_t x_vector_row_size, size_t x_vector_repeat, char* vectorname){

	printf("Input %s vector: \n", vectorname);
	printf("Row size of %s vector: %d \n", vectorname, x_vector_row_size); 
	for(size_t j = 0; j < x_vector_row_size * x_vector_repeat; j++){
		printf("%f ", x_vector->at(j));
		if(j % x_vector_row_size == (x_vector_row_size - 1) )
			printf("\n");
	}
	printf("\n");	
}

/* 
	The function for printing a raw (cleartext) matrix content for the tracking purposes  
*/
void printMatrix(std::vector<double> *matPtr, size_t mat_row_size, size_t mat_col_size, char *matname){

	printf("Input %s matrix: \n", matname);
	for(size_t j = 0; j < mat_row_size * mat_col_size; j++){
		printf("%f ", matPtr->at(j));
		if(j % mat_col_size == (mat_col_size - 1))
			printf("\n");
	}
	printf("\n");	
}

/* 
	The function for printing the application beginner banner for the tracking purposes  
*/
void printAppBeginner(){
  printf("##########################-BEGINNING-OF-OUR-APPLICATION-##########################\n");
}

/* 
	The function for printing the equation beginner banner for the  tracking purposes   
*/
void printEqBeginner(char* eqNum){
			
	printf("\n\n###############################################################################\n");
	printf("##########################-BEGINNING-OF-%s.-EQUATION-##########################\n", eqNum);
	printf("###############################################################################\n");		
	printf("###############################################################################\n");	
}

/* 
	The function for printing the equation end banner for the  tracking purposes   
*/
void printEqEnd(char* eqNum){		
	printf("\n\n###############################################################################\n");
	printf("#############################-END-OF-%s.-EQUATION-##############################\n", eqNum);
	printf("###############################################################################\n");		
	printf("###############################################################################\n");	
}
