/*
   Author				: Mestan Fırat Çeliktuğ 
   Date	 				: 27.02.2022 - 30.03.2022
   Description			: C++ class which is used to handle the cleartext (or raw) 
						  simulation data with the initializations
	Note				: Commenting and cleaning up has been on July 2022
*/

/* Import the other classes' header files*/
#include "examples.h"  // The class for the console menu and guiding the user to the preferred application   
#include "rawplain.h"  // The class which reads and stores the plain matrices
#include "printCont.h" // The class containing the printing functions for control purposes 

/* Import the important selected C libraries*/
#include <iostream>
#include <dirent.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

/* Call main namespaces */
using namespace std;
using namespace seal;

// Create Chebyshev Polynomial Coefficients Array for computing the RELU and Comparison functions for the approximation degree of twelve (12) 
double psCfAr_D12_Max []  = {15.06831662919437, 17.4154228941,  0.9390991874,  -0.729322333,  0.4948930898, -0.2756084825, 0.1017867105, 0.0117719738, -0.0666850779, 0.0769685688, -0.0619403765, 0.0388042051, -0.0175791198};
double psCfAr_D12_ISub [] = {0.15384615384615385, 0.2965730949, 0.2645311639, 0.2153445253, 0.1547681086, 0.0897951079, 0.0277609233, -0.0245940319, -0.0619810878, -0.0812285695, -0.0816694059, -0.065201076, -0.0360104934};

// Create Chebyshev Polynomial Coefficients Array for computing the RELU or Comparison functions for the approximation degree of sixteen (16) 
double psCfAr_D16_Max []  = {10.874502336837981, 13.3664957117,  1.3189604909, -0.8899067475,  0.4540466982, -0.1065472449, -0.0982366249,  0.1577428186, -0.1126917097,  0.0247397879,  0.0493727727, -0.0769482213,  0.0556748187, -0.0070457589, -0.0383123412, 0.0563087367, -0.0397478465}; 
double psCfAr_D16_ISub [] = {0.23529411764705882,  0.4294988633,  0.3187632779,  0.1715328398,  0.0299212305, -0.0694725895, -0.1074743965, -0.0873774274, -0.0315417058,  0.0287540641,  0.0659843703, 0.066545329, 0.0345932392, -0.0115915465, -0.0488053008, -0.0595871836, -0.0397989307};

/* Define constant */
#define maxlineLength 100000
char line[maxlineLength];

/*
	The function for assigning the column and row dimensions for the simulation (i.e., m and n respectively) 
*/
void assDimToSmd(struct simulationMatrixData *smd, char * DimFolderDir, int numIter){

	// Define the dimension file path in the working folder  
	char dimFileD[1000000]; 
	strcat(dimFileD, DimFolderDir);
	strcat(dimFileD, "/Dim.txt");	
	printf("Dimension file directory '%s' \n", dimFileD);			

	// Read file dimension parameters	
	FILE * file = fopen(dimFileD, "r");
	if (file == NULL) {
		perror("Failed to open general dimensions file\n");
		exit(0);
    }		
	fgets(line, 1000, file);
	fgets(line, 1000, file);
	printf("%s\n", line);	

	int *param = (int *) calloc(2, sizeof(int));
	char *token;
 	char* rest = line; 
	int i = 0;
    while ((token = strtok_r(rest, " ", &rest))){
        param[i] = atoi(token);
		i++;	
	}

	// Close the file	   	
	fclose(file);

	// Assign the zeros back to the file pointers by the memory setting 
	memset(dimFileD, 0, 100000); // Note: Without this line, yFD does not work correct although it is locally defined.

	// Assign file dimension parameters
	smd->m = param[0];
	smd->n = param[1];
	smd->tMax = numIter; 
	printf("m: %d\n", param[0]);
	printf("n: %d\n", param[1]);
}

/*
	The function for reading and filling the recorded y vector data 
*/
void readYData(struct simulationMatrixData *smd, char * directory){

	// Print the y vector records directory 	
	printf("Y-dir: '%s' \n", directory);	
	
	// Read the y-vector file    
	for(int k = 0; k < 15; k++){
		char lp[5];	
		sprintf(lp, "%d", k + 1);
		
		// Define the y file directory variable	
		char yFD[100000];
		memset(yFD, 0, 100000);
		strcat(yFD, directory);
		strcat(yFD, lp);
		strcat(yFD, ".txt");
		printf("YFD: %s \n", yFD);
		
		// Define the file pointer 
		FILE * fp;
		char * line = NULL;
		size_t len = 0;
		ssize_t read;
		fp = fopen(yFD, "r");
		if (fp == NULL)
		    exit(EXIT_FAILURE);
		
		// Read the y file and assign the y vector measurement 	
		int i = 0;
		while ((read = getline(&line, &len, fp)) != -1) {
			char *ptr;
			double sensorMeasurement;
			sensorMeasurement = strtod(line, &ptr);
			int quotient  = i/smd->n; // i/10;
			int remainder = i%smd->n; // i%10;	
			smd->yy[k * 120 + quotient][remainder] = sensorMeasurement;   
			i++;
		}
		fclose(fp);
		if(line)
		    free(line);

		// Assign the zeros back to the file pointers by the memory setting 		
		memset(yFD, 0, 100000); // Note: Without this line, yFD does not work correct although it is locally defined.
		memset(lp, 0, 5); 
	} 
}

/*
	The function for checking the content of the y vector 
*/
void readYMCheck(struct simulationMatrixData *smd, char * matname, int dim1, int dim2){
	// Print the matrix name 
	printf("Matrix '%s' read from the file shared by Luis- Sanity Check \n", matname);
	// Print the content of the y vector 	
	for (int i = 0; i < dim1; i++){
		printf("Line %d :", i + 1);
		for (int j = 0; j < dim2; j++){
			if(strcmp(matname, "y") == 0)	
				printf("%f ", smd->yy[i][j]);
		}	
		printf("\n");
	}
}

/*
	The function for checking the content of the almost each vector and matrix except the  
*/
void readMCheck(struct simulationMatrixData *smd, char * matname, int dim1, int dim2){

	// Print the beginner banner of the matrix content  	
	printf("Matrix '%s' read from the file shared by Luis- Sanity Check \n", matname);	
	cout << "Dimension-1: " << dim1 << endl;
	cout << "Dimension-2: " << dim2 << endl;

	// Print the content of the selected vector or matrix   	
	for (int i = 0; i < dim1; i++){
		for (int j = 0; j < dim2; j++){
			// Target matrix list to be read from file sent by Luis
			if(strcmp(matname, "A") == 0)
				printf("%f ", smd->AA[i][j]);
			if(strcmp(matname, "C") == 0)
				printf("%f ", smd->CC[i][j]);
			if(strcmp(matname, "L") == 0)
				printf("%f ", smd->LL[i][j]);
			if(strcmp(matname, "B") == 0)
				printf("%f ", smd->BB[i][j]);				
			if(strcmp(matname, "K") == 0)
				printf("%f ", smd->KK[i][j]);				

			// Target vector list to be read from file sent by Luis 				
			if(strcmp(matname, "v") == 0)
				printf("%f ", smd->vv[i][j]);
			if(strcmp(matname, "ur") == 0)
				printf("%f ", smd->urur[i][j]);	 
			if(strcmp(matname, "xr") == 0)
				printf("%f ", smd->xrxr[i][j]);
			if(strcmp(matname, "Tau") == 0)
				printf("%f ", smd->TAU[i][j]);  
			if(strcmp(matname, "x0") == 0)
				printf("%f ", smd->xx[i][j]);
			if(strcmp(matname, "GAMMA") == 0)
				printf("%f ", smd->GAMMA[i][j]);
			if(strcmp(matname, "ACL") == 0)
				printf("%f ", smd->ACL[i][j]);	
			if(strcmp(matname, "KG") == 0)
				printf("%f ", smd->KGKG[i][j]);	
			if(strcmp(matname, "KL") == 0)
				printf("%f ", smd->KLKL[i][j] );	
			if(strcmp(matname, "Kxug") == 0)
				printf("%f ", smd->KxuGKxuG[i][j]);
			if(strcmp(matname, "uG") == 0)
				printf("%f ", smd->uGuG[i][j]); 
			if(strcmp(matname, "xG") == 0)
				printf("%f ", smd->xGxG[i][j]);	
			if(strcmp(matname, "yNoise") == 0)
				printf("%f ", smd->yNoise[i][j]);
			if(strcmp(matname, "xNoise") == 0)
				printf("%f ", smd->xNoise[i][j]);	
		}	
		printf("\n");
	}
}

/*
	The function for reading and filling the simulation vectors and matrices  
*/
void readMatrix(struct simulationMatrixData *smd, char * dataPath, char * matname, int dim1, int dim2){  
	
	// Define the file pointer 
	FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
	
	// Open the selected file 
    fp = fopen(dataPath, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

	// Read iteratively content of the selected vector or matrix 
	for(int i = 0; i < dim1 && (read = getline(&line, &len, fp)) != -1; i++){
		char *pt;
    	pt = strtok (line,","); // Tokenize the read line
		int j = 0; 
 		while (pt != NULL) {
			double matV = atof(pt);

			// Target matrix list 	
			if(strcmp(matname, "A") == 0)
				smd->AA[i][j] = matV; 		 
			if(strcmp(matname, "C") == 0)
				smd->CC[i][j] = matV;
			if(strcmp(matname, "L") == 0)
				smd->LL[i][j] = matV;
			if(strcmp(matname, "B") == 0)
				smd->BB[i][j] = matV;
			if(strcmp(matname, "K") == 0)
				smd->KK[i][j] = matV;

			// Target vector list 				
			if(strcmp(matname, "v") == 0)
				smd->vv[i][j] = matV;
			if(strcmp(matname, "ur") == 0)
				smd->urur[i][j] = matV; 
			if(strcmp(matname, "xr") == 0)
				smd->xrxr[i][j] = matV;
			if(strcmp(matname, "Tau") == 0)
				smd->TAU[i][j]  = matV; 
			if(strcmp(matname, "x0") == 0)
				smd->xx[i][j]   = matV;	

			// Target precomputed vector and matrix list 
			if(strcmp(matname, "GAMMA") == 0)
				smd->GAMMA[i][j] = matV;	
			if(strcmp(matname, "ACL") == 0)
				smd->ACL[i][j]   = matV;	
			if(strcmp(matname, "KG") == 0)
				smd->KGKG[i][j]   = matV;	
			if(strcmp(matname, "KL") == 0)
				smd->KLKL[i][j]   	= matV;	
			if(strcmp(matname, "Kxug") == 0)
				smd->KxuGKxuG[i][j] = matV;
			if(strcmp(matname, "Kx") == 0)
				smd->KxKx[i][j]  = matV;
			if(strcmp(matname, "uG") == 0)
				smd->uGuG[i][j] = matV; 
			if(strcmp(matname, "xG") == 0)
				smd->xGxG[i][j] = matV;	
			if(strcmp(matname, "yy") == 0)
				smd->yy[i][j] = matV;	

			// Target noise vector list 	
			if(strcmp(matname, "xNoise") == 0)	
				smd->xNoise[i][j] = matV;
			if(strcmp(matname, "yNoise") == 0)	
				smd->yNoise[i][j] = matV;
        	pt = strtok (NULL, ",");
			j++;	
    	}	
	}

	// Close the selected file
	fclose(fp);
}

/*
	The function for reading a file  (Probably not-used) 
*/
void readFile(double **data, std::string input, int *size) { 
	int i = 0; 
	int j = 0; 
	fstream file; 
	string inputLine, stringinput; 
	file.open(input, ios::in); 
	while (getline(file, inputLine)) { 
		j = 0; stringstream line(inputLine); 
		while (getline(line, stringinput, ',')){ 
			data[i][j] = stod(stringinput);  
			j++; 
		} 
		i++; 
	} 
} 

/*
	The function for printing the experimental result (both the sensor and CUSUM statistics) 
*/
void printExperimentalResult(struct simulationMatrixData *smd){
		
	// Extract the alarm computation results  
	cout << "s vector result " << endl; 		
	for (int i = 0; i < smd ->tMax + 1; i++){
		for (int j = 0; j < smd ->n; j++)
			cout << smd->s_Res[i][j] << " ";  
		cout << endl;		
	}

	// Extract the sensor measurement results
	cout << "y vector result " << endl; 		
	for (int i = 0; i < smd ->tMax + 1; i++){
		for (int j = 0; j < smd ->n; j++)
			cout << smd->y_Res[i][j] << " ";  
		cout << endl;		
	}
	
}

/*
	The function for reading and assigning the values of the initial raw (cleartext) vectors and matrices 
*/
void assignValMatrixDatabyFileRead(struct simulationMatrixData *smd, char * folderPath){
	
	// Concatenate the folderPath with the specific matrix file path 
	// Create the vector and matrix file paths  	 	
	char A_p[100000],  C_p[100000],   K_p[100000],  v_p[100000],   ur_p[100000],  xr_p[100000], B_p[100000],  L_p[100000],  Tau_p[100000], x0_p[100000];
	char GAMMA_p[100000], KG_p[100000], KL_p[100000], KxuG_p[100000], Kx_p[100000], ACL_p[100000],  xG_p[100000],  uG_p[100000]; // Precomputed Matrices
	char y_p[100000], xNoise_p[100000], yNoise_p[100000];
	// Known matrices' and vectors' file paths
	strcat(A_p, folderPath);
	strcat(A_p, "/A.txt");	
	strcat(B_p, folderPath);
	strcat(B_p, "/B.txt");	
	strcat(C_p, folderPath);	
	strcat(C_p, "/C.txt");	
	strcat(K_p,  folderPath);
	strcat(K_p,  "/K.txt"); 	
	strcat(v_p,  folderPath);
	strcat(v_p,  "/Nu.txt"); 
	strcat(ur_p,  folderPath);
	strcat(ur_p,  "/ur.txt"); 
	strcat(xr_p, folderPath);
	strcat(xr_p,  "/xr.txt"); 
	strcat(L_p,  folderPath);
	strcat(L_p,  "/L.txt"); 
	strcat(Tau_p, folderPath); 
	strcat(Tau_p,  "/Tau.txt"); 
	strcat(x0_p, folderPath);
	strcat(x0_p,  "/x0.txt"); 
	strcat(y_p, folderPath);
	strcat(y_p,  "/yy.txt"); 
	// Precomputed matrices' and vectors' file paths	
	strcat(GAMMA_p, folderPath);
	strcat(GAMMA_p, "/PreCompute/Gamma.txt");
	strcat(KG_p, folderPath);
	strcat(KG_p, "/PreCompute/KG.txt");
	strcat(KL_p, folderPath);
	strcat(KL_p, "/PreCompute/KL.txt");
	strcat(KxuG_p, folderPath);
	strcat(KxuG_p, "/PreCompute/kxuG.txt");
	strcat(ACL_p, folderPath);
	strcat(ACL_p, "/PreCompute/Acl.txt");
	strcat(xG_p, folderPath);
	strcat(xG_p, "/PreCompute/xG.txt");
	strcat(uG_p, folderPath);
	strcat(uG_p, "/PreCompute/uG.txt");
	strcat(Kx_p, folderPath);
	strcat(Kx_p, "/PreCompute/Kx.txt");
	// Noise vectors' file paths	
	strcat(xNoise_p, folderPath);
	strcat(xNoise_p,  "/process_noise.txt"); 
	strcat(yNoise_p, folderPath);
	strcat(yNoise_p,  "/sensor_noise.txt"); 

	// Read the given matrices into smd struct's matrices	
	readMatrix(smd, A_p, "A", smd->n, smd->n);	
	readMatrix(smd, C_p, "C", smd->n, smd->n);
	readMatrix(smd, L_p, "L", smd->n, smd->n);
	readMatrix(smd, B_p, "B", smd->n, smd->m);
	readMatrix(smd, K_p, "K", smd->m, smd->n);
	// Read the given vectors into smd struct's vectors
	readMatrix(smd, v_p, "v" , smd->n, 1);
	readMatrix(smd, ur_p, 	"ur", smd->m, 1);
	readMatrix(smd, xr_p, 	"xr", smd->n, 1);
	readMatrix(smd, Tau_p, "Tau", smd->n, 1);
	readMatrix(smd, x0_p, 	"x0", smd->n, 1);
	// Read the precomputed vectors and matrices into smd struct's vectors
	readMatrix(smd, GAMMA_p, "GAMMA", smd->n, smd->n);	
	readMatrix(smd, ACL_p, "ACL", smd->n, smd->n);	
	readMatrix(smd, KG_p, "KG", smd->m, smd->n);	
	readMatrix(smd, KL_p, "KL", smd->m, smd->n);	
	readMatrix(smd, KxuG_p, "Kxug", smd->m, 1);	
	readMatrix(smd, Kx_p, "Kx", smd->m, 1);
	readMatrix(smd, uG_p, "uG", smd->m, 1);	
	readMatrix(smd, xG_p, "xG", smd->n, 1);
	// Read the noise vectors into smd struct's vectors
	readMatrix(smd, xNoise_p, "xNoise", 201, smd->n);
	readMatrix(smd, yNoise_p, "yNoise", 201, smd->n);

	// Assign the zeros back to the file pointers by the memory setting 
	memset(A_p, 0, 100000);
	memset(C_p, 0, 100000);
	memset(K_p, 0, 100000);
	memset(v_p, 0, 100000);
	memset(ur_p, 0, 100000);
	memset(xr_p, 0, 100000);
	memset(L_p, 0, 100000);
	memset(B_p, 0, 100000);
	memset(Tau_p, 0, 100000);
	memset(x0_p, 0, 100000);
	memset(GAMMA_p, 0, 100000);
	memset(KG_p, 0, 100000);
	memset(KL_p, 0, 100000);
	memset(Kx_p, 0, 100000);
	memset(KxuG_p, 0, 100000);
	memset(ACL_p, 0, 100000);
	memset(xG_p, 0, 100000);
	memset(uG_p, 0, 100000);
	memset(xNoise_p, 0, 100000);
	memset(yNoise_p, 0, 100000);
	memset(y_p, 0, 100000);	
}

/*
	The function for creating the simulation matrix data with the empty slots   
*/
void create_SimulationMatrixData(struct simulationMatrixData *smd){

	// Define folder paths 
	char * folderPath_y10_u2 	= "./all_data/y10_u2";
	char * folderPath_y20_u4 	= "./all_data/y20_u4";
	char * folderPath_y50_u10 	= "./all_data/y50_u10";

	// Assign m, n dimensions  	
	// assDimToSmd(smd, folderPath_y10_u2, 120); // Change here based on the used matrices
	// assDimToSmd(smd, folderPath_y20_u4, 120); // Change here based on the used matrices	
	assDimToSmd(smd, folderPath_y50_u10, 120); // Change here based on the used matrices
		
	// Assign Chebyshev Approximation Degrees 
	smd->chebDegEq8 = 16;
	smd->chebDegEq9 = 16; 	
	
	// Print the dimensions 
	printf("m: %d \n", smd->m);		
	printf("n: %d \n", smd->n);
	printf("tMax: %d \n", smd->tMax);		

	// Create each simulation matrix (that is the place where both mathematically simulation matrices and vectors are stored) first pointers
	smd->AA  		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][10]	
	smd->BB  		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: B[10][2]
	smd->CC  		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: C[10][10]
	smd->VV  		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: V[10][10]
	smd->KK  		= (double **) calloc (smd->m, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: K[2][10] 
	smd->LL  		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: L[10][10]
	smd->xrxr 		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: xr[10][1]
	smd->urur 		= (double **) calloc (smd->m, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: ur[2][1]
	smd->GAMMA  	= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][10]
	smd->xGxG 		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][1]
	smd->uGuG 		= (double **) calloc (smd->m, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [2][1]
	smd->vv 		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][1]
	smd->TAU 		= (double **) calloc (smd->n, sizeof(double *)); //	Exemplary dimensions when m = 2, n = 10: [10][1]
	smd->ALARMSYS 	= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][1]
	smd->ss 		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][1]
	smd->xx 		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][1]
	smd->xexe      	= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][1]
	smd->xpxp 		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][1] 
	smd->yy 		= (double **) calloc (2000, sizeof(double *));   // Exemplary dimensions when m = 2, n = 10: [10][1]
	
	smd->xNoise 	= (double **) calloc (201, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][1]
	smd->yNoise 	= (double **) calloc (201, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [10][1]

	smd->KGKG 		= (double **) calloc (smd->m, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: KG[2][10]  
	smd->KLKL 		= (double **) calloc (smd->m, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: KL[2][10]  
	smd->KxuGKxuG	= (double **) calloc (smd->m, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: KxuG[2][1]
	smd->KxKx 		= (double **) calloc (smd->m, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: [2][1]	
	smd->ACL 	    = (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: ACL[10][10]
		
	// Create Chebyshev Approximation Arrays
	smd->One 		= (double **) calloc (smd->n, sizeof(double *)); // Exemplary dimensions when m = 2, n = 10: Vector One[10][1]
	smd->eq8maxAppx_PS_Coeff_D12_y_10_u_2 = (double **) calloc (smd->chebDegEq8 + 1, sizeof(double *)); // First term  
	smd->eq8maxAppx_PS_FT_D12_y_10_u_2 	  = (double **) calloc (smd->n, sizeof(double *)); // Chebyshev Coefficients	
		 
	smd->eq9ISubAppx_PS_Coeff_D12_y_10_u_2 = (double **) calloc (smd->chebDegEq9 + 1, sizeof(double *)); // First term
	smd->eq9ISubAppx_PS_FT_D12_y_10_u_2    = (double **) calloc (smd->n, sizeof(double *)); // Chebyshev Coefficients

	// Experimental Result Vectors
	int numIter   = 120;
	int numStates = smd->n;   
	smd->xe_Res 	 = (double **) calloc(smd->tMax + 1, sizeof (double*)); 
	smd->u_Res  	 = (double **) calloc(smd->tMax + 1, sizeof (double*));
	smd->xp_Res      = (double **) calloc(smd->tMax + 1, sizeof (double*)); 
	smd->residue_Res = (double **) calloc(smd->tMax + 1, sizeof (double*));
	smd->sBar_Res 	 = (double **) calloc(smd->tMax + 1, sizeof (double*));
	smd->indInp_Res  = (double **) calloc(smd->tMax + 1, sizeof (double*)); 
	smd->alarm_Res 	 = (double **) calloc(smd->tMax + 1, sizeof (double*));
	smd->s_Res 		 = (double **) calloc(smd->tMax + 1, sizeof (double*));
	smd->x_Res		 = (double **) calloc(smd->tMax + 1, sizeof (double*));
	smd->y_Res 		 = (double **) calloc(smd->tMax + 1, sizeof (double*));
	
	// Create each simulation result vector second pointers
	for(int i = 0; i < smd->tMax + 1; i++){
		smd->xe_Res[i]   	= (double *) calloc (smd->n, sizeof(double)); 		 
		smd->u_Res[i]  	 	= (double *) calloc (smd->m, sizeof(double));		
		smd->xp_Res[i] 	    = (double *) calloc (smd->n, sizeof(double));		 
		smd->residue_Res[i] = (double *) calloc (smd->n, sizeof(double)); 
		smd->sBar_Res[i] 	= (double *) calloc (smd->n, sizeof(double));
		smd->indInp_Res[i]  = (double *) calloc (smd->n, sizeof(double));
		smd->alarm_Res[i] 	= (double *) calloc (smd->n, sizeof(double));
		smd->s_Res[i] 		= (double *) calloc (smd->n, sizeof(double));  
		smd->x_Res[i]		= (double *) calloc (smd->n, sizeof(double));
		smd->y_Res[i] 		= (double *) calloc (smd->n, sizeof(double)); 
	}

	// Create each simulation matrix (that is the place where both mathematically simulation matrices and vectors are stored) second pointers
	for(int i = 0; i < smd->n; i++){
		if(i < smd->m){
			smd->KK[i]   	 = (double *) calloc (smd->n, sizeof(double));
			smd->urur[i] 	 = (double *) calloc (1, sizeof(double));
			smd->uGuG[i] 	 = (double *) calloc (1, sizeof(double));

			smd->KGKG[i] 	 = (double *) calloc (smd->n, sizeof(double)); 
			smd->KLKL[i] 	 = (double *) calloc (smd->n, sizeof(double)); 
			smd->KxuGKxuG[i] = (double *) calloc (1, sizeof(double));      
			smd->KxKx[i] 	 = (double *) calloc (1, sizeof(double));      
		}
		smd->AA[i]  	 = (double *) calloc (smd->n, sizeof(double));					
		smd->BB[i]  	 = (double *) calloc (smd->m, sizeof(double)); 
		smd->CC[i]  	 = (double *) calloc (smd->n, sizeof(double)); 
		smd->LL[i]  	 = (double *) calloc (smd->n, sizeof(double)); 
		smd->xrxr[i] 	 = (double *) calloc (1, sizeof(double));      
		smd->GAMMA[i] 	 = (double *) calloc (smd->n, sizeof(double)); 
		smd->xGxG[i] 	 = (double *) calloc (1, sizeof(double));      
		smd->vv[i] 		 = (double *) calloc (1, sizeof(double));      
		smd->TAU[i] 	 = (double *) calloc (1, sizeof(double));      
		smd->ALARMSYS[i] = (double *) calloc (1, sizeof(double));      
		smd->ss[i] 		 = (double *) calloc (1, sizeof(double));      
		smd->xx[i] 		 = (double *) calloc (1, sizeof(double));      
		smd->xexe[i]     = (double *) calloc (1, sizeof(double));      
		smd->xpxp[i] 	 = (double *) calloc (1, sizeof(double));      	
		
		smd->ACL[i]      = (double *) calloc (smd->n, sizeof(double)); 
		
		smd->One[i] 	 					   = (double *) calloc (1, sizeof(double));  
		smd->eq8maxAppx_PS_FT_D12_y_10_u_2[i]  = (double *) calloc (1, sizeof(double)); // First term array of 8th equation 
		smd->eq9ISubAppx_PS_FT_D12_y_10_u_2[i] = (double *) calloc (1, sizeof(double)); // First term array of 9th equation
	}	

	// Create y sensor measurement matrix
	for(int i = 0; i < 2000; i++)  // Could be changed to 1800 afterwards
		smd->yy[i]  = (double *) calloc (smd->n, sizeof(double));

	// Create power series matrices for Chebyshev Approximation
	// Create power series matrices for Chebyshev Approximation of 8th equation
	for(int i = 0; i < smd->chebDegEq8 + 1; i++)
		smd->eq8maxAppx_PS_Coeff_D12_y_10_u_2[i]  = (double *) calloc (1, sizeof(double));  // Chebyshev Power Series coefficient array of 8th equation
	// Create power series matrices for Chebyshev Approximation of 9th equation
	for(int i = 0; i < smd->chebDegEq9 + 1; i++)
		smd->eq9ISubAppx_PS_Coeff_D12_y_10_u_2[i] = (double *) calloc (1, sizeof(double)); // Chebyshev Power Series coefficient array of 9th equation	

	// Create x and y noise vectors	
	for(int i = 0; i < 201; i++){ 
		smd->xNoise[i]  = (double *) calloc (smd->n, sizeof(double));	
		smd->yNoise[i]  = (double *) calloc (smd->n, sizeof(double));
	}
}

/*
	The function for initializing the remaining (not-read from the records received) vectors including the sensor measurement vector y 
	Status: Possibly not-used 
*/
void initRemainVec(struct simulationMatrixData *smd){

	// Prtint the beginning banner 
	cout << "IN \n" << endl;	

	// Assign the values for the remaining vectors		
	for(int i = 0; i < smd->n; i++){ 	
		smd->ALARMSYS[i][0] = 0; 			 // Alarm vector 
		smd->ss[i][0] 		= 0;			 // s vector
		smd->xexe[i][0] 	= smd->xx[i][0]; // x^p vector
		smd->xpxp[i][0] 	= smd->xx[i][0]; // x^e vector
		smd->One[i][0] 		= 1; 			 // Chebyshev Appx. One vector	
	}

	// Read measurement data - y vector 	
	char y_MeasSensDir_y10_u2[1000000]; 
	strcat(y_MeasSensDir_y10_u2, "/home/mfcustben/SEAL_1/native/examples/ZZ_Last_points/Approximation/0_PhD-SimulationData-m=2-n=10/0_Corresponding_Y_Matrix_Measurements/YMatrixContent_for_y10_u2_"); 
	char y_MeasSensDir_y50_u10[1000000];
	strcat(y_MeasSensDir_y50_u10, "/home/mfcustben/SEAL_1/native/examples/ZZ_Last_points/Approximation/0_PhD-SimulationData-m=2-n=10/0_Corresponding_Y_Matrix_Measurements/YMatrixContent_for_y50_u10_"); 
	readYData(smd, y_MeasSensDir_y10_u2);
	memset(y_MeasSensDir_y10_u2, 0, 100000);

} 

/*
	The function for assigning the Chebyshev Approximation constants and arrays for the 8th Equation (Max Function) and 9th Equation (Subtraction-based Indicator Function) 
*/
void assignCUSUMChebyshevAppxParams(struct simulationMatrixData *smd){

	/*
		Func	a 	b	Alpha			Beta +1
		Max		-4	33	0.0540540541	0.7837837838
		I-Sub	-31	2	0.0606060606	-0.8787878788
	*/

	// Print the begin banner
	// cout << "---------ININIINININ---------\n"<< endl;

	// Max(or RELU) Approximation constant variables 	
	smd->alpbetLowBouEq8 	= -5;  // Min value (namely a)
	smd->alpbetUpBouEq8  	= 25;  // Max value (namely b)
	smd->alpEq8 			= 0.0666666667;	// alpha  
	smd->betEq8 			= 0.6666666667; // beta + 1	
	// smd->alpbetLowBouEq8 = -4;  // Min value (namely a)
	// smd->alpbetUpBouEq8  = 33;  // Max value (namely b)
	// smd->alpEq8 = 0.0540540541;	// alpha  
	// smd->betEq8 = 0.7837837838; // beta + 1

	// Subtraction-based Indicator Equation constant variables
	smd->alpbetLowBouEq9 	= -31;  // Min value (namely a)
	smd->alpbetUpBouEq9  	= 5; 	 // Max value (namely b) 
	smd->alpEq9 			= 0.0555555556;	 // alpha    
	smd->betEq9 			= -0.7222222222; // beta + 1

	// Assign the Chebyshev Coefficient array of max-8th Equation  
	for(int i = 0; i < smd->chebDegEq8 + 1; i++)
	 smd->eq8maxAppx_PS_Coeff_D12_y_10_u_2[i][0]  = psCfAr_D16_Max[i];
	
	// Assign the Chebyshev Coefficient array of max-9th Equation
	for(int i = 0; i < smd->chebDegEq9 + 1; i++)
	 smd->eq9ISubAppx_PS_Coeff_D12_y_10_u_2[i][0] = psCfAr_D16_ISub[i];
	
	// Assign the content of the First Term array
	for(int i = 0; i < smd->n; i++){	
	 	smd->eq8maxAppx_PS_FT_D12_y_10_u_2[i][0]  = psCfAr_D16_Max[0];
  		smd->eq9ISubAppx_PS_FT_D12_y_10_u_2[i][0] = psCfAr_D16_ISub[0]; 
	}
}

// SENSOR MEASUREMENT UPDATES
/*
	The function for assigning the resulting x vector 
*/
void assignXResult(struct simulationMatrixData *smd, double ** result, int rowDm, int colDm){
	for(int i = 0; i < rowDm; i++)
		for(int j = 0; j < colDm; j++)
			smd->xx[i][j] = result[i][j];
}

/*
	The function for poerforming the matrix multiplication 
*/
void matrixMult(double ** result, double ** vect1, double ** vect2, int rowDm1, int colDm1, int rowDm2, int colDm2){
	
	for(int i = 0; i < rowDm1; i++){
		for(int j = 0; j < colDm2; j++){
			double sum = 0;
			for(int k = 0; k < colDm1; k++){					
				sum  = sum + vect1[i][k] * vect2[k][j];
			}			
			result[i][j] = sum;
		}
	}
}

/*
	The function for adding three matrices in a row  
*/
void matrixAdditionThree(double ** result,  double ** vect1, double ** vect2, double ** vect3, int rowDm, int colDm){

	for(int i = 0; i < rowDm; i++)
		for(int j = 0; j < colDm; j++)
			result[i][j] = vect1[i][j] + vect2[i][j] + vect3[i][j];
}

/*
	The function for copying a result vector to a target vector 
*/
void assignResult(double ** result, double ** target, int rowDm, int colDm){
	for(int i = 0; i < rowDm; i++)
		for(int j = 0; j < colDm; j++)
			target[i][j] = result[i][j];
}

/*
	The function for printing a matrix content 
*/
void printMatrix(double **matrix, int rowDm, int colDm){

	for(int i = 0; i < rowDm; i ++){
		for(int j = 0; j < colDm; j++)
			printf("%f \t", matrix[i][j]);	
		printf("\n");	
	}
}


