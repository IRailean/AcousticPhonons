#include <iostream>
#include <fstream>
#include <conio.h>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <algorithm>

#include "mkl_lapack.h" 
#include "threeLayerShear.h"


/* Constants declarations begin*/

constexpr double Pi = 3.14159265359;

/* Elastic modules */
static double  c44_Ge = 67.7;
static double c44_Si = 79.6;

/* Lattice constants */
static double  latticeConstant_Ge = 0.5658;
static double latticeConstant_Si = 0.5431;

/* Densities*/
static double density_Ge = 5.323;
static double density_Si = 2.329;

// Number of dots
// This shows how much dots we place along the axis X3 
constexpr int dotsAmount = 50;

// Structure width (nm)
static double  width = 6;

// Step size along X3 axis
static double step = width / dotsAmount;

// Wave vector number of values
// Wave vector values
static double qFirst = 0;
static double qLast = (2 * Pi) / latticeConstant_Ge;
static double qStep = (Pi / 128) / latticeConstant_Ge;
static int numOfWaveVectorValues = (qLast - qFirst) / qStep;

#define N dotsAmount
#define LDA N
#define LDVL N
#define LDVR N

/* Constants declarations end */

int solveThreeLayerShear(double width_Ge, double width_Si)
{
	/* Allocate memory for coefficients */
	//  coefficients depend of 2 type of variables - wave vector (k) [0 - numOfWaveVectorValues] and constants [0 - dotsAmount]
	double **coefU_previous = new double*[numOfWaveVectorValues];
	double **coefU_next = new double*[numOfWaveVectorValues];
	double **coefU_current = new double*[numOfWaveVectorValues];

	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		coefU_previous[k] = new double[dotsAmount];
		coefU_current[k] = new double[dotsAmount];
		coefU_next[k] = new double[dotsAmount];
	}

	/* Determine constants in three layer structure using linear function with deltaL nm transition */
	// Calculate borders position due to width of Ge and Si layers
	double deltaL = 0.5;
	double firstLeftBorder = width_Ge - deltaL / 2;
	double firstRightBorder = width_Ge + deltaL / 2;
	double secondLeftBorder = width_Si + width_Ge - deltaL / 2;
	double secondRightBorder = width_Si + width_Ge + deltaL / 2;

	// alpha is dc44/dx3, where x3 is the axis along which layers are placed
	double *c44 = new double[dotsAmount];
	double *latticeConstant = new double[dotsAmount];
	double *density = new double[dotsAmount];
	double *alpha = new double[dotsAmount];

	// Set constants in dependency of dot position, i.e where the dot is placed in structure
	for (int i = 0; i < dotsAmount; i++)
	{
		double dotPos = i * width / dotsAmount;
		if (dotPos < firstLeftBorder)
		{
			alpha[i] = 0;
			c44[i] = c44_Ge;
			density[i] = density_Ge;
			latticeConstant[i] = latticeConstant_Ge;
		}
		else if (dotPos >= firstLeftBorder && dotPos < firstRightBorder)
		{
			double deltaZ = firstRightBorder - dotPos;
			alpha[i] = (c44_Ge - c44_Si) / -deltaL;
			c44[i] = (deltaZ*c44_Ge + (deltaL - deltaZ)*c44_Si) / deltaL;
			density[i] = (deltaZ*density_Ge + (deltaL - deltaZ)*density_Si) / deltaL;
			latticeConstant[i] = (deltaZ*latticeConstant_Ge + (deltaL - deltaZ)*latticeConstant_Si) / deltaL;
		}
		else if (dotPos >= firstRightBorder && dotPos < secondLeftBorder)
		{
			alpha[i] = 0;
			c44[i] = c44_Si;
			density[i] = density_Si;
			latticeConstant[i] = latticeConstant_Si;
		}
		else if (dotPos > secondLeftBorder && dotPos < secondRightBorder)
		{
			double deltaZ = secondRightBorder - dotPos;
			alpha[i] = (c44_Si - c44_Ge) / -deltaL;
			c44[i] = (deltaZ*c44_Si + (deltaL - deltaZ)*c44_Ge) / deltaL;
			density[i] = (deltaZ*density_Si + (deltaL - deltaZ)*density_Ge) /deltaL;
			latticeConstant[i] = (deltaZ*latticeConstant_Si + (deltaL - deltaZ)*latticeConstant_Ge) / deltaL;
		}
		else if (dotPos >= secondRightBorder)
		{
			alpha[i] = 0;
			c44[i] = c44_Ge;
			density[i] = density_Ge;
			latticeConstant[i] = latticeConstant_Ge;
		}
	}
	/* Set wave vector values */
	float *qValues = new float[numOfWaveVectorValues];
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		qValues[i] = qFirst + i * qStep;
	}
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		for (int i = 0; i < dotsAmount; i++)
		{
			coefU_previous[k][i] = -(alpha[i] / (density[i] * 2 * step) + c44[i] / (density[i] * pow(step, 2)));
			coefU_current[k][i] = (2 * c44[i] / density[i]) * ((2 * pow(sin(qValues[k] * latticeConstant[i] / 2), 2) / pow(latticeConstant[i], 2)) + 1 / pow(step, 2));
			coefU_next[k][i] = -(-alpha[i] / (density[i] * 2 * step) + c44[i] / (density[i] * pow(step, 2)));
		}
	}
	
	
	/* Init matrix with zero values */
	float ***matrix = init3DMatrixWithZeros();

	/* Fill matrix with coefficients from equations */
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		for (int i = 1; i < dotsAmount - 1; i++)
		{
			matrix[k][i][i - 1] = coefU_next[k][i];
			matrix[k][i][i] = coefU_current[k][i];
			matrix[k][i][i + 1] = coefU_previous[k][i];
		}
	}
	/* Add boundary conditions*/
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		matrix[k][0][0] = coefU_current[k][0];
		matrix[k][0][1] = coefU_previous[k][0] + coefU_next[k][0];
		matrix[k][dotsAmount - 1][dotsAmount - 2] = coefU_next[k][dotsAmount - 1] + +coefU_previous[k][dotsAmount - 1];
		matrix[k][dotsAmount - 1][dotsAmount - 1] = coefU_current[k][dotsAmount - 1];
	}
	/* Calculate eigenvalues using sgeev */
	float** eigenvalues = calculateEigenvalues(matrix);

	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		eigenvalues[i][0] = qValues[i];
	}
	float** energies = eigenvaluesToEnergy(eigenvalues, numOfWaveVectorValues, dotsAmount + 1);
	/* Calculate group velocities */
	float** groupVelocities = new float*[numOfWaveVectorValues - 1];
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		groupVelocities[i] = new float[dotsAmount];
	}
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		groupVelocities[i][0] = qValues[i];
	}
	for (int i = 0; i < numOfWaveVectorValues - 1; i++)
	{
		for (int j = 1; j < dotsAmount; j++)
		{
			groupVelocities[i][j] = (energies[i + 1][j] - energies[i][j]) / (qValues[i + 1] - qValues[i]);
		}
	}
	/* Output results */
	std::string filenameEnergies = "energies.csv";
	writeEnergiesToFile(filenameEnergies, energies);
	std::string filenameGroupVelocities = "groupVelocities.csv";
	writeGroupVelocitiesToFile(filenameGroupVelocities, groupVelocities);


	std::cout << "Eigenvalues succesfully calculated, energies stored in " + filenameEnergies + " , group velocities stored in " + filenameGroupVelocities + " press any key to continue...";
	delete[] matrix;
	delete[] energies;
	delete qValues;
	_getch();
	return 0;
}
static float*** init3DMatrixWithZeros()
{
	float ***matrix = new float**[numOfWaveVectorValues];
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		matrix[k] = new float*[dotsAmount];
		for (int i = 0; i < dotsAmount; i++)
		{
			matrix[k][i] = new float[dotsAmount];
			for (int j = 0; j < dotsAmount; j++)
			{
				matrix[k][i][j] = 0;
			}
		}
	}
	return matrix;
}
static float** calculateEigenvalues(float*** matrix)
{
	int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
	float wkopt;
	float *work;
	float wr[N], wi[N], vl[N*N], vr[N*N];

	float **eigenvalues = new float*[numOfWaveVectorValues];
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		eigenvalues[i] = new float[dotsAmount + 1];
	}

	/* Calculate eigenvalues for matrices for each wave vector value */
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		lwork = -1;
		/* temp array is needed to convert from 2D matrix to 1D representation */
		float *temp = new float[dotsAmount*dotsAmount];

		for (int i = 0; i < dotsAmount*dotsAmount; i++)
		{
			temp[i] = matrix[k][i / dotsAmount][i % dotsAmount];
		}


		sgeev((char*)"Vectors", (char*)"Vectors", &n, temp, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info);
		lwork = (int)wkopt;
		work = new float[lwork];
		sgeev((char*)"Vectors", (char*)"Vectors", &n, temp, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
		if (info > 0)
		{
			std::cout << "Failed to calculate";
			exit(1);
		}
		std::sort(wr, wr + dotsAmount);
		for (int i = 1; i < dotsAmount + 1; i++)
			eigenvalues[k][i] = wr[i - 1];

		free((void*)work);
		delete temp;
	}
	return eigenvalues;
}
static void writeEnergiesToFile(std::string filename, float** energies)
{
	std::ofstream myFile;
	myFile.open(filename);

	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		for (int j = 0; j < dotsAmount + 1; j++)
		{
			myFile << energies[i][j] << ',';
		}
		myFile << '\n';
	}
	myFile.close();
}
static void writeGroupVelocitiesToFile(std::string filename, float** groupVelocities)
{
	std::ofstream myFile;
	myFile.open(filename);

	for (int i = 0; i < numOfWaveVectorValues - 1; i++)
	{
		for (int j = 0; j < dotsAmount; j++)
		{
			myFile << groupVelocities[i][j] << ',';
		}
		myFile << '\n';
	}
	myFile.close();
}
static float** eigenvaluesToEnergy(float** eigenvalues, int sizeA, int sizeB)
{
	float** energies = new float*[sizeA];
	for (int i = 0; i < sizeA; i++)
	{
		energies[i] = new float[sizeB];
	}
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		energies[i][0] = eigenvalues[i][0];
	}
	for (int i = 0; i < sizeA; i++)
	{
		for (int j = 1; j < sizeB; j++)
		{
			energies[i][j] = 0.658*sqrt(eigenvalues[i][j]);
		}
	}
	delete[] eigenvalues;
	return energies;
}