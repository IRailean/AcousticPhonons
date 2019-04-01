#include <iostream>
#include <fstream>
#include <conio.h>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include "mkl_lapack.h" 
#include "oneLayerShear.h"


/* Constants declarations begin*/

constexpr double Pi = 3.14159265359;

/* Elastic modules */
double static c44;

/* Lattice constants */
double static latticeConstant;

/* Densities*/
double static density;

// Number of dots
// This shows how much dots we place along the axis X3 
constexpr int dotsAmount = 100;

// Structure width (nm)
double static width = 6;

// Step size along X3 axis
double static step = width / dotsAmount;

// Wave vector number of values
static int numOfWaveVectorValues;

#define N dotsAmount
#define LDA N
#define LDVL N
#define LDVR N

/* Constants declarations end */

int solveOneLayerShear(std::string element)
{
	setConstants(element);

	// Wave vector values
	double qFirst = 0;
	double qLast = (2 * Pi) / latticeConstant;
	double qStep = (Pi / 128) / latticeConstant;
	numOfWaveVectorValues = (qLast - qFirst) / qStep;

	/* Calculate wave vector values and coefficients */
	double coefU_previous = -c44 / (density * pow(step, 2));
	double coefU_next = -c44 / (density * pow(step, 2));
	double *coefU_current = new double[numOfWaveVectorValues];

	float *qValues = new float[numOfWaveVectorValues];
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		qValues[i] = qFirst + i * qStep;
	}
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		coefU_current[i] = 2 * c44 * (2 * pow(sin(qValues[i] * latticeConstant / 2), 2) / pow(latticeConstant, 2) + 1 / pow(step, 2)) / density;
	}

	/* Init matrix with zero values */
	float ***matrix = init3DMatrixWithZeros();
	
	/* Fill matrix with coefficients from equations */
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		for (int i = 1; i < dotsAmount - 1; i++)
		{
			matrix[k][i][i - 1] = coefU_next;
			matrix[k][i][i] = coefU_current[k];
			matrix[k][i][i + 1] = coefU_previous;
		}
	}
	/* Add boundary conditions*/
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		matrix[k][0][0] = coefU_current[k];
		matrix[k][0][1] = coefU_previous + coefU_next;
		matrix[k][dotsAmount - 1][dotsAmount - 2] = coefU_next + coefU_previous;
		matrix[k][dotsAmount - 1][dotsAmount - 1] = coefU_current[k];
	}
	/* Calculate eigenvalues using sgeev */
	float** eigenvalues = calculateEigenvalues(matrix);

	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		eigenvalues[i][0] = qValues[i];
	}
	float** energies = eigenvaluesToEnergy(eigenvalues, numOfWaveVectorValues, dotsAmount + 1);

	/* Calculate group velocities */
	float** groupVelocities = new float*[numOfWaveVectorValues + 1];
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		groupVelocities[i] = new float[dotsAmount + 1];
	}
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		groupVelocities[i][0] = qValues[i];
	}

	for (int i = 1; i < numOfWaveVectorValues - 1; i++)
	{
		for (int j = 1; j <  dotsAmount + 1; j++)
		{
			groupVelocities[i][j] = (1.0/0.658) * (energies[i + 1][j] - energies[i - 1][j]) / (2 * qStep);
		}
	}
	for (int j = 1; j < dotsAmount + 1; j++)
	{
		groupVelocities[0][j] = (1.0 / 0.658) * (energies[1][j] - energies[0][j]) / qStep;
	}
	for (int j = 1; j < dotsAmount + 1; j++)
	{
		groupVelocities[numOfWaveVectorValues - 1][j] = (1.0 / 0.658) * (energies[numOfWaveVectorValues - 2][j] - energies[numOfWaveVectorValues - 1][j]) / qStep;
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
static void setConstants(std::string element)
{
	if (strcmp("Ge", element.c_str()) == 0)
	{
		c44 = 67.7;
		latticeConstant = 0.5658;
		density = 5.323;
	}
	else if (strcmp("Si", element.c_str()) == 0)
	{
		c44 = 79.6;
		latticeConstant = 0.5431;
		density = 2.329;
	}
	else
	{
		std::cout << "Wrong choice !";
	}
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
		for (int j = 0; j < dotsAmount + 1; j++)
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