#include <iostream>
#include <fstream>
#include <conio.h>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <algorithm>

#include "mkl_lapack.h" 
#include "oneLayerASSA.h"


/* Constants declarations begin*/

constexpr double Pi = 3.14159265359;

/* Elastic modules */
static double  c44_Ge = 67.7;
static double c44 = 79.6;
static double c11 = 166;
static double c12 = 64;

/* Lattice constants */
static double  latticeConstant_Ge = 0.5658;
static double latticeConstant = 0.5431;

/* Densities*/
static double density_Ge = 5.323;
static double density = 2.329;

// Number of dots
// This shows how much dots we place along the axis X3 
constexpr int dotsAmount = 200;

// Structure width (nm)
static double  width = 6;

// Step size along X3 axis
static double step = width / (double)dotsAmount;

// Wave vector number of values
// Wave vector values
static double divider = 256;
static double qFirst = 0;
static double qLast = (2 * Pi) / latticeConstant;
static double qStep = (Pi / divider) / latticeConstant;
static int numOfWaveVectorValues = (qLast - qFirst) / qStep;

#define N dotsAmount
#define LDA N
#define LDVL N
#define LDVR N

/* Constants declarations end */

int solveOneLayerASSA(std::string matter)
{
	/* Allocate memory for coefficients */

	double *coefU1_main_previous = new double[numOfWaveVectorValues];
	double *coefU1_main_next = new double[numOfWaveVectorValues];
	double *coefU1_main_current = new double[numOfWaveVectorValues];
	double *coefU3_secondary_previous = new double[numOfWaveVectorValues];
	double *coefU3_secondary_next = new double[numOfWaveVectorValues];

	double *coefU3_main_previous = new double[numOfWaveVectorValues];
	double *coefU3_main_next = new double[numOfWaveVectorValues];
	double *coefU3_main_current = new double[numOfWaveVectorValues];
	double *coefU1_secondary_previous = new double[numOfWaveVectorValues];
	double *coefU1_secondary_next = new double[numOfWaveVectorValues];

	/* Fill coefficients */
	divider *= 2.0;

	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		coefU1_main_current[k] = (4.0 * c11 * pow(sin(k * (Pi / divider)), 2) / (pow(latticeConstant, 2)) + 2.0 * c44 / (step * step)) * (1.0 / density);
		coefU1_main_next[k] = -c44 / (density * step * step);
		coefU1_main_previous[k] = -c44 / (density * step * step);
		coefU3_secondary_next[k] =  (-c12 * sin(k * (Pi / divider) * 2.0) / (2.0 * latticeConstant * step) - c44 * sin(k * (Pi / divider) * 2.0) / (2.0 * latticeConstant * step)) * (1.0 / density);
		coefU3_secondary_previous[k] = (c12 * sin(k * (Pi / divider) * 2.0) / (2.0 * latticeConstant * step) + c44 * sin(k * (Pi / divider) * 2.0) / (2.0 * latticeConstant * step)) * (1.0 / density);

		coefU3_main_current[k] = (4.0 * c44 * pow(sin(k * (Pi / divider)), 2.0) / pow(latticeConstant, 2.0) + 2.0 * c11 / (step * step)) * (1.0 / density);
		coefU3_main_next[k] = -c11 / (density * step * step);
		coefU3_main_previous[k] = -c11 / (density * step * step);
		coefU1_secondary_next[k] = (c12 * sin(k * (Pi / divider) * 2.0) / (2.0 * latticeConstant * step) + c44 * sin(k * (Pi / divider) * 2.0) / (2.0 * latticeConstant * step)) * (1.0 / density);
		coefU1_secondary_previous[k] = (-c12 * sin(k * (Pi / divider) * 2.0) / (2.0 * latticeConstant * step) - c44 * sin(k * (Pi / divider) * 2.0) / (2.0 * latticeConstant * step)) * (1.0 / density);
	}
	
	/* Init matrix with zero values */
	double ***matrix = init3DMatrixWithZeros();

	/* Fill matrix with coefficients from equations */
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		for (int i = 1; i < dotsAmount / 2 - 1; i++)
		{
			matrix[k][i][i - 1] = coefU1_main_previous[k];
			matrix[k][i][i] = coefU1_main_current[k];
			matrix[k][i][i + 1] = coefU1_main_next[k];

			matrix[k][i][dotsAmount / 2 + i - 1] = coefU3_secondary_previous[k];
			matrix[k][i][dotsAmount / 2 + i + 1] = coefU3_secondary_next[k];
		}
		for (int i = dotsAmount / 2 + 1; i < dotsAmount - 1; i++)
		{
			matrix[k][i][i] = coefU3_main_current[k];
			matrix[k][i][i - 1] = coefU3_main_previous[k];
			matrix[k][i][i + 1] = coefU3_main_next[k];

			matrix[k][i][i - dotsAmount / 2 - 1] = coefU1_secondary_previous[k];
			matrix[k][i][i - dotsAmount / 2 + 1] = coefU1_secondary_next[k];
		}
	}

	/* Additional coefficients */
	double *coefFirst = new double[numOfWaveVectorValues];
	double *coefSecond = new double[numOfWaveVectorValues];
	
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		coefFirst[k] = 2.0 * step * c12 / (latticeConstant * c11) * sin(k * (Pi / divider) * 2.0);
		coefSecond[k] = 2.0 * step / latticeConstant * sin(k * (Pi / divider) * 2.0);
	}

	/* Add border conditions*/
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		matrix[k][0][0] = coefU1_main_current[k] - coefU3_secondary_previous[k] * coefFirst[k];
		matrix[k][0][1] = coefU1_main_next[k] + coefU1_main_previous[k];
		matrix[k][0][dotsAmount / 2] = coefU1_main_previous[k] * coefSecond[k];
		matrix[k][0][dotsAmount / 2 + 1] = coefU3_secondary_next[k] + coefU3_secondary_previous[k];

		matrix[k][dotsAmount / 2 - 1][dotsAmount / 2 - 2] = coefU1_main_next[k] + coefU1_main_previous[k];
		matrix[k][dotsAmount / 2 - 1][dotsAmount / 2 - 1] = coefU1_main_current[k] + coefU3_secondary_next[k] * coefFirst[k];
		matrix[k][dotsAmount / 2 - 1][dotsAmount - 2] = coefU3_secondary_next[k] + coefU3_secondary_previous[k];
		matrix[k][dotsAmount / 2 - 1][dotsAmount - 1] = -coefU1_main_next[k] * coefSecond[k];

		matrix[k][dotsAmount / 2][0] = -coefU3_main_previous[k] * coefFirst[k];
		matrix[k][dotsAmount / 2][1] = coefU1_secondary_next[k] + coefU1_secondary_previous[k];
		matrix[k][dotsAmount / 2][dotsAmount / 2] = coefU3_main_current[k] + coefU1_secondary_previous[k] * coefSecond[k];
		matrix[k][dotsAmount / 2][dotsAmount / 2 + 1] = coefU3_main_next[k] + coefU3_main_previous[k];

		matrix[k][dotsAmount - 1][dotsAmount / 2 - 2] = coefU1_secondary_next[k] + coefU1_secondary_previous[k];
		matrix[k][dotsAmount - 1][dotsAmount / 2 - 1] = coefU3_main_next[k] * coefFirst[k];
		matrix[k][dotsAmount - 1][dotsAmount - 2] = coefU3_main_next[k] + coefU3_main_previous[k];
		matrix[k][dotsAmount - 1][dotsAmount - 1] = coefU3_main_current[k] - coefU1_secondary_next[k] * coefSecond[k];
	}
	/* Calculate eigenvalues using dgeev */
	double** eigenvalues = calculateEigenvalues(matrix);

	double *qValues = new double[numOfWaveVectorValues];

	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		qValues[i] = qFirst + i * qStep;
		eigenvalues[i][0] = qValues[i];
	}
	double** energies = eigenvaluesToEnergy(eigenvalues, numOfWaveVectorValues, dotsAmount + 1);

	/* Calculate group velocities */
	double** groupVelocities = new double*[numOfWaveVectorValues + 1];
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		groupVelocities[i] = new double[dotsAmount + 1];
	}
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		groupVelocities[i][0] = qValues[i];
	}
	for (int i = 1; i < numOfWaveVectorValues - 1; i++)
	{
		for (int j = 1; j < dotsAmount + 1; j++)
		{
			groupVelocities[i][j] = (1.0 / 0.658) * (energies[i + 1][j] - energies[i - 1][j]) / (2.0 * qStep);
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
static double*** init3DMatrixWithZeros()
{
	double ***matrix = new double**[numOfWaveVectorValues];
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		matrix[k] = new double*[dotsAmount];
		for (int i = 0; i < dotsAmount; i++)
		{
			matrix[k][i] = new double[dotsAmount];
			for (int j = 0; j < dotsAmount; j++)
			{
				matrix[k][i][j] = 0;
			}
		}
	}
	return matrix;
}
static double** calculateEigenvalues(double*** matrix)
{
	int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
	double wkopt;
	double *work;
	double wr[N], wi[N];
	double *vl = new double[N*N];
	double *vr = new double[N*N];

	double **eigenvalues = new double*[numOfWaveVectorValues];
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		eigenvalues[i] = new double[dotsAmount + 1];
	}

	/* Calculate eigenvalues for matrices for each wave vector value */
	for (int k = 0; k < numOfWaveVectorValues; k++)
	{
		lwork = -1;
		/* temp array is needed to convert from 2D matrix to 1D representation */
		double *temp = new double[dotsAmount*dotsAmount];

		for (int i = 0; i < dotsAmount*dotsAmount; i++)
		{
			temp[i] = matrix[k][i / dotsAmount][i % dotsAmount];
		}


		dgeev((char*)"Vectors", (char*)"Vectors", &n, temp, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info);
		lwork = (int)wkopt;
		work = new double[lwork];
		dgeev((char*)"Vectors", (char*)"Vectors", &n, temp, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
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
	delete vl;
	delete vr;
	return eigenvalues;
}
static void writeEnergiesToFile(std::string filename, double** energies)
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
static void writeGroupVelocitiesToFile(std::string filename, double** groupVelocities)
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
static double** eigenvaluesToEnergy(double** eigenvalues, int sizeA, int sizeB)
{
	double** energies = new double*[sizeA];
	for (int i = 0; i < sizeA; i++)
	{
		energies[i] = new double[sizeB];
	}
	for (int i = 0; i < numOfWaveVectorValues; i++)
	{
		energies[i][0] = eigenvalues[i][0];
	}
	for (int i = 0; i < sizeA; i++)
	{
		for (int j = 1; j < sizeB; j++)
		{
			energies[i][j] = 0.6582119514*sqrt(eigenvalues[i][j]);
		}
	}
	delete[] eigenvalues;
	return energies;
}