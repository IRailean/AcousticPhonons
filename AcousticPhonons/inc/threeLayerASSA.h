#include <string>

#ifndef THREE_LAYER_ASSA_H
#define THREE_LAYER_ASSA_H

// Code body for header file


/* Function prototypes begin */

int solveThreeLayerASSA(double, double);
static double*** init3DMatrixWithZeros();
static double** calculateEigenvalues(double***);
static void writeEnergiesToFile(std::string, double**);
static void writeGroupVelocitiesToFile(std::string, double**);
static double** eigenvaluesToEnergy(double**, int, int);

/* Function prototypes end */

#endif