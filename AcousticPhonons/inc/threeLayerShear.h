#include <string>

#ifndef THREE_LAYER_SHEAR_H
#define THREE_LAYER_SHEAR_H

// Code body for header file


/* Function prototypes begin */

int solveThreeLayerShear(double, double);
static float*** init3DMatrixWithZeros();
static float** calculateEigenvalues(float***);
static void writeEnergiesToFile(std::string, float**);
static void writeGroupVelocitiesToFile(std::string, float**);
static float** eigenvaluesToEnergy(float**, int, int);

/* Function prototypes end */

#endif