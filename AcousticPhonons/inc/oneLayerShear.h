#include <string>

#ifndef ONE_LAYER_SHEAR_H
#define ONE_LAYER_SHEAR_H

// Code body for header file


/* Function prototypes begin */

int solveOneLayerShear(std::string);
static void setConstants(std::string element);
static float*** init3DMatrixWithZeros();
static float** calculateEigenvalues(float***);
static void writeEnergiesToFile(std::string, float**);
static void writeGroupVelocitiesToFile(std::string, float**);
static float** eigenvaluesToEnergy(float**, int, int);

/* Function prototypes end */

#endif