#ifndef ONE_LAYER_SHEAR_H
#define ONE_LAYER_SHEAR_H

// Code body for header file


/* Function prototypes begin */

int solveOneLayerShear();
float*** init3DMatrixWithZeros();
float** calculateEigenvalues(float***);
void writeEnergiesToFile(std::string, float**);
void writeGroupVelocitiesToFile(std::string, float**);
float** eigenvaluesToEnergy(float**, int, int);

/* Function prototypes end */

#endif