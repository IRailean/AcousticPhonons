#include <iostream>
#include "oneLayerShear.h"
#include "threeLayerShear.h"

int main()
{
	int retVal = 0;
	//retVal = solveOneLayerShear("Si");
	solveThreeLayerShear("Ge");

	return retVal;
}