#include <iostream>
//#include "oneLayerShear.h"
//#include "threeLayerShear.h"
#include "oneLayerASSA.h"

int main()
{
	int retVal = 0;
	//retVal = solveOneLayerShear("Ge");
	solveOneLayerASSA("Si");

	return retVal;
}