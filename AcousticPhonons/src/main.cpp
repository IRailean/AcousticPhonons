#include <iostream>
//#include "oneLayerShear.h"
//#include "threeLayerShear.h"
//#include "oneLayerASSA.h"
#include "threeLayerASSA.h"

int main()
{
	int retVal = 0;
	//retVal = solveOneLayerShear("Ge");
	solveThreeLayerASSA(2, 2);

	return retVal;
}