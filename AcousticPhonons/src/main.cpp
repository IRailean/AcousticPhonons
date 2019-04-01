#include <iostream>
//#include "oneLayerShear.h"
#include "threeLayerShear.h"

int main()
{
	int retVal = 0;
	//retVal = solveOneLayerShear("Ge");
    solveThreeLayerShear(2, 2);

	return retVal;
}