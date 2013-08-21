//
//  main.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 6/21/13.
//
//

#include <iostream>
#include <cstring>

using namespace std;

#include "HEuler.hpp"
    // ============================================================
static char help[] = "Petsc Help";

int main(int argc, char **argv)
{
    
    
    PetscInitialize(&argc,&argv,PETSC_NULL,help);

    const string fileName="blabla.out";
        //    unsigned int numberOfUnknowns, unsigned int numberOfBasisFunctions, unsigned int  numberOfTimeLevels=1,
    HEulerConfigurationData config(4, 11, 1, fileName, HEulerConfigurationData::COMPRESSIBLE_PERIODIC);
    HEulerGlobalVariables globalVar;
    HEuler eulerProblem(&globalVar, &config);
    
    eulerProblem.initialiseMesh();

    
        eulerProblem.initialConditions();
    
        eulerProblem.createTheCompressibleSystem();

        eulerProblem.output();
    
        eulerProblem.solve();
    
    
        eulerProblem.output(2.0);
    
    return 0;
}