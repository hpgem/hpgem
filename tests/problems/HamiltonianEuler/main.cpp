//
//  main.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 6/21/13.
//
//

#include <iostream>
#include <string>
using namespace std;

#include <petscmat.h>
#include <petscksp.h>


#include "HEuler.hpp"
    // ============================================================
static char help[] = "Petsc Help";

int main(int argc, char **argv)
{
    
    
   PetscInitialize(&argc,&argv,PETSC_NULL,help);

    const string fileName="blabla.out";
    
        /// enum  SolutionType 		{INCOMPRESSIBLE_WALLS, INCOMPRESSIBLE_ONETHIRDPERIODIC, INCOMPRESSIBLE_PERIODIC, COMPRESSIBLE_WALLS, COMPRESSIBLE_PERIODIC};

        //    unsigned int numberOfUnknowns, unsigned int numberOfBasisFunctions, unsigned int  numberOfTimeLevels=1,
    
    HEulerConfigurationData config(4, 1, 1, fileName, HEulerConfigurationData::INCOMPRESSIBLE_WALLS);
    HEulerGlobalVariables globalVar;

    
    
    if (config.solutionType_==HEulerConfigurationData::INCOMPRESSIBLE_PERIODIC||config.solutionType_==HEulerConfigurationData::COMPRESSIBLE_PERIODIC || config.solutionType_==HEulerConfigurationData::COMPRESSIBLE_WALLS)
    {
        config.lx_ = 1;
        config.ly_ = 1;
        config.lz_ = 1;
    }
    else if (config.solutionType_==HEulerConfigurationData::INCOMPRESSIBLE_ONETHIRDPERIODIC)
    {
        config.lx_ = 2*3.14159265359;
        config.ly_ = 2*3.14159265359;
        config.lz_ = 2*3.14159265359;
    }
    else if (config.solutionType_==HEulerConfigurationData::INCOMPRESSIBLE_WALLS)
    {
        config.lx_ = 2*3.14159265359;
        config.ly_ = 3.14159265359;
        config.lz_ = 3.14159265359;
    }
    config.nx_ = 10;
    config.ny_ = 10;
    config.nz_ = 10;
    
    HEuler eulerProblem(&globalVar, &config);

    eulerProblem.initialiseMesh();

    
    eulerProblem.initialConditions();

    eulerProblem.output(0);
    
    
    
    if (config.solutionType_==HEulerConfigurationData::COMPRESSIBLE_PERIODIC || config.solutionType_==HEulerConfigurationData::COMPRESSIBLE_WALLS)
    {
         eulerProblem.createCompressibleSystem();
    }
    else
    {
        eulerProblem.createIncompressibleSystem();
    }
    
    eulerProblem.output(0.0001);
    
    eulerProblem.solve();
    
    
        //eulerProblem.output(2.0);
    
    return 0;
}