/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//this file should contain all relevant information about how the integrands look like and what problem is solved

#define _USE_MATH_DEFINES

#include "BaseExtended.h"

#include <cstdlib>
#include <iostream>
#include "math.h"
#include <ctime>

#include "DGMaxDim.h"

#include "Algorithms/DGMaxEigenValue.h"
#include "Algorithms/DGMaxHarmonic.h"
#include "Algorithms/DGMaxTimeIntegration.h"
#include "Algorithms/DivDGMaxEigenValue.h"
#include "Algorithms/DivDGMaxHarmonic.h"

#include "ProblemTypes/Harmonic/SampleHarmonicProblems.h"
#include "ProblemTypes/Time/SampleTestProblems.h"

/**
 * This class should provide problem specific information about the maxwell equations.
 */
class DGMax : public hpGemUIExtentions
{
public:
    
    DGMax(Base::GlobalData * const globalConfig, Base::ConfigurationData* elementConfig)
            : hpGemUIExtentions(globalConfig, elementConfig)
    {
    }
    
    /**
     * set up the mesh and complete initialisation of the global data and the configuration data
     */


    void createCubeMesh(std::size_t subdivisions)
    {
        Geometry::PointPhysical<DIM> bottomLeft, topRight;
        std::vector<std::size_t> numElementsOneD (DIM);
        // Configure each dimension of the unit cube/square
        for (std::size_t i = 0; i < DIM; ++i) {
            bottomLeft[i] = 0;
            topRight[i] = 1;
            numElementsOneD[i] = subdivisions;
        }

        auto mesh = new Base::MeshManipulator<DIM>(getConfigData(), Base::BoundaryType::PERIODIC,
                                                   Base::BoundaryType::PERIODIC, Base::BoundaryType::PERIODIC, getConfigData()->polynomialOrder_, 0, 2, 3, 1, 1);
        mesh->createTriangularMesh(bottomLeft, topRight, numElementsOneD);

        for (Base::MeshManipulator<DIM>::ElementIterator it = mesh->elementColBegin(Base::IteratorType::GLOBAL); it != mesh->elementColEnd(Base::IteratorType::GLOBAL); ++it)
        {
            (*it)->setUserData(new ElementInfos(**it));
        }
        addMesh(mesh);
    }
    
    void readMesh(std::string fileName)
    {
        auto mesh = new Base::MeshManipulator<DIM>(getConfigData(), Base::BoundaryType::PERIODIC,
                                                   Base::BoundaryType::PERIODIC, Base::BoundaryType::PERIODIC, getConfigData()->polynomialOrder_, 0, 2, 3, 1, 1);
        mesh->readMesh(fileName);
        addMesh(mesh);
        for (Base::MeshManipulator<DIM>::ElementIterator it = mesh->elementColBegin(Base::IteratorType::GLOBAL); it != mesh->elementColEnd(Base::IteratorType::GLOBAL); ++it)
        {
            (*it)->setUserData(new ElementInfos(**it));
        }
        
    }
};

auto& numElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", true);
auto& p = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);
auto& meshFile = Base::register_argument<std::string>('m', "meshFile", "The hpgem meshfile to use", false);

/**
 * main routine. Suggested usage: make a problem, solve it, tell the user your results.
 */
int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    std::cout<<"This is the parallel version"<<std::endl;

    logger.assert_debug( DIM >= 2 && DIM <= 3, "Can only handle 2D and 3D problems.");

    //set up timings
    time_t start, end, initialised, solved;
    time(&start);
    logger(INFO, "using % elements", std::pow(numElements.getValue(), 3) * 5);
    logger(INFO, "using polynomial order: %", p.getValue());
    logger(INFO, "Compiled to use % dimensions.", DIM);
    //set up problem and decide flux type
    //DGMax problem(new MaxwellData(numElements.getValue(), p.getValue()), new Base::ConfigurationData(DIM, 1, p.getValue(), 1), new MatrixAssemblyIP);
    const std::size_t numberOfTimeLevels = 1;
    Base::GlobalData* const globalData = new Base::GlobalData();
    globalData->numberOfTimeLevels_ = numberOfTimeLevels;
    // TODO: LC: this should be determined by the discretization, but this is
    // currently not possible yet, as we get a dependency loop (discretization
    // requires DGMax, which requires the configurationData, which would then
    // require the discretization).
    const std::size_t numberOfUnknowns = 2;
    Base::ConfigurationData* const configData = new Base::ConfigurationData(DIM, numberOfUnknowns, p.getValue(), numberOfTimeLevels);
    try
    {
        double stab = numElements.getValue() * (p.getValue() + 1) * (p.getValue() + 3);
        DivDGMaxDiscretization::Stab divStab;
        // Values from the Jelmer fix code.
        divStab.stab1 = 100 * numElements.getValue() / sqrt(2.0);
        // Note: for 2D harmonic it looks like that we need 10 instead of 0.01.
        divStab.stab2 = 0.01 / (numElements.getValue() * sqrt(2.0));
        divStab.stab3 = 10.0 * numElements.getValue() / sqrt(2.0);

        DGMax base(globalData, configData);

        if (meshFile.isUsed())
        {
            base.readMesh(meshFile.getValue());
        }
        else
        {
            // Temporary fall back for easy testing.
            base.createCubeMesh(numElements.getValue());
        }
        //base.createCentaurMesh(std::string("SmallIW_Mesh4000.hyb"));
        //base.createCentaurMesh(std::string("BoxCylinder_Mesh6000.hyb"));
        //TODO: LC: this does seem rather arbitrary and should probably be done by the solver
        base.setNumberOfTimeIntegrationVectorsGlobally(21);

        ///////////////////
        // Harmonic code //
        ///////////////////

//        DGMaxHarmonic harmonicSolver (base);
//        DivDGMaxHarmonic harmonicSolver (base);
//
//        SampleHarmonicProblems problem (SampleHarmonicProblems::SARMANY2010, 1);
//        harmonicSolver.solve(problem, divStab);
//        std::cout << "L2 error " << harmonicSolver.computeL2Error(problem) << std::endl;
//        auto errors = harmonicSolver.computeError({DGMaxDiscretization::L2, DGMaxDiscretization::HCurl}, problem);
//        std::cout << "L2 error    " << errors[DGMaxDiscretization::L2] << std::endl;
//        std::cout << "HCurl error " << errors[DGMaxDiscretization::HCurl] << std::endl;


        /////////////////////
        // Eigenvalue code //
        /////////////////////

//        DGMaxEigenValue solver (base);
        DivDGMaxEigenValue solver (base);
        EigenValueProblem input;
        solver.solve(input, divStab);

        ///////////////////////////
        // Time dependent solver //
        ///////////////////////////

//        DGMaxTimeIntegration timeSolver (base);
//        SampleTestProblems testProblem (SampleTestProblems::SARMANY2013);
//        TimeIntegrationParameters parameters;
//        parameters.stab = stab;
//        parameters.configureTraditional(DGMaxTimeIntegration::CO2, 1, p.getValue(), numElements.getValue());
//        parameters.snapshotStride = 100;
//        timeSolver.solve(testProblem, parameters);
//
//        timeSolver.writeTimeSnapshots("domokos.dat");
//        timeSolver.printErrors({DGMaxDiscretization::L2, DGMaxDiscretization::HCurl}, testProblem);


        time(&end);
        std::cout << "Spend " << (end - start) << "s" << std::endl;
    }
    catch (const char* message)
    {
        std::cout << message;
    }
    // No need to clean globalData/configData, these are cleaned by the baseAPI.
    return 0;
}
