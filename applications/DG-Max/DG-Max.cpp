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

#include "DGMaxLogger.h"

#include "Algorithms/DGMaxEigenValue.h"
#include "Algorithms/DGMaxHarmonic.h"
#include "Algorithms/DGMaxTimeIntegration.h"
#include "Algorithms/DivDGMaxEigenValue.h"
#include "Algorithms/DivDGMaxHarmonic.h"

#include "ProblemTypes/Harmonic/SampleHarmonicProblems.h"
#include "ProblemTypes/Time/SampleTestProblems.h"
#include "Utils/HomogeneousBandStructure.h"
#include "Utils/BandstructureGNUPlot.h"

/**
 * This class should provide problem specific information about the maxwell equations.
 */
template<std::size_t DIM>
class DGMax : public hpGemUIExtentions<DIM>
{
public:
    
    DGMax(Base::GlobalData * const globalConfig, Base::ConfigurationData* elementConfig)
            : hpGemUIExtentions<DIM>(globalConfig, elementConfig)
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

        auto mesh = new Base::MeshManipulator<DIM>(this->getConfigData(), Base::BoundaryType::PERIODIC,
                                                   Base::BoundaryType::PERIODIC, Base::BoundaryType::PERIODIC,
                                                   2, 3, 1, 1);
        mesh->createTriangularMesh(bottomLeft, topRight, numElementsOneD);

        for (typename Base::MeshManipulator<DIM>::ElementIterator it = mesh->elementColBegin(Base::IteratorType::GLOBAL);
                it != mesh->elementColEnd(Base::IteratorType::GLOBAL); ++it)
        {
            (*it)->setUserData(new ElementInfos(**it));
        }
        this->addMesh(mesh);
    }
    
    void readMesh(std::string fileName)
    {
        auto mesh = new Base::MeshManipulator<DIM>(this->getConfigData(), Base::BoundaryType::PERIODIC,
                                                   Base::BoundaryType::PERIODIC, Base::BoundaryType::PERIODIC,
                                                   2, 3, 1, 1);
        mesh->readMesh(fileName);
        this->addMesh(mesh);
        for (typename Base::MeshManipulator<DIM>::ElementIterator it = mesh->elementColBegin(Base::IteratorType::GLOBAL);
                it != mesh->elementColEnd(Base::IteratorType::GLOBAL); ++it)
        {
            (*it)->setUserData(new ElementInfos(**it));
        }
        
    }
};

auto& numElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", false, 0);
auto& p = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);
auto& meshFile = Base::register_argument<std::string>('m', "meshFile", "The hpgem meshfile to use", false);
auto& numEigenvalues = Base::register_argument<std::size_t>('e', "eigenvalues", "The number of eigenvalues to compute", false, 24);

//Temporary
const std::size_t DIM = 2;

void printArguments(int argc, char** argv)
{
#ifdef HPGEM_USE_MPI
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    logAll([&](){DGMaxLogger(INFO, "Proc %/%", rank, size);});
#endif
    if (!loggingSuppressed())
    {
        std::stringstream stream;
        stream << "Program arguments: " << std::endl;
        for(int i = 0; i < argc; ++i)
        {
            if (i != 0) stream << " ";
            stream << argv[i];
        }
        std::string message = stream.str();
        DGMaxLogger(INFO, message);
        DGMaxLogger(INFO, "Using DGMax for % dimensions.", DIM);
    }
}

/**
 * main routine. Suggested usage: make a problem, solve it, tell the user your results.
 */
int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    initDGMaxLogging();
    printArguments(argc, argv);


    DGMaxLogger.assert_always( DIM >= 2 && DIM <= 3, "Can only handle 2D and 3D problems.");
    DGMaxLogger.assert_always(numElements.isUsed() || meshFile.isUsed(),
            "DGMax requires either a mesh file or a number of elements");

    //set up timings
    time_t start, end, initialised, solved;
    time(&start);
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
    Base::ConfigurationData* const configData = new Base::ConfigurationData(numberOfUnknowns, numberOfTimeLevels);
    try
    {
        double stab = (p.getValue() + 1) * (p.getValue() + 3);
        DivDGMaxDiscretization<DIM>::Stab divStab;
        // Values from the Jelmer fix code.
        divStab.stab1 = 100;
        // Note: for 2D harmonic it looks like that we need 10 instead of 0.01.
        divStab.stab2 = 0.01;
        divStab.stab3 = 10.0;

        DGMax<DIM> base(globalData, configData);

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

//        DGMaxHarmonic harmonicSolver (base, p.getValue());
//        DivDGMaxHarmonic harmonicSolver (base);
//
//        SampleHarmonicProblems problem (SampleHarmonicProblems::SARMANY2010, 1);
//        harmonicSolver.solve(problem, divStab, p.getValue());
//        std::cout << "L2 error " << harmonicSolver.computeL2Error(problem) << std::endl;
//        auto errors = harmonicSolver.computeError({DGMaxDiscretization::L2, DGMaxDiscretization::HCurl}, problem);
//        std::cout << "L2 error    " << errors[DGMaxDiscretization::L2] << std::endl;
//        std::cout << "HCurl error " << errors[DGMaxDiscretization::HCurl] << std::endl;


        /////////////////////
        // Eigenvalue code //
        /////////////////////

//        DGMaxEigenValue solver (base, p.getValue());
        DivDGMaxEigenValue<DIM> solver (base);
        KSpacePath<DIM> path = KSpacePath<DIM>::cubePath(20);
        EigenValueProblem<DIM> input(path, numEigenvalues.getValue());
        DivDGMaxEigenValue<DIM>::Result result = solver.solve(input, divStab, p.getValue());
        if (Base::MPIContainer::Instance().getProcessorID() == 0)
        {
            result.printFrequencies();
        }

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


        // Quick way of plotting all the eigenvalue results against theory
        // Probably should be split off at some point to output the bandstructure
        // and do this plotting separately.
        std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors;
        for(std::size_t i = 0; i < DIM; ++i) {
            reciprocalVectors[i].set(0);
            reciprocalVectors[i][i] = 2.0*M_PI;
        }

        HomogeneousBandStructure<DIM> structure (reciprocalVectors);
        std::vector<std::string> points ({"G", "X"});
        if(DIM == 2)
        {
            points.emplace_back("M");
        }
        else
        {
            points.emplace_back("S");
            points.emplace_back("R");
        }
        BandstructureGNUPlot<DIM> output (path, points, structure, &result);
        output.plot("gnutest.in");

        time(&end);
        DGMaxLogger(INFO, "DGMax finished in %s", end - start);
    }
    catch (const char* message)
    {
        DGMaxLogger(ERROR, message);
    }
    // No need to clean globalData/configData, these are cleaned by the baseAPI.
    return 0;
}
