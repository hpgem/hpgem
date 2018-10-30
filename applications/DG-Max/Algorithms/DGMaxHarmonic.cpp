/*
This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "DGMaxHarmonic.h"

#include <petscksp.h>

#include <complex.h>
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"


DGMaxHarmonic::DGMaxHarmonic(hpGemUIExtentions &base)
        : base_ (base)
{
    discretization.initializeBasisFunctions(*(base_.getMesh(0)));
}

void DGMaxHarmonic::solve(const HarmonicProblem &harmonicProblem, double stab)
{
    PetscErrorCode error;
    std::cout << "finding a time-harmonic solution" << std::endl;
    discretization.computeElementIntegrands(
            *(base_.getMesh(0)),
            false,
            std::bind(&HarmonicProblem::sourceTerm, std::ref(harmonicProblem), std::placeholders::_1, std::placeholders::_2),
            nullptr, nullptr // No initial values
    );
    discretization.computeFaceIntegrals(
            *(base_.getMesh(0)),
            std::bind(&HarmonicProblem::boundaryCondition, std::ref(harmonicProblem), std::placeholders::_1,
                    std::placeholders::_2, std::placeholders::_3),
            stab
    );

    Utilities::GlobalPetscMatrix massMatrix(base_.getMesh(0), DGMaxDiscretization::MASS_MATRIX_ID, -1),
            stiffnessMatrix(base_.getMesh(0), DGMaxDiscretization::STIFFNESS_MATRIX_ID, DGMaxDiscretization::FACE_MATRIX_ID);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector
            resultVector(base_.getMesh(0), -1, -1), // The vector that we will use for the solution, initialize it with zeros.
            rhsVector(base_.getMesh(0), DGMaxDiscretization::SOURCE_TERM_VECTOR_ID, DGMaxDiscretization::FACE_VECTOR_ID);
    std::cout << "GlobalPetscVector initialised" << std::endl;
    resultVector.assemble();
    std::cout << "resultVector assembled" << std::endl;
    rhsVector.assemble();
    std::cout << "rhsVector assembled" << std::endl;

//    std::complex<double> I = std::complex<double>(0, 1);
    // Apply the factor i omega to the source term J.
//    error = VecScale(rhsVector, harmonicProblem.omega());
//    CHKERRABORT(PETSC_COMM_WORLD, error);

    KSP solver;
    error = KSPCreate(PETSC_COMM_WORLD, &solver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Preconditioner
    PC preconditioner;
    // The KSP solver will 'create' a valid preconditioner for us to work with.
    error = KSPGetPC(solver, &preconditioner);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = PCSetType(preconditioner, "jacobi");
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = KSPSetPC(solver, preconditioner);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    double omega2 = harmonicProblem.omega();
    omega2 *= omega2;

    error = MatAXPY(stiffnessMatrix, -1.0 * omega2, massMatrix, SUBSET_NONZERO_PATTERN);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPSetTolerances(solver, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPSetType(solver, "minres");
    CHKERRABORT(PETSC_COMM_WORLD, error);

    //everything that is set in the code, but before this line is overridden by command-line options
    error = KSPSetFromOptions(solver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPSetOperators(solver, stiffnessMatrix, stiffnessMatrix);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = KSPSetUp(solver);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = KSPSolve(solver, rhsVector, resultVector);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPDestroy(&solver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    resultVector.writeTimeIntegrationVector(0); //DOUBTFUL
}

std::map<DGMaxDiscretization::NormType, double> DGMaxHarmonic::computeError(
        const std::set<DGMaxDiscretization::NormType>& norms,
        const DGMaxDiscretization::InputFunction& exactSolution,
        const DGMaxDiscretization::InputFunction& exactSolutionCurl) const
{
    //Note, this only works by grace of distributing the solution as timeIntegrationVector
    return discretization.computeError(
            *(base_.getMesh(0)),
            0, // The abused time vector
            exactSolution,
            exactSolutionCurl,
            norms
    );
}

std::map<DGMaxDiscretization::NormType, double> DGMaxHarmonic::computeError(
        const std::set<DGMaxDiscretization::NormType> &norms, const ExactHarmonicProblem &problem) const
{
    return computeError(
            norms,
            std::bind(&ExactHarmonicProblem::exactSolution, std::ref(problem), std::placeholders::_1, std::placeholders::_2),
            std::bind(&ExactHarmonicProblem::exactSolutionCurl, std::ref(problem), std::placeholders::_1, std::placeholders::_2)
    );
}

void DGMaxHarmonic::writeTec(std::string fileName) const
{
    std::ofstream fileWriter;
    fileWriter.open(fileName);
    Output::TecplotDiscontinuousSolutionWriter<DIM> tecplotWriter(
            fileWriter, "The electric field", "012", "E0,E1,E2,H0,H1,H2");

    logger.log(Log::WARN, "The magnetic field output is incorrectly computed.");

    std::string zoneName ("field"); // Dummy
    tecplotWriter.write(base_.getMesh(0), zoneName, false,
        [&](const Base::Element* element, const Geometry::PointReference<DIM>& point, std::ostream& stream) {

        const LinearAlgebra::MiddleSizeVector coefficients = element->getTimeIntegrationVector(0);
        LinearAlgebra::SmallVector<DIM> electricField = discretization.computeField(element, point, coefficients);
        // TODO: Note that we computeCurlField already converts to real numbers, so we can not compute
        // H = i/(omega mu) curl E
        LinearAlgebra::SmallVector<DIM> curlField = discretization.computeCurlField(element, point, coefficients);
        stream
                << electricField[0] << " "
                << electricField[1] << " "
                << electricField[2] << " "
                << curlField[0] << " "
                << curlField[1] << " "
                << curlField[2] << std::endl;
    });

    fileWriter.close();
}