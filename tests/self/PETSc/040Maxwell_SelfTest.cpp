/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2015, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#include <cmath>

#include "Base/HpgemAPILinearSteadyState.h"
#include "petscksp.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "Base/HCurlConformingTransformation.h"

#include "Logger.h"
#include <CMakeDefinitions.h>


//If this test ever breaks it is not a bad thing per se. However, once this breaks a thorough convergence analysis needs to be done.
//If the results still show the theoretically optimal order of convergence, and you are convinced that your changes improved the code,
//you should change the numbers in this test to reflect the updated result. Always confer with other developers if you do this.
//Note that the nedelec basis functions only converge with order p for the l2norm. They DO converge optimally for the Hcurl norm.

/// \brief Class for solving the Poisson problem using HpgemAPILinearSteadyState.
class MaxwellTest : public Base::HpgemAPILinearSteadyState<3>
{
public:

    MaxwellTest(const std::string name, const std::size_t n, const std::size_t p, bool useNedelec) :
    Base::HpgemAPILinearSteadyState<3>(1, p, true, true),
    totalError(0)
    {
        using namespace std::string_literals;
        penalty = n * (p + 1) * (p + 3) + 1;
        readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s  + name);
        elementIntegrator_.setTransformation(std::shared_ptr<Base::CoordinateTransformation<3>>(new Base::HCurlConformingTransformation<3>));
        faceIntegrator_.setTransformation(std::shared_ptr<Base::CoordinateTransformation<3>>(new Base::HCurlConformingTransformation<3>));
        if(useNedelec)
        {
            meshes_[0]->useNedelecDGBasisFunctions(p);
        }
        else
        {
            meshes_[0]->useAinsworthCoyleDGBasisFunctions(p);
        }
    }

    ///\brief Compute the integrand for the stiffness matrix at the element.
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(Base::PhysicalElement<3> &element) override final
    {
        //Obtain the number of basisfunctions that are possibly non-zero on this element.
        const std::size_t numberOfBasisFunctions = element.getElement()->getNumberOfBasisFunctions();

        LinearAlgebra::SmallVector<3> phi_i, phi_j;

        //Create the integrandVal such that it contains as many rows and columns as
        //the number of basisfunctions.
        LinearAlgebra::MiddleSizeMatrix& integrandVal = element.getResultMatrix();

        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i)
        {
            element.basisFunction(i, phi_i);
            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j)
            {
                element.basisFunction(j, phi_j);
                //Compute the value of curl(phi_i).curl(phi_j) - phi_i.phi_j at point p and
                //store it at the appropriate place in the matrix integrandVal.
                integrandVal(i, j) = element.basisFunctionCurl(i) * element.basisFunctionCurl(j) - phi_i * phi_j;
            }
        }
        return integrandVal;
    }

    /// \brief Compute the integrand for the siffness matrix at the face.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(Base::PhysicalFace<3> &face) override final
    {
        //Get the number of basis functions, first of both sides of the face and
        //then only the basis functions associated with the left and right element.
        std::size_t numberOfBasisFunctions = face.getFace()->getNumberOfBasisFunctions();

        //Create the FaceMatrix integrandVal with the correct size.
        Base::FaceMatrix& integrandVal = face.getResultMatrix();

        //Initialize the vectors that contain gradient(phi_i), gradient(phi_j), normal_i phi_i and normal_j phi_j
        LinearAlgebra::SmallVector<3> phiNormalI, phiNormalJ, phiCurlI, phiCurlJ;

        //Transform the point from the reference value to its physical value.
        //This is necessary to check at which boundary we are if we are at a boundary face.
        const PointPhysicalT& pPhys = face.getPointPhysical();

        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i)
        {
            //normal_i phi_i is computed at point p, the result is stored in phiNormalI.
            face.basisFunctionUnitNormalCross(i, phiNormalI);
            //The gradient of basisfunction phi_i is computed at point p, the result is stored in phiDerivI.
            phiCurlI = face.basisFunctionCurl(i);

            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j)
            {
                //normal_j phi_j is computed at point p, the result is stored in phiNormalJ.
                face.basisFunctionUnitNormalCross(j, phiNormalJ);
                //The gradient of basisfunction phi_j is computed at point p, the result is stored in phiDerivJ.
                phiCurlJ = face.basisFunctionCurl(j);

                //Switch to the correct type of face, and compute the integrand accordingly
                //you could also compute the integrandVal by directly using face->basisFunctionDeriv
                //and face->basisFunctionNormal in the following lines, but this results in very long expressions
                //Internal face:
                if (face.isInternal())
                {
                    integrandVal(j, i) = -(phiNormalI * phiCurlJ + phiNormalJ * phiCurlI) / 2 + penalty * phiNormalI * phiNormalJ;
                }
                else
                {
                    integrandVal(j, i) = -(phiNormalI * phiCurlJ + phiNormalJ * phiCurlI) + penalty * phiNormalI * phiNormalJ;
                }
            }
        }

        return integrandVal;
    }

    /// \brief Define the exact solution
    std::vector<LinearAlgebra::SmallVector<3>> getExactSolutionVector(const PointPhysicalT &p)
    {
        std::vector<LinearAlgebra::SmallVector<3>> exactSolution(1);

        exactSolution[0][0]=std::sin(M_PI * p[1]) * std::sin(M_PI * p[2]);
        exactSolution[0][1]=std::sin(M_PI * p[2]) * std::sin(M_PI * p[0]);
        exactSolution[0][2]=std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);

        return exactSolution;
    }

    ///\brief Define the source term.
    std::vector<LinearAlgebra::SmallVector<3>> getSourceTermVector(const PointPhysicalT &p)
    {
        std::vector<LinearAlgebra::SmallVector<3>> sourceTerm = getExactSolutionVector(p);

        sourceTerm[0] *= -2 * M_PI * M_PI + 1;

        return sourceTerm;
    }

    std::vector<LinearAlgebra::SmallVector<3>> boundaryConditions(const PointPhysicalT &p)
    {
        return getExactSolutionVector(p);
    }

    /// \brief Compute the integrals of the right-hand side associated with faces.
    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(Base::PhysicalFace<3> &face) override final
    {
        if(face.getFace()->isInternal())
        {
            return face.getResultVector();
        }
        else
        {
            LinearAlgebra::MiddleSizeVector& result = face.getResultVector();
            PointPhysicalT pPhys = face.getPointPhysical();
            LinearAlgebra::SmallVector<3> value = boundaryConditions(pPhys)[0];
            LinearAlgebra::SmallVector<3> normalValue = LinearAlgebra::SmallMatrix<3, 2>{{face.getUnitNormalVector(), value}}.computeWedgeStuffVector();
            LinearAlgebra::SmallVector<3> phi;
            for(std::size_t i = 0; i < face.getNumberOfBasisFunctions(); ++i)
            {
                face.basisFunctionNormalCross(i, phi);
                result[i] = -face.basisFunctionCurl(i) * normalValue + penalty * phi * normalValue;
            }
            return result;
        }
    }

    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtElement(Base::PhysicalElement<3> &element, const double time, const std::size_t orderTimeDerivative) override final
    {
        // Get a reference to the result vector.
        LinearAlgebra::MiddleSizeVector &integrand = element.getResultVector();

        // Get the physical point.
        PointPhysicalT pPhys = element.getPointPhysical();

        // Compute the source term.
        std::vector<LinearAlgebra::SmallVector<3>> sourceTerm = getSourceTermVector(pPhys);

        // Get the number of basis functions.
        const std::size_t numberOfBasisFunctions = element.getElement()->getNumberOfBasisFunctions();
        LinearAlgebra::SmallVector<3> functionValue;

        // Compute the product of the source term and all test functions.
        std::size_t iVB; // indices for both variable and basis function.
        for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
        {
            iVB = element.getElement()->convertToSingleIndex(iB, 0);
            element.basisFunction(iB, functionValue);
            integrand(iVB) = functionValue * sourceTerm[0];
        }

        return integrand;
    }

    LinearAlgebra::MiddleSizeVector::type computeIntegrandTotalError(Base::PhysicalElement<3>& element, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time) override final
    {

        // Get the physical point.
        const PointPhysicalT &pPhys = element.getPointPhysical();

        // Compute the real solution.
        const std::vector<LinearAlgebra::SmallVector<3>> exactSolution = getExactSolutionVector(pPhys);
        LinearAlgebra::SmallVector<3> functionValue;

        // Get the number of basis functions.
        const std::size_t numberOfBasisFunctions = element.getElement()->getNumberOfBasisFunctions();

        // Compute the numerical solution.
        LinearAlgebra::SmallVector<3> numericalSolution;
        numericalSolution *= 0;
        for(std::size_t jB = 0; jB < numberOfBasisFunctions; jB++)
        {
            std::size_t jVB = element.getElement()->convertToSingleIndex(jB, 0);
            element.basisFunction(jB, functionValue);
            numericalSolution += functionValue * std::real(solutionCoefficients[jVB]);
        }

        // Compute the error
        const LinearAlgebra::SmallVector<3> error = exactSolution[0] - numericalSolution;

        // Compute the square of the l2 norm.
        return error * error;
    }

    LinearAlgebra::MiddleSizeVector::type getTotalError()
    {
        return this->totalError;
    }

    void solveSteadyStateWithPetsc(bool doComputeError) override final
    {
#if defined(HPGEM_USE_ANY_PETSC)
        // Create and Store things before solving the problem.
        this->tasksBeforeSolving();

        // Solve the linear problem
        Utilities::GlobalIndexing indexing (this->meshes_[0]);
        //Assemble the matrix A of the system Ax = b.
        Utilities::GlobalPetscMatrix A(indexing, this->stiffnessElementMatrixID_, this->stiffnessFaceMatrixID_);
        MatScale(A,-1);
        //Declare the vectors x and b of the system Ax = b.
        Utilities::GlobalPetscVector b(indexing, this->sourceElementVectorID_, this->sourceFaceVectorID_),
            x(indexing);

        //Assemble the vector b. This is needed because Petsc assumes you don't know
        //yet whether a vector is a variable or right-hand side the moment it is
        //declared.
        b.assemble();

        //Make the Krylov supspace method
        KSP ksp;
        PC preconditioner;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        //Tell ksp that it will solve the system Ax = b.
        KSPSetOperators(ksp, A, A);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, b, x);
        //Do PETSc magic, including solving.
        KSPConvergedReason converge;
        KSPGetConvergedReason(ksp, &converge);
        int iterations;
        KSPGetIterationNumber(ksp, &iterations);
        logger(INFO, "KSP solver ended because of % in % iterations.", KSPConvergedReasons[converge], iterations);

        x.writeTimeIntegrationVector(this->solutionVectorId_);

        if(doComputeError)
        {
            totalError = this->computeTotalError(this->solutionVectorId_, 0);
            logger(INFO, "Total error: %.", totalError);
        }

        return;
#endif
    }

private:

    LinearAlgebra::MiddleSizeVector::type totalError;
    double penalty;
};


int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    MaxwellTest test0("maxwellMesh1.hpgem", 1, 4, true);
    test0.solveSteadyStateWithPetsc(true);
    logger.assert_always(std::abs(test0.getTotalError() - 0.03355038) < 1e-8, "comparison to old results");

    MaxwellTest test1("maxwellMesh2.hpgem", 2, 3, true);
    test1.solveSteadyStateWithPetsc(true);
    logger.assert_always(std::abs(test1.getTotalError() - 0.01697762) < 1e-8, "comparison to old results");

    MaxwellTest test2("maxwellMesh3.hpgem", 4, 2, true);
    test2.solveSteadyStateWithPetsc(true);
    logger.assert_always(std::abs(test2.getTotalError() - 0.02863282) < 1e-8, "comparison to old results");

    MaxwellTest test3("maxwellMesh4.hpgem", 8, 1, true);
    test3.solveSteadyStateWithPetsc(true);
    logger.assert_always(std::abs(test3.getTotalError() - 0.1189594) < 1e-8, "comparison to old results");

    MaxwellTest test4("maxwellMesh1.hpgem", 1, 4, false);
    test4.solveSteadyStateWithPetsc(true);
    logger.assert_always(std::abs(test4.getTotalError() - 0.02266147) < 1e-8, "comparison to old results");

    MaxwellTest test5("maxwellMesh2.hpgem", 2, 3, false);
    test5.solveSteadyStateWithPetsc(true);
    logger.assert_always(std::abs(test5.getTotalError() - 0.00469313) < 1e-8, "comparison to old results");

    MaxwellTest test6("maxwellMesh3.hpgem", 4, 2, false);
    test6.solveSteadyStateWithPetsc(true);
    logger.assert_always(std::abs(test6.getTotalError() - 0.00306707) < 1e-8, "comparison to old results");

    MaxwellTest test7("maxwellMesh4.hpgem", 8, 1, false);
    test7.solveSteadyStateWithPetsc(true);
    logger.assert_always(std::abs(test7.getTotalError() - 0.00828368) < 1e-8, "comparison to old results");
}

