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

#include <cmath>

#include "Base/HpgemAPILinearSteadyState.h"
#include "petscksp.h"
#include "Utilities/BasisFunctions1DH1ConformingLine.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include "Logger.h"
#include <CMakeDefinitions.h>

//If this test ever breaks it is not a bad thing per se. However, once this breaks a thorough convergence analysis needs to be done.
//If the results still show the theoretically optimal order of convergence, and you are convinced that your changes improved the code,
//you should change the numbers in this test to reflect the updated result. Always confer with other developers if you do this.

/// \brief Class for solving the Poisson problem using HpgemAPILinearSteadyState and conforming basis functions.
template<std::size_t DIM>
class PoissonTest : public Base::HpgemAPILinearSteadyState<DIM>
{
public:
    using typename Base::HpgemAPIBase<DIM>::PointPhysicalT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceOnFaceT;

    PoissonTest(const std::string name, const std::size_t p) :
    Base::HpgemAPILinearSteadyState<DIM>(1, p, true, true),
    p_(p),
    totalError_(0)
    {
        using namespace std::string_literals;
        readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s  + name);
    }
    
    void readMesh(const std::string meshName) override final
    {
        
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numberOfElementMatrices = 2;   // Mass matrix and stiffness matrix
        std::size_t numberOfElementVectors = 1;    // Source term vector
        std::size_t numberOfFaceMatrices = 1;      // Stiffness matrix
        std::size_t numberOfFaceVectors = 1;       // Source term vector at boundary
        
        // Create mesh and set basis functions.
        this->addMesh(meshName, numberOfElementMatrices, numberOfElementVectors, numberOfFaceMatrices, numberOfFaceVectors);
        this->meshes_[0]->useDefaultConformingBasisFunctions();
        
        // Set the number of time integration vectors according to the size of the Butcher tableau.
        this->setNumberOfTimeIntegrationVectorsGlobally(this->globalNumberOfTimeIntegrationVectors_);
        
        // Plot info about the mesh
        std::size_t nElements = this->meshes_[0]->getNumberOfElements();
        logger(VERBOSE, "Total number of elements: %", nElements);
    }
    
    ///\brief Compute the integrand for the stiffness matrix at the element.
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(Base::PhysicalElement<DIM> &element) override final
    {
        //Obtain the number of basisfunctions that are possibly non-zero on this element.
        const std::size_t numberOfBasisFunctions = element.getElement()->getNrOfBasisFunctions();
        
        //Create the integrandVal such that it contains as many rows and columns as
        //the number of basisfunctions.
        LinearAlgebra::MiddleSizeMatrix& integrandVal = element.getResultMatrix();
        
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i)
        {
            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j)
            {
                //Compute the value of gradient(phi_i).gradient(phi_j) at point p and
                //store it at the appropriate place in the matrix integrandVal.
                integrandVal(j, i) = element.basisFunctionDeriv(i) * element.basisFunctionDeriv(j);
            }
        }
        
        return integrandVal;
    }

    //the default hpGEM solver expects to have to construct a face matrix, just give it the default one
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(Base::PhysicalFace<DIM> &face) override final
    {
        return face.getResultMatrix();
    }

    //the default hpGEM solver expects to have to construct a face vector, just give it the default one
    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(Base::PhysicalFace<DIM> &face) override final
    {
        return face.getResultVector();
    }

    /// \brief Define the exact solution
    LinearAlgebra::MiddleSizeVector getExactSolution(const PointPhysicalT &p) override final
    {
        LinearAlgebra::MiddleSizeVector exactSolution(1);
        
        double ret = std::sin(2 * M_PI * p[0]);
        if (p.size() > 1)
        {
            ret *= std::cos(2 * M_PI * p[1]) / 2.;
        }
        if (p.size() > 2)
        {
            ret *= std::cos(2 * M_PI * p[2]) * 2.;
        }
        
        exactSolution[0] = ret;
        return exactSolution;
    }
    
    ///\brief Define the source term.
    LinearAlgebra::MiddleSizeVector getSourceTerm(const PointPhysicalT &p) override final
    {
        LinearAlgebra::MiddleSizeVector sourceTerm(1);
        
        double ret = -std::sin(2 * M_PI * p[0]) * (4 * M_PI * M_PI);
        if (DIM > 1)
        {
            ret *= std::cos(2 * M_PI * p[1]);
        }
        if (DIM > 2)
        {
            ret *= std::cos(2 * M_PI * p[2]) * 3;
        }
        
        sourceTerm[0] = ret;
        return sourceTerm;
    }
    
    //This routine alters the matrix such that it can deal with conforming boundaries. It assumes correct boundary values are provided in
    //its third argument (the rest of the vector can be garbage) and that the second vector will be used as the RHS of a linear system solve
    //it clears the rows corresponding to the boudary nodes to have only a 1 on the diagonal and sets the RHS to the appropriate value
    void insertDirichletBoundary(Utilities::GlobalPetscMatrix& A, Utilities::GlobalPetscVector& b, Utilities::GlobalPetscVector& x)
    {
        std::size_t numberOfRows(0);
        std::vector<int> rows(0);
        Geometry::PointPhysical<DIM> pPhys;
        for (Base::Face* face : this->meshes_[0]->getFacesList())
        {
            const PointReferenceOnFaceT& center = face->getReferenceGeometry()->getCenter();
            pPhys = face->referenceToPhysical(center);
            //if the face is on a dirichlet boundary
            //if(face->faceType_=(...))
            if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1) < 1e-9)
            {
                //fetch the row numbers
                A.getMatrixBCEntries(face, numberOfRows, rows);
            }
        }
        int ierr = MatZeroRows(A, numberOfRows, &rows[0], 1.0, x, b);
        CHKERRV(ierr);
    }
    
    void solveSteadyStateWithPetsc(bool doComputeError) override final
    {
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
        // Create and Store things before solving the problem.
        this->tasksBeforeSolving();
        
        // Solve the linear problem
        //Assemble the matrix A of the system Ax = b.
        //The special value -1 is used to indicate there is no face matrix (since there is no flux in the conforming case)
        Utilities::GlobalPetscMatrix A(this->meshes_[0], this->stiffnessElementMatrixID_, -1);
        MatScale(A,-1);
        //Declare the vectors x and b of the system Ax = b.
        //The special value -1 is used to indicate there is no face matrix (since there is no flux in the conforming case)
        Utilities::GlobalPetscVector b(this->meshes_[0], this->sourceElementVectorID_, -1), x(this->meshes_[0]);
        
        //Assemble the vector b. This is needed because Petsc assumes you don't know
        //yet whether a vector is a variable or right-hand side the moment it is
        //declared.
        b.assemble();
        VecSet(x, 0);
        insertDirichletBoundary(A, b, x);
        
        //Make the Krylov supspace method
        KSP ksp;
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
            LinearAlgebra::MiddleSizeVector::type totalError = this->computeTotalError(this->solutionVectorId_, 0);
            totalError_ = totalError;
            logger(INFO, "Total error: %.", totalError);
            LinearAlgebra::MiddleSizeVector maxError = this->computeMaxError(this->solutionVectorId_, 0);
            logger.assert_debug(maxError.size() == this->configData_->numberOfUnknowns_, "Size of maxError (%) not equal to the number of variables (%)",
                                maxError.size(), this->configData_->numberOfUnknowns_);
            for(std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; iV ++)
            {
                logger(INFO, "Maximum error %: %", this->variableNames_[iV], maxError(iV));
            }
        }
        
        return;
#endif
    }
    
    LinearAlgebra::MiddleSizeVector::type getTotalError()
    {
        return totalError_;
    }
    
private:
    
    ///number of elements per cardinal direction
    int n_;
    
    ///polynomial order of the approximation
    int p_;
    
    /// Weighted L2 norm of the error
    LinearAlgebra::MiddleSizeVector::type totalError_;
};

 
int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    
    PoissonTest<1> test0("poissonMesh1.hpgem", 2);
    test0.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test0.getTotalError() - 0.48478776) < 1e-8), "comparison to old results");
    
    PoissonTest<1> test1("poissonMesh2.hpgem", 3);
    test1.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test1.getTotalError() - 0.02225892) < 1e-8), "comparison to old results");
    
    PoissonTest<1> test2("poissonMesh3.hpgem", 4);
    test2.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test2.getTotalError() - 0.00008248) < 1e-8), "comparison to old results");
    
    PoissonTest<1> test3("poissonMesh4.hpgem", 5);
    test3.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test3.getTotalError() - 0.00000008) < 1e-8), "comparison to old results");
    
    PoissonTest<1> test4("poissonMesh5.hpgem", 1);
    test4.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test4.getTotalError() - 0.00904309) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test5("poissonMesh6.hpgem", 2);
    test5.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test5.getTotalError() - 0.21870166) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test6("poissonMesh7.hpgem", 3);
    test6.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test6.getTotalError() - 0.02345377) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test7("poissonMesh8.hpgem", 4);
    test7.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test7.getTotalError() - 0.00039351) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test8("poissonMesh9.hpgem", 5);
    test8.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test8.getTotalError() - 0.00000066) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test9("poissonMesh10.hpgem", 1);
    test9.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test9.getTotalError() - 0.00911139) < 1e-8), "comparison to old results");
    
    PoissonTest<3> test10("poissonMesh11.hpgem", 5);
    test10.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test10.getTotalError() - 0.0121718) < 1e-8), "comparison to old results");

    PoissonTest<3> test11("poissonMesh12.hpgem", 4);
    test11.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test11.getTotalError() - 0.00216138) < 1e-8), "comparison to old results");

    PoissonTest<3> test12("poissonMesh13.hpgem", 3);
    test12.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test12.getTotalError() - 0.000906667) < 1e-8), "comparison to old results");

    PoissonTest<3> test13("poissonMesh14.hpgem", 2);
    test13.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test13.getTotalError() - 0.00138552) < 1e-8), "comparison to old results");

    PoissonTest<3> test14("poissonMesh15.hpgem", 1);
    test14.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test14.getTotalError() - 0.00454794) < 1e-8), "comparison to old results");
     
    return 0;
}

