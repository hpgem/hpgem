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

//#define hpGEM_INCLUDE_PETSC_SUPPORT//temporarily activating this definition makes development easier on some IDEs

#include <cmath>

#include "Base/HpgemAPILinearSteadyState.h"
#include "petscksp.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include "Logger.h"

//If this test ever breaks it is not a bad thing per se. However, once this breaks a thorough convergence analysis needs to be done.
//If the results still show the theoretically optimal order of convergence, and you are convinced that your changes improved the code,
//you should change the numbers in this test to reflect the updated result. Always confer with other developers if you do this.

/// \brief Class for solving the Poisson problem using HpgemAPILinearSteadyState.
class PoissonTest : public Base::HpgemAPILinearSteadyState
{
public:
    PoissonTest(const std::size_t n, const std::size_t p, const std::size_t dimension, const Base::MeshType meshType) :
    HpgemAPILinearSteadyState(dimension, 1, p, true, true),
    n_(n),
    p_(p),
    DIM_(dimension),
    totalError_(0)
    {
        penalty_ = 3 * n_ * p_ * (p_ + DIM_ - 1) + 1;
        createMesh(n_, meshType);
    }
    
    ///\brief set up the mesh
    Base::RectangularMeshDescriptor createMeshDescription(const std::size_t numOfElementPerDirection) override final
    {
        //describes a rectangular domain
        Base::RectangularMeshDescriptor description(DIM_);
        
        for (std::size_t i = 0; i < DIM_; ++i)
        {
            //define the value of the bottom left corner in each dimension
            description.bottomLeft_[i] = 0;
            //define the value of the top right corner in each dimension
            description.topRight_[i] = 1;
            //define how many elements there should be in the direction of dimension
            //At this stage, the mesh first consists of n^2 squares, and later these
            //squares can be divided in two triangles each if a triangular mesh is desired.
            description.numElementsInDIM_[i] = n_;
            //define whether you have periodic boundary conditions or a solid wall in this direction.
            description.boundaryConditions_[i] = Base::BoundaryType::SOLID_WALL;
        }
        
        return description;
    }
    
    ///\brief Compute the integrand for the stiffness matrix at the element.
    LinearAlgebra::Matrix computeIntegrandStiffnessMatrixAtElement(const Base::Element *element, const PointReferenceT &point) override final
    {
        //Obtain the number of basisfunctions that are possibly non-zero on this element.
        const std::size_t numBasisFunctions = element->getNrOfBasisFunctions();
        
        //Create the integrandVal such that it contains as many rows and columns as
        //the number of basisfunctions.
        LinearAlgebra::Matrix integrandVal(numBasisFunctions, numBasisFunctions);
        
        for (std::size_t i = 0; i < numBasisFunctions; ++i)
        {
            for (std::size_t j = 0; j < numBasisFunctions; ++j)
            {
                //Compute the value of gradient(phi_i).gradient(phi_j) at point p and
                //store it at the appropriate place in the matrix integrandVal.
                integrandVal(j, i) = element->basisFunctionDeriv(i, point) * element->basisFunctionDeriv(j, point);
            }
        }
        
        return integrandVal;
    }
    
    /// \brief Compute the integrand for the siffness matrix at the face.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(const Base::Face* face, const LinearAlgebra::MiddleSizeVector& normal, const PointReferenceT& p) override final
    {
        //Get the number of basis functions, first of both sides of the face and
        //then only the basis functions associated with the left and right element.
        std::size_t numBasisFunctions = face->getNrOfBasisFunctions();
        std::size_t nLeft = face->getPtrElementLeft()->getNrOfBasisFunctions();
        std::size_t nRight = 0;
        if(face->isInternal())
        {
            nRight = face->getPtrElementLeft()->getNrOfBasisFunctions();
        }
        
        //Create the FaceMatrix integrandVal with the correct size.
        Base::FaceMatrix integrandVal(nLeft, nRight);
        
        //Initialize the vectors that contain gradient(phi_i), gradient(phi_j), normal_i phi_i and normal_j phi_j
        LinearAlgebra::MiddleSizeVector phiNormalI(DIM_), phiNormalJ(DIM_), phiDerivI(DIM_), phiDerivJ(DIM_);
        
        //Transform the point from the reference value to its physical value.
        //This is necessary to check at which boundary we are if we are at a boundary face.
        PointPhysicalT pPhys = face->referenceToPhysical(p);
        
        for (int i = 0; i < numBasisFunctions; ++i)
        {
            //normal_i phi_i is computed at point p, the result is stored in phiNormalI.
            phiNormalI = face->basisFunctionNormal(i, normal, p);
            //The gradient of basisfunction phi_i is computed at point p, the result is stored in phiDerivI.
            phiDerivI = face->basisFunctionDeriv(i, p);
            
            for (int j = 0; j < numBasisFunctions; ++j)
            {
                //normal_j phi_j is computed at point p, the result is stored in phiNormalJ.
                phiNormalJ = face->basisFunctionNormal(j, normal, p);
                //The gradient of basisfunction phi_j is computed at point p, the result is stored in phiDerivJ.
                phiDerivJ = face->basisFunctionDeriv(j, p);
                
                //Switch to the correct type of face, and compute the integrand accordingly
                //you could also compute the integrandVal by directly using face->basisFunctionDeriv
                //and face->basisFunctionNormal in the following lines, but this results in very long expressions
                //Internal face:
                if (face->isInternal())
                {
                    integrandVal(j, i) = -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) / 2 + penalty_ * phiNormalI * phiNormalJ;
                }
                //Boundary face with Dirichlet boundary conditions:
                else if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1.) < 1e-9)
                {
                    integrandVal(j, i) = -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) + penalty_ * phiNormalI * phiNormalJ * 2;
                }
                //Boundary face with homogeneous Neumann boundary conditions:
                else
                {
                    integrandVal(j, i) = 0;
                }
            }
        }
        
        return integrandVal;
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
        if (DIM_ > 1)
        {
            ret *= std::cos(2 * M_PI * p[1]);
        }
        if (DIM_ > 2)
        {
            ret *= std::cos(2 * M_PI * p[2]) * 3;
        }
        
        sourceTerm[0] = ret;
        return sourceTerm;
    }
    
    /// \brief Compute the integrals of the right-hand side associated with faces.
    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(const Base::Face* face, const LinearAlgebra::MiddleSizeVector& normal, const PointReferenceT& p) override final
    {
        //Obtain the number of basisfunctions that are possibly non-zero
        const std::size_t numBasisFunctions = face->getNrOfBasisFunctions();
        //Resize the integrandVal such that it contains as many rows as
        //the number of basisfunctions.
        LinearAlgebra::MiddleSizeVector integrandVal(numBasisFunctions);
        
        //Compute the value of the integrand
        //We have no rhs face integrals, so this is just 0.
        for (std::size_t i = 0; i < numBasisFunctions; ++i)
        {
            integrandVal[i] = 0;
        }
        
        return integrandVal;
    }
    
    void solveSteadyStateWithPetsc(bool doComputeError) override final
    {
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
        // Create and Store things before solving the problem.
        tasksBeforeSolving();
        
        // Solve the linear problem
        //Assemble the matrix A of the system Ax = b.
        Utilities::GlobalPetscMatrix A(HpgemAPIBase::meshes_[0], stiffnessElementMatrixID_, stiffnessFaceMatrixID_);
        MatScale(A,-1);
        //Declare the vectors x and b of the system Ax = b.
        Utilities::GlobalPetscVector b(HpgemAPIBase::meshes_[0], sourceElementVectorID_, sourceFaceVectorID_), x(HpgemAPIBase::meshes_[0]);
        
        //Assemble the vector b. This is needed because Petsc assumes you don't know
        //yet whether a vector is a variable or right-hand side the moment it is
        //declared.
        b.assemble();
        
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
        
        x.writeTimeLevelData(solutionTimeLevel_);
        
        if(doComputeError)
        {
            double totalError = computeTotalError(solutionTimeLevel_, 0);
            totalError_ = totalError;
            logger(INFO, "Total error: %.", totalError);
            LinearAlgebra::MiddleSizeVector maxError = computeMaxError(solutionTimeLevel_, 0);
            logger.assert(maxError.size() == configData_->numberOfUnknowns_, "Size of maxError (%) not equal to the number of variables (%)", maxError.size(), configData_->numberOfUnknowns_);
            for(std::size_t iV = 0; iV < configData_->numberOfUnknowns_; iV ++)
            {
                logger(INFO, "Maximum error %: %", variableNames_[iV], maxError(iV));
            }
        }
        
        return;
#endif
    }
    
    double getTotalError()
    {
        return totalError_;
    }
        
private:
    
    ///number of elements per cardinal direction
    int n_;
    
    ///polynomial order of the approximation
    int p_;
    
    ///Dimension of the domain, in this case 2
    int DIM_;
    
    ///\brief Penalty parameter
    ///
    ///Penalty parameter that is associated with the interior penalty discontinuous
    ///Galerkin method. This parameter is initialized in the constructor, and has
    ///to be greater than 3 * n_ * p_ * (p_ + DIM - 1) in order for the method to be stable.
    double penalty_;
    
    /// Weighted L2 norm of the error
    double totalError_;
};


int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    //no 3D testing due to speed related issues
    
    PoissonTest test0(1, 2, 1, Base::MeshType::RECTANGULAR);
    test0.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test0.getTotalError() - 0.35188045) < 1e-8), "comparison to old results");
    
    PoissonTest test1(2, 3, 1, Base::MeshType::RECTANGULAR);
    test1.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test1.getTotalError() - 0.01607749) < 1e-8), "comparison to old results");
    
    PoissonTest test2(4, 4, 1, Base::MeshType::RECTANGULAR);
    test2.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test2.getTotalError() - 0.00007200) < 1e-8), "comparison to old results");
    
    PoissonTest test3(8, 5, 1, Base::MeshType::RECTANGULAR);
    test3.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test3.getTotalError() - 0.00000008) < 1e-8), "comparison to old results");
    
    PoissonTest test4(16, 1, 1, Base::MeshType::RECTANGULAR);
    test4.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test4.getTotalError() - 0.00880380) < 1e-8), "comparison to old results");
    
    PoissonTest test5(1, 2, 2, Base::MeshType::TRIANGULAR);
    test5.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test5.getTotalError() - 0.17226144) < 1e-8), "comparison to old results");
    
    PoissonTest test6(2, 3, 2, Base::MeshType::TRIANGULAR);
    test6.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test6.getTotalError() - 0.01782337) < 1e-8), "comparison to old results");
    
    PoissonTest test7(4, 4, 2, Base::MeshType::TRIANGULAR);
    test7.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test7.getTotalError() - 0.00035302) < 1e-8), "comparison to old results");
    
    PoissonTest test8(8, 5, 2, Base::MeshType::TRIANGULAR);
    test8.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test8.getTotalError() - 0.00000061) < 1e-8), "comparison to old results");
    
    PoissonTest test9(16, 1, 2, Base::MeshType::TRIANGULAR);
    test9.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test9.getTotalError() - 0.00448270) < 1e-8), "comparison to old results");
    
    return 0;
}
