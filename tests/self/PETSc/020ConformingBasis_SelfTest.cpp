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
#include "Utilities/BasisFunctions1DH1ConformingLine.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include "Logger.h"

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

    PoissonTest(const std::size_t n, const std::size_t p, const std::size_t dimension, const Base::MeshType meshType) :
    Base::HpgemAPILinearSteadyState<DIM>(dimension, 1, p, true, true),
    n_(n),
    p_(p),
    DIM_(dimension),
    totalError_(0)
    {
        penalty_ = 3 * n_ * p_ * (p_ + DIM_ - 1) + 1;
        createMesh(n_, meshType);
    }
    
    ///\brief set up the mesh
    Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numOfElementPerDirection) override final
    {
        //describes a rectangular domain
        Base::RectangularMeshDescriptor<DIM> description;
        
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
    
    void createMesh(const std::size_t numOfElementsPerDirection, const Base::MeshType meshType) override final
    {
        const Base::RectangularMeshDescriptor<DIM> description = createMeshDescription(numOfElementsPerDirection);
        
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numOfElementMatrices = 2;   // Mass matrix and stiffness matrix
        std::size_t numOfElementVectors = 1;    // Source term vector
        std::size_t numOfFaceMatrices = 1;      // Stiffness matrix
        std::size_t numOfFaceVectors = 1;       // Source term vector at boundary
        
        // Create mesh and set basis functions.
        this->addMesh(description, meshType, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
        this->meshes_[0]->useDefaultConformingBasisFunctions();
        
        std::size_t nElements = this->meshes_[0]->getNumberOfElements();
        logger(VERBOSE, "Total number of elements: %", nElements);
    }
    
    ///\brief Compute the integrand for the stiffness matrix at the element.
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(const Base::Element *element, const PointReferenceT &point) override final
    {
        //Obtain the number of basisfunctions that are possibly non-zero on this element.
        const std::size_t numBasisFunctions = element->getNrOfBasisFunctions();
        
        //Create the integrandVal such that it contains as many rows and columns as
        //the number of basisfunctions.
        LinearAlgebra::MiddleSizeMatrix integrandVal(numBasisFunctions, numBasisFunctions);
        
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
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(const Base::Face* face, const LinearAlgebra::SmallVector<DIM>& normal, const PointReferenceOnFaceT& p) override final
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
        LinearAlgebra::SmallVector<DIM> phiNormalI, phiNormalJ, phiDerivI, phiDerivJ;
        
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
    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(const Base::Face* face, const LinearAlgebra::SmallVector<DIM>& normal, const PointReferenceOnFaceT& p) override final
    {
        //Obtain the number of basisfunctions that are possibly non-zero
        const std::size_t numBasisFunctions = face->getNrOfBasisFunctions();
        //Resize the integrandVal such that it contains as many rows as
        //the number of basisfunctions.
        LinearAlgebra::MiddleSizeVector integrandVal(numBasisFunctions);
        
        PointPhysicalT pPhys = face->referenceToPhysical(p);
        if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1) < 1e-9)
        { //Dirichlet
            LinearAlgebra::SmallVector<DIM> phiDeriv;
            for (std::size_t i = 0; i < numBasisFunctions; ++i)
            {
                phiDeriv = face->basisFunctionDeriv(i, p);
                integrandVal[i] = (-normal * phiDeriv / Base::L2Norm(normal) + penalty_ * face->basisFunction(i, p)) * 0;
            }
        }
        else
        {
            for (std::size_t i = 0; i < numBasisFunctions; ++i)
            {
                integrandVal[i] = 0;
            }
        }
        
        return integrandVal;
    }
    
    void insertDirichletBoundary(Utilities::GlobalPetscMatrix& A, Utilities::GlobalPetscVector& b, Utilities::GlobalPetscVector& x)
    {
        std::size_t numberOfRows(0);
        std::vector<int> rows(0);
        Geometry::PointPhysical<DIM> pPhys;
        for (Base::Face* face : this->meshes_[0]->getFacesList())
        {
            const PointReferenceOnFaceT& center = face->getReferenceGeometry()->getCenter();
            pPhys = face->referenceToPhysical(center);
            //if(face->faceType_=(...))
            if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1) < 1e-9)
            {
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
        Utilities::GlobalPetscMatrix A(this->meshes_[0], this->stiffnessElementMatrixID_, this->stiffnessFaceMatrixID_);
        MatScale(A,-1);
        //Declare the vectors x and b of the system Ax = b.
        Utilities::GlobalPetscVector b(this->meshes_[0], this->sourceElementVectorID_, this->sourceFaceVectorID_), x(this->meshes_[0]);
        
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
        
        x.writeTimeLevelData(this->solutionTimeLevel_);
        
        if(doComputeError)
        {
            LinearAlgebra::MiddleSizeVector::type totalError = this->computeTotalError(this->solutionTimeLevel_, 0);
            totalError_ = totalError;
            logger(INFO, "Total error: %.", totalError);
            LinearAlgebra::MiddleSizeVector maxError = this->computeMaxError(this->solutionTimeLevel_, 0);
            logger.assert(maxError.size() == this->configData_->numberOfUnknowns_, "Size of maxError (%) not equal to the number of variables (%)", maxError.size(), this->configData_->numberOfUnknowns_);
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
    
    ///Dimension of the domain, in this case 1 or 2
    int DIM_;
    
    ///\brief Penalty parameter
    ///
    ///Penalty parameter that is associated with the interior penalty discontinuous
    ///Galerkin method. This parameter is initialized in the constructor, and has
    ///to be greater than 3 * n_ * p_ * (p_ + DIM - 1) in order for the method to be stable.
    double penalty_;
    
    /// Weighted L2 norm of the error
    LinearAlgebra::MiddleSizeVector::type totalError_;
};

 
int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    
    PoissonTest<1> test0(1, 2, 1, Base::MeshType::RECTANGULAR);
    test0.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test0.getTotalError() - 0.48478776) < 1e-8), "comparison to old results");
    
    PoissonTest<1> test1(2, 3, 1, Base::MeshType::RECTANGULAR);
    test1.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test1.getTotalError() - 0.02225892) < 1e-8), "comparison to old results");
    
    PoissonTest<1> test2(4, 4, 1, Base::MeshType::RECTANGULAR);
    test2.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test2.getTotalError() - 0.00008248) < 1e-8), "comparison to old results");
    
    PoissonTest<1> test3(8, 5, 1, Base::MeshType::RECTANGULAR);
    test3.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test3.getTotalError() - 0.00000008) < 1e-8), "comparison to old results");
    
    PoissonTest<1> test4(16, 1, 1, Base::MeshType::RECTANGULAR);
    test4.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test4.getTotalError() - 0.00904309) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test5(1, 2, 2, Base::MeshType::TRIANGULAR);
    test5.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test5.getTotalError() - 0.21870166) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test6(2, 3, 2, Base::MeshType::TRIANGULAR);
    test6.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test6.getTotalError() - 0.02345377) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test7(4, 4, 2, Base::MeshType::TRIANGULAR);
    test7.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test7.getTotalError() - 0.00039351) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test8(8, 5, 2, Base::MeshType::TRIANGULAR);
    test8.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test8.getTotalError() - 0.00000066) < 1e-8), "comparison to old results");
    
    PoissonTest<2> test9(16, 1, 2, Base::MeshType::TRIANGULAR);
    test9.solveSteadyStateWithPetsc(true);
    logger.assert_always((std::abs(test9.getTotalError() - 0.00911139) < 1e-8), "comparison to old results");
    
    //no 3D testing due to speed related issues
    /*
    Laplace test0(1, 2, 1, Base::MeshType::RECTANGULAR);
    test0.initialise();
    std::cout.precision(10);
    std::cout << test0.solve() << std::endl;
    logger.assert_always((std::abs(test0.solve() - 0.48478776) < 1e-8), "comparison to old results");
    Laplace test1(2, 3, 1, Base::MeshType::RECTANGULAR);
    test1.initialise();
    std::cout << test1.solve() << std::endl;
    logger.assert_always((std::abs(test1.solve() - 0.02225892) < 1e-8), "comparison to old results");
    Laplace test2(4, 4, 1, Base::MeshType::RECTANGULAR);
    test2.initialise();
    std::cout << test2.solve() << std::endl;
    logger.assert_always((std::abs(test2.solve() - 0.00008248) < 1e-8), "comparison to old results");
    Laplace test3(8, 5, 1, Base::MeshType::RECTANGULAR);
    test3.initialise();
    std::cout << test3.solve() << std::endl;
    logger.assert_always((std::abs(test3.solve() - 0.00000008) < 1e-8), "comparison to old results");
    Laplace test4(16, 1, 1, Base::MeshType::RECTANGULAR);
    test4.initialise();
    std::cout << test4.solve() << std::endl;
    logger.assert_always((std::abs(test4.solve() - 0.00904309) < 1e-8), "comparison to old results");
    Laplace test5(1, 2, 2, Base::MeshType::TRIANGULAR);
    test5.initialise();
    std::cout << test5.solve() << std::endl;
    logger.assert_always((std::abs(test5.solve() - 0.21870166) < 1e-8), "comparison to old results");
    Laplace test6(2, 3, 2, Base::MeshType::TRIANGULAR);
    test6.initialise();
    std::cout << test6.solve() << std::endl;
    logger.assert_always((std::abs(test6.solve() - 0.02345377) < 1e-8), "comparison to old results");
    Laplace test7(4, 4, 2, Base::MeshType::TRIANGULAR);
    test7.initialise();
    std::cout << test7.solve() << std::endl;
    logger.assert_always((std::abs(test7.solve() - 0.00039351) < 1e-8), "comparison to old results");
    Laplace test8(8, 5, 2, Base::MeshType::TRIANGULAR);
    test8.initialise();
    std::cout << test8.solve() << std::endl;
    logger.assert_always((std::abs(test8.solve() - 0.00000066) < 1e-8), "comparison to old results");
    Laplace test9(16, 1, 2, Base::MeshType::TRIANGULAR);
    test9.initialise();
    std::cout << test9.solve() << std::endl;
    logger.assert_always((std::abs(test9.solve() - 0.00911139) < 1e-8), "comparison to old results");
     */
     
    return 0;
}

