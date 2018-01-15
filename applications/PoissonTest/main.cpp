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

#include <fstream>
#include <iostream>
#include <chrono>

#include "Base/HpgemAPILinearSteadyState.h"
#include "Base/ElementCacheData.h"
#include "Base/FaceCacheData.h"
#include "Output/GNUPlotDiscontinuousSolutionWriter.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "petscksp.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Geometry/ReferenceTetrahedron.h"
#include "Base/Element.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Integration/ElementIntegral.h"
#include "Base/CommandLineOptions.h"

    // Choose the dimension (1 or 2 or 3)
const static std::size_t DIM = 3;

///\brief Test application that solves the Poisson equation.
///
///In here, all methods that you use
///for the computations and output are defined.
///We are going to solve the problem
/// delta u(x,y) = g(x,y)     on [0,1]^2
/// u(x,y) = 0                if x = 0 or x = 1
/// du(x,y)/dn = 0            if y = 0 or y = 1
/// g(x,y) = -8 pi^2 sin(2pi x)cos(2pi y)
///The analytical solution of this problem is u(x,y) = sin(2pi x)cos(2pi y)
///The problem is discretised with the Interior Penalty Discontinuous Galerkin method.
///
class PoissonTest : public Base::HpgemAPILinearSteadyState<DIM>
{
public:
    ///Constructor: assign the dimension, number of elements and maximum order
    ///of basisfunctions. Furthermore, assign a value to the penalty parameter.
    ///n stands for number of elements, p stands for polynomial order.
    ///Lastly, construct the mesh with this number of elements and polynomial order.
    PoissonTest(const std::size_t dimension, const std::size_t n, const std::size_t p) :
    Base::HpgemAPILinearSteadyState<DIM>(1, p, true, true),
    n_(n),
    p_(p)
    {
        penalty_ = 3 * n_ * p_ * (p_ + DIM - 1) + 1;
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
    ///
    ///For every element, we want to compute the integral of grad(phi_i).grad(phi_j) for all combinations of i and j.
    ///This function will compute the values of gradient(phi_i).gradient(phi_j) for all combinations of i and j for one element
    ///and one reference point p. Then the integral can later be computed with appropriate (Gauss-)quadrature rules.
    ///The resulting matrix of values is then given in the matrix integrandVal, which we return.
    ///Please note that you pass a reference point to the basisfunctions and the
    ///transformations are done internally.
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(Base::PhysicalElement<DIM>& element) override final
    {
        //Obtain the number of basisfunctions that are possibly non-zero on this element.
        const std::size_t numberOfBasisFunctions = element.getElement()->getNumberOfBasisFunctions();
        
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
    
    /// \brief Compute the integrand for the siffness matrix at the face.
    ///
    ///For every internal face, we want to compute the integral of
    /// -1/2 (normal_i phi_i * gradient phi_j) - 1/2 (normal_j phi_j * gradient phi_i) + penalty * (normal_i phi_i * normal_j phi_j)
    ///for all basisfunctions phi_i and phi_j that are non-zero on that face.
    ///For boundary faces, similar expressions can be obtained depending of the type of boundary condition.
    ///This function will compute these integrands for all basisfunctions phi_i and phi_j
    ///on a certain face at a reference point p. Then the integral can later be computed with appropriate (Gauss-)quadrature rules.
    ///The resulting matrix of values is then given in the matrix integrandVal, which is returned.
    ///Please note that you pass a reference point to the basisfunctions and the
    ///transformations are done internally.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(Base::PhysicalFace<DIM>& face) override final
    {
        //Get the number of basis functions, first of both sides of the face and
        //then only the basis functions associated with the left and right element.
        std::size_t numberOfBasisFunctions = face.getFace()->getNumberOfBasisFunctions();
        
        //Create the FaceMatrix integrandVal with the correct size.
        Base::FaceMatrix& integrandVal = face.getResultMatrix();
        
        //Initialize the vectors that contain gradient(phi_i), gradient(phi_j), normal_i phi_i and normal_j phi_j
        LinearAlgebra::SmallVector<DIM> phiNormalI, phiNormalJ, phiDerivI, phiDerivJ;
        
        //Transform the point from the reference value to its physical value.
        //This is necessary to check at which boundary we are if we are at a boundary face.
        const PointPhysicalT& pPhys = face.getPointPhysical();
        
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i)
        {
            //normal_i phi_i is computed at point p, the result is stored in phiNormalI.
            phiNormalI = face.basisFunctionUnitNormal(i);
            //The gradient of basisfunction phi_i is computed at point p, the result is stored in phiDerivI.
            phiDerivI = face.basisFunctionDeriv(i);
            
            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j)
            {
                //normal_j phi_j is computed at point p, the result is stored in phiNormalJ.
                phiNormalJ = face.basisFunctionUnitNormal(j);
                //The gradient of basisfunction phi_j is computed at point p, the result is stored in phiDerivJ.
                phiDerivJ = face.basisFunctionDeriv(j);
                
                //Switch to the correct type of face, and compute the integrand accordingly
                //you could also compute the integrandVal by directly using face->basisFunctionDeriv
                //and face->basisFunctionNormal in the following lines, but this results in very long expressions
                //Internal face:
                if (face.isInternal())
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
    /// \details In this case the exact solution is u(x,y) = sin(2pi x) * cos(2pi y).
    LinearAlgebra::MiddleSizeVector getExactSolution(const PointPhysicalT &p) override final
    {
        LinearAlgebra::MiddleSizeVector exactSolution(1);
        exactSolution[0] = std::sin(2 * M_PI * p[0]) * std::cos(2 * M_PI * p[1]) * std::cos(2 * M_PI * p[2]);
        return exactSolution;
    }
    
    ///\brief Define the source term.
    ///
    ///Define the source, which is the right hand side of laplacian(u) = f(x,y).
    ///Here: f(x,y) = -8pi^2 * sin(2pi x) * sin(2pi y).
    LinearAlgebra::MiddleSizeVector getSourceTerm(const PointPhysicalT &p) override final
    {
        LinearAlgebra::MiddleSizeVector sourceTerm(1);
        sourceTerm[0] = (-12 * M_PI * M_PI) * std::sin(2 * M_PI * p[0]) * std::cos(2 * M_PI * p[1]) * std::cos(2 * M_PI * p[2]);
        //sourceTerm[0] = 10 * std::exp(-((p[0]-0.5) * (p[0]-0.5) + (p[1]-0.5) * (p[1]-0.5))/0.02);
        return sourceTerm;
    }
    
    //For the right hand side, we also need to integrate over elements and faces.
    //This will be done in the two functions below.
    //After they are computed, we have a vector, each element representing one basisfunction.
    
    /// \brief Compute the integrals of the right-hand side associated with faces.
    ///
    ///To implement the boundary conditions, one often needs to compute the face
    ///integral on the right-hand side. However, in our application we do not have
    ///contributions for the boundary conditions, so the vector has only zeroes.
    ///The input/output structure is the same as the other faceIntegrand function.
    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(Base::PhysicalFace<DIM>& face) override final
    {
        //Obtain the number of basisfunctions that are possibly non-zero
        const std::size_t numberOfBasisFunctions = face.getFace()->getNumberOfBasisFunctions();
        auto pPhys = face.getPointPhysical();
        //Resize the integrandVal such that it contains as many rows as
        //the number of basisfunctions.
        LinearAlgebra::MiddleSizeVector integrandVal(numberOfBasisFunctions);
        
        //Compute the value of the integrand
        //We have no rhs face integrals, so this is just 0.
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i)
        {
            if (std::abs(pPhys[0] + 1.e3) < 1e-9)
            {
                integrandVal[i] = -1.;
            }
            else if (std::abs(pPhys[0] - 1.e3) < 1e-9)
            {
                integrandVal[i] = 1.;
            }
            else
            {
                integrandVal[i] = 0;
            }
        }
        
        return integrandVal;
    }

    void tasksBeforeSolving() override final
    {
        Base::HpgemAPILinearSteadyState<DIM>::tasksBeforeSolving();
        std::function<LinearAlgebra::SmallVector<DIM>(Base::Element*, const Geometry::PointReference<DIM>&, std::size_t)> newFunction =
                [](Base::Element* element, const Geometry::PointReference<DIM>& point, std::size_t timeLevel)->LinearAlgebra::SmallVector<DIM>
                {
                    if(std::isnan(Base::L2Norm(element->getSolutionGradient(timeLevel, point)[0])))
                    {
                        logger(INFO, "%", element->calcJacobian(point));
                        logger(ERROR, "nan gradient detected at element % (type %) at point %", element->getID(), element->getReferenceGeometry()->getName(), point);
                    }
                    return element->getSolutionGradient(timeLevel, point)[0];
                };
        registerVTKWriteFunction(newFunction, "gradient");
    }
    
    ///\brief Write the output to an outputstream (file)
    ///
    ///This function is used by other functions to write the tecplot file, so it is
    ///not a function to write the whole file.
    ///The only thing this has to write in the file is the value of the solution.
    void writeToTecplotFile(const Base::Element* element, const PointReferenceT& point, std::ostream& out) override final
    {
        LinearAlgebra::MiddleSizeVector value(1);
        value = element->getSolution(0, point);
        out << value[0];
    }

        
private:
    
    ///number of elements per cardinal direction
    std::size_t n_;
    
    ///polynomial order of the approximation
    std::size_t p_;
    
    ///\brief Penalty parameter
    ///
    ///Penalty parameter that is associated with the interior penalty discontinuous
    ///Galerkin method. This parameter is initialized in the constructor, and has
    ///to be greater than 3 * n_ * p_ * (p_ + DIM - 1) in order for the method to be stable.
    double penalty_;
};

auto& meshName = Base::register_argument<std::string>('\0', "meshName", "name of the Mesh file", true);
auto& numberOfElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", true);
auto& p = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);
///Example of using the Laplace class.
///This implementation asks for commandline input arguments for the number of elements
///in each direction and the polynomial order. Then make the object test, initialise
///it with the input parameters and solve it. PetscInitialize and PetscFinalize are
///necessary since we need to use the library PETSc in the solve routine.
int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);

    // Choose a mesh type (e.g. TRIANGULAR, RECTANGULAR).
    const Base::MeshType meshType = Base::MeshType::TRIANGULAR;

    // Choose variable name(s). Since we have a scalar function, we only need to chooes one name.
    std::vector<std::string> variableNames;
    variableNames.push_back("u");

    //Make the object test with n elements in each direction and polynomial order p.
    PoissonTest test(DIM, numberOfElements.getValue(), p.getValue());

    //Create the mesh
    test.readMesh(meshName.getValue());

    // Set the names for the output file
    test.setOutputNames("output", "PoissonTest", "PoissonTest", variableNames);

    auto startTime = std::chrono::steady_clock::now();

    //Solve the system.
    test.solveSteadyStateWithPetsc(true);

    auto endTime = std::chrono::steady_clock::now();

    logger(INFO, "Simulation took %ms.", std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count());

    return 0;
}
        
