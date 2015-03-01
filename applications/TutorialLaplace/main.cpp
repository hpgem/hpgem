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

#define hpGEM_INCLUDE_PETSC_SUPPORT


#include <fstream>
#include <iostream>

#include "Base/HpgemUISimplified.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"
#include "Output/GNUPlotDiscontinuousSolutionWriter.hpp"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.hpp"
#include "Utilities/BasisFunctions2DH1ConformingSquare.hpp"
#include "Utilities/GlobalMatrix.hpp"
#include "Utilities/GlobalVector.hpp"
#include "petscksp.h"

#include "Base/Norm2.hpp"
#include "Output/TecplotSingleElementWriter.hpp"
#include "Geometry/ReferenceTetrahedron.hpp"
#include "Base/Element.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Base/CommandLineOptions.hpp"

///\brief Make the encapsulation for the DG problem. 
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
class TutorialLaplace : public Base::HpgemUISimplified, Output::SingleElementWriter
{
public:
    ///Constructor: assign the dimension, number of elements and maximum order 
    ///of basisfunctions. Furthermore, assign a value to the penalty parameter.
    ///n stands for number of elements, p stands for polynomial order.
    ///Lastly, construct the mesh with this number of elements and polynomial order.
    TutorialLaplace(int n, int p) : HpgemUISimplified(2, p), n_(n), p_(p)
    {
        DIM_ = 2;
        penalty_ = 3 * n_ * p_ * (p_ + DIM_ - 1) + 1;
        initialise();
    }

    ///\brief set up the mesh  
    ///
    ///In this example, we are going to make a triangular mesh on the domain [0,1]^2
    ///First define the domain, number of elements in each direction and whether or
    ///no there are periodic boundary conditions. Then make a triangular mesh and 
    ///generate the basisfunctions on the reference domain.
    bool initialise() override
    {
        //Make the framework for the mesh
        RectangularMeshDescriptorT description(DIM_);

        for (int i = 0; i < DIM_; ++i)
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
            description.boundaryConditions_[i] = RectangularMeshDescriptorT::SOLID_WALL;
        }
        //Make a triangular mesh from the mesh descriptor. 
        addMesh(description, Base::TRIANGULAR, 1, 1, 1, 1);
        //On this mesh, make the default basisfunctions the standard DG basisfunctions on a triangle.
        meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Triangle(p_));
        return true;
    }
 
    ///\brief Compute the integrals of the left-hand side associated with elements
    ///
    ///For every element, we want to compute the integral of grad(phi_i).grad(phi_j) for all combinations of i and j.
    ///This function will compute the values of gradient(phi_i).gradient(phi_j) for all combinations of i and j for one element
    ///and one reference point p. Then the integral can later be computed with appropriate (Gauss-)quadrature rules.
    ///The resulting matrix of values is then given in the matrix integrandVal, to which we passed a reference when calling it.
    ///Please note that you pass a reference point to the basisfunctions and the 
    ///transformations are done internally.
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& integrandVal) override
    {
        //Obtain the number of basisfunctions that are possibly non-zero on this element.
        const std::size_t numBasisFunctions = element->getNrOfBasisFunctions();

        //Resize the integrandVal such that it contains as many rows and columns as 
        //the number of basisfunctions.
        integrandVal.resize(numBasisFunctions, numBasisFunctions);

        //Initialize the vectors that contain gradient(phi_i) and gradient(phi_j)
        LinearAlgebra::NumericalVector phiDerivI(DIM_), phiDerivJ(DIM_);

        for (std::size_t i = 0; i < numBasisFunctions; ++i)
        {
            //The gradient of basisfunction phi_i is computed at point p, the result is stored in phiDerivI.
            element->basisFunctionDeriv(i, p, phiDerivI);

            for (std::size_t j = 0; j < numBasisFunctions; ++j)
            {
                //The gradient of basisfunction phi_j is computed at point p, the result is stored in phiDerivJ.
                element->basisFunctionDeriv(j, p, phiDerivJ);
                //Compute the value of gradient(phi_i).gradient(phi_j) at point p and 
                //store it at the appropriate place in the matrix integrandVal.
                integrandVal(j, i) = phiDerivI * phiDerivJ;
            }
        }
    }

    /// \brief Compute the integrals of the left-hand side associated with faces.
    ///
    ///For every internal face, we want to compute the integral of 
    /// -1/2 (normal_i phi_i * gradient phi_j) - 1/2 (normal_j phi_j * gradient phi_i) + penalty * (normal_i phi_i * normal_j phi_j)
    ///for all basisfunctions phi_i and phi_j that are non-zero on that face.
    ///For boundary faces, similar expressions can be obtained depending of the type of boundary condition.
    ///This function will compute these integrands for all basisfunctions phi_i and phi_j
    ///on a certain face at a reference point p. Then the integral can later be computed with appropriate (Gauss-)quadrature rules.  
    ///The resulting matrix of values is then given in the matrix integrandVal, to which we passed a reference when calling it.
    ///Please note that you pass a reference point to the basisfunctions and the 
    ///transformations are done internally. If you expect 4 matrices here, 
    ///you can assume that integrandVal is block structured with 4 blocks in total such
    ///that basisfunctions belonging to the left element are on the left and top.
    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p, LinearAlgebra::Matrix& integrandVal) override
    {
        //Obtain the number of basisfunctions that are possibly non-zero at this face.
        const std::size_t numBasisFunctions = face->getNrOfBasisFunctions();

        //Resize the integrandVal such that it contains as many rows and columns as 
        //the number of basisfunctions.
        integrandVal.resize(numBasisFunctions, numBasisFunctions);

        //Initialize the vectors that contain gradient(phi_i), gradient(phi_j), normal_i phi_i and normal_j phi_j
        LinearAlgebra::NumericalVector phiNormalI(DIM_), phiNormalJ(DIM_), phiDerivI(DIM_), phiDerivJ(DIM_);

        //Transform the point from the reference value to its physical value.
        //This is necessary to check at which boundary we are if we are at a boundary face.
        PointPhysicalT pPhys(DIM_);
        face->referenceToPhysical(p, pPhys);

        for (int i = 0; i < numBasisFunctions; ++i)
        {
            //normal_i phi_i is computed at point p, the result is stored in phiNormalI.
            face->basisFunctionNormal(i, normal, p, phiNormalI);
            //The gradient of basisfunction phi_i is computed at point p, the result is stored in phiDerivI.
            face->basisFunctionDeriv(i, p, phiDerivI);

            for (int j = 0; j < numBasisFunctions; ++j)
            {
                //normal_j phi_j is computed at point p, the result is stored in phiNormalJ.
                face->basisFunctionNormal(j, normal, p, phiNormalJ);
                //The gradient of basisfunction phi_j is computed at point p, the result is stored in phiDerivJ.
                face->basisFunctionDeriv(j, p, phiDerivJ);

                //Switch to the correct type of face, and compute the integrand accordingly
                //Internal face:
                if (face->isInternal())
                {
                    integrandVal(j, i) = -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) / 2
                        + penalty_ * phiNormalI * phiNormalJ;
                }
                //Boundary face with Dirichlet boundary conditions:
                else if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1.) < 1e-9)
                {
                    integrandVal(j, i) = -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI)
                        + penalty_ * phiNormalI * phiNormalJ;
                }
                //Boundary face with homogeneous Neumann boundary conditions:
                else
                {
                    integrandVal(j, i) = 0;
                }
            }
        }
    }

    ///\brief Define initial conditions.
    ///
    ///In a non-steady state problem, the initial conditions can be inserted here.
    ///Since we have a steady state problem, there is no initial condition
    double initialConditions(const PointPhysicalT& p) override
    {
        return 0;
    }

    ///\brief Define the source term.
    ///
    ///Define the source, which is the right hand side of laplacian(u) = f(x,y).
    ///Here: f(x,y) = -8pi^2 * sin(2pi x) * sin(2pi y).
    double source(const PointPhysicalT& p)
    {
        return (-8 * M_PI * M_PI) * sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
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
    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p, LinearAlgebra::NumericalVector& integrandVal) override
    {
        //Obtain the number of basisfunctions that are possibly non-zero
        const std::size_t numBasisFunctions = face->getNrOfBasisFunctions();
        //Resize the integrandVal such that it contains as many rows as 
        //the number of basisfunctions.
        integrandVal.resize(numBasisFunctions);

        //Compute the value of the integrand
        //We have no rhs face integrals, so this is just 0.
        for (std::size_t i = 0; i < numBasisFunctions; ++i)
        {
            integrandVal[i] = 0;
        }
    }

    ///\brief Compute the integrals of the right-hand side associated with elements, 
    /// which is in this case interpolating the source term.
    ///
    ///To implement the source term, one needs to integrate over the basisfunction 
    ///times the source term for all basisfunctions. The set of all basisfunctions
    ///is the same as the set of all basisfunctions per element for all elements,
    ///so we can just integrate this for all basisfunctions for all elements.
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::NumericalVector& integrandVal) override
    {
        //Obtain the number of basisfunctions that are possibly non-zero
        const std::size_t numBasisFunctions = element->getNrOfBasisFunctions();
        //Transform reference point p to physical point pPhys to evaluate the source term
        PointPhysicalT pPhys(DIM_);
        element->referenceToPhysical(p, pPhys);

        //Resize the right-hand side vector to get the same amount of rows as the 
        //number of basisFunctions.
        integrandVal.resize(numBasisFunctions);

        //compute value of the integrand, namely basisfunction times source
        for (std::size_t i = 0; i < numBasisFunctions; ++i)
        {
            integrandVal[i] = - element->basisFunction(i, p) * source(pPhys);
        }
    }
    
    ///\brief Write the output to an outputstream (file)
    ///
    ///This function is used by other functions to write the tecplot file, so it is
    ///not a function to write the whole file.
    ///The only thing this has to write in the file is the value of the solution.
    void writeOutput(const ElementT* element, const PointReferenceT& p, std::ostream& out) override
    {
        LinearAlgebra::NumericalVector value(1);
        element->getSolution(0, p, value);
        out << value[0];
    }
    
    /// \brief Solve the system
    ///
    /// This function contains a lot of PETSc code and should be automated in the future.
    /// First compute all integrals and assemble the matrix and right-hand side vector.
    /// Then solve this system with a krylov subspace method from the PETSc package.
    /// Finally, write the solution and mesh to the Tecplot file output.dat.
    bool solve()
    {
        //Compute all element integrals.
        doAllElementIntegration();
        //Compute all face integrals.
        doAllFaceIntegration();
        //Assemble the matrix A of the system Ax = b.
        Utilities::GlobalPetscMatrix A(HpgemUI::meshes_[0], 0, 0);
        //Declare the vectors x and b of the system Ax = b.
        Utilities::GlobalPetscVector b(HpgemUI::meshes_[0], 0, 0), x(HpgemUI::meshes_[0]);

        //Assemble the vector b. This is needed because Petsc assumes you don't know
        //yet whether a vector is a variable or right-hand side the moment it is 
        //declared.
        b.assemble();

        //Make the Krylov supspace method
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        //Tell ksp that it will solve the system Ax = b.
        KSPSetOperators(ksp, A, A);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, b, x);
        //Do PETSc magic, including solving.
        KSPConvergedReason conferge;
        KSPGetConvergedReason(ksp, &conferge);
        int iterations;
        KSPGetIterationNumber(ksp, &iterations);
        std::cout << "KSP solver ended because of " << KSPConvergedReasons[conferge] <<
            " in " << iterations << " iterations." << std::endl;

        x.writeTimeLevelData(0);

        //Write solution to file (can be opened with GNUPlot).
        std::ofstream outFile("output.dat");
        Output::GNUPlotDiscontinuousSolutionWriter writeFunc(outFile, "title", "01", "value");
        writeFunc.write(meshes_[0], this);
        return true;
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
};

auto& numBasisFuns = Base::register_argument<int>('n', "numElems", "number of elements per dimension", true);
auto& p = Base::register_argument<int>('p', "order", "polynomial order of the solution", true);
///Example of using the Laplace class. 
///This implementation asks for commandline input arguments for the number of elements
///in each direction and the polynomial order. Then make the object test, initialise
///it with the input parameters and solve it. PetscInitialize and PetscFinalize are
///necessary since we need to use the library PETSc in the solve routine.
int main(int argc, char **argv)
{
    try
    {
        //read the number of elements and polynomial order from the command line.
        //this will also take care of initialisation and finalisation of required external packages
        Base::parse_options(argc, argv);
        
        //Make the object test with n elements in each direction and polynomial order p.
        TutorialLaplace test(numBasisFuns.getValue(), p.getValue());
        
        //Solve the system.
        test.solve();
        
        return 0;
    }
    catch (const char* e)
    {
        std::cout << e;
    }
}

